"""
Experiment 7 — Sub-algebra Gate: group-Lasso recovery of su(3) ⊂ su(4)
========================================================================

QUESTION
--------
Does the combination of Killing crystallization (T_K) + orthogonal confinement
(T_ortho) + group-Lasso gate spontaneously identify the su(3) sub-algebra
embedded in a larger su(4) input space?

TASK DESIGN
-----------
Input  : x ∈ R^{15}  (full su(4) GML coordinates, 15 = 4²-1 generators)
Target : y = sign(Re(Tr(M_sub³)))
         M_sub = Σ_{a ∈ SU3_GENS} x_a · T_a^{su(4)}
         — the su(3)-cubic-sign task, but embedded in the larger su(4) space.

Generator index map  (su(4) GML basis, 0-indexed, from _su_basis(4)):
  Group 1 — Symmetric off-diagonal (i < j):
    (0,1)=0    (0,2)=1    (0,3)=2*   (1,2)=3    (1,3)=4*   (2,3)=5*
  Group 2 — Antisymmetric off-diagonal (i < j):
    (0,1)=6    (0,2)=7    (0,3)=8*   (1,2)=9    (1,3)=10*  (2,3)=11*
  Group 3 — Diagonal Cartan:
    H_1=12     H_2=13     H_3=14*
  (* generator involves row/col 3 → NOT in the su(3) subalgebra)

  su(3) signal generators  : {0, 1, 3, 6, 7, 9, 12, 13}   — 8 total
  su(4)-only noise gens    : {2, 4, 5, 8, 10, 11, 14}       — 7 total

  Verification: T_a^{su4} for a ∈ {0,1,3,6,7,9,12,13} are exactly the
  su(3) Gell-Mann generators padded with zeros in row/col 3.  Their 3×3
  upper-left blocks are identical to _su_basis(3) in the order
  [0,1,2,3,4,5,6,7] → [0,1,3,6,7,9,12,13].

MODELS
------
  1. MLP            : standard 15→hidden→2, no algebraic structure
  2. AlgNet-NoGate  : T_K + T_ortho (λ_K=2x, λ_o=1.0), no gate, 1500 ep
  3. AlgNet-Gate-lo : T_K + T_ortho + group-Lasso gate (λ_gate=0.05)
  4. AlgNet-Gate-med: T_K + T_ortho + group-Lasso gate (λ_gate=0.15)
  5. Oracle         : theta=I (fixed), gate=su3_mask (fixed), head trains only

GROUP-LASSO GATE
----------------
  gate(g_raw) = softplus(g_raw)   —  smooth, always positive
  T_gate = Σ_a gate_a             —  L1-type sparsity penalty
  Loss = L_CE + λ_K T_K + λ_o T_ortho + λ_gate T_gate

  Initialisation: g_raw = GATE_INIT (scalar) → gate starts at softplus(GATE_INIT)
  At optimum: gate a ≈ 0 for noise gens (task gradient is zero, L1 pushes down)
              gate a > 0 for signal gens (task gradient resists collapse)

NEW DIAGNOSTICS
---------------
  GateSparsity   : fraction of gate values < 0.1
  GateAlignScore : |top-8 by gate ∩ SU3_GENS_IDX| / 8   (1.0 = perfect)
  GateSignal     : mean gate value at su(3) positions
  GateNoise      : mean gate value at su(4)-only positions
  SignalNoiseRatio: GateSignal / max(GateNoise, 1e-6)

PHYSICAL FRAMING
----------------
  The group-Lasso gate is the neural analog of the 'isospin projection'
  in su(2) ⊂ su(3): a variational handle on which generators carry the
  task-relevant invariant.  SS-PEG postulates that Killing crystallization
  provides exactly this handle because T_K = 0 selects a valid algebra basis
  within which the task-relevant subalgebra can be isolated.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
#  Path setup
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from killing_utils import (
    killing_tension,
    _su_basis,
    params_to_su_generators,
)

# ---------------------------------------------------------------------------
#  Devices
# ---------------------------------------------------------------------------
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# ---------------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------------
D_FULL    = 4
N_GEN     = D_FULL * D_FULL - 1          # 15
D_SUB     = 3
N_GEN_SUB = D_SUB * D_SUB - 1            # 8

# su(3) generator indices inside _su_basis(4) — see docstring for derivation
SU3_GENS_IDX   = [0, 1, 3, 6, 7, 9, 12, 13]
SU4_EXTRA_IDX  = [2, 4, 5, 8, 10, 11, 14]

LAMBDA_K_BASE  = 2.0
N_GEN_REF      = 8
LK_EFF         = LAMBDA_K_BASE * N_GEN_REF / N_GEN   # ≈ 1.0667
LK_2X          = 2.0 * LK_EFF                         # ≈ 2.1333

LAMBDA_ORTHO   = 1.0      # single-phase, per Exp6c finding
LAMBDA_GATE_LO = 0.05
LAMBDA_GATE_MED= 0.15
GATE_INIT      = 1.5      # init: gate = softplus(1.5) ≈ 1.69

N_EPOCHS       = 1500
REPORT_EVERY   = 300
BATCH_SIZE     = 256
N_TRAIN        = 6_000
N_VAL          = 1_500
SEED           = 42
LR_MAIN        = 1e-3

OUT_DIR = os.path.join(HERE, "exp7_outputs")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
#  Tension functions
# ---------------------------------------------------------------------------

def ortho_tension(theta: torch.Tensor) -> torch.Tensor:
    """T_ortho = ||theta @ theta.T - I||_F^2  (manifold confinement)."""
    gram = theta @ theta.T
    I_n  = torch.eye(gram.shape[0], dtype=theta.dtype, device=theta.device)
    return ((gram - I_n) ** 2).sum()


def lasso_tension(gate_values: torch.Tensor) -> torch.Tensor:
    """L1 gate penalty: Σ_a gate_a  (promotes sparsity)."""
    return gate_values.sum()


# ---------------------------------------------------------------------------
#  Dataset
# ---------------------------------------------------------------------------

def generate_su4_su3task_dataset(n: int, seed: int = SEED):
    """
    x ∈ R^{15}  (su(4) GML coordinates).
    y = sign(Re(Tr(M_sub^3)))  where M_sub uses only su(3)-block generators.

    Coordinates in SU4_EXTRA_IDX are drawn from the same distribution as
    SU3_GENS_IDX but are entirely irrelevant to the target — pure noise.
    """
    rng  = np.random.default_rng(seed)
    X    = rng.standard_normal((n, N_GEN)).astype(np.float64)

    basis_su4 = _su_basis(D_FULL, "cpu").numpy()              # (15, 4, 4) complex
    basis_sub = basis_su4[SU3_GENS_IDX]                        # (8, 4, 4) su(3) block

    # Project only the su(3) coordinates
    X_sub = X[:, SU3_GENS_IDX]                                 # (n, 8)
    M_sub = np.einsum("na,aij->nij", X_sub, basis_sub)         # (n, 4, 4) complex
    M3    = M_sub @ M_sub @ M_sub
    tr3   = np.trace(M3, axis1=-2, axis2=-1).real

    y   = (tr3 > 0).astype(np.int64)
    pos = y.sum()
    print(f"  su({D_SUB})⊂su({D_FULL}) dataset  n={n}  "
          f"balance: {n-pos} neg / {pos} pos  (frac={pos/n:.3f})")
    return (
        torch.tensor(X,  dtype=torch.float64, device=DEVICE),
        torch.tensor(y,  dtype=torch.long,    device=DEVICE),
    )


def make_loaders():
    X_tr, y_tr = generate_su4_su3task_dataset(N_TRAIN, seed=SEED)
    X_va, y_va = generate_su4_su3task_dataset(N_VAL,   seed=SEED + 1)
    tr = DataLoader(TensorDataset(X_tr, y_tr), batch_size=BATCH_SIZE, shuffle=True)
    va = DataLoader(TensorDataset(X_va, y_va), batch_size=512,         shuffle=False)
    return tr, va


# ---------------------------------------------------------------------------
#  Common head
# ---------------------------------------------------------------------------

def make_head(n_in: int) -> nn.Sequential:
    hidden = min(4 * n_in, 128)
    return nn.Sequential(
        nn.ReLU(),
        nn.Linear(n_in,   hidden, dtype=torch.float64),
        nn.ReLU(),
        nn.Linear(hidden, 2,     dtype=torch.float64),
    )


# ---------------------------------------------------------------------------
#  Models
# ---------------------------------------------------------------------------

class AlgebraicProjection(nn.Module):
    """
    Linear algebraic projection layer.
    G_a = Σ_k theta_{ak} T_k   (linear combination of GML generators)
    Forward: f_a = Re(Tr(G_a M(x)))  for each sample x.

    Regularization (T_K, T_ortho) is applied externally in the training loop.
    """

    def __init__(self, d: int, fixed: bool = False):
        super().__init__()
        self.d          = d
        self.n_gen_full = d * d - 1
        self.fixed      = fixed

        init = torch.eye(self.n_gen_full, dtype=torch.float64)
        if fixed:
            self.register_buffer("theta_buf", init)
        else:
            noise = torch.randn_like(init) * 0.05
            self.theta = nn.Parameter(init + noise)
        self.bias = nn.Parameter(torch.zeros(self.n_gen_full, dtype=torch.float64))

    def _theta(self) -> torch.Tensor:
        return self.theta_buf if self.fixed else self.theta

    def get_generators(self) -> torch.Tensor:
        """(n_gen, d, d) complex generators."""
        return params_to_su_generators(self._theta(), self.d)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        dev   = str(self._theta().device)
        basis = _su_basis(self.d, dev)
        M_x   = torch.einsum("...b,bij->...ij", x.to(torch.complex128), basis)
        G     = self.get_generators()
        f     = torch.einsum("aij,...ji->...a", G, M_x).real
        return f + self.bias


class AlgebraicNetGated(nn.Module):
    """
    AlgebraicProjection + learnable element-wise gate + classification head.

    gate(g_raw) = softplus(g_raw) applied element-wise to the feature vector.
    Sparsity promotion: add  λ_gate * gate.sum()  to the loss externally.

    If fixed_gate is provided (a buffer Tensor), the gate is fixed (Oracle mode).
    """

    def __init__(self, d: int, fixed: bool = False, fixed_gate: torch.Tensor | None = None):
        super().__init__()
        self.d          = d
        self.n_gen      = d * d - 1
        self.alg_layer  = AlgebraicProjection(d=d, fixed=fixed)
        self.head       = make_head(self.n_gen)

        if fixed_gate is not None:
            self.register_buffer("gate_buf", fixed_gate.to(torch.float64))
            self._has_learnable_gate = False
        else:
            self.gate_raw = nn.Parameter(
                torch.full((self.n_gen,), GATE_INIT, dtype=torch.float64)
            )
            self._has_learnable_gate = True

    def gate(self) -> torch.Tensor:
        if self._has_learnable_gate:
            return F.softplus(self.gate_raw)
        return self.gate_buf

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        f      = self.alg_layer(x)          # (batch, n_gen)
        gated  = f * self.gate()            # element-wise
        return self.head(gated)

    def get_generators(self) -> torch.Tensor:
        with torch.no_grad():
            return self.alg_layer.get_generators()


class StandardMLP(nn.Module):
    """Baseline: linear layer (no algebraic structure) + head."""

    def __init__(self):
        super().__init__()
        self.layer1 = nn.Linear(N_GEN, N_GEN, dtype=torch.float64)
        self.head   = make_head(N_GEN)

    def forward(self, x):
        return self.head(self.layer1(x))


# ---------------------------------------------------------------------------
#  Evaluation
# ---------------------------------------------------------------------------

def _evaluate(model: nn.Module, val_loader) -> float:
    model.eval()
    correct, total = 0, 0
    with torch.no_grad():
        for xb, yb in val_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            correct += (model(xb).argmax(1) == yb).sum().item()
            total   += xb.size(0)
    return correct / total


def _full_diagnostics(model: nn.Module):
    """
    Returns (T_K, OrthoErr, GateSparsity, GateAlignScore, GateSignal, GateNoise).
    nan for fields not applicable to the given model.
    """
    # Algebraic (T_K, OrthoErr) — only for non-fixed AlgebraicProjection
    has_alg = hasattr(model, "alg_layer") and not model.alg_layer.fixed
    if has_alg:
        with torch.no_grad():
            G     = model.alg_layer.get_generators()
            theta = model.alg_layer.theta
            t_k   = killing_tension(G).item()
            th    = theta.detach().cpu().double().numpy()
            o_err = float(np.linalg.norm(th @ th.T - np.eye(th.shape[0]), "fro"))
    else:
        t_k, o_err = float("nan"), float("nan")

    # Gate diagnostics
    has_gate = hasattr(model, "_has_learnable_gate")
    if has_gate:
        with torch.no_grad():
            g = model.gate().cpu().numpy()
        sparsity  = float((g < 0.1).mean())
        top_k     = set(np.argsort(-g)[:N_GEN_SUB])
        align     = len(top_k & set(SU3_GENS_IDX)) / N_GEN_SUB
        g_signal  = float(g[SU3_GENS_IDX].mean())
        g_noise   = float(g[SU4_EXTRA_IDX].mean())
    else:
        sparsity, align, g_signal, g_noise = (float("nan"),) * 4

    return t_k, o_err, sparsity, align, g_signal, g_noise


# ---------------------------------------------------------------------------
#  Training
# ---------------------------------------------------------------------------

def train_one_model(
    model        : nn.Module,
    train_loader,
    val_loader,
    name         : str,
    lambda_k     : float = LK_2X,
    lambda_ortho : float = LAMBDA_ORTHO,
    lambda_gate  : float = 0.0,
    n_epochs     : int   = N_EPOCHS,
    report_every : int   = REPORT_EVERY,
) -> dict:
    """
    Single-phase training with optional T_K, T_ortho, and gate-Lasso penalties.
    All three penalties support the same warmup (30 ep) for numerical stability.
    """
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=LR_MAIN,
    )
    WARMUP = 30

    hist = {
        "val_acc"   : [],
        "T_K"       : [],
        "ortho_err" : [],
        "sparsity"  : [],
        "align"     : [],
        "g_signal"  : [],
        "g_noise"   : [],
    }
    best_val = 0.0

    print(f"\n  >> {name} ({n_epochs}ep  lK={lambda_k:.4f}  lo={lambda_ortho:.4f}"
          f"  lg={lambda_gate:.4f})")

    for epoch in range(n_epochs):
        warmup_scale = min(1.0, (epoch + 1) / WARMUP)
        lk_ep = lambda_k     * warmup_scale
        lo_ep = lambda_ortho * warmup_scale
        lg_ep = lambda_gate  * warmup_scale

        model.train()
        for xb, yb in train_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            optimizer.zero_grad()
            logits = model(xb)
            loss   = criterion(logits, yb)

            # Algebraic regularization (applicable only to non-fixed AlgNetGated)
            if hasattr(model, "alg_layer") and not model.alg_layer.fixed:
                G     = model.alg_layer.get_generators()
                theta = model.alg_layer.theta
                if lk_ep > 0:
                    loss = loss + lk_ep * killing_tension(G)
                if lo_ep > 0:
                    loss = loss + lo_ep * ortho_tension(theta)

            # Gate Lasso penalty (only when gate is learnable)
            if lg_ep > 0 and hasattr(model, "gate_raw"):
                loss = loss + lg_ep * lasso_tension(F.softplus(model.gate_raw))

            loss.backward()
            optimizer.step()

        val_acc  = _evaluate(model, val_loader)
        best_val = max(best_val, val_acc)
        t_k, o_err, spars, align, g_sig, g_noi = _full_diagnostics(model)

        hist["val_acc"].append(val_acc)
        hist["T_K"].append(t_k)
        hist["ortho_err"].append(o_err)
        hist["sparsity"].append(spars)
        hist["align"].append(align)
        hist["g_signal"].append(g_sig)
        hist["g_noise"].append(g_noi)

        if (epoch + 1) % report_every == 0:
            align_str = f"  Align={align:.3f}" if not np.isnan(align) else ""
            spars_str = f"  Spar={spars:.3f}" if not np.isnan(spars) else ""
            tk_str    = f"  T_K={t_k:.5f}" if not np.isnan(t_k) else ""
            print(f"    [{name} ep{epoch+1:4d}]  val={val_acc:.4f}"
                  f"{tk_str:14s}{spars_str}{align_str}")

    hist["best_val"]      = best_val
    hist["total_epochs"]  = n_epochs
    return hist


# ---------------------------------------------------------------------------
#  Main experiment runner
# ---------------------------------------------------------------------------

def run_exp7(tr, va) -> dict:
    """Train all 5 Exp7 models and return results dict."""
    results = {}

    # ------------------------------------------------------------------
    # 1. MLP baseline
    # ------------------------------------------------------------------
    model = StandardMLP().to(DEVICE)
    results["MLP"] = train_one_model(
        model, tr, va, name="MLP",
        lambda_k=0.0, lambda_ortho=0.0, lambda_gate=0.0,
    )

    # ------------------------------------------------------------------
    # 2. AlgNet-NoGate  (T_K_2x + T_ortho, no gate)
    # ------------------------------------------------------------------
    model = AlgebraicNetGated(d=D_FULL, fixed=False).to(DEVICE)
    # Disable the gate penalty by setting lambda_gate=0 (gate is still present
    # in the architecture but not penalized — its gradient comes from task loss only)
    results["AlgNet-NoGate"] = train_one_model(
        model, tr, va, name="AlgNet-NoGate",
        lambda_k=LK_2X, lambda_ortho=LAMBDA_ORTHO, lambda_gate=0.0,
    )

    # ------------------------------------------------------------------
    # 3. AlgNet-Gate-lo  (λ_gate=0.05)
    # ------------------------------------------------------------------
    model = AlgebraicNetGated(d=D_FULL, fixed=False).to(DEVICE)
    results["AlgNet-Gate-lo"] = train_one_model(
        model, tr, va, name="AlgNet-Gate-lo",
        lambda_k=LK_2X, lambda_ortho=LAMBDA_ORTHO, lambda_gate=LAMBDA_GATE_LO,
    )

    # ------------------------------------------------------------------
    # 4. AlgNet-Gate-med  (λ_gate=0.15)
    # ------------------------------------------------------------------
    model = AlgebraicNetGated(d=D_FULL, fixed=False).to(DEVICE)
    results["AlgNet-Gate-med"] = train_one_model(
        model, tr, va, name="AlgNet-Gate-med",
        lambda_k=LK_2X, lambda_ortho=LAMBDA_ORTHO, lambda_gate=LAMBDA_GATE_MED,
    )

    # ------------------------------------------------------------------
    # 5. Oracle  (theta=I fixed, gate=su3_mask fixed, head trains)
    # ------------------------------------------------------------------
    oracle_gate = torch.zeros(N_GEN, dtype=torch.float64)
    oracle_gate[SU3_GENS_IDX] = 1.0
    model = AlgebraicNetGated(
        d=D_FULL, fixed=True, fixed_gate=oracle_gate
    ).to(DEVICE)
    results["Oracle"] = train_one_model(
        model, tr, va, name="Oracle",
        lambda_k=0.0, lambda_ortho=0.0, lambda_gate=0.0,
    )

    return results


# ---------------------------------------------------------------------------
#  Summary printing
# ---------------------------------------------------------------------------

def print_summary(results: dict) -> None:
    oracle_best = results.get("Oracle", {}).get("best_val", float("nan"))

    print(f"\n{'='*80}")
    print(f"  EXP7 SUMMARY  (d={D_FULL}, su({D_SUB})⊂su({D_FULL}) task)")
    print(f"{'='*80}")
    print(f"  Oracle best: {oracle_best:.4f}")
    print()
    header = (
        f"  {'Model':<22}  {'Best':>6}  {'Final':>6}  "
        f"{'T_K':>9}  {'OrthoErr':>8}  {'Sparsity':>8}  "
        f"{'AlignScore':>10}  {'Sig/Noi':>8}  {'Δ_Oracle':>8}"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))

    for name, h in results.items():
        best  = h["best_val"]
        final = h["val_acc"][-1]
        t_k   = h["T_K"][-1]
        o_err = h["ortho_err"][-1]
        spars = h["sparsity"][-1]
        align = h["align"][-1]
        g_sig = h["g_signal"][-1]
        g_noi = h["g_noise"][-1]
        delta = best - oracle_best

        def _fmt(v, fmt):
            return f"{v:{fmt}}" if not (isinstance(v, float) and np.isnan(v)) else "     --"

        ratio_str = (
            f"{g_sig / max(g_noi, 1e-6):7.2f}x"
            if not (np.isnan(g_sig) or np.isnan(g_noi))
            else "      --"
        )
        print(f"  {name:<22}  {best:.4f}  {final:.4f}  "
              f"{_fmt(t_k, '9.5f')}  {_fmt(o_err, '8.3f')}  "
              f"{_fmt(spars, '8.3f')}  {_fmt(align, '10.3f')}  "
              f"{ratio_str}  {delta:+.4f}")


# ---------------------------------------------------------------------------
#  Plotting
# ---------------------------------------------------------------------------

def plot_results(results: dict) -> None:
    eps     = np.arange(1, N_EPOCHS + 1)
    COLORS  = {
        "MLP"            : "tab:blue",
        "AlgNet-NoGate"  : "tab:orange",
        "AlgNet-Gate-lo" : "tab:green",
        "AlgNet-Gate-med": "tab:red",
        "Oracle"         : "tab:purple",
    }
    DASHES  = {
        "MLP"            : "--",
        "AlgNet-NoGate"  : "-.",
        "AlgNet-Gate-lo" : "-",
        "AlgNet-Gate-med": "-",
        "Oracle"         : ":",
    }

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle(
        f"Exp7 — Sub-algebra Gate: su({D_SUB})⊂su({D_FULL}) recovery\n"
        f"d={D_FULL}, n_gen={N_GEN}, signal gens={N_GEN_SUB}  "
        f"λ_K={LK_2X:.3f}  λ_o={LAMBDA_ORTHO}  "
        f"λ_gate_lo={LAMBDA_GATE_LO}  λ_gate_med={LAMBDA_GATE_MED}",
        fontsize=11,
    )

    def _safe(arr):
        return [np.nan if isinstance(v, float) and np.isnan(v) else v for v in arr]

    # ── (0,0)  Validation accuracy ──────────────────────────────────────────
    ax = axes[0, 0]
    for name, h in results.items():
        ax.plot(eps, _safe(h["val_acc"]),
                color=COLORS[name], ls=DASHES[name], label=name, alpha=0.85)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Val accuracy")
    ax.set_title("Validation accuracy")
    ax.legend(fontsize=7)
    ax.set_ylim(0.5, 0.85)
    ax.grid(True, alpha=0.3)

    # ── (0,1)  Killing tension T_K ──────────────────────────────────────────
    ax = axes[0, 1]
    for name, h in results.items():
        tk = _safe(h["T_K"])
        if all(np.isnan(v) for v in tk):
            continue
        ax.semilogy(eps, tk, color=COLORS[name], ls=DASHES[name], label=name, alpha=0.85)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("T_K  (log scale)")
    ax.set_title("Killing tension T_K")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # ── (0,2)  Ortho error ──────────────────────────────────────────────────
    ax = axes[0, 2]
    for name, h in results.items():
        oe = _safe(h["ortho_err"])
        if all(np.isnan(v) for v in oe):
            continue
        ax.plot(eps, oe, color=COLORS[name], ls=DASHES[name], label=name, alpha=0.85)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("OrthoErr  ||θθᵀ - I||_F")
    ax.set_title("Orthogonal manifold error")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # ── (1,0)  Gate alignment score ─────────────────────────────────────────
    ax = axes[1, 0]
    ax.axhline(1.0, color="k", ls=":", lw=1, alpha=0.5, label="Perfect align")
    for name, h in results.items():
        al = _safe(h["align"])
        if all(np.isnan(v) for v in al):
            continue
        ax.plot(eps, al, color=COLORS[name], ls=DASHES[name], label=name, alpha=0.85)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("GateAlignScore")
    ax.set_title(f"Gate alignment (top-{N_GEN_SUB} vs true su({D_SUB}) gens)")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # ── (1,1)  Gate sparsity ─────────────────────────────────────────────────
    ax = axes[1, 1]
    ax.axhline(7 / N_GEN, color="k", ls=":", lw=1, alpha=0.5, label="Target (7/15)")
    for name, h in results.items():
        sp = _safe(h["sparsity"])
        if all(np.isnan(v) for v in sp):
            continue
        ax.plot(eps, sp, color=COLORS[name], ls=DASHES[name], label=name, alpha=0.85)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("GateSparsity  (frac < 0.1)")
    ax.set_title("Gate sparsity (target = 7/15 ≈ 0.47 noise gens)")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # ── (1,2)  Final gate values bar chart ──────────────────────────────────
    ax = axes[1, 2]
    x_pos = np.arange(N_GEN)
    offset_map = {
        "AlgNet-NoGate"  : -1.5 * 0.15,
        "AlgNet-Gate-lo" : -0.5 * 0.15,
        "AlgNet-Gate-med":  0.5 * 0.15,
        "Oracle"         :  1.5 * 0.15,
    }
    width = 0.14
    for name, h in results.items():
        if name not in offset_map:
            continue
        model_results = h
        # reconstruct final gate from last stored align/sparsity (approximate);
        # we store values per epoch, so retrieve from the last evaluation
        # (actual gate values are NOT stored per epoch; use the model directly)
        # → We record g_signal and g_noise, but not per-generator values.
        # For the bar chart, we use the "align" trace end + signal/noise means.
        # (Full per-generator bars would require storing the full gate each epoch.)
        # Approximate display: shade bars by whether they are in SU3_GENS_IDX.
        pass

    # Simplified bar chart: show SU3 vs noise generator positions with
    # GateSignal / GateNoise means for each model at final epoch.
    bar_data = []
    for name in ["AlgNet-NoGate", "AlgNet-Gate-lo", "AlgNet-Gate-med", "Oracle"]:
        if name not in results:
            continue
        h      = results[name]
        g_sig  = h["g_signal"][-1]
        g_noi  = h["g_noise"][-1]
        bar_data.append((name, g_sig, g_noi))

    categories = [b[0] for b in bar_data]
    sig_vals   = [b[1] for b in bar_data]
    noi_vals   = [b[2] for b in bar_data]
    x = np.arange(len(categories))
    ax.bar(x - 0.2, sig_vals, 0.4, label=f"su({D_SUB}) signal gens mean", color="tab:green",  alpha=0.8)
    ax.bar(x + 0.2, noi_vals, 0.4, label="su(4)-only noise gens mean",    color="tab:red",    alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([b.replace("AlgNet-", "") for b in categories], fontsize=8, rotation=15)
    ax.set_ylabel("Mean gate value")
    ax.set_title(f"Final gate: signal vs noise  (ideal: signal≫noise)")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3, axis="y")

    fig.tight_layout()
    out_path = os.path.join(OUT_DIR, "exp7_results.png")
    fig.savefig(out_path, dpi=150)
    print(f"\n  Figure saved -> {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
#  Save numerics
# ---------------------------------------------------------------------------

def save_npz(results: dict) -> None:
    npz_data = {}
    for name, h in results.items():
        key = name.replace("-", "_").replace(" ", "_")
        for metric, arr in h.items():
            if isinstance(arr, list):
                npz_data[f"{key}__{metric}"] = np.array(arr, dtype=np.float64)
            else:
                npz_data[f"{key}__{metric}"] = np.float64(arr)
    out_path = os.path.join(OUT_DIR, "exp7_results.npz")
    np.savez(out_path, **npz_data)
    print(f"  Numerics saved  -> {out_path}")


# ---------------------------------------------------------------------------
#  Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    t0 = time.time()

    gpu_str = (
        f"GPU    : {torch.cuda.get_device_name(0)}"
        if torch.cuda.is_available() else "CPU mode"
    )
    print(f"""
#################################################################
  Experiment 7 — Sub-algebra Gate  su({D_SUB})⊂su({D_FULL})
  Device : {DEVICE}
  {gpu_str}
  n_gen_full={N_GEN}  n_gen_sub={N_GEN_SUB}
  Signal gens : {SU3_GENS_IDX}
  Noise  gens : {SU4_EXTRA_IDX}
  λ_K_eff={LK_EFF:.4f}  λ_K_2x={LK_2X:.4f}
  λ_ortho={LAMBDA_ORTHO}  λ_gate_lo={LAMBDA_GATE_LO}  λ_gate_med={LAMBDA_GATE_MED}
  Epochs : {N_EPOCHS}
#################################################################
""")

    print("  Preparing datasets...")
    tr, va = make_loaders()

    print(f"\n{'='*60}")
    print(f"  EXP 7 — Training models")
    print(f"{'='*60}")

    results = run_exp7(tr, va)

    print_summary(results)
    plot_results(results)
    save_npz(results)

    print(f"\n  Total wall time: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
