"""
Experiment 6 -- T_ortho Annealing Schedule
============================================================
Physical analogy: G2/F4 closure-tension confinement -> neural manifold confinement.

MOTIVATION
----------
Exp5E revealed two persistent failure modes for d > 3:
  d=4 : AlgNet-TK underperforms MLP despite T_K -> 0
        (root cause: unnormalized lambda_K, fixed in Part E, advantage = +0.009)
  d=5 : T_K = 0.0019,  OrthoErr = 54.4
        theta pulled off O(n_gen) manifold by task gradient -> unstable crystallization

Physical reference: SS-PEG verification of G2 and F4 required TWO terms beyond T_K:
  (a) Closure tension  T_closure: penalizes [G_a, G_b] escaping span(generators)
  (b) Annealing schedule: gamma_1 = 50  (Phase 1), gamma_2 = 200 (Phase 2)

Neural structural analog:
  T_ortho = || theta @ theta.T - I_{n_gen} ||_F^2

  T_ortho confines theta to O(n_gen): the manifold of orthonormal rotations of
  the GML basis. Any theta in O(n_gen) gives a valid su(d) frame with T_K = 0.
  When the task gradient pulls theta off O(n_gen), the Killing metric ceases to
  correctly reflect the algebra structure -> T_K cannot relax to zero.

TRAINING PROTOCOL (3-phase)
---------------------------
Phase 1 (600ep, lr=1e-3):  lambda_K = lambda_eff(d),     lambda_ortho = 0.1
                            weak manifold confinement while task is learned

Phase 2 (200ep, lr=1e-4):  lambda_K = 5 * lambda_eff(d), lambda_ortho = 2.0
                            annealing -- freeze head, strong confinement
                            -> theta retracts to O(n_gen)

Phase 3 (200ep, lr=1e-3):  lambda_K = lambda_eff(d),     lambda_ortho = 0.5
                            fine-tune with head unfrozen at moderate confinement

MODELS COMPARED (d = 4 and d = 5)
-----------------------------------
  1. GML-Fixed       : oracle (theta = I, no training)
  2. AlgNet-NoOrtho  : Part E baseline (T_K only, normalized lambda_K, 1000ep)
  3. AlgNet-Ortho1   : T_K + T_ortho = 1.0 (single phase, 1000ep, no annealing)
  4. AlgNet-OrthoAnn : 3-phase training (main Exp6 contribution)
  5. MLP             : standard linear + MLP (no algebraic structure)

lambda_K normalization (Part E fix):
  lambda_K_eff = LAMBDA_K_BASE * N_GEN_REF / n_gen   with LAMBDA_K_BASE = 2.0, N_GEN_REF = 8
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
#  Path setup: killing_utils lives in the same directory
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from killing_utils import (
    killing_tension,
    closure_tension,
    analyze_generators,
    print_algebra_report,
    _su_basis,
    params_to_su_generators,
)

# ---------------------------------------------------------------------------
#  Global config
# ---------------------------------------------------------------------------
DEVICE          = torch.device("cuda" if torch.cuda.is_available() else "cpu")

LAMBDA_K_BASE   = 2.0      # Killing pressure base (validated in Exp4)
N_GEN_REF       = 8        # Reference n_gen for normalization (su(3))

LR_MAIN         = 1e-3     # Default learning rate
LR_ANNEAL       = 1e-4     # Phase 2 annealing learning rate
BATCH_SIZE      = 256
N_TRAIN         = 6_000
N_VAL           = 1_500
SEED            = 42

# Phase epochs
PH1_EPOCHS      = 600
PH2_EPOCHS      = 200
PH3_EPOCHS      = 200
SINGLE_EPOCHS   = 1_000    # Baselines (AlgNet-NoOrtho, AlgNet-Ortho1, MLP)

# Phase lambda_ortho values
PH1_LO          = 0.1
PH2_LO          = 2.0
PH3_LO          = 0.5
SINGLE_LO       = 1.0    # AlgNet-Ortho1

D_VALUES        = [4, 5]

OUT_DIR = os.path.join(HERE, "exp6_outputs")
os.makedirs(OUT_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
#  Helper: normalized lambda_K
# ---------------------------------------------------------------------------

def lambda_k_eff(d: int) -> float:
    """
    Normalized Killing pressure: constant effective pressure per generator.
    From Part E: lambda_K_eff = LAMBDA_K_BASE * N_GEN_REF / n_gen
    """
    n_gen = d * d - 1
    return LAMBDA_K_BASE * N_GEN_REF / n_gen


# ============================================================
#  T_ortho -- orthogonality loss (the neural T_closure analog)
# ============================================================

def ortho_tension(theta: torch.Tensor) -> torch.Tensor:
    """
    T_ortho = || theta @ theta.T - I_{n_gen} ||_F^2

    theta : (n_gen_active, n_gen_full)  real, learnable

    Penalizes deviation of theta from O(n_gen_active): the manifold of
    orthonormal row bases in R^{n_gen_full}.  When theta in O(n_gen), each
    row selects an orthogonal direction in the GML basis space, giving maximal
    independence of the learned generators.

    Physical analog: gamma * T_closure in SS-PEG confines commutators to
    span(generators). T_ortho confines theta to the generator manifold itself.
    """
    gram   = theta @ theta.T                                       # (n, n)
    I_n    = torch.eye(gram.shape[0], dtype=theta.dtype, device=theta.device)
    return ((gram - I_n) ** 2).sum()


# ============================================================
#  Dataset
# ============================================================

def generate_su_d_cubic_dataset(
    d: int, n: int, seed: int = SEED,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    x in R^{d^2-1} (GML coordinates), y = sign(Re(Tr(M(x)^3))).
    Exact same generation as Exp5 for comparability.
    """
    n_gen = d * d - 1
    rng   = np.random.default_rng(seed)
    X     = rng.standard_normal((n, n_gen)).astype(np.float64)

    basis = _su_basis(d, "cpu").numpy()
    M_c   = np.einsum("na,aij->nij", X, basis)
    M3_c  = M_c @ M_c @ M_c
    tr3   = np.trace(M3_c, axis1=-2, axis2=-1).real

    y = (tr3 > 0).astype(np.int64)
    pos = y.sum()
    print(f"  su({d}) dataset  n={n}  balance: {n-pos} neg / {pos} pos  (frac={pos/n:.3f})")
    return (
        torch.tensor(X, dtype=torch.float64, device=DEVICE),
        torch.tensor(y, dtype=torch.long,    device=DEVICE),
    )


def make_loaders(d: int, seed: int = SEED):
    X_tr, y_tr = generate_su_d_cubic_dataset(d, N_TRAIN, seed=seed)
    X_va, y_va = generate_su_d_cubic_dataset(d, N_VAL,   seed=seed + 1)
    tr = DataLoader(TensorDataset(X_tr, y_tr), batch_size=BATCH_SIZE, shuffle=True)
    va = DataLoader(TensorDataset(X_va, y_va), batch_size=512,         shuffle=False)
    return tr, va


# ============================================================
#  Models
# ============================================================

def make_head(n_gen: int) -> nn.Sequential:
    hidden = min(4 * n_gen, 128)
    return nn.Sequential(
        nn.ReLU(),
        nn.Linear(n_gen, hidden, dtype=torch.float64),
        nn.ReLU(),
        nn.Linear(hidden, 2,      dtype=torch.float64),
    )


class AlgebraicProjection(nn.Module):
    """
    Algebraic first layer: G_a = sum_k theta_{ak} T_k, where T_k is the GML basis.
    Regularization is computed EXTERNALLY (in training loop) for phase-dependent control.
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
        """(n_gen, d, d) complex generators -- gradient-tracked."""
        return params_to_su_generators(self._theta(), self.d)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        dev   = str(self._theta().device)
        basis = _su_basis(self.d, dev)
        M_x   = torch.einsum("...b,bij->...ij", x.to(torch.complex128), basis)
        G     = self.get_generators()
        f     = torch.einsum("aij,...ji->...a", G, M_x).real
        return f + self.bias


class AlgebraicNet(nn.Module):
    """AlgebraicProjection + classification head."""

    def __init__(self, d: int, fixed: bool = False):
        super().__init__()
        self.alg_layer = AlgebraicProjection(d=d, fixed=fixed)
        self.head      = make_head(d * d - 1)

    def forward(self, x):
        return self.head(self.alg_layer(x))

    def get_generators(self):
        with torch.no_grad():
            return self.alg_layer.get_generators()


class StandardMLP(nn.Module):
    """Baseline: linear layer (no algebraic structure) + same head."""

    def __init__(self, d: int):
        super().__init__()
        n = d * d - 1
        self.layer1 = nn.Linear(n, n, dtype=torch.float64)
        self.head   = make_head(n)

    def forward(self, x):
        return self.head(self.layer1(x))


# ============================================================
#  Training loop (single phase, external regularization)
# ============================================================

def _evaluate(model: nn.Module, val_loader) -> float:
    model.eval()
    correct, total = 0, 0
    with torch.no_grad():
        for xb, yb in val_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            correct += (model(xb).argmax(1) == yb).sum().item()
            total   += xb.size(0)
    return correct / total


def _algebraic_diagnostics(model: nn.Module) -> tuple[float, float, float]:
    """Returns (T_K, T_ortho, OrthoErr) with no_grad."""
    if not isinstance(model, AlgebraicNet) or model.alg_layer.fixed:
        return float("nan"), float("nan"), float("nan")
    with torch.no_grad():
        G     = model.alg_layer.get_generators()
        theta = model.alg_layer.theta
        t_k   = killing_tension(G).item()
        t_o   = ortho_tension(theta).item()
        th    = theta.detach().cpu().double().numpy()
        o_err = float(np.linalg.norm(th @ th.T - np.eye(th.shape[0]), "fro"))
    return t_k, t_o, o_err


def train_phase(
    model       : nn.Module,
    train_loader,
    val_loader,
    n_epochs    : int,
    lambda_k    : float,
    lambda_ortho: float,
    lr          : float   = LR_MAIN,
    freeze_head : bool    = False,
    phase_id    : int     = 1,
    warmup      : int     = 30,
    quiet       : bool    = False,
    report_every: int     = 100,
) -> dict:
    """
    Train one phase with fixed hyperparameters.
    Regularization is applied externally to allow per-phase lambda changes.
    """
    criterion = nn.CrossEntropyLoss()

    # Freeze/unfreeze head parameters
    if hasattr(model, "head"):
        for p in model.head.parameters():
            p.requires_grad = not freeze_head

    optimizer = optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=lr,
    )

    hist = {"val_acc": [], "T_K": [], "T_ortho": [], "ortho_err": []}
    best_val = 0.0

    if not quiet:
        fh_str = " [head FROZEN]" if freeze_head else ""
        print(f"    Phase {phase_id}{fh_str}: {n_epochs}ep  "
              f"lr={lr}  lK={lambda_k:.4f}  lO={lambda_ortho:.4f}")

    for epoch in range(n_epochs):
        # Linear warmup within this phase
        warmup_scale = min(1.0, (epoch + 1) / max(warmup, 1))
        lk_ep = lambda_k    * warmup_scale
        lo_ep = lambda_ortho * warmup_scale

        model.train()
        for xb, yb in train_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            optimizer.zero_grad()
            logits    = model(xb)
            loss      = criterion(logits, yb)

            if isinstance(model, AlgebraicNet) and not model.alg_layer.fixed:
                G     = model.alg_layer.get_generators()
                theta = model.alg_layer.theta
                if lk_ep > 0:
                    loss = loss + lk_ep * killing_tension(G)
                if lo_ep > 0:
                    loss = loss + lo_ep * ortho_tension(theta)

            loss.backward()
            optimizer.step()

        val_acc  = _evaluate(model, val_loader)
        best_val = max(best_val, val_acc)
        t_k, t_o, o_err = _algebraic_diagnostics(model)

        hist["val_acc"].append(val_acc)
        hist["T_K"].append(t_k)
        hist["T_ortho"].append(t_o)
        hist["ortho_err"].append(o_err)

        if not quiet and (epoch + 1) % report_every == 0:
            print(f"      Ph{phase_id} ep{epoch+1:4d}/{n_epochs}  "
                  f"val={val_acc:.4f}  T_K={t_k:.4f}  T_o={t_o:.6f}  OErr={o_err:.2f}")

    hist["best_val"] = best_val
    return hist


# ============================================================
#  Multi-phase training (Exp6 main contribution)
# ============================================================

def train_multiphase(
    model       : nn.Module,
    train_loader,
    val_loader,
    d           : int,
    name        : str  = "AlgNet-OrthoAnn",
    quiet       : bool = False,
) -> dict:
    """
    3-phase training schedule -- analog to F4 Stage 1 (gamma=50) + Stage 2 (gamma=200).

    Phase 1: task + weak ortho confinement (learn task geometry, stay near manifold)
    Phase 2: annealing -- freeze head, strong Killing + ortho pressure
             -> theta retracts to O(n_gen), generators crystallize
    Phase 3: fine-tune -- unfreeze head, moderate ortho
             -> task performance recovers with crystallized generators
    """
    lk = lambda_k_eff(d)

    print(f"\n  [{name}]  3-phase training  lambda_K_eff={lk:.4f}")
    print(f"    Ph1: {PH1_EPOCHS}ep  lK={lk:.4f}  lO={PH1_LO}")
    print(f"    Ph2: {PH2_EPOCHS}ep  lK={5*lk:.4f}  lO={PH2_LO}  [freeze head, lr=1e-4]")
    print(f"    Ph3: {PH3_EPOCHS}ep  lK={lk:.4f}  lO={PH3_LO}  [unfreeze]")

    h1 = train_phase(
        model, train_loader, val_loader,
        n_epochs=PH1_EPOCHS, lambda_k=lk, lambda_ortho=PH1_LO,
        lr=LR_MAIN, freeze_head=False, phase_id=1,
        quiet=quiet, report_every=150,
    )
    h2 = train_phase(
        model, train_loader, val_loader,
        n_epochs=PH2_EPOCHS, lambda_k=5 * lk, lambda_ortho=PH2_LO,
        lr=LR_ANNEAL, freeze_head=True, phase_id=2,
        quiet=quiet, report_every=50,
    )
    # Ensure head is unfrozen for Phase 3
    if hasattr(model, "head"):
        for p in model.head.parameters():
            p.requires_grad = True

    h3 = train_phase(
        model, train_loader, val_loader,
        n_epochs=PH3_EPOCHS, lambda_k=lk, lambda_ortho=PH3_LO,
        lr=LR_MAIN, freeze_head=False, phase_id=3,
        quiet=quiet, report_every=100,
    )

    total_ep = PH1_EPOCHS + PH2_EPOCHS + PH3_EPOCHS
    combined = {
        "val_acc"   : h1["val_acc"]    + h2["val_acc"]    + h3["val_acc"],
        "T_K"       : h1["T_K"]        + h2["T_K"]        + h3["T_K"],
        "T_ortho"   : h1["T_ortho"]    + h2["T_ortho"]    + h3["T_ortho"],
        "ortho_err" : h1["ortho_err"]  + h2["ortho_err"]  + h3["ortho_err"],
        "best_val"  : max(h1["best_val"], h2["best_val"], h3["best_val"]),
        "phase_boundaries": [PH1_EPOCHS, PH1_EPOCHS + PH2_EPOCHS],
        "total_epochs": total_ep,
    }
    return combined


# ============================================================
#  Single-phase training wrapper (for baselines)
# ============================================================

def train_single(
    model       : nn.Module,
    train_loader,
    val_loader,
    d           : int,
    name        : str   = "Model",
    lambda_ortho: float = 0.0,
    n_epochs    : int   = SINGLE_EPOCHS,
    quiet       : bool  = False,
) -> dict:
    """Single-phase training: T_K (normalized) + optional T_ortho, no annealing."""
    lk = lambda_k_eff(d) if isinstance(model, AlgebraicNet) and not model.alg_layer.fixed else 0.0
    print(f"\n  [{name}]  {n_epochs} epochs  lK={lk:.4f}  lO={lambda_ortho:.4f}")

    h = train_phase(
        model, train_loader, val_loader,
        n_epochs=n_epochs, lambda_k=lk, lambda_ortho=lambda_ortho,
        lr=LR_MAIN, freeze_head=False, phase_id=1,
        quiet=quiet, report_every=200,
    )
    h["phase_boundaries"] = []
    h["total_epochs"]     = n_epochs
    return h


# ============================================================
#  Run experiment for one d value
# ============================================================

def run_d(d: int) -> dict:
    n_gen = d * d - 1
    lk    = lambda_k_eff(d)
    print(f"\n{'='*70}")
    print(f"  Exp6  d={d}  su({d})  n_gen={n_gen}  lambda_K_eff={lk:.4f}")
    print(f"{'='*70}")

    tr, va = make_loaders(d)

    results: dict[str, dict] = {}

    # ---------------------------------------------------------
    # 1. GML-Fixed (oracle)
    # ---------------------------------------------------------
    model = AlgebraicNet(d=d, fixed=True).to(DEVICE)
    results["GML-Fixed"] = train_single(
        model, tr, va, d=d, name="GML-Fixed",
        lambda_ortho=0.0, n_epochs=SINGLE_EPOCHS,
    )

    # ---------------------------------------------------------
    # 2. MLP baseline
    # ---------------------------------------------------------
    model = StandardMLP(d=d).to(DEVICE)
    results["MLP"] = train_single(
        model, tr, va, d=d, name="MLP",
        lambda_ortho=0.0, n_epochs=SINGLE_EPOCHS,
    )

    # ---------------------------------------------------------
    # 3. AlgNet-NoOrtho  (Part E baseline: T_K only, normalized lambda)
    # ---------------------------------------------------------
    model = AlgebraicNet(d=d, fixed=False).to(DEVICE)
    results["AlgNet-NoOrtho"] = train_single(
        model, tr, va, d=d, name="AlgNet-NoOrtho",
        lambda_ortho=0.0, n_epochs=SINGLE_EPOCHS,
    )

    # ---------------------------------------------------------
    # 4. AlgNet-Ortho1  (single phase, T_K + T_ortho=1.0)
    # ---------------------------------------------------------
    model = AlgebraicNet(d=d, fixed=False).to(DEVICE)
    results["AlgNet-Ortho1"] = train_single(
        model, tr, va, d=d, name="AlgNet-Ortho1",
        lambda_ortho=SINGLE_LO, n_epochs=SINGLE_EPOCHS,
    )

    # ---------------------------------------------------------
    # 5. AlgNet-OrthoAnn  (3-phase annealing -- main Exp6 model)
    # ---------------------------------------------------------
    model = AlgebraicNet(d=d, fixed=False).to(DEVICE)
    results["AlgNet-OrthoAnn"] = train_multiphase(
        model, tr, va, d=d, name="AlgNet-OrthoAnn",
    )

    # ---------------------------------------------------------
    #  Summary table
    # ---------------------------------------------------------
    print(f"\n  EXP6 SUMMARY -- d={d}  (su({d}), n_gen={n_gen}, lambda_K_eff={lk:.4f})")
    header = (f"  {'Model':<20}  {'Best acc':>10}  {'Final acc':>10}  "
              f"{'T_K final':>10}  {'T_ortho':>9}  {'OrthoErr':>9}  {'Eps':>5}")
    print(header)
    print("  " + "-" * (len(header) - 2))

    mlp_best  = results["MLP"]["best_val"]
    gml_best  = results["GML-Fixed"]["best_val"]

    for name, h in results.items():
        best   = h["best_val"]
        final  = h["val_acc"][-1]
        t_k    = h["T_K"][-1]
        t_o    = h["T_ortho"][-1]
        o_err  = h["ortho_err"][-1]
        ep     = h["total_epochs"]
        delta  = best - mlp_best

        tk_s  = f"{t_k:10.4f}"  if not (isinstance(t_k,  float) and np.isnan(t_k))  else f"{'--':>10}"
        to_s  = f"{t_o:9.4f}"   if not (isinstance(t_o,  float) and np.isnan(t_o))  else f"{'--':>9}"
        oe_s  = f"{o_err:9.2f}" if not (isinstance(o_err, float) and np.isnan(o_err)) else f"{'--':>9}"
        print(
            f"  {name:<20}  {best:>10.4f}  {final:>10.4f}  "
            f"{tk_s}  {to_s}  {oe_s}  {ep:>5}"
            f"  [Delta_MLP={delta:+.4f}]"
        )

    print(f"\n  GML-Fixed (oracle): {gml_best:.4f}")
    print(f"  MLP baseline:       {mlp_best:.4f}")

    # ---------------------------------------------------------
    #  Algebraic report for trained models
    # ---------------------------------------------------------
    print(f"\n  Exp6 algebraic fingerprints (d={d}):")
    for name in ("AlgNet-NoOrtho", "AlgNet-Ortho1", "AlgNet-OrthoAnn"):
        # Retrieve model -- need to reconstruct since we discarded it
        # (We track diagnostics in history, full report here uses final snapshot)
        pass  # diagnostics already in hist["T_K"][-1] etc.

    return results


# ============================================================
#  Plotting
# ============================================================

def plot_results(all_results: dict[int, dict[str, dict]], out_path: str) -> None:
    """
    2 rows (one per d) x 3 cols (val_acc, T_K, OrthoErr).
    Phase boundaries shown as vertical dashed lines.
    """
    n_d = len(all_results)
    fig, axes = plt.subplots(n_d, 3, figsize=(16, 5 * n_d))
    if n_d == 1:
        axes = axes[None, :]   # 2D indexing

    color_map = {
        "GML-Fixed"       : ("tab:green",  "--"),
        "MLP"             : ("tab:blue",   "-"),
        "AlgNet-NoOrtho"  : ("tab:orange", "-"),
        "AlgNet-Ortho1"   : ("tab:red",    "-"),
        "AlgNet-OrthoAnn" : ("tab:purple", "-"),
    }

    for row_idx, (d, results) in enumerate(all_results.items()):
        # Determine total epochs and phase boundaries from OrthoAnn
        ann_h = results.get("AlgNet-OrthoAnn", {})
        boundaries = ann_h.get("phase_boundaries", [])

        ax_acc, ax_tk, ax_oe = axes[row_idx]
        ax_acc.set_title(f"d={d} — Validation Accuracy")
        ax_tk .set_title(f"d={d} — Killing Tension T_K")
        ax_oe .set_title(f"d={d} — Ortho Error ||theta*theta.T - I||")

        for name, h in results.items():
            c, ls = color_map.get(name, ("black", "-"))
            eps = list(range(1, len(h["val_acc"]) + 1))

            ax_acc.plot(eps, h["val_acc"],  color=c, ls=ls, lw=1.8, label=name)

            tk_vals = [max(v, 1e-6) if not np.isnan(v) else np.nan for v in h["T_K"]]
            if not all(np.isnan(v) for v in tk_vals):
                valid = [(e, v) for e, v in zip(eps, tk_vals) if not np.isnan(v)]
                if valid:
                    ev, vv = zip(*valid)
                    ax_tk.semilogy(ev, vv, color=c, ls=ls, lw=1.8, label=name)

            oe_vals = [v if not np.isnan(v) else np.nan for v in h["ortho_err"]]
            if not all(np.isnan(v) for v in oe_vals):
                valid = [(e, v) for e, v in zip(eps, oe_vals) if not np.isnan(v)]
                if valid:
                    ev, vv = zip(*valid)
                    ax_oe.plot(ev, vv, color=c, ls=ls, lw=1.8, label=name)

        # Phase boundary lines for OrthoAnn
        for ax in (ax_acc, ax_tk, ax_oe):
            for bnd in boundaries:
                ax.axvline(bnd, color="gray", ls=":", lw=1.2, alpha=0.6)

        for ax in (ax_acc, ax_tk, ax_oe):
            ax.set_xlabel("Epoch")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

        ax_acc.set_ylabel("Accuracy")
        ax_tk .set_ylabel("T_K (log scale)")
        ax_oe .set_ylabel("OrthoErr (Frobenius)")

    fig.suptitle(
        "Exp6 — T_ortho Annealing  (G2/F4 closure-confinement analog)\n"
        "Gray dotted lines: phase boundaries (Ph1|Ph2=600, Ph2|Ph3=800)",
        fontsize=11,
    )
    plt.tight_layout()
    plt.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close()
    print(f"\n  Figure saved -> {out_path}")


# ============================================================
#  Main
# ============================================================

def main():
    print(f"\n{'#'*70}")
    print(f"  Experiment 6 -- T_ortho Annealing (physical confinement analog)")
    print(f"  Device: {DEVICE}")
    if DEVICE.type == "cuda":
        print(f"  GPU: {torch.cuda.get_device_name(0)}")
    print(f"  d values: {D_VALUES}")
    print(f"  Total epochs per model: single={SINGLE_EPOCHS},"
          f"  multiphase={PH1_EPOCHS+PH2_EPOCHS+PH3_EPOCHS}")
    print(f"  lambda_K_base={LAMBDA_K_BASE}  N_gen_ref={N_GEN_REF}")
    print(f"  Phase lambdas: Ph1_ortho={PH1_LO}, Ph2_ortho={PH2_LO}, Ph3_ortho={PH3_LO}")
    print(f"{'#'*70}")

    all_results: dict[int, dict] = {}
    t0_total = time.perf_counter()

    for d in D_VALUES:
        t0 = time.perf_counter()
        all_results[d] = run_d(d)
        elapsed = time.perf_counter() - t0
        print(f"\n  d={d} done in {elapsed:.1f}s")

    # Final cross-d summary
    print(f"\n\n{'='*70}")
    print(f"  EXP6 CROSS-d SUMMARY")
    print(f"{'='*70}")
    print(f"  {'d':>3}  {'n_gen':>5}  {'lK_eff':>8}  {'GML-Fixed':>10}  "
          f"{'NoOrtho':>9}  {'Ortho1':>8}  {'OrthoAnn':>10}  {'MLP':>8}")
    print(f"  " + "-" * 80)
    for d, results in all_results.items():
        n_gen = d * d - 1
        lk    = lambda_k_eff(d)
        vals  = {k: results[k]["best_val"] for k in results}
        print(
            f"  {d:>3}  {n_gen:>5}  {lk:>8.4f}  "
            f"{vals.get('GML-Fixed', float('nan')):>10.4f}  "
            f"{vals.get('AlgNet-NoOrtho', float('nan')):>9.4f}  "
            f"{vals.get('AlgNet-Ortho1', float('nan')):>8.4f}  "
            f"{vals.get('AlgNet-OrthoAnn', float('nan')):>10.4f}  "
            f"{vals.get('MLP', float('nan')):>8.4f}"
        )
    print(f"\n  Deltas (OrthoAnn - MLP):")
    for d, results in all_results.items():
        ann   = results.get("AlgNet-OrthoAnn", {}).get("best_val", float("nan"))
        mlp   = results.get("MLP",             {}).get("best_val", float("nan"))
        noort = results.get("AlgNet-NoOrtho",  {}).get("best_val", float("nan"))
        print(
            f"    d={d}:  OrthoAnn - MLP = {ann - mlp:+.4f}  "
            f"|  NoOrtho - MLP = {noort - mlp:+.4f}  "
            f"|  OrthoAnn vs NoOrtho = {ann - noort:+.4f}"
        )

    total_time = time.perf_counter() - t0_total
    print(f"\n  Total time: {total_time:.1f}s")

    # Save plot
    plot_results(
        all_results,
        out_path=os.path.join(OUT_DIR, "exp6_ortho_anneal.png"),
    )

    # Save numeric results
    save_dict = {}
    for d, results in all_results.items():
        for name, h in results.items():
            key = f"d{d}_{name.replace('-', '_')}"
            save_dict[key + "_val_acc"]   = np.array(h["val_acc"])
            save_dict[key + "_T_K"]       = np.array(h["T_K"])
            save_dict[key + "_T_ortho"]   = np.array(h["T_ortho"])
            save_dict[key + "_ortho_err"] = np.array(h["ortho_err"])
            save_dict[key + "_best"]      = np.array([h["best_val"]])

    npz_path = os.path.join(OUT_DIR, "exp6_results.npz")
    np.savez(npz_path, **save_dict)
    print(f"  Numeric results saved -> {npz_path}")


if __name__ == "__main__":
    main()
