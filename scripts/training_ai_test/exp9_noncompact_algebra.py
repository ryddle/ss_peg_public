"""
Exp9: Non-Compact Algebra — Killing Regularization on sl(3,ℝ)
==============================================================

Question:
  Does algebraic regularization (T_K + T_ortho) generalize to non-compact
  Lie algebras?  Specifically, does the COMPACT Killing tension (which
  targets K ∝ −I, i.e., negative-definite) harm performance on a task whose
  underlying symmetry is the NON-COMPACT algebra sl(3,ℝ)?

Design:
  Task      : sign(Tr(M³)),  M = Σ_a x_a B_a ∈ sl(3,ℝ)
              sl(3,ℝ) has n_gen = 8 real 3×3 traceless generators —
              the same count as su(3) but with indefinite Killing form.
  n_gen     : 8,  n_dim = 3  (real matrices)
  n_train   : 6 000  (fixed)
  n_val     : 1 500
  N_reps    : 3 seeds
  Epochs    : 1 500  (same as Exp8 n=6000)
  LR        : 1e-3  (Adam, grad-clip max_norm=1.0)

Models:
  1. MLP           — no algebraic structure (baseline)
  2. AlgNet-Compact — T_K_compact (targets K ∝ −I, WRONG prior) + T_ortho
  3. AlgNet-NC      — T_K_noncompact (targets K_ref of sl(3,ℝ), CORRECT prior) + T_ortho
  4. AlgNet-NoTK    — T_ortho only, no T_K (pure manifold confinement)
  5. Oracle         — θ = I fixed (exact sl(3,ℝ) generators), head trains only

Diagnosis:
  For each model, report the eigenvalue signature of K(θ) at convergence.
  For sl(3,ℝ), K_ref has 2 positive + 4 zero + 2 negative eigenvalues
  (due to the mixed-signature Killing form of a non-compact algebra).

Hypothesis (G13):
  a) AlgNet-NC ≥ AlgNet-Compact  (correct structural prior helps)
  b) AlgNet-Compact may be BELOW MLP if compact T_K actively fights sl(3,ℝ) structure
  c) AlgNet-NoTK provides a lower bound for the benefit of Killing geometry
  d) Oracle sets the upper bound

Outputs:
  exp9_outputs/exp9_results.png
  exp9_outputs/exp9_results.npz
"""

# ---------------------------------------------------------------------------
#  Paths and imports
# ---------------------------------------------------------------------------
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import time
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from killing_utils import killing_metric, killing_tension

# ---------------------------------------------------------------------------
#  Device
# ---------------------------------------------------------------------------
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {DEVICE}")

# ---------------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------------
N_GEN       = 8     # sl(3,ℝ): n = 3^2 - 1
N_DIM       = 3     # generator dimension

LAMBDA_K_BASE = 2.0
N_GEN_REF     = 8
LK_EFF        = LAMBDA_K_BASE * N_GEN_REF / N_GEN   # 2.0  (N_GEN=8 → LK_EFF=2.0)
LK_2X         = 2.0 * LK_EFF                         # 4.0
# Use same conservative scaling as exp7/8 (8-generator formula)
LK_FINAL      = 2.0 * 8 / N_GEN * 2.0               # ≈ 4.0 → use exp7 scale
LAMBDA_K      = LK_EFF                                # 2.0  (single factor)
LAMBDA_ORTHO  = 1.0

BATCH_SIZE  = 256
LR_MAIN     = 1e-3
N_EPOCHS    = 1500
REPORT_EVERY= 300
N_TRAIN     = 6_000
N_VAL       = 1_500
BASE_SEEDS  = [42, 43, 44]
N_REPS      = 3

OUT_DIR = os.path.join(HERE, "exp9_outputs")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
#  sl(3,ℝ) basis
# ---------------------------------------------------------------------------

def _sl3r_basis(device_str: str = "cpu") -> torch.Tensor:
    """
    8 real traceless 3×3 generators for sl(3,ℝ) in the Chevalley-Serre basis.

    Ordering:
      0: H1 = diag(1, -1,  0)              — Cartan
      1: H2 = diag(0,  1, -1)              — Cartan
      2: E01 [entry (0,1) = 1]             — raising
      3: E02 [entry (0,2) = 1]             — raising
      4: E12 [entry (1,2) = 1]             — raising
      5: F01 [entry (1,0) = 1]             — lowering
      6: F02 [entry (2,0) = 1]             — lowering
      7: F12 [entry (2,1) = 1]             — lowering

    These satisfy: Tr(B_a B_b) = non-trivial indefinite metric.
    The Killing form K_ab = Tr(ad_a ad_b) is indefinite for sl(3,ℝ).
    """
    device = torch.device(device_str)
    d      = torch.float64
    Z      = lambda: torch.zeros(3, 3, dtype=d, device=device)

    gens = []
    # Cartan
    H1 = Z(); H1[0,0] =  1.0; H1[1,1] = -1.0
    H2 = Z(); H2[1,1] =  1.0; H2[2,2] = -1.0
    gens += [H1, H2]
    # Raising (E_ij = e_ij, i < j)
    for (i, j) in [(0,1), (0,2), (1,2)]:
        E = Z(); E[i, j] = 1.0; gens.append(E)
    # Lowering (F_ij = e_ji, i < j)
    for (i, j) in [(0,1), (0,2), (1,2)]:
        F = Z(); F[j, i] = 1.0; gens.append(F)

    return torch.stack(gens)   # (8, 3, 3) real float64


# Pre-compute reference Killing form for sl(3,ℝ) (used by nc_killing_tension)
def _sl3r_killing_ref(device_str: str = "cpu") -> torch.Tensor:
    """
    Normalised reference Killing form K_ref for sl(3,ℝ).
    K_ref[a,b] = Tr(ad_a ∘ ad_b)  normalised so mean |diag| = 1.
    """
    basis = _sl3r_basis(device_str)
    K = killing_metric(basis)                                  # (8, 8) real
    scale = K.diagonal().abs().mean().clamp(min=1e-8)
    return K / scale                                           # (8, 8) real, normalised


# ---------------------------------------------------------------------------
#  Tension helpers
# ---------------------------------------------------------------------------

def ortho_tension(theta: torch.Tensor) -> torch.Tensor:
    """T_ortho = ||theta @ theta.T − I||_F^2  (manifold confinement)."""
    gram = theta @ theta.T
    I_n  = torch.eye(gram.shape[0], dtype=theta.dtype, device=theta.device)
    return ((gram - I_n) ** 2).sum()


def nc_killing_tension(G: torch.Tensor, K_ref_norm: torch.Tensor) -> torch.Tensor:
    """
    Non-compact Killing tension.
    Penalises deviation of K(G) from the reference indefinite K_ref of sl(3,ℝ).

    G         : (n_gen, 3, 3) real learnable generators
    K_ref_norm: (n_gen, n_gen) normalised reference Killing matrix (precomputed)
    Returns   : scalar ≥ 0
    """
    K = killing_metric(G)                                      # (8, 8) real
    scale = K.abs().mean().clamp(min=1e-8)
    K_n   = K / scale
    return ((K_n - K_ref_norm) ** 2).sum()


# ---------------------------------------------------------------------------
#  Dataset
# ---------------------------------------------------------------------------

def generate_sl3r_cubic_dataset(n: int, seed: int):
    """
    x ∈ R^8  (sl(3,ℝ) coordinates w.r.t. the Chevalley-Serre basis).
    M = Σ_a x_a B_a ∈ sl(3,ℝ)   (real 3×3 traceless)
    y = sign(Tr(M³))

    By antisymmetry of the cubic polynomial under x → −x (and N(0,1) symmetry),
    P(y=1) = 0.5 exactly — class-balanced without any correction.
    """
    rng   = np.random.default_rng(seed)
    X     = rng.standard_normal((n, N_GEN)).astype(np.float64)
    basis = _sl3r_basis("cpu").numpy()                          # (8, 3, 3)
    M     = np.einsum("na,aij->nij", X, basis)                  # (n, 3, 3)
    M2    = np.einsum("nij,njk->nik",  M,  M)
    M3    = np.einsum("nij,njk->nik",  M2, M)
    tr3   = np.trace(M3, axis1=-2, axis2=-1)                    # (n,)
    y     = (tr3 > 0).astype(np.int64)
    pos   = y.sum()
    print(f"  sl(3,ℝ) cubic dataset  n={n}  "
          f"balance: {n-pos} neg / {pos} pos  (frac={pos/n:.3f})")
    return (
        torch.tensor(X, dtype=torch.float64, device=DEVICE),
        torch.tensor(y, dtype=torch.long,    device=DEVICE),
    )


# ---------------------------------------------------------------------------
#  Models
# ---------------------------------------------------------------------------

def _make_head(n_in: int) -> nn.Sequential:
    hidden = min(4 * n_in, 128)
    return nn.Sequential(
        nn.ReLU(),
        nn.Linear(n_in,    hidden, dtype=torch.float64),
        nn.ReLU(),
        nn.Linear(hidden, 2,      dtype=torch.float64),
    )


class AlgebraicNetReal(nn.Module):
    """
    Algebraic projection layer for real Lie algebras.

    Generators G_a = Σ_k theta_{ak} B_k  (real linear combination of basis).
    Feature f_a = Tr(G_a @ M(x))  for each sample x.

    If fixed=True, theta is set to I and frozen (Oracle mode).
    """

    def __init__(self, fixed: bool = False):
        super().__init__()
        init  = torch.eye(N_GEN, dtype=torch.float64)
        if fixed:
            self.register_buffer("theta_buf", init)
        else:
            noise = torch.randn_like(init) * 0.05
            self.theta = nn.Parameter(init + noise)
        self.bias   = nn.Parameter(torch.zeros(N_GEN, dtype=torch.float64))
        self.head   = _make_head(N_GEN)
        self.fixed  = fixed

    def _theta(self) -> torch.Tensor:
        return self.theta_buf if self.fixed else self.theta

    def _get_generators(self) -> torch.Tensor:
        dev   = str(self._theta().device)
        basis = _sl3r_basis(dev)
        return torch.einsum("ab,bij->aij", self._theta(), basis)   # (8, 3, 3)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        dev   = str(self._theta().device)
        basis = _sl3r_basis(dev)
        # M(x) = Σ_a x_a B_a
        M_x   = torch.einsum("...a,aij->...ij", x, basis)           # (..., 3, 3)
        G     = self._get_generators()                               # (8, 3, 3)
        # f_a = Tr(G_a @ M(x)) = Σ_{ij} G_a[i,j] M_x[j,i]
        f     = torch.einsum("aij,...ji->...a", G, M_x)              # (..., 8)
        return self.head(f + self.bias)


class StandardMLP(nn.Module):
    """Baseline MLP: linear layer + head."""

    def __init__(self):
        super().__init__()
        self.layer1 = nn.Linear(N_GEN, N_GEN, dtype=torch.float64)
        self.head   = _make_head(N_GEN)

    def forward(self, x):
        return self.head(self.layer1(x))


# ---------------------------------------------------------------------------
#  K-signature diagnostic
# ---------------------------------------------------------------------------

def _k_signature(model: nn.Module) -> str:
    """
    Returns the eigenvalue signature of the Killing form as a string '+p -m ~z'.
    Only meaningful for AlgebraicNetReal (non-Oracle).
    """
    if not isinstance(model, AlgebraicNetReal) or model.fixed:
        return "n/a"
    with torch.no_grad():
        G = model._get_generators()
        K = killing_metric(G).cpu().double().numpy()
        eigs = np.linalg.eigvalsh(K)
        pos  = int((eigs >  0.1 * abs(eigs).max()).sum())
        neg  = int((eigs < -0.1 * abs(eigs).max()).sum())
        zer  = N_GEN - pos - neg
    return f"+{pos} -{neg} ~{zer}"


# ---------------------------------------------------------------------------
#  Evaluation
# ---------------------------------------------------------------------------

def _evaluate(model: nn.Module, loader) -> float:
    model.eval()
    correct, total = 0, 0
    with torch.no_grad():
        for xb, yb in loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            correct += (model(xb).argmax(1) == yb).sum().item()
            total   += xb.size(0)
    return correct / total


# ---------------------------------------------------------------------------
#  Training
# ---------------------------------------------------------------------------

def train_one(
    model         : nn.Module,
    train_loader,
    val_loader,
    name          : str,
    lambda_k_comp : float,      # compact T_K coefficient
    lambda_k_nc   : float,      # non-compact T_K coefficient
    lambda_ortho  : float,
    K_ref_norm    : torch.Tensor,
    n_epochs      : int   = N_EPOCHS,
    report_every  : int   = REPORT_EVERY,
) -> dict:
    """
    Single-phase training with one of three regularization modes:
      compact:     lambda_k_comp > 0, lambda_k_nc = 0
      noncompact:  lambda_k_nc   > 0, lambda_k_comp = 0
      notk:        both = 0, lambda_ortho > 0
      oracle/mlp:  all = 0
    """
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=LR_MAIN,
    )
    WARMUP = 30
    best_val = 0.0
    hist = {"val_acc": [], "T_K_comp": [], "T_K_nc": []}

    print(f"\n  >> {name} ({n_epochs}ep  lKc={lambda_k_comp:.3f}"
          f"  lKnc={lambda_k_nc:.3f}  lo={lambda_ortho:.3f})")

    for epoch in range(n_epochs):
        scale  = min(1.0, (epoch + 1) / WARMUP)
        lkc_ep = lambda_k_comp * scale
        lknc_ep= lambda_k_nc   * scale
        lo_ep  = lambda_ortho  * scale

        model.train()
        for xb, yb in train_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            optimizer.zero_grad()
            loss = criterion(model(xb), yb)

            if isinstance(model, AlgebraicNetReal) and not model.fixed:
                G = model._get_generators()
                if lkc_ep > 0:
                    loss = loss + lkc_ep * killing_tension(G)
                if lknc_ep > 0:
                    loss = loss + lknc_ep * nc_killing_tension(G, K_ref_norm)
                if lo_ep > 0:
                    loss = loss + lo_ep  * ortho_tension(model.theta)

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

        val_acc  = _evaluate(model, val_loader)
        best_val = max(best_val, val_acc)

        # Compute diagnostic tensions
        if isinstance(model, AlgebraicNetReal) and not model.fixed:
            with torch.no_grad():
                G   = model._get_generators()
                tkc = killing_tension(G).item()
                tknc= nc_killing_tension(G, K_ref_norm).item()
        else:
            tkc = tknc = float("nan")

        hist["val_acc"].append(val_acc)
        hist["T_K_comp"].append(tkc)
        hist["T_K_nc"].append(tknc)

        if (epoch + 1) % report_every == 0:
            print(f"    [{name} ep{epoch+1:4d}]  val={val_acc:.4f}  best={best_val:.4f}")

    hist["best_val"] = best_val
    hist["k_sig"]    = _k_signature(model)
    return hist


# ---------------------------------------------------------------------------
#  Experiment runner
# ---------------------------------------------------------------------------

def run_exp9():
    """Outer loop: seed × 5 models."""
    X_val, y_val = generate_sl3r_cubic_dataset(N_VAL, seed=999)
    val_loader   = DataLoader(TensorDataset(X_val, y_val),
                              batch_size=512, shuffle=False)

    # Precompute reference K for non-compact tension (constant)
    K_ref_norm = _sl3r_killing_ref(str(DEVICE)).to(DEVICE)

    # Print K_ref signature for reference
    eigs_ref = np.linalg.eigvalsh(K_ref_norm.cpu().numpy())
    pos = int((eigs_ref >  0.1 * abs(eigs_ref).max()).sum())
    neg = int((eigs_ref < -0.1 * abs(eigs_ref).max()).sum())
    zer = N_GEN - pos - neg
    print(f"\n  sl(3,ℝ) K_ref eigenvalue signature: +{pos} -{neg} ~{zer}")
    print(f"  K_ref eigenvalues: {np.round(eigs_ref, 3)}")

    results = {name: [] for name in
               ["MLP", "AlgNet-Compact", "AlgNet-NC", "AlgNet-NoTK", "Oracle"]}

    for seed_idx, seed in enumerate(BASE_SEEDS):
        print(f"\n{'='*70}")
        print(f"  SEED {seed}  ({seed_idx+1}/{N_REPS})")
        print(f"{'='*70}")

        X_tr, y_tr = generate_sl3r_cubic_dataset(N_TRAIN, seed=seed)
        tr_loader  = DataLoader(TensorDataset(X_tr, y_tr),
                                batch_size=BATCH_SIZE, shuffle=True)

        # 1. MLP
        model = StandardMLP().to(DEVICE)
        h = train_one(model, tr_loader, val_loader, name=f"MLP-s{seed}",
                      lambda_k_comp=0.0, lambda_k_nc=0.0, lambda_ortho=0.0,
                      K_ref_norm=K_ref_norm)
        results["MLP"].append(h)
        print(f"    MLP           best={h['best_val']:.4f}  K-sig={h['k_sig']}")

        # 2. AlgNet-Compact (compact T_K — WRONG prior for sl(3,ℝ))
        model = AlgebraicNetReal(fixed=False).to(DEVICE)
        h = train_one(model, tr_loader, val_loader, name=f"AlgNet-Compact-s{seed}",
                      lambda_k_comp=LAMBDA_K, lambda_k_nc=0.0, lambda_ortho=LAMBDA_ORTHO,
                      K_ref_norm=K_ref_norm)
        results["AlgNet-Compact"].append(h)
        print(f"    AlgNet-Compact best={h['best_val']:.4f}  K-sig={h['k_sig']}")

        # 3. AlgNet-NC (non-compact T_K — CORRECT prior)
        model = AlgebraicNetReal(fixed=False).to(DEVICE)
        h = train_one(model, tr_loader, val_loader, name=f"AlgNet-NC-s{seed}",
                      lambda_k_comp=0.0, lambda_k_nc=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
                      K_ref_norm=K_ref_norm)
        results["AlgNet-NC"].append(h)
        print(f"    AlgNet-NC      best={h['best_val']:.4f}  K-sig={h['k_sig']}")

        # 4. AlgNet-NoTK (T_ortho only, no T_K)
        model = AlgebraicNetReal(fixed=False).to(DEVICE)
        h = train_one(model, tr_loader, val_loader, name=f"AlgNet-NoTK-s{seed}",
                      lambda_k_comp=0.0, lambda_k_nc=0.0, lambda_ortho=LAMBDA_ORTHO,
                      K_ref_norm=K_ref_norm)
        results["AlgNet-NoTK"].append(h)
        print(f"    AlgNet-NoTK    best={h['best_val']:.4f}  K-sig={h['k_sig']}")

        # 5. Oracle (fixed exact sl(3,ℝ) generators)
        model = AlgebraicNetReal(fixed=True).to(DEVICE)
        h = train_one(model, tr_loader, val_loader, name=f"Oracle-s{seed}",
                      lambda_k_comp=0.0, lambda_k_nc=0.0, lambda_ortho=0.0,
                      K_ref_norm=K_ref_norm)
        results["Oracle"].append(h)
        print(f"    Oracle         best={h['best_val']:.4f}  K-sig={h['k_sig']}")

    return results


# ---------------------------------------------------------------------------
#  Summary printing
# ---------------------------------------------------------------------------

def print_summary(results: dict) -> None:
    from scipy import stats as sp_stats

    oracle_mean = float(np.mean([h["best_val"] for h in results["Oracle"]]))
    print(f"\n{'='*80}")
    print(f"  EXP9 SUMMARY  — sl(3,ℝ) Non-Compact Algebra  (n_train={N_TRAIN})")
    print(f"{'='*80}")
    print(f"  Oracle mean best: {oracle_mean:.4f}")
    print()
    header = (
        f"  {'Model':<22}  {'Mean':>6}  {'Std':>6}  "
        f"{'K-sig (final)':>14}  {'Δ_Oracle':>9}"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))

    for name, hist_list in results.items():
        vals  = [h["best_val"] for h in hist_list]
        k_sigs= [h["k_sig"]   for h in hist_list]
        m, s  = float(np.mean(vals)), float(np.std(vals))
        delta = m - oracle_mean
        k_str = k_sigs[-1]  # last seed's signature
        print(f"  {name:<22}  {m:.4f}  {s:.4f}  {k_str:>14}  {delta:+.4f}")

    # Pairwise t-tests vs MLP
    print()
    print(f"  {'Comparison':<35}  {'Δ':>7}  {'p-val':>7}")
    print("  " + "-" * 55)
    mlp_vals = [h["best_val"] for h in results["MLP"]]
    for cmp_name in ["AlgNet-Compact", "AlgNet-NC", "AlgNet-NoTK"]:
        cmp_vals = [h["best_val"] for h in results[cmp_name]]
        delta    = float(np.mean(cmp_vals)) - float(np.mean(mlp_vals))
        try:
            _, pv = sp_stats.ttest_ind(cmp_vals, mlp_vals, equal_var=False)
        except Exception:
            pv = float("nan")
        print(f"  {cmp_name} vs MLP          {delta:+.4f}  {pv:.4f}")


# ---------------------------------------------------------------------------
#  Plotting
# ---------------------------------------------------------------------------

def plot_results(results: dict, path: str) -> None:
    model_names = ["MLP", "AlgNet-Compact", "AlgNet-NC", "AlgNet-NoTK", "Oracle"]
    colors      = ["royalblue", "tomato", "forestgreen", "darkorange", "purple"]

    means = [float(np.mean([h["best_val"] for h in results[n]])) for n in model_names]
    stds  = [float(np.std( [h["best_val"] for h in results[n]])) for n in model_names]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Bar chart
    x = np.arange(len(model_names))
    bars = ax1.bar(x, means, yerr=stds, capsize=5, color=colors, alpha=0.8, width=0.6)
    ax1.axhline(0.5, color="black", linestyle="--", linewidth=0.8, label="chance")
    ax1.set_xticks(x)
    ax1.set_xticklabels(model_names, rotation=15, ha="right")
    ax1.set_ylabel("Best val accuracy  (mean ± std, 3 seeds)")
    ax1.set_title(f"Exp9 — sl(3,ℝ) cubic task  (n_train={N_TRAIN})")
    ax1.legend()
    ax1.grid(True, axis="y", alpha=0.3)
    ax1.set_ylim(0.48, min(1.0, max(means) + 0.05))

    # Training curves (mean over 3 seeds)
    ep = np.arange(1, N_EPOCHS + 1)
    for name, color in zip(model_names, colors):
        curves = np.array([[v for v in h["val_acc"]] for h in results[name]])
        m      = curves.mean(axis=0)
        s      = curves.std(axis=0)
        ax2.plot(ep, m, color=color, label=name, linewidth=1.5)
        ax2.fill_between(ep, m - s, m + s, color=color, alpha=0.15)
    ax2.set_xlabel("Epoch")
    ax2.set_ylabel("Val accuracy")
    ax2.set_title(f"Exp9 — Training curves  (mean ± std over {N_REPS} seeds)")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    print(f"\n  Figure saved: {path}")


# ---------------------------------------------------------------------------
#  Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    t0 = time.time()

    print(__doc__)
    print(f"LAMBDA_K={LAMBDA_K:.4f}  LAMBDA_ORTHO={LAMBDA_ORTHO:.4f}")
    print(f"n_train={N_TRAIN}  n_val={N_VAL}  n_epochs={N_EPOCHS}  seeds={BASE_SEEDS}")

    results = run_exp9()
    print_summary(results)

    fig_path = os.path.join(OUT_DIR, "exp9_results.png")
    plot_results(results, fig_path)

    # Save raw results
    npz_path = os.path.join(OUT_DIR, "exp9_results.npz")
    save_dict = {}
    for name, hist_list in results.items():
        key = name.replace("-", "_").replace("(", "").replace(")", "")
        save_dict[f"{key}_bests"] = np.array([h["best_val"] for h in hist_list])
        for i, h in enumerate(hist_list):
            save_dict[f"{key}_s{i}_val_acc"]   = np.array(h["val_acc"])
            save_dict[f"{key}_s{i}_TK_comp"]   = np.array(h["T_K_comp"])
            save_dict[f"{key}_s{i}_TK_nc"]     = np.array(h["T_K_nc"])
    np.savez(npz_path, **save_dict)
    print(f"  Data saved:   {npz_path}")

    elapsed = time.time() - t0
    print(f"\n  Wall time: {elapsed:.1f}s  ({elapsed/60:.1f}min)")
    print("  Exp9 complete.")
