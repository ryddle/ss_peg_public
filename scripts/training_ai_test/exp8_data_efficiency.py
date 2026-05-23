"""
Exp8: Data-Efficiency Sweep
============================

Question:
  Does the algebraic inductive bias (Killing crystallization + T_ortho)
  provide a larger accuracy advantage over MLP when training data is scarce?

Design:
  Task      : sign(Re(Tr(M^3))),  M ∈ su(4)  [full algebra, all 15 generators]
  d         : 4  (n_gen = 15)
  n_train   : {100, 300, 1 000, 3 000, 6 000}
  n_val     : 1 500  (fixed — large enough for stable accuracy estimates)
  Models    : MLP  vs  AlgNet (T_K_2x + T_ortho, single-phase, no gate)
  N_reps    : 3 seeds per (model, n_train)
  Epochs    : adaptive  {100: 4000, 300: 3000, 1000: 2000, 3000: 1500, 6000: 1500}
  LR        : 1e-3  (Adam)

Hypothesis (G12):
  Δ(AlgNet − MLP) increases as n_train decreases.
  The structural prior is most valuable in the data-scarce regime.
  Possible crossover: AlgNet may overtake MLP at small n.

Outputs:
  exp8_outputs/exp8_results.png  — two-panel figure
  exp8_outputs/exp8_results.npz  — raw arrays for post-analysis
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

from killing_utils import killing_tension, _su_basis, params_to_su_generators

# ---------------------------------------------------------------------------
#  Device
# ---------------------------------------------------------------------------
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {DEVICE}")

# ---------------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------------
D         = 4
N_GEN     = D * D - 1          # 15

# Killing / ortho lambda — same normalisation as exp7 (proven numerically stable).
# N_GEN_REF=8 keeps per-generator penalty scale identical to the exp7 oracle run.
LAMBDA_K_BASE = 2.0
N_GEN_REF     = 8
LK_EFF        = LAMBDA_K_BASE * N_GEN_REF / N_GEN   # ≈ 1.0667  (exp7 formula)
LK_2X         = 2.0 * LK_EFF                         # ≈ 2.1333  (2× schedule)
LAMBDA_ORTHO  = 1.0

BATCH_SIZE  = 256
LR_MAIN     = 1e-3

N_VAL       = 1500
N_REPS      = 3
N_TRAIN_VALS = [100, 300, 1_000, 3_000, 6_000]
EPOCH_MAP    = {100: 4000, 300: 3000, 1_000: 2000, 3_000: 1500, 6_000: 1500}
REPORT_EVERY = 500

BASE_SEEDS = [42, 43, 44]

OUT_DIR = os.path.join(HERE, "exp8_outputs")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
#  Dataset
# ---------------------------------------------------------------------------

def generate_full_su4_cubic_dataset(n: int, seed: int):
    """
    x ∈ R^{15}  (su(4) GML coordinates, all generators).
    y = sign(Re(Tr(M^3)))  where M = Σ_a x_a T_a uses ALL 15 su(4) generators.
    """
    rng   = np.random.default_rng(seed)
    X     = rng.standard_normal((n, N_GEN)).astype(np.float64)
    basis = _su_basis(D, "cpu").numpy()                     # (15, 4, 4) complex
    M     = np.einsum("na,aij->nij", X, basis)              # (n, 4, 4) complex
    M3    = M @ M @ M
    tr3   = np.trace(M3, axis1=-2, axis2=-1).real
    y     = (tr3 > 0).astype(np.int64)
    pos   = y.sum()
    print(f"  su({D}) full dataset  n={n}  "
          f"balance: {n-pos} neg / {pos} pos  (frac={pos/n:.3f})")
    return (
        torch.tensor(X, dtype=torch.float64, device=DEVICE),
        torch.tensor(y, dtype=torch.long,    device=DEVICE),
    )


# ---------------------------------------------------------------------------
#  Tension helpers
# ---------------------------------------------------------------------------

def ortho_tension(theta: torch.Tensor) -> torch.Tensor:
    """T_ortho = ||theta @ theta.T - I||_F^2"""
    gram = theta @ theta.T
    I_n  = torch.eye(gram.shape[0], dtype=theta.dtype, device=theta.device)
    return ((gram - I_n) ** 2).sum()


# ---------------------------------------------------------------------------
#  Models
# ---------------------------------------------------------------------------

def _make_head(n_in: int) -> nn.Sequential:
    hidden = min(4 * n_in, 128)
    return nn.Sequential(
        nn.ReLU(),
        nn.Linear(n_in,   hidden, dtype=torch.float64),
        nn.ReLU(),
        nn.Linear(hidden, 2,     dtype=torch.float64),
    )


class AlgebraicNet(nn.Module):
    """
    AlgebraicProjection + classification head, no gate.
    Forward: f_a = Re(Tr(G_a M(x)))  for each a ∈ {0…n_gen-1}.
    T_K and T_ortho are applied externally during training.
    """

    def __init__(self):
        super().__init__()
        init       = torch.eye(N_GEN, dtype=torch.float64)
        noise      = torch.randn_like(init) * 0.05
        self.theta = nn.Parameter(init + noise)
        self.bias  = nn.Parameter(torch.zeros(N_GEN, dtype=torch.float64))
        self.head  = _make_head(N_GEN)

    def _get_generators(self) -> torch.Tensor:
        return params_to_su_generators(self.theta, D)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        dev  = str(self.theta.device)
        basis = _su_basis(D, dev)
        M_x  = torch.einsum("...b,bij->...ij", x.to(torch.complex128), basis)
        G    = self._get_generators()
        f    = torch.einsum("aij,...ji->...a", G, M_x).real
        return self.head(f + self.bias)


class StandardMLP(nn.Module):
    """Baseline MLP: linear layer + head (no algebraic structure)."""

    def __init__(self):
        super().__init__()
        self.layer1 = nn.Linear(N_GEN, N_GEN, dtype=torch.float64)
        self.head   = _make_head(N_GEN)

    def forward(self, x):
        return self.head(self.layer1(x))


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
    model        : nn.Module,
    train_loader,
    val_loader,
    name         : str,
    lambda_k     : float,
    lambda_ortho : float,
    n_epochs     : int,
    report_every : int = REPORT_EVERY,
) -> float:
    """
    Single-phase training. Returns best validation accuracy.
    """
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=LR_MAIN,
    )
    WARMUP   = 30
    best_val = 0.0

    for epoch in range(n_epochs):
        scale  = min(1.0, (epoch + 1) / WARMUP)
        lk_ep  = lambda_k     * scale
        lo_ep  = lambda_ortho * scale

        model.train()
        for xb, yb in train_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            optimizer.zero_grad()
            loss = criterion(model(xb), yb)

            if isinstance(model, AlgebraicNet):
                G = model._get_generators()
                if lk_ep > 0:
                    loss = loss + lk_ep * killing_tension(G)
                if lo_ep > 0:
                    loss = loss + lo_ep * ortho_tension(model.theta)

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

        val_acc  = _evaluate(model, val_loader)
        best_val = max(best_val, val_acc)

        if (epoch + 1) % report_every == 0:
            print(f"      [{name} ep{epoch+1:4d}]  val={val_acc:.4f}  best={best_val:.4f}")

    return best_val


# ---------------------------------------------------------------------------
#  Experiment runner
# ---------------------------------------------------------------------------

def run_exp8() -> dict:
    """
    Outer loop: n_train × model × seed.
    Returns nested dict results[n_train][model_name] = list-of-best_val.
    """
    # Fixed validation set (large — does not change across n_train sweep)
    X_val, y_val = generate_full_su4_cubic_dataset(N_VAL, seed=999)
    val_loader   = DataLoader(
        TensorDataset(X_val, y_val), batch_size=512, shuffle=False
    )

    results = {n: {"MLP": [], "AlgNet": []} for n in N_TRAIN_VALS}

    total_cells = len(N_TRAIN_VALS) * 2 * N_REPS
    cell = 0

    for n_train in N_TRAIN_VALS:
        n_epochs = EPOCH_MAP[n_train]
        print(f"\n{'='*70}")
        print(f"  n_train={n_train}  epochs={n_epochs}")
        print(f"{'='*70}")

        for rep, seed in enumerate(BASE_SEEDS):
            cell += 2
            print(f"\n  --- rep {rep+1}/{N_REPS}  seed={seed}  "
                  f"[cell {cell-1}&{cell}/{total_cells}] ---")

            # Training data for this (n_train, seed)
            X_tr, y_tr = generate_full_su4_cubic_dataset(n_train, seed=seed)
            tr_loader  = DataLoader(
                TensorDataset(X_tr, y_tr),
                batch_size=min(BATCH_SIZE, n_train),
                shuffle=True,
            )

            # MLP
            mlp   = StandardMLP().to(DEVICE)
            best  = train_one(
                mlp, tr_loader, val_loader, name=f"MLP-n{n_train}-s{seed}",
                lambda_k=0.0, lambda_ortho=0.0, n_epochs=n_epochs,
            )
            results[n_train]["MLP"].append(best)
            print(f"    MLP      best={best:.4f}")

            # AlgNet (T_K_2x + T_ortho)
            algnet = AlgebraicNet().to(DEVICE)
            best   = train_one(
                algnet, tr_loader, val_loader, name=f"AlgNet-n{n_train}-s{seed}",
                lambda_k=LK_2X, lambda_ortho=LAMBDA_ORTHO, n_epochs=n_epochs,
            )
            results[n_train]["AlgNet"].append(best)
            print(f"    AlgNet   best={best:.4f}")

    return results


# ---------------------------------------------------------------------------
#  Summary printing
# ---------------------------------------------------------------------------

def print_summary(results: dict) -> None:
    from scipy import stats as sp_stats

    print(f"\n{'='*80}")
    print(f"  EXP8 SUMMARY  — Data-Efficiency Sweep  (d={D}, full su({D}) cubic task)")
    print(f"{'='*80}")
    header = (
        f"  {'n_train':>8}  {'MLP mean±std':>14}  {'AlgNet mean±std':>16}"
        f"  {'Δ':>7}  {'p-val':>7}  {'winner':>8}"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))

    for n_train in N_TRAIN_VALS:
        mlp_vals    = results[n_train]["MLP"]
        alg_vals    = results[n_train]["AlgNet"]
        mlp_m, mlp_s   = float(np.mean(mlp_vals)), float(np.std(mlp_vals))
        alg_m, alg_s   = float(np.mean(alg_vals)), float(np.std(alg_vals))
        delta          = alg_m - mlp_m
        try:
            _, pv = sp_stats.ttest_ind(alg_vals, mlp_vals, equal_var=False)
        except Exception:
            pv = float("nan")
        winner = "AlgNet" if delta > 0 else "MLP   "
        print(f"  {n_train:>8}  {mlp_m:.4f}±{mlp_s:.4f}  "
              f"{alg_m:.4f}±{alg_s:.4f}  {delta:+.4f}  {pv:.4f}  {winner}")


# ---------------------------------------------------------------------------
#  Plotting
# ---------------------------------------------------------------------------

def plot_results(results: dict, path: str) -> None:
    n_trains = N_TRAIN_VALS

    mlp_means,    mlp_stds    = [], []
    alg_means,    alg_stds    = [], []
    deltas,       delta_errs  = [], []

    for n in n_trains:
        mlp_v = np.array(results[n]["MLP"])
        alg_v = np.array(results[n]["AlgNet"])
        mlp_m, mlp_s = mlp_v.mean(), mlp_v.std()
        alg_m, alg_s = alg_v.mean(), alg_v.std()
        mlp_means.append(mlp_m);  mlp_stds.append(mlp_s)
        alg_means.append(alg_m);  alg_stds.append(alg_s)
        deltas.append(alg_m - mlp_m)
        delta_errs.append(np.sqrt(mlp_s**2 + alg_s**2))  # propagated std

    x      = np.array(n_trains, dtype=float)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # --- Panel 1: accuracy vs n_train ---
    ax1.errorbar(x, mlp_means, yerr=mlp_stds,
                 fmt="o-", color="royalblue",  label="MLP",    capsize=4)
    ax1.errorbar(x, alg_means, yerr=alg_stds,
                 fmt="s-", color="tomato",     label="AlgNet", capsize=4)
    ax1.set_xscale("log")
    ax1.set_xlabel("n_train")
    ax1.set_ylabel("Best val accuracy  (mean ± std, 3 seeds)")
    ax1.set_title(f"Exp8 — Accuracy vs n_train  (d={D})")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(n_trains)
    ax1.set_xticklabels([str(n) for n in n_trains])

    # --- Panel 2: Δ = AlgNet − MLP ---
    colors = ["tomato" if d > 0 else "royalblue" for d in deltas]
    ax2.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax2.errorbar(x, deltas, yerr=delta_errs,
                 fmt="D-", color="purple", capsize=4, label="Δ (AlgNet − MLP)")
    for xi, yi, c in zip(x, deltas, colors):
        ax2.scatter([xi], [yi], color=c, zorder=5)
    ax2.set_xscale("log")
    ax2.set_xlabel("n_train")
    ax2.set_ylabel("Δ accuracy  (AlgNet − MLP)")
    ax2.set_title(f"Exp8 — Algebraic Advantage vs n_train  (d={D})")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(n_trains)
    ax2.set_xticklabels([str(n) for n in n_trains])

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    print(f"\n  Figure saved: {path}")


# ---------------------------------------------------------------------------
#  Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    t0 = time.time()

    print(__doc__)
    print(f"LK_2X={LK_2X:.4f}  LAMBDA_ORTHO={LAMBDA_ORTHO:.4f}")
    print(f"n_train sweep: {N_TRAIN_VALS}")
    print(f"Seeds: {BASE_SEEDS}  (N_reps={N_REPS})")
    print(f"n_val: {N_VAL}")

    results = run_exp8()
    print_summary(results)

    fig_path = os.path.join(OUT_DIR, "exp8_results.png")
    plot_results(results, fig_path)

    # Save raw results
    npz_path = os.path.join(OUT_DIR, "exp8_results.npz")
    save_dict = {}
    for n in N_TRAIN_VALS:
        save_dict[f"mlp_{n}"]    = np.array(results[n]["MLP"])
        save_dict[f"algnet_{n}"] = np.array(results[n]["AlgNet"])
    np.savez(npz_path, **save_dict)
    print(f"  Data saved:   {npz_path}")

    elapsed = time.time() - t0
    print(f"\n  Wall time: {elapsed:.1f}s  ({elapsed/60:.1f}min)")
    print("  Exp8 complete.")
