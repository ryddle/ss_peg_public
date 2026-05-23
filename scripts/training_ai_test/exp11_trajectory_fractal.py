"""
Exp11: Optimizer Trajectory Fractal Dimension (Grassberger-Procaccia D₂)
=========================================================================

Question:
  Does the AlgNet optimizer trajectory in θ-space live on a fractal attractor?
  If the correlation dimension D₂ ∉ ℤ, the algebraic constraint creates genuine
  fractal geometry in parameter space — connecting to Julia sets in ℂ² and
  percolation critical clusters.

Design:
  Task      : sign(Re(Tr(M^3))),  M = Σ_a x_a T_a ∈ su(3)
              N_GEN = 8, N_DIM = 3, complex anti-Hermitian generators
  n_train   : 6 000  (same as Exp9/10)
  n_val     : 1 500
  N_REPS    : 5 seeds
  Epochs    : 3 000
  LR        : 1e-3, grad-clip max_norm=1.0
  λ_K       : 2.0,  λ_ortho : 1.0
  RECORD_EVERY : 5 epochs  → 600 trajectory points per run

Trajectory recording:
  AlgNet-1L : θ_t ∈ ℝ^64  (8×8 rotation parameter matrix)
  MLP-1L    : w_t ∈ ℝ^64  (first 8×8 weight matrix, for fair comparison)
  Both      : full trajectory + final 50% (convergence phase)

Grassberger-Procaccia correlation dimension D₂:
  1. Apply PCA to trajectory: project to K_PCA=16 principal components
  2. Compute pairwise distances || p_i − p_j ||₂  (N×N matrix)
  3. Correlation sum: C(r, N) = (2 / N(N-1)) Σ_{i<j} θ(r − dist(i,j))
  4. Scaling region: r ∈ [r_min, r_max] where C(r) is power-law
  5. D₂ = slope of log C(r) vs log r  (linear regression in scaling region)
  6. Bootstrap CI on D₂

Comparisons:
  - AlgNet D₂ (full trajectory) vs MLP D₂ (full trajectory)
  - AlgNet D₂ (convergence phase) vs MLP D₂ (convergence phase)
  - Hypothesis: AlgNet trajectory confined to algebraic manifold → lower D₂

Secondary analysis:
  - PCA variance explained vs component index
  - 2D PCA projection of trajectories (visual)
  - Compare trajectory length ||Δθ||₂ accumulated

Outputs:
  exp11_outputs/exp11_results.png     — correlation function C(r), D₂ vs model
  exp11_outputs/exp11_trajectory.png  — PCA projections + length
  exp11_outputs/exp11_results.npz     — full trajectory arrays and D₂ estimates
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

from killing_utils import killing_metric, killing_tension, _su_basis

# ---------------------------------------------------------------------------
#  Device
# ---------------------------------------------------------------------------
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {DEVICE}")

# ---------------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------------
N_DIM           = 3
N_GEN           = N_DIM ** 2 - 1      # 8 for su(3)

LAMBDA_K        = 2.0
LAMBDA_ORTHO    = 1.0
WARMUP_EP       = 30

BATCH_SIZE      = 256
LR_MAIN         = 1e-3
N_EPOCHS        = 3_000
RECORD_EVERY    = 5                    # record θ every this many epochs
REPORT_EVERY    = 500

N_TRAIN         = 6_000
N_VAL           = 1_500
N_REPS          = 5
BASE_SEEDS      = [42, 43, 44, 45, 46]

# Grassberger-Procaccia settings
K_PCA           = 16        # PCA components before GP
N_R_POINTS      = 60        # number of r values for C(r)
SCALING_FRAC    = (0.1, 0.6)  # fraction of r range used for slope fit
N_BOOTSTRAP     = 200       # bootstrap samples for D₂ CI

OUT_DIR = os.path.join(HERE, "exp11_outputs")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
#  Dataset: su(3) cubic trace
# ---------------------------------------------------------------------------

def generate_su3_cubic_dataset(n: int, seed: int):
    rng   = np.random.default_rng(seed)
    X     = rng.standard_normal((n, N_GEN)).astype(np.float64)
    basis = _su_basis(N_DIM, "cpu").numpy()
    M     = np.einsum("na,aij->nij", X, basis)
    M3    = M @ M @ M
    tr3   = np.trace(M3, axis1=-2, axis2=-1).real
    y     = (tr3 > 0).astype(np.int64)
    return (
        torch.tensor(X, dtype=torch.float64, device=DEVICE),
        torch.tensor(y, dtype=torch.long,    device=DEVICE),
    )


# ---------------------------------------------------------------------------
#  Models
# ---------------------------------------------------------------------------

class AlgebraicLayer(nn.Module):
    def __init__(self):
        super().__init__()
        init_t = torch.eye(N_GEN, dtype=torch.float64)
        self.theta = nn.Parameter(init_t + torch.randn_like(init_t) * 0.05)
        self.bias  = nn.Parameter(torch.zeros(N_GEN, dtype=torch.float64))

    def _get_generators(self, basis):
        return torch.einsum("ab,bij->aij", self.theta.to(basis.dtype), basis)

    def forward(self, x, basis):
        G   = self._get_generators(basis)
        M_x = torch.einsum("...a,aij->...ij", x.to(basis.dtype), basis)
        f   = torch.einsum("aij,...ji->...a", G, M_x)
        return f.real.to(torch.float64) + self.bias

    def t_ortho(self):
        gram = self.theta @ self.theta.T
        return ((gram - torch.eye(N_GEN, dtype=torch.float64, device=self.theta.device))**2).sum()

    def t_k(self, basis):
        return killing_tension(self._get_generators(basis))


class AlgNet1L(nn.Module):
    def __init__(self):
        super().__init__()
        self.layer = AlgebraicLayer()
        self.relu  = nn.ReLU()
        self.head  = nn.Linear(N_GEN, 2, dtype=torch.float64)

    def forward(self, x, basis):
        return self.head(self.relu(self.layer(x, basis)))

    def total_reg(self, basis, scale):
        return (LAMBDA_K * self.layer.t_k(basis) + LAMBDA_ORTHO * self.layer.t_ortho()) * scale

    def get_trajectory_point(self) -> np.ndarray:
        """Return θ as a flat CPU array for trajectory recording."""
        return self.layer.theta.detach().cpu().numpy().ravel().copy()


class MLP1L(nn.Module):
    def __init__(self):
        super().__init__()
        self.fc1  = nn.Linear(N_GEN, N_GEN, dtype=torch.float64)
        self.relu = nn.ReLU()
        self.head = nn.Linear(N_GEN, 2, dtype=torch.float64)

    def forward(self, x, basis=None):
        return self.head(self.relu(self.fc1(x)))

    def get_trajectory_point(self) -> np.ndarray:
        """Return W₁ (first layer weights) as a flat CPU array."""
        return self.fc1.weight.detach().cpu().numpy().ravel().copy()


# ---------------------------------------------------------------------------
#  Evaluation
# ---------------------------------------------------------------------------

def evaluate(model, loader, basis) -> float:
    model.eval()
    correct, total = 0, 0
    with torch.no_grad():
        for xb, yb in loader:
            out = model(xb, basis)
            correct += (out.argmax(1) == yb).sum().item()
            total   += xb.size(0)
    return correct / total


# ---------------------------------------------------------------------------
#  Training with trajectory recording
# ---------------------------------------------------------------------------

def train_and_record(model, train_loader, val_loader, basis,
                     name: str, is_algnet: bool) -> dict:
    """
    Train model for N_EPOCHS, recording trajectory every RECORD_EVERY epochs.
    Returns trajectory array (n_points, n_params) and accuracy history.
    """
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()), lr=LR_MAIN
    )

    trajectory = []
    val_accs   = []
    best_val   = 0.0

    for ep in range(1, N_EPOCHS + 1):
        warmup_scale = min(1.0, ep / WARMUP_EP) if is_algnet else 0.0

        model.train()
        for xb, yb in train_loader:
            optimizer.zero_grad()
            logits = model(xb, basis)
            loss   = criterion(logits, yb)
            if is_algnet:
                loss = loss + model.total_reg(basis, warmup_scale)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

        if ep % RECORD_EVERY == 0:
            pt = model.get_trajectory_point()
            trajectory.append(pt)

        val_acc = evaluate(model, val_loader, basis)
        val_accs.append(val_acc)
        if val_acc > best_val:
            best_val = val_acc

        if ep % REPORT_EVERY == 0 or ep == N_EPOCHS:
            print(f"    ep {ep:4d}  val={val_acc:.4f}  best={best_val:.4f}")

    trajectory = np.array(trajectory)   # (n_points, n_params)
    print(f"  {name}: trajectory shape = {trajectory.shape}, "
          f"total displacement = {_total_displacement(trajectory):.4f}")

    return {
        "name":       name,
        "trajectory": trajectory,       # (n_points, n_params=64)
        "val_accs":   np.array(val_accs),
        "best_val":   best_val,
    }


def _total_displacement(traj: np.ndarray) -> float:
    """Sum of step-by-step L2 distances along trajectory."""
    if len(traj) < 2:
        return 0.0
    steps = np.diff(traj, axis=0)
    return float(np.linalg.norm(steps, axis=1).sum())


# ---------------------------------------------------------------------------
#  Grassberger-Procaccia algorithm
# ---------------------------------------------------------------------------

def pca_project(traj: np.ndarray, k: int) -> np.ndarray:
    """Zero-mean PCA projection of trajectory onto top-k principal components."""
    mu    = traj.mean(axis=0)
    X     = traj - mu
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    return X @ Vt[:k].T    # (n_points, k)


def correlation_sum(pts: np.ndarray, r_values: np.ndarray) -> np.ndarray:
    """
    C(r) = (2 / N(N-1)) Σ_{i<j} θ(r - ||p_i − p_j||₂)
    Computed efficiently with batched distance matrix.
    pts      : (N, k)
    r_values : (M,)
    Returns  : (M,) array of C(r)
    """
    N = len(pts)
    # Pairwise distances — full matrix (upper triangle only)
    diff   = pts[:, None, :] - pts[None, :, :]          # (N, N, k)
    dists  = np.sqrt((diff**2).sum(axis=-1))             # (N, N)
    # Extract upper triangle
    idx    = np.triu_indices(N, k=1)
    d_triu = dists[idx]                                  # (N*(N-1)/2,)
    n_pairs = len(d_triu)

    C = np.array([np.sum(d_triu < r) / n_pairs for r in r_values])
    return C


def estimate_d2(traj: np.ndarray, label: str, phase: str = "full") -> dict:
    """
    Estimate correlation dimension D₂ via Grassberger-Procaccia on 'traj'.
    traj  : (T, P) trajectory
    Returns dict with D2, CI low/high, r values, C values, PCA variance.
    """
    T, P = traj.shape
    k    = min(K_PCA, T - 1, P)

    # PCA projection
    mu    = traj.mean(axis=0)
    X     = traj - mu
    _, sv, Vt = np.linalg.svd(X, full_matrices=False)
    var_explained = (sv**2) / (sv**2).sum()
    pts   = X @ Vt[:k].T         # (T, k)

    # r range from 5th–95th percentile of pairwise distances (sample 500 pairs for speed)
    rng   = np.random.default_rng(0)
    idx_s = rng.choice(T, size=min(T, 200), replace=False)
    pts_s = pts[idx_s]
    diff_s  = pts_s[:, None, :] - pts_s[None, :, :]
    d_s     = np.sqrt((diff_s**2).sum(axis=-1))
    d_flat  = d_s[np.triu_indices(len(idx_s), k=1)]
    d_flat  = d_flat[d_flat > 1e-15]
    r_min   = np.percentile(d_flat, 5)
    r_max   = np.percentile(d_flat, 95)
    r_vals  = np.logspace(np.log10(r_min), np.log10(r_max), N_R_POINTS)

    # Full correlation sum
    C_vals  = correlation_sum(pts, r_vals)

    # Find scaling region: exclude C~0 and C~1
    valid   = (C_vals > 0.02) & (C_vals < 0.98)
    lo_frac, hi_frac = SCALING_FRAC
    n_v     = valid.sum()
    if n_v < 5:
        print(f"  WARNING: {label} [{phase}]: only {n_v} valid C(r) points — D₂ unreliable")
        lo_idx = max(0, int(len(r_vals)*lo_frac))
        hi_idx = max(lo_idx+4, int(len(r_vals)*hi_frac))
        valid_idx = np.arange(lo_idx, min(hi_idx+1, len(r_vals)))
    else:
        valid_idx = np.where(valid)[0]
        lo_idx = valid_idx[int(len(valid_idx)*lo_frac)]
        hi_idx = valid_idx[int(len(valid_idx)*hi_frac)]
        valid_idx = np.arange(lo_idx, hi_idx+1)

    log_r = np.log(r_vals[valid_idx])
    log_C = np.log(np.clip(C_vals[valid_idx], 1e-15, None))

    # Ordinary least squares slope
    A  = np.column_stack([log_r, np.ones_like(log_r)])
    if len(log_r) < 2:
        d2 = float("nan")
        d2_std = float("nan")
    else:
        coef, *_ = np.linalg.lstsq(A, log_C, rcond=None)
        d2 = coef[0]

        # Bootstrap CI
        d2_boots = []
        for _ in range(N_BOOTSTRAP):
            bi = rng.choice(len(log_r), size=len(log_r), replace=True)
            if len(np.unique(bi)) < 2:
                continue
            c_, *_ = np.linalg.lstsq(A[bi], log_C[bi], rcond=None)
            d2_boots.append(c_[0])
        d2_boots = np.array(d2_boots)
        d2_std = d2_boots.std() if len(d2_boots) > 5 else float("nan")

    print(f"  {label} [{phase}]: T={T}, k={k}, D₂={d2:.3f}±{d2_std:.3f}  "
          f"valid_r=[{r_vals[valid_idx[0]]:.3e}, {r_vals[valid_idx[-1]]:.3e}]"
          f"  var_expl(PC1-3)={var_explained[:3].sum():.3f}")

    return {
        "d2":           d2,
        "d2_std":       d2_std,
        "r_vals":       r_vals,
        "C_vals":       C_vals,
        "valid_idx":    valid_idx,
        "var_explained": var_explained,
        "label":        label,
        "phase":        phase,
        "T":            T,
    }


# ---------------------------------------------------------------------------
#  Run one seed
# ---------------------------------------------------------------------------

def run_one_seed(seed: int):
    print(f"\n{'='*60}")
    print(f"Seed {seed}")
    print(f"{'='*60}")

    X_tr, y_tr = generate_su3_cubic_dataset(N_TRAIN, seed=seed)
    X_vl, y_vl = generate_su3_cubic_dataset(N_VAL,   seed=seed + 100)
    basis       = _su_basis(N_DIM, DEVICE)

    tr_loader = DataLoader(TensorDataset(X_tr, y_tr), batch_size=BATCH_SIZE, shuffle=True)
    vl_loader = DataLoader(TensorDataset(X_vl, y_vl), batch_size=512)

    torch.manual_seed(seed * 1000 + 1)
    algnet = AlgNet1L().to(DEVICE)
    torch.manual_seed(seed * 1000 + 2)
    mlp    = MLP1L().to(DEVICE)

    print("  --- AlgNet ---")
    res_alg = train_and_record(algnet, tr_loader, vl_loader, basis,
                                "AlgNet-1L", is_algnet=True)
    print("  --- MLP ---")
    res_mlp = train_and_record(mlp, tr_loader, vl_loader, basis,
                                "MLP-1L", is_algnet=False)

    # D₂ computation: full trajectory and convergence phase (last 50%)
    n_half = len(res_alg["trajectory"]) // 2

    d2_alg_full = estimate_d2(res_alg["trajectory"],          "AlgNet", "full")
    d2_alg_conv = estimate_d2(res_alg["trajectory"][n_half:], "AlgNet", "convergence")
    d2_mlp_full = estimate_d2(res_mlp["trajectory"],          "MLP",    "full")
    d2_mlp_conv = estimate_d2(res_mlp["trajectory"][n_half:], "MLP",    "convergence")

    return {
        "seed":       seed,
        "alg":        res_alg,
        "mlp":        res_mlp,
        "d2_alg_full": d2_alg_full,
        "d2_alg_conv": d2_alg_conv,
        "d2_mlp_full": d2_mlp_full,
        "d2_mlp_conv": d2_mlp_conv,
    }


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def main():
    t_start = time.time()
    seed_results = []

    for seed in BASE_SEEDS:
        seed_results.append(run_one_seed(seed))

    # Aggregate D₂ across seeds
    d2_alg_full_l = [r["d2_alg_full"]["d2"] for r in seed_results]
    d2_alg_conv_l = [r["d2_alg_conv"]["d2"] for r in seed_results]
    d2_mlp_full_l = [r["d2_mlp_full"]["d2"] for r in seed_results]
    d2_mlp_conv_l = [r["d2_mlp_conv"]["d2"] for r in seed_results]

    def _fmt(lst):
        a = np.array([x for x in lst if not np.isnan(x)])
        if len(a) == 0:
            return "n/a"
        return f"{a.mean():.3f} ± {a.std(ddof=1):.3f}" if len(a) > 1 else f"{a[0]:.3f}"

    print("\n" + "="*60)
    print("  D₂ Summary (mean ± std across seeds)")
    print("  ----------------------------------------")
    print(f"  AlgNet  full       : {_fmt(d2_alg_full_l)}")
    print(f"  AlgNet  convergence: {_fmt(d2_alg_conv_l)}")
    print(f"  MLP     full       : {_fmt(d2_mlp_full_l)}")
    print(f"  MLP     convergence: {_fmt(d2_mlp_conv_l)}")

    # Save NPZ
    npz_path = os.path.join(OUT_DIR, "exp11_results.npz")
    np.savez(
        npz_path,
        d2_alg_full = np.array(d2_alg_full_l),
        d2_alg_conv = np.array(d2_alg_conv_l),
        d2_mlp_full = np.array(d2_mlp_full_l),
        d2_mlp_conv = np.array(d2_mlp_conv_l),
        # Store one representative trajectory per model (seed 0)
        traj_alg_0  = seed_results[0]["alg"]["trajectory"],
        traj_mlp_0  = seed_results[0]["mlp"]["trajectory"],
    )
    print(f"Saved {npz_path}")

    _plot_results(seed_results)
    _plot_trajectories(seed_results)

    elapsed = time.time() - t_start
    print(f"\nTotal runtime: {elapsed/60:.1f} min")


# ---------------------------------------------------------------------------
#  Plotting
# ---------------------------------------------------------------------------

def _plot_results(seed_results: list):
    """Plot C(r) curves and D₂ bar chart."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Left & Center: C(r) curves for seed 0 (full, conv)
    for ax, phase_key, title in [
        (axes[0], "d2_alg_full", "AlgNet  full trajectory"),
        (axes[1], "d2_mlp_full", "MLP  full trajectory"),
    ]:
        res0 = seed_results[0][phase_key]
        ax.loglog(res0["r_vals"], np.clip(res0["C_vals"], 1e-10, None),
                  "o-", markersize=3, color="steelblue" if "alg" in phase_key.lower() else "coral")
        # Scaling-region markers
        vi = res0["valid_idx"]
        ax.loglog(res0["r_vals"][vi], np.clip(res0["C_vals"][vi], 1e-10, None),
                  "x", markersize=6, color="black", label=f"D₂={res0['d2']:.2f}")
        ax.set_xlabel("r")
        ax.set_ylabel("C(r)")
        ax.set_title(title)
        ax.legend()
        ax.grid(True, which="both", alpha=0.3)

    # Right: D₂ bar chart (full trajectory)
    ax = axes[2]
    d2_alg = np.array([r["d2_alg_full"]["d2"] for r in seed_results])
    d2_mlp = np.array([r["d2_mlp_full"]["d2"] for r in seed_results])
    x = np.arange(len(BASE_SEEDS))
    ax.bar(x - 0.2, d2_alg, 0.4, label="AlgNet", color="steelblue", alpha=0.8)
    ax.bar(x + 0.2, d2_mlp, 0.4, label="MLP",    color="coral",     alpha=0.8)

    # Mean lines
    ax.axhline(np.nanmean(d2_alg), color="steelblue", linestyle="--",
               linewidth=1.5, label=f"AlgNet mean={np.nanmean(d2_alg):.2f}")
    ax.axhline(np.nanmean(d2_mlp), color="coral",     linestyle="--",
               linewidth=1.5, label=f"MLP mean={np.nanmean(d2_mlp):.2f}")

    # Integer reference lines
    for n in range(1, 6):
        ax.axhline(n, color="gray", linestyle=":", alpha=0.4, linewidth=0.8)
        ax.text(len(BASE_SEEDS)-0.3, n+0.05, f"D={n}", fontsize=8, color="gray")

    ax.set_xlabel("Seed")
    ax.set_ylabel("Correlation dimension D₂")
    ax.set_title("D₂ (full trajectory)")
    ax.set_xticks(x)
    ax.set_xticklabels([str(s) for s in BASE_SEEDS])
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    fig.suptitle("Exp11: Grassberger-Procaccia Correlation Dimension D₂", fontsize=13, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "exp11_results.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved {path}")
    plt.close(fig)


def _plot_trajectories(seed_results: list):
    """PCA projections and trajectory displacement."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Row 0: PCA projection (seed 0)
    for ax, model_key, color, model_name in [
        (axes[0, 0], "alg", "steelblue", "AlgNet-1L"),
        (axes[0, 1], "mlp", "coral",     "MLP-1L"),
    ]:
        traj = seed_results[0][model_key]["trajectory"]
        mu   = traj.mean(axis=0)
        X    = traj - mu
        _, sv, Vt = np.linalg.svd(X, full_matrices=False)
        pts2 = X @ Vt[:2].T
        sc   = ax.scatter(pts2[:, 0], pts2[:, 1],
                          c=np.arange(len(pts2)), cmap="viridis", s=5, alpha=0.7)
        plt.colorbar(sc, ax=ax, label="epoch / RECORD_EVERY")
        ax.set_title(f"{model_name} — PCA projection (seed 0)")
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.grid(True, alpha=0.3)

    # PCA variance explained
    ax = axes[0, 2]
    for model_key, color, label in [("alg", "steelblue", "AlgNet"), ("mlp", "coral", "MLP")]:
        traj = seed_results[0][model_key]["trajectory"]
        mu   = traj.mean(axis=0)
        X    = traj - mu
        _, sv, _ = np.linalg.svd(X, full_matrices=False)
        var_exp = (sv**2) / (sv**2).sum()
        ax.plot(np.arange(1, len(var_exp)+1), np.cumsum(var_exp), marker=".", color=color, label=label)
    ax.set_xlabel("PCA component")
    ax.set_ylabel("Cumulative variance explained")
    ax.set_title("PCA variance (seed 0)")
    ax.axhline(0.95, color="gray", linestyle="--", alpha=0.7, label="95%")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 30)

    # Row 1: step-wise displacement
    ax = axes[1, 0]
    for model_key, color, label in [("alg", "steelblue", "AlgNet"), ("mlp", "coral", "MLP")]:
        traj   = seed_results[0][model_key]["trajectory"]
        steps  = np.linalg.norm(np.diff(traj, axis=0), axis=1)
        epochs = (np.arange(len(steps)) + 1) * RECORD_EVERY
        ax.semilogy(epochs, steps, color=color, alpha=0.7, label=label, linewidth=0.8)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Step displacement ||Δθ||₂")
    ax.set_title("Trajectory velocity (seed 0)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # D₂ full vs convergence comparison
    ax = axes[1, 1]
    d2_alg_f = np.nanmean([r["d2_alg_full"]["d2"] for r in seed_results])
    d2_alg_c = np.nanmean([r["d2_alg_conv"]["d2"] for r in seed_results])
    d2_mlp_f = np.nanmean([r["d2_mlp_full"]["d2"] for r in seed_results])
    d2_mlp_c = np.nanmean([r["d2_mlp_conv"]["d2"] for r in seed_results])
    labels_   = ["AlgNet\nfull", "AlgNet\nconv", "MLP\nfull", "MLP\nconv"]
    vals_     = [d2_alg_f, d2_alg_c, d2_mlp_f, d2_mlp_c]
    colors_   = ["steelblue", "steelblue", "coral", "coral"]
    alphas_   = [0.8, 0.4, 0.8, 0.4]
    for xi, (lab, val, col, alph) in enumerate(zip(labels_, vals_, colors_, alphas_)):
        ax.bar(xi, val, color=col, alpha=alph, edgecolor="black", linewidth=0.5)
        ax.text(xi, val + 0.05, f"{val:.2f}", ha="center", fontsize=9)
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels_, fontsize=8)
    ax.set_ylabel("Mean D₂")
    ax.set_title("D₂: full vs convergence phase")
    ax.set_ylim(0, max(vals_)*1.3 + 0.5)
    for n in range(1, 5):
        ax.axhline(n, color="gray", linestyle=":", alpha=0.4, linewidth=0.8)
    ax.grid(True, alpha=0.3, axis="y")

    # Val accuracy comparison
    ax = axes[1, 2]
    epochs_ax = np.arange(1, N_EPOCHS + 1)
    for model_key, color, label in [("alg", "steelblue", "AlgNet"), ("mlp", "coral", "MLP")]:
        acc_matrix = np.array([r[model_key]["val_accs"] for r in seed_results])
        acc_mean   = acc_matrix.mean(axis=0)
        acc_std    = acc_matrix.std(axis=0)
        ax.plot(epochs_ax, acc_mean, color=color, label=label, linewidth=1)
        ax.fill_between(epochs_ax, acc_mean - acc_std, acc_mean + acc_std,
                        alpha=0.15, color=color)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Validation accuracy")
    ax.set_title("Learning curves (mean ± std, N=5 seeds)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.suptitle("Exp11: Trajectories, PCA, and Fractal Dimension", fontsize=13, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "exp11_trajectory.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved {path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
#  Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    main()
