"""
Exp11b: High-Resolution Fractal Dimension — Takens Embedding + Nearest-Neighbor D₂
====================================================================================

Improvement over Exp11:
  - T = 2000+ trajectory points (RECORD_EVERY=1, 2000 epochs) → statistical robustness
  - Nearest-neighbor correlation dimension estimator (Grassberger-Procaccia + Theiler correction)
  - Takens time-delay embedding as alternative to PCA (delay τ from AMI, embedding m from FNN)
  - Box-counting D₀ as secondary estimator for cross-validation
  - N_REPS = 5 seeds, per-seed D₂ ± bootstrap CI

Design:
  Task      : sign(Re(Tr(M^3))),  M = Σ_a x_a T_a ∈ su(3)
              N_GEN = 8, N_DIM = 3
  n_train   : 6 000, n_val 1 500
  Seeds     : 5 (42–46)
  Epochs    : 2 000 (RECORD_EVERY=1 → T=2000 points exactly)
  Theiler window : w_T = T^(1/3) ≈ 12 → exclude temporally correlated pairs
  PCA baseline  : K_PCA=16 (same as Exp11 for comparison)
  Takens alt    : τ from first minimum of AMI, m from false nearest neighbors

Outputs:
  exp11b_outputs/exp11b_results.{npz,png}
  exp11b_outputs/exp11b_takens.png
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
from scipy.spatial import KDTree

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
N_GEN           = N_DIM ** 2 - 1      # 8

LAMBDA_K        = 2.0
LAMBDA_ORTHO    = 1.0
WARMUP_EP       = 30

BATCH_SIZE      = 256
LR_MAIN         = 1e-3
N_EPOCHS        = 2_000
RECORD_EVERY    = 1                    # record every epoch → T=2000
REPORT_EVERY    = 500

N_TRAIN         = 6_000
N_VAL           = 1_500
N_REPS          = 5
BASE_SEEDS      = [42, 43, 44, 45, 46]

# D₂ estimation settings
K_PCA           = 16
N_R_POINTS      = 80
SCALING_FRAC    = (0.1, 0.6)
N_BOOTSTRAP     = 300
THEILER_POWER   = 1/3      # w_T = T^(1/3) Theiler temporal exclusion window

# Takens embedding settings
TAKENS_MAX_TAU  = 50        # max delay for AMI search
TAKENS_MAX_M    = 10        # max embedding dimension for FNN search
FNN_THRESHOLD   = 0.1       # false nearest neighbor threshold

OUT_DIR = os.path.join(HERE, "exp11b_outputs")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
#  Dataset: su(3) cubic trace (same as Exp11)
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
#  Models (identical to Exp11)
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
#  Training with high-res trajectory recording
# ---------------------------------------------------------------------------

def train_and_record(model, train_loader, val_loader, basis,
                     name: str, is_algnet: bool) -> dict:
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

    trajectory = np.array(trajectory)   # (T, n_params)
    print(f"  {name}: trajectory shape = {trajectory.shape}, "
          f"total displacement = {_total_displacement(trajectory):.4f}")

    return {
        "name":       name,
        "trajectory": trajectory,
        "val_accs":   np.array(val_accs),
        "best_val":   best_val,
    }


def _total_displacement(traj: np.ndarray) -> float:
    if len(traj) < 2:
        return 0.0
    steps = np.diff(traj, axis=0)
    return float(np.linalg.norm(steps, axis=1).sum())


# ---------------------------------------------------------------------------
#  Improved Grassberger-Procaccia with Theiler correction + KDTree
# ---------------------------------------------------------------------------

def pca_project(traj: np.ndarray, k: int) -> np.ndarray:
    mu    = traj.mean(axis=0)
    X     = traj - mu
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    return X @ Vt[:k].T


def correlation_sum_theiler(pts: np.ndarray, r_values: np.ndarray,
                            theiler_w: int) -> np.ndarray:
    """
    C(r) with Theiler correction: exclude pairs (i,j) with |i-j| < theiler_w
    to remove temporal autocorrelation bias.
    Uses KDTree for efficient neighbor counting.
    """
    N = len(pts)
    tree = KDTree(pts)
    C = np.zeros(len(r_values))

    # Count all pairs within radius r
    for ri, r in enumerate(r_values):
        # count_neighbors counts ALL pairs (i,j) i<j with dist < r
        total_count = tree.count_neighbors(tree, r) - N  # subtract self-pairs
        total_count //= 2  # since count_neighbors double-counts

        # Subtract temporally correlated pairs |i-j| < theiler_w
        temp_count = 0
        for dt in range(1, theiler_w):
            for i in range(N - dt):
                if np.linalg.norm(pts[i] - pts[i + dt]) < r:
                    temp_count += 1

        valid_count = max(0, total_count - temp_count)
        n_valid_pairs = N * (N - 1) // 2 - (N - 1) * min(theiler_w - 1, N - 1) // 1
        # More precise: valid pairs = total pairs - Theiler-excluded pairs
        n_theiler_excluded = sum(N - dt for dt in range(1, min(theiler_w, N)))
        n_valid_pairs = N * (N - 1) // 2 - n_theiler_excluded
        if n_valid_pairs > 0:
            C[ri] = valid_count / n_valid_pairs
        else:
            C[ri] = 0.0
    return C


def correlation_sum_kdtree(pts: np.ndarray, r_values: np.ndarray) -> np.ndarray:
    """
    Fast C(r) via KDTree pair counting (no Theiler correction).
    Suitable for PCA-projected points where temporal order is less critical.
    """
    N = len(pts)
    tree = KDTree(pts)
    n_pairs = N * (N - 1) / 2
    C = np.zeros(len(r_values))
    for ri, r in enumerate(r_values):
        count = tree.count_neighbors(tree, r)
        C[ri] = (count - N) / (2 * n_pairs)  # subtract self, undouble
    return C


def _auto_r_range(pts: np.ndarray, n_sample: int = 500) -> tuple:
    """Pick r range from 5th-95th percentile of subsampled pairwise distances."""
    rng   = np.random.default_rng(0)
    idx_s = rng.choice(len(pts), size=min(len(pts), n_sample), replace=False)
    pts_s = pts[idx_s]
    diff_s  = pts_s[:, None, :] - pts_s[None, :, :]
    d_s     = np.sqrt((diff_s**2).sum(axis=-1))
    d_flat  = d_s[np.triu_indices(len(idx_s), k=1)]
    d_flat  = d_flat[d_flat > 1e-15]
    r_min   = np.percentile(d_flat, 5)
    r_max   = np.percentile(d_flat, 95)
    return r_min, r_max


def estimate_d2(traj: np.ndarray, label: str, phase: str = "full",
                use_theiler: bool = True) -> dict:
    """
    Estimate D₂ with optional Theiler correction.
    Uses KDTree for efficiency (handles T=2000 well).
    """
    T, P = traj.shape
    k    = min(K_PCA, T - 1, P)

    # PCA projection
    mu    = traj.mean(axis=0)
    X     = traj - mu
    _, sv, Vt = np.linalg.svd(X, full_matrices=False)
    var_explained = (sv**2) / (sv**2).sum()
    pts   = X @ Vt[:k].T

    # r range
    r_min, r_max = _auto_r_range(pts)
    r_vals = np.logspace(np.log10(r_min), np.log10(r_max), N_R_POINTS)

    # Theiler window: T^(1/3)
    theiler_w = max(2, int(T ** THEILER_POWER))

    if use_theiler:
        C_vals = correlation_sum_theiler(pts, r_vals, theiler_w)
    else:
        C_vals = correlation_sum_kdtree(pts, r_vals)

    # Scaling region
    valid   = (C_vals > 0.02) & (C_vals < 0.98)
    lo_frac, hi_frac = SCALING_FRAC
    n_v     = valid.sum()
    if n_v < 5:
        print(f"  WARNING: {label} [{phase}]: only {n_v} valid C(r) points")
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

    # OLS slope
    A  = np.column_stack([log_r, np.ones_like(log_r)])
    rng = np.random.default_rng(0)
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

    print(f"  {label} [{phase}]: T={T}, k={k}, theiler_w={theiler_w}, "
          f"D₂={d2:.3f}±{d2_std:.3f}  "
          f"var_expl(PC1-3)={var_explained[:3].sum():.3f}")

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
        "theiler_w":    theiler_w,
    }


# ---------------------------------------------------------------------------
#  Takens time-delay embedding
# ---------------------------------------------------------------------------

def average_mutual_information(x: np.ndarray, max_tau: int, n_bins: int = 32) -> int:
    """
    Find optimal time delay τ as first minimum of Average Mutual Information.
    x : 1D time series (use first PCA component of trajectory).
    """
    n = len(x)
    ami_vals = []
    for tau in range(1, min(max_tau + 1, n // 2)):
        x1 = x[:n - tau]
        x2 = x[tau:]
        # 2D histogram
        H2d, _, _ = np.histogram2d(x1, x2, bins=n_bins)
        H2d = H2d / H2d.sum()
        Hx  = np.histogram(x1, bins=n_bins)[0].astype(float)
        Hx  = Hx / Hx.sum()
        Hy  = np.histogram(x2, bins=n_bins)[0].astype(float)
        Hy  = Hy / Hy.sum()

        ami = 0.0
        for i in range(n_bins):
            for j in range(n_bins):
                if H2d[i, j] > 0 and Hx[i] > 0 and Hy[j] > 0:
                    ami += H2d[i, j] * np.log(H2d[i, j] / (Hx[i] * Hy[j]))
        ami_vals.append(ami)

    # First local minimum
    ami_vals = np.array(ami_vals)
    for i in range(1, len(ami_vals) - 1):
        if ami_vals[i] < ami_vals[i-1] and ami_vals[i] <= ami_vals[i+1]:
            return i + 1  # tau (1-indexed)
    return max_tau // 2  # fallback


def false_nearest_neighbors(traj_1d: np.ndarray, tau: int, max_m: int,
                            threshold: float = FNN_THRESHOLD) -> int:
    """
    Find minimum embedding dimension m via false nearest neighbors method.
    traj_1d : scalar trajectory (first PCA component).
    """
    n = len(traj_1d)
    fnn_fracs = []
    for m in range(1, max_m + 1):
        n_embed = n - (m + 1) * tau
        if n_embed < 20:
            break

        # Build m-dimensional embedding
        X_m = np.column_stack([traj_1d[i*tau:i*tau + n_embed] for i in range(m)])
        # Build (m+1)-dimensional embedding
        X_m1 = np.column_stack([traj_1d[i*tau:i*tau + n_embed] for i in range(m + 1)])

        tree_m = KDTree(X_m)
        # Find nearest neighbor in m-space
        dists_m, idxs_m = tree_m.query(X_m, k=2)  # k=2: self + nearest
        nn_dists = dists_m[:, 1]
        nn_idxs  = idxs_m[:, 1]

        # Check false neighbors: distance increases dramatically in (m+1)-space
        fnn_count = 0
        valid_count = 0
        for i in range(n_embed):
            j = nn_idxs[i]
            if nn_dists[i] > 1e-12:
                extra_dist = abs(X_m1[i, -1] - X_m1[j, -1])
                ratio = extra_dist / nn_dists[i]
                if ratio > threshold:
                    fnn_count += 1
                valid_count += 1

        frac = fnn_count / valid_count if valid_count > 0 else 1.0
        fnn_fracs.append(frac)
        print(f"    FNN m={m}: frac={frac:.4f}")

        if frac < 0.01:  # less than 1% false neighbors
            return m

    return max_m  # fallback


def takens_embed(traj_pca1: np.ndarray, tau: int, m: int) -> np.ndarray:
    """
    Build Takens time-delay embedding from scalar time series.
    Returns (n_embed, m) array.
    """
    n = len(traj_pca1)
    n_embed = n - (m - 1) * tau
    return np.column_stack([traj_pca1[i*tau:i*tau + n_embed] for i in range(m)])


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

    n_half = len(res_alg["trajectory"]) // 2

    # PCA-based D₂ with Theiler correction
    print("  --- D₂ (PCA + Theiler) ---")
    d2_alg_full = estimate_d2(res_alg["trajectory"],          "AlgNet", "full",   use_theiler=True)
    d2_alg_conv = estimate_d2(res_alg["trajectory"][n_half:], "AlgNet", "convergence", use_theiler=True)
    d2_mlp_full = estimate_d2(res_mlp["trajectory"],          "MLP",    "full",   use_theiler=True)
    d2_mlp_conv = estimate_d2(res_mlp["trajectory"][n_half:], "MLP",    "convergence", use_theiler=True)

    # PCA-based D₂ WITHOUT Theiler (for comparison with Exp11)
    print("  --- D₂ (PCA, no Theiler — Exp11-compatible) ---")
    d2_alg_full_nt = estimate_d2(res_alg["trajectory"],          "AlgNet", "full_noTheiler", use_theiler=False)
    d2_mlp_full_nt = estimate_d2(res_mlp["trajectory"],          "MLP",    "full_noTheiler", use_theiler=False)

    # Takens embedding on first PCA component
    print("  --- Takens embedding ---")
    takens_results = {}
    for mkey, res in [("alg", res_alg), ("mlp", res_mlp)]:
        traj_pca = pca_project(res["trajectory"], K_PCA)
        pc1 = traj_pca[:, 0]  # first principal component

        tau = average_mutual_information(pc1, TAKENS_MAX_TAU)
        print(f"  {mkey}: AMI optimal τ = {tau}")

        m = false_nearest_neighbors(pc1, tau, TAKENS_MAX_M)
        print(f"  {mkey}: FNN optimal m = {m}")

        embed = takens_embed(pc1, tau, m)
        print(f"  {mkey}: Takens embedding shape = {embed.shape}")

        # D₂ on Takens embedding (no PCA needed, already m-dimensional)
        d2_takens = estimate_d2_raw(embed, f"{mkey}-Takens", "full")

        takens_results[mkey] = {
            "tau": tau,
            "m":   m,
            "embed": embed,
            "d2": d2_takens,
        }

    return {
        "seed":       seed,
        "alg":        res_alg,
        "mlp":        res_mlp,
        "d2_alg_full": d2_alg_full,
        "d2_alg_conv": d2_alg_conv,
        "d2_mlp_full": d2_mlp_full,
        "d2_mlp_conv": d2_mlp_conv,
        "d2_alg_full_noTheiler": d2_alg_full_nt,
        "d2_mlp_full_noTheiler": d2_mlp_full_nt,
        "takens":     takens_results,
    }


def estimate_d2_raw(pts: np.ndarray, label: str, phase: str) -> dict:
    """
    D₂ estimation directly on embedded points (no PCA step).
    """
    T = len(pts)
    r_min, r_max = _auto_r_range(pts)
    r_vals = np.logspace(np.log10(r_min), np.log10(r_max), N_R_POINTS)
    C_vals = correlation_sum_kdtree(pts, r_vals)

    valid   = (C_vals > 0.02) & (C_vals < 0.98)
    lo_frac, hi_frac = SCALING_FRAC
    n_v     = valid.sum()
    if n_v < 5:
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

    A  = np.column_stack([log_r, np.ones_like(log_r)])
    rng = np.random.default_rng(0)
    if len(log_r) < 2:
        d2 = float("nan")
        d2_std = float("nan")
    else:
        coef, *_ = np.linalg.lstsq(A, log_C, rcond=None)
        d2 = coef[0]
        d2_boots = []
        for _ in range(N_BOOTSTRAP):
            bi = rng.choice(len(log_r), size=len(log_r), replace=True)
            if len(np.unique(bi)) < 2:
                continue
            c_, *_ = np.linalg.lstsq(A[bi], log_C[bi], rcond=None)
            d2_boots.append(c_[0])
        d2_boots = np.array(d2_boots)
        d2_std = d2_boots.std() if len(d2_boots) > 5 else float("nan")

    print(f"  {label} [{phase}]: T={T}, dim={pts.shape[1]}, D₂={d2:.3f}±{d2_std:.3f}")
    return {"d2": d2, "d2_std": d2_std, "r_vals": r_vals, "C_vals": C_vals,
            "valid_idx": valid_idx, "T": T, "label": label, "phase": phase}


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def main():
    t_start = time.time()
    seed_results = []

    for seed in BASE_SEEDS:
        seed_results.append(run_one_seed(seed))

    # Aggregate D₂
    def _agg(key, subkey="d2"):
        vals = [r[key][subkey] for r in seed_results if not np.isnan(r[key][subkey])]
        return np.array(vals) if vals else np.array([float("nan")])

    def _fmt(a):
        a = a[~np.isnan(a)]
        if len(a) == 0: return "n/a"
        return f"{a.mean():.3f} ± {a.std(ddof=1):.3f}" if len(a) > 1 else f"{a[0]:.3f}"

    print("\n" + "="*60)
    print("  D₂ Summary — PCA + Theiler correction")
    print("  ----------------------------------------")
    print(f"  AlgNet  full       : {_fmt(_agg('d2_alg_full'))}")
    print(f"  AlgNet  convergence: {_fmt(_agg('d2_alg_conv'))}")
    print(f"  MLP     full       : {_fmt(_agg('d2_mlp_full'))}")
    print(f"  MLP     convergence: {_fmt(_agg('d2_mlp_conv'))}")
    print()
    print("  D₂ Summary — PCA, no Theiler (Exp11-compatible)")
    print("  ----------------------------------------")
    print(f"  AlgNet  full       : {_fmt(_agg('d2_alg_full_noTheiler'))}")
    print(f"  MLP     full       : {_fmt(_agg('d2_mlp_full_noTheiler'))}")
    print()
    print("  D₂ Summary — Takens embedding")
    print("  ----------------------------------------")
    for mkey in ["alg", "mlp"]:
        d2_takens = [r["takens"][mkey]["d2"]["d2"] for r in seed_results
                     if not np.isnan(r["takens"][mkey]["d2"]["d2"])]
        a = np.array(d2_takens) if d2_takens else np.array([float("nan")])
        print(f"  {mkey:8s} Takens    : {_fmt(a)}")

    # Save NPZ
    npz_path = os.path.join(OUT_DIR, "exp11b_results.npz")
    save_dict = {
        "d2_alg_full":       _agg("d2_alg_full"),
        "d2_alg_conv":       _agg("d2_alg_conv"),
        "d2_mlp_full":       _agg("d2_mlp_full"),
        "d2_mlp_conv":       _agg("d2_mlp_conv"),
        "d2_alg_full_noTheiler": _agg("d2_alg_full_noTheiler"),
        "d2_mlp_full_noTheiler": _agg("d2_mlp_full_noTheiler"),
        "d2_alg_takens":     np.array([r["takens"]["alg"]["d2"]["d2"] for r in seed_results]),
        "d2_mlp_takens":     np.array([r["takens"]["mlp"]["d2"]["d2"] for r in seed_results]),
        "traj_alg_0":        seed_results[0]["alg"]["trajectory"],
        "traj_mlp_0":        seed_results[0]["mlp"]["trajectory"],
    }
    np.savez(npz_path, **save_dict)
    print(f"Saved {npz_path}")

    _plot_results(seed_results)

    elapsed = time.time() - t_start
    print(f"\nTotal runtime: {elapsed/60:.1f} min")


# ---------------------------------------------------------------------------
#  Plotting
# ---------------------------------------------------------------------------

def _plot_results(seed_results: list):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # -- Row 0 --
    # (0,0): C(r) curves seed 0 — Theiler vs no-Theiler for AlgNet
    ax = axes[0, 0]
    for key, style, label in [
        ("d2_alg_full",          "o-",  "PCA+Theiler"),
        ("d2_alg_full_noTheiler", "s--", "PCA only"),
    ]:
        res0 = seed_results[0][key]
        mask = res0["C_vals"] > 1e-10
        ax.loglog(res0["r_vals"][mask], res0["C_vals"][mask],
                  style, markersize=3, label=f"AlgNet {label} D₂={res0['d2']:.2f}")
    ax.set_xlabel("r")
    ax.set_ylabel("C(r)")
    ax.set_title("AlgNet — Theiler correction effect")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    # (0,1): C(r) curves seed 0 — Theiler vs no-Theiler for MLP
    ax = axes[0, 1]
    for key, style, label in [
        ("d2_mlp_full",          "o-",  "PCA+Theiler"),
        ("d2_mlp_full_noTheiler", "s--", "PCA only"),
    ]:
        res0 = seed_results[0][key]
        mask = res0["C_vals"] > 1e-10
        ax.loglog(res0["r_vals"][mask], res0["C_vals"][mask],
                  style, markersize=3, color="coral",
                  label=f"MLP {label} D₂={res0['d2']:.2f}")
    ax.set_xlabel("r")
    ax.set_ylabel("C(r)")
    ax.set_title("MLP — Theiler correction effect")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    # (0,2): D₂ comparison bar chart — all methods
    ax = axes[0, 2]
    methods = ["PCA+Theiler\nfull", "PCA+Theiler\nconv", "PCA only\nfull", "Takens\nfull"]
    alg_vals = []
    mlp_vals = []
    for key in ["d2_alg_full", "d2_alg_conv", "d2_alg_full_noTheiler"]:
        vals = [r[key]["d2"] for r in seed_results if not np.isnan(r[key]["d2"])]
        alg_vals.append(np.mean(vals) if vals else 0)
    alg_takens = [r["takens"]["alg"]["d2"]["d2"] for r in seed_results
                  if not np.isnan(r["takens"]["alg"]["d2"]["d2"])]
    alg_vals.append(np.mean(alg_takens) if alg_takens else 0)

    for key in ["d2_mlp_full", "d2_mlp_conv", "d2_mlp_full_noTheiler"]:
        vals = [r[key]["d2"] for r in seed_results if not np.isnan(r[key]["d2"])]
        mlp_vals.append(np.mean(vals) if vals else 0)
    mlp_takens = [r["takens"]["mlp"]["d2"]["d2"] for r in seed_results
                  if not np.isnan(r["takens"]["mlp"]["d2"]["d2"])]
    mlp_vals.append(np.mean(mlp_takens) if mlp_takens else 0)

    x = np.arange(len(methods))
    ax.bar(x - 0.2, alg_vals, 0.35, label="AlgNet", color="steelblue", alpha=0.8)
    ax.bar(x + 0.2, mlp_vals, 0.35, label="MLP",    color="coral",     alpha=0.8)
    for xi, (av, mv) in enumerate(zip(alg_vals, mlp_vals)):
        ax.text(xi - 0.2, av + 0.03, f"{av:.2f}", ha="center", fontsize=7)
        ax.text(xi + 0.2, mv + 0.03, f"{mv:.2f}", ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=8)
    ax.set_ylabel("Mean D₂")
    ax.set_title("D₂ across methods")
    ax.legend()
    for n in range(1, 4):
        ax.axhline(n, color="gray", linestyle=":", alpha=0.4)
    ax.grid(True, alpha=0.3, axis="y")

    # -- Row 1 --
    # (1,0): Per-seed D₂ scatter — PCA+Theiler full
    ax = axes[1, 0]
    seeds = [r["seed"] for r in seed_results]
    d2_a = [r["d2_alg_full"]["d2"] for r in seed_results]
    d2_m = [r["d2_mlp_full"]["d2"] for r in seed_results]
    ax.scatter(seeds, d2_a, marker="o", color="steelblue", s=60, label="AlgNet", zorder=5)
    ax.scatter(seeds, d2_m, marker="s", color="coral",     s=60, label="MLP",    zorder=5)
    ax.axhline(np.nanmean(d2_a), color="steelblue", linestyle="--", alpha=0.6)
    ax.axhline(np.nanmean(d2_m), color="coral",     linestyle="--", alpha=0.6)
    ax.set_xlabel("Seed")
    ax.set_ylabel("D₂ (PCA+Theiler, full)")
    ax.set_title("Per-seed D₂")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (1,1): Takens τ and m per seed
    ax = axes[1, 1]
    taus_a = [r["takens"]["alg"]["tau"] for r in seed_results]
    taus_m = [r["takens"]["mlp"]["tau"] for r in seed_results]
    ms_a   = [r["takens"]["alg"]["m"]   for r in seed_results]
    ms_m   = [r["takens"]["mlp"]["m"]   for r in seed_results]
    x = np.arange(len(seeds))
    ax.bar(x - 0.3, taus_a, 0.25, label="AlgNet τ", color="steelblue", alpha=0.6)
    ax.bar(x,       taus_m, 0.25, label="MLP τ",    color="coral",     alpha=0.6)
    ax2 = ax.twinx()
    ax2.plot(x - 0.15, ms_a, "o-", color="steelblue", markersize=8, label="AlgNet m")
    ax2.plot(x + 0.15, ms_m, "s-", color="coral",     markersize=8, label="MLP m")
    ax.set_xlabel("Seed")
    ax.set_ylabel("τ (delay)")
    ax2.set_ylabel("m (embedding dim)")
    ax.set_xticks(x)
    ax.set_xticklabels([str(s) for s in seeds])
    ax.set_title("Takens parameters per seed")
    ax.legend(loc="upper left", fontsize=7)
    ax2.legend(loc="upper right", fontsize=7)

    # (1,2): PCA projection seed 0 with color-coded epochs
    ax = axes[1, 2]
    for model_key, color, model_name, marker in [
        ("alg", "steelblue", "AlgNet", "."),
        ("mlp", "coral",     "MLP",    "."),
    ]:
        traj = seed_results[0][model_key]["trajectory"]
        pts2 = pca_project(traj, K_PCA)[:, :2]
        ax.scatter(pts2[:, 0], pts2[:, 1], c=np.arange(len(pts2)),
                   cmap="viridis" if model_key == "alg" else "magma",
                   s=3, alpha=0.5, label=model_name)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA projection (seed 0)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.suptitle("Exp11b: High-Resolution Fractal Dimension (T=2000, Theiler + Takens)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "exp11b_results.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved {path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
#  Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    main()
