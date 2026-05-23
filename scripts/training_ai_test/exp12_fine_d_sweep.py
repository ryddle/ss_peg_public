"""
Exp12: Fine d-Sweep — Algebraic Advantage Δ(d) with Error Bars
===============================================================

Question:
  Where is the true peak of algebraic advantage Δ(d) = acc(AlgNet) − acc(MLP)?
  Current data (d=3,4,5,6,8 from Exp5) shows non-monotonic behaviour.
  A fine sweep d=2..12 with N=10 reps per point resolves the peak shape
  and distinguishes smooth interpolation from possible fractal structure.

Design:
  Task      : sign(Re(Tr(M^3))),  M = Σ_a x_a T_a ∈ su(d)
              N_GEN = d^2 - 1,  complex anti-Hermitian generators
  d range   : 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
  n_train   : 6 000  (fixed across all d for comparability)
  n_val     : 1 500
  N_reps    : 10 seeds per (d, model)
  Epochs    : max(500, 4000 // d)   — adaptive to GPU budget
  LR        : 1e-3, grad-clip max_norm=1.0
  λ_K       : 2.0  (compact killing_tension, targets K ∝ −I)
  λ_ortho   : 1.0  (T_ortho = ||θθᵀ − I||_F², manifold confinement)
  Warmup    : 30 epochs for regularisation scaling

Architecture (same as Exp10, depth=1):
  AlgNet-1L : AlgebraicLayer(T_K + T_ortho) + ReLU → Linear(N_GEN, 2)
  MLP-1L    : Linear(N_GEN, N_GEN) + ReLU → Linear(N_GEN, 2)

Primary metric:
  Δ(d) = mean(AlgNet acc) − mean(MLP acc),  error bars = pooled SEM
  p-value from unpaired two-sided t-test (α=0.05)

Secondary metrics:
  T_K final,  θ orthogonality error,  K-signature string

Outputs:
  exp12_outputs/exp12_results.png    — Δ(d) plot with error bars
  exp12_outputs/exp12_results.npz    — full results array
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
from scipy import stats

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
D_VALUES        = list(range(2, 13))        # d = 2, 3, ..., 12
LAMBDA_K        = 2.0
LAMBDA_ORTHO    = 1.0
WARMUP_EP       = 30

BATCH_SIZE      = 256
LR_MAIN         = 1e-3
REPORT_EVERY    = 200

N_TRAIN         = 6_000
N_VAL           = 1_500
N_REPS          = 10
BASE_SEEDS      = list(range(N_REPS))       # seeds 0..9

OUT_DIR = os.path.join(HERE, "exp12_outputs")
os.makedirs(OUT_DIR, exist_ok=True)


def adaptive_epochs(d: int) -> int:
    """Epoch budget scaled to keep run-time roughly constant across d."""
    return max(500, 4_000 // d)


# ---------------------------------------------------------------------------
#  Dataset: su(d) cubic trace
# ---------------------------------------------------------------------------

def generate_su_cubic_dataset(d: int, n: int, seed: int, device):
    """
    y = sign(Re(Tr(M^3))),  M = Σ_a x_a T_a ∈ su(d),  x ~ N(0,1)^(d²-1).
    Returns balanced-class X, y on 'device'.
    """
    n_gen = d * d - 1
    rng   = np.random.default_rng(seed)
    X     = rng.standard_normal((n, n_gen)).astype(np.float64)
    basis = _su_basis(d, "cpu").numpy()             # (n_gen, d, d) complex128
    M     = np.einsum("na,aij->nij", X, basis)      # (n, d, d)
    M3    = M @ M @ M
    tr3   = np.trace(M3, axis1=-2, axis2=-1).real   # (n,)
    y     = (tr3 > 0).astype(np.int64)
    pos   = y.sum()
    print(f"  su({d}) cubic  n={n}  {n-pos} neg / {pos} pos  (frac={pos/n:.3f})")
    return (
        torch.tensor(X,  dtype=torch.float64, device=device),
        torch.tensor(y,  dtype=torch.long,    device=device),
    )


# ---------------------------------------------------------------------------
#  Models
# ---------------------------------------------------------------------------

class AlgebraicLayer(nn.Module):
    """
    One algebraic projection layer for su(d) with n_gen = d²-1 generators.
    Generator : G_a = Σ_b theta_ab B_b   (B = su(d) basis, complex).
    Forward   : h_a = Re(Tr(G_a M(x))) + bias_a
    Regs      : T_ortho = ||theta theta^T − I||_F^2
                T_K     = compact killing_tension(G)
    """

    def __init__(self, d: int):
        super().__init__()
        self.d     = d
        self.n_gen = d * d - 1
        init_t = torch.eye(self.n_gen, dtype=torch.float64)
        self.theta = nn.Parameter(init_t + torch.randn_like(init_t) * 0.05)
        self.bias  = nn.Parameter(torch.zeros(self.n_gen, dtype=torch.float64))

    def _get_generators(self, basis: torch.Tensor) -> torch.Tensor:
        th = self.theta.to(basis.dtype)
        return torch.einsum("ab,bij->aij", th, basis)

    def forward(self, x: torch.Tensor, basis: torch.Tensor) -> torch.Tensor:
        G   = self._get_generators(basis)
        M_x = torch.einsum("...a,aij->...ij", x.to(basis.dtype), basis)
        f   = torch.einsum("aij,...ji->...a", G, M_x)
        return f.real.to(torch.float64) + self.bias

    def t_ortho(self) -> torch.Tensor:
        gram = self.theta @ self.theta.T
        I_n  = torch.eye(self.n_gen, dtype=torch.float64, device=self.theta.device)
        return ((gram - I_n) ** 2).sum()

    def t_k(self, basis: torch.Tensor) -> torch.Tensor:
        return killing_tension(self._get_generators(basis))

    def k_signature(self, basis: torch.Tensor) -> str:
        with torch.no_grad():
            K    = killing_metric(self._get_generators(basis)).cpu().numpy()
            if np.iscomplexobj(K):
                K = K.real
            eigs = np.linalg.eigvalsh(K)
            thr  = 0.1 * abs(eigs).max() if abs(eigs).max() > 1e-10 else 1e-8
            pos  = int((eigs >  thr).sum())
            neg  = int((eigs < -thr).sum())
            zer  = self.n_gen - pos - neg
        return f"+{pos} -{neg} ~{zer}"


class AlgNet1L(nn.Module):
    def __init__(self, d: int):
        super().__init__()
        n_gen = d * d - 1
        self.layer = AlgebraicLayer(d)
        self.relu  = nn.ReLU()
        self.head  = nn.Linear(n_gen, 2, dtype=torch.float64)

    def forward(self, x, basis):
        return self.head(self.relu(self.layer(x, basis)))

    def total_reg(self, basis, scale: float):
        return (
            LAMBDA_K     * self.layer.t_k(basis)
            + LAMBDA_ORTHO * self.layer.t_ortho()
        ) * scale


class MLP1L(nn.Module):
    def __init__(self, d: int):
        super().__init__()
        n_gen = d * d - 1
        self.net  = nn.Sequential(
            nn.Linear(n_gen, n_gen, dtype=torch.float64),
            nn.ReLU(),
        )
        self.head = nn.Linear(n_gen, 2, dtype=torch.float64)

    def forward(self, x, basis=None):
        return self.head(self.net(x))


# ---------------------------------------------------------------------------
#  Evaluation
# ---------------------------------------------------------------------------

def evaluate(model, loader, basis) -> float:
    model.eval()
    correct, total = 0, 0
    with torch.no_grad():
        for xb, yb in loader:
            out  = model(xb, basis)
            correct += (out.argmax(1) == yb).sum().item()
            total   += xb.size(0)
    return correct / total


# ---------------------------------------------------------------------------
#  Training loop for one (model, d, seed)
# ---------------------------------------------------------------------------

def train_one(model, train_loader, val_loader, basis, name: str,
              is_algnet: bool, n_epochs: int) -> dict:
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=LR_MAIN
    )

    best_val = 0.0
    final_tk = 0.0

    for ep in range(1, n_epochs + 1):
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

        val_acc = evaluate(model, val_loader, basis)
        if val_acc > best_val:
            best_val = val_acc

        if ep % REPORT_EVERY == 0 or ep == n_epochs:
            if is_algnet:
                with torch.no_grad():
                    final_tk = model.layer.t_k(basis).item()
            print(f"    ep {ep:4d}/{n_epochs}  val={val_acc:.4f}  best={best_val:.4f}"
                  + (f"  T_K={final_tk:.5f}" if is_algnet else ""))

    k_sig = model.layer.k_signature(basis) if is_algnet else "n/a"
    ortho_err = _ortho_error(model) if is_algnet else float("nan")
    return {
        "name":      name,
        "best_val":  best_val,
        "final_val": val_acc,
        "final_tk":  final_tk,
        "k_sig":     k_sig,
        "ortho_err": ortho_err,
    }


def _ortho_error(model) -> float:
    """||theta theta^T − I||_F."""
    with torch.no_grad():
        t = model.layer.theta
        gram = (t @ t.T).cpu().numpy()
        return float(np.linalg.norm(gram - np.eye(gram.shape[0])))


# ---------------------------------------------------------------------------
#  Run one (d, seed) pair
# ---------------------------------------------------------------------------

def run_single(d: int, seed: int, n_epochs: int):
    print(f"  d={d}  seed={seed}  epochs={n_epochs}")

    # Shared dataset
    X_tr, y_tr = generate_su_cubic_dataset(d, N_TRAIN, seed=seed,     device=DEVICE)
    X_vl, y_vl = generate_su_cubic_dataset(d, N_VAL,   seed=seed+100, device=DEVICE)
    basis       = _su_basis(d, DEVICE)

    tr_loader = DataLoader(TensorDataset(X_tr, y_tr), batch_size=BATCH_SIZE, shuffle=True)
    vl_loader = DataLoader(TensorDataset(X_vl, y_vl), batch_size=512)

    torch.manual_seed(seed * 1000 + 1)
    algnet = AlgNet1L(d).to(DEVICE)
    torch.manual_seed(seed * 1000 + 2)
    mlp    = MLP1L(d).to(DEVICE)

    print("  --- AlgNet ---")
    res_alg = train_one(algnet, tr_loader, vl_loader, basis,
                        f"AlgNet d={d}", is_algnet=True, n_epochs=n_epochs)
    print("  --- MLP ---")
    res_mlp = train_one(mlp, tr_loader, vl_loader, basis,
                        f"MLP d={d}",    is_algnet=False, n_epochs=n_epochs)
    return res_alg, res_mlp


# ---------------------------------------------------------------------------
#  Main: sweep over d and N_REPS seeds
# ---------------------------------------------------------------------------

def main():
    t_start = time.time()

    # results[d][seed] = (res_alg, res_mlp)
    all_results = {}

    for d in D_VALUES:
        n_epochs = adaptive_epochs(d)
        print(f"\n{'='*60}")
        print(f"  d = {d}  (N_GEN={d*d-1}, epochs={n_epochs})")
        print(f"{'='*60}")
        alg_accs = []
        mlp_accs = []
        for seed in BASE_SEEDS:
            res_alg, res_mlp = run_single(d, seed, n_epochs)
            alg_accs.append(res_alg["best_val"])
            mlp_accs.append(res_mlp["best_val"])
            print(f"    seed={seed}  AlgNet={res_alg['best_val']:.4f}"
                  f"  MLP={res_mlp['best_val']:.4f}"
                  f"  Δ={res_alg['best_val']-res_mlp['best_val']:+.4f}")
        all_results[d] = {"alg": alg_accs, "mlp": mlp_accs}

        delta_mean = np.mean(alg_accs) - np.mean(mlp_accs)
        t_stat, p_val = stats.ttest_ind(alg_accs, mlp_accs)
        print(f"  >>> d={d}  Δ_mean={delta_mean:+.4f}  p={p_val:.4f}")

    # Summary table
    print("\n" + "="*70)
    print(f"  {'d':>4}  {'AlgNet':>8}  {'MLP':>8}  {'Δ':>8}  {'p':>8}  {'sig':>5}")
    print("  " + "-"*60)
    for d in D_VALUES:
        alg = all_results[d]["alg"]
        mlp = all_results[d]["mlp"]
        dm  = np.mean(alg) - np.mean(mlp)
        _, p = stats.ttest_ind(alg, mlp)
        sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "n.s."))
        print(f"  {d:>4}  {np.mean(alg):>8.4f}  {np.mean(mlp):>8.4f}"
              f"  {dm:>+8.4f}  {p:>8.4f}  {sig:>5}")

    # Save NPZ
    npz_path = os.path.join(OUT_DIR, "exp12_results.npz")
    np.savez(
        npz_path,
        d_values   = np.array(D_VALUES),
        alg_matrix = np.array([all_results[d]["alg"] for d in D_VALUES]),  # (11, 10)
        mlp_matrix = np.array([all_results[d]["mlp"] for d in D_VALUES]),  # (11, 10)
    )
    print(f"\nSaved {npz_path}")

    # Plot
    _plot(all_results)

    elapsed = time.time() - t_start
    print(f"\nTotal runtime: {elapsed/60:.1f} min")


def _plot(all_results: dict):
    d_arr  = np.array(D_VALUES)
    alg_m  = np.array([np.mean(all_results[d]["alg"]) for d in D_VALUES])
    alg_s  = np.array([np.std(all_results[d]["alg"], ddof=1) / np.sqrt(N_REPS) for d in D_VALUES])
    mlp_m  = np.array([np.mean(all_results[d]["mlp"]) for d in D_VALUES])
    mlp_s  = np.array([np.std(all_results[d]["mlp"], ddof=1) / np.sqrt(N_REPS) for d in D_VALUES])
    delta  = alg_m - mlp_m
    delta_s = np.sqrt(alg_s**2 + mlp_s**2)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: raw accuracy vs d
    ax = axes[0]
    ax.errorbar(d_arr, alg_m, yerr=alg_s, marker="o", label="AlgNet-1L", color="steelblue", capsize=4)
    ax.errorbar(d_arr, mlp_m, yerr=mlp_s, marker="s", label="MLP-1L",    color="coral",     capsize=4)
    ax.set_xlabel("d  (su(d) dimension)")
    ax.set_ylabel("Best validation accuracy")
    ax.set_title("Accuracy vs d  (N=10 seeds)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xticks(d_arr)

    # Right: Δ = AlgNet − MLP vs d
    ax = axes[1]
    bar_colors = ["steelblue" if v >= 0 else "coral" for v in delta]
    ax.bar(d_arr, delta, color=bar_colors, alpha=0.7, edgecolor="black", linewidth=0.5)
    ax.errorbar(d_arr, delta, yerr=delta_s, fmt="none", color="black", capsize=4)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel("d  (su(d) dimension)")
    ax.set_ylabel("Δ = AlgNet − MLP  (best val acc)")
    ax.set_title("Algebraic advantage Δ(d)  (N=10)")
    ax.grid(True, alpha=0.3, axis="y")
    ax.set_xticks(d_arr)

    # Annotate peak
    peak_d = d_arr[np.argmax(delta)]
    ax.annotate(f"peak d={peak_d}", xy=(peak_d, delta.max()),
                xytext=(peak_d + 0.5, delta.max() + 0.003),
                fontsize=9, arrowprops=dict(arrowstyle="->"))

    fig.suptitle("Exp12: Fine d-Sweep — Algebraic Advantage Δ(d)", fontsize=13, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "exp12_results.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved {path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
#  Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    main()
