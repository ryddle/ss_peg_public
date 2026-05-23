"""
Exp14: Head Width Sweep — Stiefel Bottleneck Threshold h*
==========================================================

Question:
  What is the minimum head width h* at which AlgNet overcomes the
  Stiefel bottleneck (Δ(h) transitions from negative to non-negative)?
  Is the bottleneck caused by dimensionality or nonlinear expressiveness?

Design:
  Task      : sign(Re(Tr(M^3))),  M = Σ_a x_a T_a ∈ su(3)
              N_GEN = 8, N_DIM = 3, complex anti-Hermitian generators
  n_train   : 6 000, n_val 1 500
  N_reps    : 10 seeds per (model, h, head_type) combination
  Epochs    : 1 500, LR 1e-3, grad-clip max_norm=1.0, warmup 30 ep
  λ_K       : 2.0 (compact), λ_ortho 1.0

Head architectures (two variants):
  A) 1-layer: Linear(N_GEN, h) → Linear(h, 2)           [capacity only]
  B) 2-layer: Linear(N_GEN, h) → ReLU → Linear(h, 2)    [capacity + nonlinearity]
  h ∈ {2, 4, 8, 16, 32}

Models per (h, head_type):
  AlgNet-1L  — AlgebraicLayer(T_K+T_ortho) → ReLU → head_{A or B}
  MLP-1L     — Linear(8,8) → ReLU → head_{A or B}

Hypothesis:
  ∃ h* such that Δ(h) = AlgNet − MLP changes sign from negative to ≥0.
  If h*_A > h*_B, the bottleneck is nonlinear expressiveness (ReLU in head
  helps AlgNet more than extra width), confirming that the Stiefel constraint
  limits the representable function class, not just output dimensionality.

Outputs:
  exp14_outputs/exp14_results.{npz,png}
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
N_DIM       = 3
N_GEN       = N_DIM ** 2 - 1   # 8

LAMBDA_K    = 2.0
LAMBDA_ORTHO = 1.0
WARMUP_EP   = 30

BATCH_SIZE  = 256
LR_MAIN     = 1e-3
N_EPOCHS    = 1_500
REPORT_EVERY = 500

N_TRAIN     = 6_000
N_VAL       = 1_500
N_REPS      = 10
BASE_SEED   = 42

HEAD_WIDTHS    = [2, 4, 8, 16, 32]
HEAD_TYPES     = ["1L", "2L"]   # 1-layer (linear), 2-layer (linear+ReLU+linear)

OUT_DIR = os.path.join(HERE, "exp14_outputs")
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
#  Head modules
# ---------------------------------------------------------------------------

def make_head_1L(n_in: int, h: int) -> nn.Sequential:
    """1-layer head: Linear(n_in, h) → Linear(h, 2). No intermediate nonlinearity."""
    return nn.Sequential(
        nn.Linear(n_in, h, dtype=torch.float64),
        nn.Linear(h, 2, dtype=torch.float64),
    )


def make_head_2L(n_in: int, h: int) -> nn.Sequential:
    """2-layer head: Linear(n_in, h) → ReLU → Linear(h, 2)."""
    return nn.Sequential(
        nn.Linear(n_in, h, dtype=torch.float64),
        nn.ReLU(),
        nn.Linear(h, 2, dtype=torch.float64),
    )


def make_head(n_in: int, h: int, head_type: str) -> nn.Sequential:
    if head_type == "1L":
        return make_head_1L(n_in, h)
    elif head_type == "2L":
        return make_head_2L(n_in, h)
    else:
        raise ValueError(f"Unknown head_type: {head_type}")


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


class AlgNet(nn.Module):
    def __init__(self, h: int, head_type: str):
        super().__init__()
        self.layer = AlgebraicLayer()
        self.relu  = nn.ReLU()
        self.head  = make_head(N_GEN, h, head_type)

    def forward(self, x, basis):
        feat = self.relu(self.layer(x, basis))
        return self.head(feat)

    def total_reg(self, basis, scale):
        return (LAMBDA_K * self.layer.t_k(basis) + LAMBDA_ORTHO * self.layer.t_ortho()) * scale


class MLPNet(nn.Module):
    def __init__(self, h: int, head_type: str):
        super().__init__()
        self.fc1  = nn.Linear(N_GEN, N_GEN, dtype=torch.float64)
        self.relu = nn.ReLU()
        self.head = make_head(N_GEN, h, head_type)

    def forward(self, x, basis=None):
        feat = self.relu(self.fc1(x))
        return self.head(feat)


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
#  Training
# ---------------------------------------------------------------------------

def train_model(model, train_loader, val_loader, basis,
                is_algnet: bool, label: str) -> float:
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()), lr=LR_MAIN
    )

    best_val = 0.0
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

        val_acc = evaluate(model, val_loader, basis)
        if val_acc > best_val:
            best_val = val_acc

        if ep % REPORT_EVERY == 0 or ep == N_EPOCHS:
            print(f"      ep {ep:4d}  val={val_acc:.4f}  best={best_val:.4f}")

    return best_val


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def main():
    t_start = time.time()

    # Pre-generate shared datasets (different seed for data vs model init)
    basis = _su_basis(N_DIM, DEVICE)
    X_tr, y_tr = generate_su3_cubic_dataset(N_TRAIN, seed=999)
    X_vl, y_vl = generate_su3_cubic_dataset(N_VAL,   seed=1999)
    tr_loader = DataLoader(TensorDataset(X_tr, y_tr), batch_size=BATCH_SIZE, shuffle=True)
    vl_loader = DataLoader(TensorDataset(X_vl, y_vl), batch_size=512)

    # Results storage: results[head_type][h] = {"alg": [...], "mlp": [...]}
    results = {}

    for head_type in HEAD_TYPES:
        results[head_type] = {}
        for h in HEAD_WIDTHS:
            results[head_type][h] = {"alg": [], "mlp": []}
            print(f"\n{'='*60}")
            print(f"  Head type={head_type}, width h={h}")
            print(f"{'='*60}")

            for rep in range(N_REPS):
                seed = BASE_SEED + rep
                print(f"\n  --- Rep {rep+1}/{N_REPS} (seed={seed}) ---")

                # AlgNet
                torch.manual_seed(seed * 1000 + 1)
                algnet = AlgNet(h, head_type).to(DEVICE)
                n_params_alg = sum(p.numel() for p in algnet.parameters())
                if rep == 0:
                    print(f"    AlgNet params: {n_params_alg}")

                best_alg = train_model(algnet, tr_loader, vl_loader, basis,
                                       is_algnet=True, label=f"AlgNet-{head_type}-h{h}")
                results[head_type][h]["alg"].append(best_alg)

                # MLP
                torch.manual_seed(seed * 1000 + 2)
                mlp = MLPNet(h, head_type).to(DEVICE)
                n_params_mlp = sum(p.numel() for p in mlp.parameters())
                if rep == 0:
                    print(f"    MLP params: {n_params_mlp}")

                best_mlp = train_model(mlp, tr_loader, vl_loader, basis,
                                       is_algnet=False, label=f"MLP-{head_type}-h{h}")
                results[head_type][h]["mlp"].append(best_mlp)

                print(f"    AlgNet={best_alg:.4f}  MLP={best_mlp:.4f}  "
                      f"Δ={best_alg - best_mlp:+.4f}")

    # -----------------------------------------------------------------------
    #  Print summary tables
    # -----------------------------------------------------------------------
    print("\n" + "="*70)
    print("  RESULTS SUMMARY")
    print("="*70)

    for head_type in HEAD_TYPES:
        print(f"\n  Head type: {head_type} "
              f"({'Linear→Linear' if head_type == '1L' else 'Linear→ReLU→Linear'})")
        print(f"  {'h':>4s}  {'AlgNet':>14s}  {'MLP':>14s}  {'Δ':>8s}  {'p':>8s}  {'sig':>4s}")
        print(f"  {'-'*4}  {'-'*14}  {'-'*14}  {'-'*8}  {'-'*8}  {'-'*4}")

        for h in HEAD_WIDTHS:
            alg_vals = np.array(results[head_type][h]["alg"])
            mlp_vals = np.array(results[head_type][h]["mlp"])
            delta = alg_vals.mean() - mlp_vals.mean()
            t_stat, p_val = stats.ttest_ind(alg_vals, mlp_vals, equal_var=False)
            sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "n.s."
            print(f"  {h:4d}  {alg_vals.mean():.4f}±{alg_vals.std():.4f}  "
                  f"{mlp_vals.mean():.4f}±{mlp_vals.std():.4f}  "
                  f"{delta:+.4f}  {p_val:.4f}  {sig}")

    # -----------------------------------------------------------------------
    #  Determine h* for each head type
    # -----------------------------------------------------------------------
    h_star = {}
    for head_type in HEAD_TYPES:
        h_star[head_type] = None
        for h in HEAD_WIDTHS:
            alg_vals = np.array(results[head_type][h]["alg"])
            mlp_vals = np.array(results[head_type][h]["mlp"])
            delta = alg_vals.mean() - mlp_vals.mean()
            if delta >= 0:
                h_star[head_type] = h
                break
        print(f"\n  h* ({head_type}): {h_star[head_type] if h_star[head_type] else '>32 (not found)'}")

    if h_star["1L"] and h_star["2L"]:
        if h_star["1L"] > h_star["2L"]:
            print("  → h*_1L > h*_2L: bottleneck is NONLINEAR expressiveness")
        elif h_star["1L"] == h_star["2L"]:
            print("  → h*_1L = h*_2L: bottleneck is DIMENSIONAL (width-limited)")
        else:
            print("  → h*_1L < h*_2L: unexpected — ReLU in head HURTS AlgNet advantage")

    # -----------------------------------------------------------------------
    #  Save NPZ
    # -----------------------------------------------------------------------
    npz_path = os.path.join(OUT_DIR, "exp14_results.npz")
    save_dict = {}
    for head_type in HEAD_TYPES:
        for h in HEAD_WIDTHS:
            save_dict[f"alg_{head_type}_h{h}"] = np.array(results[head_type][h]["alg"])
            save_dict[f"mlp_{head_type}_h{h}"] = np.array(results[head_type][h]["mlp"])
    np.savez(npz_path, **save_dict)
    print(f"\nSaved {npz_path}")

    # -----------------------------------------------------------------------
    #  Plotting
    # -----------------------------------------------------------------------
    _plot_results(results, h_star)

    elapsed = time.time() - t_start
    print(f"\nTotal runtime: {elapsed/60:.1f} min")


# ---------------------------------------------------------------------------
#  Plotting
# ---------------------------------------------------------------------------

def _plot_results(results: dict, h_star: dict):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Top row: accuracy vs h for each head type
    for col, head_type in enumerate(HEAD_TYPES):
        ax = axes[0, col]
        h_arr = np.array(HEAD_WIDTHS)
        alg_means = []
        alg_stds  = []
        mlp_means = []
        mlp_stds  = []
        for h in HEAD_WIDTHS:
            a = np.array(results[head_type][h]["alg"])
            m = np.array(results[head_type][h]["mlp"])
            alg_means.append(a.mean())
            alg_stds.append(a.std(ddof=1))
            mlp_means.append(m.mean())
            mlp_stds.append(m.std(ddof=1))

        alg_means = np.array(alg_means)
        alg_stds  = np.array(alg_stds)
        mlp_means = np.array(mlp_means)
        mlp_stds  = np.array(mlp_stds)

        ax.errorbar(h_arr, alg_means, yerr=alg_stds, marker="o", capsize=4,
                    color="steelblue", label="AlgNet", linewidth=1.5)
        ax.errorbar(h_arr, mlp_means, yerr=mlp_stds, marker="s", capsize=4,
                    color="coral", label="MLP", linewidth=1.5)

        if h_star[head_type]:
            ax.axvline(h_star[head_type], color="green", linestyle="--", alpha=0.7,
                       label=f"h*={h_star[head_type]}")

        head_label = "Linear→Linear" if head_type == "1L" else "Linear→ReLU→Linear"
        ax.set_xlabel("Head width h")
        ax.set_ylabel("Best validation accuracy")
        ax.set_title(f"Head type {head_type}: {head_label}")
        ax.set_xscale("log", base=2)
        ax.set_xticks(HEAD_WIDTHS)
        ax.set_xticklabels([str(h) for h in HEAD_WIDTHS])
        ax.legend()
        ax.grid(True, alpha=0.3)

    # Bottom-left: Δ(h) for both head types
    ax = axes[1, 0]
    for head_type, color, style in [("1L", "purple", "o-"), ("2L", "darkgreen", "s--")]:
        deltas = []
        delta_cis = []
        for h in HEAD_WIDTHS:
            a = np.array(results[head_type][h]["alg"])
            m = np.array(results[head_type][h]["mlp"])
            d = a.mean() - m.mean()
            # SE of difference
            se = np.sqrt(a.var(ddof=1)/len(a) + m.var(ddof=1)/len(m))
            deltas.append(d)
            delta_cis.append(1.96 * se)
        deltas = np.array(deltas)
        delta_cis = np.array(delta_cis)
        head_label = "1-layer" if head_type == "1L" else "2-layer"
        ax.errorbar(HEAD_WIDTHS, deltas, yerr=delta_cis, fmt=style,
                    color=color, capsize=4, linewidth=1.5, label=f"Δ ({head_label} head)")

    ax.axhline(0, color="black", linewidth=1, linestyle="-")
    ax.set_xlabel("Head width h")
    ax.set_ylabel("Δ = AlgNet − MLP")
    ax.set_title("Δ(h) — sign transition identifies h*")
    ax.set_xscale("log", base=2)
    ax.set_xticks(HEAD_WIDTHS)
    ax.set_xticklabels([str(h) for h in HEAD_WIDTHS])
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Bottom-right: p-value heatmap
    ax = axes[1, 1]
    p_matrix = np.zeros((len(HEAD_TYPES), len(HEAD_WIDTHS)))
    for ri, head_type in enumerate(HEAD_TYPES):
        for ci, h in enumerate(HEAD_WIDTHS):
            a = np.array(results[head_type][h]["alg"])
            m = np.array(results[head_type][h]["mlp"])
            _, p_val = stats.ttest_ind(a, m, equal_var=False)
            p_matrix[ri, ci] = p_val

    im = ax.imshow(p_matrix, cmap="RdYlGn", vmin=0, vmax=0.1, aspect="auto")
    ax.set_xticks(range(len(HEAD_WIDTHS)))
    ax.set_xticklabels([str(h) for h in HEAD_WIDTHS])
    ax.set_yticks(range(len(HEAD_TYPES)))
    ax.set_yticklabels([f"{ht} head" for ht in HEAD_TYPES])
    ax.set_xlabel("Head width h")
    ax.set_title("p-value (Welch t-test, N=10)")
    for ri in range(len(HEAD_TYPES)):
        for ci in range(len(HEAD_WIDTHS)):
            txt = f"{p_matrix[ri, ci]:.3f}"
            color = "white" if p_matrix[ri, ci] < 0.03 else "black"
            ax.text(ci, ri, txt, ha="center", va="center", fontsize=9, color=color)
    fig.colorbar(im, ax=ax, label="p-value")

    fig.suptitle("Exp14: Head Width Sweep — Stiefel Bottleneck h*\n"
                 "su(3) cubic, N_GEN=8, depth-1 AlgNet vs MLP, N=10 reps",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "exp14_results.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved {path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
#  Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    main()
