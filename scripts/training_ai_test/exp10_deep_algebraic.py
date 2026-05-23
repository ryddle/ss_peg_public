"""
Exp10: Deep Algebraic Network — Stacked Killing Layers on su(3)
================================================================

Question:
  Can multiple algebraic layers (each regularized by T_K + T_ortho) be stacked?
  Does depth amplify the accuracy advantage of algebraic structure over plain MLP?
  Is the Killing signature preserved at every layer by T_ortho (Sylvester invariance)?

Design:
  Task      : sign(Re(Tr(M^3))),  M = Σ_a x_a T_a ∈ su(3)
              N_GEN = 8, N_DIM = 3, complex anti-Hermitian generators
  n_train   : 6 000  (same as Exp9)
  n_val     : 1 500
  N_reps    : 3 seeds per model
  Epochs    : 1 500, LR 1e-3, grad-clip max_norm=1.0
  λ_K       : 2.0 per layer  (compact killing_tension, targets K ∝ −I)
  λ_ortho   : 1.0 per layer

Architecture:
  Feature layers : depth k of (AlgebraicLayer or Linear(8,8)) + ReLU
  Head           : Linear(8, 2) [same for ALL models]

Models (7 total):
  MLP-1L     — 1 × Linear(8,8)+ReLU → head
  MLP-2L     — 2 × Linear(8,8)+ReLU → head    (stacked)
  MLP-3L     — 3 × Linear(8,8)+ReLU → head    (stacked)
  AlgNet-1L  — 1 × AlgLayer(T_K+T_ortho)+ReLU → head
  AlgNet-2L  — 2 × AlgLayer(T_K+T_ortho)+ReLU → head  (stacked)
  AlgNet-3L  — 3 × AlgLayer(T_K+T_ortho)+ReLU → head  (stacked)
  Oracle     — 3 × AlgLayer(theta=I fixed)   +ReLU → head

Hypothesis (G14):
  a) AlgNet-kL ≥ MLP-kL at each depth k (algebraic bias helps at every depth)
  b) Δ(AlgNet-kL, MLP-kL) non-decreasing in k (depth amplifies structural advantage)
  c) K-signature +0 -8 ~0 preserved at each layer by T_ortho (Sylvester argument:
     theta ∈ O(8) → K(G) = theta K_ref theta^T has same signature as K_ref)
  d) Oracle upper-bounds AlgNet-3L (or is exceeded as in Exp9)

Outputs:
  exp10_outputs/exp10_results.png
  exp10_outputs/exp10_results.npz
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
N_GEN       = N_DIM ** 2 - 1   # 8 generators for su(3)
DEPTHS      = [1, 2, 3]

LAMBDA_K      = 2.0   # compact T_K per layer  (su(3) ∝ −I target)
LAMBDA_ORTHO  = 1.0   # T_ortho per layer
WARMUP_EP     = 30

BATCH_SIZE    = 256
LR_MAIN       = 1e-3
N_EPOCHS      = 1500
REPORT_EVERY  = 300

N_TRAIN     = 6_000
N_VAL       = 1_500
BASE_SEEDS  = [42, 43, 44]
N_REPS      = 3

OUT_DIR = os.path.join(HERE, "exp10_outputs")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
#  Dataset: su(3) cubic trace
# ---------------------------------------------------------------------------

def generate_su3_cubic_dataset(n: int, seed: int):
    """
    y = sign(Re(Tr(M^3))),  M = Σ_a x_a T_a ∈ su(3),  x ~ N(0,1)^8.
    Class-balanced by antisymmetry under x → −x.
    """
    rng   = np.random.default_rng(seed)
    X     = rng.standard_normal((n, N_GEN)).astype(np.float64)
    basis = _su_basis(N_DIM, "cpu").numpy()        # (8, 3, 3) complex128
    M     = np.einsum("na,aij->nij", X, basis)     # (n, 3, 3) complex
    M3    = M @ M @ M
    tr3   = np.trace(M3, axis1=-2, axis2=-1).real  # (n,)
    y     = (tr3 > 0).astype(np.int64)
    pos   = y.sum()
    print(f"  su(3) cubic  n={n}  {n-pos} neg / {pos} pos  (frac={pos/n:.3f})")
    return (
        torch.tensor(X, dtype=torch.float64, device=DEVICE),
        torch.tensor(y, dtype=torch.long,    device=DEVICE),
    )

# ---------------------------------------------------------------------------
#  Algebraic Layer
# ---------------------------------------------------------------------------

class AlgebraicLayer(nn.Module):
    """
    One real-parametrized algebraic projection layer.

    Generators  G_a = Σ_b theta_ab B_b   (B = su(3) basis, complex).
    Forward     h_a = Re(Tr(G_a M(x))) + bias_a
                where M(x) = Σ_a x_a B_a  is the Lie-algebra element.
    T_ortho     ||theta theta^T − I||_F^2        (manifold confinement)
    T_K         compact killing_tension(G)        (targets K ∝ −I)

    If fixed=True, theta is frozen to I (Oracle mode: exact su(3) generators).
    """

    def __init__(self, fixed: bool = False):
        super().__init__()
        init_t = torch.eye(N_GEN, dtype=torch.float64)
        if fixed:
            self.register_buffer("theta", init_t)
        else:
            self.theta = nn.Parameter(init_t + torch.randn_like(init_t) * 0.05)
        self.bias  = nn.Parameter(torch.zeros(N_GEN, dtype=torch.float64))
        self.fixed = fixed

    def _get_generators(self, basis: torch.Tensor) -> torch.Tensor:
        """G = theta @ basis  →  (N_GEN, N_DIM, N_DIM) complex."""
        th = self.theta.to(basis.dtype)            # float64 → complex128 (zero imag)
        return torch.einsum("ab,bij->aij", th, basis)

    def forward(self, x: torch.Tensor, basis: torch.Tensor) -> torch.Tensor:
        """
        x     : (..., N_GEN) real — Lie algebra coordinates
        basis : (N_GEN, N_DIM, N_DIM) complex — fixed reference generators
        """
        G   = self._get_generators(basis)                             # (8, 3, 3) complex
        M_x = torch.einsum("...a,aij->...ij", x.to(basis.dtype), basis)  # (..., 3, 3)
        f   = torch.einsum("aij,...ji->...a", G, M_x)                 # (..., 8) complex
        return f.real.to(torch.float64) + self.bias                   # (..., 8) real

    def t_ortho(self) -> torch.Tensor:
        if self.fixed:
            return self.bias.new_zeros(())
        gram = self.theta @ self.theta.T
        I_n  = torch.eye(N_GEN, dtype=torch.float64, device=self.theta.device)
        return ((gram - I_n) ** 2).sum()

    def t_k(self, basis: torch.Tensor) -> torch.Tensor:
        if self.fixed:
            return self.bias.new_zeros(())
        G = self._get_generators(basis)
        return killing_tension(G)

    def k_signature(self, basis: torch.Tensor) -> str:
        """Eigenvalue signature string of the Killing form K(G)."""
        if self.fixed:
            return "n/a"
        with torch.no_grad():
            G    = self._get_generators(basis)
            K    = killing_metric(G).cpu().numpy()
            if np.iscomplexobj(K):
                K = K.real
            eigs = np.linalg.eigvalsh(K)
            thr  = 0.1 * abs(eigs).max() if abs(eigs).max() > 1e-10 else 1e-8
            pos  = int((eigs >  thr).sum())
            neg  = int((eigs < -thr).sum())
            zer  = N_GEN - pos - neg
        return f"+{pos} -{neg} ~{zer}"

# ---------------------------------------------------------------------------
#  Models
# ---------------------------------------------------------------------------

class AlgNetDeep(nn.Module):
    """
    Depth-k algebraic network.
    k AlgebraicLayers each regularized by T_K + T_ortho, interleaved with ReLU,
    followed by a linear classifier head.
    """

    def __init__(self, depth: int, fixed: bool = False):
        super().__init__()
        self.alg_layers = nn.ModuleList(
            [AlgebraicLayer(fixed=fixed) for _ in range(depth)]
        )
        self.relu   = nn.ReLU()
        self.head   = nn.Linear(N_GEN, 2, dtype=torch.float64)
        self.depth  = depth

    def forward(self, x: torch.Tensor, basis: torch.Tensor) -> torch.Tensor:
        h = x
        for layer in self.alg_layers:
            h = self.relu(layer(h, basis))
        return self.head(h)

    def total_reg(self, basis: torch.Tensor, scale: float) -> torch.Tensor:
        """Sum of (lambda_K * T_K + lambda_ortho * T_ortho) over all layers."""
        total = torch.zeros((), dtype=torch.float64, device=self.head.weight.device)
        for layer in self.alg_layers:
            total = total + (
                LAMBDA_K    * layer.t_k(basis)
                + LAMBDA_ORTHO * layer.t_ortho()
            ) * scale
        return total

    def k_signatures(self, basis: torch.Tensor) -> list:
        return [lyr.k_signature(basis) for lyr in self.alg_layers]


class MLPDeep(nn.Module):
    """
    Depth-k plain MLP.
    k × Linear(8,8)+ReLU layers then a linear classifier head.
    Identical parameter count structure to AlgNetDeep(k).
    """

    def __init__(self, depth: int):
        super().__init__()
        blocks = []
        for _ in range(depth):
            blocks += [nn.Linear(N_GEN, N_GEN, dtype=torch.float64), nn.ReLU()]
        self.net   = nn.Sequential(*blocks)
        self.head  = nn.Linear(N_GEN, 2, dtype=torch.float64)
        self.depth = depth

    def forward(self, x: torch.Tensor, basis: torch.Tensor = None) -> torch.Tensor:
        return self.head(self.net(x))

# ---------------------------------------------------------------------------
#  Evaluation
# ---------------------------------------------------------------------------

def _evaluate(model: nn.Module, loader, basis: torch.Tensor) -> float:
    model.eval()
    correct, total = 0, 0
    with torch.no_grad():
        for xb, yb in loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            correct += (model(xb, basis).argmax(1) == yb).sum().item()
            total   += xb.size(0)
    return correct / total

# ---------------------------------------------------------------------------
#  Training
# ---------------------------------------------------------------------------

def train_one(
    model    : nn.Module,
    t_loader,
    v_loader,
    basis    : torch.Tensor,
    name     : str,
    is_algnet: bool,
) -> dict:
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=LR_MAIN,
    )
    best_val = 0.0
    hist_val = []

    print(f"\n  >> {name}  (depth={model.depth}, {N_EPOCHS}ep"
          f"  λ_K={LAMBDA_K if is_algnet else 0:.1f}"
          f"  λ_o={LAMBDA_ORTHO if is_algnet else 0:.1f})")

    for epoch in range(N_EPOCHS):
        scale = min(1.0, (epoch + 1) / WARMUP_EP)
        model.train()

        for xb, yb in t_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            optimizer.zero_grad()
            loss = criterion(model(xb, basis), yb)
            if is_algnet:
                loss = loss + model.total_reg(basis, scale)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

        val_acc  = _evaluate(model, v_loader, basis)
        best_val = max(best_val, val_acc)
        hist_val.append(val_acc)

        if (epoch + 1) % REPORT_EVERY == 0:
            print(f"     ep{epoch+1:4d}  val={val_acc:.4f}  best={best_val:.4f}")

    sigs = model.k_signatures(basis) if is_algnet else ["n/a"] * model.depth
    print(f"     FINAL  best={best_val:.4f}  K-sigs={sigs}")
    return {"best": best_val, "hist": hist_val, "sigs": sigs}

# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def main():
    t0 = time.time()

    # Pre-compute su(3) basis on device (complex float64)
    basis = _su_basis(N_DIM, "cpu").to(DEVICE)

    # Shared validation set (deterministic seed)
    X_val, y_val = generate_su3_cubic_dataset(N_VAL, seed=999)
    val_loader   = DataLoader(TensorDataset(X_val, y_val),
                              batch_size=512, shuffle=False)

    # Print su(3) K_ref signature (theoretical: +0 -8 ~0)
    K_ref = killing_metric(basis)
    if K_ref.is_complex():
        K_ref = K_ref.real
    eigs_ref = np.linalg.eigvalsh(K_ref.cpu().numpy())
    pos_r = int((eigs_ref >  0.1 * abs(eigs_ref).max()).sum())
    neg_r = int((eigs_ref < -0.1 * abs(eigs_ref).max()).sum())
    zer_r = N_GEN - pos_r - neg_r
    print(f"\n  su(3) K_ref eigenvalue signature: +{pos_r} -{neg_r} ~{zer_r}")
    print(f"  K_ref eigenvalues: {np.round(eigs_ref, 2)}")

    # (name, ModelClass, depth, is_algnet, fixed)
    configs = [
        ("MLP-1L",    MLPDeep,    1, False, False),
        ("MLP-2L",    MLPDeep,    2, False, False),
        ("MLP-3L",    MLPDeep,    3, False, False),
        ("AlgNet-1L", AlgNetDeep, 1, True,  False),
        ("AlgNet-2L", AlgNetDeep, 2, True,  False),
        ("AlgNet-3L", AlgNetDeep, 3, True,  False),
        ("Oracle",    AlgNetDeep, 3, True,  True),
    ]

    results = {cfg[0]: [] for cfg in configs}

    for name, Cls, depth, is_algnet, fixed in configs:
        print(f"\n{'='*70}")
        print(f"  MODEL: {name}")
        print(f"{'='*70}")

        for seed in BASE_SEEDS:
            torch.manual_seed(seed)
            np.random.seed(seed)

            if Cls is AlgNetDeep:
                model = AlgNetDeep(depth=depth, fixed=fixed).to(DEVICE)
            else:
                model = MLPDeep(depth=depth).to(DEVICE)

            X_tr, y_tr = generate_su3_cubic_dataset(N_TRAIN, seed=seed)
            t_loader   = DataLoader(TensorDataset(X_tr, y_tr),
                                    batch_size=BATCH_SIZE, shuffle=True)

            run = train_one(model, t_loader, val_loader, basis,
                            f"{name}[s{seed}]", is_algnet)
            run["seed"] = seed
            results[name].append(run)

    # --------------- summary ---------------
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    summary = {}
    for name in results:
        bests = [r["best"] for r in results[name]]
        mean  = float(np.mean(bests))
        std   = float(np.std(bests, ddof=1))
        summary[name] = {"mean": mean, "std": std, "bests": bests}
        sigs_last = results[name][-1]["sigs"]
        print(f"  {name:12s}  mean={mean:.4f}±{std:.4f}  "
              f"bests={[round(b, 4) for b in bests]}")
        for idx, sig in enumerate(sigs_last):
            if sig != "n/a":
                print(f"     Layer {idx + 1} K-sig: {sig}")

    print("\nPairwise AlgNet vs MLP at each depth (Welch t-test):")
    for d in DEPTHS:
        alg_key, mlp_key = f"AlgNet-{d}L", f"MLP-{d}L"
        ab = summary[alg_key]["bests"]
        mb = summary[mlp_key]["bests"]
        delta  = np.mean(ab) - np.mean(mb)
        _, pval = stats.ttest_ind(ab, mb)
        print(f"  depth={d}  AlgNet={summary[alg_key]['mean']:.4f}  "
              f"MLP={summary[mlp_key]['mean']:.4f}  "
              f"Δ={delta:+.4f}  p={pval:.4f}")

    oracle_mean = summary["Oracle"]["mean"]
    alg3_mean   = summary["AlgNet-3L"]["mean"]
    print(f"\n  Oracle={oracle_mean:.4f}  AlgNet-3L={alg3_mean:.4f}  "
          f"Δ(AlgNet3 − Oracle)={alg3_mean - oracle_mean:+.4f}")

    elapsed = time.time() - t0
    print(f"\nWall time: {elapsed:.1f}s ({elapsed / 60:.1f} min)")

    # --------------- per-seed table ---------------
    print("\nPer-seed results:")
    header = f"  {'Model':12s}  " + "  ".join(f"s{s}" for s in BASE_SEEDS)
    print(header)
    for name in results:
        row = f"  {name:12s}  " + "  ".join(
            f"{r['best']:.4f}" for r in results[name]
        )
        print(row)

    # --------------- save ---------------
    model_names = list(summary.keys())
    np.savez(
        os.path.join(OUT_DIR, "exp10_results.npz"),
        model_names=np.array(model_names),
        means =np.array([summary[n]["mean"]  for n in model_names]),
        stds  =np.array([summary[n]["std"]   for n in model_names]),
        bests =np.array([summary[n]["bests"] for n in model_names]),
    )

    # --------------- plot ---------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    means_arr = [summary[n]["mean"] for n in model_names]
    stds_arr  = [summary[n]["std"]  for n in model_names]
    colors    = [
        "#bbbbbb", "#888888", "#444444",    # MLP greys (1L, 2L, 3L)
        "#ffaaaa", "#ff4444", "#cc0000",    # AlgNet reds (1L, 2L, 3L)
        "#3399ff",                          # Oracle blue
    ]
    xpos = np.arange(len(model_names))

    axes[0].bar(xpos, means_arr, yerr=stds_arr, color=colors, alpha=0.85, capsize=5)
    axes[0].set_xticks(xpos)
    axes[0].set_xticklabels(model_names, rotation=30, ha="right", fontsize=9)
    axes[0].set_ylabel("Best val accuracy")
    axes[0].set_title("Exp10: Depth × Algebraic bias  [su(3) cubic]")
    axes[0].set_ylim(0.5, 1.0)
    axes[0].axhline(oracle_mean, ls="--", color="#3399ff", alpha=0.5,
                    label=f"Oracle={oracle_mean:.3f}")
    axes[0].legend(fontsize=8)

    # Accuracy vs depth
    for label, marker, keys, color in [
        ("MLP",    "o", ["MLP-1L",    "MLP-2L",    "MLP-3L"],    "#666666"),
        ("AlgNet", "s", ["AlgNet-1L", "AlgNet-2L", "AlgNet-3L"], "#cc0000"),
    ]:
        ym = [summary[k]["mean"] for k in keys]
        ye = [summary[k]["std"]  for k in keys]
        axes[1].errorbar(DEPTHS, ym, yerr=ye, marker=marker, color=color,
                         capsize=5, label=label, linewidth=1.8)
    axes[1].axhline(oracle_mean, ls="--", color="#3399ff", alpha=0.5, label="Oracle")
    axes[1].set_xlabel("Depth  (# algebraic / linear layers before head)")
    axes[1].set_ylabel("Best val accuracy")
    axes[1].set_title("Accuracy vs Depth")
    axes[1].set_xticks(DEPTHS)
    axes[1].legend()

    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "exp10_results.png"), dpi=150)
    plt.close(fig)

    print(f"\nSaved: {OUT_DIR}/exp10_results.{{png,npz}}")


if __name__ == "__main__":
    main()
