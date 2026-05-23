"""
Experiment 4: SU(3)-Symmetric Cubic Invariant Task
=====================================================

Key changes from Exp3:
  Exp3 input:  x ∈ R^9 — arbitrary real 3×3 matrix (NOT an su(3) element)
  Exp4 input:  x ∈ R^8 — GML coordinates of Hermitian traceless M = Σ_a x_a T_a ∈ su(3)

  Exp3 task:  sign(Tr(M^3)) for arbitrary real M
  Exp4 task:  sign(Tr(M(x)^3)) where M(x) IS in su(3)  ← genuine su(3) cubic invariant

Mathematical structure:
  The su(3) cubic invariant is d_{abc} x_a x_b x_c ∝ Tr(M^3) where
  d_{abc} = (1/4) Tr({T_a, T_b} T_c) = the fully symmetric rank-3 su(3) d-tensor.

  Under any SU(3) adjoint action M → U M U†:
    Tr(M^3) = Tr((UMU†)^3) = Tr(M^3)   [invariant!]

  This means the task is genuinely symmetric under SU(3) ⊂ O(8):
  any rotation in the SU(3) orbit yields the same label.

  AlgebraicProjection2 forward:
    f_a = Re(Tr(G_a @ M(x)))   where M(x) = Σ_b x_b T_b

    = Re(Σ_b x_b Tr(G_a T_b))  [bilinear, linear in x]
    = Σ_b x_b · Re(Σ_k θ_{ak} Tr(T_k T_b))  [G_a = Σ_k θ_{ak} T_k]
    = Σ_b x_b · 2 θ_{ab}       [Tr(T_k T_b) = 2δ_{kb}]
    = 2 (θ @ x^T)^T

  → When θ = I (GML-Fixed oracle): f_a = 2 x_a   [EXACT coordinate recovery]
  → When T_K → 0:  θ ∈ O(8)  →  first layer is an isometric rotation of R^8
  → Without T_K:   θ may degenerate, losing rank, destroying information

Models:
  1. MLP-8     : free Linear(8→8) → ReLU → Linear(8→32) → ReLU → Linear(32→2)
  2. AlgNet-Free: AlgebraicProjection2 (θ free, no T_K) + same head
  3. AlgNet-TK  : AlgebraicProjection2 + T_K loss (λ_K = 2.0) + same head
  4. GML-Fixed  : FIXED θ = I (first layer not trained) + same head  ← oracle

Scientific hypothesis:
  GML-Fixed ≥ AlgNet-TK > AlgNet-Free ≈ MLP-8
  (with enough training, crystallized AlgNet approaches oracle)

Hardware: NVIDIA RTX 5070 Ti
Output  : exp4_results.npz + exp4_results.png + exp4_generators_*.png
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader, TensorDataset

from killing_utils import (
    DEVICE,
    killing_metric,
    killing_tension,
    closure_tension,
    analyze_generators,
    print_algebra_report,
    params_to_su_generators,
    _su_basis,
)

torch.set_default_dtype(torch.float64)

# ============================================================
#  HYPER-PARAMETERS
# ============================================================

D          = 3              # su(3)
N_GEN      = D**2 - 1      # 8 generators
N_IN       = N_GEN          # 8 input features (GML coordinates)

N_EPOCHS   = 1000
BATCH_SIZE = 256
LR         = 1e-3
N_SAMPLES  = 20_000

LAMBDA_K   = 2.0            # Killing tension (AlgNet-TK)
WARMUP     = 50             # Longer warmup — give task loss time to orient generators


def lam(epoch: int, lam_max: float, warmup: int = WARMUP) -> float:
    """Linear warmup schedule."""
    return lam_max * min(1.0, epoch / max(warmup, 1))


# ============================================================
#  DATASET: su(3)-symmetric cubic invariant
# ============================================================

def generate_su3_cubic_dataset(
    n:    int = N_SAMPLES,
    seed: int = 42,
) -> tuple:
    """
    SU(3)-symmetric classification task.

    Input:  x ∈ R^8  (coordinates in GML basis: M = Σ_a x_a T_a ∈ su(3))
    Label:  sign(Tr(M(x)^3)) — the su(3) cubic invariant d_{abc} x_a x_b x_c

    Key property:  Tr(M^3) ∈ R for Hermitian M  (verifiable since M†=M → M^3 Hermitian-ish)
    Note:          sign is +/- balanced because the d-tensor changes sign under some reflections.

    The task is genuinely invariant under SU(3) adjoint action on x.
    """
    rng = np.random.default_rng(seed)
    X   = rng.standard_normal((n, N_GEN)).astype(np.float64)  # (N, 8) GML coordinates

    # Build the basis on CPU (numpy)
    basis_t = _su_basis(D, "cpu")                          # (8, 3, 3) complex128 tensor
    basis   = basis_t.numpy()                              # numpy for fast batched ops

    # M = Σ_a x_a T_a  →  compute M^3  →  extract Tr(M^3)
    M_c  = np.einsum('na,aij->nij', X, basis)             # (N, 3, 3) complex
    M3_c = np.einsum('nij,njk,nkl->nil', M_c, M_c, M_c)  # (N, 3, 3) complex
    tr3  = np.trace(M3_c, axis1=-2, axis2=-1).real        # (N,) — imaginary part numerically 0

    y    = (tr3 > 0).astype(np.int64)

    # Sanity: print balance
    pos = y.sum(); neg = n - pos
    print(f"  Dataset balance: {neg} neg / {pos} pos  (frac={pos/n:.3f})")

    # Sanity: max imaginary residual should be tiny
    max_imag = np.abs(np.trace(M3_c, axis1=-2, axis2=-1).imag).max()
    if max_imag > 1e-10:
        print(f"  WARNING: max imag(Tr(M^3)) = {max_imag:.2e}  (should be ~0)")

    return (
        torch.tensor(X,   dtype=torch.float64, device=DEVICE),
        torch.tensor(y,   dtype=torch.long,    device=DEVICE),
    )


# ============================================================
#  MODELS
# ============================================================

def _make_head() -> nn.Sequential:
    """Shared head: 8 features → 2-class logits. Same for all models."""
    return nn.Sequential(
        nn.ReLU(),
        nn.Linear(N_GEN, 32),
        nn.ReLU(),
        nn.Linear(32, 2),
    )


class StandardMLP8(nn.Module):
    """
    Baseline: free Linear(8→8) + head.

    No algebraic constraint on the first layer — just learns whatever 8×8 matrix
    minimises the task loss. Can represent any linear combination of input coords.
    """

    def __init__(self):
        super().__init__()
        self.layer1 = nn.Linear(N_IN, N_GEN)
        self.head   = _make_head()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.head(self.layer1(x))


class AlgebraicProjection2(nn.Module):
    """
    Algebraic first layer: x ∈ R^8 → features f ∈ R^8 via su(3) inner products.

    M(x) = Σ_b x_b T_b        (always a valid su(3) Hermitian traceless element)
    f_a  = Re(Tr(G_a @ M(x))) (projection onto learned generator G_a)
         = 2 Σ_b θ_{ab} x_b   (reduces to a linear map with weight W = 2θ)

    The algebraic structure CONSTRAINS W: it must come from valid su(3) generators.
    T_K pressure drives θ → O(8) → first layer becomes an isometry (no info loss).

    When θ = I (oracle/GML-Fixed): f_a = 2 x_a  — exact coordinate recovery.
    """

    def __init__(self, fixed: bool = False, lambda_k: float = 0.0):
        super().__init__()
        self.d        = D
        self.n_gen    = N_GEN
        self.fixed    = fixed
        self.lambda_k = lambda_k

        # θ: coordinates of each learned generator in the GML basis
        # Init: identity = GML basis (T_K=0 starting point), small noise for free models
        init = torch.eye(N_GEN, dtype=torch.float64)
        if fixed:
            # Oracle: θ = I, frozen. Use register_buffer so it moves with .to(device)
            self.register_buffer("theta_buf", init)
        else:
            self.theta = nn.Parameter(
                init + torch.randn(N_GEN, N_GEN, dtype=torch.float64) * 0.05
            )

        self.bias = nn.Parameter(torch.zeros(N_GEN, dtype=torch.float64))

    def _theta(self) -> torch.Tensor:
        """Return θ regardless of fixed/free mode."""
        return self.theta_buf if self.fixed else self.theta

    def get_generators(self) -> torch.Tensor:
        """(n_gen, d, d) complex — always hermitian traceless.
        Uses the device of θ so the result is always on the correct device."""
        return params_to_su_generators(self._theta(), self.d)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        x:  (B, 8)  →  f: (B, 8)
        Computes f_a = Re(Tr(G_a @ M(x))) where M(x) = Σ_b x_b T_b.

        Device-agnostic: basis and generators are placed on the same device as θ.
        """
        dev   = str(self._theta().device)
        basis = _su_basis(self.d, dev)                                  # (8, 3, 3) complex

        # M(x) = Σ_b x_b T_b:  (B, 3, 3) complex
        M_x   = torch.einsum('...b,bij->...ij',
                              x.to(basis.device).to(torch.complex128), basis)

        # G_a ∈ (8, 3, 3) complex  (on device of θ)
        G     = self.get_generators()

        # f_a = Re(Tr(G_a @ M(x)))  →  shortcut:
        #   Tr(G_a M) = Σ_{ij} G_a[i,j] M[j,i] = einsum('aij,...ji->...a', G, M)
        #   NOTE: no conjugate on G — Tr(G_a M), not <G_a, M>_F
        f     = torch.einsum('aij,...ji->...a', G, M_x).real            # (B, 8)
        return f + self.bias

    def reg_loss(self, epoch: int = 0) -> torch.Tensor:
        """Killing tension with warmup schedule."""
        if self.fixed or self.lambda_k == 0.0:
            dev = str(self._theta().device)
            return torch.zeros(1, dtype=torch.float64, device=dev).squeeze()
        lk = lam(epoch, self.lambda_k)
        return lk * killing_tension(self.get_generators())


class AlgebraicNet2(nn.Module):
    """Neural net with AlgebraicProjection2 as Layer 1 + shared head."""

    def __init__(self, fixed: bool = False, lambda_k: float = 0.0):
        super().__init__()
        self.alg_layer = AlgebraicProjection2(fixed=fixed, lambda_k=lambda_k)
        self.head      = _make_head()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.head(self.alg_layer(x))

    def reg_loss(self, epoch: int = 0) -> torch.Tensor:
        return self.alg_layer.reg_loss(epoch)

    def get_generators(self) -> torch.Tensor:
        return self.alg_layer.get_generators().detach()


# ============================================================
#  TRAINING
# ============================================================

def train(
    model:        nn.Module,
    train_loader: DataLoader,
    val_loader:   DataLoader,
    name:         str = "Model",
    n_epochs:     int = N_EPOCHS,
) -> dict:
    model     = model.to(DEVICE)

    # GML-Fixed first layer has no parameters — only head will be optimised
    optimizer = optim.Adam(filter(lambda p: p.requires_grad, model.parameters()), lr=LR)
    criterion = nn.CrossEntropyLoss()

    hist = {
        "train_acc":  [], "val_acc":   [],
        "train_loss": [], "val_loss":  [],
        "T_K":        [], "T_closure": [],
    }

    best_val   = 0.0
    best_epoch = 0
    print(f"\n  [{name}]  {n_epochs} epochs ...")

    for epoch in range(n_epochs):
        model.train()
        t_loss, correct, total = 0.0, 0, 0
        for xb, yb in train_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            optimizer.zero_grad()
            logits    = model(xb)
            task_loss = criterion(logits, yb)
            loss      = task_loss
            if isinstance(model, AlgebraicNet2):
                loss = loss + model.reg_loss(epoch)
            loss.backward()
            optimizer.step()
            t_loss  += task_loss.item() * xb.size(0)
            correct += (logits.argmax(1) == yb).sum().item()
            total   += xb.size(0)

        model.eval()
        v_loss, v_correct, v_total = 0.0, 0, 0
        with torch.no_grad():
            for xb, yb in val_loader:
                xb, yb = xb.to(DEVICE), yb.to(DEVICE)
                out      = model(xb)
                v_loss   += criterion(out, yb).item() * xb.size(0)
                v_correct += (out.argmax(1) == yb).sum().item()
                v_total   += xb.size(0)

        val_acc = v_correct / v_total
        if val_acc > best_val:
            best_val   = val_acc
            best_epoch = epoch + 1

        # Algebraic diagnostics
        if isinstance(model, AlgebraicNet2):
            with torch.no_grad():
                G   = model.alg_layer.get_generators()
                t_k = killing_tension(G).item()
                t_c = closure_tension(G).item()
        else:
            t_k = t_c = float("nan")

        hist["train_loss"].append(t_loss / total)
        hist["train_acc"].append(correct / total)
        hist["val_loss"].append(v_loss / v_total)
        hist["val_acc"].append(val_acc)
        hist["T_K"].append(t_k)
        hist["T_closure"].append(t_c)

        if (epoch + 1) % 100 == 0:
            print(
                f"    epoch {epoch+1:4d}/{n_epochs}"
                f"  val_acc={val_acc:.4f}"
                f"  T_K={t_k:.4f}  T_cl={t_c:.4f}"
            )

    hist["best_val"]   = best_val
    hist["best_epoch"] = best_epoch
    return hist


# ============================================================
#  GENERATOR + THETA VISUALIZATION
# ============================================================

def plot_theta(theta: torch.Tensor, title: str, out_path: str) -> None:
    """Visualize the θ matrix (8×8 real) as a heatmap.
    Off-diagonal entries should be small when crystallized (θ ≈ I or orthogonal).
    """
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    mat = theta.cpu().float().numpy()
    im = ax.imshow(mat, cmap="RdBu", vmin=-1.5, vmax=1.5)
    plt.colorbar(im, ax=ax)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("GML basis index k"); ax.set_ylabel("Generator index a")
    plt.tight_layout()
    plt.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close()
    print(f"  Theta heatmap → {out_path}")


def plot_generators(generators: torch.Tensor, title: str, out_path: str) -> None:
    """Visualize the 8 su(3) generators as 3×3 heatmaps (real and imaginary parts)."""
    G   = generators.detach().cpu()
    fig, axes = plt.subplots(4, 4, figsize=(12, 12))
    fig.suptitle(title, fontsize=12)
    for a in range(N_GEN):
        ax = axes[a // 2][a % 2 * 2]
        ax.imshow(G[a].real.numpy(), cmap='RdBu', vmin=-1, vmax=1)
        ax.set_title(f"G_{a+1} Re", fontsize=9)
        ax.axis("off")
        ax = axes[a // 2][a % 2 * 2 + 1]
        ax.imshow(G[a].imag.numpy(), cmap='RdBu', vmin=-1, vmax=1)
        ax.set_title(f"G_{a+1} Im", fontsize=9)
        ax.axis("off")
    plt.tight_layout()
    plt.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close()
    print(f"  Generator figure → {out_path}")


# ============================================================
#  RESULTS PLOT
# ============================================================

def plot_results(results: dict, out_path: str) -> None:
    """5-panel figure: val accuracy, loss, T_K, T_closure, θ-orthogonality."""
    epochs     = range(1, N_EPOCHS + 1)
    model_keys = list(results.keys())

    colors = {
        "MLP-8":       "tab:blue",
        "AlgNet-Free": "tab:gray",
        "AlgNet-TK":   "tab:orange",
        "GML-Fixed":   "tab:green",
    }
    ls_map = {
        "MLP-8":       "-",
        "AlgNet-Free": "--",
        "AlgNet-TK":   "-",
        "GML-Fixed":   "-.",
    }

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(
        f"Experiment 4 — SU(3) Cubic Invariant Task  "
        f"(d={D}, n_gen={N_GEN},  λ_K={LAMBDA_K},  input=x∈R^8 GML coords)",
        fontsize=11,
    )

    # Val accuracy
    ax = axes[0, 0]
    for name in model_keys:
        h = results[name]
        ax.plot(epochs, h["val_acc"], label=name,
                color=colors.get(name, "black"), ls=ls_map.get(name, "-"), lw=1.8)
    ax.set_title("Validation Accuracy")
    ax.set_xlabel("Epoch"); ax.set_ylabel("Accuracy")
    ax.legend(); ax.grid(True, alpha=0.3)

    # Val loss
    ax = axes[0, 1]
    for name in model_keys:
        h = results[name]
        ax.plot(epochs, h["val_loss"], label=name,
                color=colors.get(name, "black"), ls=ls_map.get(name, "-"), lw=1.8)
    ax.set_title("Validation Loss (Cross-Entropy)")
    ax.set_xlabel("Epoch"); ax.set_ylabel("Loss")
    ax.legend(); ax.grid(True, alpha=0.3)

    # T_K (log scale)
    ax = axes[1, 0]
    for name in model_keys:
        h = results[name]
        vals = h["T_K"]
        if not all(np.isnan(v) for v in vals):
            ax.semilogy(list(epochs), vals, label=name,
                        color=colors.get(name, "black"), ls=ls_map.get(name, "-"), lw=1.8)
    ax.axhline(1e-3, ls=":", color="gray", alpha=0.5, label="T_K=0.001")
    ax.axvline(WARMUP, ls="--", color="gray", alpha=0.4, label=f"warmup={WARMUP}")
    ax.set_title("Killing Tension (log)  → 0 = Lie algebra crystallized")
    ax.set_xlabel("Epoch"); ax.set_ylabel("T_K")
    ax.legend(); ax.grid(True, alpha=0.3)

    # T_closure (log scale)
    ax = axes[1, 1]
    for name in model_keys:
        h = results[name]
        vals = h["T_closure"]
        if not all(np.isnan(v) for v in vals):
            ax.semilogy(list(epochs), [max(v, 1e-12) for v in vals], label=name,
                        color=colors.get(name, "black"), ls=ls_map.get(name, "-"), lw=1.8)
    ax.set_title("Closure Fraction (log)  → 0 = closed under commutator")
    ax.set_xlabel("Epoch"); ax.set_ylabel("T_closure")
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close()
    print(f"  Figure saved → {out_path}")


# ============================================================
#  MAIN
# ============================================================

if __name__ == "__main__":

    OUT_DIR   = os.path.dirname(os.path.abspath(__file__))
    BASE_NAME = os.path.join(OUT_DIR, "exp4")

    print("=" * 68)
    print("Experiment 4: SU(3)-Symmetric Cubic Invariant Task")
    print("=" * 68)
    print(f"  Device      : {DEVICE}")
    print(f"  GPU         : {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'}")
    print(f"  d={D}  n_gen={N_GEN}  n_in={N_IN}  epochs={N_EPOCHS}")
    print(f"  λ_K={LAMBDA_K}  warmup={WARMUP} epochs")
    print(f"  INPUT CHANGE vs Exp3: x ∈ R^8 = GML coords of su(3) element")
    print(f"  ORACLE MODEL: GML-Fixed (θ=I, first layer frozen → f_a = 2x_a)")

    # ----------------------------------------------------------
    # [1] Dataset
    # ----------------------------------------------------------
    print("\n[1/5] Generating su(3) cubic invariant dataset ...")
    print(f"  N={N_SAMPLES}, task: sign(Tr(M(x)^3))  where M(x) = Σ_a x_a T_a ∈ su(3)")
    X, y = generate_su3_cubic_dataset(N_SAMPLES, seed=42)
    print(f"  Input shape: {X.shape}  (R^8 GML coordinates)")

    n_train = int(0.8 * N_SAMPLES)
    X_tr, y_tr = X[:n_train],  y[:n_train]
    X_val, y_val = X[n_train:], y[n_train:]

    train_ds     = TensorDataset(X_tr, y_tr)
    val_ds       = TensorDataset(X_val, y_val)
    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
    val_loader   = DataLoader(val_ds,   batch_size=BATCH_SIZE)

    # ----------------------------------------------------------
    # [2] Verify GML-Fixed oracle analytically
    # ----------------------------------------------------------
    print("\n[2/5] Oracle sanity check ...")
    gml_proj = AlgebraicProjection2(fixed=True).to(DEVICE)
    gml_proj.eval()
    # For a sample batch: confirm f_a = 2 x_a (within numerical precision)
    x_test = X[:10]
    with torch.no_grad():
        f_test = gml_proj(x_test)  # (10, 8) on DEVICE
    diff_oracle = (f_test - 2 * x_test).abs().max().item()
    print(f"  max |f_a - 2 x_a| (oracle check) : {diff_oracle:.2e}  "
          f"{'✓' if diff_oracle < 1e-8 else '✗ WARNING'}")

    # Also check free model (θ ≈ I + noise): verify forward pass runs on DEVICE
    free_proj = AlgebraicProjection2(fixed=False, lambda_k=0.0).to(DEVICE)
    free_proj.eval()
    with torch.no_grad():
        f_free = free_proj(x_test)
    print(f"  Initial free θ deviation from I  : "
          f"{(free_proj.theta - torch.eye(N_GEN, device=DEVICE)).abs().max().item():.2e}")

    # Check that GML-Fixed has only bias as trainable parameter in alg layer
    n_fixed_params = sum(
        p.numel() for name, p in gml_proj.named_parameters() if "bias" not in name
    )
    n_bias_params = sum(
        p.numel() for name, p in gml_proj.named_parameters() if "bias" in name
    )
    print(f"  GML-Fixed non-bias trainable params: {n_fixed_params}  (expected 0)")

    # ----------------------------------------------------------
    # [3] Build models
    # ----------------------------------------------------------
    models_to_train = {
        "MLP-8":       StandardMLP8(),
        "AlgNet-Free": AlgebraicNet2(fixed=False, lambda_k=0.0),
        "AlgNet-TK":   AlgebraicNet2(fixed=False, lambda_k=LAMBDA_K),
        "GML-Fixed":   AlgebraicNet2(fixed=True,  lambda_k=0.0),
    }

    # Print parameter counts
    print("\n  Model parameter counts:")
    for name, m in models_to_train.items():
        n_params = sum(p.numel() for p in m.parameters())
        print(f"    {name:<14}  {n_params:5d} params")

    # ----------------------------------------------------------
    # [4] Train all models sequentially
    # ----------------------------------------------------------
    results = {}
    for i, (name, model) in enumerate(models_to_train.items(), 1):
        print(f"\n[{i+2}/5] Training {name} ...")
        results[name] = train(model, train_loader, val_loader, name)

    # ----------------------------------------------------------
    # [5] Save results and generate plots FIRST (before analysis that might fail)
    # ----------------------------------------------------------

    # Summary table
    print("\n" + "=" * 68)
    print("SUMMARY")
    print("=" * 68)
    print(f"  {'Model':<16}  {'Best acc':>9}  {'@ epoch':>8}  {'Final acc':>10}  {'T_K final':>10}")
    print(f"  {'-' * 62}")

    for name, h in results.items():
        best   = h.get("best_val",   max(h["val_acc"]))
        bepoch = h.get("best_epoch", int(np.argmax(h["val_acc"])) + 1)
        final  = h["val_acc"][-1]
        t_k    = h["T_K"][-1]
        t_k_s  = f"{t_k:.4f}" if not np.isnan(t_k) else "     —"
        print(f"  {name:<16}  {best:>9.4f}  {bepoch:>8d}  {final:>10.4f}  {t_k_s:>10}")

    ref = results["GML-Fixed"]["val_acc"][-1]
    print(f"\n  Oracle GML-Fixed final acc: {ref:.4f}")
    for name in ("MLP-8", "AlgNet-Free", "AlgNet-TK"):
        delta = results[name]["val_acc"][-1] - ref
        print(f"  Δ({name} vs GML-Fixed) : {delta:+.4f}")

    print("\n  Killing tension convergence (AlgebraicNet models):")
    for name in ("AlgNet-Free", "AlgNet-TK"):
        h  = results[name]
        t0 = h["T_K"][0]
        tb = min(h["T_K"])
        tf = h["T_K"][-1]
        print(f"    {name:<16}  initial={t0:.4f}  best={tb:.6f}  final={tf:.4f}")

    # Plots
    plot_results(results, BASE_NAME + "_results.png")

    # Numeric results
    save_data = {}
    for name, h in results.items():
        key = name.lower().replace(" ", "_").replace("-", "")
        for metric, vals in h.items():
            if isinstance(vals, list):
                save_data[f"{key}_{metric}"] = np.array(vals)
            else:
                save_data[f"{key}_{metric}_scalar"] = np.array(vals)
    np.savez(BASE_NAME + "_results.npz", **save_data)
    print(f"  Numeric results → {BASE_NAME}_results.npz")

    # ----------------------------------------------------------
    # [6] Post-training algebraic analysis (after saves — safe to fail)
    # ----------------------------------------------------------
    print("\n[6/6] Post-training algebraic analysis ...")

    for name in ("AlgNet-Free", "AlgNet-TK", "GML-Fixed"):
        model = models_to_train[name]
        model.eval()

        if isinstance(model, AlgebraicNet2) and isinstance(model.alg_layer, AlgebraicProjection2):
            G = model.alg_layer.get_generators()
            try:
                from killing_utils import analyze_generators
                metrics = analyze_generators(G)
                print_algebra_report(metrics, name)
            except Exception as exc:
                print(f"  [{name}] analysis failed: {exc}")

            # Theta orthogonality: ||θᵀθ - I||_F
            if not model.alg_layer.fixed:
                theta     = model.alg_layer.theta.detach().cpu()
                orth_err  = (theta.T @ theta - torch.eye(N_GEN)).norm().item()
                print(f"  θ orthogonality error ||θᵀθ - I||_F : {orth_err:.4f}")
                plot_theta(
                    theta, f"θ matrix — {name}",
                    BASE_NAME + f"_theta_{name.lower().replace('-','').replace(' ','')}.png"
                )

            # Generator heatmap
            slug = name.lower().replace("-", "").replace(" ", "")
            plot_generators(
                G, f"Learned generators — {name}",
                BASE_NAME + f"_generators_{slug}.png"
            )

    print("\nExperiment 4 complete.")
