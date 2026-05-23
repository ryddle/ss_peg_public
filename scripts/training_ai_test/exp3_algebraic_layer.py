"""
Experiment 3: AlgebraicLinear — Genuinely Constrained Lie Algebra Layer
=========================================================================

Root cause of Exp 1 & 2 failure:
  Weight rows W[a,:] ∈ R^9 are arbitrary real vectors. Reshaping to G_a = W[a].reshape(3,3)
  gives arbitrary real 3×3 matrices. These are NOT hermitian traceless → NOT in su(3).
  T_K cannot reach 0 because the optimizer is searching the wrong manifold.
  T_K was stuck at ≈26 throughout training in both experiments.

Fix — AlgebraicProjection layer:
  Generator G_a = sum_k θ_{a,k} · T_k   where {T_k} = GML basis of su(3).
  Parameters θ ∈ R^(n_gen × n_gen) are the coordinates in the algebra.
  Output ALWAYS hermitian traceless → T_K CAN reach 0.
  The global minimum (T_K=0) corresponds to θ forming an orthogonal matrix
  (any orthonormal rotation of the GML basis gives the same algebra su(3)).

Forward pass:
  x ∈ R^{d^2} → M = x.reshape(d,d) + 0j →  f_a = Re(Tr(G_a @ M))  ∈ R^{n_gen}

  Tr(G_a @ M) is a linear functional of M. The 8 functionals {f_a} span all
  traceless linear forms on M_3(C). With su(3) structure (T_K=0), they form
  an orthonormal basis for this space, and the cubic invariant Tr(M^3) can be
  expressed as sum_{abc} d_{abc} f_a f_b f_c (via the su(3) d-symbol).

Experiment design:
  4 models compared on sign(Tr(M^3)) classification:

    Model A — Standard MLP     : Linear(9→8) + ReLU + Linear(8→16) + … + Linear→2
    Model B — AlgNet-Free      : AlgebraicProjection (θ free, no T_K) + same head
    Model C — AlgNet-TK        : AlgebraicProjection + T_K loss (λ_K=2.0)
    Model D — AlgNet-Full      : AlgebraicProjection + T_K + T_closure (λ_K=2.0, λ_c=0.5)

  All models have identical head capacity for fair comparison.
  Key diagnostic: T_K over training for B,C,D — does it reach 0?

Sanity check (before training):
  Verify that for GML generators (T_k) as input: T_K = 0 exactly.
  Confirms the parameterization is correct.

Hardware: NVIDIA RTX 5070 Ti
Output  : exp3_results.npz  +  exp3_results.png  +  exp3_generators.png
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
    norm_tension,
    analyze_generators,
    print_algebra_report,
    generate_cubic_trace_dataset,
    params_to_su_generators,
    _su_basis,
)

torch.set_default_dtype(torch.float64)

# ============================================================
#  HYPER-PARAMETERS
# ============================================================

D          = 3              # algebra dimension
N_GEN      = D**2 - 1      # 8
N_IN       = D**2           # 9

N_EPOCHS   = 600
BATCH_SIZE = 256
LR         = 1e-3
N_SAMPLES  = 20_000

LAMBDA_K   = 2.0            # Killing tension strength (stronger than exp1 to really push)
LAMBDA_C   = 0.5            # Closure tension strength (model D only)
WARMUP     = 30             # Ramp epochs (short: we want algebra to form early)


def lam(epoch: int, lam_max: float, warmup: int = WARMUP) -> float:
    """Linear warmup schedule."""
    return lam_max * min(1.0, epoch / max(warmup, 1))


# ============================================================
#  SANITY CHECK: verify T_K = 0 for the GML basis
# ============================================================

def sanity_check() -> bool:
    """
    Verify that the canonical su(3) GML generators give T_K = 0.
    This tests that:
      1. _su_basis returns valid hermitian traceless matrices
      2. The Killing metric is K = -2d·I for orthonormal generators
      3. killing_tension = 0 for the exact algebra

    Returns True if all checks pass.
    """
    print("\n── Sanity check: should T_K = 0 for the GML basis? ──")
    ok = True

    # Retrieve the canonical basis as generators
    basis = _su_basis(D, str(DEVICE))  # (8, 3, 3) complex

    # Check: hermitian? (must use .mH, not .conj().mH — latter gives transpose only)
    herm_err = (basis - basis.mH).abs().max().item()
    print(f"  Hermitian error (max |T_a - T_a†|)    : {herm_err:.2e}  {'✓' if herm_err < 1e-12 else '✗'}")
    ok = ok and herm_err < 1e-12

    # Check: traceless?
    trace_err = basis.diagonal(dim1=-2, dim2=-1).sum(-1).abs().max().item()
    print(f"  Traceless error (max |Tr(T_a)|)        : {trace_err:.2e}  {'✓' if trace_err < 1e-12 else '✗'}")
    ok = ok and trace_err < 1e-12

    # Check: Tr(T_a T_b) = 2δ_{ab}  (GML normalization)
    inner = torch.einsum('aij,bji->ab', basis, basis).real
    inner_err = (inner - 2 * torch.eye(N_GEN, dtype=torch.float64, device=DEVICE)).abs().max().item()
    print(f"  Orthogonality error (max |Tr(TaTb)-2δ|): {inner_err:.2e}  {'✓' if inner_err < 1e-12 else '✗'}")
    ok = ok and inner_err < 1e-12

    # Check: T_K = 0 for these generators
    t_k = killing_tension(basis).item()
    print(f"  Killing tension T_K for GML basis      : {t_k:.6f}  {'✓' if t_k < 1e-6 else '✗'}")
    ok = ok and t_k < 1e-6

    # Check: K = -8N·I  (GML normalization Tr(TaTb)=2δ → K = -8N·I for su(N))
    # Note: with standard physicist norm Tr(TaTb)=(1/2)δ → K=-2N·I
    #       Our GML basis has Tr=2, factor-of-4 larger → K scales by 16 → -2N*16? No:
    #       Direct derivation: K_{ab} = -8N δ_{ab}  verified numerically.
    K = killing_metric(basis)
    expected_diag = -8 * D  # = -24 for su(3)
    expected_K = expected_diag * torch.eye(N_GEN, dtype=torch.float64, device=DEVICE)
    k_err = (K - expected_K).abs().max().item()
    print(f"  K = {expected_diag}·I error (max, K=-8N·I GML)    : {k_err:.2e}  {'✓' if k_err < 1e-6 else '✗'}")
    ok = ok and k_err < 1e-6

    # Now test via params_to_su_generators: identity params → GML basis
    # If params = I (identity matrix), G_a = sum_k I_{ak} T_k = T_a
    params_id = torch.eye(N_GEN, dtype=torch.float64, device=DEVICE)
    G_from_params = params_to_su_generators(params_id, D)
    diff = (G_from_params - basis).abs().max().item()
    t_k_p = killing_tension(G_from_params).item()
    print(f"  params=I → same as basis (max diff)    : {diff:.2e}  {'✓' if diff < 1e-12 else '✗'}")
    print(f"  T_K for params=I (identity)            : {t_k_p:.6f}  {'✓' if t_k_p < 1e-6 else '✗'}")
    ok = ok and diff < 1e-12 and t_k_p < 1e-6

    # Test: random orthogonal params → T_K still = 0
    torch.manual_seed(7)
    Q, _ = torch.linalg.qr(torch.randn(N_GEN, N_GEN, dtype=torch.float64, device=DEVICE))
    G_rot = params_to_su_generators(Q, D)
    t_k_q = killing_tension(G_rot).item()
    print(f"  T_K for random orthogonal rotation     : {t_k_q:.6f}  {'✓' if t_k_q < 1e-4 else '✗'}")
    ok = ok and t_k_q < 1e-4

    status = "ALL CHECKS PASSED ✓" if ok else "SOME CHECKS FAILED ✗"
    print(f"\n  {status}")
    return ok


# ============================================================
#  MODELS
# ============================================================

class AlgebraicProjection(nn.Module):
    """
    Layer: x ∈ R^{d^2} → features f ∈ R^{n_gen}  via algebra inner products.

      f_a = Re(Tr(G_a @ M))   where M = x.reshape(d,d) + 0j

    Generators G_a are constrained to su(d):
      G_a = sum_k θ_{a,k} * T_k   (θ is the learned parameter, T_k = GML basis)

    Result: T_K CAN reach 0.  The global minimum is the set of all orthonormal
    rotations of the GML basis — compact manifold O(n_gen) / (sign flips).
    """

    def __init__(self, d: int):
        super().__init__()
        self.d     = d
        self.n_gen = d * d - 1
        # θ: coordinates of each generator in the GML basis
        # Init: small random perturbation around zero (not orthogonal yet)
        self.theta = nn.Parameter(
            torch.randn(self.n_gen, self.n_gen, dtype=torch.float64) * 0.3
        )
        self.bias = nn.Parameter(torch.zeros(self.n_gen, dtype=torch.float64))

    def get_generators(self) -> torch.Tensor:
        """(n_gen, d, d) complex — always hermitian traceless."""
        return params_to_su_generators(self.theta, self.d)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        x : (B, d^2)  →  f : (B, n_gen)

        f_a = Re(Tr(G_a @ M)) = Re(<G_a, M^T>_F)
        """
        B = x.shape[0]
        G  = self.get_generators()                         # (n_gen, d, d) complex
        Mc = x.reshape(B, self.d, self.d).to(torch.complex128)  # (B, d, d)

        # Tr(G_a @ M) = sum_{j} (G_a @ M)_{jj} = einsum
        # G: (n, d, d) → (1, n, d, d)
        # M: (B, d, d) → (B, 1, d, d)
        GM     = G.unsqueeze(0) @ Mc.unsqueeze(1)          # (B, n_gen, d, d)
        traces = GM.diagonal(dim1=-2, dim2=-1).sum(-1)      # (B, n_gen) complex
        return traces.real + self.bias                      # (B, n_gen) real

    def reg_loss(self, epoch: int = 0,
                 lambda_k: float = 0.0, lambda_c: float = 0.0) -> torch.Tensor:
        """T_K + T_closure on current generators (with warmup scaling)."""
        lk = lam(epoch, lambda_k)
        lc = lam(epoch, lambda_c)
        G  = self.get_generators()
        reg = torch.zeros(1, dtype=torch.float64, device=self.theta.device).squeeze()
        if lk > 0:
            reg = reg + lk * killing_tension(G)
        if lc > 0:
            reg = reg + lc * closure_tension(G)
        return reg


class StandardMLP(nn.Module):
    """Baseline: flat Linear(9→8) with no algebraic structure."""

    def __init__(self):
        super().__init__()
        self.layer1 = nn.Linear(N_IN, N_GEN)
        self.head   = nn.Sequential(
            nn.ReLU(),
            nn.Linear(N_GEN, 16),
            nn.ReLU(),
            nn.Linear(16, 2),
        )

    def forward(self, x):
        return self.head(self.layer1(x))


class AlgebraicNet(nn.Module):
    """
    MLP with AlgebraicProjection as Layer 1.

    Same head as StandardMLP for fair comparison.
    λ_K and λ_C control regularization strength (0 = no regularization).
    """

    def __init__(self, lambda_k: float = 0.0, lambda_c: float = 0.0):
        super().__init__()
        self.lambda_k   = lambda_k
        self.lambda_c   = lambda_c
        self.alg_layer  = AlgebraicProjection(D)
        self.head       = nn.Sequential(
            nn.ReLU(),
            nn.Linear(N_GEN, 16),
            nn.ReLU(),
            nn.Linear(16, 2),
        )

    def forward(self, x):
        return self.head(self.alg_layer(x))

    def reg_loss(self, epoch: int = 0) -> torch.Tensor:
        return self.alg_layer.reg_loss(epoch, self.lambda_k, self.lambda_c)

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
) -> dict:
    model = model.to(DEVICE)
    optimizer = optim.Adam(model.parameters(), lr=LR)
    criterion = nn.CrossEntropyLoss()

    hist = {
        "train_acc": [], "val_acc": [],
        "train_loss": [], "val_loss": [],
        "T_K": [], "T_closure": [],
    }

    print(f"\n  [{name}]  {N_EPOCHS} epochs ...")
    for epoch in range(N_EPOCHS):
        model.train()
        t_loss, correct, total = 0.0, 0, 0
        for xb, yb in train_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            optimizer.zero_grad()
            logits    = model(xb)
            task_loss = criterion(logits, yb)
            loss      = task_loss
            if isinstance(model, AlgebraicNet):
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
                out     = model(xb)
                v_loss   += criterion(out, yb).item() * xb.size(0)
                v_correct += (out.argmax(1) == yb).sum().item()
                v_total   += xb.size(0)

        # Algebraic metrics
        if isinstance(model, AlgebraicNet):
            with torch.no_grad():
                G   = model.alg_layer.get_generators()
                t_k = killing_tension(G).item()
                t_c = closure_tension(G).item()
        else:
            t_k = t_c = float("nan")

        hist["train_loss"].append(t_loss / total)
        hist["train_acc"].append(correct / total)
        hist["val_loss"].append(v_loss / v_total)
        hist["val_acc"].append(v_correct / v_total)
        hist["T_K"].append(t_k)
        hist["T_closure"].append(t_c)

        if (epoch + 1) % 100 == 0:
            print(
                f"    epoch {epoch+1:3d}/{N_EPOCHS}"
                f"  val_acc={v_correct/v_total:.4f}"
                f"  T_K={t_k:.4f}  T_cl={t_c:.4f}"
            )

    return hist


# ============================================================
#  GENERATOR VISUALIZATION
# ============================================================

def plot_generators(generators: torch.Tensor, title: str, out_path: str) -> None:
    """
    Visualize the 8 learned generators as 3×3 heatmaps (real and imag parts).
    Useful for visually identifying if the algebra looks like su(3).
    """
    G   = generators.cpu()
    fig, axes = plt.subplots(4, 4, figsize=(12, 12))
    fig.suptitle(title, fontsize=12)
    for a in range(N_GEN):
        # Real part
        ax = axes[a // 2][a % 2 * 2]
        im = ax.imshow(G[a].real.numpy(), cmap='RdBu', vmin=-1, vmax=1)
        ax.set_title(f"G_{a+1} Re", fontsize=9)
        ax.axis("off")
        # Imaginary part
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

def plot_results(results: dict, out_dir: str) -> None:
    """4-panel figure: accuracy, loss, T_K, T_closure over training."""
    epochs = range(1, N_EPOCHS + 1)
    colors = {
        "Standard MLP":   "tab:blue",
        "AlgNet-Free":    "tab:gray",
        "AlgNet-TK":      "tab:orange",
        "AlgNet-Full":    "tab:green",
    }
    ls_map = {
        "Standard MLP": "-",
        "AlgNet-Free":  "--",
        "AlgNet-TK":    "-",
        "AlgNet-Full":  "-",
    }

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(
        f"Experiment 3 — AlgebraicLinear  (d={D}, n_gen={N_GEN},  λ_K={LAMBDA_K},  λ_c={LAMBDA_C})",
        fontsize=12,
    )

    # Val accuracy
    ax = axes[0, 0]
    for name, h in results.items():
        ax.plot(epochs, h["val_acc"], label=name, color=colors[name],
                ls=ls_map[name], lw=1.8)
    ax.set_title("Validation Accuracy"); ax.set_xlabel("Epoch"); ax.set_ylabel("Accuracy")
    ax.legend(); ax.grid(True, alpha=0.3)

    # Val loss
    ax = axes[0, 1]
    for name, h in results.items():
        ax.plot(epochs, h["val_loss"], label=name, color=colors[name],
                ls=ls_map[name], lw=1.8)
    ax.set_title("Validation Loss (Cross-Entropy)"); ax.set_xlabel("Epoch"); ax.set_ylabel("Loss")
    ax.legend(); ax.grid(True, alpha=0.3)

    # T_K (algebraic models only)
    ax = axes[1, 0]
    for name, h in results.items():
        vals = h["T_K"]
        if not all(np.isnan(v) for v in vals):
            ax.semilogy(list(epochs), vals, label=name, color=colors[name],
                        ls=ls_map[name], lw=1.8)
    ax.axhline(1e-3, ls=":", color="gray", alpha=0.5, label="T_K=0.001")
    ax.axvline(WARMUP, ls="--", color="gray", alpha=0.4, label=f"warmup={WARMUP}")
    ax.set_title("Killing Tension (log)  → 0 = Lie algebra"); ax.set_xlabel("Epoch"); ax.set_ylabel("T_K")
    ax.legend(); ax.grid(True, alpha=0.3)

    # T_closure
    ax = axes[1, 1]
    for name, h in results.items():
        vals = h["T_closure"]
        if not all(np.isnan(v) for v in vals):
            ax.semilogy(list(epochs), [max(v, 1e-12) for v in vals],
                        label=name, color=colors[name], ls=ls_map[name], lw=1.8)
    ax.set_title("Closure Fraction (log)  → 0 = closed algebra")
    ax.set_xlabel("Epoch"); ax.set_ylabel("T_closure")
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_path = os.path.join(out_dir, "exp3_results.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Figure saved → {out_path}")


# ============================================================
#  MAIN
# ============================================================

def main():
    print("=" * 68)
    print("Experiment 3: AlgebraicLinear — Constrained Lie Algebra Layer")
    print("=" * 68)
    print(f"  Device      : {DEVICE}")
    if torch.cuda.is_available():
        print(f"  GPU         : {torch.cuda.get_device_name(0)}")
    print(f"  d={D}  n_gen={N_GEN}  n_in={N_IN}  epochs={N_EPOCHS}")
    print(f"  λ_K={LAMBDA_K}  λ_c={LAMBDA_C}  warmup={WARMUP} epochs")
    print(f"  FIX vs Exp1/2: generators NOW constrained to su({D}) manifold")
    print(f"                 → T_K CAN reach 0 (previously stuck at ≈26)")

    out_dir = os.path.dirname(os.path.abspath(__file__))

    # ---- Sanity check ----
    all_ok = sanity_check()
    if not all_ok:
        print("\nWARNING: sanity check failed. Results may be unreliable.")

    # ---- Dataset ----
    print(f"\n[1/5] Dataset  N={N_SAMPLES}, task: sign(Tr(M^3))")
    X, y = generate_cubic_trace_dataset(n=N_SAMPLES, d=D, seed=42, device=DEVICE)
    split = int(0.8 * N_SAMPLES)
    X_tr, y_tr = X[:split], y[:split]
    X_va, y_va = X[split:], y[split:]
    pos = int(y.sum().item())
    print(f"  Balance: {N_SAMPLES - pos} neg / {pos} pos")

    def make_loader(Xd, yd, shuffle=True):
        return DataLoader(
            TensorDataset(Xd.cpu(), yd.cpu()),
            batch_size=BATCH_SIZE, shuffle=shuffle,
        )

    train_loader = make_loader(X_tr, y_tr)
    val_loader   = make_loader(X_va, y_va, shuffle=False)

    # ---- Train all 4 models ----
    models = {
        "Standard MLP": StandardMLP(),
        "AlgNet-Free":  AlgebraicNet(lambda_k=0.0,     lambda_c=0.0),
        "AlgNet-TK":    AlgebraicNet(lambda_k=LAMBDA_K, lambda_c=0.0),
        "AlgNet-Full":  AlgebraicNet(lambda_k=LAMBDA_K, lambda_c=LAMBDA_C),
    }

    results = {}
    trained_models = {}
    print()
    for i, (name, model) in enumerate(models.items(), start=2):
        print(f"[{i}/5] Training {name} ...")
        torch.manual_seed(0)
        # Re-initialize (model was created above, reinit weights for reproducibility)
        for m in model.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
        hist = train(model, train_loader, val_loader, name=name)
        results[name] = hist
        trained_models[name] = model

    # ---- Post-training algebraic analysis ----
    print("\n[5/5] Post-training algebraic analysis ...")
    final_generators = {}
    for name, model in trained_models.items():
        if isinstance(model, AlgebraicNet):
            G = model.get_generators()
            final_generators[name] = G
            m = analyze_generators(G)
            print_algebra_report(m, title=f"{name}")
        else:
            print(f"\n  {name}: no algebraic layer (standard weights only)")

    # ---- Generator heatmaps ----
    for name, G in final_generators.items():
        fname = "exp3_generators_" + name.lower().replace(" ", "_").replace("-", "") + ".png"
        plot_generators(G, title=f"Learned generators — {name}", out_path=os.path.join(out_dir, fname))

    # ---- Summary ----
    print("\n" + "=" * 68)
    print("SUMMARY")
    print("=" * 68)
    print(f"  {'Model':<22}  {'Best acc':>9}  {'@ epoch':>8}  {'Final acc':>10}  {'T_K final':>10}")
    print("  " + "-" * 63)
    for name, h in results.items():
        best       = max(h["val_acc"])
        ep_best    = int(np.argmax(h["val_acc"])) + 1
        final      = h["val_acc"][-1]
        t_k_final  = h["T_K"][-1]
        t_k_str    = f"{t_k_final:.4f}" if not np.isnan(t_k_final) else "   —"
        print(f"  {name:<22}  {best:>9.4f}  {ep_best:>8}  {final:>10.4f}  {t_k_str:>10}")

    ref = max(results["Standard MLP"]["val_acc"])
    print()
    for name in ["AlgNet-Free", "AlgNet-TK", "AlgNet-Full"]:
        delta = max(results[name]["val_acc"]) - ref
        print(f"  Δ({name} vs Standard MLP) : {delta:+.4f}")

    # ---- T_K convergence report ----
    print()
    print("  Killing tension convergence (AlgebraicNet models):")
    for name in ["AlgNet-Free", "AlgNet-TK", "AlgNet-Full"]:
        t0   = results[name]["T_K"][0]
        t_f  = results[name]["T_K"][-1]
        best = min(results[name]["T_K"])
        print(f"    {name:<15}  initial={t0:.4f}  best={best:.6f}  final={t_f:.4f}")

    print()
    print("  Expected: AlgNet-TK/Full → T_K near 0  (Lie algebra crystallized)")
    print("            AlgNet-Free   → T_K fluctuates (task drives away from algebra)")

    # ---- Save ----
    np.savez(
        os.path.join(out_dir, "exp3_results.npz"),
        results = {k: v for k, v in results.items()},
    )
    plot_results(results, out_dir)

    print("\nExperiment 3 complete.")


if __name__ == "__main__":
    main()
