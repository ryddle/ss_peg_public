"""
Experiment 1: Killing-Regularized MLP vs Standard MLP
=======================================================

Hypothesis:
  Adding T_K (Killing tension) as regularization to a neural network layer whose
  weight rows are d×d matrices drives those matrices toward a compact Lie algebra.
  On a task with hidden algebraic symmetry, this leads to:
    (a) structured weight matrices (verifiable post-training)
    (b) potentially improved generalization / sample efficiency

Experiment design:
  Task   : Binary classification — sign(Tr(M^3)) > 0
           where M ∈ R^(3×3) is the input reshaped from x ∈ R^9.
           Tr(M^3) is a degree-3 polynomial conjugacy invariant,
           related to the third-order Casimir of sl(3).

  Key layer: Linear(9 → 8)  [8 = 3^2 - 1 = dim su(3)]
    Weight W ∈ R^(8 × 9): each row w_a ∈ R^9 → reshape to G_a = w_a.reshape(3,3)
    Killing regularization:
      T_K   = || K_normalized + I_8 ||^2   (drives K → -I: compact algebra)
    Full loss: L = CE(logits, y) + λ_K * T_K

  Models compared:
    (A) StandardMLP  — cross-entropy only
    (B) KillingMLP   — cross-entropy + λ_K * T_K on layer 1
    (C) KillingMLP+C — cross-entropy + λ_K * T_K + λ_c * T_closure on layer 1

  Metrics:
    - Validation accuracy over training (does T_K help?)
    - T_K value over training           (does algebra crystallize?)
    - Killing eigenvalues at end        (are they proportional to -1?)
    - Closure fraction at end           (is the algebra self-consistent?)

Hardware: NVIDIA RTX 5070 Ti — GPU-accelerated via PyTorch CUDA
Output  : exp1_results.npz  +  exp1_results.png
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
    killing_regularization,
    analyze_generators,
    print_algebra_report,
    extract_generators,
    generate_cubic_trace_dataset,
)

torch.set_default_dtype(torch.float64)

# ============================================================
#  HYPER-PARAMETERS
# ============================================================

D          = 3          # algebra representation dimension
N_GEN      = D**2 - 1  # = 8 = dim su(3)
N_IN       = D**2      # = 9 = input dimension (rows of W reshape to d×d)
N_SAMPLES  = 20_000
N_EPOCHS   = 400
BATCH_SIZE = 256
LR         = 1e-3

# Lambda schedule: ramp from 0 → λ_max over WARMUP epochs, then hold
LAMBDA_K_MAX  = 0.05    # Killing regularization strength
LAMBDA_C_MAX  = 0.0     # Closure regularization (off in exp1, used in KillingMLP+C)
WARMUP_EPOCHS = 50      # Ramp duration (let task loss stabilize first)


def lambda_schedule(epoch: int, lam_max: float, warmup: int) -> float:
    """Linear ramp-up from 0 to lam_max over warmup epochs, then constant."""
    if epoch < warmup:
        return lam_max * epoch / warmup
    return lam_max


# ============================================================
#  MODELS
# ============================================================

class StandardMLP(nn.Module):
    """Baseline MLP — no algebraic regularization."""

    def __init__(self):
        super().__init__()
        self.layer1 = nn.Linear(N_IN, N_GEN)
        self.rest   = nn.Sequential(
            nn.ReLU(),
            nn.Linear(N_GEN, 16),
            nn.ReLU(),
            nn.Linear(16, 2),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.rest(self.layer1(x))


class KillingMLP(nn.Module):
    """
    MLP with Killing regularization on the first layer.

    Layer 1: Linear(d^2, d^2-1)
      Weight W ∈ R^(n_gen × d^2): rows w_a → G_a = w_a.reshape(d, d)
      T_K({ G_a }) is computed and returned via reg_loss() for the training loop.
    """

    def __init__(self, lambda_killing: float = LAMBDA_K_MAX,
                 lambda_closure: float = 0.0):
        super().__init__()
        self.d               = D
        self.lambda_killing  = lambda_killing
        self.lambda_closure  = lambda_closure

        self.layer1 = nn.Linear(N_IN, N_GEN)
        self.rest   = nn.Sequential(
            nn.ReLU(),
            nn.Linear(N_GEN, 16),
            nn.ReLU(),
            nn.Linear(16, 2),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.rest(self.layer1(x))

    def reg_loss(self, epoch: int = 0) -> torch.Tensor:
        """Compute Killing (+ closure) regularization for current weights."""
        lk = lambda_schedule(epoch, self.lambda_killing, WARMUP_EPOCHS)
        lc = lambda_schedule(epoch, self.lambda_closure, WARMUP_EPOCHS)
        return killing_regularization(
            self.layer1.weight,
            d               = self.d,
            lambda_killing  = lk,
            lambda_closure  = lc,
            lambda_norm     = 0.0,    # norm stabilizer off for task training
        )

    def get_generators(self) -> torch.Tensor:
        return extract_generators(self.layer1.weight, self.d).detach()


# ============================================================
#  TRAINING LOOP
# ============================================================

def train(
    model:       nn.Module,
    train_loader: DataLoader,
    val_loader:   DataLoader,
    name:         str = "Model",
) -> dict:
    """Train model and record history every epoch."""
    model = model.to(DEVICE)
    optimizer = optim.Adam(model.parameters(), lr=LR)
    criterion = nn.CrossEntropyLoss()

    hist = {
        "train_acc": [], "val_acc": [],
        "train_loss": [], "val_loss": [],
        "T_K": [], "T_closure": [],
    }

    print(f"\n  Training {name} for {N_EPOCHS} epochs ...")
    for epoch in range(N_EPOCHS):
        # ---- train ----
        model.train()
        t_loss, correct, total = 0.0, 0, 0
        for xb, yb in train_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            optimizer.zero_grad()
            logits    = model(xb)
            task_loss = criterion(logits, yb)
            loss      = task_loss
            if isinstance(model, KillingMLP):
                loss = loss + model.reg_loss(epoch)
            loss.backward()
            optimizer.step()
            t_loss  += task_loss.item() * xb.size(0)
            correct += (logits.argmax(1) == yb).sum().item()
            total   += xb.size(0)

        # ---- validate ----
        model.eval()
        v_loss, v_correct, v_total = 0.0, 0, 0
        with torch.no_grad():
            for xb, yb in val_loader:
                xb, yb = xb.to(DEVICE), yb.to(DEVICE)
                logits    = model(xb)
                v_loss   += criterion(logits, yb).item() * xb.size(0)
                v_correct += (logits.argmax(1) == yb).sum().item()
                v_total   += xb.size(0)

        # ---- algebraic metrics ----
        if isinstance(model, KillingMLP):
            with torch.no_grad():
                G = extract_generators(model.layer1.weight.double(), D)
                from killing_utils import killing_tension, closure_tension
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

        if (epoch + 1) % 50 == 0:
            print(
                f"    [{name}]  epoch {epoch+1:3d}/{N_EPOCHS}"
                f"  val_acc={v_correct/v_total:.3f}"
                f"  T_K={t_k:.4f}  T_cl={t_c:.4f}"
            )

    return hist


# ============================================================
#  PLOTTING
# ============================================================

def plot_results(results: dict, out_path: str) -> None:
    """4-panel figure: accuracy, loss, T_K, T_closure."""
    epochs = range(1, N_EPOCHS + 1)
    colors = {"Standard MLP": "tab:blue", "Killing MLP": "tab:orange",
              "Killing MLP+C": "tab:green"}

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(
        f"Experiment 1 — Killing-Regularized MLP  (d={D}, n_gen={N_GEN}, λ_K={LAMBDA_K_MAX})",
        fontsize=13,
    )

    # Validation accuracy
    ax = axes[0, 0]
    for name, hist in results.items():
        ax.plot(epochs, hist["val_acc"], label=name, color=colors[name], lw=1.8)
    ax.set_xlabel("Epoch"); ax.set_ylabel("Validation Accuracy")
    ax.set_title("Validation Accuracy"); ax.legend(); ax.grid(True, alpha=0.3)

    # Validation loss
    ax = axes[0, 1]
    for name, hist in results.items():
        ax.plot(epochs, hist["val_loss"], label=name, color=colors[name], lw=1.8)
    ax.set_xlabel("Epoch"); ax.set_ylabel("Cross-Entropy Loss")
    ax.set_title("Validation Loss"); ax.legend(); ax.grid(True, alpha=0.3)

    # T_K over training (Killing models only)
    ax = axes[1, 0]
    for name, hist in results.items():
        t_k = hist["T_K"]
        if not all(np.isnan(v) for v in t_k):
            ax.plot(epochs, t_k, label=name, color=colors[name], lw=1.8)
    ax.axvline(WARMUP_EPOCHS, ls="--", color="gray", alpha=0.5, label="warmup end")
    ax.set_xlabel("Epoch"); ax.set_ylabel("T_K")
    ax.set_title("Killing Tension in Layer 1  (→ 0 = Lie algebra)"); ax.legend()
    ax.grid(True, alpha=0.3)
    try:
        ax.set_yscale("log")
    except Exception:
        pass

    # T_closure over training
    ax = axes[1, 1]
    for name, hist in results.items():
        t_c = hist["T_closure"]
        if not all(np.isnan(v) for v in t_c):
            ax.plot(epochs, t_c, label=name, color=colors[name], lw=1.8)
    ax.set_xlabel("Epoch"); ax.set_ylabel("T_closure")
    ax.set_title("Closure Fraction in Layer 1  (→ 0 = closed algebra)")
    ax.legend(); ax.grid(True, alpha=0.3)
    try:
        ax.set_yscale("log")
    except Exception:
        pass

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Figure saved → {out_path}")


# ============================================================
#  MAIN
# ============================================================

def main():
    print("=" * 60)
    print("Experiment 1: Killing-Regularized MLP vs Standard MLP")
    print("=" * 60)
    print(f"  Device  : {DEVICE}")
    if torch.cuda.is_available():
        print(f"  GPU     : {torch.cuda.get_device_name(0)}")
    print(f"  d       : {D}   →  n_gen = {N_GEN}  (=dim su({D}))")
    print(f"  n_in    : {N_IN}")
    print(f"  λ_K max : {LAMBDA_K_MAX}   (warmup = {WARMUP_EPOCHS} epochs)")

    # ---- 1. Dataset ----
    print(f"\n[1/5] Generating dataset  (N={N_SAMPLES}, task: sign(Tr(M^3)))...")
    X, y = generate_cubic_trace_dataset(n=N_SAMPLES, d=D, seed=42, device=DEVICE)
    split = int(0.8 * N_SAMPLES)
    X_tr, y_tr = X[:split], y[:split]
    X_va, y_va = X[split:], y[split:]
    counts = y.cpu().numpy()
    pos    = int(counts.sum())
    print(f"  Balance: {N_SAMPLES - pos} neg / {pos} pos  ({100*pos/N_SAMPLES:.1f}% positive)")

    def make_loader(Xd, yd, shuffle=True):
        return DataLoader(
            TensorDataset(Xd.cpu(), yd.cpu()),
            batch_size=BATCH_SIZE, shuffle=shuffle,
        )

    train_loader = make_loader(X_tr, y_tr, shuffle=True)
    val_loader   = make_loader(X_va, y_va, shuffle=False)

    # ---- 2. Standard MLP ----
    print("\n[2/5] Training Standard MLP ...")
    torch.manual_seed(0)
    std_mlp = StandardMLP()
    hist_std = train(std_mlp, train_loader, val_loader, name="Standard MLP")

    # ---- 3. Killing MLP (T_K only) ----
    print("\n[3/5] Training Killing MLP  (λ_K=%g) ..." % LAMBDA_K_MAX)
    torch.manual_seed(0)
    kil_mlp = KillingMLP(lambda_killing=LAMBDA_K_MAX, lambda_closure=0.0)
    hist_kil = train(kil_mlp, train_loader, val_loader, name="Killing MLP")

    # ---- 4. Killing MLP + Closure ----
    print("\n[4/5] Training Killing MLP+C  (λ_K=%g, λ_c=0.02) ..." % LAMBDA_K_MAX)
    torch.manual_seed(0)
    kil_c_mlp = KillingMLP(lambda_killing=LAMBDA_K_MAX, lambda_closure=0.02)
    hist_kil_c = train(kil_c_mlp, train_loader, val_loader, name="Killing MLP+C")

    # ---- 5. Analysis ----
    print("\n[5/5] Post-training algebraic analysis ...")

    for name, model in [("Killing MLP", kil_mlp), ("Killing MLP+C", kil_c_mlp)]:
        G = extract_generators(model.layer1.weight.double(), D)
        m = analyze_generators(G)
        print_algebra_report(m, title=f"Layer 1 — {name}")

    # Standard MLP: just measure its weight structure for comparison
    G_std = extract_generators(std_mlp.layer1.weight.double(), D)
    m_std = analyze_generators(G_std)
    print_algebra_report(m_std, title="Layer 1 — Standard MLP  (no regularization)")

    # ---- Summary ----
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for name, h in [("Standard MLP", hist_std), ("Killing MLP", hist_kil),
                    ("Killing MLP+C", hist_kil_c)]:
        best = max(h["val_acc"])
        final = h["val_acc"][-1]
        print(f"  {name:20s}  best={best:.4f}  final={final:.4f}")
    print()
    b_std = max(hist_std["val_acc"])
    b_kil = max(hist_kil["val_acc"])
    b_klc = max(hist_kil_c["val_acc"])
    print(f"  Δ(Killing vs Standard)   : {b_kil - b_std:+.4f}")
    print(f"  Δ(Killing+C vs Standard) : {b_klc - b_std:+.4f}")

    # ---- Save ----
    out_dir = os.path.dirname(os.path.abspath(__file__))
    np.savez(
        os.path.join(out_dir, "exp1_results.npz"),
        hist_standard   = hist_std,
        hist_killing    = hist_kil,
        hist_killing_c  = hist_kil_c,
    )

    results = {
        "Standard MLP":   hist_std,
        "Killing MLP":    hist_kil,
        "Killing MLP+C":  hist_kil_c,
    }
    plot_results(results, os.path.join(out_dir, "exp1_results.png"))

    print("\nExperiment 1 complete.")


if __name__ == "__main__":
    main()
