"""
Experiment 2: Killing Pretraining — Symmetry as Initialization
===============================================================

Two-phase experiment. Tests whether using T_K as a *self-supervised pretraining
objective* gives a better initialization for task learning.

Motivation:
  In the SS-PEG verification, the optimizer finds exact Lie algebra generators
  by minimizing T_K alone — without any task-specific data.  If we pretrain a
  network layer this way (pure T_K, no labels), do the resulting weights form a
  'better prior' that accelerates learning on a downstream task?

Phase 1 — Killing Pretraining (unsupervised, no labels):
  Minimize only T_K({G_a}) on layer 1 of the network.
  No data needed — pure geometric optimization.
  Expected result: layer 1 weights converge toward su(d) generators (T_K → 0).

Phase 2 — Task Fine-tuning (supervised):
  Freeze or unfreeze layer 1, then train on the classification task.

Baselines:
  (A) Random init → full training from scratch  [standard baseline]
  (B) Killing pretrain → fine-tune all layers   [main hypothesis]
  (C) Killing pretrain → fine-tune rest only     [frozen algebra layer]

Task:
  Same as Exp 1: sign(Tr(M^3)) classification on R^9 input.
  Architecture: Linear(9→8) + ReLU + Linear(8→16) + ReLU + Linear(16→2)

Key question:
  Does starting with algebraically structured weights (T_K ≈ 0) lead to
  faster convergence and/or higher final accuracy?

Hardware: NVIDIA RTX 5070 Ti
Output  : exp2_results.npz  +  exp2_results.png
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
import copy

from killing_utils import (
    DEVICE,
    killing_tension,
    closure_tension,
    norm_tension,
    extract_generators,
    analyze_generators,
    print_algebra_report,
    generate_cubic_trace_dataset,
)

torch.set_default_dtype(torch.float64)

# ============================================================
#  HYPER-PARAMETERS
# ============================================================

D           = 3             # algebra dimension
N_GEN       = D**2 - 1     # 8
N_IN        = D**2          # 9

# Pretraining
PRETRAIN_STEPS  = 2_000     # pure T_K optimization steps (no data)
PRETRAIN_LR     = 5e-3
LAMBDA_K_PRE    = 5.0       # aggressive: same scale as verification/* tests
LAMBDA_N_PRE    = 1.0       # norm stabilizer during pretraining

# Fine-tuning / full training
FINETUNE_EPOCHS = 300
FULL_EPOCHS     = 300       # baseline trains for same total budget
BATCH_SIZE      = 256
LR              = 1e-3

N_SAMPLES       = 20_000


# ============================================================
#  SHARED NETWORK DEFINITION
# ============================================================

class AlgebraNet(nn.Module):
    """
    Simple MLP with an identifiable 'algebraic layer' (layer 1).
    n_in = d^2,  n_hidden = d^2 - 1.
    """

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

    def generators(self) -> torch.Tensor:
        return extract_generators(self.layer1.weight.double(), D).detach()


# ============================================================
#  PHASE 1: KILLING PRETRAINING (pure geometry, no data)
# ============================================================

def killing_pretrain(model: AlgebraNet) -> dict:
    """
    Minimize T_K + T_norm on layer 1 weights only.
    No data. Direct algebraic optimization — mirrors the SS-PEG verification tests.

    Returns history dict with T_K and T_norm per step.
    """
    # Only optimize layer1 parameters
    optimizer = optim.Adam(model.layer1.parameters(), lr=PRETRAIN_LR)

    hist = {"T_K": [], "T_norm": [], "step": []}

    model.train()
    print(f"\n  [Pretraining] {PRETRAIN_STEPS} steps, λ_K={LAMBDA_K_PRE}, λ_n={LAMBDA_N_PRE}")

    for step in range(PRETRAIN_STEPS):
        optimizer.zero_grad()

        W = model.layer1.weight.double()
        G = extract_generators(W, D)

        t_k   = killing_tension(G)
        t_n   = norm_tension(G)
        loss  = LAMBDA_K_PRE * t_k + LAMBDA_N_PRE * t_n

        loss.backward()
        optimizer.step()

        hist["T_K"].append(t_k.item())
        hist["T_norm"].append(t_n.item())
        hist["step"].append(step)

        if (step + 1) % 500 == 0:
            print(
                f"    step {step+1:5d}/{PRETRAIN_STEPS}"
                f"  T_K={t_k.item():.4f}  T_norm={t_n.item():.4f}"
            )

    print(f"  Pretraining done.  Final T_K = {hist['T_K'][-1]:.6f}")
    return hist


# ============================================================
#  PHASE 2: TASK FINE-TUNING / FULL TRAINING
# ============================================================

def train_task(
    model:        AlgebraNet,
    train_loader: DataLoader,
    val_loader:   DataLoader,
    n_epochs:     int   = FINETUNE_EPOCHS,
    freeze_l1:    bool  = False,
    name:         str   = "Model",
) -> dict:
    """
    Train on the classification task.
    If freeze_l1=True, layer 1 weights are frozen (only rest is updated).
    """
    if freeze_l1:
        for p in model.layer1.parameters():
            p.requires_grad_(False)
    else:
        for p in model.layer1.parameters():
            p.requires_grad_(True)

    model = model.to(DEVICE)
    optimizer = optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=LR,
    )
    criterion = nn.CrossEntropyLoss()

    hist = {
        "train_acc": [], "val_acc": [],
        "train_loss": [], "val_loss": [],
        "T_K": [], "T_closure": [],
    }

    print(f"\n  Training {name} for {n_epochs} epochs "
          f"{'(L1 frozen)' if freeze_l1 else '(all layers)'}...")

    for epoch in range(n_epochs):
        model.train()
        t_loss, correct, total = 0.0, 0, 0
        for xb, yb in train_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            optimizer.zero_grad()
            loss = criterion(model(xb), yb)
            loss.backward()
            optimizer.step()
            t_loss  += loss.item() * xb.size(0)
            correct += (model(xb).argmax(1) == yb).sum().item()
            total   += xb.size(0)

        model.eval()
        v_loss, v_correct, v_total = 0.0, 0, 0
        with torch.no_grad():
            for xb, yb in val_loader:
                xb, yb = xb.to(DEVICE), yb.to(DEVICE)
                out = model(xb)
                v_loss    += criterion(out, yb).item() * xb.size(0)
                v_correct += (out.argmax(1) == yb).sum().item()
                v_total   += xb.size(0)

        # Algebraic metrics on layer 1
        with torch.no_grad():
            G   = extract_generators(model.layer1.weight.double(), D)
            t_k = killing_tension(G).item()
            t_c = closure_tension(G).item()

        hist["train_loss"].append(t_loss / total)
        hist["train_acc"].append(correct / total)
        hist["val_loss"].append(v_loss / v_total)
        hist["val_acc"].append(v_correct / v_total)
        hist["T_K"].append(t_k)
        hist["T_closure"].append(t_c)

        if (epoch + 1) % 50 == 0:
            print(
                f"    [{name}]  epoch {epoch+1:3d}/{n_epochs}"
                f"  val_acc={v_correct/v_total:.3f}"
                f"  T_K={t_k:.4f}  T_cl={t_c:.4f}"
            )

    return hist


# ============================================================
#  PLOTTING
# ============================================================

def plot_pretraining(pre_hist: dict, out_path: str) -> None:
    """Plot T_K convergence during the pretraining phase."""
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.semilogy(pre_hist["step"], pre_hist["T_K"], color="tab:red", lw=1.8)
    ax.set_xlabel("Optimization step"); ax.set_ylabel("T_K  (log scale)")
    ax.set_title(f"Phase 1 — Killing Pretraining  (d={D}, n_gen={N_GEN})")
    ax.grid(True, alpha=0.3)
    ax.axhline(0.1, ls="--", color="gray", alpha=0.5, label="T_K = 0.1 threshold")
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Figure saved → {out_path}")


def plot_finetuning(results: dict, out_path: str) -> None:
    """3-panel figure: accuracy, T_K, T_closure over fine-tuning epochs."""
    colors = {
        "Random init":           "tab:blue",
        "Killing pretrain (all)":"tab:orange",
        "Killing pretrain (L1 frozen)": "tab:green",
    }

    epochs = range(1, FINETUNE_EPOCHS + 1)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        f"Experiment 2 — Killing Pretraining  (d={D}, pretrain={PRETRAIN_STEPS} steps)",
        fontsize=13,
    )

    ax = axes[0]
    for name, h in results.items():
        ax.plot(epochs, h["val_acc"], label=name, color=colors[name], lw=1.8)
    ax.set_xlabel("Epoch"); ax.set_ylabel("Val Accuracy")
    ax.set_title("Downstream Task Accuracy"); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1]
    for name, h in results.items():
        t_k = [v for v in h["T_K"] if not np.isnan(v)]
        if t_k:
            ax.semilogy(list(epochs)[:len(t_k)], t_k,
                        label=name, color=colors[name], lw=1.8)
    ax.set_xlabel("Epoch"); ax.set_ylabel("T_K  (log)")
    ax.set_title("Killing Tension  (→ 0 = Lie algebra)"); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[2]
    for name, h in results.items():
        ax.plot(epochs, h["T_closure"], label=name, color=colors[name], lw=1.8)
    ax.set_xlabel("Epoch"); ax.set_ylabel("T_closure")
    ax.set_title("Closure Fraction  (→ 0 = closed)"); ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Figure saved → {out_path}")


# ============================================================
#  MAIN
# ============================================================

def main():
    print("=" * 62)
    print("Experiment 2: Killing Pretraining as Algebraic Initialization")
    print("=" * 62)
    print(f"  Device  : {DEVICE}")
    if torch.cuda.is_available():
        print(f"  GPU     : {torch.cuda.get_device_name(0)}")
    print(f"  d={D}  n_gen={N_GEN}  n_in={N_IN}")
    print(f"  Pretrain: {PRETRAIN_STEPS} steps (T_K only, no data)")
    print(f"  Finetune: {FINETUNE_EPOCHS} epochs")

    out_dir = os.path.dirname(os.path.abspath(__file__))

    # ---- Dataset ----
    print(f"\n[1/6] Generating dataset  (N={N_SAMPLES})...")
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

    # ---- Phase 1: Killing pretraining ----
    print("\n[2/6] Phase 1 — Killing Pretraining (no data, pure T_K) ...")
    torch.manual_seed(0)
    pretrained_net = AlgebraNet()
    pre_hist = killing_pretrain(pretrained_net)

    # Post-pretraining analysis
    G_pre = pretrained_net.generators()
    m_pre = analyze_generators(G_pre)
    print_algebra_report(m_pre, title="Layer 1 after Killing Pretraining")

    plot_pretraining(pre_hist, os.path.join(out_dir, "exp2_pretraining.png"))

    # ---- Baselines ----
    # (A) Random init → full training
    print("\n[3/6] Baseline (A): Random init → full training ...")
    torch.manual_seed(0)
    rand_net = AlgebraNet()
    hist_rand = train_task(
        rand_net, train_loader, val_loader,
        n_epochs=FINETUNE_EPOCHS, freeze_l1=False, name="Random init",
    )

    # (B) Killing pretrain → fine-tune all layers
    print("\n[4/6] Baseline (B): Killing pretrain → fine-tune all layers ...")
    net_b = copy.deepcopy(pretrained_net)
    hist_ft_all = train_task(
        net_b, train_loader, val_loader,
        n_epochs=FINETUNE_EPOCHS, freeze_l1=False, name="Killing pretrain (all)",
    )

    # (C) Killing pretrain → fine-tune rest only (layer 1 frozen)
    print("\n[5/6] Baseline (C): Killing pretrain → layer 1 frozen ...")
    net_c = copy.deepcopy(pretrained_net)
    hist_ft_frozen = train_task(
        net_c, train_loader, val_loader,
        n_epochs=FINETUNE_EPOCHS, freeze_l1=True, name="Killing pretrain (L1 frozen)",
    )

    # ---- Post-training analysis ----
    print("\n[6/6] Post-training algebraic analysis ...")
    for name, net in [
        ("Random init (after full training)",       rand_net),
        ("Killing pretrain (after fine-tune all)",  net_b),
        ("Killing pretrain (after L1 frozen)",      net_c),
    ]:
        G = net.generators()
        m = analyze_generators(G)
        print_algebra_report(m, title=f"Layer 1 — {name}")

    # ---- Summary ----
    print("\n" + "=" * 62)
    print("SUMMARY")
    print("=" * 62)
    results = {
        "Random init":                  hist_rand,
        "Killing pretrain (all)":       hist_ft_all,
        "Killing pretrain (L1 frozen)": hist_ft_frozen,
    }
    for name, h in results.items():
        best  = max(h["val_acc"])
        final = h["val_acc"][-1]
        step_best = int(np.argmax(h["val_acc"])) + 1
        print(f"  {name:36s}  best={best:.4f} @ ep {step_best:3d}  final={final:.4f}")

    b_rand = max(hist_rand["val_acc"])
    b_ft   = max(hist_ft_all["val_acc"])
    b_fr   = max(hist_ft_frozen["val_acc"])
    print(f"\n  Δ(pretrain-all  vs random) : {b_ft - b_rand:+.4f}")
    print(f"  Δ(pretrain-frz  vs random) : {b_fr - b_rand:+.4f}")
    print(f"\n  T_K before pretraining     : {pre_hist['T_K'][0]:.4f}")
    print(f"  T_K after  pretraining     : {pre_hist['T_K'][-1]:.4f}")

    # ---- Save ----
    np.savez(
        os.path.join(out_dir, "exp2_results.npz"),
        pre_hist        = pre_hist,
        hist_rand       = hist_rand,
        hist_ft_all     = hist_ft_all,
        hist_ft_frozen  = hist_ft_frozen,
    )
    plot_finetuning(results, os.path.join(out_dir, "exp2_results.png"))

    print("\nExperiment 2 complete.")


if __name__ == "__main__":
    main()
