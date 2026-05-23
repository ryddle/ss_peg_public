"""
Experiment 5 -- Scaling Study: Variable n_gen and d-sweep
=========================================================

Two complementary experiments in one script:

Part A  (d=3, n_gen sweep):
    Compare models with different numbers of active generators for su(3):
      1. AlgNet-8   : 8 generators (full su(3) basis), λ_K=2.0
      2. AlgNet-5   : 5 generators explicitly, λ_K=2.0
      3. AlgNet-Gate: 8 generators + per-generator soft gate (group-Lasso)
         -> automatic rank discovery: starts with 8, sparsity prunes to k active
      4. GML-Fixed  : oracle with full 8 GML generators (upper bound)

Part B  (d-sweep, d=2..8):
    For each d: fit GML-Fixed, AlgNet-TK (n_gen=d²-1), MLP-baseline
    Task: sign(Tr(M(x)^3))  --  generalizes cleanly to any su(d)
    Metrics recorded per d:
      - val accuracy (best + final)
      - T_K final
      - θ orthogonality error
      - Killing eigenvalue uniformity
      - effective n_gen (from gate or rank)
      - GPU peak memory (MB)

Design notes:
  - input dim  = d²-1  (GML coordinates, one per generator)
  - head scales proportionally: Linear(n_gen -> min(4*n_gen, 128))
  - epochs scale inversely with d: epochs = max(400, 1500 // d)  (GPU budget)
  - λ_K = 2.0 throughout (validated in Exp4)
  - λ_gate = 0.05  (soft-gate L1 strength, tuned empirically for d=3)
"""

from __future__ import annotations

import os
import sys
import time
import warnings
import argparse

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ---------------------------------------------------------------------------
#  Path setup -- pick up killing_utils from the same directory
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from killing_utils import (
    killing_tension,
    closure_tension,
    analyze_generators,
    print_algebra_report,
    _su_basis,
    params_to_su_generators,
)

# ---------------------------------------------------------------------------
#  Global config
# ---------------------------------------------------------------------------
DEVICE     = torch.device("cuda" if torch.cuda.is_available() else "cpu")
LAMBDA_K   = 2.0       # Killing pressure (same as Exp4)
LAMBDA_GATE= 0.05      # Group-Lasso for gated model
LR         = 1e-3
BATCH_SIZE = 256
N_TRAIN    = 6_000
N_VAL      = 1_500
SEED       = 42

# Part A config
PART_A_D       = 3
PART_A_N_EPOCHS= 1000

# Part B config
D_VALUES       = [2, 3, 4, 5, 6, 8]  # d-sweep (remove 7 to save time; add if desired)
PART_B_N_EPOCHS_BASE = 1500           # scaled as max(400, base // d)

OUT_DIR = os.path.join(HERE, "exp5_outputs")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
#  Schedule helpers
# ---------------------------------------------------------------------------

def lam_schedule(epoch: int, lambda_k: float, warmup: int = 50) -> float:
    """Linear warmup then constant."""
    if epoch < warmup:
        return lambda_k * (epoch / warmup)
    return lambda_k


def epochs_for_d(d: int) -> int:
    return max(400, PART_B_N_EPOCHS_BASE // d)


# ============================================================
#  DATASET – generalised su(d) cubic invariant
# ============================================================

def generate_su_d_cubic_dataset(
    d: int,
    n: int,
    seed: int = SEED,
    device: torch.device = DEVICE,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Generate n samples of the su(d) cubic invariant classification task.

    x    ∈ R^(d²-1)   (GML coordinates, i.i.d. standard normal)
    M(x) = Σ_a x_a T_a  ∈ su(d)  (Hermitian traceless)
    y    = sign( Re(Tr(M(x)^3)) )   ∈ {0, 1}

    Re(Tr(M³)) is real for Hermitian M (since M†=M ⟹ Tr(M³)=Tr((M†)³)*=Tr(M³̄) ).
    The cubic invariant is genuinely SU(d)-symmetric.
    """
    n_gen  = d * d - 1
    rng    = np.random.default_rng(seed)
    X      = rng.standard_normal((n, n_gen)).astype(np.float64)

    basis_t = _su_basis(d, "cpu")
    basis   = basis_t.numpy()                              # (n_gen, d, d) complex128

    M_c  = np.einsum("na,aij->nij", X, basis)             # (N, d, d) complex
    M3_c = M_c @ M_c @ M_c                                 # uses matmul broadcasting
    tr3  = np.trace(M3_c, axis1=-2, axis2=-1).real        # (N,) real

    y = (tr3 > 0).astype(np.int64)
    pos = y.sum()
    print(f"  su({d}) dataset  n={n}  balance: {n-pos} neg / {pos} pos  "
          f"(frac={pos/n:.3f})")
    max_imag = np.abs(np.trace(M3_c, axis1=-2, axis2=-1).imag).max()
    if max_imag > 1e-8:
        print(f"  WARNING: max|Im(Tr(M³))| = {max_imag:.2e}")

    return (
        torch.tensor(X, dtype=torch.float64, device=device),
        torch.tensor(y, dtype=torch.long,    device=device),
    )


def make_loaders(d: int, seed: int = SEED):
    X_tr, y_tr = generate_su_d_cubic_dataset(d, N_TRAIN, seed=seed)
    X_va, y_va = generate_su_d_cubic_dataset(d, N_VAL,   seed=seed + 1)
    tr = DataLoader(TensorDataset(X_tr, y_tr), batch_size=BATCH_SIZE, shuffle=True)
    va = DataLoader(TensorDataset(X_va, y_va), batch_size=512,         shuffle=False)
    return tr, va


# ============================================================
#  HEAD FACTORY (scales with n_gen)
# ============================================================

def make_head(n_gen: int) -> nn.Sequential:
    hidden = min(4 * n_gen, 128)
    return nn.Sequential(
        nn.ReLU(),
        nn.Linear(n_gen, hidden, dtype=torch.float64),
        nn.ReLU(),
        nn.Linear(hidden, 2, dtype=torch.float64),
    )


# ============================================================
#  ALGEBRAIC PROJECTION  (generalised to arbitrary d, n_gen_active)
# ============================================================

class AlgebraicProjection(nn.Module):
    """
    Algebraic first layer for su(d) with configurable number of active generators.

    G_a = Σ_k θ_{ak} T_k   (a=0..n_gen_active-1,  k=0..d²-2)

    f_a = Re(Tr(G_a @ M(x)))  =  2 Σ_k θ_{ak} x_k   (linear map with W=2θ)

    Parameters
    ----------
    d            : algebra dimension
    n_gen_active : how many generators to learn  (≤ d²-1)
    fixed        : if True, θ = I[:n_gen_active, :] fixed (oracle)
    lambda_k     : Killing pressure coefficient
    """

    def __init__(
        self,
        d           : int,
        n_gen_active: int,
        fixed       : bool  = False,
        lambda_k    : float = 0.0,
    ):
        super().__init__()
        self.d            = d
        self.n_gen_full   = d * d - 1          # full basis size
        self.n_gen_active = n_gen_active
        self.fixed        = fixed
        self.lambda_k     = lambda_k

        # θ shape: (n_gen_active, n_gen_full)
        # Oracle init = I[:n_gen_active, :]  ->  first n_gen_active GML generators
        init = torch.eye(n_gen_active, self.n_gen_full, dtype=torch.float64)
        if fixed:
            self.register_buffer("theta_buf", init)
        else:
            noise = torch.randn_like(init) * 0.05
            self.theta = nn.Parameter(init + noise)

        self.bias = nn.Parameter(torch.zeros(n_gen_active, dtype=torch.float64))

    def _theta(self) -> torch.Tensor:
        return self.theta_buf if self.fixed else self.theta

    def get_generators(self) -> torch.Tensor:
        """(n_gen_active, d, d) complex on same device as θ."""
        return params_to_su_generators(self._theta(), self.d)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        dev   = str(self._theta().device)
        basis = _su_basis(self.d, dev)                       # (n_gen_full, d, d) complex
        M_x   = torch.einsum(
            "...b,bij->...ij",
            x.to(torch.complex128),
            basis,
        )                                                    # (B, d, d) complex
        G  = self.get_generators()                           # (n_active, d, d) complex
        f  = torch.einsum("aij,...ji->...a", G, M_x).real   # (B, n_active) -- no conj on G
        return f + self.bias

    def reg_loss(self, epoch: int = 0) -> torch.Tensor:
        dev = str(self._theta().device)
        if self.fixed or self.lambda_k == 0.0:
            return torch.zeros(1, dtype=torch.float64, device=dev).squeeze()
        lk = lam_schedule(epoch, self.lambda_k)
        return lk * killing_tension(self.get_generators())


# ============================================================
#  GATED ALGEBRAIC PROJECTION  (automatic rank discovery)
# ============================================================

class GatedAlgebraicProjection(nn.Module):
    """
    Algebraic layer with per-generator learnable scale gates.

    G_active_a = s_a * G_a    where  s_a = softplus(log_s_a)  ∈ (0, ∞)

    Group-Lasso regularization:   λ_gate * Σ_a s_a
    drives unused generators toward s_a≈0 -> effective rank collapses automatically.

    Monitor effective count:  (s > threshold).sum()
    """

    def __init__(
        self,
        d          : int,
        lambda_k   : float = LAMBDA_K,
        lambda_gate: float = LAMBDA_GATE,
    ):
        super().__init__()
        self.d           = d
        self.n_gen_full  = d * d - 1
        self.lambda_k    = lambda_k
        self.lambda_gate = lambda_gate

        init  = torch.eye(self.n_gen_full, dtype=torch.float64)
        noise = torch.randn_like(init) * 0.05
        self.theta   = nn.Parameter(init + noise)
        self.log_s   = nn.Parameter(torch.zeros(self.n_gen_full, dtype=torch.float64))
        self.bias    = nn.Parameter(torch.zeros(self.n_gen_full, dtype=torch.float64))

    def scales(self) -> torch.Tensor:
        """Positive per-generator scales via softplus."""
        return nn.functional.softplus(self.log_s)          # (n_gen,) > 0

    def get_generators(self) -> torch.Tensor:
        """(n_gen, d, d) scaled generators for Killing analysis."""
        G_raw = params_to_su_generators(self.theta, self.d)  # (n, d, d)
        s     = self.scales()
        return G_raw * s[:, None, None]

    def effective_n_gen(self, threshold: float = 0.1) -> int:
        with torch.no_grad():
            return int((self.scales() > threshold).sum().item())

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        dev   = str(self.theta.device)
        basis = _su_basis(self.d, dev)
        M_x   = torch.einsum(
            "...b,bij->...ij",
            x.to(torch.complex128),
            basis,
        )
        G  = self.get_generators()                           # scaled generators
        f  = torch.einsum("aij,...ji->...a", G, M_x).real
        return f + self.bias

    def reg_loss(self, epoch: int = 0) -> torch.Tensor:
        dev   = str(self.theta.device)
        lk    = lam_schedule(epoch, self.lambda_k)
        # Killing pressure on scaled generators
        tk    = lk * killing_tension(self.get_generators())
        # Group-Lasso on scales (drives unused generators to zero)
        gl    = self.lambda_gate * self.scales().sum()
        return tk + gl


# ============================================================
#  FULL MODELS
# ============================================================

class AlgebraicNet(nn.Module):
    def __init__(
        self,
        d           : int,
        n_gen_active: int,
        fixed       : bool  = False,
        lambda_k    : float = 0.0,
    ):
        super().__init__()
        self.alg_layer = AlgebraicProjection(
            d=d, n_gen_active=n_gen_active, fixed=fixed, lambda_k=lambda_k
        )
        self.head = make_head(n_gen_active)

    def forward(self, x):
        return self.head(self.alg_layer(x))

    def reg_loss(self, epoch=0):
        return self.alg_layer.reg_loss(epoch)

    def get_generators(self):
        with torch.no_grad():
            return self.alg_layer.get_generators()


class GatedAlgebraicNet(nn.Module):
    def __init__(self, d: int, lambda_k=LAMBDA_K, lambda_gate=LAMBDA_GATE):
        super().__init__()
        n_gen = d * d - 1
        self.alg_layer = GatedAlgebraicProjection(
            d=d, lambda_k=lambda_k, lambda_gate=lambda_gate
        )
        self.head = make_head(n_gen)

    def forward(self, x):
        return self.head(self.alg_layer(x))

    def reg_loss(self, epoch=0):
        return self.alg_layer.reg_loss(epoch)

    def get_generators(self):
        with torch.no_grad():
            return self.alg_layer.get_generators()

    def effective_n_gen(self, threshold=0.1):
        return self.alg_layer.effective_n_gen(threshold)


class StandardMLP(nn.Module):
    def __init__(self, d: int):
        super().__init__()
        n_gen = d * d - 1
        self.layer1 = nn.Linear(n_gen, n_gen, dtype=torch.float64)
        self.head   = make_head(n_gen)

    def forward(self, x):
        return self.head(self.layer1(x))


# ============================================================
#  TRAINING LOOP
# ============================================================

def train(
    model      : nn.Module,
    train_loader,
    val_loader,
    name       : str   = "Model",
    n_epochs   : int   = 1000,
    quiet      : bool  = False,
    report_every: int  = 100,
) -> dict:
    model = model.to(DEVICE)
    optimizer = optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=LR,
    )
    criterion = nn.CrossEntropyLoss()

    hist = {
        "train_acc": [], "val_acc":  [],
        "train_loss": [], "val_loss": [],
        "T_K": [], "T_closure": [],
        "eff_n_gen": [],
    }
    best_val   = 0.0
    best_epoch = 0

    if not quiet:
        print(f"\n  [{name}]  {n_epochs} epochs ...")

    for epoch in range(n_epochs):
        model.train()
        t_loss, correct, total = 0.0, 0, 0
        for xb, yb in train_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            optimizer.zero_grad()
            logits    = model(xb)
            task_loss = criterion(logits, yb)
            loss = task_loss
            if hasattr(model, "reg_loss"):
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
                out    = model(xb)
                v_loss += criterion(out, yb).item() * xb.size(0)
                v_correct += (out.argmax(1) == yb).sum().item()
                v_total   += xb.size(0)

        val_acc = v_correct / v_total
        if val_acc > best_val:
            best_val   = val_acc
            best_epoch = epoch + 1

        # Algebraic diagnostics
        t_k = t_c = float("nan")
        eff  = -1
        if hasattr(model, "get_generators"):
            with torch.no_grad():
                G   = model.get_generators()
                t_k = killing_tension(G).item()
                t_c = closure_tension(G).item()
        if isinstance(model, GatedAlgebraicNet):
            eff = model.effective_n_gen()

        hist["train_loss"].append(t_loss / total)
        hist["train_acc"].append(correct / total)
        hist["val_loss"].append(v_loss / v_total)
        hist["val_acc"].append(val_acc)
        hist["T_K"].append(t_k)
        hist["T_closure"].append(t_c)
        hist["eff_n_gen"].append(eff)

        if not quiet and (epoch + 1) % report_every == 0:
            eff_str = f"  eff_n={eff}" if eff >= 0 else ""
            print(
                f"    epoch {epoch+1:4d}/{n_epochs}"
                f"  val_acc={val_acc:.4f}"
                f"  T_K={t_k:.4f}  T_cl={t_c:.4f}{eff_str}"
            )

    hist["best_val"]   = best_val
    hist["best_epoch"] = best_epoch
    return hist


# ============================================================
#  PART A -- n_gen sweep for su(3)
# ============================================================

def run_part_a():
    print("\n" + "=" * 70)
    print("  Part A: Variable n_gen for su(3)")
    print("=" * 70)

    d      = PART_A_D
    n_gen  = d * d - 1   # 8 for su(3)
    epochs = PART_A_N_EPOCHS
    tr, va = make_loaders(d)

    models = {
        "GML-Fixed-8" : AlgebraicNet(d=d, n_gen_active=n_gen, fixed=True,  lambda_k=0.0),
        "AlgNet-TK-8" : AlgebraicNet(d=d, n_gen_active=n_gen, fixed=False, lambda_k=LAMBDA_K),
        "AlgNet-TK-5" : AlgebraicNet(d=d, n_gen_active=5,     fixed=False, lambda_k=LAMBDA_K),
        "GML-Fixed-5" : AlgebraicNet(d=d, n_gen_active=5,     fixed=True,  lambda_k=0.0),
        "AlgNet-Gated": GatedAlgebraicNet(d=d, lambda_k=LAMBDA_K, lambda_gate=LAMBDA_GATE),
        "MLP-8"       : StandardMLP(d=d),
    }

    results = {}
    for name, model in models.items():
        torch.cuda.reset_peak_memory_stats() if DEVICE.type == "cuda" else None
        t0 = time.perf_counter()
        results[name] = train(model, tr, va, name=name, n_epochs=epochs, report_every=200)
        results[name]["time_s"]   = time.perf_counter() - t0
        results[name]["peak_mem_mb"] = (
            torch.cuda.max_memory_allocated() / 1e6 if DEVICE.type == "cuda" else 0.0
        )

    # Print summary
    print("\n  PART A SUMMARY")
    print(f"  {'Model':<20} {'Best acc':>10} {'@epoch':>8} {'Final acc':>10} "
          f"{'T_K final':>12} {'eff_n_gen':>10}")
    print("  " + "-" * 75)
    for name, h in results.items():
        tk  = h["T_K"][-1]
        eff = h["eff_n_gen"][-1]
        eff_str = f"{eff:>10d}" if eff >= 0 else f"{'--':>10}"
        tk_str  = f"{tk:12.4f}" if not np.isnan(tk) else f"{'--':>12}"
        print(
            f"  {name:<20} {h['best_val']:>10.4f}  {h['best_epoch']:>7d}  "
            f"{h['val_acc'][-1]:>10.4f}  {tk_str}{eff_str}"
        )

    # Plot Part A
    _plot_part_a(results, epochs, out_path=os.path.join(OUT_DIR, "exp5_part_a.png"))

    # Algebraic report for each AlgNet model
    print("\n  PART A -- Algebraic analysis:")
    for name, model in models.items():
        if hasattr(model, "get_generators"):
            G = model.get_generators().to(DEVICE)
            metrics = analyze_generators(G)
            print_algebra_report(metrics, title=name)

    # Save scales for gated model
    gated_model = models["AlgNet-Gated"]
    with torch.no_grad():
        scales = gated_model.alg_layer.scales().cpu().numpy()
    print(f"\n  AlgNet-Gated final gate scales: {np.round(scales, 3)}")
    print(f"  Effective n_gen (s>0.1): {(scales > 0.1).sum()}")
    print(f"  Effective n_gen (s>0.3): {(scales > 0.3).sum()}")

    # Save numeric results
    np.savez(
        os.path.join(OUT_DIR, "exp5_part_a.npz"),
        model_names=list(results.keys()),
        **{k.replace("-", "_") + "_hist": np.array(v["val_acc"])
           for k, v in results.items()},
        gate_scales=scales,
    )
    print(f"  Saved -> {OUT_DIR}/exp5_part_a.npz")
    return results


def _plot_part_a(results: dict, n_epochs: int, out_path: str) -> None:
    epochs = range(1, n_epochs + 1)

    color_map = {
        "GML-Fixed-8" : ("tab:green",  "-."),
        "GML-Fixed-5" : ("tab:olive",  ":"),
        "AlgNet-TK-8" : ("tab:orange", "-"),
        "AlgNet-TK-5" : ("tab:red",    "-"),
        "AlgNet-Gated": ("tab:purple", "--"),
        "MLP-8"       : ("tab:blue",   "-"),
    }

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("Exp 5A -- Variable n_gen, su(3)", fontsize=12)

    # Val accuracy
    ax = axes[0]
    for name, h in results.items():
        c, ls = color_map.get(name, ("black", "-"))
        ax.plot(epochs, h["val_acc"], label=name, color=c, ls=ls, lw=1.8)
    ax.set_title("Validation Accuracy")
    ax.set_xlabel("Epoch"); ax.set_ylabel("Accuracy")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # T_K
    ax = axes[1]
    for name, h in results.items():
        c, ls = color_map.get(name, ("black", "-"))
        vals = h["T_K"]
        if not all(np.isnan(v) for v in vals):
            y_plot = [max(v, 1e-5) for v in vals]
            ax.semilogy(list(epochs), y_plot, label=name, color=c, ls=ls, lw=1.8)
    ax.set_title("Killing Tension T_K")
    ax.set_xlabel("Epoch"); ax.set_ylabel("T_K (log)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Effective n_gen for gated model
    ax = axes[2]
    if "AlgNet-Gated" in results:
        eff = results["AlgNet-Gated"]["eff_n_gen"]
        ax.plot(list(epochs), eff, color="tab:purple", lw=2, label="AlgNet-Gated eff_n")
        ax.set_ylabel("Effective n_gen (scale > 0.1)")
    ax.set_title("Automatic Generator Discovery (Gated)")
    ax.set_xlabel("Epoch")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close()
    print(f"\n  Part A figure -> {out_path}")


# ============================================================
#  PART B -- d-scaling sweep
# ============================================================

def run_part_b():
    print("\n" + "=" * 70)
    print("  Part B: d-scaling sweep  (d ∈ {})".format(D_VALUES))
    print("=" * 70)

    scaling_results = {}

    for d in D_VALUES:
        n_gen  = d * d - 1
        epochs = epochs_for_d(d)

        print(f"\n  --- d={d}  (su({d}), n_gen={n_gen}, epochs={epochs}) ---")
        tr, va = make_loaders(d, seed=SEED)

        models_d = {
            f"GML-Fixed-d{d}" : AlgebraicNet(d=d, n_gen_active=n_gen, fixed=True,  lambda_k=0.0),
            f"AlgNet-TK-d{d}" : AlgebraicNet(d=d, n_gen_active=n_gen, fixed=False, lambda_k=LAMBDA_K),
            f"MLP-d{d}"       : StandardMLP(d=d),
        }

        d_res = {}
        for name, model in models_d.items():
            torch.cuda.reset_peak_memory_stats() if DEVICE.type == "cuda" else None
            t0 = time.perf_counter()
            h  = train(
                model, tr, va,
                name=name, n_epochs=epochs,
                quiet=False, report_every=max(100, epochs // 5),
            )
            h["time_s"]      = time.perf_counter() - t0
            h["peak_mem_mb"] = (
                torch.cuda.max_memory_allocated() / 1e6 if DEVICE.type == "cuda" else 0.0
            )

            # Final algebraic analysis
            if hasattr(model, "get_generators"):
                G = model.get_generators().to(DEVICE)
                metrics = analyze_generators(G)
                h["ortho_error"]    = float(
                    np.linalg.norm(
                        model.alg_layer._theta().detach().cpu().numpy()
                        @ model.alg_layer._theta().detach().cpu().numpy().T
                        - np.eye(n_gen)
                    )
                    if not model.alg_layer.fixed else 0.0
                )
                h["eigenvalue_uniformity"] = metrics["eigenvalue_uniformity"]
                h["killing_tension_final"] = metrics["killing_tension"]
            else:
                h["ortho_error"]           = float("nan")
                h["eigenvalue_uniformity"] = float("nan")
                h["killing_tension_final"] = float("nan")

            d_res[name] = h

        scaling_results[d] = d_res

        # Brief summary for this d
        print(f"\n  d={d}  results:")
        for name, h in d_res.items():
            tk = h["T_K"][-1]
            print(f"    {name:<22}  best={h['best_val']:.4f}  "
                  f"T_K={tk:.4f}  mem={h['peak_mem_mb']:.1f}MB  "
                  f"time={h['time_s']:.1f}s")

    # Print big summary table
    print("\n\n  PART B SCALING SUMMARY")
    print(f"  {'d':>4}  {'Model':<22}  {'Best acc':>10}  {'T_K final':>12}  "
          f"{'OrthoErr':>10}  {'EigUnif':>9}  {'Mem MB':>8}")
    print("  " + "-" * 85)
    for d, d_res in scaling_results.items():
        n_gen = d * d - 1
        oracle_acc  = d_res.get(f"GML-Fixed-d{d}", {}).get("best_val", float("nan"))
        algnet_acc  = d_res.get(f"AlgNet-TK-d{d}", {}).get("best_val", float("nan"))
        mlp_acc     = d_res.get(f"MLP-d{d}", {}).get("best_val", float("nan"))
        gap_alg_oracle = algnet_acc - oracle_acc
        gap_alg_mlp    = algnet_acc - mlp_acc
        for name, h in d_res.items():
            tk = h.get("killing_tension_final", h["T_K"][-1])
            print(
                f"  {d:>4}  {name:<22}  {h['best_val']:>10.4f}  "
                f"{tk:>12.4f}  {h.get('ortho_error', float('nan')):>10.4f}  "
                f"{h.get('eigenvalue_uniformity', float('nan')):>9.4f}  "
                f"{h['peak_mem_mb']:>8.1f}"
            )
        print(f"  {d:>4}  {'  Δ(TK - Oracle)':>22}  {gap_alg_oracle:>+10.4f}")
        print(f"  {d:>4}  {'  Δ(TK - MLP)':>22}  {gap_alg_mlp:>+10.4f}")
        print()

    # Save numeric
    save_dict = {}
    for d, d_res in scaling_results.items():
        for name, h in d_res.items():
            key = name.replace("-", "_")
            save_dict[key + "_val_acc"] = np.array(h["val_acc"])
            save_dict[key + "_T_K"]     = np.array(h["T_K"])
    save_dict["d_values"] = np.array(D_VALUES)
    np.savez(os.path.join(OUT_DIR, "exp5_part_b.npz"), **save_dict)
    print(f"  Saved -> {OUT_DIR}/exp5_part_b.npz")

    # Plot Part B
    _plot_part_b(scaling_results, out_path=os.path.join(OUT_DIR, "exp5_part_b.png"))

    return scaling_results


def _plot_part_b(scaling_results: dict, out_path: str) -> None:
    """
    Multi-panel scaling figure:
      Row 1: val-accuracy curves per d (one sub-panel per d)
      Row 2: scalar metrics vs d  (accuracy gap, T_K, ortho error, eigenvalue uniformity)
    """
    d_list  = sorted(scaling_results.keys())
    n_d     = len(d_list)

    fig = plt.figure(figsize=(4 * n_d, 14))
    gs  = fig.add_gridspec(3, n_d, hspace=0.35, wspace=0.3)

    c_map = {
        "GML"  : "tab:green",
        "Alg"  : "tab:orange",
        "MLP"  : "tab:blue",
    }

    # --- Row 0: Validation accuracy curves per d ---
    for col, d in enumerate(d_list):
        ax = fig.add_subplot(gs[0, col])
        d_res   = scaling_results[d]
        epochs  = range(1, len(next(iter(d_res.values()))["val_acc"]) + 1)
        for name, h in d_res.items():
            prefix = "GML" if "GML" in name else ("Alg" if "Alg" in name else "MLP")
            ax.plot(list(epochs), h["val_acc"],
                    color=c_map[prefix], lw=1.5,
                    label=prefix)
        ax.set_title(f"su({d})  n_gen={d*d-1}", fontsize=9)
        ax.set_xlabel("Epoch", fontsize=8)
        ax.set_ylabel("Val Acc", fontsize=8)
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # --- Row 1: T_K curves per d ---
    for col, d in enumerate(d_list):
        ax = fig.add_subplot(gs[1, col])
        d_res   = scaling_results[d]
        epochs  = range(1, len(next(iter(d_res.values()))["T_K"]) + 1)
        for name, h in d_res.items():
            vals = h["T_K"]
            if not all(np.isnan(v) for v in vals):
                prefix = "GML" if "GML" in name else ("Alg" if "Alg" in name else "MLP")
                y_plot = [max(v, 1e-5) if not np.isnan(v) else 1e-5 for v in vals]
                ax.semilogy(list(epochs), y_plot,
                            color=c_map[prefix], lw=1.5, label=prefix)
        ax.set_title(f"T_K  su({d})", fontsize=9)
        ax.set_xlabel("Epoch", fontsize=8)
        ax.set_ylabel("T_K (log)", fontsize=8)
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # --- Row 2: scalar metrics vs d ---
    d_arr          = np.array(d_list)
    oracle_accs    = []
    algnet_accs    = []
    mlp_accs       = []
    gap_alg_oracle = []
    gap_alg_mlp    = []
    tk_finals      = []
    ortho_errs     = []
    eig_unis       = []

    for d in d_list:
        d_res = scaling_results[d]
        gml_h = d_res.get(f"GML-Fixed-d{d}", {})
        alg_h = d_res.get(f"AlgNet-TK-d{d}", {})
        mlp_h = d_res.get(f"MLP-d{d}", {})

        oa = gml_h.get("best_val", float("nan"))
        aa = alg_h.get("best_val", float("nan"))
        ma = mlp_h.get("best_val", float("nan"))

        oracle_accs.append(oa)
        algnet_accs.append(aa)
        mlp_accs.append(ma)
        gap_alg_oracle.append(aa - oa)
        gap_alg_mlp.append(aa - ma)
        tk_finals.append(alg_h.get("killing_tension_final", alg_h["T_K"][-1] if alg_h else float("nan")))
        ortho_errs.append(alg_h.get("ortho_error", float("nan")))
        eig_unis.append(alg_h.get("eigenvalue_uniformity", float("nan")))

    # accuracy vs d
    ax = fig.add_subplot(gs[2, :2])
    ax.plot(d_arr, oracle_accs, "s-",  color="tab:green",  lw=2, label="GML-Fixed (oracle)")
    ax.plot(d_arr, algnet_accs, "o-",  color="tab:orange", lw=2, label="AlgNet-TK")
    ax.plot(d_arr, mlp_accs,    "^-",  color="tab:blue",   lw=2, label="MLP baseline")
    ax.set_title("Best Accuracy vs d", fontsize=10)
    ax.set_xlabel("d"); ax.set_ylabel("Best Val Accuracy")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    # accuracy gap vs d
    ax2 = ax.twinx()
    ax2.plot(d_arr, gap_alg_mlp, "D--", color="tab:red", lw=1.5, label="Δ(TK−MLP)")
    ax2.axhline(0, color="gray", lw=0.8, ls=":")
    ax2.set_ylabel("Accuracy Gap (AlgNet−MLP)", color="tab:red", fontsize=8)
    ax2.legend(loc="lower right", fontsize=8)

    # T_K + ortho error vs d
    ax = fig.add_subplot(gs[2, 2:4])
    ax.semilogy(d_arr, np.maximum(tk_finals, 1e-6), "o-",  color="tab:orange", lw=2, label="T_K final")
    ax.set_title("T_K final & Ortho Error vs d", fontsize=10)
    ax.set_xlabel("d"); ax.set_ylabel("T_K (log)", color="tab:orange", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(d_arr, ortho_errs, "s--", color="tab:purple", lw=1.8, label="θ ortho err")
    ax2.set_ylabel("θ ortho error", color="tab:purple", fontsize=9)
    ax2.legend(fontsize=8)
    ax.legend(fontsize=8, loc="upper left")

    if n_d > 4:
        ax = fig.add_subplot(gs[2, 4:])
        ax.plot(d_arr, eig_unis, "v-", color="tab:cyan", lw=2, label="Eig uniformity")
        ax.set_title("Killing Eigenvalue Uniformity vs d", fontsize=10)
        ax.set_xlabel("d"); ax.set_ylabel("Uniformity (lower=better)")
        ax.grid(True, alpha=0.3); ax.legend(fontsize=8)

    fig.suptitle(
        "Exp 5B -- Scaling Study: su(d) Cubic Invariant  "
        f"(d ∈ {d_list}, λ_K={LAMBDA_K})",
        fontsize=12, y=1.01,
    )
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"\n  Part B figure -> {out_path}")


# ============================================================
#  PART C -- lambda_gate sweep (automatic rank discovery tuning)
# ============================================================

# lambda_gate values to compare; 0.05 is the reference from Part A
PART_C_LAMBDAS = [0.005, 0.01, 0.02, 0.05]
PART_C_N_EPOCHS = 1000


def run_part_c():
    print("\n" + "=" * 70)
    print("  Part C: lambda_gate sweep for GatedAlgebraicNet, su(3)")
    print("  Testing whether smaller lambda_gate enables real rank selection")
    print("=" * 70)

    d      = PART_A_D           # 3
    n_gen  = d * d - 1          # 8
    epochs = PART_C_N_EPOCHS
    tr, va = make_loaders(d)

    # Oracle reference (fixed, no gate)
    oracle = AlgebraicNet(d=d, n_gen_active=n_gen, fixed=True, lambda_k=0.0)
    print(f"\n  [GML-Fixed-8  (oracle reference)]  {epochs} epochs ...")
    torch.cuda.reset_peak_memory_stats() if DEVICE.type == "cuda" else None
    oracle_hist = train(oracle, tr, va, name="GML-Fixed-8", n_epochs=epochs, report_every=200)

    results   = {"GML-Fixed-8": oracle_hist}
    gate_models = {}

    for lg in PART_C_LAMBDAS:
        name = f"Gated-lg{lg}"
        model = GatedAlgebraicNet(d=d, lambda_k=LAMBDA_K, lambda_gate=lg)
        torch.cuda.reset_peak_memory_stats() if DEVICE.type == "cuda" else None
        t0  = time.perf_counter()
        h   = train(model, tr, va, name=name, n_epochs=epochs, report_every=200)
        h["time_s"]      = time.perf_counter() - t0
        h["peak_mem_mb"] = torch.cuda.max_memory_allocated() / 1e6 if DEVICE.type == "cuda" else 0.0
        results[name]    = h
        gate_models[name] = model

    # Summary table
    print("\n  PART C SUMMARY")
    print(f"  {'Model':<20} {'Best acc':>10} {'@epoch':>8} {'Final acc':>10} "
          f"{'T_K final':>12} {'eff_n(>0.1)':>12} {'eff_n(>0.3)':>12}")
    print("  " + "-" * 88)

    for name, h in results.items():
        tk  = h["T_K"][-1]
        eff = h["eff_n_gen"][-1]
        eff_str1 = f"{eff:>12d}" if eff >= 0 else f"{'--':>12}"
        tk_str   = f"{tk:12.4f}" if not np.isnan(tk) else f"{'--':>12}"
        # secondary threshold 0.3
        if name in gate_models:
            with torch.no_grad():
                s  = gate_models[name].alg_layer.scales().cpu().numpy()
                e3 = int((s > 0.3).sum())
        else:
            e3 = -1
        eff_str2 = f"{e3:>12d}" if e3 >= 0 else f"{'--':>12}"
        print(
            f"  {name:<20} {h['best_val']:>10.4f}  {h['best_epoch']:>7d}  "
            f"{h['val_acc'][-1]:>10.4f}  {tk_str}{eff_str1}{eff_str2}"
        )

    # Gate scale trajectories + accuracy plot
    _plot_part_c(results, gate_models, epochs,
                 out_path=os.path.join(OUT_DIR, "exp5_part_c.png"))

    # Per-model scale histograms
    print("\n  Final gate scales per model:")
    for name, model in gate_models.items():
        with torch.no_grad():
            s = model.alg_layer.scales().cpu().numpy()
        print(f"    {name:<20}  scales={np.round(s, 3)}  "
              f"eff>0.1={int((s>0.1).sum())}  eff>0.3={int((s>0.3).sum())}")

    np.savez(
        os.path.join(OUT_DIR, "exp5_part_c.npz"),
        model_names=list(results.keys()),
        lambda_gate_values=np.array(PART_C_LAMBDAS),
        **{k.replace("-", "_").replace(".", "p") + "_val_acc": np.array(v["val_acc"])
           for k, v in results.items()},
        **{k.replace("-", "_").replace(".", "p") + "_eff_n": np.array(v["eff_n_gen"])
           for k, v in results.items() if k != "GML-Fixed-8"},
    )
    print(f"  Saved -> {OUT_DIR}/exp5_part_c.npz")
    return results


def _plot_part_c(results: dict, gate_models: dict, n_epochs: int, out_path: str) -> None:
    epochs   = list(range(1, n_epochs + 1))
    lg_vals  = PART_C_LAMBDAS
    cmap     = plt.get_cmap("plasma")
    colors   = {f"Gated-lg{lg}": cmap(i / max(len(lg_vals) - 1, 1))
                for i, lg in enumerate(lg_vals)}
    oracle_c = "tab:green"

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    fig.suptitle(
        f"Exp 5C -- lambda_gate sweep, GatedAlgebraicNet, su(3), {n_epochs} epochs",
        fontsize=12,
    )

    # Val accuracy
    ax = axes[0]
    ax.plot(epochs, results["GML-Fixed-8"]["val_acc"],
            color=oracle_c, ls="-.", lw=2, label="GML-Fixed-8 (oracle)")
    for lg in lg_vals:
        name = f"Gated-lg{lg}"
        ax.plot(epochs, results[name]["val_acc"],
                color=colors[name], lw=1.8, label=f"lg={lg}")
    ax.set_title("Validation Accuracy")
    ax.set_xlabel("Epoch"); ax.set_ylabel("Accuracy")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    # Effective n_gen over time
    ax = axes[1]
    for lg in lg_vals:
        name = f"Gated-lg{lg}"
        eff  = results[name]["eff_n_gen"]
        ax.plot(epochs, eff, color=colors[name], lw=1.8, label=f"lg={lg}")
    ax.set_title("Effective n_gen (scale > 0.1)")
    ax.set_xlabel("Epoch"); ax.set_ylabel("Active generators")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    # Final scale distributions
    ax = axes[2]
    x_pos = np.arange(PART_A_D * PART_A_D - 1)  # 0..7
    for i, lg in enumerate(lg_vals):
        name = f"Gated-lg{lg}"
        if name in gate_models:
            with torch.no_grad():
                s = gate_models[name].alg_layer.scales().cpu().numpy()
            ax.bar(x_pos + i * 0.18 - 0.27, s, width=0.16,
                   color=colors[name], alpha=0.8, label=f"lg={lg}")
    ax.axhline(0.1, color="gray", ls="--", lw=1, label="threshold 0.1")
    ax.axhline(0.3, color="black", ls=":",  lw=1, label="threshold 0.3")
    ax.set_title("Final Gate Scales per Generator")
    ax.set_xlabel("Generator index"); ax.set_ylabel("Scale s_a")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close()
    print(f"\n  Part C figure -> {out_path}")


# ============================================================
#  PART D -- extended epochs for d >= 4
# ============================================================

PART_D_D_VALUES  = [4, 5, 6, 8]
PART_D_N_EPOCHS  = 1000


def run_part_d():
    print("\n" + "=" * 70)
    print(f"  Part D: Extended epochs ({PART_D_N_EPOCHS}) for d >= 4")
    print(f"  d values: {PART_D_D_VALUES}")
    print("=" * 70)

    scaling_results = {}

    for d in PART_D_D_VALUES:
        n_gen  = d * d - 1
        epochs = PART_D_N_EPOCHS

        print(f"\n  --- d={d}  (su({d}), n_gen={n_gen}, epochs={epochs}) ---")
        tr, va = make_loaders(d, seed=SEED)

        models_d = {
            f"GML-Fixed-d{d}" : AlgebraicNet(d=d, n_gen_active=n_gen, fixed=True,  lambda_k=0.0),
            f"AlgNet-TK-d{d}" : AlgebraicNet(d=d, n_gen_active=n_gen, fixed=False, lambda_k=LAMBDA_K),
            f"MLP-d{d}"       : StandardMLP(d=d),
        }

        d_res = {}
        for name, model in models_d.items():
            torch.cuda.reset_peak_memory_stats() if DEVICE.type == "cuda" else None
            t0 = time.perf_counter()
            h  = train(
                model, tr, va,
                name=name, n_epochs=epochs,
                quiet=False, report_every=max(100, epochs // 5),
            )
            h["time_s"]      = time.perf_counter() - t0
            h["peak_mem_mb"] = (
                torch.cuda.max_memory_allocated() / 1e6 if DEVICE.type == "cuda" else 0.0
            )

            if hasattr(model, "get_generators") and not getattr(
                getattr(model, "alg_layer", None), "fixed", False
            ):
                alg_layer = model.alg_layer
                theta_np  = alg_layer._theta().detach().cpu().numpy()
                h["ortho_error"] = float(
                    np.linalg.norm(theta_np @ theta_np.T - np.eye(n_gen))
                )
            else:
                h["ortho_error"] = float("nan")

            d_res[name] = h

        scaling_results[d] = d_res

        print(f"\n  d={d}  results:")
        for name, h in d_res.items():
            tk = h["T_K"][-1]
            tk_str = f"{tk:.4f}" if not np.isnan(tk) else "nan"
            print(f"    {name:<22}  best={h['best_val']:.4f}  "
                  f"T_K={tk_str}  mem={h['peak_mem_mb']:.1f}MB  "
                  f"time={h['time_s']:.1f}s")

    # Summary
    print("\n\n  PART D SCALING SUMMARY (1000 epochs)")
    print(f"  {'d':>4}  {'Model':<22}  {'Best acc':>10}  {'T_K final':>12}  "
          f"{'OrthoErr':>10}  {'Mem MB':>8}")
    print("  " + "-" * 72)
    for d, d_res in scaling_results.items():
        n_gen = d * d - 1
        for name, h in d_res.items():
            tk = h["T_K"][-1]
            tk_str = f"{tk:12.4f}" if not np.isnan(tk) else f"{'nan':>12}"
            print(
                f"  {d:>4}  {name:<22}  {h['best_val']:>10.4f}  "
                f"{tk_str}  {h.get('ortho_error', float('nan')):>10.4f}  "
                f"{h['peak_mem_mb']:>8.1f}"
            )
        alg_acc = d_res.get(f"AlgNet-TK-d{d}", {}).get("best_val", float("nan"))
        mlp_acc = d_res.get(f"MLP-d{d}", {}).get("best_val", float("nan"))
        oracle_acc = d_res.get(f"GML-Fixed-d{d}", {}).get("best_val", float("nan"))
        print(f"  {d:>4}  {'  delta(TK - Oracle)':>22}  {alg_acc - oracle_acc:>+10.4f}")
        print(f"  {d:>4}  {'  delta(TK - MLP)':>22}  {alg_acc - mlp_acc:>+10.4f}")
        print()

    # Save + plot
    save_dict = {"d_values": np.array(PART_D_D_VALUES)}
    for d, d_res in scaling_results.items():
        for name, h in d_res.items():
            key = name.replace("-", "_")
            save_dict[key + "_val_acc"] = np.array(h["val_acc"])
            save_dict[key + "_T_K"]     = np.array(h["T_K"])
    np.savez(os.path.join(OUT_DIR, "exp5_part_d.npz"), **save_dict)
    print(f"  Saved -> {OUT_DIR}/exp5_part_d.npz")

    _plot_part_d(scaling_results, out_path=os.path.join(OUT_DIR, "exp5_part_d.png"))
    return scaling_results


def _plot_part_d(scaling_results: dict, out_path: str) -> None:
    d_list = sorted(scaling_results.keys())
    n_d    = len(d_list)
    cmap   = {"GML": "tab:green", "Alg": "tab:orange", "MLP": "tab:blue"}
    ls_map = {"GML": "-.",        "Alg": "-",           "MLP": "--"}

    fig, axes = plt.subplots(2, n_d, figsize=(4.5 * n_d, 9))
    if n_d == 1:
        axes = axes.reshape(2, 1)
    fig.suptitle(
        f"Exp 5D -- Extended {PART_D_N_EPOCHS} epochs, d in {d_list}",
        fontsize=12,
    )

    for col, d in enumerate(d_list):
        d_res  = scaling_results[d]
        n_ep   = len(next(iter(d_res.values()))["val_acc"])
        epochs = list(range(1, n_ep + 1))

        # accuracy curves
        ax = axes[0, col]
        for name, h in d_res.items():
            prefix = "GML" if "GML" in name else ("Alg" if "Alg" in name else "MLP")
            ax.plot(epochs, h["val_acc"],
                    color=cmap[prefix], ls=ls_map[prefix], lw=1.8, label=prefix)
        ax.set_title(f"su({d})  n_gen={d*d-1}", fontsize=10)
        ax.set_xlabel("Epoch"); ax.set_ylabel("Val Acc")
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

        # T_K curves
        ax = axes[1, col]
        for name, h in d_res.items():
            vals = h["T_K"]
            if not all(np.isnan(v) for v in vals):
                prefix = "GML" if "GML" in name else ("Alg" if "Alg" in name else "MLP")
                y_plot = [max(v, 1e-5) if not np.isnan(v) else float("nan") for v in vals]
                ax.semilogy(epochs, y_plot,
                            color=cmap[prefix], ls=ls_map[prefix], lw=1.8, label=prefix)
        ax.set_title(f"T_K  su({d})", fontsize=10)
        ax.set_xlabel("Epoch"); ax.set_ylabel("T_K (log)")
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"\n  Part D figure -> {out_path}")


# ============================================================
#  PART E -- normalized lambda_K (lambda_K / n_gen) across d
# ============================================================

PART_E_D_VALUES  = [3, 4, 5, 6, 8]
PART_E_N_EPOCHS  = 1000
PART_E_NGEN_REF  = 8   # reference n_gen for d=3; lambda_K is calibrated here


def lambda_k_for_d(d: int) -> float:
    """Scale Killing pressure so effective gradient magnitude is constant across d.

    Motivation: killing_tension() sums over n_gen^2 elements, so the raw
    gradient ~O(n_gen * lambda_K).  Normalizing by the d=3 reference keeps
    per-generator pressure identical regardless of algebra size.
    """
    n_gen = d * d - 1
    return LAMBDA_K * PART_E_NGEN_REF / n_gen


def run_part_e():
    n_gen_ref = PART_E_NGEN_REF
    print("\n" + "=" * 70)
    print(f"  Part E: Normalized lambda_K = {LAMBDA_K} * {n_gen_ref} / n_gen")
    print(f"  d values: {PART_E_D_VALUES}   epochs: {PART_E_N_EPOCHS}")
    print("=" * 70)

    scaling_results = {}

    for d in PART_E_D_VALUES:
        n_gen   = d * d - 1
        epochs  = PART_E_N_EPOCHS
        lk_eff  = lambda_k_for_d(d)

        print(f"\n  --- d={d}  (su({d}), n_gen={n_gen}, lk_eff={lk_eff:.4f}, epochs={epochs}) ---")
        tr, va = make_loaders(d, seed=SEED)

        models_d = {
            f"GML-Fixed-d{d}" : AlgebraicNet(d=d, n_gen_active=n_gen, fixed=True,  lambda_k=0.0),
            f"AlgNet-TK-d{d}" : AlgebraicNet(d=d, n_gen_active=n_gen, fixed=False, lambda_k=lk_eff),
            f"MLP-d{d}"       : StandardMLP(d=d),
        }

        d_res = {}
        for name, model in models_d.items():
            torch.cuda.reset_peak_memory_stats() if DEVICE.type == "cuda" else None
            t0 = time.perf_counter()
            h  = train(
                model, tr, va,
                name=name, n_epochs=epochs,
                quiet=False, report_every=max(100, epochs // 5),
            )
            h["time_s"]      = time.perf_counter() - t0
            h["peak_mem_mb"] = (
                torch.cuda.max_memory_allocated() / 1e6 if DEVICE.type == "cuda" else 0.0
            )

            if hasattr(model, "get_generators") and not getattr(
                getattr(model, "alg_layer", None), "fixed", False
            ):
                alg_layer = model.alg_layer
                theta_np  = alg_layer._theta().detach().cpu().numpy()
                h["ortho_error"] = float(
                    np.linalg.norm(theta_np @ theta_np.T - np.eye(n_gen))
                )
            else:
                h["ortho_error"] = float("nan")

            d_res[name] = h

        scaling_results[d] = d_res

        print(f"\n  d={d}  results (lambda_K_eff={lk_eff:.4f}):")
        for name, h in d_res.items():
            tk = h["T_K"][-1]
            tk_str = f"{tk:.4f}" if not np.isnan(tk) else "nan"
            print(f"    {name:<22}  best={h['best_val']:.4f}  "
                  f"T_K={tk_str}  mem={h['peak_mem_mb']:.1f}MB  "
                  f"time={h['time_s']:.1f}s")

    # Summary table
    print("\n\n  PART E SCALING SUMMARY (norm. lambda_K, 1000 epochs)")
    print(f"  {'d':>4}  {'Model':<22}  {'Best acc':>10}  {'T_K final':>12}  "
          f"{'lk_eff':>8}  {'OrthoErr':>10}  {'Mem MB':>8}")
    print("  " + "-" * 82)
    for d, d_res in scaling_results.items():
        n_gen  = d * d - 1
        lk_eff = lambda_k_for_d(d)
        for name, h in d_res.items():
            tk = h["T_K"][-1]
            tk_str = f"{tk:12.4f}" if not np.isnan(tk) else f"{'nan':>12}"
            print(
                f"  {d:>4}  {name:<22}  {h['best_val']:>10.4f}  "
                f"{tk_str}  {lk_eff:>8.4f}  "
                f"{h.get('ortho_error', float('nan')):>10.4f}  "
                f"{h['peak_mem_mb']:>8.1f}"
            )
        alg_acc    = d_res.get(f"AlgNet-TK-d{d}", {}).get("best_val", float("nan"))
        mlp_acc    = d_res.get(f"MLP-d{d}",       {}).get("best_val", float("nan"))
        oracle_acc = d_res.get(f"GML-Fixed-d{d}", {}).get("best_val", float("nan"))
        print(f"  {d:>4}  {'  delta(TK - Oracle)':>22}  {alg_acc - oracle_acc:>+10.4f}")
        print(f"  {d:>4}  {'  delta(TK - MLP)':>22}  {alg_acc - mlp_acc:>+10.4f}")
        print()

    # Save
    save_dict = {"d_values": np.array(PART_E_D_VALUES)}
    for d, d_res in scaling_results.items():
        for name, h in d_res.items():
            key = name.replace("-", "_")
            save_dict[key + "_val_acc"] = np.array(h["val_acc"])
            save_dict[key + "_T_K"]     = np.array(h["T_K"])
    np.savez(os.path.join(OUT_DIR, "exp5_part_e.npz"), **save_dict)
    print(f"  Saved -> {OUT_DIR}/exp5_part_e.npz")

    _plot_part_e(scaling_results, out_path=os.path.join(OUT_DIR, "exp5_part_e.png"))
    return scaling_results


def _plot_part_e(scaling_results: dict, out_path: str) -> None:
    d_list = sorted(scaling_results.keys())
    n_d    = len(d_list)
    cmap   = {"GML": "tab:green", "Alg": "tab:orange", "MLP": "tab:blue"}
    ls_map = {"GML": "-.",        "Alg": "-",           "MLP": "--"}

    # 3 rows: accuracy, T_K log, delta bars
    fig, axes = plt.subplots(3, n_d, figsize=(4.5 * n_d, 12))
    if n_d == 1:
        axes = axes.reshape(3, 1)
    fig.suptitle(
        f"Exp 5E -- Normalized lambda_K, {PART_E_N_EPOCHS} epochs",
        fontsize=12,
    )

    deltas_alg_mlp  = []
    d_labels        = []

    for col, d in enumerate(d_list):
        d_res  = scaling_results[d]
        n_ep   = len(next(iter(d_res.values()))["val_acc"])
        epochs = list(range(1, n_ep + 1))
        lk_eff = lambda_k_for_d(d)

        # Row 0: accuracy
        ax = axes[0, col]
        for name, h in d_res.items():
            prefix = "GML" if "GML" in name else ("Alg" if "Alg" in name else "MLP")
            ax.plot(epochs, h["val_acc"],
                    color=cmap[prefix], ls=ls_map[prefix], lw=1.8, label=prefix)
        ax.set_title(f"su({d})  n_gen={d*d-1}  lk={lk_eff:.3f}", fontsize=9)
        ax.set_xlabel("Epoch"); ax.set_ylabel("Val Acc")
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

        # Row 1: T_K log
        ax = axes[1, col]
        for name, h in d_res.items():
            vals = h["T_K"]
            if not all(np.isnan(v) for v in vals):
                prefix = "GML" if "GML" in name else ("Alg" if "Alg" in name else "MLP")
                y_plot = [max(v, 1e-5) if not np.isnan(v) else float("nan") for v in vals]
                ax.semilogy(epochs, y_plot,
                            color=cmap[prefix], ls=ls_map[prefix], lw=1.8, label=prefix)
        ax.set_title(f"T_K  su({d})", fontsize=9)
        ax.set_xlabel("Epoch"); ax.set_ylabel("T_K (log)")
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

        # Collect delta for bar chart
        alg_acc = d_res.get(f"AlgNet-TK-d{d}", {}).get("best_val", float("nan"))
        mlp_acc = d_res.get(f"MLP-d{d}",       {}).get("best_val", float("nan"))
        deltas_alg_mlp.append(alg_acc - mlp_acc)
        d_labels.append(f"d={d}")

        # Row 2: bar comparing Alg vs oracle vs MLP
        ax = axes[2, col]
        names_bar = ["GML-Fixed", "AlgNet-TK", "MLP"]
        accs_bar  = [
            d_res.get(f"GML-Fixed-d{d}", {}).get("best_val", 0.0),
            d_res.get(f"AlgNet-TK-d{d}", {}).get("best_val", 0.0),
            d_res.get(f"MLP-d{d}",       {}).get("best_val", 0.0),
        ]
        colors_bar = [cmap["GML"], cmap["Alg"], cmap["MLP"]]
        bars = ax.bar(names_bar, accs_bar, color=colors_bar, alpha=0.8)
        ax.set_ylim(min(accs_bar) - 0.02, max(accs_bar) + 0.02)
        for bar, acc in zip(bars, accs_bar):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.001,
                    f"{acc:.4f}", ha="center", va="bottom", fontsize=7)
        ax.set_title(f"Best acc  su({d})", fontsize=9)
        ax.set_ylabel("Val Acc"); ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"\n  Part E figure -> {out_path}")

    # Print delta summary
    print(f"\n  Delta(AlgNet-TK - MLP) with normalized lambda_K:")
    for d, delta in zip(d_list, deltas_alg_mlp):
        bar = "#" * int(abs(delta) * 200)
        sign = "+" if delta >= 0 else ""
        print(f"    d={d}: {sign}{delta:+.4f}  {bar}")


# ============================================================
#  ENTRY POINT
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Exp5: scaling study")
    parser.add_argument(
        "--part", choices=["a", "b", "c", "d", "e", "all"], default="all",
        help=("Which part to run: 'a' (n_gen sweep), 'b' (d-sweep), "
              "'c' (lambda_gate sweep), 'd' (extended epochs d>=4), "
              "'e' (normalized lambda_K), or 'all'"),
    )
    parser.add_argument(
        "--d-values", nargs="+", type=int, default=None,
        help="Override d values for Part B, e.g. --d-values 2 3 4 5",
    )
    parser.add_argument(
        "--epochs-a", type=int, default=None,
        help="Epochs for Part A (default: use script constant)",
    )
    args = parser.parse_args()

    global D_VALUES, PART_A_N_EPOCHS
    if args.d_values is not None:
        D_VALUES = sorted(args.d_values)
    if args.epochs_a is not None:
        PART_A_N_EPOCHS = args.epochs_a

    print(f"\nDevice      : {DEVICE}")
    if DEVICE.type == "cuda":
        print(f"GPU         : {torch.cuda.get_device_name(DEVICE)}")
        print(f"VRAM total  : {torch.cuda.get_device_properties(DEVICE).total_memory / 1e9:.1f} GB")
    print(f"Output dir  : {OUT_DIR}\n")

    if args.part in ("a", "all"):
        run_part_a()
    if args.part in ("b", "all"):
        run_part_b()
    if args.part in ("c", "all"):
        run_part_c()
    if args.part in ("d", "all"):
        run_part_d()
    if args.part in ("e", "all"):
        run_part_e()

    print("\n  All done.")


if __name__ == "__main__":
    main()
