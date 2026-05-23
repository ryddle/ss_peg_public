"""
exp6bc_timing_pressure.py
=========================
Mini-Exp6b  —  Phase-timing ablation for d=5
    Hypothesis: Phase1=800ep lets the task head mature (val>0.58) before the
    orthogonal-freeze in Phase2, fixing the d=5 OrthoAnn collapse seen in Exp6.
    Protocol variants compared (all 1000 total epochs):
        OrthoAnn-orig  : Ph1=600 + Ph2=200 + Ph3=200   (Exp6 baseline)
        OrthoAnn-6b    : Ph1=800 + Ph2=100 + Ph3=100   (new timing)
        Ortho1         : single phase, lambda_o=1.0
        MLP            : no algebraic structure

Mini-Exp6c  —  Killing-pressure doubling for d=5
    Hypothesis: lambda_K_eff*2 forces tighter Lie-algebra confinement,
    pushing T_K from ~0.0012 below the crystallisation threshold 0.0010.
    Protocols compared (1000 epochs, single phase):
        Ortho1-base    : lambda_K = lambda_K_eff(5) = 0.6667, lambda_o=1.0
        Ortho1-2x      : lambda_K = 2 * lambda_K_eff(5) = 1.3333, lambda_o=1.0
        MLP            : baseline (shared with Exp6b run)

Physical analogy
----------------
Phase timing ~ confinement timescale: in QCD, colour charges must first
reach a bound state (task-relevant representation) before the confining
flux tube (T_ortho pressure) is fully applied.  Applying confinement too
early on a disordered state binds random colour configurations.
"""

import os
import sys
import time

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ─── add parent dir so killing_utils is importable ───────────────────────────
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from training_ai_test.killing_utils import (
    killing_tension,
    params_to_su_generators,
)

# ─── output directory ─────────────────────────────────────────────────────────
OUT_DIR = os.path.join(os.path.dirname(__file__), "exp6bc_outputs")
os.makedirs(OUT_DIR, exist_ok=True)

# ─── device ──────────────────────────────────────────────────────────────────
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# ─── fixed experiment dimension ───────────────────────────────────────────────
D = 5
N_GEN = D * D - 1        # 24 for su(5)

# ─── data config ─────────────────────────────────────────────────────────────
N_TRAIN        = 6_000
N_VAL          = 1_500
BATCH_SIZE     = 256
DATA_SEED      = 42        # same as Exp6 for comparability

# ─── training config ──────────────────────────────────────────────────────────
LR_MAIN        = 1e-3
LR_ANNEAL      = 1e-4

# Killing-tension normalisation (same formula as Exp6)
LAMBDA_K_BASE  = 2.0
N_GEN_REF      = 8
def lambda_k_eff(d: int) -> float:
    return LAMBDA_K_BASE * N_GEN_REF / (d * d - 1)

LK5            = lambda_k_eff(D)       # 16/24 ≈ 0.6667
LK5_2X         = 2.0 * LK5            # ≈ 1.3333

SINGLE_EPOCHS  = 1_000

# ─── Mini-Exp6b phase configs ─────────────────────────────────────────────────
# PH1_LO / PH2_LO / PH3_LO: lambda_ortho in each phase
PH1_LO = 0.1
PH2_LO = 2.0
PH3_LO = 0.5

CONFIGS_6B = {
    "OrthoAnn-orig": {"ph1": 600, "ph2": 200, "ph3": 200},   # Exp6 baseline
    "OrthoAnn-6b"  : {"ph1": 800, "ph2": 100, "ph3": 100},   # new timing
}

REPORT_EVERY = 200


# ============================================================
#  Dataset
# ============================================================

def generate_su_d_cubic_dataset(d: int, n: int, seed: int = 0):
    """
    Generate n samples of SU(d) cubic-invariant classification.
    Features: real Killing-metric coordinates (n_gen dims).
    Label: sign(Re(Tr(M^3))) where M = sum_a x_a * T_a in su(d).
    """
    from training_ai_test.killing_utils import _su_basis
    rng   = np.random.default_rng(seed)
    n_gen = d * d - 1
    basis = _su_basis(d, "cpu").numpy()    # (n_gen, d, d) complex128

    X   = rng.standard_normal((n, n_gen)).astype(np.float64)
    M_c = np.einsum("na,aij->nij", X, basis)     # (n, d, d) complex
    tr3 = np.trace(M_c @ M_c @ M_c, axis1=-2, axis2=-1).real
    y   = (tr3 > 0).astype(np.int64)

    return torch.tensor(X, dtype=torch.float64), torch.tensor(y, dtype=torch.long)


def make_loaders(d: int, seed: int = DATA_SEED):
    X_tr, y_tr = generate_su_d_cubic_dataset(d, N_TRAIN, seed=seed)
    X_va, y_va = generate_su_d_cubic_dataset(d, N_VAL,   seed=seed + 1000)
    tr = DataLoader(TensorDataset(X_tr, y_tr), batch_size=BATCH_SIZE, shuffle=True)
    va = DataLoader(TensorDataset(X_va, y_va), batch_size=BATCH_SIZE)
    return tr, va


# ============================================================
#  Models
# ============================================================

def make_head(n_gen: int) -> nn.Sequential:
    return nn.Sequential(
        nn.ReLU(),
        nn.Linear(n_gen, 4 * n_gen),
        nn.ReLU(),
        nn.Linear(4 * n_gen, 2),
    )


class AlgebraicProjection(nn.Module):
    """Learnable linear projection onto su(d) coordinates via generators θ."""

    def __init__(self, d: int, fixed: bool = False):
        super().__init__()
        self.d     = d
        self.n_gen = d * d - 1
        theta_init = torch.eye(self.n_gen, dtype=torch.float64)
        self.theta = nn.Parameter(theta_init, requires_grad=not fixed)

    def get_generators(self):
        return params_to_su_generators(self.theta, self.d)

    def forward(self, x):
        return x @ self.theta.T


class AlgebraicNet(nn.Module):
    def __init__(self, d: int, fixed: bool = False):
        super().__init__()
        self.proj = AlgebraicProjection(d, fixed)
        self.head = make_head(d * d - 1)

    def forward(self, x):
        return self.head(self.proj(x))


class StandardMLP(nn.Module):
    def __init__(self, d: int):
        super().__init__()
        n = d * d - 1
        self.net = nn.Sequential(
            nn.Linear(n, n, dtype=torch.float64),
            *make_head(n),
        )

    def forward(self, x):
        return self.net(x)


# ============================================================
#  Evaluation & diagnostics
# ============================================================

@torch.no_grad()
def _evaluate(model, val_loader) -> float:
    model.eval()
    correct = total = 0
    for xb, yb in val_loader:
        xb, yb = xb.to(DEVICE), yb.to(DEVICE)
        pred = model(xb).argmax(dim=1)
        correct += (pred == yb).sum().item()
        total   += len(yb)
    return correct / total


@torch.no_grad()
def _algebraic_diagnostics(model):
    """Returns (T_K, T_ortho, OrthoErr)  — float('nan') for non-algebraic models."""
    if not hasattr(model, "proj"):
        return (float("nan"), float("nan"), float("nan"))
    theta = model.proj.theta
    gens  = model.proj.get_generators()
    T_K   = killing_tension(gens).item()
    eye   = torch.eye(theta.shape[0], dtype=theta.dtype, device=theta.device)
    diff  = theta @ theta.T - eye
    T_o   = (diff ** 2).sum().item()
    O_err = diff.norm(p="fro").item()
    return (T_K, T_o, O_err)


# ============================================================
#  Training primitives
# ============================================================

def train_phase(
    model,
    train_loader,
    val_loader,
    *,
    n_epochs: int,
    lambda_k: float,
    lambda_ortho: float,
    lr: float,
    freeze_head: bool,
    phase_id: str,
    warmup: int = 30,
    report_every: int = REPORT_EVERY,
    quiet: bool = False,
) -> dict:
    """Train one phase; returns history dict."""
    model.train()
    if freeze_head and hasattr(model, "head"):
        for p in model.head.parameters():
            p.requires_grad_(False)
    else:
        if hasattr(model, "head"):
            for p in model.head.parameters():
                p.requires_grad_(True)

    opt = torch.optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()), lr=lr
    )
    ce  = nn.CrossEntropyLoss()

    hist = {"val_acc": [], "T_K": [], "T_ortho": [], "ortho_err": [], "best_val": 0.0}

    for ep in range(1, n_epochs + 1):
        model.train()
        for xb, yb in train_loader:
            xb, yb = xb.to(DEVICE), yb.to(DEVICE)
            logits = model(xb)
            loss   = ce(logits, yb)

            # Geometric regularisation (algebraic models only)
            if hasattr(model, "proj") and (lambda_k > 0 or lambda_ortho > 0):
                gens = model.proj.get_generators()
                if lambda_k > 0:
                    loss = loss + lambda_k * killing_tension(gens)
                if lambda_ortho > 0:
                    theta = model.proj.theta
                    eye   = torch.eye(theta.shape[0], dtype=theta.dtype, device=theta.device)
                    diff  = theta @ theta.T - eye
                    # warmup ramp for lambda_ortho
                    ramp = min(1.0, ep / warmup) if warmup > 0 else 1.0
                    loss = loss + lambda_ortho * ramp * (diff ** 2).sum()

            opt.zero_grad()
            loss.backward()
            opt.step()

        val  = _evaluate(model, val_loader)
        T_K, T_o, O_err = _algebraic_diagnostics(model)
        hist["val_acc"].append(val)
        hist["T_K"].append(T_K)
        hist["T_ortho"].append(T_o)
        hist["ortho_err"].append(O_err)
        if val > hist["best_val"]:
            hist["best_val"] = val

        if not quiet and ep % report_every == 0:
            print(
                f"    [{phase_id} ep{ep:>4}]  val={val:.4f}  "
                f"T_K={T_K:.5f}  T_o={T_o:.5f}  OErr={O_err:.3f}"
            )

    # Unfreeze head for downstream phases
    if hasattr(model, "head"):
        for p in model.head.parameters():
            p.requires_grad_(True)

    return hist


def _concat_hist(h1: dict, h2: dict, h3: dict | None = None) -> dict:
    keys = ["val_acc", "T_K", "T_ortho", "ortho_err"]
    out  = {k: h1[k] + h2[k] + (h3[k] if h3 else []) for k in keys}
    out["best_val"]    = max(h1["best_val"], h2["best_val"], *([] if h3 is None else [h3["best_val"]]))
    out["total_epochs"] = len(out["val_acc"])
    return out


def train_multiphase(model, tr, va, *, ph1: int, ph2: int, ph3: int, lk: float, name: str) -> dict:
    """
    3-phase annealing:
        Phase1: task learning, mild ortho  (lr=LR_MAIN)
        Phase2: crystallisation, head frozen, strong ortho  (lr=LR_ANNEAL)
        Phase3: fine-tune, ortho relaxed  (lr=LR_MAIN)
    """
    print(f"\n  >> {name} (3-phase {ph1}+{ph2}+{ph3}ep)  lk={lk:.4f}")

    h1 = train_phase(model, tr, va, n_epochs=ph1, lambda_k=lk, lambda_ortho=PH1_LO,
                     lr=LR_MAIN, freeze_head=False, phase_id=f"{name}/Ph1",
                     warmup=30, report_every=ph1 // 4)

    print(f"    -- Phase1 done: best_val={h1['best_val']:.4f}  "
          f"T_K={h1['T_K'][-1]:.5f}  T_o={h1['T_ortho'][-1]:.5f}")

    h2 = train_phase(model, tr, va, n_epochs=ph2, lambda_k=lk, lambda_ortho=PH2_LO,
                     lr=LR_ANNEAL, freeze_head=True, phase_id=f"{name}/Ph2",
                     warmup=10, report_every=max(10, ph2 // 4))

    print(f"    -- Phase2 done: best_val={h2['best_val']:.4f}  "
          f"T_K={h2['T_K'][-1]:.5f}  T_o={h2['T_ortho'][-1]:.5f}")

    h3 = train_phase(model, tr, va, n_epochs=ph3, lambda_k=lk, lambda_ortho=PH3_LO,
                     lr=LR_MAIN, freeze_head=False, phase_id=f"{name}/Ph3",
                     warmup=0, report_every=max(10, ph3 // 4))

    print(f"    -- Phase3 done: best_val={h3['best_val']:.4f}  "
          f"T_K={h3['T_K'][-1]:.5f}  T_o={h3['T_ortho'][-1]:.5f}")

    combined             = _concat_hist(h1, h2, h3)
    combined["phase_boundaries"] = [ph1, ph1 + ph2]
    return combined


def train_single(model, tr, va, *, lk: float, lambda_ortho: float,
                 n_epochs: int, name: str) -> dict:
    """Single-phase training."""
    print(f"\n  >> {name} (single {n_epochs}ep  lk={lk:.4f}  lo={lambda_ortho:.4f})")
    h = train_phase(model, tr, va, n_epochs=n_epochs, lambda_k=lk,
                    lambda_ortho=lambda_ortho, lr=LR_MAIN,
                    freeze_head=False, phase_id=name,
                    warmup=30, report_every=REPORT_EVERY)
    h["total_epochs"]    = n_epochs
    h["phase_boundaries"] = []
    return h


# ============================================================
#  Experiment runners
# ============================================================

def run_mini_6b(tr, va) -> dict:
    """
    Mini-Exp6b — Phase-timing ablation (d=5).
    Models: MLP, Ortho1, OrthoAnn-orig (Exp6), OrthoAnn-6b (delayed anneal).
    """
    print(f"\n{'='*60}")
    print(f"  MINI-EXP 6b — Phase timing ablation  (d={D})")
    print(f"{'='*60}")
    results = {}

    # MLP baseline
    model = StandardMLP(D).double().to(DEVICE)
    results["MLP"] = train_single(model, tr, va, lk=0.0, lambda_ortho=0.0,
                                  n_epochs=SINGLE_EPOCHS, name="MLP")

    # Ortho1: single phase, lambda_o=1.0
    model = AlgebraicNet(D, fixed=False).double().to(DEVICE)
    results["Ortho1"] = train_single(model, tr, va, lk=LK5, lambda_ortho=1.0,
                                     n_epochs=SINGLE_EPOCHS, name="Ortho1")

    # OrthoAnn-orig: Exp6 baseline timings
    cfg = CONFIGS_6B["OrthoAnn-orig"]
    model = AlgebraicNet(D, fixed=False).double().to(DEVICE)
    results["OrthoAnn-orig"] = train_multiphase(
        model, tr, va, ph1=cfg["ph1"], ph2=cfg["ph2"], ph3=cfg["ph3"],
        lk=LK5, name="OrthoAnn-orig",
    )

    # OrthoAnn-6b: delayed crystallisation
    cfg = CONFIGS_6B["OrthoAnn-6b"]
    model = AlgebraicNet(D, fixed=False).double().to(DEVICE)
    results["OrthoAnn-6b"] = train_multiphase(
        model, tr, va, ph1=cfg["ph1"], ph2=cfg["ph2"], ph3=cfg["ph3"],
        lk=LK5, name="OrthoAnn-6b",
    )

    _print_summary(results, "Mini-Exp6b")
    return results


def run_mini_6c(tr, va, mlp_hist: dict, ortho1_hist: dict) -> dict:
    """
    Mini-Exp6c — Killing-pressure doubling (d=5).
    MLP and Ortho1-base reused from Exp6b run.
    New model: Ortho1-2x with lk=2*lambda_K_eff(5).
    """
    print(f"\n{'='*60}")
    print(f"  MINI-EXP 6c — Killing pressure doubling  (d={D})")
    print(f"{'='*60}")
    results = {}

    # Reuse shared runs from 6b
    results["MLP"]        = mlp_hist
    results["Ortho1-base"] = ortho1_hist   # same as Ortho1 from 6b

    # Ortho1-2x: doubled lambda_K
    model = AlgebraicNet(D, fixed=False).double().to(DEVICE)
    results["Ortho1-2x"] = train_single(model, tr, va, lk=LK5_2X, lambda_ortho=1.0,
                                        n_epochs=SINGLE_EPOCHS, name="Ortho1-2x")

    _print_summary(results, "Mini-Exp6c")
    return results


def _print_summary(results: dict, label: str):
    mlp_best = results.get("MLP", {}).get("best_val", float("nan"))
    print(f"\n  {label} SUMMARY  (d={D})")
    header = f"  {'Model':<20}  {'Best':>8}  {'Final':>8}  {'T_K':>9}  {'T_o':>9}  {'OrthoErr':>9}  {'Ep':>5}"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for name, h in results.items():
        best   = h["best_val"]
        final  = h["val_acc"][-1]
        T_K    = h["T_K"][-1]
        T_o    = h["T_ortho"][-1]
        O_err  = h["ortho_err"][-1]
        ep     = h.get("total_epochs", len(h["val_acc"]))
        delta  = best - mlp_best

        def _fmt(v, w, d):
            return f"{v:{w}.{d}f}" if not (isinstance(v, float) and np.isnan(v)) else f"{'--':>{w}}"

        print(
            f"  {name:<20}  {best:>8.4f}  {final:>8.4f}  "
            f"{_fmt(T_K,9,5)}  {_fmt(T_o,9,5)}  {_fmt(O_err,9,3)}  {ep:>5}"
            f"  [Δ_MLP={delta:+.4f}]"
        )


# ============================================================
#  Plotting
# ============================================================

_COLOR = {
    "MLP"          : ("tab:blue",   "-",  1.5),
    "Ortho1"       : ("tab:red",    "-",  1.8),
    "Ortho1-base"  : ("tab:red",    "-",  1.8),
    "Ortho1-2x"    : ("tab:orange", "-",  2.0),
    "OrthoAnn-orig": ("tab:purple", "--", 1.8),
    "OrthoAnn-6b"  : ("tab:green",  "-",  2.2),
}


def _plot_panel(ax_acc, ax_tk, ax_oe, results: dict, title: str):
    for name, h in results.items():
        c, ls, lw = _COLOR.get(name, ("black", "-", 1.5))
        eps = list(range(1, len(h["val_acc"]) + 1))
        ax_acc.plot(eps, h["val_acc"], color=c, ls=ls, lw=lw, label=name)

        tk_vals = [max(v, 1e-7) if not np.isnan(v) else np.nan for v in h["T_K"]]
        valid = [(e, v) for e, v in zip(eps, tk_vals) if not np.isnan(v)]
        if valid:
            ev, vv = zip(*valid)
            ax_tk.semilogy(ev, vv, color=c, ls=ls, lw=lw, label=name)

        oe_vals = [v if not np.isnan(v) else np.nan for v in h["ortho_err"]]
        valid   = [(e, v) for e, v in zip(eps, oe_vals) if not np.isnan(v)]
        if valid:
            ev, vv = zip(*valid)
            ax_oe.plot(ev, vv, color=c, ls=ls, lw=lw, label=name)

    # Phase boundaries (only annealing models have them)
    for name, h in results.items():
        for bnd in h.get("phase_boundaries", []):
            for ax in (ax_acc, ax_tk, ax_oe):
                ax.axvline(bnd, color="gray", ls=":", lw=1.0, alpha=0.5)

    ax_acc.set_title(f"{title} — Val Accuracy")
    ax_tk.set_title(f"{title} — T_K (log)")
    ax_oe.set_title(f"{title} — OrthoErr")
    for ax in (ax_acc, ax_tk, ax_oe):
        ax.set_xlabel("Epoch")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    ax_acc.set_ylabel("Accuracy")
    ax_tk.set_ylabel("T_K (log)")
    ax_oe.set_ylabel("OrthoErr (Frobenius)")

    # crystallisation threshold line on T_K panel
    ax_tk.axhline(1e-3, color="red", ls=":", lw=1.0, alpha=0.4, label="T_K=1e-3")


def plot_all(results_6b: dict, results_6c: dict):
    fig, axes = plt.subplots(2, 3, figsize=(17, 9))
    _plot_panel(*axes[0], results_6b,
                title=f"Mini-Exp6b d={D}  (OrthoAnn timing)")
    _plot_panel(*axes[1], results_6c,
                title=f"Mini-Exp6c d={D}  (λ_K pressure)")
    fig.suptitle(
        f"Exp6b+6c — Phase-timing & Killing-pressure ablation  (su({D}), n_gen={N_GEN})\n"
        f"λ_K_eff={LK5:.4f}  λ_K_2x={LK5_2X:.4f}  |  "
        f"Gray dotted = phase boundaries",
        fontsize=10,
    )
    plt.tight_layout()
    out = os.path.join(OUT_DIR, "exp6bc_results.png")
    plt.savefig(out, dpi=130, bbox_inches="tight")
    plt.close()
    print(f"\n  Figure saved -> {out}")


# ============================================================
#  Main
# ============================================================

def main():
    print(f"\n{'#'*65}")
    print(f"  Mini-Exp 6b & 6c  — Phase timing & Killing-pressure ablation")
    print(f"  Device : {DEVICE}")
    if DEVICE.type == "cuda":
        print(f"  GPU    : {torch.cuda.get_device_name(0)}")
    print(f"  d={D}  n_gen={N_GEN}  λ_K_eff={LK5:.4f}  λ_K_2x={LK5_2X:.4f}")
    print(f"  Epochs : {SINGLE_EPOCHS} (single) / {max(c['ph1']+c['ph2']+c['ph3'] for c in CONFIGS_6B.values())} (multiphase)")
    print(f"{'#'*65}\n")

    tr, va = make_loaders(D)

    t0 = time.perf_counter()

    # ── Mini-Exp6b ──
    results_6b = run_mini_6b(tr, va)

    # ── Mini-Exp6c: reuse MLP and Ortho1 from 6b ──
    results_6c = run_mini_6c(tr, va,
                              mlp_hist=results_6b["MLP"],
                              ortho1_hist=results_6b["Ortho1"])

    elapsed = time.perf_counter() - t0

    # ── Cross-mini summary ──
    print(f"\n\n{'='*65}")
    print(f"  EXP6B+6C  CROSS-EXPERIMENT SUMMARY  (d={D})")
    print(f"{'='*65}")
    mlp_best = results_6b["MLP"]["best_val"]
    print(f"  MLP baseline: {mlp_best:.4f}")

    for exp_label, results in [("6b", results_6b), ("6c", results_6c)]:
        print(f"\n  -- {exp_label} --")
        for name, h in results.items():
            if name == "MLP":
                continue
            print(f"    {name:<20}  best={h['best_val']:.4f}  "
                  f"T_K_final={h['T_K'][-1]:.5f}  "
                  f"OrthoErr_final={h['ortho_err'][-1]:.3f}  "
                  f"Δ_MLP={h['best_val']-mlp_best:+.4f}")

    print(f"\n  Total wall time: {elapsed:.1f}s")

    # ── Plots ──
    plot_all(results_6b, results_6c)

    # ── Save numerics ──
    save_dict = {}
    for exp_label, results in [("6b", results_6b), ("6c", results_6c)]:
        for name, h in results.items():
            key = f"{exp_label}_{name.replace('-', '_')}"
            for field in ("val_acc", "T_K", "T_ortho", "ortho_err"):
                save_dict[f"{key}_{field}"] = np.array(h[field])
            save_dict[f"{key}_best"] = np.array([h["best_val"]])

    npz = os.path.join(OUT_DIR, "exp6bc_results.npz")
    np.savez(npz, **save_dict)
    print(f"  Numerics saved  -> {npz}")


if __name__ == "__main__":
    main()
