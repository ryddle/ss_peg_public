"""
Exp13: Lyapunov Exponents of Adam Dynamics in Parameter Space
=============================================================

Question:
  Is the Adam optimization trajectory in θ-space deterministically chaotic?
  Positive maximum Lyapunov exponent λ_max → chaos → fractal attractor
  with Kaplan-Yorke dimension d_KY = j + (λ₁+...+λ_j) / |λ_{j+1}|.

Context (§11g of EXPERIMENTAL_SUMMARY.md):
  In the F₄ algebra discovery experiment, a Lyapunov bifurcation was found at
  step ~1200 with λ ≈ 0.061/step — ε~10⁻¹³ CUDA perturbations were amplified
  by 13 orders of magnitude in ~500 steps. This experiment measures the same
  phenomenon in the neural-network (AlgNet on su(3)) context.

Design:
  Task      : sign(Re(Tr(M^3))),  M = Σ_a x_a T_a ∈ su(3)
              N_GEN = 8, N_DIM = 3, complex anti-Hermitian generators
  n_train   : 6 000,  n_val : 1 500
  N_REPS    : 5 seeds
  Epochs    : 3 000
  LR        : 1e-3, grad-clip max_norm=1.0
  Models    : AlgNet-1L (θ ∈ ℝ^64) and MLP-1L (W₁ ∈ ℝ^64)

Lyapunov measurement — twin-trajectory method:
  At each PROBE_INTERVAL epochs, inject a small parameter perturbation
  δθ (Gaussian, normalised to ε₀) into a COPY of the model with an
  INDEPENDENT copy of the Adam state (m, v, step).
  Run both primary and perturbed trajectories for DIVERGE_STEPS more steps.
  Measure finite-time Lyapunov exponent:
      λ(t) = (1 / DIVERGE_STEPS) × log( ||Δθ_final|| / ||Δθ_initial|| )
  where ||Δθ|| is the L2 norm of the parameter difference.

  This captures the maximum Lyapunov exponent λ_max of the extended-phase-space
  map F: (θ, m, v) → (θ', m', v') induced by one Adam update step.

Extended phase space connection (§11g):
  Experiments use both:
  (a) θ-only norm: ||Δθ||                     (parameter space only)
  (b) full-state norm: ||(Δθ, Δm, Δv)||       (extended phase space, z as in §11g)
  The §11g finding: F₄ basin exists in (θ, m, v) space — the Adam momentum
  state is the causal factor. This experiment tests whether the same extended-
  space sensitivity appears in the AlgNet algebraic classification task.

Outputs per probe:
  λ_theta(t)     — Lyapunov in θ-space alone
  λ_extended(t)  — Lyapunov in (θ,m,v) extended phase space
  ||Δθ||_t       — divergence curve (primary + perturbed parameter distance)

Summary:
  λ(t) vs training epoch for AlgNet vs MLP
  Bifurcation onset (epoch where λ transitions positive → negative or vice versa)
  Kaplan-Yorke dimension d_KY (estimated from λ(t) near convergence)

Outputs:
  exp13_outputs/exp13_results.png  — λ(t) traces + bifurcation analysis
  exp13_outputs/exp13_results.npz  — full λ arrays
"""

# ---------------------------------------------------------------------------
#  Paths and imports
# ---------------------------------------------------------------------------
import os
import sys
import copy

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import time
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

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
N_GEN           = N_DIM ** 2 - 1

LAMBDA_K        = 2.0
LAMBDA_ORTHO    = 1.0
WARMUP_EP       = 30

BATCH_SIZE      = 256
LR_MAIN         = 1e-3
N_EPOCHS        = 3_000
REPORT_EVERY    = 500

N_TRAIN         = 6_000
N_VAL           = 1_500
N_REPS          = 5
BASE_SEEDS      = [42, 43, 44, 45, 46]

# Lyapunov probe settings
PROBE_INTERVAL  = 100      # inject perturbation every this many epochs
DIVERGE_STEPS   = 50       # steps to run each clone pair (in epochs)
EPSILON_0       = 1e-6     # initial perturbation magnitude (⟨||δθ||⟩ = ε₀)

OUT_DIR = os.path.join(HERE, "exp13_outputs")
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
        return (torch.einsum("aij,...ji->...a", G, M_x)).real.to(torch.float64) + self.bias

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


class MLP1L(nn.Module):
    def __init__(self):
        super().__init__()
        self.fc1  = nn.Linear(N_GEN, N_GEN, dtype=torch.float64)
        self.relu = nn.ReLU()
        self.head = nn.Linear(N_GEN, 2, dtype=torch.float64)

    def forward(self, x, basis=None):
        return self.head(self.relu(self.fc1(x)))


# ---------------------------------------------------------------------------
#  Utilities
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


def get_param_vector(model: nn.Module) -> torch.Tensor:
    """Return all learnable parameters as a single flat vector."""
    return torch.cat([p.data.ravel() for p in model.parameters() if p.requires_grad])


def set_param_vector(model: nn.Module, vec: torch.Tensor):
    """Set all learnable parameters from flat vector."""
    offset = 0
    for p in model.parameters():
        if p.requires_grad:
            n = p.numel()
            p.data.copy_(vec[offset:offset + n].reshape(p.shape))
            offset += n


def get_adam_state_vectors(optimizer: optim.Adam) -> tuple[torch.Tensor, torch.Tensor]:
    """Return (m, v) as flat vectors from all Adam state groups."""
    m_parts, v_parts = [], []
    for group in optimizer.param_groups:
        for p in group["params"]:
            if p.requires_grad:
                state = optimizer.state.get(p, {})
                if "exp_avg" in state:
                    m_parts.append(state["exp_avg"].data.ravel())
                    v_parts.append(state["exp_avg_sq"].data.ravel())
                else:
                    m_parts.append(torch.zeros(p.numel(), dtype=p.dtype, device=p.device))
                    v_parts.append(torch.zeros(p.numel(), dtype=p.dtype, device=p.device))
    return (torch.cat(m_parts) if m_parts else torch.zeros(1),
            torch.cat(v_parts) if v_parts else torch.zeros(1))


def count_params(model: nn.Module) -> int:
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


# ---------------------------------------------------------------------------
#  One optimisation step (for the clone trajectories)
# ---------------------------------------------------------------------------

def one_step(model, optimizer, data_iter, train_loader, basis,
             criterion, is_algnet: bool, warmup_scale: float):
    """One epoch of training."""
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


# ---------------------------------------------------------------------------
#  Lyapunov probe: measure λ at a given checkpoint
# ---------------------------------------------------------------------------

def probe_lyapunov(
    primary_model:     nn.Module,
    primary_optim:     optim.Adam,
    train_loader,
    basis:             torch.Tensor,
    criterion,
    is_algnet:         bool,
    current_ep:        int,
) -> dict:
    """
    Clone the primary (model + Adam state) and add a small perturbation to θ.
    Run both primary and clone for DIVERGE_STEPS epochs.
    Measure λ in both θ-space and extended phase space (θ, m, v).

    Returns: {"lambda_theta": float, "lambda_extended": float,
              "delta_theta_curve": np.ndarray}
    """
    # Clone primary model and optimiser
    clone_model = copy.deepcopy(primary_model)
    clone_optim = optim.Adam(
        filter(lambda p: p.requires_grad, clone_model.parameters()),
        lr=LR_MAIN
    )
    # Mirror Adam state from primary
    primary_state = primary_optim.state_dict()
    clone_optim.load_state_dict(copy.deepcopy(primary_state))

    # Inject normalised perturbation into clone θ
    rng_probe = np.random.default_rng(current_ep)
    theta_vec  = get_param_vector(clone_model).cpu().numpy()
    noise      = rng_probe.standard_normal(len(theta_vec))
    noise      = noise / (np.linalg.norm(noise) + 1e-30) * EPSILON_0
    perturbed  = theta_vec + noise
    set_param_vector(clone_model, torch.tensor(perturbed, dtype=torch.float64, device=DEVICE))

    # Initial separation
    delta_theta_0 = EPSILON_0  # by construction

    # Also record initial extended-space separation (m and v are shared at t=0, Δm=Δv=0)
    # The extended-space divergence starts purely in θ direction
    delta_extended_0 = EPSILON_0

    delta_theta_curve = [delta_theta_0]
    primary_param_curve  = [get_param_vector(primary_model).cpu().numpy().copy()]
    clone_param_curve    = [get_param_vector(clone_model  ).cpu().numpy().copy()]

    # Run both trajectories forward
    for step in range(DIVERGE_STEPS):
        ep_eff = current_ep + step
        ws     = min(1.0, ep_eff / WARMUP_EP) if is_algnet else 0.0

        one_step(primary_model, primary_optim, None, train_loader,
                 basis, criterion, is_algnet, ws)
        one_step(clone_model,   clone_optim,   None, train_loader,
                 basis, criterion, is_algnet, ws)

        # θ separation
        pv  = get_param_vector(primary_model).cpu().numpy()
        cv  = get_param_vector(clone_model).cpu().numpy()
        dt  = float(np.linalg.norm(pv - cv))
        delta_theta_curve.append(dt)

    delta_theta_final = delta_theta_curve[-1]

    # λ_theta (finite-time, DIVERGE_STEPS)
    with np.errstate(divide="ignore", invalid="ignore"):
        if delta_theta_final > 1e-30 and delta_theta_0 > 1e-30:
            lambda_theta = np.log(delta_theta_final / delta_theta_0) / DIVERGE_STEPS
        else:
            lambda_theta = float("nan")

    # Extended-space: get final (m,v) separation
    m_p, v_p = get_adam_state_vectors(primary_optim)
    m_c, v_c = get_adam_state_vectors(clone_optim)
    dm = (m_p - m_c).cpu().numpy()
    dv = (v_p - v_c).cpu().numpy()
    dt_ext = float(np.linalg.norm(np.concatenate([pv - cv, dm, dv])))
    dt_ext_0 = delta_extended_0  # initial Δm=Δv=0, so same as δθ

    if dt_ext > 1e-30 and dt_ext_0 > 1e-30:
        lambda_extended = np.log(dt_ext / dt_ext_0) / DIVERGE_STEPS
    else:
        lambda_extended = float("nan")

    return {
        "lambda_theta":      lambda_theta,
        "lambda_extended":   lambda_extended,
        "delta_theta_curve": np.array(delta_theta_curve),
    }


# ---------------------------------------------------------------------------
#  Full training run with Lyapunov probes
# ---------------------------------------------------------------------------

def run_with_lyapunov_probes(model, train_loader, val_loader,
                              basis, name: str, is_algnet: bool):
    """
    Train model for N_EPOCHS. Every PROBE_INTERVAL epochs, measure λ.
    Returns dict with λ history and accuracy history.
    """
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=LR_MAIN
    )

    lambda_theta_history    = []
    lambda_extended_history = []
    probe_epochs            = []
    val_accs                = []
    best_val                = 0.0
    delta_curves_sample     = {}  # store delta_theta_curve for a few probes (viz)

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
        val_accs.append(val_acc)
        if val_acc > best_val:
            best_val = val_acc

        # Lyapunov probe
        if ep % PROBE_INTERVAL == 0:
            probe = probe_lyapunov(
                primary_model = model,
                primary_optim = optimizer,
                train_loader  = train_loader,
                basis         = basis,
                criterion     = criterion,
                is_algnet     = is_algnet,
                current_ep    = ep,
            )
            lambda_theta_history.append(probe["lambda_theta"])
            lambda_extended_history.append(probe["lambda_extended"])
            probe_epochs.append(ep)
            if ep in [100, 500, 1000, 2000, 3000]:
                delta_curves_sample[ep] = probe["delta_theta_curve"]

            if ep % REPORT_EVERY == 0:
                print(f"    ep {ep:4d}  val={val_acc:.4f}  best={best_val:.4f}"
                      f"  λ_θ={probe['lambda_theta']:+.4f}"
                      f"  λ_ext={probe['lambda_extended']:+.4f}")
        elif ep % REPORT_EVERY == 0:
            print(f"    ep {ep:4d}  val={val_acc:.4f}  best={best_val:.4f}")

    return {
        "name":                    name,
        "probe_epochs":            np.array(probe_epochs),
        "lambda_theta":            np.array(lambda_theta_history),
        "lambda_extended":         np.array(lambda_extended_history),
        "val_accs":                np.array(val_accs),
        "best_val":                best_val,
        "delta_curves_sample":     delta_curves_sample,
    }


# ---------------------------------------------------------------------------
#  Kaplan-Yorke dimension estimate
# ---------------------------------------------------------------------------

def kaplan_yorke_dim(lambda_series: np.ndarray) -> float:
    """
    Approximate d_KY from a sorted Lyapunov spectrum λ₁ ≥ λ₂ ≥ ...
    d_KY = j + (λ₁ + ... + λ_j) / |λ_{j+1}|
    where j is the largest index s.t. Σ_{k=1}^j λ_k > 0.

    Here we only have λ_max (scalar), so if λ_max > 0, d_KY > 1.
    For a single scalar estimate: d_KY_lower = 1 + λ_max / |λ_max|
    is trivially 2, so instead we report λ_max as the dynamical chaos indicator.
    A proper d_KY requires the full Lyapunov spectrum (QR method), which would
    require computing the full Jacobian — that is left as an open extension.
    """
    finite = lambda_series[np.isfinite(lambda_series)]
    if len(finite) == 0:
        return float("nan")
    return float(np.mean(finite))


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def main():
    t_start = time.time()

    X_tr_all = {}
    X_vl_all = {}
    for seed in BASE_SEEDS:
        X_tr_all[seed] = generate_su3_cubic_dataset(N_TRAIN, seed=seed)
        X_vl_all[seed] = generate_su3_cubic_dataset(N_VAL,   seed=seed + 100)

    basis = _su_basis(N_DIM, DEVICE)
    seed_results = []

    for seed in BASE_SEEDS:
        print(f"\n{'='*60}")
        print(f"Seed {seed}")
        print(f"{'='*60}")

        X_tr, y_tr = X_tr_all[seed]
        X_vl, y_vl = X_vl_all[seed]
        tr_loader = DataLoader(TensorDataset(X_tr, y_tr), batch_size=BATCH_SIZE, shuffle=True)
        vl_loader = DataLoader(TensorDataset(X_vl, y_vl), batch_size=512)

        torch.manual_seed(seed * 1000 + 1)
        algnet = AlgNet1L().to(DEVICE)
        torch.manual_seed(seed * 1000 + 2)
        mlp    = MLP1L().to(DEVICE)

        print("  --- AlgNet-1L ---")
        res_alg = run_with_lyapunov_probes(algnet, tr_loader, vl_loader, basis,
                                            "AlgNet-1L", is_algnet=True)
        print("  --- MLP-1L ---")
        res_mlp = run_with_lyapunov_probes(mlp, tr_loader, vl_loader, basis,
                                            "MLP-1L", is_algnet=False)
        seed_results.append({"seed": seed, "alg": res_alg, "mlp": res_mlp})

        # Per-seed summary
        for mkey, label in [("alg", "AlgNet"), ("mlp", "MLP")]:
            r = seed_results[-1][mkey]
            lam = r["lambda_theta"]
            print(f"  {label}: best_val={r['best_val']:.4f}  "
                  f"λ_θ early (med 10ep)={np.nanmedian(lam[:10]):+.4f}  "
                  f"λ_θ late (med 10ep)={np.nanmedian(lam[-10:]):+.4f}")

    # Aggregate summary
    print("\n" + "="*70)
    print("  Lyapunov summary  (median ± IQR/2  across seeds)")
    print("  " + "-"*60)
    for mkey, label in [("alg", "AlgNet-1L"), ("mlp", "MLP-1L")]:
        all_lam_early = np.array([np.nanmedian(r[mkey]["lambda_theta"][:10])
                                   for r in seed_results])
        all_lam_late  = np.array([np.nanmedian(r[mkey]["lambda_theta"][-10:])
                                   for r in seed_results])
        all_lam_ext_e = np.array([np.nanmedian(r[mkey]["lambda_extended"][:10])
                                   for r in seed_results])
        all_lam_ext_l = np.array([np.nanmedian(r[mkey]["lambda_extended"][-10:])
                                   for r in seed_results])
        print(f"  {label:<12}  λ_θ early={np.nanmedian(all_lam_early):+.4f}"
              f"  λ_θ late={np.nanmedian(all_lam_late):+.4f}"
              f"  λ_ext early={np.nanmedian(all_lam_ext_e):+.4f}"
              f"  λ_ext late={np.nanmedian(all_lam_ext_l):+.4f}")

    # Save NPZ
    npz_path = os.path.join(OUT_DIR, "exp13_results.npz")
    save_dict = {}
    for seed_r in seed_results:
        s = seed_r["seed"]
        save_dict[f"probe_epochs_s{s}"]   = seed_r["alg"]["probe_epochs"]
        save_dict[f"lam_theta_alg_s{s}"]  = seed_r["alg"]["lambda_theta"]
        save_dict[f"lam_ext_alg_s{s}"]    = seed_r["alg"]["lambda_extended"]
        save_dict[f"lam_theta_mlp_s{s}"]  = seed_r["mlp"]["lambda_theta"]
        save_dict[f"lam_ext_mlp_s{s}"]    = seed_r["mlp"]["lambda_extended"]
    np.savez(npz_path, **save_dict)
    print(f"\nSaved {npz_path}")

    _plot_results(seed_results)

    elapsed = time.time() - t_start
    print(f"\nTotal runtime: {elapsed/60:.1f} min")


# ---------------------------------------------------------------------------
#  Plotting
# ---------------------------------------------------------------------------

def _plot_results(seed_results: list):
    probes = seed_results[0]["alg"]["probe_epochs"]
    n_p    = len(probes)

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # Helper: aggregate λ across seeds
    def _agg(mkey: str, lam_key: str):
        arr = np.array([r[mkey][lam_key] for r in seed_results],
                       dtype=float)   # (N_REPS, n_probes)
        return np.nanmedian(arr, axis=0), np.nanpercentile(arr, 25, axis=0), np.nanpercentile(arr, 75, axis=0)

    # --- Row 0: λ_theta(t) ---
    ax = axes[0, 0]
    for mkey, color, label in [("alg", "steelblue", "AlgNet-1L"), ("mlp", "coral", "MLP-1L")]:
        med, lo, hi = _agg(mkey, "lambda_theta")
        ax.plot(probes, med, color=color, label=label, linewidth=1.5)
        ax.fill_between(probes, lo, hi, alpha=0.2, color=color)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel("Training epoch")
    ax.set_ylabel("λ_θ (per epoch)")
    ax.set_title("Max Lyapunov (θ-space)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # --- Row 0: λ_extended(t) ---
    ax = axes[0, 1]
    for mkey, color, label in [("alg", "steelblue", "AlgNet-1L"), ("mlp", "coral", "MLP-1L")]:
        med, lo, hi = _agg(mkey, "lambda_extended")
        ax.plot(probes, med, color=color, label=label, linewidth=1.5)
        ax.fill_between(probes, lo, hi, alpha=0.2, color=color)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    # Mark §11g reference λ ≈ 0.061/step (F₄ discovery) for context
    ax.axhline(0.061, color="purple", linewidth=0.8, linestyle=":",
               label="§11g F₄ λ≈0.061")
    ax.set_xlabel("Training epoch")
    ax.set_ylabel("λ_ext (per epoch)")
    ax.set_title("Max Lyapunov (extended: θ,m,v)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Row 0: λ_theta AlgNet vs MLP at early/late phases ---
    ax = axes[0, 2]
    n_early, n_late = 5, 5
    for xi, (mkey, color, label) in enumerate([("alg", "steelblue", "AlgNet"), ("mlp", "coral", "MLP")]):
        all_lam_e = [np.nanmedian(r[mkey]["lambda_theta"][:n_early]) for r in seed_results]
        all_lam_l = [np.nanmedian(r[mkey]["lambda_theta"][-n_late:])  for r in seed_results]
        ax.bar(xi*2,     np.nanmean(all_lam_e), color=color, alpha=0.8, width=0.7,
               label=f"{label} early",  edgecolor="black", linewidth=0.5)
        ax.bar(xi*2+0.8, np.nanmean(all_lam_l), color=color, alpha=0.4, width=0.7,
               label=f"{label} late",   edgecolor="black", linewidth=0.5)
        # Error bars
        ax.errorbar(xi*2,     np.nanmean(all_lam_e), yerr=np.nanstd(all_lam_e)/max(1,np.sqrt(len(all_lam_e))),
                    fmt="none", color="black", capsize=4)
        ax.errorbar(xi*2+0.8, np.nanmean(all_lam_l), yerr=np.nanstd(all_lam_l)/max(1,np.sqrt(len(all_lam_l))),
                    fmt="none", color="black", capsize=4)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xticks([0.35, 2.35])
    ax.set_xticklabels(["AlgNet", "MLP"])
    ax.set_ylabel("Median λ_θ")
    ax.set_title("Early vs Late λ_θ")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    # --- Row 1: divergence curve ||Δθ||_t at selected epochs (seed 0, AlgNet) ---
    ax = axes[1, 0]
    delta_alg = seed_results[0]["alg"]["delta_curves_sample"]
    colors_ep = plt.cm.viridis(np.linspace(0, 1, len(delta_alg)))
    for (ep, curve), c in zip(sorted(delta_alg.items()), colors_ep):
        steps = np.arange(len(curve))
        ax.semilogy(steps, curve, color=c, label=f"ep {ep}", linewidth=1.2)
    ax.axhline(EPSILON_0, color="black", linewidth=0.8, linestyle="--",
               label=f"ε₀={EPSILON_0:.0e}")
    ax.set_xlabel("Steps after injection")
    ax.set_ylabel("||Δθ||")
    ax.set_title("AlgNet: divergence curves (seed 0)")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # --- Row 1: same for MLP ---
    ax = axes[1, 1]
    delta_mlp = seed_results[0]["mlp"]["delta_curves_sample"]
    for (ep, curve), c in zip(sorted(delta_mlp.items()), colors_ep):
        steps = np.arange(len(curve))
        ax.semilogy(steps, curve, color=c, label=f"ep {ep}", linewidth=1.2)
    ax.axhline(EPSILON_0, color="black", linewidth=0.8, linestyle="--",
               label=f"ε₀={EPSILON_0:.0e}")
    ax.set_xlabel("Steps after injection")
    ax.set_ylabel("||Δθ||")
    ax.set_title("MLP: divergence curves (seed 0)")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # --- Row 1: validation accuracy comparison ---
    ax = axes[1, 2]
    ep_ax = np.arange(1, len(seed_results[0]["alg"]["val_accs"]) + 1)
    for mkey, color, label in [("alg", "steelblue", "AlgNet-1L"), ("mlp", "coral", "MLP-1L")]:
        acc_m = np.array([r[mkey]["val_accs"] for r in seed_results])
        ax.plot(ep_ax, acc_m.mean(axis=0), color=color, label=label, linewidth=1)
        ax.fill_between(ep_ax, acc_m.mean(0)-acc_m.std(0), acc_m.mean(0)+acc_m.std(0),
                        alpha=0.15, color=color)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Validation accuracy")
    ax.set_title("Learning curves (N=5)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.suptitle("Exp13: Lyapunov Exponents of Adam Dynamics  (AlgNet vs MLP)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "exp13_results.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved {path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
#  Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    main()
