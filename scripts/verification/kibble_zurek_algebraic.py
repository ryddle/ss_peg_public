"""
Task 3 — Kibble-Zurek Algebraic Phase Transition
═════════════════════════════════════════════════

Grid (N×N, periodic BC) of coupled algebraic optimizer instances.
Each cell: 3 generators in SU(3) hermitian-traceless parameterization.
→ cells crystallize as su(2) ⊂ su(3); the *orientation* of the
  su(2) embedding varies, creating domain structure.

Coupling
  Fingerprint  F_i = Σ_a  T_a^(i) ⊗ T_a^(i)   (frame operator).
  ||F_i − F_j||²  penalises orientation mismatch between neighbours.
  (Killing coupling cannot distinguish orientations of the same algebra.)

Quench protocol (crystallise → ramp coupling)
  Phase 1 — Independent crystallisation (PRE_STEPS):
    α = α_max, λ_c = 0.  Each cell crystallises to su(2) with
    a random orientation determined by its initial seed.
  Phase 2 — Linear λ-ramp (τ_Q steps):
    α = α_max, λ(t) = λ_max · (t − t₀ + 1) / τ_Q (linear ramp 0 → λ_max).
    Gradual onset lets orientations adjust before coupling saturates,
    reducing geometric frustration in CP² orientation space.
    Shorter τ_Q → steeper ramp → less coarsening → more frozen defects.
    Analogous to KZ freeze-out with adiabatic approach to transition.

Measurements
  · Mean domain size ⟨s⟩ vs τ_Q
  · Defect density ρ_d vs τ_Q   (fraction of misaligned neighbour pairs)
  · Mean fingerprint distance ⟨d⟩ vs τ_Q
  · KZ power-law fit: ⟨s⟩ ∝ τ_Q^β

Output
  output/kz/kz_results.md        — full report
  output/kz/kz_scaling.png       — power-law plot
  output/kz/kz_grids.png         — domain maps
  output/kz/kz_raw_data.json     — raw data

Usage
  python verification/kibble_zurek_algebraic.py              # full experiment
  python verification/kibble_zurek_algebraic.py --quick      # quick test (4×4)
  python verification/kibble_zurek_algebraic.py --grid 12    # larger grid
"""
import sys, os, time, json, argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import torch
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import Counter

from shared_utils_gpu import (
    params_to_generators,
    generators_to_params,
    random_params,
    batched_killing_metric,
    batched_non_commutativity_tension,
    batched_norm_tension,
    killing_metric_np,
    check_in_span_np,
    commutator_np,
    get_device_info,
    DEVICE, DTYPE, FDTYPE,
)

# ═══════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(WORKSPACE_ROOT, "output", "kz")

DEFAULT_CONFIG = dict(
    N_GRID=20,
    N_GEN=3,            # su(2) generators
    DIM=3,              # in SU(3) parameterization  (ρ = 3/8)
    ALPHA_MAX=5.0,      # fixed α, well into crystallisation regime
    BETA=1.0,
    LAMBDA_C=0.1,       # fingerprint coupling (domain alignment)
    KILLING_TARGET=-8.0,# su(2) Killing in dim=3 rep
    PRE_STEPS=15_000,   # Phase 1: independent crystallisation
    LR=0.001,
    TAU_Q_VALUES=[500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000],
    N_SEEDS=3,
    DOMAIN_THRESHOLD=0.10,
)

QUICK_CONFIG = dict(
    N_GRID=20,
    N_GEN=3,
    DIM=3,
    ALPHA_MAX=5.0,
    BETA=1.0,
    LAMBDA_C=0.1,
    KILLING_TARGET=-8.0,
    PRE_STEPS=15_000,
    LR=0.001,
    TAU_Q_VALUES=[500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000],
    N_SEEDS=1,
    DOMAIN_THRESHOLD=0.10,
)


# ═══════════════════════════════════════════════════════════════
# GRID UTILITIES
# ═══════════════════════════════════════════════════════════════
def build_neighbor_pairs(N):
    """Neighbor pairs for N×N grid with periodic BC (torus).

    Returns list of (cell_i, cell_j) for right and down neighbors.
    Total pairs = 2·N² (each cell contributes right + down, wrapping).
    """
    pairs = []
    for y in range(N):
        for x in range(N):
            cell = y * N + x
            right = y * N + ((x + 1) % N)
            down = ((y + 1) % N) * N + x
            pairs.append((cell, right))
            pairs.append((cell, down))
    return pairs


# ═══════════════════════════════════════════════════════════════
# FINGERPRINT — permutation-invariant, orientation-sensitive
# ═══════════════════════════════════════════════════════════════
def compute_fingerprint(generators):
    """Frame operator  F_i = Σ_a  T_a^(i) ⊗ T_a^(i).

    Parameters
    ----------
    generators : (B, n_gen, dim, dim) complex tensor

    Returns
    -------
    F : (B, dim², dim²) complex tensor
    """
    B, n_gen, dim, _ = generators.shape
    G_flat = generators.reshape(B, n_gen, dim * dim)
    return torch.einsum("bai,baj->bij", G_flat.conj(), G_flat)


# ═══════════════════════════════════════════════════════════════
# OPTIMIZATION LOOP — two-phase: crystallise then couple
# ═══════════════════════════════════════════════════════════════
def run_single_quench(tau_Q, seed, cfg):
    """Quench protocol: crystallise then ramp coupling.

    Phase 1 — Independent crystallisation (PRE_STEPS):
      α = α_max, λ_c = 0.  Each cell crystallises to su(2)
      with a random orientation set by its initial seed.

    Phase 2 — Linear λ-ramp (tau_Q steps):
      α = α_max, λ(t) = λ_max · (t − pre_steps + 1) / τ_Q (0 → λ_max).
      Gradual onset reduces geometric frustration in CP² orientation space.
      Shorter τ_Q → steeper ramp → less coarsening → more frozen defects.

    Returns
    -------
    final_np : ndarray (n_cells, n_gen, dim, dim) complex
    fingerprints_np : ndarray (n_cells, dim², dim²) complex
    history : dict
    """
    N = cfg["N_GRID"]
    n_cells = N * N
    n_gen = cfg["N_GEN"]
    dim = cfg["DIM"]
    pre_steps = cfg["PRE_STEPS"]
    lr = cfg["LR"]

    # Neighbor index tensors
    pairs = build_neighbor_pairs(N)
    i_idx = torch.tensor([p[0] for p in pairs], device=DEVICE)
    j_idx = torch.tensor([p[1] for p in pairs], device=DEVICE)

    # Random initial parameters — distinct seed per cell
    cell_seeds = [seed * 10_000 + i for i in range(n_cells)]
    params = random_params(n_cells, n_gen, dim, device=DEVICE, seeds=cell_seeds)

    # Normalize initial generators
    with torch.no_grad():
        gens = params_to_generators(params, n_gen, dim=dim)
        norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
        gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
        params = generators_to_params(gens)
    params = params.clone().detach().requires_grad_(True)

    optimizer = torch.optim.Adam([params], lr=lr)
    alpha = cfg["ALPHA_MAX"]
    K_target_mat = cfg["KILLING_TARGET"] * torch.eye(
        n_gen, dtype=FDTYPE, device=DEVICE
    )

    total_steps = pre_steps + tau_Q
    history = {
        "step": [], "total_loss": [], "cell_loss": [],
        "coupling_loss": [], "phase": [],
    }

    for step in range(total_steps):
        optimizer.zero_grad()
        generators = params_to_generators(params, n_gen, dim=dim)

        # Per-cell tensions
        T_nc = batched_non_commutativity_tension(generators)
        K = batched_killing_metric(generators)
        T_k = ((K - K_target_mat.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
        T_norm = batched_norm_tension(generators)
        cell_loss = (T_nc + alpha * T_k + cfg["BETA"] * T_norm).sum()

        # Phase 2: fingerprint coupling (linear λ-ramp: 0 → λ_max over τ_Q steps)
        if step >= pre_steps:
            t_phase2 = step - pre_steps
            lambda_t = cfg["LAMBDA_C"] * (t_phase2 + 1) / tau_Q
            F = compute_fingerprint(generators)
            F_diff = F[i_idx] - F[j_idx]
            coupling_loss = lambda_t * (
                (F_diff.real ** 2 + F_diff.imag ** 2).sum(dim=(-2, -1)).sum()
            )
        else:
            coupling_loss = torch.tensor(0.0, device=DEVICE)

        total_loss = cell_loss + coupling_loss
        total_loss.backward()

        # Clip gradients to prevent NaN (safety)
        torch.nn.utils.clip_grad_norm_([params], max_norm=100.0)

        optimizer.step()

        # Reset optimizer at phase boundary → fresh Adam for Phase 2
        if step == pre_steps - 1:
            optimizer = torch.optim.Adam([params], lr=lr)

        if step % 2000 == 0 or step == total_steps - 1 or step == pre_steps:
            phase = "crystal" if step < pre_steps else "ramp"
            with torch.no_grad():
                cl = coupling_loss.item()
                history["step"].append(step)
                history["total_loss"].append(total_loss.item())
                history["cell_loss"].append(cell_loss.item())
                history["coupling_loss"].append(cl)
                history["phase"].append(phase)
                print(
                    f"    step {step:>5d}/{total_steps}  [{phase:>7s}]  "
                    f"cell={cell_loss.item():>10.2f}  "
                    f"coup={cl:>8.2f}  "
                    f"total={total_loss.item():>10.2f}"
                )

    with torch.no_grad():
        final_gens = params_to_generators(params, n_gen, dim=dim)
        final_F = compute_fingerprint(final_gens)
        final_np = final_gens.cpu().numpy()
        fp_np = final_F.cpu().numpy()

    return final_np, fp_np, history


# ═══════════════════════════════════════════════════════════════
# CELL CLASSIFICATION (diagnostic — confirms su(2) crystallisation)
# ═══════════════════════════════════════════════════════════════
def classify_cell(gens_single):
    """Classify emerged algebra from one cell's generators (numpy)."""
    generators = list(gens_single)
    n = len(generators)

    K = killing_metric_np(generators)
    K_eigs = np.linalg.eigvalsh(K)
    K_diag_mean = float(np.mean(np.diag(K)))

    total_pairs = n * (n - 1) // 2
    closes = 0
    for i in range(n):
        for j in range(i + 1, n):
            C = commutator_np(generators[i], generators[j])
            if np.linalg.norm(C) < 1e-6:
                closes += 1
            elif check_in_span_np(C, generators, tol=0.3):
                closes += 1
    closure = closes / total_pairs if total_pairs > 0 else 0.0

    eig_scale = max(float(np.abs(K_eigs).max()), 1e-8)
    tol = 0.01 * eig_scale
    n_neg = int(np.sum(K_eigs < -tol))
    n_pos = int(np.sum(K_eigs > tol))
    n_zero = n - n_neg - n_pos

    kappa = float(np.mean([np.trace(g @ g).real for g in generators]))
    k_ratio = abs(K_diag_mean) / abs(kappa) ** 2 if abs(kappa) > 1e-10 else 0.0

    return {
        "closure": float(closure),
        "K_eigs": K_eigs.tolist(),
        "K_diag_mean": K_diag_mean,
        "n_neg": n_neg,
        "n_pos": n_pos,
        "n_zero": n_zero,
        "k_ratio": float(k_ratio),
        "h_dual": float(k_ratio),
    }


def algebra_label(cls):
    """Human-readable algebra label for reporting."""
    if cls["closure"] < 0.5:
        return "disordered"
    # su(2) in dim=3: n_neg=3, h_dual≈4 (Killing=-8 with unit-norm gens)
    if cls["n_neg"] == 3 and 3.0 < cls["h_dual"] < 5.0:
        return "su2"
    return f"n{cls['n_neg']}_h{round(cls['h_dual'])}"


# ═══════════════════════════════════════════════════════════════
# DOMAIN DETECTION — fingerprint-based (orientation-sensitive)
# ═══════════════════════════════════════════════════════════════
def fingerprint_distance(fp_i, fp_j):
    """Relative Frobenius distance between two fingerprints."""
    diff_norm = np.linalg.norm(fp_i - fp_j)
    ref_norm = max(np.linalg.norm(fp_i), np.linalg.norm(fp_j), 1e-12)
    return diff_norm / ref_norm


def detect_domains_fp(fingerprints, N, threshold):
    """Connected components on N×N periodic grid based on fingerprint similarity.

    Two neighbouring cells are in the same domain iff
    their relative fingerprint distance < threshold.
    """
    n_cells = N * N
    visited = [False] * n_cells
    domains = []

    for start in range(n_cells):
        if visited[start]:
            continue
        cells = []
        stack = [start]
        while stack:
            cell = stack.pop()
            if visited[cell]:
                continue
            visited[cell] = True
            cells.append(cell)
            x, y = cell % N, cell // N
            for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                nx, ny = (x + dx) % N, (y + dy) % N
                nc = ny * N + nx
                if not visited[nc]:
                    d = fingerprint_distance(fingerprints[cell], fingerprints[nc])
                    if d < threshold:
                        stack.append(nc)
        domains.append({"cells": cells, "size": len(cells)})

    return domains


def compute_defect_density_fp(fingerprints, N, threshold):
    """Fraction of neighbour pairs whose fingerprints differ (ρ_d)."""
    n_mismatch = 0
    n_total = 0
    for y in range(N):
        for x in range(N):
            cell = y * N + x
            right = y * N + ((x + 1) % N)
            down = ((y + 1) % N) * N + x
            if fingerprint_distance(fingerprints[cell], fingerprints[right]) >= threshold:
                n_mismatch += 1
            if fingerprint_distance(fingerprints[cell], fingerprints[down]) >= threshold:
                n_mismatch += 1
            n_total += 2
    return n_mismatch / n_total if n_total > 0 else 0.0


def mean_fp_distance(fingerprints, N):
    """Mean neighbour fingerprint distance across the whole grid."""
    dists = []
    for y in range(N):
        for x in range(N):
            cell = y * N + x
            right = y * N + ((x + 1) % N)
            down = ((y + 1) % N) * N + x
            dists.append(fingerprint_distance(fingerprints[cell], fingerprints[right]))
            dists.append(fingerprint_distance(fingerprints[cell], fingerprints[down]))
    return float(np.mean(dists))


# ═══════════════════════════════════════════════════════════════
# MAIN EXPERIMENT LOOP
# ═══════════════════════════════════════════════════════════════
def run_experiment(cfg):
    N = cfg["N_GRID"]
    n_cells = N * N
    tau_values = cfg["TAU_Q_VALUES"]
    n_seeds = cfg["N_SEEDS"]
    total_runs = len(tau_values) * n_seeds
    threshold = cfg["DOMAIN_THRESHOLD"]

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 72)
    print("  TASK 3 \u2014 KIBBLE-ZUREK ALGEBRAIC PHASE TRANSITION")
    print("=" * 72)
    print(f"  Device:       {get_device_info()}")
    print(f"  Grid:         {N}\u00d7{N} = {n_cells} cells (periodic BC)")
    print(f"  Per cell:     {cfg['N_GEN']} generators in SU({cfg['DIM']})")
    print(f"  \u03b1:            {cfg['ALPHA_MAX']}")
    print(f"  λ_coupling:   {cfg['LAMBDA_C']}  (fingerprint, linear ramp 0→λ_max)")
    print(f"  K_target:     {cfg['KILLING_TARGET']}")
    print(f"  Pre-steps:    {cfg['PRE_STEPS']}  (independent crystallisation)")
    print(f"  LR:           {cfg['LR']}")
    print(f"  τ_Q values:   {tau_values}  (coupling duration)")
    print(f"  Seeds:        {n_seeds}")
    print(f"  Total runs:   {total_runs}")
    print(f"  Domain thr:   {threshold}")
    print()

    all_results = []
    t_start = time.time()

    for ti, tau_Q in enumerate(tau_values):
        for si in range(n_seeds):
            seed = si + 1
            run_num = ti * n_seeds + si + 1
            print(
                f"  \u2500\u2500 Run {run_num}/{total_runs}: "
                f"\u03c4_Q={tau_Q}, seed={seed} \u2500\u2500"
            )

            t0 = time.time()
            final_gens, fp_np, history = run_single_quench(tau_Q, seed, cfg)
            gpu_time = time.time() - t0

            # classify (diagnostic)
            print(f"    Classifying {n_cells} cells ...")
            t1 = time.time()
            classifications = []
            for c in range(n_cells):
                cls = classify_cell(final_gens[c])
                classifications.append(cls)
            classify_time = time.time() - t1

            # domain detection (fingerprint-based)
            domains = detect_domains_fp(fp_np, N, threshold)
            defect_density = compute_defect_density_fp(fp_np, N, threshold)
            mean_d = mean_fp_distance(fp_np, N)

            # build per-cell domain ID array (for grid plot)
            cell_domain_id = np.full(n_cells, -1, dtype=int)
            for dom_idx, dom in enumerate(domains):
                for c in dom["cells"]:
                    cell_domain_id[c] = dom_idx

            mean_closure = float(np.mean([c["closure"] for c in classifications]))
            domain_sizes = [d["size"] for d in domains]
            mean_domain_size = float(np.mean(domain_sizes))
            label_counts = Counter(algebra_label(c) for c in classifications)

            result = {
                "tau_Q": tau_Q,
                "seed": seed,
                "n_domains": len(domains),
                "mean_domain_size": mean_domain_size,
                "max_domain_size": int(max(domain_sizes)),
                "defect_density": defect_density,
                "mean_fp_distance": mean_d,
                "mean_closure": mean_closure,
                "label_distribution": dict(label_counts),
                "domain_sizes": domain_sizes,
                "cell_domain_ids": cell_domain_id.tolist(),
                "gpu_time": gpu_time,
                "classify_time": classify_time,
                "history": history,
                "classifications": [
                    {k: v for k, v in c.items() if k != "K_eigs"}
                    for c in classifications
                ],
                "fp_norms": [float(np.linalg.norm(fp_np[c])) for c in range(n_cells)],
            }
            all_results.append(result)

            print(
                f"    \u2713 {len(domains)} domains, \u27e8s\u27e9={mean_domain_size:.1f}, "
                f"\u03c1_d={defect_density:.3f}, \u27e8d\u27e9={mean_d:.4f}, "
                f"\u27e8closure\u27e9={mean_closure:.2f}"
            )
            print(f"      Labels: {dict(label_counts)}")
            print(
                f"      Time: GPU={gpu_time:.1f}s, "
                f"classify={classify_time:.1f}s"
            )
            print()

            torch.cuda.empty_cache()

    total_time = time.time() - t_start
    print(f"  Total experiment time: {total_time:.0f}s ({total_time / 60:.1f} min)")

    return all_results


# ═══════════════════════════════════════════════════════════════
# ANALYSIS AND OUTPUT
# ═══════════════════════════════════════════════════════════════
def analyze_and_output(all_results, cfg):
    N = cfg["N_GRID"]
    tau_values = cfg["TAU_Q_VALUES"]

    # ── aggregate by τ_Q ──
    agg = {}
    for r in all_results:
        tq = r["tau_Q"]
        if tq not in agg:
            agg[tq] = {
                "domain_sizes": [], "defect_densities": [],
                "n_domains": [], "closures": [], "fp_dists": [],
            }
        agg[tq]["domain_sizes"].append(r["mean_domain_size"])
        agg[tq]["defect_densities"].append(r["defect_density"])
        agg[tq]["n_domains"].append(r["n_domains"])
        agg[tq]["closures"].append(r["mean_closure"])
        agg[tq]["fp_dists"].append(r["mean_fp_distance"])

    tau_arr = np.array(sorted(agg.keys()), dtype=float)
    mean_s = np.array([np.mean(agg[t]["domain_sizes"]) for t in tau_arr])
    std_s = np.array([np.std(agg[t]["domain_sizes"]) for t in tau_arr])
    mean_rho = np.array([np.mean(agg[t]["defect_densities"]) for t in tau_arr])
    std_rho = np.array([np.std(agg[t]["defect_densities"]) for t in tau_arr])
    mean_nd = np.array([np.mean(agg[t]["n_domains"]) for t in tau_arr])
    mean_cl = np.array([np.mean(agg[t]["closures"]) for t in tau_arr])
    mean_fp = np.array([np.mean(agg[t]["fp_dists"]) for t in tau_arr])

    # ── power-law fits (exclude saturated points) ──
    n_cells = N * N
    # Fit regime: 1.05 < ⟨s⟩ < n_cells/2  (pre-saturation, some merging)
    fit_mask = (mean_s > 1.05) & (mean_s < n_cells / 2)
    log_tau = np.log(tau_arr)
    log_s = np.log(np.clip(mean_s, 1e-6, None))

    if fit_mask.sum() >= 2:
        fit_s = np.polyfit(log_tau[fit_mask], log_s[fit_mask], 1)
        beta_kz = float(fit_s[0])
    else:
        fit_s = np.polyfit(log_tau, log_s, 1) if len(log_tau) >= 2 else np.array([0.0, 0.0])
        beta_kz = float(fit_s[0])
        fit_mask = np.ones(len(tau_arr), dtype=bool)

    rho_fit_mask = fit_mask & (mean_rho > 1e-6)
    log_rho = np.log(np.clip(mean_rho, 1e-6, None))
    if rho_fit_mask.sum() >= 2:
        fit_rho = np.polyfit(log_tau[rho_fit_mask], log_rho[rho_fit_mask], 1)
        delta_kz = float(-fit_rho[0])
    else:
        fit_rho = np.array([0.0, 0.0])
        delta_kz = 0.0

    n_fit_pts = int(fit_mask.sum())
    fit_range = f"τ_Q ∈ [{tau_arr[fit_mask][0]:.0f}, {tau_arr[fit_mask][-1]:.0f}]" if n_fit_pts > 0 else "N/A"

    # ── Plot 1: scaling laws ──
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    ax = axes[0]
    ax.errorbar(
        tau_arr, mean_s, yerr=std_s,
        fmt="o-", color="#2196F3", capsize=4, linewidth=2, markersize=8,
        label="Data",
    )
    # Show fit only in the fit range
    fit_x = np.linspace(log_tau[fit_mask].min(), log_tau[fit_mask].max(), 100) if n_fit_pts >= 2 else log_tau
    ax.plot(
        np.exp(fit_x), np.exp(fit_s[0] * fit_x + fit_s[1]),
        "--", color="#F44336", linewidth=1.5,
        label=f"Fit: ⟨s⟩ ∝ τ_Q^{{{beta_kz:.2f}}} ({n_fit_pts} pts)",
    )
    # Mark saturated / non-scaling points
    if fit_mask.sum() < len(tau_arr):
        ax.scatter(tau_arr[~fit_mask], mean_s[~fit_mask], s=120,
                   facecolors="none", edgecolors="#999", linewidths=2, zorder=5,
                   label="Excluded (saturation)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("τ_Q (quench time)", fontsize=12)
    ax.set_ylabel("⟨s⟩ (mean domain size)", fontsize=12)
    ax.set_title("Domain Size Scaling", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10); ax.grid(True, alpha=0.3)

    ax = axes[1]
    rho_pos = mean_rho > 0
    if rho_pos.sum() > 0:
        ax.errorbar(
            tau_arr[rho_pos], mean_rho[rho_pos], yerr=std_rho[rho_pos],
            fmt="s-", color="#FF9800", capsize=4, linewidth=2, markersize=8,
            label="Data",
        )
        if rho_fit_mask.sum() >= 2:
            rho_fit_x = np.linspace(log_tau[rho_fit_mask].min(), log_tau[rho_fit_mask].max(), 100)
            ax.plot(
                np.exp(rho_fit_x), np.exp(fit_rho[0] * rho_fit_x + fit_rho[1]),
                "--", color="#9C27B0", linewidth=1.5,
                label=f"Fit: ρ_d ∝ τ_Q^{{{fit_rho[0]:.2f}}}",
            )
        ax.set_xscale("log"); ax.set_yscale("log")
    else:
        ax.plot(tau_arr, mean_rho, "s-", color="#FF9800", linewidth=2, markersize=8)
        ax.set_xscale("log")
    ax.set_xlabel("τ_Q", fontsize=12)
    ax.set_ylabel("ρ_d (defect density)", fontsize=12)
    ax.set_title("Defect Density Scaling", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10); ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.plot(tau_arr, mean_fp, "D-", color="#4CAF50", linewidth=2, markersize=8)
    ax.set_xscale("log")
    ax.set_xlabel("\u03c4_Q", fontsize=12)
    ax.set_ylabel("\u27e8d\u27e9 (mean FP distance)", fontsize=12)
    ax.set_title("Mean Fingerprint Distance", fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        f"Kibble-Zurek Algebraic Phase Transition \u2014 "
        f"{N}\u00d7{N} grid, SU({cfg['DIM']}) / {cfg['N_GEN']} gens",
        fontsize=14, fontweight="bold",
    )
    plt.tight_layout()
    plt.savefig(
        os.path.join(OUTPUT_DIR, "kz_scaling.png"), dpi=150, bbox_inches="tight"
    )
    plt.close()

    # ── Plot 2: grid domain maps ──
    n_tau = len(tau_values)
    n_seeds_plot = cfg["N_SEEDS"]

    fig, axes = plt.subplots(
        n_seeds_plot, n_tau, figsize=(3 * n_tau, 3 * n_seeds_plot), squeeze=False,
    )
    for r in all_results:
        tau_idx = tau_values.index(r["tau_Q"])
        seed_idx = r["seed"] - 1
        ax = axes[seed_idx][tau_idx]

        grid = np.array(r["cell_domain_ids"]).reshape(N, N)
        n_dom = r["n_domains"]
        ax.imshow(
            grid, cmap=plt.cm.tab20, vmin=0, vmax=max(n_dom - 1, 1),
            interpolation="nearest",
        )
        ax.set_title(
            f"\u03c4_Q={r['tau_Q']}, s={r['seed']}\n"
            f"{n_dom} dom, \u03c1={r['defect_density']:.2f}",
            fontsize=8,
        )
        ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle(
        "Domain Maps (colour \u2248 fingerprint norm bin)", fontsize=13, fontweight="bold",
    )
    plt.tight_layout()
    plt.savefig(
        os.path.join(OUTPUT_DIR, "kz_grids.png"), dpi=150, bbox_inches="tight"
    )
    plt.close()

    # ── Markdown report ──
    lines = []
    lines.append("# Task 3 \u2014 Kibble-Zurek Algebraic Phase Transition: Results\n")
    lines.append(f"**Date**: {time.strftime('%Y-%m-%d %H:%M')}\n")
    lines.append(f"**Grid**: {N}\u00d7{N} = {N*N} cells (periodic BC)  ")
    lines.append(
        f"**Parameterization**: {cfg['N_GEN']} generators in SU({cfg['DIM']})  "
    )
    lines.append(f"**\u03b1**: {cfg['ALPHA_MAX']}  ")
    lines.append(f"**K_target**: {cfg['KILLING_TARGET']}  ")
    lines.append(f"**Coupling**: λ_c = {cfg['LAMBDA_C']}  (fingerprint, linear ramp 0→λ_max)  ")
    lines.append(f"**Pre-steps**: {cfg['PRE_STEPS']}  (independent crystallisation)  ")
    lines.append(f"**LR**: {cfg['LR']}  ")
    lines.append(f"**τ_Q values**: {tau_values}  (coupling duration)  ")
    lines.append(f"**Seeds per \u03c4_Q**: {cfg['N_SEEDS']}  ")
    lines.append(f"**Domain threshold**: {cfg['DOMAIN_THRESHOLD']}\n")

    lines.append("## 1. Per-Run Summary\n")
    lines.append(
        "| \u03c4_Q | Seed | Domains | \u27e8s\u27e9 | Max s | \u03c1_d | \u27e8d\u27e9 | \u27e8closure\u27e9 | Labels |"
    )
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for r in all_results:
        lab = ", ".join(
            f"{k}:{v}" for k, v in sorted(r["label_distribution"].items())
        )
        lines.append(
            f"| {r['tau_Q']} | {r['seed']} | {r['n_domains']} | "
            f"{r['mean_domain_size']:.1f} | {r['max_domain_size']} | "
            f"{r['defect_density']:.3f} | {r['mean_fp_distance']:.4f} | "
            f"{r['mean_closure']:.3f} | {lab} |"
        )
    lines.append("")

    lines.append("## 2. Aggregated by \u03c4_Q\n")
    lines.append(
        "| \u03c4_Q | \u27e8s\u27e9 (mean\u00b1std) | \u27e8\u03c1_d\u27e9 (mean\u00b1std) | \u27e8d\u27e9 | \u27e8N_dom\u27e9 | \u27e8closure\u27e9 |"
    )
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for i, t in enumerate(tau_arr):
        lines.append(
            f"| {int(t)} | {mean_s[i]:.2f}\u00b1{std_s[i]:.2f} | "
            f"{mean_rho[i]:.4f}\u00b1{std_rho[i]:.4f} | "
            f"{mean_fp[i]:.4f} | "
            f"{mean_nd[i]:.1f} | {mean_cl[i]:.3f} |"
        )
    lines.append("")

    lines.append("## 3. Power-Law Fits\n")
    lines.append(
        f"**Domain size**: ⟨s⟩ ∝ τ_Q^β with **β = {beta_kz:.4f}**  "
    )
    lines.append(
        f"**Defect density**: ρ_d ∝ τ_Q^(−δ) with **δ = {delta_kz:.4f}**  "
    )
    lines.append(f"**Fit range**: {fit_range} ({n_fit_pts} data points)  ")
    lines.append(f"**Excluded**: {len(tau_arr) - n_fit_pts} points (⟨s⟩≤1.05 or ⟨s⟩≥N²/2 → saturation)")
    lines.append("")
    if abs(beta_kz) > 0.01:
        lines.append(
            f"KZ interpretation: \u03b2 = 2\u03bd/(1+z\u03bd) \u2192 measured = {beta_kz:.4f}"
        )
    else:
        lines.append(
            "**Note**: Exponent near zero \u2014 no clear power law detected.  "
        )
        lines.append(
            "Consider adjusting \u03bb_c, \u03b1 range, or grid size."
        )
    lines.append("")

    lines.append("## 4. Plots\n")
    lines.append("![Scaling](kz_scaling.png)\n")
    lines.append("![Grids](kz_grids.png)\n")

    report_path = os.path.join(OUTPUT_DIR, "kz_results.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    # ── JSON ──
    json_data = {
        "config": cfg,
        "aggregated": {
            "tau_Q": tau_arr.tolist(),
            "mean_domain_size": mean_s.tolist(),
            "std_domain_size": std_s.tolist(),
            "mean_defect_density": mean_rho.tolist(),
            "std_defect_density": std_rho.tolist(),
            "mean_fp_distance": mean_fp.tolist(),
            "beta_kz": beta_kz,
            "delta_kz": delta_kz,
        },
        "runs": [
            {k: v for k, v in r.items() if k != "classifications"}
            for r in all_results
        ],
    }
    json_path = os.path.join(OUTPUT_DIR, "kz_raw_data.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(json_data, f, indent=2, default=str)

    print(f"\n  Output:")
    print(f"    {report_path}")
    print(f"    {os.path.join(OUTPUT_DIR, 'kz_scaling.png')}")
    print(f"    {os.path.join(OUTPUT_DIR, 'kz_grids.png')}")
    print(f"    {json_path}")

    return beta_kz, delta_kz


# ═══════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(
        description="Task 3: Kibble-Zurek Algebraic Phase Transition"
    )
    parser.add_argument(
        "--quick", action="store_true",
        help="Quick test: 4\u00d74 grid, 2 \u03c4_Q values, 1 seed",
    )
    parser.add_argument("--grid", type=int, help="Override grid size N")
    parser.add_argument("--pre-steps", type=int, help="Override pre-crystallisation steps")
    parser.add_argument("--seeds", type=int, help="Override number of seeds")
    parser.add_argument(
        "--lambda-c", type=float, help="Override coupling strength"
    )
    parser.add_argument(
        "--threshold", type=float, help="Override domain detection threshold"
    )
    args = parser.parse_args()

    cfg = QUICK_CONFIG.copy() if args.quick else DEFAULT_CONFIG.copy()
    if args.grid:
        cfg["N_GRID"] = args.grid
    if args.pre_steps:
        cfg["PRE_STEPS"] = args.pre_steps
    if args.seeds:
        cfg["N_SEEDS"] = args.seeds
    if args.lambda_c is not None:
        cfg["LAMBDA_C"] = args.lambda_c
    if args.threshold is not None:
        cfg["DOMAIN_THRESHOLD"] = args.threshold

    all_results = run_experiment(cfg)
    beta_kz, delta_kz = analyze_and_output(all_results, cfg)

    print(f"\n  {'=' * 72}")
    print(f"  RESULT: \u03b2_KZ = {beta_kz:.4f}, \u03b4_KZ = {delta_kz:.4f}")
    print(f"  {'=' * 72}")


if __name__ == "__main__":
    main()
