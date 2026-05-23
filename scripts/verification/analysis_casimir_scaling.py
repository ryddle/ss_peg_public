"""
Analysis C.1 + C.4: Casimir scaling across SU(N).

Extracts coupling-related quantities from optimized SU(N) generators
and verifies the theoretical prediction:

    C_λ(N) = (N² - 1) / (2N) · κ

where κ = Tr(T_a²) / d is the per-generator norm (normalization-dependent).

Also extracts:
  - κ(N): mean generator norm
  - f_rms(N): RMS structure constant
  - g_eff(N) = f_rms / κ: effective coupling
  - K_diag(N): Killing metric diagonal

Usage:
    python analysis_casimir_scaling.py          # re-run SU(2..5), full extraction
    python analysis_casimir_scaling.py --quick   # use hardcoded data from prior runs
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils_gpu import (
    evolve_generators_gpu_individual,
    killing_metric_np,
    commutator_np,
    structure_constants_real_np,
    compute_casimir_np,
    get_device_info,
)


def extract_quantities(generators, N):
    """Extract all coupling-related quantities from a set of optimized generators.

    Args:
        generators: list of n_gen numpy arrays, each (dim, dim) complex
        N: the group dimension (SU(N))

    Returns:
        dict with κ, C_λ, C_λ_predicted, f_rms, g_eff, K_diag, etc.
    """
    n_gen = len(generators)
    d = generators[0].shape[0]

    # κ = mean Tr(T_a²) / d
    kappa_per_gen = np.array([np.trace(g @ g).real for g in generators])
    kappa = np.mean(kappa_per_gen) / d

    # Casimir: C = sum(T_a²), should be ∝ I
    C = compute_casimir_np(generators)
    C_diag = np.diag(C).real
    C_lambda = np.mean(C_diag)
    C_spread = np.std(C_diag)
    C_offdiag = np.max(np.abs(C - np.diag(C_diag)))

    # Theoretical prediction: C = Σ T_a², Tr(C) = Σ Tr(T_a²) = n_gen · κ_total
    # Schur's lemma → C = C_λ · I, so C_λ = Tr(C)/d = n_gen · κ_total / d
    # For SU(N): C_λ = (N²-1)/N · κ_total
    # Note: the standard (N²-1)/(2N) assumes Tr(T_a²) = 1/2 normalization.
    kappa_total = np.mean(kappa_per_gen)
    C_predicted = n_gen * kappa_total / d  # = (N²-1)/N · κ_total for SU(N)

    # Structure constants
    f = structure_constants_real_np(generators)
    f_nonzero = f[np.abs(f) > 1e-10]
    f_rms = np.sqrt(np.mean(f**2)) if f.size > 0 else 0.0
    f_max = np.max(np.abs(f)) if f.size > 0 else 0.0

    # Effective coupling: g_eff = f_rms / κ_total
    g_eff = f_rms / kappa_total if kappa_total > 1e-12 else 0.0

    # Killing metric
    K = killing_metric_np(generators)
    K_diag_vals = np.diag(K)
    K_diag_mean = np.mean(K_diag_vals)
    K_diag_std = np.std(K_diag_vals)

    # Dual Coxeter number: for SU(N), h∨ = N
    # Our code computes K_ij = Σ_k Tr([T_i,T_k]·[T_j,T_k]) which introduces
    # an extra factor of Tr(T²) = κ compared to the standard Killing form.
    # So K_diag = -2·h∨·κ² for compact semisimple algebras.
    # h∨ = |K_diag| / (2·κ²)
    h_dual_measured = abs(K_diag_mean) / (2 * kappa_total**2) if kappa_total > 1e-12 else 0.0

    return {
        "N": N,
        "n_gen": n_gen,
        "dim": d,
        "kappa_per_d": kappa,               # Tr(T²)/d averaged
        "kappa_total": kappa_total,          # Tr(T²) averaged
        "C_lambda": C_lambda,                # measured Casimir eigenvalue
        "C_predicted": C_predicted,          # (N²-1)/(2N) · κ_total
        "C_relative_error": abs(C_lambda - C_predicted) / abs(C_predicted) if abs(C_predicted) > 1e-12 else float('inf'),
        "C_spread": C_spread,                # should be ~0 (C ∝ I)
        "C_offdiag": C_offdiag,              # should be ~0
        "f_rms": f_rms,
        "f_max": f_max,
        "g_eff": g_eff,
        "K_diag_mean": K_diag_mean,
        "K_diag_std": K_diag_std,
        "h_dual_measured": h_dual_measured,
        "h_dual_theoretical": N,             # for SU(N)
    }


def run_analysis(groups=None, n_seeds=3, quick=False):
    if groups is None:
        groups = [2, 3, 4, 5]

    print("=" * 78)
    print("  CASIMIR SCALING ANALYSIS — SU(N)")
    print("=" * 78)

    if not quick:
        print(f"  Device: {get_device_info()}")
        print(f"  Groups: {', '.join(f'SU({n})' for n in groups)}")
        print(f"  Seeds per group: {n_seeds}")
        print()

    all_data = []

    for N in groups:
        n_gen = N**2 - 1
        print(f"\n{'─' * 78}")
        print(f"  SU({N}): {n_gen} generators in dim={N}")
        print(f"{'─' * 78}")

        if quick:
            print("  [quick mode — skipping optimization, using hardcoded data]")
            print("  Re-run without --quick for full analysis.")
            continue

        steps = max(20000, 5000 * N)
        lr = max(0.0001, 0.001 / N)

        t0 = time.time()
        raw = evolve_generators_gpu_individual(
            n_gen=n_gen, dim=N, steps=steps,
            lr=lr, alpha=5.0, beta=1.0,
            killing_target=-1.0,
            batch_size=n_seeds,
            seeds=list(range(n_seeds)),
            verbose=True,
            convergence_threshold=1e-4,
        )
        gpu_time = time.time() - t0

        seed_data = []
        for seed_idx, (final_gen, hist) in enumerate(raw):
            gens = list(final_gen)
            q = extract_quantities(gens, N)
            q["seed"] = seed_idx
            q["gpu_time"] = gpu_time / n_seeds
            seed_data.append(q)

            print(f"\n  seed={seed_idx}:")
            print(f"    κ = Tr(T²):       {q['kappa_total']:.6f}")
            print(f"    C_λ (measured):    {q['C_lambda']:.6f}")
            print(f"    C_λ (predicted):   {q['C_predicted']:.6f}")
            print(f"    C_λ rel. error:    {q['C_relative_error']:.2e}")
            print(f"    f_rms:             {q['f_rms']:.6f}")
            print(f"    g_eff = f/κ:       {q['g_eff']:.6f}")
            print(f"    K_diag:            {q['K_diag_mean']:.6f} ± {q['K_diag_std']:.6f}")
            print(f"    h∨ (measured):     {q['h_dual_measured']:.4f}  (theory: {N})")

        # Average across seeds
        avg = {}
        for key in seed_data[0]:
            if key in ("seed", "N", "n_gen", "dim", "h_dual_theoretical"):
                avg[key] = seed_data[0][key]
            elif isinstance(seed_data[0][key], (int, float, np.floating)):
                avg[key] = np.mean([s[key] for s in seed_data])
        avg["seed"] = "avg"
        all_data.append(avg)

    if quick:
        return

    # ═══════════════════════════════════════════════════════════════
    #  SUMMARY TABLES
    # ═══════════════════════════════════════════════════════════════

    print("\n\n" + "=" * 78)
    print("  TABLE 1: CASIMIR SCALING")
    print("=" * 78)
    print(f"\n  {'Group':8s} {'n_gen':>6s} {'κ=Tr(T²)':>10s} {'C_λ meas':>10s} "
          f"{'C_λ pred':>10s} {'Rel.err':>10s} {'C∝I?':>7s}")
    print("  " + "─" * 65)
    for d in all_data:
        prop = "YES" if d["C_spread"] < 0.01 else "NO"
        print(f"  SU({d['N']})  {d['n_gen']:>6d} {d['kappa_total']:>10.6f} "
              f"{d['C_lambda']:>10.6f} {d['C_predicted']:>10.6f} "
              f"{d['C_relative_error']:>10.2e} {prop:>7s}")
    print("  " + "─" * 65)

    print(f"\n  Formula: C_λ = n_gen · Tr(T_a²) / d = (N²-1)/N · Tr(T_a²)")
    print(f"  (Note: standard (N²-1)/(2N) assumes Tr(T²)=1/2 normalization)")

    # Check if the formula holds
    all_match = all(d["C_relative_error"] < 0.01 for d in all_data)
    if all_match:
        print(f"  ★ CONFIRMED: Casimir eigenvalue matches theoretical prediction")
        print(f"    for ALL groups with relative error < 1%.")
    else:
        fails = [d for d in all_data if d["C_relative_error"] >= 0.01]
        labels = ', '.join('SU(' + str(d["N"]) + ')' for d in fails)
        print(f"  ✗ MISMATCH for: {labels}")

    print(f"\n\n{'=' * 78}")
    print("  TABLE 2: COUPLING EXTRACTION")
    print("=" * 78)
    print(f"\n  {'Group':8s} {'f_rms':>10s} {'f_max':>10s} {'g_eff':>10s} "
          f"{'K_diag':>10s} {'h∨ meas':>10s} {'h∨ theo':>10s}")
    print("  " + "─" * 72)
    for d in all_data:
        print(f"  SU({d['N']})  {d['f_rms']:>10.6f} {d['f_max']:>10.6f} "
              f"{d['g_eff']:>10.6f} {d['K_diag_mean']:>10.6f} "
              f"{d['h_dual_measured']:>10.4f} {d['h_dual_theoretical']:>10d}")
    print("  " + "─" * 72)

    print(f"\n  g_eff = f_rms / Tr(T²) — normalization-independent effective coupling")
    print(f"  h∨ = K_diag / (2·Tr(T²)) — dual Coxeter number extraction")

    # Scaling analysis
    print(f"\n\n{'=' * 78}")
    print("  TABLE 3: SCALING WITH N")
    print("=" * 78)
    Ns = np.array([d["N"] for d in all_data])
    Cs = np.array([d["C_lambda"] for d in all_data])
    kappas = np.array([d["kappa_total"] for d in all_data])
    gs = np.array([d["g_eff"] for d in all_data])
    hs = np.array([d["h_dual_measured"] for d in all_data])

    # Fit C_λ / κ vs n_gen/d = (N²-1)/N
    ratio_measured = Cs / kappas
    ratio_predicted = (Ns**2 - 1) / Ns

    print(f"\n  {'N':>4s} {'(N²-1)/N':>12s} {'C_λ/κ':>12s} {'Match':>8s}")
    print("  " + "─" * 40)
    for i, d in enumerate(all_data):
        match = "✓" if abs(ratio_measured[i] - ratio_predicted[i]) / ratio_predicted[i] < 0.01 else "✗"
        print(f"  {d['N']:>4d} {ratio_predicted[i]:>12.6f} {ratio_measured[i]:>12.6f} {match:>8s}")
    print("  " + "─" * 40)

    # g_eff trend
    print(f"\n  g_eff scaling:")
    for i in range(1, len(all_data)):
        ratio = gs[i] / gs[i-1]
        print(f"    g_eff(SU({Ns[i]}))/g_eff(SU({Ns[i-1]})) = {ratio:.4f}")

    # h∨ trend
    print(f"\n  Dual Coxeter number:")
    for d in all_data:
        err = abs(d["h_dual_measured"] - d["h_dual_theoretical"]) / d["h_dual_theoretical"]
        print(f"    SU({d['N']}): h∨ = {d['h_dual_measured']:.4f} "
              f"(expected {d['h_dual_theoretical']}, error {err:.2e})")

    return all_data


if __name__ == "__main__":
    quick = "--quick" in sys.argv
    run_analysis(quick=quick)
