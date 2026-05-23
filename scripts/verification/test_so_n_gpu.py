"""
Test A.2 + A.3: SO(n) Lie algebra emergence from Killing metric optimization.

Tests:
  - SO(3): 3 generators in dim=3 (cross-check with SU(2) ≅ SO(3))
  - SO(4): 6 generators in dim=4
  - SO(5): 10 generators in dim=5
  - SO(6): 15 generators in dim=6 (cross-check with SU(4) — same n_gen)

Usage:
    python test_so_n_gpu.py                 # run SO(3..6) sweep
    python test_so_n_gpu.py --groups 3      # only SO(3)
    python test_so_n_gpu.py --groups 3,6    # SO(3) and SO(6)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils_gpu import (
    params_to_antisymmetric,
    antisymmetric_to_params,
    random_params_so,
    evolve_generators_generic,
    killing_metric_np,
    structure_constants_general_np,
    compute_casimir_np,
    check_in_span_np,
    commutator_np,
    get_device_info,
)


def n_generators_so(n):
    """Number of generators for SO(n)."""
    return n * (n - 1) // 2


def reference_so_generators(n):
    """Build standard antisymmetric basis for so(n).

    E_{ij} has 1 at (i,j), -1 at (j,i), zeros elsewhere.
    Returns list of n(n-1)/2 numpy arrays, each (n, n) complex.
    """
    gens = []
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n), dtype=complex)
            E[i, j] = 1.0
            E[j, i] = -1.0
            gens.append(E)
    return gens


def analyze_so_result(generators, n, seed_idx, verbose=True):
    """Full analysis of an optimized SO(n) generator set.

    Returns dict with all measured quantities.
    """
    n_gen = len(generators)
    dim = generators[0].shape[0]
    expected_n_gen = n_generators_so(n)

    # Antisymmetry check: G + G^T = 0
    antisym_err = max(np.linalg.norm(g + g.T) for g in generators)

    # Reality check: Im(G) ≈ 0 (SO generators should be real)
    reality_err = max(np.max(np.abs(g.imag)) for g in generators)

    # Killing metric
    K = killing_metric_np(generators)
    K_diag = np.diag(K)
    K_diag_mean = np.mean(K_diag)
    K_diag_std = np.std(K_diag)
    K_offdiag = np.abs(K - np.diag(K_diag)).max()

    # Closure check: [T_i, T_j] ∈ span({T_k}) for all i,j
    closes_count = 0
    total_pairs = n_gen * (n_gen - 1) // 2
    max_closure_err = 0
    for i in range(n_gen):
        for j in range(i + 1, n_gen):
            C = commutator_np(generators[i], generators[j])
            if np.linalg.norm(C) < 1e-10:
                closes_count += 1
            elif check_in_span_np(C, generators, tol=1e-4):
                closes_count += 1
            else:
                residual = np.linalg.norm(C)
                max_closure_err = max(max_closure_err, residual)
    closure_frac = closes_count / total_pairs if total_pairs > 0 else 1.0

    # Structure constants (general convention — no factor of i for SO)
    f = structure_constants_general_np(generators)
    f_rms = np.sqrt(np.mean(f**2))
    f_max = np.max(np.abs(f))

    # Jacobi identity: J_{abc} = f_{abd}f_{dce} + cyclic
    jacobi_max = 0.0
    for a in range(min(n_gen, 6)):
        for b in range(a + 1, min(n_gen, 6)):
            for c in range(b + 1, min(n_gen, 6)):
                j_val = sum(
                    f[a, b, d] * f[d, c, :].sum() +
                    f[b, c, d] * f[d, a, :].sum() +
                    f[c, a, d] * f[d, b, :].sum()
                    for d in range(n_gen)
                )
                jacobi_max = max(jacobi_max, abs(j_val))

    # Casimir
    C_mat = compute_casimir_np(generators)
    C_diag = np.diag(C_mat).real
    C_lambda = np.mean(C_diag)
    C_spread = np.std(C_diag)
    C_offdiag = np.max(np.abs(C_mat - np.diag(C_diag)))

    # Generator norm κ = Tr(T²) — negative for antisymmetric generators
    kappa = np.mean([np.trace(g @ g).real for g in generators])
    kappa_abs = abs(kappa)

    # Observable ratio: |K_diag|/κ² — representation-dependent
    # For SO(n) in vector rep: should give exactly (n-2) (empirical)
    # For SU(N) in fundamental: gives 2N
    k_over_kappa2 = abs(K_diag_mean) / kappa_abs**2 if kappa_abs > 1e-12 else 0.0
    h_dual_theoretical = n - 2

    result = {
        "N": n,
        "n_gen": n_gen,
        "n_gen_expected": expected_n_gen,
        "dim": dim,
        "seed": seed_idx,
        "antisym_err": antisym_err,
        "reality_err": reality_err,
        "K_diag_mean": K_diag_mean,
        "K_diag_std": K_diag_std,
        "K_offdiag": K_offdiag,
        "closure_frac": closure_frac,
        "closes_count": closes_count,
        "total_pairs": total_pairs,
        "max_closure_err": max_closure_err,
        "f_rms": f_rms,
        "f_max": f_max,
        "jacobi_max": jacobi_max,
        "C_lambda": C_lambda,
        "C_spread": C_spread,
        "C_offdiag": C_offdiag,
        "kappa": kappa,
        "k_over_kappa2": k_over_kappa2,
        "h_dual_theoretical": h_dual_theoretical,
    }

    if verbose:
        tag = "PASS" if closure_frac > 0.99 else "PARTIAL" if closure_frac > 0.5 else "FAIL"
        print(f"\n  seed={seed_idx} [{tag}]:")
        print(f"    Antisymmetry err: {antisym_err:.2e}")
        print(f"    Reality err:      {reality_err:.2e}")
        print(f"    K_diag:           {K_diag_mean:.6f} +/- {K_diag_std:.6f}")
        print(f"    K_offdiag max:    {K_offdiag:.2e}")
        print(f"    Closure:          {closes_count}/{total_pairs} ({closure_frac:.1%})")
        print(f"    f_rms:            {f_rms:.6f}")
        print(f"    Jacobi max:       {jacobi_max:.2e}")
        print(f"    C_lambda:         {C_lambda:.6f} (spread: {C_spread:.2e})")
        print(f"    kappa=Tr(T^2):    {kappa:.6f}")
        print(f"    |K|/kappa^2:      {k_over_kappa2:.4f} (h_dual_theo: {h_dual_theoretical})")

    return result


def run_test_so(n, n_seeds=3, steps=None, lr=None):
    """Run SO(n) optimizer and analyze results."""
    n_gen = n_generators_so(n)
    dim = n

    if steps is None:
        steps = max(20000, 5000 * n)
    if lr is None:
        lr = max(0.0001, 0.001 / n)

    print(f"\n{'=' * 78}")
    print(f"  SO({n}): {n_gen} generators in dim={dim}")
    print(f"  steps={steps}, lr={lr}, seeds={n_seeds}")
    print(f"{'=' * 78}")

    t0 = time.time()
    raw = evolve_generators_generic(
        n_gen=n_gen,
        dim=dim,
        steps=steps,
        lr=lr,
        alpha=5.0,
        beta=1.0,
        killing_target=-1.0,
        batch_size=n_seeds,
        seeds=list(range(n_seeds)),
        verbose=True,
        convergence_threshold=1e-4,
        param_to_gen_fn=params_to_antisymmetric,
        gen_to_param_fn=antisymmetric_to_params,
        random_param_fn=lambda bs, ng, device, seeds: random_params_so(bs, ng, dim, device, seeds),
        param_to_gen_kwargs={"dim": dim},
    )
    gpu_time = time.time() - t0

    results = []
    for seed_idx, (gens_np, hist) in enumerate(raw):
        gen_list = list(gens_np)
        r = analyze_so_result(gen_list, n, seed_idx)
        r["gpu_time"] = gpu_time / n_seeds
        results.append(r)

    # Summary
    any_pass = any(r["closure_frac"] > 0.99 for r in results)
    avg_closure = np.mean([r["closure_frac"] for r in results])
    avg_K = np.mean([r["K_diag_mean"] for r in results])
    avg_ratio = np.mean([r["k_over_kappa2"] for r in results])
    h_theo = n - 2

    print(f"\n  SUMMARY SO({n}):")
    print(f"    GPU time:          {gpu_time:.1f}s ({gpu_time/n_seeds:.1f}s/seed)")
    print(f"    Avg closure:       {avg_closure:.1%}")
    print(f"    Avg K_diag:        {avg_K:.6f}")
    print(f"    |K|/kappa^2:       {avg_ratio:.4f}  (h_dual = {h_theo})")
    print(f"    Verdict:           {'PASS' if any_pass else 'FAIL'}")

    return {
        "group": f"SO({n})",
        "N": n,
        "n_gen": n_gen,
        "results": results,
        "gpu_time": gpu_time,
        "verdict": "PASS" if any_pass else "FAIL",
    }


def run_test_so_sweep(groups=None, n_seeds=3):
    """Run SO(n) sweep across multiple groups."""
    if groups is None:
        groups = [3, 4, 5, 6]

    print("=" * 78)
    print("  SO(n) SWEEP — Killing Metric Optimization")
    print("=" * 78)
    print(f"  Device: {get_device_info()}")
    print(f"  Groups: {', '.join(f'SO({n})' for n in groups)}")
    print(f"  Seeds per group: {n_seeds}")

    all_results = []
    t_total = time.time()

    for n in groups:
        r = run_test_so(n, n_seeds=n_seeds)
        all_results.append(r)

    t_total = time.time() - t_total

    # Comparative summary
    print(f"\n\n{'=' * 78}")
    print(f"  COMPARATIVE SUMMARY — SO(n) sweep")
    print(f"{'=' * 78}")
    print(f"\n  {'Group':8s} {'n_gen':>6s} {'Closure':>9s} {'K_diag':>10s} "
          f"{'|K|/k^2':>10s} {'h_theo':>10s} {'C_lam':>10s} {'Time':>8s} {'Verdict':>8s}")
    print("  " + "─" * 85)
    for r in all_results:
        avg = lambda key: np.mean([s[key] for s in r["results"]])
        print(f"  {r['group']:8s} {r['n_gen']:>6d} {avg('closure_frac'):>8.1%} "
              f"{avg('K_diag_mean'):>10.6f} {avg('k_over_kappa2'):>10.4f} "
              f"{int(avg('h_dual_theoretical')):>10d} {avg('C_lambda'):>10.6f} "
              f"{r['gpu_time']:>7.1f}s {r['verdict']:>8s}")
    print("  " + "─" * 85)
    print(f"\n  Total wall clock: {t_total:.1f}s")

    # Cross-check: SO(3) ≅ SU(2)?
    so3 = next((r for r in all_results if r["N"] == 3), None)
    if so3:
        avg_k = np.mean([s["K_diag_mean"] for s in so3["results"]])
        avg_c = np.mean([s["C_lambda"] for s in so3["results"]])
        print(f"\n  Cross-check SO(3) ≅ SU(2):")
        print(f"    SO(3) K_diag = {avg_k:.6f}")
        print(f"    SU(2) K_diag = -0.960925  (from prior runs)")
        print(f"    SO(3) C_λ    = {avg_c:.6f}")
        print(f"    SU(2) C_λ    = 0.735201   (from prior runs)")

    # Cross-check: SO(6) vs SU(4) — same n_gen=15
    so6 = next((r for r in all_results if r["N"] == 6), None)
    if so6:
        avg_k = np.mean([s["K_diag_mean"] for s in so6["results"]])
        avg_cl = np.mean([s["closure_frac"] for s in so6["results"]])
        print(f"\n  Cross-check SO(6) ~ SU(4) (same n_gen=15, different dim):")
        print(f"    SO(6) K_diag = {avg_k:.6f}  (dim=6)")
        print(f"    SU(4) K_diag = -0.962575   (dim=4, from prior runs)")

    return all_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--groups", default="3,4,5,6",
                        help="Comma-separated list of SO(n) groups")
    parser.add_argument("--seeds", type=int, default=3)
    args = parser.parse_args()

    groups = [int(x) for x in args.groups.split(",")]
    run_test_so_sweep(groups=groups, n_seeds=args.seeds)
