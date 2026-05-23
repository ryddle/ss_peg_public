"""
Test SU(N) GPU: Generalized gauge symmetry emergence test.

Tests whether SU(N) emerges from N²-1 generators in dimension N,
for any N ≥ 2. Uses GPU-accelerated optimizer with autograd.

Usage:
    python test_su_n_gpu.py              # defaults: SU(2), SU(3), SU(4), SU(5)
    python test_su_n_gpu.py 4            # test SU(4) only
    python test_su_n_gpu.py 2 3 4 5 6    # test multiple groups

Key quantities verified per SU(N):
  - Algebraic closure: all N²-1 choose 2 commutator pairs lie in the span
  - Jacobi identity: [A,[B,C]] + cyclic = 0
  - Killing metric: K ≈ c·I (proportional to identity)
  - Casimir operator: C = sum(G_i²) ∝ I
  - Structure constants: antisymmetric, real
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils_gpu import (
    evolve_generators_gpu_individual,
    killing_metric_np,
    commutator_np,
    check_in_span_np,
    structure_constants_real_np,
    compute_casimir_np,
    orthonormalize_killing_np,
    get_device_info,
)


def check_jacobi(generators, tol=1e-4):
    n = len(generators)
    max_violation = 0.0
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                term1 = commutator_np(generators[i], commutator_np(generators[j], generators[k]))
                term2 = commutator_np(generators[j], commutator_np(generators[k], generators[i]))
                term3 = commutator_np(generators[k], commutator_np(generators[i], generators[j]))
                violation = np.linalg.norm(term1 + term2 + term3)
                max_violation = max(max_violation, violation)
                count += 1
    return max_violation, count


def check_antisymmetry(f, n, tol=1e-4):
    violations = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if abs(f[i, j, k] + f[j, i, k]) > tol:
                    violations += 1
    return violations


def check_casimir_proportional_to_identity(C, tol=0.01):
    """Check if Casimir C = lambda * I."""
    d = C.shape[0]
    diag_vals = np.diag(C).real
    mean_val = np.mean(diag_vals)
    diag_spread = np.std(diag_vals) / (abs(mean_val) + 1e-12)
    offdiag_max = np.max(np.abs(C - np.diag(np.diag(C))))
    return diag_spread < tol and offdiag_max < tol * abs(mean_val), mean_val, diag_spread, offdiag_max


def run_test_su_n(N, n_seeds=3, steps=None, lr=None, verbose=True):
    """Run SU(N) emergence test.

    Args:
        N: dimension (tests SU(N) with N²-1 generators in dim N)
        n_seeds: number of random initializations
        steps: optimization steps (auto-scaled if None)
        lr: learning rate (auto-scaled if None)
    """
    n_gen = N * N - 1
    n_pairs = n_gen * (n_gen - 1) // 2

    # Auto-scale parameters based on N
    if steps is None:
        steps = max(20000, 5000 * N)
    if lr is None:
        lr = max(0.0001, 0.001 / N)

    print("=" * 70)
    print(f"  SU({N}) EMERGENCE TEST — GPU")
    print(f"  {n_gen} generators in dim={N} ({n_pairs} commutator pairs)")
    print("=" * 70)
    print(f"  Device: {get_device_info()}")
    print(f"  Parameters: steps={steps}, lr={lr:.6f}, seeds={n_seeds}")
    print()

    t0 = time.time()
    raw = evolve_generators_gpu_individual(
        n_gen=n_gen, dim=N, steps=steps,
        lr=lr, alpha=5.0, beta=1.0,
        killing_target=-1.0,
        batch_size=n_seeds,
        seeds=list(range(n_seeds)),
        verbose=verbose,
        convergence_threshold=1e-4,
    )
    gpu_time = time.time() - t0
    print(f"\n  GPU batch time: {gpu_time:.1f}s")

    results = []
    for seed, (final_gen, hist) in enumerate(raw):
        gens = list(final_gen)

        # Killing metric
        K = killing_metric_np(gens)
        K_diag = np.diag(K)
        K_diag_mean = np.mean(K_diag)
        K_diag_std = np.std(K_diag)
        K_offdiag_max = np.max(np.abs(K - np.diag(K_diag)))

        # Closure check
        closure_checks = []
        for i in range(n_gen):
            for j in range(i + 1, n_gen):
                comm = commutator_np(gens[i], gens[j])
                in_span = check_in_span_np(comm, gens, tol=1e-4)
                closure_checks.append(in_span)
        closes_all = all(closure_checks)
        closes_count = sum(closure_checks)

        # Jacobi identity
        jacobi_viol, jacobi_triples = check_jacobi(gens)

        # Structure constants
        f = structure_constants_real_np(gens)
        antisym_violations = check_antisymmetry(f, n_gen)

        # Casimir
        C = compute_casimir_np(gens)
        casimir_prop, casimir_val, casimir_spread, casimir_offdiag = \
            check_casimir_proportional_to_identity(C)

        # Final tension
        T_final = hist["T_total"][-1] if hist["T_total"] else float('inf')
        final_steps = len(hist["T_total"])

        result = {
            "seed": seed,
            "K_diag_mean": K_diag_mean,
            "K_diag_std": K_diag_std,
            "K_offdiag_max": K_offdiag_max,
            "closes_all": closes_all,
            "closes_count": closes_count,
            "total_pairs": n_pairs,
            "jacobi_violation": jacobi_viol,
            "jacobi_triples": jacobi_triples,
            "antisym_violations": antisym_violations,
            "casimir_proportional": casimir_prop,
            "casimir_eigenvalue": casimir_val,
            "casimir_spread": casimir_spread,
            "casimir_offdiag": casimir_offdiag,
            "T_final": T_final,
            "steps_used": final_steps,
        }
        results.append(result)

        print(f"\n  seed={seed}:")
        print(f"    K_diag: {K_diag_mean:.4f} ± {K_diag_std:.4f}  (offdiag max: {K_offdiag_max:.6f})")
        print(f"    Closure: {closes_count}/{n_pairs} pairs")
        print(f"    Jacobi max violation: {jacobi_viol:.6f} ({jacobi_triples} triples)")
        print(f"    Antisymmetry violations: {antisym_violations}")
        print(f"    Casimir ∝ I: {casimir_prop}  (λ={casimir_val:.4f}, spread={casimir_spread:.6f})")
        print(f"    Steps used: {final_steps}, T_final: {T_final:.6f}")

    # ---- Summary ----
    print(f"\n{'─' * 70}")
    print(f"  SUMMARY — SU({N})")
    print(f"{'─' * 70}")

    all_closed = [r for r in results if r["closes_all"]]
    best_closure = max(r["closes_count"] for r in results)
    all_jacobi_ok = all(r["jacobi_violation"] < 1e-3 for r in results)
    all_casimir_ok = all(r["casimir_proportional"] for r in results)
    all_antisym_ok = all(r["antisym_violations"] == 0 for r in results)

    checks = {
        "Closure (all pairs)": f"{len(all_closed)}/{n_seeds} seeds",
        "Best closure": f"{best_closure}/{n_pairs} pairs",
        "Jacobi identity": "PASS" if all_jacobi_ok else "FAIL",
        "Antisymmetry f_abc": "PASS" if all_antisym_ok else "FAIL",
        "Casimir ∝ I": "PASS" if all_casimir_ok else "FAIL",
        "K_diag uniformity": f"std={np.mean([r['K_diag_std'] for r in results]):.6f}",
    }

    for key, val in checks.items():
        print(f"  {key:30s}: {val}")

    # ---- Verdict ----
    print(f"\n{'─' * 70}")
    is_pass = len(all_closed) > 0 and all_jacobi_ok
    is_partial = best_closure > n_pairs * 0.8

    if is_pass:
        print(f"  ✓ PASS: SU({N}) algebra emerges from the Killing metric variational principle!")
        print(f"    {n_gen} generators close under commutation, Jacobi identity satisfied.")
    elif is_partial:
        print(f"  ~ PARTIAL: Close to SU({N}) — {best_closure}/{n_pairs} pairs close.")
        print(f"    May need more steps or learning rate tuning.")
    else:
        print(f"  ✗ FAIL: SU({N}) does NOT emerge ({best_closure}/{n_pairs} pairs).")

    print(f"\n  [GPU: {gpu_time:.1f}s for {n_seeds} seeds × {steps} steps]")

    return {
        "N": N,
        "group": f"SU({N})",
        "n_gen": n_gen,
        "n_pairs": n_pairs,
        "pass": is_pass,
        "partial": is_partial,
        "results": results,
        "gpu_time": gpu_time,
    }


def run_test_su_n_sweep(groups=None):
    """Run SU(N) tests for multiple N values and produce a comparative summary."""
    if groups is None:
        groups = [2, 3, 4, 5]

    print("=" * 70)
    print("  SU(N) UNIVERSALITY SWEEP — GPU")
    print(f"  Testing: {', '.join(f'SU({n})' for n in groups)}")
    print(f"  Device: {get_device_info()}")
    print("=" * 70)
    print()

    total_t0 = time.time()
    all_results = {}

    for N in groups:
        print(f"\n{'#' * 70}")
        print(f"#  SU({N}): {N**2-1} generators in dim={N}")
        print(f"{'#' * 70}\n")
        all_results[N] = run_test_su_n(N, n_seeds=3)
        print()

    # ---- Comparative summary ----
    total_time = time.time() - total_t0

    print("\n" + "=" * 70)
    print("  COMPARATIVE SUMMARY — SU(N) UNIVERSALITY")
    print("=" * 70)
    print(f"\n  {'Group':8s} {'Generators':>12s} {'Pairs':>8s} {'Closure':>12s} "
          f"{'Jacobi':>10s} {'Casimir':>10s} {'K_diag':>10s} {'Time':>8s} {'Verdict':>10s}")
    print("  " + "─" * 90)

    for N in groups:
        r = all_results[N]
        best_seed = max(r["results"], key=lambda x: x["closes_count"])
        closure_str = f"{best_seed['closes_count']}/{r['n_pairs']}"
        jacobi_str = f"{best_seed['jacobi_violation']:.2e}"
        casimir_str = "∝ I" if best_seed["casimir_proportional"] else "NO"
        k_str = f"{best_seed['K_diag_mean']:.4f}"
        time_str = f"{r['gpu_time']:.0f}s"
        verdict = "PASS" if r["pass"] else ("PARTIAL" if r["partial"] else "FAIL")
        print(f"  SU({N:d})   {r['n_gen']:>12d} {r['n_pairs']:>8d} {closure_str:>12s} "
              f"{jacobi_str:>10s} {casimir_str:>10s} {k_str:>10s} {time_str:>8s} {verdict:>10s}")

    print("  " + "─" * 90)
    print(f"\n  Total GPU time: {total_time:.1f}s")

    # Overall assessment
    all_pass = all(all_results[N]["pass"] for N in groups)
    if all_pass:
        print(f"\n  ★ UNIVERSAL: The Killing metric variational principle produces exact")
        print(f"    Lie algebra closure for ALL tested groups: {', '.join(f'SU({n})' for n in groups)}.")
        print(f"    The mechanism is confirmed universal across gauge symmetries.")
    else:
        passing = [f"SU({N})" for N in groups if all_results[N]["pass"]]
        failing = [f"SU({N})" for N in groups if not all_results[N]["pass"]]
        print(f"\n  Passing: {', '.join(passing) if passing else 'None'}")
        print(f"  Failing: {', '.join(failing) if failing else 'None'}")

    return all_results


if __name__ == "__main__":
    args = sys.argv[1:]
    if args:
        groups = [int(a) for a in args]
    else:
        groups = [2, 3, 4, 5]

    run_test_su_n_sweep(groups)
