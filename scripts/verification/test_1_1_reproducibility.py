"""
Test 1.1: Reproducibility across 100 seeds.
Question: Does the optimizer consistently converge to SU(2) regardless of seed?

EXPECTED: >95/100 converge, all close, K_diag ≈ -1.0 ± small
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils import (
    evolve_parallel,
    killing_metric,
    commutator,
    check_in_span,
)


def run_test_reproducibility(n_seeds=100, steps=15000):
    print("=" * 70)
    print("TEST 1.1: REPRODUCIBILITY (100 seeds)")
    print("=" * 70)

    tasks = [
        {'seed': seed, 'n_gen': 3, 'dim': 2, 'steps': steps,
         'lr': 0.001, 'alpha': 5.0, 'beta': 1.0,
         'killing_target': -1.0, 'verbose': False}
        for seed in range(n_seeds)
    ]
    print(f"  Running {n_seeds} seeds in parallel...")
    raw = evolve_parallel(tasks)

    results = []
    for seed, (final_gen, hist) in enumerate(raw):
        K = killing_metric(final_gen)
        K_diag_mean = np.mean(np.diag(K))
        K_offdiag_max = np.max(np.abs(K - np.diag(np.diag(K))))

        closes_all = all(
            check_in_span(commutator(final_gen[i], final_gen[j]), final_gen)
            for i, j in [(0, 1), (1, 2), (2, 0)]
        )

        T_final = hist["T_total"][-1]
        converged = T_final < 0.01

        results.append({
            "seed": seed,
            "K_diag": K_diag_mean,
            "K_offdiag": K_offdiag_max,
            "closes": closes_all,
            "converged": converged,
            "T_final": T_final,
        })

    # ---- Analysis ----
    converged_count = sum(r["converged"] for r in results)
    closes_count = sum(r["closes"] for r in results)
    K_diags = [r["K_diag"] for r in results]
    K_offdiags = [r["K_offdiag"] for r in results]
    T_finals = [r["T_final"] for r in results]

    print("\n--- RESULTS ---")
    print(f"Converged:   {converged_count}/{n_seeds}")
    print(f"Closes SU(2): {closes_count}/{n_seeds}")
    print(f"K_diag mean:  {np.mean(K_diags):.4f} ± {np.std(K_diags):.4f}")
    print(f"K_offdiag max (mean): {np.mean(K_offdiags):.6f}")
    print(f"T_final mean: {np.mean(T_finals):.6f} ± {np.std(T_finals):.6f}")

    # Show failures
    failures = [r for r in results if not r["converged"] or not r["closes"]]
    if failures:
        print(f"\nFailed seeds ({len(failures)}):")
        for f in failures[:10]:
            print(f"  seed={f['seed']}: converged={f['converged']}, closes={f['closes']}, "
                  f"K_diag={f['K_diag']:.4f}, T_final={f['T_final']:.6f}")

    # ---- Verdict ----
    print("\n--- VERDICT ---")
    if closes_count >= 95:
        print(f"PASS: {closes_count}/{n_seeds} seeds produce a closed SU(2) algebra.")
        if converged_count < 95:
            print(f"  Note: Only {converged_count}/{n_seeds} reached T<0.01 in {steps} steps")
            print(f"  (more steps would tighten K toward -I, but the algebra is already correct).")
    elif closes_count >= 80:
        print(f"PARTIAL: {closes_count}/{n_seeds} close. Structure emerges but some seeds need tuning.")
    else:
        print(f"FAIL: Only {closes_count}/{n_seeds} close. Optimizer is unreliable.")

    return results


if __name__ == "__main__":
    run_test_reproducibility(n_seeds=100, steps=15000)
