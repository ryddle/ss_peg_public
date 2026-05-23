"""
Test 2.1: Does the mechanism work in dim=3 (3x3 matrices)?
Question: With 3 generators in dim=3, does SU(2) emerge embedded in 3x3?

EXPECTED: SU(2) sub-algebra emerges (spin-1 representation) embedded in 3x3.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils import (
    evolve_parallel,
    killing_metric,
    commutator,
    check_in_span,
    structure_constants_real,
    compute_casimir,
)


def run_test_dim3(n_seeds=10, steps=15000):
    print("=" * 70)
    print("TEST 2.1: 3 GENERATORS IN dim=3 (3x3 MATRICES)")
    print("=" * 70)

    tasks = [
        {'seed': seed, 'n_gen': 3, 'dim': 3, 'steps': steps,
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

        f = structure_constants_real(final_gen)
        f_012 = f[0, 1, 2]
        f_120 = f[1, 2, 0]
        f_201 = f[2, 0, 1]

        C = compute_casimir(final_gen)
        C_diag = np.diag(C).real

        results.append({
            "seed": seed,
            "K_diag": K_diag_mean,
            "K_offdiag": K_offdiag_max,
            "closes": closes_all,
            "f_012": f_012,
            "f_120": f_120,
            "f_201": f_201,
            "casimir_diag": C_diag,
            "T_final": hist["T_total"][-1],
        })

        print(f"  seed={seed}: K_diag={K_diag_mean:.4f}, closes={closes_all}, "
              f"f_{{012}}={f_012:.3f}, Casimir_diag={np.round(C_diag, 3)}")

    # ---- Analysis ----
    closes_count = sum(r["closes"] for r in results)

    print("\n--- RESULTS ---")
    print(f"Closes SU(2): {closes_count}/{n_seeds}")
    print(f"K_diag mean: {np.mean([r['K_diag'] for r in results]):.4f}")

    # Compare with spin-1 representation of SU(2)
    # Spin-1 Casimir = j(j+1) = 1*2 = 2, so C = 2*I_3
    print(f"\nCasimir comparison:")
    print(f"  Expected for spin-1 (j=1): C = 2*I_3 (diagonal = [2, 2, 2])")
    print(f"  Actual Casimir diagonals (first 3 seeds):")
    for r in results[:3]:
        print(f"    seed={r['seed']}: {np.round(r['casimir_diag'], 3)}")

    # ---- Verdict ----
    print("\n--- VERDICT ---")
    if closes_count >= n_seeds - 1:
        print("PASS: SU(2) emerges in dim=3 (likely spin-1 representation).")
    else:
        print("PARTIAL/FAIL: Closure not consistent in dim=3.")

    return results


if __name__ == "__main__":
    run_test_dim3(n_seeds=10, steps=15000)
