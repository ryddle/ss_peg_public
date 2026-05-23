"""
Test 3.2: Is T_Killing necessary?
Evolve with ONLY T_nc + T_norm (alpha=0, no Killing regularization).

Question: Without Killing, does it converge to SU(2) or something else?

EXPECTED: Does NOT converge well — needs Killing to regularize.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils import (
    evolve_generators,
    evolve_parallel,
    killing_metric,
    commutator,
    check_in_span,
    non_commutativity_tension,
)


def run_test_only_nc(n_seeds=10, steps=15000):
    print("=" * 70)
    print("TEST 3.2: EVOLVE WITH ONLY T_nc (no T_Killing, alpha=0)")
    print("=" * 70)

    tasks = [
        {'seed': seed, 'n_gen': 3, 'dim': 2, 'steps': steps,
         'lr': 0.001, 'alpha': 0.0, 'beta': 1.0,
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

        T_nc = non_commutativity_tension(final_gen)

        K_is_proportional = K_offdiag_max < 0.1

        results.append({
            "seed": seed,
            "K_diag": K_diag_mean,
            "K_offdiag": K_offdiag_max,
            "closes": closes_all,
            "T_nc": T_nc,
            "K_proportional_to_I": K_is_proportional,
        })

        print(f"  seed={seed:2d}: K_diag={K_diag_mean:.4f}, K_offdiag={K_offdiag_max:.4f}, "
              f"closes={closes_all}, T_nc={T_nc:.4f}")

    # ---- Analysis ----
    closes_count = sum(r["closes"] for r in results)
    K_prop_count = sum(r["K_proportional_to_I"] for r in results)

    print("\n--- RESULTS ---")
    print(f"Closes SU(2): {closes_count}/{n_seeds}")
    print(f"K proportional to I: {K_prop_count}/{n_seeds}")
    print(f"K_diag mean: {np.mean([r['K_diag'] for r in results]):.4f} ± "
          f"{np.std([r['K_diag'] for r in results]):.4f}")
    print(f"K_offdiag mean: {np.mean([r['K_offdiag'] for r in results]):.4f}")

    # ---- Compare with full optimization ----
    print("\n--- Comparison: full vs no-Killing ---")
    np.random.seed(42)
    full_gen, _ = evolve_generators(
        n_gen=3, dim=2, steps=steps, lr=0.001,
        alpha=5.0, beta=1.0, killing_target=-1.0, verbose=False,
    )
    K_full = killing_metric(full_gen)
    print(f"  Full (alpha=5): K_diag={np.mean(np.diag(K_full)):.4f}, "
          f"K_offdiag={np.max(np.abs(K_full - np.diag(np.diag(K_full)))):.6f}")
    print(f"  No-Killing (alpha=0): K_diag={results[0]['K_diag']:.4f}, "
          f"K_offdiag={results[0]['K_offdiag']:.4f}")

    # ---- Verdict ----
    print("\n--- VERDICT ---")
    if closes_count < n_seeds * 0.5:
        print("PASS (as expected): Without Killing regularization, SU(2) does NOT emerge reliably.")
        print("  Killing metric is necessary to guide convergence.")
    elif closes_count == n_seeds and K_prop_count == n_seeds:
        print("SURPRISING: SU(2) emerges even without Killing!")
        print("  T_nc alone is sufficient. Killing is just acceleration.")
    else:
        print(f"PARTIAL: SU(2) emerges sometimes ({closes_count}/{n_seeds}) but Killing helps.")

    return results


if __name__ == "__main__":
    run_test_only_nc(n_seeds=10, steps=15000)
