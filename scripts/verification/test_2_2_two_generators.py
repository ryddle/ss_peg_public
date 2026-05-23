"""
Test 2.2: What happens with only 2 generators?
Question: With 2 generators in dim=2, does U(1) (commutative) emerge?

EXPECTED: Generators should commute — no room for SU(2) with only 2 generators.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils import (
    evolve_parallel,
    killing_metric,
    commutator,
)


def run_test_2_generators(n_seeds=20, steps=15000):
    print("=" * 70)
    print("TEST 2.2: 2 GENERATORS IN dim=2")
    print("=" * 70)

    tasks = [
        {'seed': seed, 'n_gen': 2, 'dim': 2, 'steps': steps,
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

        comm = commutator(final_gen[0], final_gen[1])
        comm_norm = np.linalg.norm(comm)

        T_nc = hist["T_nc"][-1]

        results.append({
            "seed": seed,
            "K_diag": K_diag_mean,
            "comm_norm": comm_norm,
            "T_nc": T_nc,
            "T_final": hist["T_total"][-1],
        })

        print(f"  seed={seed:2d}: ||[G0,G1]||={comm_norm:.6f}, "
              f"K_diag={K_diag_mean:.4f}, T_nc={T_nc:.6f}")

    # ---- Analysis ----
    comm_norms = [r["comm_norm"] for r in results]
    commute_count = sum(1 for c in comm_norms if c < 1e-4)

    print("\n--- RESULTS ---")
    print(f"Commutator norms: mean={np.mean(comm_norms):.6f}, max={np.max(comm_norms):.6f}")
    print(f"Effectively commute (<1e-4): {commute_count}/{n_seeds}")

    # ---- Verdict ----
    print("\n--- VERDICT ---")
    if commute_count == n_seeds:
        print("PASS: 2 generators always commute — U(1)×U(1) or trivial, no room for SU(2).")
    elif np.mean(comm_norms) > 0.1:
        print("INTERESTING: 2 generators do NOT commute. They may form a partial SU(2) sub-algebra.")
        print("  This would mean the commutator [G0,G1] generates a third direction implicitly.")
    else:
        print("PARTIAL: Mixed behavior. Some seeds commute, some don't.")

    return results


if __name__ == "__main__":
    run_test_2_generators(n_seeds=20, steps=15000)
