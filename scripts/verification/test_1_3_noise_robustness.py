"""
Test 1.3: Robustness to perturbations.
Question: Starting from Pauli + noise, does the optimizer return to SU(2)?

EXPECTED: Always returns to SU(2), even with large noise.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils import (
    pauli_matrices,
    random_hermitian,
    evolve_parallel,
    killing_metric,
    commutator,
    check_in_span,
    structure_constants_real,
)


def measure_distance_to_pauli(generators):
    """
    Measure distance to Pauli/2 basis (up to unitary rotation).
    Uses the Killing metric as a basis-independent invariant.
    """
    pauli = [s / 2 for s in pauli_matrices()]
    K_test = killing_metric(generators)
    K_pauli = killing_metric(pauli)
    return np.linalg.norm(K_test - K_pauli)


def run_test_noise_robustness(steps=10000):
    print("=" * 70)
    print("TEST 1.3: ROBUSTNESS TO PERTURBATIONS")
    print("=" * 70)

    pauli_gens = [s / 2 for s in pauli_matrices()]
    noise_levels = [0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0]

    # Pre-create noisy initial generators for each noise level
    tasks = []
    for noise in noise_levels:
        np.random.seed(42)
        noisy_gens = [G + noise * random_hermitian(2) for G in pauli_gens]
        for i in range(len(noisy_gens)):
            noisy_gens[i] = (noisy_gens[i] + noisy_gens[i].conj().T) / 2
            noisy_gens[i] -= np.trace(noisy_gens[i]) * np.eye(2) / 2
        tasks.append({
            'initial_generators': noisy_gens,
            'steps': steps, 'lr': 0.001, 'alpha': 5.0, 'beta': 1.0,
            'killing_target': -1.0, 'verbose': False,
        })

    print(f"  Running {len(noise_levels)} noise levels in parallel...")
    raw = evolve_parallel(tasks)

    results = []
    for noise, (final_gen, hist) in zip(noise_levels, raw):
        K = killing_metric(final_gen)
        K_diag_mean = np.mean(np.diag(K))
        K_offdiag_max = np.max(np.abs(K - np.diag(np.diag(K))))

        closes_all = all(
            check_in_span(commutator(final_gen[i], final_gen[j]), final_gen)
            for i, j in [(0, 1), (1, 2), (2, 0)]
        )

        dist_pauli = measure_distance_to_pauli(final_gen)

        results.append({
            "noise": noise,
            "K_diag": K_diag_mean,
            "K_offdiag": K_offdiag_max,
            "closes": closes_all,
            "dist_pauli": dist_pauli,
            "T_final": hist["T_total"][-1],
        })

        print(f"  noise={noise:.2f}: K_diag={K_diag_mean:.4f}, closes={closes_all}, "
              f"dist_pauli={dist_pauli:.4f}")

    # ---- Verdict ----
    print("\n--- VERDICT ---")
    all_close = all(r["closes"] for r in results)
    if all_close:
        print("PASS: Always returns to SU(2) even with large noise.")
    else:
        failed_noise = [r["noise"] for r in results if not r["closes"]]
        print(f"FAIL: Does NOT return to SU(2) for noise levels: {failed_noise}")

    return results


if __name__ == "__main__":
    run_test_noise_robustness(steps=10000)
