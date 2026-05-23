"""
Test 3.1: Is T_nc (non-commutativity tension) necessary?
Evolve with ONLY T_Killing (alpha >> 0, no T_nc contribution).

Question: Without T_nc, does it collapse to commutative (trivial)?

EXPECTED: Without T_nc, collapses to commuting operators — trivial algebra.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils import (
    evolve_generators,
    killing_metric,
    commutator,
    check_in_span,
    non_commutativity_tension,
)


def evolve_only_killing(n_gen=3, dim=2, steps=10000, lr=0.001, killing_target=-1.0, seed=42):
    """
    Custom optimizer: minimize ONLY T_Killing + T_norm (no T_nc).
    alpha for Killing, beta for norm, but T_nc is explicitly excluded.
    """
    np.random.seed(seed)
    from shared_utils import random_generators

    generators = random_generators(n=n_gen, dim=dim)
    for i in range(len(generators)):
        norm = np.linalg.norm(generators[i])
        if norm > 1e-10:
            generators[i] = generators[i] / norm

    K_target = killing_target * np.eye(n_gen)
    alpha = 5.0
    beta = 1.0
    history = {"T_k": [], "T_norm": [], "T_nc_observed": []}

    for step in range(steps):
        K = killing_metric(generators)
        T_k = np.linalg.norm(K - K_target) ** 2
        T_norm = sum((np.linalg.norm(G) - 1.0) ** 2 for G in generators)

        # Observe T_nc but do NOT include in loss
        T_nc = non_commutativity_tension(generators)

        T_total = alpha * T_k + beta * T_norm  # <-- NO T_nc

        history["T_k"].append(T_k)
        history["T_norm"].append(T_norm)
        history["T_nc_observed"].append(T_nc)

        # Numerical gradient of T_total (without T_nc)
        grad = []
        epsilon = 1e-6
        for i, G in enumerate(generators):
            grad_G = np.zeros_like(G)
            for r in range(G.shape[0]):
                for c in range(G.shape[1]):
                    G_plus = G.copy()
                    G_plus[r, c] += epsilon
                    gens_plus = [g.copy() for g in generators]
                    gens_plus[i] = G_plus
                    K_p = killing_metric(gens_plus)
                    T_p = (alpha * np.linalg.norm(K_p - K_target) ** 2
                           + beta * sum((np.linalg.norm(Gp) - 1.0) ** 2 for Gp in gens_plus))

                    G_minus = G.copy()
                    G_minus[r, c] -= epsilon
                    gens_minus = [g.copy() for g in generators]
                    gens_minus[i] = G_minus
                    K_m = killing_metric(gens_minus)
                    T_m = (alpha * np.linalg.norm(K_m - K_target) ** 2
                           + beta * sum((np.linalg.norm(Gm) - 1.0) ** 2 for Gm in gens_minus))

                    grad_G[r, c] = (T_p - T_m) / (2 * epsilon)
            grad.append(grad_G)

        for i in range(len(generators)):
            generators[i] = generators[i] - lr * grad[i]
            generators[i] = (generators[i] + generators[i].conj().T) / 2
            generators[i] -= np.trace(generators[i]) * np.eye(dim) / dim

        if step % 2000 == 0:
            print(f"  Step {step}: T_k={T_k:.6f}, T_norm={T_norm:.6f}, "
                  f"T_nc(observed)={T_nc:.6f}")

    return generators, history


def _killing_worker(args):
    """Worker for parallel evolve_only_killing."""
    seed, steps = args
    return evolve_only_killing(steps=steps, seed=seed)


def run_test_only_killing(n_seeds=5, steps=10000):
    print("=" * 70)
    print("TEST 3.1: EVOLVE WITH ONLY T_Killing (no T_nc)")
    print("=" * 70)

    from concurrent.futures import ProcessPoolExecutor
    import os

    tasks = [(seed, steps) for seed in range(n_seeds)]
    max_workers = min(n_seeds, os.cpu_count() or 4)
    print(f"  Running {n_seeds} seeds in parallel...")

    thread_vars = ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS')
    old_env = {v: os.environ.get(v) for v in thread_vars}
    for v in thread_vars:
        os.environ[v] = '1'

    try:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            raw = list(executor.map(_killing_worker, tasks))
    finally:
        for v, val in old_env.items():
            if val is None:
                os.environ.pop(v, None)
            else:
                os.environ[v] = val

    results = []
    for seed, (final_gen, hist) in enumerate(raw):
        K = killing_metric(final_gen)
        T_nc_final = non_commutativity_tension(final_gen)

        closes_all = all(
            check_in_span(commutator(final_gen[i], final_gen[j]), final_gen)
            for i, j in [(0, 1), (1, 2), (2, 0)]
        )

        comm_norms = []
        for i, j in [(0, 1), (1, 2), (2, 0)]:
            comm_norms.append(np.linalg.norm(commutator(final_gen[i], final_gen[j])))

        results.append({
            "seed": seed,
            "K_diag": np.mean(np.diag(K)),
            "T_nc_final": T_nc_final,
            "closes": closes_all,
            "comm_norms": comm_norms,
        })

        print(f"  seed={seed}: K_diag={np.mean(np.diag(K)):.4f}, T_nc={T_nc_final:.6f}, "
              f"comm_norms={[f'{c:.4f}' for c in comm_norms]}")

    # ---- Verdict ----
    print("\n--- VERDICT ---")
    avg_T_nc = np.mean([r["T_nc_final"] for r in results])
    has_su2 = any(r["closes"] and r["T_nc_final"] > 0.01 for r in results)

    if avg_T_nc < 1e-4:
        print("PASS (as expected): Without T_nc, generators collapse to commutative.")
        print("  This confirms T_nc is necessary to drive non-commutativity.")
    elif has_su2:
        print("SURPRISING: SU(2) still emerges even without T_nc!")
        print("  This would mean Killing metric alone forces non-commutativity.")
    else:
        print(f"PARTIAL: T_nc is reduced but not zero (avg={avg_T_nc:.6f}).")

    return results


if __name__ == "__main__":
    run_test_only_killing(n_seeds=5, steps=10000)
