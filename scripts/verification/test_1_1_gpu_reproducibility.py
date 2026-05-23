"""
Test 1.1 GPU: Reproducibility across 100 seeds — batched on GPU.
Question: Does the optimizer consistently converge to SU(2) regardless of seed?

GPU version: all 100 seeds run simultaneously in a single batched tensor
using torch.autograd for analytical gradients.

EXPECTED: >=95/100 converge, all close, K_diag ≈ -1.0 ± small
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils_gpu import (
    evolve_generators_gpu_individual,
    killing_metric_np,
    commutator_np,
    check_in_span_np,
    get_device_info,
)


def run_test_reproducibility_gpu(n_seeds=100, steps=15000):
    print("=" * 70)
    print("TEST 1.1 GPU: REPRODUCIBILITY (100 seeds, batched on GPU)")
    print("=" * 70)
    print(f"  Device: {get_device_info()}")
    print(f"  Running {n_seeds} seeds in a single GPU batch...")
    print()

    t0 = time.time()
    raw = evolve_generators_gpu_individual(
        n_gen=3, dim=2, steps=steps,
        lr=0.001, alpha=5.0, beta=1.0,
        killing_target=-1.0,
        batch_size=n_seeds,
        seeds=list(range(n_seeds)),
        verbose=True,
    )
    gpu_time = time.time() - t0
    print(f"\n  GPU batch time: {gpu_time:.1f}s")

    results = []
    for seed, (final_gen, hist) in enumerate(raw):
        gens = list(final_gen)  # list of (dim, dim) arrays
        K = killing_metric_np(gens)
        K_diag_mean = np.mean(np.diag(K))
        K_offdiag_max = np.max(np.abs(K - np.diag(np.diag(K))))

        closes_all = all(
            check_in_span_np(commutator_np(gens[i], gens[j]), gens)
            for i, j in [(0, 1), (1, 2), (2, 0)]
        )

        T_final = hist["T_total"][-1] if hist["T_total"] else float('inf')

        results.append({
            "seed": seed,
            "K_diag": K_diag_mean,
            "K_offdiag": K_offdiag_max,
            "closes": closes_all,
            "T_final": T_final,
        })

    # ---- Analysis ----
    closes_count = sum(r["closes"] for r in results)
    K_diags = [r["K_diag"] for r in results]
    K_offdiags = [r["K_offdiag"] for r in results]
    T_finals = [r["T_final"] for r in results]

    print("\n--- RESULTS ---")
    print(f"  Closes SU(2): {closes_count}/{n_seeds}")
    print(f"  K_diag mean:  {np.mean(K_diags):.4f} ± {np.std(K_diags):.4f}")
    print(f"  K_offdiag max (mean): {np.mean(K_offdiags):.6f}")
    print(f"  T_final mean: {np.mean(T_finals):.6f} ± {np.std(T_finals):.6f}")

    failures = [r for r in results if not r["closes"]]
    if failures:
        print(f"\n  Failed seeds ({len(failures)}):")
        for f in failures[:10]:
            print(f"    seed={f['seed']}: closes={f['closes']}, "
                  f"K_diag={f['K_diag']:.4f}, T_final={f['T_final']:.6f}")

    # ---- Verdict ----
    print("\n--- VERDICT ---")
    if closes_count >= 95:
        print(f"  PASS: {closes_count}/{n_seeds} seeds produce a closed SU(2) algebra.")
    elif closes_count >= 80:
        print(f"  PARTIAL: {closes_count}/{n_seeds} close. Structure emerges but some seeds need tuning.")
    else:
        print(f"  FAIL: Only {closes_count}/{n_seeds} close. Optimizer is unreliable.")

    print(f"\n  [GPU batch: {gpu_time:.1f}s for {n_seeds} seeds × {steps} steps]")
    return results


if __name__ == "__main__":
    run_test_reproducibility_gpu(n_seeds=100, steps=15000)
