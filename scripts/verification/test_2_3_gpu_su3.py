"""
Test 2.3 GPU: 8 generators in dim=3 → SU(3) emergence.
Question: With GPU autograd, does SU(3) emerge more efficiently?

GPU version: Uses torch.autograd for analytical gradients (vs CPU finite differences).
8 generators × 3×3 = 72 complex parameters → huge speedup from autograd.

EXPECTED: All 28 commutator pairs close, Jacobi identity holds.
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
    gell_mann_matrices_torch,
    get_device_info,
)


def check_antisymmetry(f, n, tol=1e-4):
    violations = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if abs(f[i, j, k] + f[j, i, k]) > tol:
                    violations += 1
    return violations


def check_jacobi(generators, tol=1e-4):
    n = len(generators)
    max_violation = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                term1 = commutator_np(generators[i], commutator_np(generators[j], generators[k]))
                term2 = commutator_np(generators[j], commutator_np(generators[k], generators[i]))
                term3 = commutator_np(generators[k], commutator_np(generators[i], generators[j]))
                violation = np.linalg.norm(term1 + term2 + term3)
                max_violation = max(max_violation, violation)
    return max_violation


def run_test_su3_gpu(n_seeds=5, steps=20000):
    print("=" * 70)
    print("TEST 2.3 GPU: 8 GENERATORS IN dim=3 — SU(3) CANDIDATE")
    print("=" * 70)
    print(f"  Device: {get_device_info()}")
    print()

    # Reference: Gell-Mann/2
    gm_torch = gell_mann_matrices_torch()
    gm = [g.cpu().numpy() / 2 for g in gm_torch]
    K_gm = killing_metric_np(gm)
    C_gm = compute_casimir_np(gm)
    print("  Reference Gell-Mann/2:")
    print(f"    K_gm diagonal: {np.round(np.diag(K_gm), 3)}")
    print(f"    Casimir diagonal: {np.round(np.diag(C_gm).real, 3)}")
    print()

    print(f"  Running {n_seeds} seeds in GPU batch...")
    t0 = time.time()
    raw = evolve_generators_gpu_individual(
        n_gen=8, dim=3, steps=steps,
        lr=0.0005, alpha=5.0, beta=1.0,
        killing_target=-1.0,
        batch_size=n_seeds,
        seeds=list(range(n_seeds)),
        verbose=True,
        convergence_threshold=1e-4,
    )
    gpu_time = time.time() - t0
    print(f"\n  GPU batch time: {gpu_time:.1f}s")

    results = []
    for seed, (final_gen, hist) in enumerate(raw):
        gens = list(final_gen)
        K = killing_metric_np(gens)
        K_diag = np.diag(K)
        K_offdiag_max = np.max(np.abs(K - np.diag(K_diag)))

        closure_checks = []
        for i in range(8):
            for j in range(i + 1, 8):
                comm = commutator_np(gens[i], gens[j])
                in_span = check_in_span_np(comm, gens, tol=1e-4)
                closure_checks.append(in_span)
        closes_all = all(closure_checks)
        closes_count = sum(closure_checks)

        jacobi_viol = check_jacobi(gens)

        f = structure_constants_real_np(gens)
        antisym_violations = check_antisymmetry(f, 8)

        C = compute_casimir_np(gens)
        C_diag = np.diag(C).real

        results.append({
            "seed": seed,
            "K_diag_mean": np.mean(K_diag),
            "K_diag_std": np.std(K_diag),
            "K_offdiag_max": K_offdiag_max,
            "closes_all": closes_all,
            "closes_count": closes_count,
            "total_pairs": len(closure_checks),
            "jacobi_violation": jacobi_viol,
            "antisym_violations": antisym_violations,
            "casimir_diag": C_diag,
            "T_final": hist["T_total"][-1] if hist["T_total"] else float('inf'),
        })

        print(f"\n  seed={seed}: K_diag mean={np.mean(K_diag):.4f} ± {np.std(K_diag):.4f}")
        print(f"    Closure: {closes_count}/{len(closure_checks)} pairs")
        print(f"    Jacobi violation: {jacobi_viol:.6f}")
        print(f"    Casimir diagonal: {np.round(C_diag, 3)}")

    # ---- Verdict ----
    print("\n--- VERDICT ---")
    any_closed = any(r["closes_all"] for r in results)
    best_closure = max(r["closes_count"] for r in results)
    total_pairs = results[0]["total_pairs"]

    if any_closed:
        print(f"  PASS: SU(3) algebra emerges! Full closure achieved.")
    elif best_closure > total_pairs * 0.8:
        print(f"  PARTIAL: Close to SU(3) — {best_closure}/{total_pairs} pairs close.")
    else:
        print(f"  FAIL: SU(3) does NOT emerge reliably ({best_closure}/{total_pairs} pairs).")

    print(f"\n  [GPU batch: {gpu_time:.1f}s for {n_seeds} seeds × {steps} steps]")
    return results


if __name__ == "__main__":
    run_test_su3_gpu(n_seeds=5, steps=20000)
