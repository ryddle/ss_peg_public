"""
Test 2.3: With 8 generators in dim=3, does SU(3) emerge?
Question: If the mechanism is universal, 8 generators + 3x3 matrices → SU(3)?

EXPECTED: Converges to SU(3) canonical (Gell-Mann matrices up to rotation).
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
    gell_mann_matrices,
)


def check_antisymmetry(f, n, tol=1e-4):
    """Check that f_ijk = -f_jik (antisymmetry in first two indices)."""
    violations = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if abs(f[i, j, k] + f[j, i, k]) > tol:
                    violations += 1
    return violations


def check_jacobi(generators, tol=1e-4):
    """Check Jacobi identity: [A,[B,C]] + [B,[C,A]] + [C,[A,B]] = 0"""
    n = len(generators)
    max_violation = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                term1 = commutator(generators[i], commutator(generators[j], generators[k]))
                term2 = commutator(generators[j], commutator(generators[k], generators[i]))
                term3 = commutator(generators[k], commutator(generators[i], generators[j]))
                violation = np.linalg.norm(term1 + term2 + term3)
                max_violation = max(max_violation, violation)
    return max_violation


def run_test_su3(n_seeds=5, steps=20000):
    print("=" * 70)
    print("TEST 2.3: 8 GENERATORS IN dim=3 — SU(3) CANDIDATE")
    print("=" * 70)
    print("(This test is computationally expensive due to 8 generators × 3x3)")
    print()

    # First: reference values from Gell-Mann
    gm = [l / 2 for l in gell_mann_matrices()]
    K_gm = killing_metric(gm)
    C_gm = compute_casimir(gm)
    print("Reference Gell-Mann/2:")
    print(f"  K_gm diagonal: {np.round(np.diag(K_gm), 3)}")
    print(f"  Casimir diagonal: {np.round(np.diag(C_gm).real, 3)}")
    print()

    tasks = [
        {'seed': seed, 'n_gen': 8, 'dim': 3, 'steps': steps,
         'lr': 0.0005, 'alpha': 5.0, 'beta': 1.0,
         'killing_target': -1.0, 'verbose': False,
         'convergence_threshold': 1e-4}
        for seed in range(n_seeds)
    ]
    print(f"  Running {n_seeds} seeds in parallel...")
    raw = evolve_parallel(tasks)

    results = []
    for seed, (final_gen, hist) in enumerate(raw):
        K = killing_metric(final_gen)
        K_diag = np.diag(K)
        K_offdiag_max = np.max(np.abs(K - np.diag(K_diag)))

        closure_checks = []
        for i in range(8):
            for j in range(i + 1, 8):
                comm = commutator(final_gen[i], final_gen[j])
                in_span = check_in_span(comm, final_gen, tol=1e-4)
                closure_checks.append(in_span)
        closes_all = all(closure_checks)
        closes_count = sum(closure_checks)

        jacobi_viol = check_jacobi(final_gen, tol=1e-4)

        f = structure_constants_real(final_gen)
        antisym_violations = check_antisymmetry(f, 8)

        C = compute_casimir(final_gen)
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
            "T_final": hist["T_total"][-1],
        })

        print(f"  seed={seed}: K_diag mean={np.mean(K_diag):.4f} ± {np.std(K_diag):.4f}")
        print(f"    Closure: {closes_count}/{len(closure_checks)} pairs")
        print(f"    Jacobi violation: {jacobi_viol:.6f}")
        print(f"    Casimir diagonal: {np.round(C_diag, 3)}")
        print()

    # ---- Verdict ----
    print("\n--- VERDICT ---")
    any_closed = any(r["closes_all"] for r in results)
    best_closure = max(r["closes_count"] for r in results)
    total_pairs = results[0]["total_pairs"]

    if any_closed:
        print(f"PASS: SU(3) algebra emerges! Full closure achieved.")
    elif best_closure > total_pairs * 0.8:
        print(f"PARTIAL: Close to SU(3) — {best_closure}/{total_pairs} pairs close. "
              f"May need more steps or LR tuning.")
    else:
        print(f"FAIL: SU(3) does NOT emerge reliably ({best_closure}/{total_pairs} pairs).")
        print("  Possible reasons: more steps needed, different LR, or mechanism is specific to SU(2).")

    return results


if __name__ == "__main__":
    run_test_su3(n_seeds=3, steps=20000)
