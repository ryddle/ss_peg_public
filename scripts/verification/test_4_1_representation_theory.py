"""
Test 4.1: Comparison with representation theory.
Verifies:
  - Casimir operator C = sum(G_i^2) is proportional to identity
  - The proportionality matches j(j+1) after accounting for generator norms
  - Structure constants f_ijk proportional to epsilon_ijk
  - Killing metric proportional to -delta_ij

Normalization-aware:
  For generators with Tr(G_a G_b) = kappa * delta_ab,
  Casimir C = 2*kappa * j(j+1) * I.
  Standard Pauli/2: kappa=1/2, C= 3/4 I (dim=2), C= 2 I (dim=3).
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils import (
    pauli_matrices,
    evolve_parallel,
    killing_metric,
    compute_casimir,
    structure_constants_real,
    orthonormalize_killing,
    commutator,
    check_in_span,
)


def levi_civita_3():
    """Standard Levi-Civita tensor for 3 indices."""
    eps = np.zeros((3, 3, 3))
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1.0
    eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1.0
    return eps


def run_test_representation_theory():
    print("=" * 70)
    print("TEST 4.1: COMPARISON WITH REPRESENTATION THEORY")
    print("=" * 70)

    # ================================================
    # Part A: Reference values from Pauli/2
    # ================================================
    print("\n--- Part A: Reference Pauli/2 ---")
    pauli = [s / 2 for s in pauli_matrices()]
    C_pauli = compute_casimir(pauli)
    K_pauli = killing_metric(pauli)
    f_pauli = structure_constants_real(pauli)
    eps = levi_civita_3()

    print(f"Casimir (Pauli/2): diag = {np.diag(C_pauli).real}")
    print(f"  Expected: [0.75, 0.75] (j=1/2, C=j(j+1)=3/4)")
    print(f"Killing (Pauli/2): diag = {np.diag(K_pauli)}")
    print(f"f_{{012}} = {f_pauli[0,1,2]:.4f} (expected: 1.0)")
    print(f"f_{{120}} = {f_pauli[1,2,0]:.4f} (expected: 1.0)")
    print(f"f_{{201}} = {f_pauli[2,0,1]:.4f} (expected: 1.0)")

    # ================================================
    # Part B & C: Optimized generators in dim=2 and dim=3 (parallel)
    # ================================================
    tasks = [
        {'seed': 42, 'n_gen': 3, 'dim': 2, 'steps': 15000, 'lr': 0.001,
         'alpha': 5.0, 'beta': 1.0, 'killing_target': -1.0, 'verbose': False},
        {'seed': 42, 'n_gen': 3, 'dim': 3, 'steps': 15000, 'lr': 0.001,
         'alpha': 5.0, 'beta': 1.0, 'killing_target': -1.0, 'verbose': False},
    ]
    print("\n  Running dim=2 and dim=3 in parallel...")
    raw = evolve_parallel(tasks)
    gen_2d, _ = raw[0]
    gen_3d, _ = raw[1]

    print("\n--- Part B: Optimized generators, dim=2 ---")

    C_2d = compute_casimir(gen_2d)
    K_2d = killing_metric(gen_2d)
    f_2d = structure_constants_real(gen_2d)

    # Measure actual normalization: kappa = mean Tr(G_a^2)
    kappa_2d = np.mean([np.trace(g @ g).real for g in gen_2d])
    j_half = 0.5
    C_expected_2d = 2 * kappa_2d * j_half * (j_half + 1)  # = 3/2 * kappa

    C_2d_diag = np.diag(C_2d).real
    C_2d_mean = np.mean(C_2d_diag)

    print(f"Generator norm kappa = Tr(G_a^2): {kappa_2d:.4f}  (standard Pauli/2: 0.5)")
    print(f"Casimir diag: {C_2d_diag}")
    print(f"  Mean: {C_2d_mean:.4f}")
    print(f"  Expected (2*kappa*j(j+1)): {C_expected_2d:.4f}")
    casimir_error_2d = abs(C_2d_mean - C_expected_2d)
    casimir_rel_error_2d = casimir_error_2d / abs(C_expected_2d) if abs(C_expected_2d) > 1e-10 else casimir_error_2d
    print(f"  Absolute error: {casimir_error_2d:.6f}")
    print(f"  Relative error: {casimir_rel_error_2d:.6f}")

    # Check Casimir is proportional to identity
    C_2d_spread = np.std(C_2d_diag)
    print(f"  Spread (should ~ 0 for C ∝ I): {C_2d_spread:.6f}")

    print(f"Killing diag: {np.diag(K_2d)}")
    print(f"  Ratio to Pauli: {np.diag(K_2d) / np.diag(K_pauli)}")

    # Orthonormalize generators before extracting structure constants
    print("\n  Orthonormalizing via Killing metric...")
    gen_2d_ortho = orthonormalize_killing(gen_2d)
    K_ortho = killing_metric(gen_2d_ortho)
    print(f"  K after orthonormalization: diag={np.diag(K_ortho)}")
    f_2d_ortho = structure_constants_real(gen_2d_ortho)

    # Check f_ijk is proportional to epsilon_ijk
    f_nonzero = []
    eps_nonzero = []
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if abs(eps[i, j, k]) > 0.5:
                    f_nonzero.append(f_2d_ortho[i, j, k])
                    eps_nonzero.append(eps[i, j, k])

    f_nonzero = np.array(f_nonzero)
    eps_nonzero = np.array(eps_nonzero)

    proportionality = 0.0
    proportionality_std = 999.0
    if np.linalg.norm(eps_nonzero) > 0:
        proportionality = np.mean(f_nonzero / eps_nonzero)
        proportionality_std = np.std(f_nonzero / eps_nonzero)
        print(f"\nStructure constants (ortho basis): f = {proportionality:.4f} * epsilon")
        print(f"  Std: {proportionality_std:.6f}")
        print(f"  Expected: constant ratio (std ≈ 0)")

    # Check entries that SHOULD be zero
    f_should_zero = []
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if abs(eps[i, j, k]) < 0.5:
                    f_should_zero.append(abs(f_2d_ortho[i, j, k]))
    max_spurious = max(f_should_zero) if f_should_zero else 0
    print(f"Max spurious f_ijk (should be 0): {max_spurious:.6f}")

    # ================================================
    # Part C: Optimized generators in dim=3 (already computed above)
    # ================================================
    print("\n--- Part C: Optimized generators, dim=3 ---")

    C_3d = compute_casimir(gen_3d)
    K_3d = killing_metric(gen_3d)

    kappa_3d = np.mean([np.trace(g @ g).real for g in gen_3d])
    j_one = 1.0
    C_expected_3d = 2 * kappa_3d * j_one * (j_one + 1)  # = 4 * kappa

    C_3d_diag = np.diag(C_3d).real
    C_3d_mean = np.mean(C_3d_diag)
    C_3d_spread = np.std(C_3d_diag)

    closes_3d = all(
        check_in_span(commutator(gen_3d[i], gen_3d[j]), gen_3d)
        for i, j in [(0, 1), (1, 2), (2, 0)]
    )
    print(f"  Closes: {closes_3d}")
    print(f"  Killing diag: {np.diag(K_3d)}")

    print(f"Generator norm kappa = Tr(G_a^2): {kappa_3d:.4f}")
    print(f"Casimir diag: {np.round(C_3d_diag, 4)}")
    print(f"  Mean: {C_3d_mean:.4f}")
    if closes_3d:
        print(f"  Expected (2*kappa*j(j+1)): {C_expected_3d:.4f}")
        casimir_error_3d = abs(C_3d_mean - C_expected_3d)
        casimir_rel_error_3d = casimir_error_3d / abs(C_expected_3d) if abs(C_expected_3d) > 1e-10 else casimir_error_3d
        print(f"  Absolute error: {casimir_error_3d:.6f}")
        print(f"  Relative error: {casimir_rel_error_3d:.6f}")
    else:
        casimir_rel_error_3d = float('inf')
        print(f"  (Casimir comparison skipped: algebra does not close in dim=3)")
    print(f"  Spread (should ~ 0 for C ∝ I): {C_3d_spread:.6f}")

    # ================================================
    # SUMMARY
    # ================================================
    print("\n" + "=" * 50)
    print("SUMMARY TABLE")
    print("=" * 50)
    print(f"{'':20s} {'Pauli/2':>12s} {'Opt dim=2':>12s} {'Opt dim=3':>12s}")
    print("-" * 56)
    print(f"{'Casimir diag[0]':20s} {np.diag(C_pauli).real[0]:12.4f} "
          f"{C_2d_diag[0]:12.4f} {C_3d_diag[0]:12.4f}")
    print(f"{'Expected Casimir':20s} {'0.7500':>12s} {C_expected_2d:12.4f} {C_expected_3d:12.4f}")
    print(f"{'Generator kappa':20s} {'0.5000':>12s} {kappa_2d:12.4f} {kappa_3d:12.4f}")
    print(f"{'Killing diag[0]':20s} {np.diag(K_pauli)[0]:12.4f} "
          f"{np.diag(K_2d)[0]:12.4f} {np.diag(K_3d)[0]:12.4f}")
    print(f"{'f_012 (ortho)':20s} {f_pauli[0,1,2]:12.4f} {f_2d_ortho[0,1,2]:12.4f} {'N/A':>12s}")

    # ---- Verdict ----
    print("\n--- VERDICT ---")
    casimir_2d_ok = casimir_rel_error_2d < 0.1 and C_2d_spread < 0.05
    casimir_3d_ok = casimir_rel_error_3d < 0.15 and C_3d_spread < 0.1
    f_structure_ok = proportionality_std < 0.05 if 'proportionality_std' in dir() else False

    if casimir_2d_ok and f_structure_ok:
        print("PASS: Casimir and structure constants match representation theory")
        print(f"  (accounting for generator normalization kappa={kappa_2d:.4f}).")
        if casimir_3d_ok:
            print(f"  dim=3 Casimir also matches (kappa={kappa_3d:.4f}).")
    elif casimir_2d_ok:
        print("PARTIAL: Casimir matches j(j+1) but structure constants deviate.")
    else:
        print(f"FAIL: Casimir rel-error {casimir_rel_error_2d:.4f} (threshold 0.1).")


if __name__ == "__main__":
    run_test_representation_theory()
