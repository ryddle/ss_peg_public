"""
Test 3.3: Energy landscape scan.
Question: Is SU(2) at the global minimum of the energy landscape?

Samples many random operator sets and plots T_nc vs T_k, marking
which ones have SU(2) structure.

EXPECTED: SU(2) configurations cluster at the global minimum.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils import (
    random_generators,
    killing_metric,
    non_commutativity_tension,
    commutator,
    check_in_span,
    pauli_matrices,
)


def check_is_su2_structure(generators, tol=0.1):
    """Heuristic: closure + K proportional to I."""
    if len(generators) != 3:
        return False

    # Check closure
    for i, j in [(0, 1), (1, 2), (2, 0)]:
        comm = commutator(generators[i], generators[j])
        if not check_in_span(comm, generators, tol=tol):
            return False

    # Check K proportional to identity
    K = killing_metric(generators)
    K_offdiag = np.max(np.abs(K - np.diag(np.diag(K))))
    if K_offdiag > tol:
        return False

    return True


def run_test_landscape(n_samples=2000):
    print("=" * 70)
    print("TEST 3.3: ENERGY LANDSCAPE SCAN")
    print("=" * 70)

    results = []
    np.random.seed(0)

    for i in range(n_samples):
        gens = random_generators(n=3, dim=2)

        # Normalize to unit norm
        for j in range(len(gens)):
            norm = np.linalg.norm(gens[j])
            if norm > 1e-10:
                gens[j] = gens[j] / norm

        T_nc = non_commutativity_tension(gens)
        K = killing_metric(gens)
        K_target = -1.0 * np.eye(3)
        T_k = np.linalg.norm(K - K_target) ** 2
        T_total = T_nc + 5.0 * T_k

        is_su2 = check_is_su2_structure(gens, tol=0.15)

        results.append({
            "T_nc": T_nc,
            "T_k": T_k,
            "T_total": T_total,
            "is_su2": is_su2,
            "K_diag": np.mean(np.diag(K)),
        })

        if (i + 1) % 500 == 0:
            print(f"  Sampled {i + 1}/{n_samples}...")

    # ---- Include reference points ----
    # Pauli/2 (perfect SU(2))
    pauli_gens = [s / 2 for s in pauli_matrices()]
    T_nc_pauli = non_commutativity_tension(pauli_gens)
    K_pauli = killing_metric(pauli_gens)
    T_k_pauli = np.linalg.norm(K_pauli - (-1.0 * np.eye(3))) ** 2
    T_total_pauli = T_nc_pauli + 5.0 * T_k_pauli

    # ---- Analysis ----
    su2_points = [r for r in results if r["is_su2"]]
    non_su2_points = [r for r in results if not r["is_su2"]]

    print(f"\n--- RESULTS ---")
    print(f"Total samples: {n_samples}")
    print(f"SU(2)-like: {len(su2_points)}")
    print(f"T_total range: [{min(r['T_total'] for r in results):.4f}, "
          f"{max(r['T_total'] for r in results):.4f}]")
    print(f"Pauli reference: T_nc={T_nc_pauli:.4f}, T_k={T_k_pauli:.4f}, "
          f"T_total={T_total_pauli:.4f}")

    if su2_points:
        print(f"SU(2) T_total range: [{min(r['T_total'] for r in su2_points):.4f}, "
              f"{max(r['T_total'] for r in su2_points):.4f}]")

    # ---- Plot ----
    try:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Plot 1: T_nc vs T_k
        ax = axes[0]
        ax.scatter(
            [r["T_nc"] for r in non_su2_points],
            [r["T_k"] for r in non_su2_points],
            alpha=0.2, s=10, label="Non-SU(2)", c="gray",
        )
        if su2_points:
            ax.scatter(
                [r["T_nc"] for r in su2_points],
                [r["T_k"] for r in su2_points],
                alpha=0.8, s=50, label="SU(2)", c="red",
            )
        ax.scatter([T_nc_pauli], [T_k_pauli], s=200, c="gold", marker="*",
                   label="Pauli/2", zorder=10, edgecolors="black")
        ax.set_xlabel("T_nc (non-commutativity)")
        ax.set_ylabel("T_k (Killing deviation)")
        ax.legend()
        ax.set_title("Energy Landscape: T_nc vs T_k")
        ax.grid(True, alpha=0.3)

        # Plot 2: Histogram of T_total
        ax = axes[1]
        ax.hist([r["T_total"] for r in results], bins=50, alpha=0.7, label="All samples")
        ax.axvline(T_total_pauli, color="red", linestyle="--", linewidth=2, label="Pauli/2")
        ax.set_xlabel("T_total")
        ax.set_ylabel("Count")
        ax.set_title("Distribution of Total Tension")
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(os.path.dirname(__file__), "test_3_3_energy_landscape.png"), dpi=150)
        print("\nPlot saved: test_3_3_energy_landscape.png")
    except ImportError:
        print("(matplotlib not available, skipping plot)")

    # ---- Verdict ----
    print("\n--- VERDICT ---")
    # Check if Pauli is at/near the minimum
    min_T_total = min(r["T_total"] for r in results)
    if T_total_pauli <= min_T_total * 1.1:
        print("PASS: Pauli/SU(2) is at the global minimum of the energy landscape.")
    else:
        print(f"FAIL: Pauli T_total={T_total_pauli:.4f} but sampled min={min_T_total:.4f}.")
        print("  SU(2) is NOT at the global minimum.")

    return results


if __name__ == "__main__":
    run_test_landscape(n_samples=2000)
