"""
Ramanujan Check for SS-PEG Emergent Algebras
==============================================

For each Lie algebra that emerges from the SS-PEG variational principle,
we construct the "adjoint graph" from the structure constants and compute:

  1. Adjacency matrix A from the structure constants f_{abc}
  2. Eigenvalues of A (spectrum)
  3. Ihara zeta function Z(u) and its poles
  4. Ramanujan property: are all non-trivial eigenvalues ≤ 2√(q)?
  5. Spectral gap analysis
  6. Connection to Killing eigenvalue classification (cyclic/free/decoupled)

The adjoint graph of a Lie algebra:
  - Vertices = generators {T_a}
  - Edge (a,b) exists if [T_a, T_b] ≠ 0 (they interact via commutator)
  - Edge weight = Σ_c |f_{abc}|² (strength of interaction)

The Ramanujan property would mean the algebra's interaction graph has
OPTIMAL spectral gap — maximum expansion/mixing with minimal number of edges.

Key prediction: the density transition at ρ=1 (from SYNTHESIS.md) should
correspond to the graph becoming Ramanujan or near-Ramanujan.

Usage:
    python ramanujan_check.py [--save]
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import argparse
import sys
import os

# Add verification/ to path for shared_utils
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'verification'))
from shared_utils import (
    pauli_matrices, gell_mann_matrices, commutator,
    killing_metric, structure_constants
)


# ============================================================
# Lie algebra generators catalog
# ============================================================

def su2_generators():
    """su(2): 3 generators (Pauli/2)."""
    paulis = pauli_matrices()
    return [p / 2 for p in paulis]


def su3_generators():
    """su(3): 8 Gell-Mann matrices (normalized)."""
    gm = gell_mann_matrices()
    return [g / 2 for g in gm]


def so3_generators():
    """so(3) ≅ su(2): 3 generators."""
    J1 = np.array([[0,0,0],[0,0,-1],[0,1,0]], dtype=np.complex128)
    J2 = np.array([[0,0,1],[0,0,0],[-1,0,0]], dtype=np.complex128)
    J3 = np.array([[0,-1,0],[1,0,0],[0,0,0]], dtype=np.complex128)
    return [J1, J2, J3]


def so4_generators():
    """so(4): 6 generators (antisymmetric 4×4)."""
    gens = []
    dim = 4
    for i in range(dim):
        for j in range(i+1, dim):
            E = np.zeros((dim, dim), dtype=np.complex128)
            E[i, j] = 1
            E[j, i] = -1
            gens.append(E)
    return gens


def so31_generators():
    """so(3,1) Lorentz: 6 generators (3 rotations + 3 boosts)."""
    # Rotation generators (compact)
    J1 = np.zeros((4,4), dtype=np.complex128); J1[1,2] = -1; J1[2,1] = 1
    J2 = np.zeros((4,4), dtype=np.complex128); J2[0,2] = 1; J2[2,0] = -1
    J3 = np.zeros((4,4), dtype=np.complex128); J3[0,1] = -1; J3[1,0] = 1

    # Boost generators (non-compact): K_i with η = diag(-1,1,1,1)
    K1 = np.zeros((4,4), dtype=np.complex128); K1[0,1] = 1; K1[1,0] = 1
    K2 = np.zeros((4,4), dtype=np.complex128); K2[0,2] = 1; K2[2,0] = 1
    K3 = np.zeros((4,4), dtype=np.complex128); K3[0,3] = 1; K3[3,0] = 1

    return [J1, J2, J3, K1, K2, K3]


def sp4_generators():
    """sp(4): 10 generators. X such that X^T Ω + Ω X = 0, Ω = [[0,I],[-I,0]]."""
    gens = []
    dim = 4
    Omega = np.zeros((dim, dim), dtype=np.complex128)
    half = dim // 2
    Omega[:half, half:] = np.eye(half)
    Omega[half:, :half] = -np.eye(half)

    # Basis: matrices X satisfying X^T Ω + Ω X = 0
    for i in range(dim):
        for j in range(dim):
            E = np.zeros((dim, dim), dtype=np.complex128)
            E[i, j] = 1
            # Symplectic condition: X^T Ω + Ω X = 0
            X = E
            cond = X.T @ Omega + Omega @ X
            if np.linalg.norm(cond) < 1e-10:
                gens.append(X)

    # If brute force doesn't give clean basis, construct explicitly
    if len(gens) != 10:
        gens = []
        # sp(4) has basis: symmetric + antisymmetric blocks
        # Block form: X = [[A, B], [C, -A^T]] with B=B^T, C=C^T
        for i in range(2):
            for j in range(2):
                # A block
                E = np.zeros((4,4), dtype=np.complex128)
                E[i, j] = 1
                E[j+2, i+2] = -1
                if np.linalg.norm(E) > 0.5:
                    gens.append(E)

        # B block (symmetric)
        for i in range(2):
            for j in range(i, 2):
                E = np.zeros((4,4), dtype=np.complex128)
                E[i, j+2] = 1
                if i != j:
                    E[j, i+2] = 1
                else:
                    pass  # diagonal
                gens.append(E)

        # C block (symmetric)
        for i in range(2):
            for j in range(i, 2):
                E = np.zeros((4,4), dtype=np.complex128)
                E[i+2, j] = 1
                if i != j:
                    E[j+2, i] = 1
                gens.append(E)

    # Deduplicate and ensure we have 10
    unique_gens = []
    for g in gens:
        is_dup = False
        for ug in unique_gens:
            if np.linalg.norm(g - ug) < 1e-10:
                is_dup = True
                break
        if not is_dup and np.linalg.norm(g) > 1e-10:
            unique_gens.append(g)

    return unique_gens[:10]


def u2_generators():
    """u(2) = su(2) ⊕ u(1): 4 generators."""
    paulis = pauli_matrices()
    gens = [p / 2 for p in paulis]
    gens.append(np.eye(2, dtype=np.complex128) / 2)  # U(1) part
    return gens


# ============================================================
# Adjoint graph construction
# ============================================================

def build_adjoint_graph(generators, threshold=1e-8):
    """
    Build the adjoint interaction graph of a Lie algebra.

    Vertices = generators
    Edge (a,b) iff [T_a, T_b] ≠ 0
    Weight(a,b) = ||[T_a, T_b]||²

    Returns:
      adjacency: (n, n) binary adjacency matrix
      weighted: (n, n) weighted adjacency matrix
      degree: vertex degrees
    """
    n = len(generators)
    adjacency = np.zeros((n, n))
    weighted = np.zeros((n, n))

    for i in range(n):
        for j in range(i+1, n):
            C = commutator(generators[i], generators[j])
            norm_sq = np.linalg.norm(C)**2
            if norm_sq > threshold:
                adjacency[i, j] = 1
                adjacency[j, i] = 1
                weighted[i, j] = norm_sq
                weighted[j, i] = norm_sq

    degree = adjacency.sum(axis=1).astype(int)
    return adjacency, weighted, degree


# ============================================================
# Spectral analysis
# ============================================================

def spectral_analysis(adjacency):
    """
    Compute eigenvalues of adjacency matrix and check Ramanujan property.

    For a (q+1)-regular graph, Ramanujan ⟺ all |λ| ≤ 2√q for
    non-trivial eigenvalues (those ≠ ±(q+1)).

    For non-regular graphs, we use the general bound and report
    the spectral gap.
    """
    n = adjacency.shape[0]
    eigenvalues = np.linalg.eigvalsh(adjacency)  # symmetric → real eigenvalues
    eigenvalues = np.sort(eigenvalues)[::-1]  # descending

    # Check regularity
    degrees = adjacency.sum(axis=1)
    is_regular = np.allclose(degrees, degrees[0])
    q_plus_1 = degrees[0] if is_regular else np.mean(degrees)
    q = q_plus_1 - 1

    # Trivial eigenvalue is the largest: should be q+1 for regular
    lambda_max = eigenvalues[0]
    lambda_2 = eigenvalues[1] if n > 1 else 0

    # Ramanujan bound
    ramanujan_bound = 2 * np.sqrt(max(q, 0))

    # Non-trivial eigenvalues (exclude largest and smallest if bipartite)
    nontrivial = eigenvalues[1:-1] if n > 2 else eigenvalues[1:]
    max_nontrivial = np.max(np.abs(nontrivial)) if len(nontrivial) > 0 else 0

    is_ramanujan = max_nontrivial <= ramanujan_bound + 1e-10

    # Spectral gap
    spectral_gap = lambda_max - lambda_2

    return {
        'eigenvalues': eigenvalues,
        'is_regular': is_regular,
        'degree': q_plus_1,
        'q': q,
        'lambda_max': lambda_max,
        'lambda_2': lambda_2,
        'max_nontrivial': max_nontrivial,
        'ramanujan_bound': ramanujan_bound,
        'is_ramanujan': is_ramanujan,
        'spectral_gap': spectral_gap,
        'ratio': max_nontrivial / ramanujan_bound if ramanujan_bound > 0 else float('inf')
    }


# ============================================================
# Ihara zeta for the adjoint graph
# ============================================================

def ihara_zeta_poles(adjacency):
    """
    For a graph with adjacency matrix A, the Ihara zeta poles come from:
      det(I - uA + u²Q) = 0
    where Q = D - I (D = degree matrix).

    For regular graphs: Q = qI, so poles at u = 1/λ where λ satisfies
    1 - uλ + qu² = 0  →  u = (λ ± √(λ²-4q)) / (2q)
    """
    n = adjacency.shape[0]
    D = np.diag(adjacency.sum(axis=1))
    Q = D - np.eye(n)
    I = np.eye(n)

    # Find u where det(I - uA + u²Q) = 0
    # This is a polynomial eigenvalue problem. We convert to a companion form.
    # det(I - uA + u²Q) = det(Q) · det(Q^{-1} - uQ^{-1}A + u²I)
    # = det(Q) · det(u²I - u(Q^{-1}A) + Q^{-1})

    # Alternative: linearize as a 2n×2n eigenvalue problem
    # [I  -A] [ψ ]     [0  -Q] [ψ ]
    # [0   I] [uψ]  = u[I   0] [uψ]
    # → generalized eigenvalue: (M - u N)v = 0

    Z = np.zeros((n, n))
    M = np.block([[I, -adjacency], [Z, I]])
    N = np.block([[Z, -Q], [I, Z]])

    # Solve generalized eigenvalue problem
    try:
        from numpy.linalg import solve
        # M v = u N v → N^{-1} M v = u v, if N is invertible
        if abs(np.linalg.det(N)) > 1e-15:
            eigs = np.linalg.eigvals(np.linalg.solve(N, M))
        else:
            # Fallback: direct from adjacency eigenvalues (regular case)
            degrees = adjacency.sum(axis=1)
            q = degrees[0] - 1 if np.allclose(degrees, degrees[0]) else np.mean(degrees) - 1
            eigs_A = np.linalg.eigvals(adjacency)
            eigs = []
            for lam in eigs_A:
                disc = lam**2 - 4*q + 0j
                sqrt_d = np.sqrt(disc)
                if abs(q) > 1e-10:
                    eigs.extend([(lam + sqrt_d)/(2*q), (lam - sqrt_d)/(2*q)])
                else:
                    if abs(lam) > 1e-10:
                        eigs.append(1.0/lam)
            eigs = np.array(eigs)
    except Exception:
        eigs = np.array([])

    return eigs


# ============================================================
# Killing classification overlay
# ============================================================

def killing_edge_classification(generators):
    """
    Classify edges according to Killing eigenvalues:
      K_λ < 0 → cyclic (compact)
      K_λ > 0 → free (non-compact)
      K_λ = 0 → decoupled (abelian)

    Returns counts and the Killing eigenvalue spectrum.
    """
    K = killing_metric(generators)
    eigenvalues = np.linalg.eigvalsh(K)

    n_cyclic = np.sum(eigenvalues < -1e-10)
    n_free = np.sum(eigenvalues > 1e-10)
    n_decoupled = np.sum(np.abs(eigenvalues) <= 1e-10)

    return {
        'killing_eigenvalues': eigenvalues,
        'n_cyclic': int(n_cyclic),
        'n_free': int(n_free),
        'n_decoupled': int(n_decoupled),
        'signature': f"({n_free}+, {n_decoupled}₀, {n_cyclic}−)"
    }


# ============================================================
# Full analysis for all algebras
# ============================================================

ALGEBRA_CATALOG = {
    'su(2)': su2_generators,
    'su(3)': su3_generators,
    'so(3)': so3_generators,
    'so(4)': so4_generators,
    'so(3,1)': so31_generators,
    'u(2)': u2_generators,
}


def analyze_all():
    """Run full analysis on all algebras in catalog."""
    results = {}

    for name, gen_fn in ALGEBRA_CATALOG.items():
        print(f"\n{'='*60}")
        print(f"  {name}")
        print(f"{'='*60}")

        generators = gen_fn()
        n_gen = len(generators)
        dim = generators[0].shape[0]

        print(f"  dim(algebra) = {n_gen}, representation dim = {dim}")

        # Adjoint graph
        adj, weighted, degrees = build_adjoint_graph(generators)
        print(f"  Adjacency: {int(adj.sum()//2)} edges, degrees = {degrees}")

        # Spectral analysis
        spec = spectral_analysis(adj)
        print(f"  Regular: {spec['is_regular']}, degree = {spec['degree']:.0f}")
        print(f"  Eigenvalues: {spec['eigenvalues']}")
        print(f"  λ_max = {spec['lambda_max']:.4f}, λ_2 = {spec['lambda_2']:.4f}")
        print(f"  Spectral gap = {spec['spectral_gap']:.4f}")
        print(f"  Ramanujan bound 2√q = {spec['ramanujan_bound']:.4f}")
        print(f"  max|λ_nontrivial| = {spec['max_nontrivial']:.4f}")
        print(f"  Ramanujan ratio = {spec['ratio']:.4f}")
        print(f"  IS RAMANUJAN: {'YES ✓' if spec['is_ramanujan'] else 'NO ✗'}")

        # Killing classification
        killing = killing_edge_classification(generators)
        print(f"  Killing signature: {killing['signature']}")
        print(f"  Killing eigenvalues: {killing['killing_eigenvalues']}")

        # Ihara zeta poles
        poles = ihara_zeta_poles(adj)
        if len(poles) > 0:
            print(f"  Ihara poles (|u| sorted): {np.sort(np.abs(poles))[:6]}...")
        else:
            print(f"  Ihara poles: could not compute")

        results[name] = {
            'generators': generators,
            'adjacency': adj,
            'weighted': weighted,
            'degrees': degrees,
            'spectral': spec,
            'killing': killing,
            'ihara_poles': poles
        }

    return results


# ============================================================
# Visualization
# ============================================================

def plot_all(results, save=False):
    algebras = list(results.keys())
    n_alg = len(algebras)

    fig = plt.figure(figsize=(20, 16))
    fig.suptitle(
        'SS-PEG: Ramanujan Analysis of Emergent Algebra Adjoint Graphs',
        fontsize=16, fontweight='bold'
    )
    gs = GridSpec(3, n_alg, figure=fig, hspace=0.4, wspace=0.35)

    # Row 1: Adjacency eigenvalue spectra
    for i, name in enumerate(algebras):
        ax = fig.add_subplot(gs[0, i])
        spec = results[name]['spectral']
        eigs = spec['eigenvalues']
        colors = ['red' if abs(e - spec['lambda_max']) < 0.01 else
                  'orange' if abs(e) > spec['ramanujan_bound'] else
                  'steelblue' for e in eigs]
        ax.bar(range(len(eigs)), eigs, color=colors, alpha=0.7, edgecolor='black',
              linewidth=0.5)
        ax.axhline(spec['ramanujan_bound'], color='green', linestyle='--',
                  alpha=0.7, label=f'2√q={spec["ramanujan_bound"]:.2f}')
        ax.axhline(-spec['ramanujan_bound'], color='green', linestyle='--', alpha=0.7)
        ram_str = '✓' if spec['is_ramanujan'] else '✗'
        ax.set_title(f'{name}\nRamanujan: {ram_str} (r={spec["ratio"]:.2f})',
                    fontsize=10)
        ax.set_xlabel('Index')
        if i == 0:
            ax.set_ylabel('Eigenvalue')
        ax.legend(fontsize=6)
        ax.grid(True, alpha=0.3)

    # Row 2: Ihara zeta poles in complex plane
    for i, name in enumerate(algebras):
        ax = fig.add_subplot(gs[1, i])
        poles = results[name]['ihara_poles']
        if len(poles) > 0:
            ax.scatter(poles.real, poles.imag, s=40, c='red', marker='x',
                      linewidths=1.5, alpha=0.7)

        theta = np.linspace(0, 2*np.pi, 100)
        ax.plot(np.cos(theta), np.sin(theta), 'k--', alpha=0.2, label='|u|=1')

        spec = results[name]['spectral']
        q = max(spec['q'], 0.01)
        r_ram = 1.0 / np.sqrt(q) if q > 0 else 2.0
        if r_ram < 3:
            ax.plot(r_ram*np.cos(theta), r_ram*np.sin(theta), 'g--',
                   alpha=0.4, label=f'|u|=1/√q={r_ram:.2f}')

        ax.set_title(f'{name} — Ihara poles', fontsize=10)
        ax.set_aspect('equal')
        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(-2.5, 2.5)
        ax.legend(fontsize=6)
        ax.grid(True, alpha=0.3)

    # Row 3: Killing eigenvalue spectrum + edge classification
    for i, name in enumerate(algebras):
        ax = fig.add_subplot(gs[2, i])
        killing = results[name]['killing']
        k_eigs = killing['killing_eigenvalues']
        colors = ['green' if e < -1e-10 else 'red' if e > 1e-10 else 'grey'
                  for e in k_eigs]
        ax.bar(range(len(k_eigs)), k_eigs, color=colors, alpha=0.7,
              edgecolor='black', linewidth=0.5)
        ax.axhline(0, color='black', linewidth=0.5)
        ax.set_title(f'{name}\nKilling: {killing["signature"]}', fontsize=10)
        ax.set_xlabel('Generator index')
        if i == 0:
            ax.set_ylabel('K_λ')
        ax.grid(True, alpha=0.3)

        # Legend for edge types
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='green', alpha=0.7, label=f'Cyclic ({killing["n_cyclic"]})'),
            Patch(facecolor='red', alpha=0.7, label=f'Free ({killing["n_free"]})'),
            Patch(facecolor='grey', alpha=0.7, label=f'Decoupled ({killing["n_decoupled"]})'),
        ]
        ax.legend(handles=legend_elements, fontsize=6, loc='upper right')

    plt.tight_layout()

    if save:
        out_path = 'ramanujan_check_results.png'
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        print(f"\nSaved to {out_path}")

    plt.show()


def print_summary_table(results):
    """Print a compact summary table."""
    print("\n" + "=" * 100)
    print("SUMMARY: Ramanujan Property of SS-PEG Adjoint Graphs")
    print("=" * 100)
    print(f"{'Algebra':>10} | {'dim':>4} | {'Edges':>5} | {'Regular':>7} | "
          f"{'λ_max':>6} | {'λ_2':>6} | {'2√q':>6} | {'max|λ_nt|':>9} | "
          f"{'Ratio':>6} | {'Ramanujan':>9} | {'Killing sig':>15}")
    print("-" * 100)

    for name, res in results.items():
        spec = res['spectral']
        killing = res['killing']
        n_gen = len(res['generators'])
        n_edges = int(res['adjacency'].sum() // 2)
        ram_str = 'YES ✓' if spec['is_ramanujan'] else 'NO ✗'

        print(f"{name:>10} | {n_gen:4d} | {n_edges:5d} | "
              f"{'yes' if spec['is_regular'] else 'no':>7} | "
              f"{spec['lambda_max']:6.2f} | {spec['lambda_2']:6.2f} | "
              f"{spec['ramanujan_bound']:6.2f} | {spec['max_nontrivial']:9.4f} | "
              f"{spec['ratio']:6.3f} | {ram_str:>9} | {killing['signature']:>15}")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Ramanujan check for SS-PEG emergent algebra adjoint graphs'
    )
    parser.add_argument('--save', action='store_true', help='Save plots to file')
    args = parser.parse_args()

    results = analyze_all()
    print_summary_table(results)
    plot_all(results, save=args.save)
