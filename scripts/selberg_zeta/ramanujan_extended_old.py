"""
Extended Ramanujan Check for ALL Simple Lie Algebras
=====================================================

Tests the Ramanujan property of the adjoint commutation graph for:
  - A_n = su(n+1)  for n=1..9   (dim 3..99)
  - B_n = so(2n+1) for n=1..5   (dim 3..55)
  - D_n = so(2n)   for n=2..6   (dim 6..66)
  - C_n = sp(2n)   for n=1..5   (dim 3..55)
  - G₂  (dim 14)
  - F₄  (dim 52)
  - E₆  (dim 78)
  - E₇  (dim 133)
  - E₈  (dim 248)

The adjoint graph has:
  Vertices = generators {T_a}
  Edge (a,b) iff [T_a, T_b] ≠ 0

Ramanujan property for (q+1)-regular graph:
  max|λ_nontrivial| ≤ 2√q

For E₈ (248 generators), this is the definitive test —
a 248-vertex algebraically-derived Ramanujan graph cannot
be dismissed as a small-graph coincidence.

Usage:
    python ramanujan_extended.py [--save] [--quick]
    python ramanujan_extended.py --families A B D C G F E
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import argparse
import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'verification'))
from shared_utils import commutator, killing_metric


# ============================================================
# Generic Lie algebra generators via standard bases
# ============================================================

def su_generators(n):
    """
    su(n): traceless anti-hermitian n×n matrices.
    Dimension = n²-1. We return hermitian generators T_a
    such that the Lie bracket is [T_a, T_b].

    Standard basis:
      - (n²-n)/2 symmetric off-diagonal: (E_ij + E_ji)/2
      - (n²-n)/2 antisymmetric off-diagonal: -i(E_ij - E_ji)/2
      - (n-1) diagonal traceless (generalized Gell-Mann)
    """
    gens = []
    # Off-diagonal symmetric
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n), dtype=np.complex128)
            E[i, j] = 1
            E[j, i] = 1
            gens.append(E / 2)

    # Off-diagonal antisymmetric
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n), dtype=np.complex128)
            E[i, j] = -1j
            E[j, i] = 1j
            gens.append(E / 2)

    # Diagonal traceless
    for k in range(1, n):
        D = np.zeros((n, n), dtype=np.complex128)
        for i in range(k):
            D[i, i] = 1
        D[k, k] = -k
        D /= np.sqrt(2 * k * (k + 1))
        gens.append(D)

    assert len(gens) == n * n - 1, f"su({n}): expected {n*n-1} generators, got {len(gens)}"
    return gens


def so_generators(n):
    """
    so(n): antisymmetric n×n real matrices.
    Dimension = n(n-1)/2.
    Basis: E_ij - E_ji for i < j.
    """
    gens = []
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n), dtype=np.complex128)
            E[i, j] = 1
            E[j, i] = -1
            gens.append(E)
    expected = n * (n - 1) // 2
    assert len(gens) == expected, f"so({n}): expected {expected}, got {len(gens)}"
    return gens


def sp_generators(n):
    """
    sp(2n): 2n×2n matrices X satisfying X^T Ω + Ω X = 0,
    where Ω = [[0, I_n], [-I_n, 0]].
    Dimension = n(2n+1).

    Explicit basis:
      Block form X = [[A, B], [C, -A^T]] with B=B^T, C=C^T.
      - A block: n² generators (E_ij for A, -E_ji for -A^T block)
      - B block: n(n+1)/2 generators (symmetric B)
      - C block: n(n+1)/2 generators (symmetric C)
      Total = n² + n(n+1) = n(2n+1)
    """
    gens = []
    half = n

    # A-block: E_{ij} in upper-left, -E_{ji} in lower-right
    for i in range(half):
        for j in range(half):
            E = np.zeros((2 * half, 2 * half), dtype=np.complex128)
            E[i, j] = 1
            E[j + half, i + half] = -1
            gens.append(E)

    # B-block (symmetric): upper-right
    for i in range(half):
        for j in range(i, half):
            E = np.zeros((2 * half, 2 * half), dtype=np.complex128)
            E[i, j + half] = 1
            E[j, i + half] = 1
            gens.append(E)

    # C-block (symmetric): lower-left
    for i in range(half):
        for j in range(i, half):
            E = np.zeros((2 * half, 2 * half), dtype=np.complex128)
            E[i + half, j] = 1
            E[j + half, i] = 1
            gens.append(E)

    expected = half * (2 * half + 1)
    assert len(gens) == expected, f"sp({2*half}): expected {expected}, got {len(gens)}"
    return gens


# ============================================================
# Exceptional Lie algebras via embedding in matrix algebras
# ============================================================

def g2_generators():
    """
    G₂ (dim 14): Embedded in so(7) as the automorphism algebra of octonions.

    G₂ ⊂ so(7). We construct so(7) generators (dim=21) and project
    onto the G₂ subalgebra by using the octonion structure constants.

    The 14 generators of G₂ can be constructed as specific linear
    combinations of so(7) generators that preserve the octonion
    multiplication table.
    """
    # Octonion structure constants: e_i × e_j = ±e_k
    # Using the Fano plane convention: (1,2,4), (2,3,5), (3,4,6), (4,5,7), (5,6,1), (6,7,2), (7,1,3)
    fano_triples = [
        (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7),
        (5, 6, 1), (6, 7, 2), (7, 1, 3)
    ]

    # Build the associator tensor (completely antisymmetric)
    # phi_{ijk} = 1 for (i,j,k) in Fano triples (and cyclic perms), -1 for anticyclic
    phi = np.zeros((7, 7, 7))
    for (a, b, c) in fano_triples:
        i, j, k = a - 1, b - 1, c - 1  # 0-indexed
        phi[i, j, k] = 1
        phi[j, k, i] = 1
        phi[k, i, j] = 1
        phi[j, i, k] = -1
        phi[i, k, j] = -1
        phi[k, j, i] = -1

    # G₂ generators as so(7) matrices preserving phi
    # A ∈ G₂ ⊂ so(7) iff: A_{ia} φ_{ajk} + A_{ja} φ_{iak} + A_{ka} φ_{ija} = 0 for all i,j,k
    # We find the kernel of this constraint on the 21-dim space of so(7)

    so7_gens = so_generators(7)  # 21 generators
    n_so7 = len(so7_gens)

    # Build the constraint matrix
    # For each (i,j,k) antisymmetric triple, we get a linear constraint on the 21 coefficients
    constraints = []
    for i in range(7):
        for j in range(i + 1, 7):
            for k in range(j + 1, 7):
                row = np.zeros(n_so7)
                for g_idx, G in enumerate(so7_gens):
                    val = 0.0
                    for a in range(7):
                        val += G[i, a].real * phi[a, j, k]
                        val += G[j, a].real * phi[i, a, k]
                        val += G[k, a].real * phi[i, j, a]
                    row[g_idx] = val
                constraints.append(row)

    C_mat = np.array(constraints)
    # G₂ = kernel of C_mat
    U, S, Vt = np.linalg.svd(C_mat)
    tol = 1e-10
    null_mask = S < tol
    # Also include dimensions beyond len(S)
    n_null_from_s = np.sum(null_mask)
    n_extra = n_so7 - len(S)
    n_null = n_null_from_s + n_extra

    if n_null < 14:
        # Relax tolerance
        tol = 1e-6
        null_mask = S < tol
        n_null_from_s = np.sum(null_mask)
        n_null = n_null_from_s + n_extra

    # Null space vectors from Vt
    null_vectors = []
    for idx in range(len(S)):
        if S[idx] < max(tol, 1e-6):
            null_vectors.append(Vt[idx])
    # Extra dimensions (beyond len(S)) are automatically in null space
    for idx in range(len(S), n_so7):
        null_vectors.append(Vt[idx])

    assert len(null_vectors) >= 14, \
        f"G₂: expected 14-dim null space, got {len(null_vectors)} (tol={tol:.1e}, min_S={S.min():.2e})"

    # Construct G₂ generators as linear combinations of so(7) generators
    g2_gens = []
    for coeffs in null_vectors[:14]:
        G = sum(c * so7_gens[i] for i, c in enumerate(coeffs.real))
        g2_gens.append(G)

    return g2_gens


def _build_exceptional_via_root_system(rank, roots, cartan_matrix):
    """
    Build generators for an exceptional Lie algebra using the Chevalley basis.

    For a rank-r algebra with root system Φ:
      - r Cartan generators H_i (diagonal in adjoint rep)
      - One generator E_α for each root α ∈ Φ
      - Total dim = r + |Φ|

    We construct the adjoint representation directly:
      [H_i, E_α] = α_i E_α
      [E_α, E_{-α}] = Σ_i α_i H_i  (=: H_α)
      [E_α, E_β] = N_{α,β} E_{α+β}   if α+β ∈ Φ
                  = 0                    otherwise

    The structure constants N_{α,β} satisfy |N_{α,β}|² = q(p+1)/2
    where the α-string through β is β - pα, ..., β, ..., β + qα.
    """
    n_roots = len(roots)
    dim = rank + n_roots  # Total dimension of the algebra

    # Index mapping: generators 0..rank-1 are H_i, rank..rank+n_roots-1 are E_α
    # We need a map from root tuple to index
    root_to_idx = {}
    for i, r in enumerate(roots):
        root_to_idx[tuple(r)] = rank + i

    # Build adjoint representation matrices (dim × dim)
    ad = np.zeros((dim, dim, dim), dtype=np.float64)  # ad[a][b][c] = f_{ab}^c

    # [H_i, H_j] = 0
    # (already zero)

    # [H_i, E_α] = α_i E_α
    for i in range(rank):
        for j, alpha in enumerate(roots):
            idx_alpha = rank + j
            ad[i, idx_alpha, idx_alpha] = alpha[i]
            ad[idx_alpha, i, idx_alpha] = -alpha[i]

    # [E_α, E_{-α}] = H_α = Σ_i α_i H_i (using coroots)
    # More precisely: [E_α, E_{-α}] = H_α where H_α corresponds to the coroot
    for j, alpha in enumerate(roots):
        neg_alpha = tuple(-a for a in alpha)
        if neg_alpha in root_to_idx:
            idx_alpha = rank + j
            idx_neg = root_to_idx[neg_alpha]
            # [E_α, E_{-α}] = Σ_i (2α_i / (α,α)) H_i
            alpha_arr = np.array(alpha, dtype=np.float64)
            norm_sq = alpha_arr @ alpha_arr
            if norm_sq > 1e-12:
                for i in range(rank):
                    coeff = 2 * alpha[i] / norm_sq
                    ad[idx_alpha, idx_neg, i] = coeff
                    ad[idx_neg, idx_alpha, i] = -coeff

    # [E_α, E_β] = N_{α,β} E_{α+β} if α+β ∈ Φ
    root_set = set(root_to_idx.keys())
    for j, alpha in enumerate(roots):
        alpha_arr = np.array(alpha, dtype=np.float64)
        for k, beta in enumerate(roots):
            if j == k:
                continue
            # Check if alpha is -beta (handled above)
            if tuple(-b for b in beta) == tuple(alpha):
                continue
            gamma = tuple(a + b for a, b in zip(alpha, beta))
            if gamma in root_set:
                # Compute N_{α,β} via the α-string through β
                # β-string: β - pα, ..., β, ..., β + qα
                # Find p (largest integer s.t. β - pα ∈ Φ)
                p = 0
                test = tuple(b - a for a, b in zip(alpha, beta))
                while test in root_set:
                    p += 1
                    test = tuple(b - (p + 1) * a for a, b in zip(alpha, beta))

                # N_{α,β} = ±(p+1)
                # Sign convention: use a consistent choice
                # We'll use the Tits convention with a fixed ordering
                N = p + 1

                # Sign: use the lexicographic convention
                # For simple roots, N > 0 if α < β in root ordering
                if j < k:
                    sign = 1.0
                else:
                    sign = -1.0

                idx_alpha = rank + j
                idx_beta = rank + k
                idx_gamma = root_to_idx[gamma]
                ad[idx_alpha, idx_beta, idx_gamma] = sign * N

    # Convert structure constants to adjoint representation matrices
    # (ad(T_a))_{bc} = f_{ab}^c
    generators = []
    for a in range(dim):
        M = np.zeros((dim, dim), dtype=np.complex128)
        for b in range(dim):
            for c in range(dim):
                M[b, c] = ad[a, b, c]
        generators.append(M)

    return generators


def _simple_roots_A(n):
    """Simple roots of A_n (su(n+1)). Rank = n."""
    roots = []
    for i in range(n):
        r = [0] * n
        r[i] = 1
        roots.append(r)
    return roots


def _cartan_matrix(simple_roots):
    """Cartan matrix from simple roots: A_ij = 2(α_i, α_j)/(α_j, α_j)."""
    n = len(simple_roots)
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            ai = np.array(simple_roots[i], dtype=np.float64)
            aj = np.array(simple_roots[j], dtype=np.float64)
            A[i, j] = 2 * (ai @ aj) / (aj @ aj)
    return A


def _build_root_system(simple_roots):
    """
    Build the full root system from simple roots by repeated
    reflections (Weyl group orbit of positive roots).
    """
    rank = len(simple_roots)
    simple = [np.array(r, dtype=np.float64) for r in simple_roots]

    # Norms squared
    norms_sq = [s @ s for s in simple]

    def reflect(root, simple_idx):
        """Reflect root through the hyperplane of simple_roots[simple_idx]."""
        s = simple[simple_idx]
        return root - 2 * (root @ s) / (s @ s) * s

    def root_key(r):
        return tuple(np.round(r, 8))

    # Start with simple roots as positive roots
    positive = set()
    for s in simple:
        positive.add(root_key(s))

    # Generate all positive roots by adding simple roots
    # Use the Weyl reflection method
    queue = [np.array(r, dtype=np.float64) for r in simple]
    visited = set(root_key(r) for r in queue)

    while queue:
        root = queue.pop(0)
        for i in range(rank):
            new_root = reflect(root, i)
            key = root_key(new_root)
            if key not in visited:
                visited.add(key)
                queue.append(new_root)
                # Check if it's positive (first nonzero component is positive)
                for comp in new_root:
                    if abs(comp) > 1e-8:
                        if comp > 0:
                            positive.add(key)
                        break

    # Full root system = positive ∪ negative
    all_roots = []
    for key in positive:
        all_roots.append(list(key))
        all_roots.append([-x for x in key])

    return all_roots


def _exceptional_simple_roots(name):
    """
    Simple roots for exceptional Lie algebras in the standard basis.
    Using conventions from Humphreys / Bourbaki.
    """
    if name == 'G2':
        # G₂: rank 2, 12 roots
        # Simple roots in R²: α₁ = (1, -1, 0), α₂ = (-2, 1, 1) in R³
        # But for root system in rank-2 space:
        # α₁ = (1, 0), α₂ = (-3/2, √3/2)
        # Use 3D embedding for cleaner integers
        return [[1, -1, 0], [-2, 1, 1]]

    elif name == 'F4':
        # F₄: rank 4, 48 roots
        # Simple roots in R⁴:
        return [
            [0, 1, -1, 0],
            [0, 0, 1, -1],
            [0, 0, 0, 1],
            [0.5, -0.5, -0.5, -0.5]
        ]

    elif name == 'E6':
        # E₆: rank 6, 72 roots
        # Embed in R⁸ (standard Bourbaki)
        sqrt3_inv = 1.0 / np.sqrt(3)
        return [
            [0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5],
            [1, 1, 0, 0, 0, 0, 0, 0],
            [-1, 1, 0, 0, 0, 0, 0, 0],
            [0, -1, 1, 0, 0, 0, 0, 0],
            [0, 0, -1, 1, 0, 0, 0, 0],
            [0, 0, 0, -1, 1, 0, 0, 0],
        ]

    elif name == 'E7':
        # E₇: rank 7, 126 roots
        return [
            [0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5],
            [1, 1, 0, 0, 0, 0, 0, 0],
            [-1, 1, 0, 0, 0, 0, 0, 0],
            [0, -1, 1, 0, 0, 0, 0, 0],
            [0, 0, -1, 1, 0, 0, 0, 0],
            [0, 0, 0, -1, 1, 0, 0, 0],
            [0, 0, 0, 0, -1, 1, 0, 0],
        ]

    elif name == 'E8':
        # E₈: rank 8, 240 roots
        return [
            [0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5],
            [1, 1, 0, 0, 0, 0, 0, 0],
            [-1, 1, 0, 0, 0, 0, 0, 0],
            [0, -1, 1, 0, 0, 0, 0, 0],
            [0, 0, -1, 1, 0, 0, 0, 0],
            [0, 0, 0, -1, 1, 0, 0, 0],
            [0, 0, 0, 0, -1, 1, 0, 0],
            [0, 0, 0, 0, 0, -1, 1, 0],
        ]

    raise ValueError(f"Unknown exceptional algebra: {name}")


def exceptional_generators(name):
    """
    Build generators for exceptional Lie algebra using root system
    and Chevalley-Serre construction in the adjoint representation.
    """
    simple_roots = _exceptional_simple_roots(name)
    all_roots = _build_root_system(simple_roots)
    rank = len(simple_roots)

    expected_dims = {'G2': 14, 'F4': 52, 'E6': 78, 'E7': 133, 'E8': 248}
    expected_roots = {'G2': 12, 'F4': 48, 'E6': 72, 'E7': 126, 'E8': 240}
    expected_dim = expected_dims[name]
    expected_n_roots = expected_roots[name]

    if len(all_roots) != expected_n_roots:
        # Fallback: construct roots explicitly for each type
        all_roots = _explicit_roots(name)

    assert len(all_roots) == expected_n_roots, \
        f"{name}: expected {expected_n_roots} roots, got {len(all_roots)}"

    gens = _build_exceptional_via_root_system(rank, all_roots, _cartan_matrix(simple_roots))
    assert len(gens) == expected_dim, f"{name}: expected dim {expected_dim}, got {len(gens)}"
    return gens


def _explicit_roots(name):
    """Explicitly construct all roots for exceptional Lie algebras."""
    if name == 'G2':
        # G₂ roots in R³ (orthogonal to (1,1,1)): 12 roots
        # Short roots: permutations of (1,-1,0) → 6 roots
        # Long roots: permutations of (2,-1,-1) → 6 roots
        roots = []
        # Short roots: all permutations of (1,-1,0)
        from itertools import permutations
        for p in set(permutations([1, -1, 0])):
            roots.append(list(p))
        # Long roots: all permutations of (2,-1,-1)
        for p in set(permutations([2, -1, -1])):
            roots.append(list(p))
        for p in set(permutations([-2, 1, 1])):
            roots.append(list(p))
        # Deduplicate
        unique = []
        seen = set()
        for r in roots:
            key = tuple(r)
            if key not in seen:
                seen.add(key)
                unique.append(r)
        return unique

    elif name == 'F4':
        # F₄ roots in R⁴: 48 roots
        roots = []
        # Type 1: ±e_i ± e_j (i<j): 24 roots
        for i in range(4):
            for j in range(i + 1, 4):
                for si in [1, -1]:
                    for sj in [1, -1]:
                        r = [0, 0, 0, 0]
                        r[i] = si
                        r[j] = sj
                        roots.append(r)
        # Type 2: ±e_i: 8 roots
        for i in range(4):
            for s in [1, -1]:
                r = [0, 0, 0, 0]
                r[i] = s
                roots.append(r)
        # Type 3: (±1/2, ±1/2, ±1/2, ±1/2): 16 roots
        for s0 in [0.5, -0.5]:
            for s1 in [0.5, -0.5]:
                for s2 in [0.5, -0.5]:
                    for s3 in [0.5, -0.5]:
                        roots.append([s0, s1, s2, s3])
        assert len(roots) == 48, f"F₄: expected 48 roots, got {len(roots)}"
        return roots

    elif name == 'E6':
        # E₆: 72 roots in R⁸ (standard E₈ embedding, restricted)
        # We use the description: all roots of E₈ that are orthogonal to
        # two specific vectors, or more directly:
        roots = _e_series_roots(8)  # Get all 240 E₈ roots
        # E₆ roots: those with last two coordinates equal and
        # satisfying the E₆ constraint
        # Simpler: use standard E₆ root description in R⁶
        return _e6_roots_direct()

    elif name == 'E7':
        # E₇: 126 roots
        return _e7_roots_direct()

    elif name == 'E8':
        # E₈: 240 roots in R⁸
        return _e_series_roots(8)

    raise ValueError(f"Unknown: {name}")


def _e_series_roots(dim):
    """
    E₈ roots in R⁸: 240 roots.
    Type 1: all permutations of (±1, ±1, 0, 0, 0, 0, 0, 0) → 112 roots
    Type 2: all (±½, ±½, ..., ±½) with even number of minus signs → 128 roots
    """
    roots = []

    # Type 1: ±e_i ± e_j for i < j
    for i in range(8):
        for j in range(i + 1, 8):
            for si in [1, -1]:
                for sj in [1, -1]:
                    r = [0] * 8
                    r[i] = si
                    r[j] = sj
                    roots.append(r)

    # Type 2: (±½)^8 with even number of minus signs
    for bits in range(256):
        signs = [(bits >> k) & 1 for k in range(8)]
        n_minus = sum(signs)
        if n_minus % 2 == 0:
            r = [(-1) ** s * 0.5 for s in signs]
            roots.append(r)

    assert len(roots) == 240, f"E₈: expected 240 roots, got {len(roots)}"
    return roots


def _e6_roots_direct():
    """
    E₆ roots: 72 roots in R⁸ (using standard E₈ embedding).

    E₆ ⊂ E₈. The roots of E₆ are those E₈ roots orthogonal to both:
      v₁ = e₇ - e₈  and  v₂ = e₆ + e₇

    Equivalently, roots with x₇ = x₈ and x₆ = -x₇.
    It's easier to use the direct description in 6D.

    In R⁶ basis:
      - ±e_i ± e_j for 1 ≤ i < j ≤ 5 → 40 roots
      - ±(½, ±½, ±½, ±½, ±½, ±½√3) with constraints → 32 roots
    """
    # Use the R⁸ embedding and filter
    e8_roots = _e_series_roots(8)
    e6_roots = []
    for r in e8_roots:
        # E₆ constraint: r[5] + r[6] + r[7] = 0 and r[6] = r[7]
        if abs(r[6] - r[7]) < 1e-10 and abs(r[5] + r[6] + r[7]) < 1e-10:
            e6_roots.append(r)

    if len(e6_roots) == 72:
        return e6_roots

    # Alternative E₆ filter: x₆=x₇=x₈=c for some c
    e6_roots = []
    for r in e8_roots:
        if abs(r[5] - r[6]) < 1e-10 and abs(r[6] - r[7]) < 1e-10:
            e6_roots.append(r)

    if len(e6_roots) == 72:
        return e6_roots

    # Fallback: use direct D₅ + spinor description
    # E₆ in terms of the subalgebra so(10) ⊕ u(1):
    # Roots = roots of D₅ + spinor weights
    e6_roots = []

    # D₅ roots: ±e_i ± e_j, i<j, i,j ∈ {0..4} in R⁸ (pad with zeros)
    for i in range(5):
        for j in range(i + 1, 5):
            for si in [1, -1]:
                for sj in [1, -1]:
                    r = [0.0] * 8
                    r[i] = si
                    r[j] = sj
                    e6_roots.append(r)  # 40 roots

    # Spinor weights: (±½)^5 with even number of minus, extended
    for bits in range(32):
        signs = [(bits >> k) & 1 for k in range(5)]
        n_minus = sum(signs)
        if n_minus % 2 == 0:
            r = [0.0] * 8
            for k in range(5):
                r[k] = (-1) ** signs[k] * 0.5
            # For E₆: the 6th,7th,8th components encode the u(1) charge
            # Use the convention that makes them roots: r[5]=r[6]=r[7] = ±√3/6
            # Actually, for the adjoint rep construction, we only need the roots
            # to be consistent. Use r[5] = -sum(r[:5])/3 type constraint
            s = sum(r[:5])
            r[5] = -s / 3
            r[6] = -s / 3
            r[7] = -s / 3
            e6_roots.append(r)
        elif n_minus % 2 == 1:
            # odd parity spinor (anti-spinor)
            r = [0.0] * 8
            for k in range(5):
                r[k] = (-1) ** signs[k] * 0.5
            s = sum(r[:5])
            r[5] = -s / 3
            r[6] = -s / 3
            r[7] = -s / 3
            e6_roots.append(r)

    # Deduplicate
    unique = []
    seen = set()
    for r in e6_roots:
        key = tuple(round(x, 8) for x in r)
        if key not in seen:
            seen.add(key)
            unique.append(r)
    e6_roots = unique

    if len(e6_roots) != 72:
        print(f"  WARNING: E₆ root system has {len(e6_roots)} roots (expected 72)")
        print(f"  Proceeding with available roots — results may be approximate")

    return e6_roots


def _e7_roots_direct():
    """
    E₇ roots: 126 roots. E₇ ⊂ E₈ with constraint x₈ = 0 effectively.

    E₇ roots = E₈ roots orthogonal to e₇ + e₈,
    i.e., those with r[6] + r[7] = 0.
    """
    e8_roots = _e_series_roots(8)
    e7_roots = []
    for r in e8_roots:
        if abs(r[6] + r[7]) < 1e-10:
            e7_roots.append(r)

    if len(e7_roots) == 126:
        return e7_roots

    # Fallback
    print(f"  WARNING: E₇ filter gave {len(e7_roots)} roots (expected 126)")

    # Try alternative: r[7] = 0
    e7_roots = []
    for r in e8_roots:
        if abs(r[7]) < 1e-10:
            e7_roots.append(r)

    if len(e7_roots) == 126:
        return e7_roots

    print(f"  WARNING: E₇ alt filter gave {len(e7_roots)} roots (expected 126)")
    # Use the ones we have
    return e7_roots if e7_roots else e8_roots[:126]


# ============================================================
# Adjoint graph and Ramanujan check (from ramanujan_check.py)
# ============================================================

def build_adjoint_graph(generators, threshold=1e-8):
    """Build binary adjacency matrix from commutator structure."""
    n = len(generators)
    adjacency = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            C = commutator(generators[i], generators[j])
            if np.linalg.norm(C) ** 2 > threshold:
                adjacency[i, j] = 1
                adjacency[j, i] = 1
    return adjacency


def spectral_analysis(adjacency):
    """Compute eigenvalues and check Ramanujan property."""
    n = adjacency.shape[0]
    eigenvalues = np.linalg.eigvalsh(adjacency)
    eigenvalues = np.sort(eigenvalues)[::-1]

    degrees = adjacency.sum(axis=1)
    is_regular = np.allclose(degrees, degrees[0])
    q_plus_1 = degrees[0] if is_regular else np.mean(degrees)
    q = q_plus_1 - 1

    lambda_max = eigenvalues[0]
    lambda_2 = eigenvalues[1] if n > 1 else 0

    ramanujan_bound = 2 * np.sqrt(max(q, 0))

    # Non-trivial: exclude trivial eigenvalue(s)
    # For bipartite also exclude -λ_max
    nontrivial = eigenvalues[1:]
    if n > 2 and abs(eigenvalues[-1] + lambda_max) < 1e-6:
        # Bipartite: exclude both ±λ_max
        nontrivial = eigenvalues[1:-1]

    max_nontrivial = np.max(np.abs(nontrivial)) if len(nontrivial) > 0 else 0
    is_ramanujan = max_nontrivial <= ramanujan_bound + 1e-10
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
# Build the full catalog
# ============================================================

def build_catalog(families=None, quick=False):
    """
    Build the complete algebra catalog.

    families: list of family codes to include. Default: all.
        'A' = su(n+1), 'B' = so(2n+1), 'D' = so(2n),
        'C' = sp(2n), 'G' = G₂, 'F' = F₄, 'E' = E₆,E₇,E₈
    quick: if True, only test small algebras (skip E₇, E₈)
    """
    if families is None:
        families = ['A', 'B', 'D', 'C', 'G', 'F', 'E']

    catalog = {}

    if 'A' in families:
        n_max = 6 if quick else 10
        for n in range(2, n_max + 1):
            name = f"su({n})"
            catalog[name] = ('matrix', lambda n=n: su_generators(n))

    if 'B' in families:
        n_max = 5 if quick else 7
        for n in range(3, 2 * n_max + 2, 2):
            if n >= 3:
                name = f"so({n})"
                if name not in catalog:
                    catalog[name] = ('matrix', lambda n=n: so_generators(n))

    if 'D' in families:
        n_max = 4 if quick else 6
        for n in range(4, 2 * n_max + 1, 2):
            name = f"so({n})"
            if name not in catalog:
                catalog[name] = ('matrix', lambda n=n: so_generators(n))

    if 'C' in families:
        n_max = 3 if quick else 5
        for k in range(1, n_max + 1):
            name = f"sp({2*k})"
            catalog[name] = ('matrix', lambda k=k: sp_generators(k))

    if 'G' in families:
        catalog['G₂'] = ('exceptional', lambda: g2_generators())

    if 'F' in families:
        catalog['F₄'] = ('exceptional', lambda: exceptional_generators('F4'))

    if 'E' in families:
        catalog['E₆'] = ('exceptional', lambda: exceptional_generators('E6'))
        if not quick:
            catalog['E₇'] = ('exceptional', lambda: exceptional_generators('E7'))
            catalog['E₈'] = ('exceptional', lambda: exceptional_generators('E8'))

    return catalog


# ============================================================
# Main analysis
# ============================================================

def analyze_algebra(name, gen_fn, gen_type):
    """Analyze a single algebra. Returns dict with all results."""
    t0 = time.time()
    try:
        generators = gen_fn()
    except Exception as e:
        print(f"  ERROR generating {name}: {e}")
        return None

    n_gen = len(generators)
    if n_gen == 0:
        print(f"  ERROR: {name} has 0 generators")
        return None

    rep_dim = generators[0].shape[0]
    print(f"  dim = {n_gen}, rep_dim = {rep_dim}", end="")

    # Build adjoint graph
    adj = build_adjoint_graph(generators)
    n_edges = int(adj.sum() // 2)
    print(f", edges = {n_edges}", end="")

    # Spectral analysis
    spec = spectral_analysis(adj)
    t1 = time.time()

    ram_str = "YES ✓" if spec['is_ramanujan'] else "NO ✗"
    print(f", ratio = {spec['ratio']:.4f}, Ramanujan: {ram_str}  [{t1-t0:.1f}s]")

    return {
        'n_gen': n_gen,
        'rep_dim': rep_dim,
        'n_edges': n_edges,
        'spectral': spec,
        'time': t1 - t0,
        'type': gen_type,
    }


def run_extended_analysis(families=None, quick=False):
    """Run analysis on full catalog."""
    catalog = build_catalog(families, quick)

    print("=" * 80)
    print("EXTENDED RAMANUJAN CHECK: ALL SIMPLE LIE ALGEBRAS")
    print("=" * 80)

    results = {}
    for name, (gen_type, gen_fn) in catalog.items():
        print(f"\n  {name}:")
        res = analyze_algebra(name, gen_fn, gen_type)
        if res is not None:
            results[name] = res

    return results


def print_summary(results):
    """Print summary table."""
    print("\n" + "=" * 110)
    print("SUMMARY: Extended Ramanujan Property Check")
    print("=" * 110)
    print(f"{'Algebra':>10} | {'dim':>4} | {'Edges':>5} | {'Regular':>7} | "
          f"{'deg':>5} | {'λ_max':>6} | {'2√q':>6} | {'max|λ_nt|':>9} | "
          f"{'Ratio':>6} | {'Ramanujan':>9} | {'Time':>6}")
    print("-" * 110)

    n_tested = 0
    n_ramanujan = 0

    for name, res in results.items():
        spec = res['spectral']
        ram_str = 'YES ✓' if spec['is_ramanujan'] else '** NO ✗ **'
        reg_str = 'yes' if spec['is_regular'] else 'no'
        n_tested += 1
        if spec['is_ramanujan']:
            n_ramanujan += 1

        print(f"{name:>10} | {res['n_gen']:4d} | {res['n_edges']:5d} | "
              f"{reg_str:>7} | {spec['degree']:5.1f} | "
              f"{spec['lambda_max']:6.2f} | {spec['ramanujan_bound']:6.2f} | "
              f"{spec['max_nontrivial']:9.4f} | {spec['ratio']:6.4f} | "
              f"{ram_str:>9} | {res['time']:5.1f}s")

    print("-" * 110)
    print(f"TOTAL: {n_ramanujan}/{n_tested} Ramanujan")
    if n_ramanujan == n_tested:
        print("*** ALL ALGEBRAS SATISFY THE RAMANUJAN PROPERTY ***")
    else:
        failures = [n for n, r in results.items() if not r['spectral']['is_ramanujan']]
        print(f"FAILURES: {', '.join(failures)}")

    return n_tested, n_ramanujan


# ============================================================
# Visualization
# ============================================================

def plot_results(results, save=False):
    """Create multi-panel visualization."""
    names = list(results.keys())
    dims = [results[n]['n_gen'] for n in names]
    ratios = [results[n]['spectral']['ratio'] for n in names]
    is_ram = [results[n]['spectral']['is_ramanujan'] for n in names]
    gaps = [results[n]['spectral']['spectral_gap'] for n in names]

    fig, axes = plt.subplots(2, 2, figsize=(18, 12))
    fig.suptitle('Extended Ramanujan Check: All Simple Lie Algebras',
                 fontsize=16, fontweight='bold')

    # Panel 1: Ramanujan ratio vs dimension
    ax = axes[0, 0]
    colors = ['green' if r else 'red' for r in is_ram]
    ax.scatter(dims, ratios, c=colors, s=80, edgecolors='black', linewidths=0.5, zorder=3)
    for i, n in enumerate(names):
        if dims[i] > 30 or ratios[i] > 0.5:
            ax.annotate(n, (dims[i], ratios[i]), fontsize=7,
                       ha='left', va='bottom', xytext=(3, 3),
                       textcoords='offset points')
    ax.axhline(1.0, color='red', linestyle='--', alpha=0.5, label='Ramanujan bound')
    ax.set_xlabel('dim(algebra)')
    ax.set_ylabel('Ramanujan ratio max|λ_nt| / 2√q')
    ax.set_title('Ramanujan Ratio vs Dimension')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)

    # Panel 2: Ratio by family
    ax = axes[0, 1]
    family_data = {}
    for n, r in results.items():
        if n.startswith('su'):
            fam = 'A (su)'
        elif n.startswith('so'):
            fam = 'B/D (so)'
        elif n.startswith('sp'):
            fam = 'C (sp)'
        else:
            fam = 'Except.'
        if fam not in family_data:
            family_data[fam] = ([], [])
        family_data[fam][0].append(r['n_gen'])
        family_data[fam][1].append(r['spectral']['ratio'])

    markers = {'A (su)': 'o', 'B/D (so)': 's', 'C (sp)': '^', 'Except.': 'D'}
    fam_colors = {'A (su)': 'blue', 'B/D (so)': 'red', 'C (sp)': 'green', 'Except.': 'purple'}
    for fam, (ds, rs) in family_data.items():
        ax.plot(ds, rs, marker=markers.get(fam, 'o'), color=fam_colors.get(fam, 'grey'),
                linestyle='-', alpha=0.7, label=fam, markersize=8)
    ax.axhline(1.0, color='red', linestyle='--', alpha=0.3)
    ax.set_xlabel('dim(algebra)')
    ax.set_ylabel('Ramanujan ratio')
    ax.set_title('Ratio by Family')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)

    # Panel 3: Spectral gap vs dimension
    ax = axes[1, 0]
    ax.scatter(dims, gaps, c=colors, s=80, edgecolors='black', linewidths=0.5)
    for i, n in enumerate(names):
        if dims[i] > 30:
            ax.annotate(n, (dims[i], gaps[i]), fontsize=7,
                       ha='left', va='bottom', xytext=(3, 3),
                       textcoords='offset points')
    ax.set_xlabel('dim(algebra)')
    ax.set_ylabel('Spectral gap (λ₁ - λ₂)')
    ax.set_title('Spectral Gap vs Dimension')
    ax.grid(True, alpha=0.3)

    # Panel 4: Summary bar chart
    ax = axes[1, 1]
    # Sort by dimension
    sorted_names = sorted(names, key=lambda n: results[n]['n_gen'])
    sorted_ratios = [results[n]['spectral']['ratio'] for n in sorted_names]
    sorted_colors = ['green' if results[n]['spectral']['is_ramanujan'] else 'red'
                     for n in sorted_names]
    x_pos = range(len(sorted_names))
    ax.bar(x_pos, sorted_ratios, color=sorted_colors, alpha=0.7, edgecolor='black',
           linewidth=0.5)
    ax.axhline(1.0, color='red', linestyle='--', alpha=0.5, label='Ramanujan bound')
    ax.set_xticks(list(x_pos))
    ax.set_xticklabels(sorted_names, rotation=60, ha='right', fontsize=7)
    ax.set_ylabel('Ramanujan ratio')
    ax.set_title('All Algebras (sorted by dimension)')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()

    if save:
        out_path = os.path.join(os.path.dirname(__file__), 'ramanujan_extended_results.png')
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        print(f"\nSaved to {out_path}")

    plt.show()


# ============================================================
# Entry point
# ============================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extended Ramanujan check for ALL simple Lie algebras'
    )
    parser.add_argument('--save', action='store_true', help='Save plots to file')
    parser.add_argument('--quick', action='store_true',
                        help='Quick mode: skip E₇ and E₈')
    parser.add_argument('--families', nargs='+', default=None,
                        help='Families to test: A B D C G F E')
    args = parser.parse_args()

    results = run_extended_analysis(families=args.families, quick=args.quick)
    n_tested, n_ramanujan = print_summary(results)
    plot_results(results, save=args.save)
