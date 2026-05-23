"""Quick test: Ramanujan in physics basis for su(n), so(n), sp(2n)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'verification'))
from shared_utils import commutator
from ramanujan_extended import (su_generators, adjacency_from_matrices,
                                spectral_analysis, roots_B, roots_C, roots_D,
                                adjacency_from_roots)


def so_generators(n):
    """so(n) generators: antisymmetric n×n matrices."""
    gens = []
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n), dtype=np.complex128)
            E[i, j] = 1
            E[j, i] = -1
            gens.append(E)
    return gens


def sp_generators(n):
    """sp(2n) generators: 2n×2n matrices X s.t. X^T J + J X = 0, J=[[0,I],[-I,0]]."""
    gens = []
    # Symmetric: e_ij + e_ji blocks
    for i in range(n):
        for j in range(i, n):
            # Upper-right block (symmetric n×n)
            E = np.zeros((2 * n, 2 * n), dtype=np.complex128)
            E[i, n + j] = 1
            E[j, n + i] = 1
            gens.append(E)
            # Lower-left block (symmetric n×n)
            E2 = np.zeros((2 * n, 2 * n), dtype=np.complex128)
            E2[n + i, j] = 1
            E2[n + j, i] = 1
            gens.append(E2)
    # Diagonal blocks: e_ij - e_{j+n,i+n} for the gl(n) part
    for i in range(n):
        for j in range(n):
            E = np.zeros((2 * n, 2 * n), dtype=np.complex128)
            E[i, j] = 1
            E[n + j, n + i] = -1
            if np.linalg.norm(E) > 1e-10:
                gens.append(E)
    return gens


def test_family(name_fmt, gen_fn, dim_fn, params, label):
    """Test Ramanujan for a family of algebras in physics basis."""
    print(f"\n{label}:")
    for p in params:
        gens = gen_fn(p)
        adj = adjacency_from_matrices(gens)
        s = spectral_analysis(adj)
        dim = len(gens)
        ram = "YES" if s["is_ramanujan"] else "NO"
        print(f"  {name_fmt + '(' + str(p) + ')':>8}: dim={dim:3d}, deg={s['degree']:5.1f}, "
              f"ratio={s['ratio']:.4f}, Ramanujan: {ram}")


# Extended Cartan-Weyl test for higher ranks
def test_cartan_weyl_extended():
    """Test B, C, D families at higher ranks in Cartan-Weyl basis."""
    print("\nCartan-Weyl basis - extended ranks:")
    for k in range(2, 9):
        rank, roots, simple = roots_B(k)
        adj = adjacency_from_roots(rank, roots, simple)
        s = spectral_analysis(adj)
        dim = rank + len(roots)
        ram = "YES" if s["is_ramanujan"] else "NO"
        print(f"  so({2*k+1:2d}): dim={dim:3d}, deg={s['degree']:5.1f}, "
              f"ratio={s['ratio']:.4f}, Ramanujan: {ram}")

    for k in range(2, 9):
        rank, roots, simple = roots_C(k)
        adj = adjacency_from_roots(rank, roots, simple)
        s = spectral_analysis(adj)
        dim = rank + len(roots)
        ram = "YES" if s["is_ramanujan"] else "NO"
        print(f"  sp({2*k:2d}): dim={dim:3d}, deg={s['degree']:5.1f}, "
              f"ratio={s['ratio']:.4f}, Ramanujan: {ram}")

    for k in range(3, 10):
        rank, roots, simple = roots_D(k)
        adj = adjacency_from_roots(rank, roots, simple)
        s = spectral_analysis(adj)
        dim = rank + len(roots)
        ram = "YES" if s["is_ramanujan"] else "NO"
        print(f"  so({2*k:2d}): dim={dim:3d}, deg={s['degree']:5.1f}, "
              f"ratio={s['ratio']:.4f}, Ramanujan: {ram}")


# Run tests
print("=" * 70)
print("  PHYSICS BASIS vs CARTAN-WEYL: Ramanujan comparison")
print("=" * 70)

test_family("su", su_generators, lambda n: n*n-1, range(2, 11),
            "su(n) - Physics basis (Gell-Mann)")

test_family("so", so_generators, lambda n: n*(n-1)//2, range(3, 16),
            "so(n) - Physics basis (antisymmetric)")

test_cartan_weyl_extended()
