"""
Extended Ramanujan Check for ALL Simple Lie Algebras
=====================================================

Tests the Ramanujan property of the adjoint commutation graph for:
  - A_n = su(n+1)  for n=1..9   (dim 3..99)
  - B_n = so(2n+1) for n=2..5   (dim 10..55)
  - D_n = so(2n)   for n=3..6   (dim 15..66)
  - C_n = sp(2n)   for n=2..5   (dim 10..55)
  - G2  (dim 14)
  - F4  (dim 52)
  - E6  (dim 78)
  - E7  (dim 133)
  - E8  (dim 248)

KEY INSIGHT: The adjacency matrix of the adjoint graph depends ONLY on
the root system. No matrix representations needed.

  Vertices = {H_1,...,H_r, E_{a1},...,E_{a|Phi|}}
  H_i -- H_j : NEVER (Cartan subalgebra is abelian)
  H_i -- E_a : iff (a, simple_root_i) != 0
  E_a -- E_b : iff a + b in Phi  OR  b = -a

This makes E8 (248 generators, 240 roots) trivial to compute.

Usage:
    python ramanujan_extended.py [--save] [--quick]
    python ramanujan_extended.py --families A B D C G F E
    python ramanujan_extended.py --families E  (exceptionals only)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import sys
import os
import time
from collections import OrderedDict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'verification'))
from shared_utils import commutator


# ============================================================
# Root systems for all simple Lie algebras
# ============================================================

def roots_A(n):
    """A_n root system: su(n+1). Roots = e_i - e_j in R^{n+1}."""
    rank = n
    dim_ambient = n + 1
    roots = []
    for i in range(dim_ambient):
        for j in range(dim_ambient):
            if i != j:
                r = [0] * dim_ambient
                r[i] = 1
                r[j] = -1
                roots.append(tuple(r))
    simple = []
    for i in range(n):
        r = [0] * dim_ambient
        r[i] = 1
        r[i + 1] = -1
        simple.append(tuple(r))
    return rank, roots, simple


def roots_B(n):
    """B_n root system: so(2n+1). Roots = +-e_i+-e_j (i<j) and +-e_i."""
    rank = n
    roots = []
    for i in range(n):
        for j in range(i + 1, n):
            for si in [1, -1]:
                for sj in [1, -1]:
                    r = [0] * n
                    r[i] = si
                    r[j] = sj
                    roots.append(tuple(r))
    for i in range(n):
        for s in [1, -1]:
            r = [0] * n
            r[i] = s
            roots.append(tuple(r))
    simple = []
    for i in range(n - 1):
        r = [0] * n
        r[i] = 1
        r[i + 1] = -1
        simple.append(tuple(r))
    r = [0] * n
    r[n - 1] = 1
    simple.append(tuple(r))
    return rank, roots, simple


def roots_C(n):
    """C_n root system: sp(2n). Roots = +-e_i+-e_j (i<j) and +-2e_i."""
    rank = n
    roots = []
    for i in range(n):
        for j in range(i + 1, n):
            for si in [1, -1]:
                for sj in [1, -1]:
                    r = [0] * n
                    r[i] = si
                    r[j] = sj
                    roots.append(tuple(r))
    for i in range(n):
        for s in [2, -2]:
            r = [0] * n
            r[i] = s
            roots.append(tuple(r))
    simple = []
    for i in range(n - 1):
        r = [0] * n
        r[i] = 1
        r[i + 1] = -1
        simple.append(tuple(r))
    r = [0] * n
    r[n - 1] = 2
    simple.append(tuple(r))
    return rank, roots, simple


def roots_D(n):
    """D_n root system: so(2n). Roots = +-e_i+-e_j (i<j). Requires n >= 2."""
    rank = n
    roots = []
    for i in range(n):
        for j in range(i + 1, n):
            for si in [1, -1]:
                for sj in [1, -1]:
                    r = [0] * n
                    r[i] = si
                    r[j] = sj
                    roots.append(tuple(r))
    simple = []
    for i in range(n - 1):
        r = [0] * n
        r[i] = 1
        r[i + 1] = -1
        simple.append(tuple(r))
    r = [0] * n
    r[n - 2] = 1
    r[n - 1] = 1
    simple.append(tuple(r))
    return rank, roots, simple


def roots_G2():
    """G2 root system: 14-dimensional algebra, 12 roots, rank 2."""
    from itertools import permutations
    roots = []
    seen = set()
    for p in permutations([1, -1, 0]):
        t = tuple(p)
        if t not in seen:
            seen.add(t)
            roots.append(t)
    for p in permutations([2, -1, -1]):
        t = tuple(p)
        if t not in seen:
            seen.add(t)
            roots.append(t)
    for p in permutations([-2, 1, 1]):
        t = tuple(p)
        if t not in seen:
            seen.add(t)
            roots.append(t)
    assert len(roots) == 12, f"G2: expected 12 roots, got {len(roots)}"
    simple = [(1, -1, 0), (-2, 1, 1)]
    return 2, roots, simple


def roots_F4():
    """F4 root system: 52-dimensional algebra, 48 roots, rank 4."""
    roots = []
    for i in range(4):
        for j in range(i + 1, 4):
            for si in [1, -1]:
                for sj in [1, -1]:
                    r = [0, 0, 0, 0]
                    r[i] = si
                    r[j] = sj
                    roots.append(tuple(r))
    for i in range(4):
        for s in [1, -1]:
            r = [0, 0, 0, 0]
            r[i] = s
            roots.append(tuple(r))
    for s0 in [0.5, -0.5]:
        for s1 in [0.5, -0.5]:
            for s2 in [0.5, -0.5]:
                for s3 in [0.5, -0.5]:
                    roots.append((s0, s1, s2, s3))
    assert len(roots) == 48, f"F4: expected 48 roots, got {len(roots)}"
    simple = [
        (0, 1, -1, 0),
        (0, 0, 1, -1),
        (0, 0, 0, 1),
        (0.5, -0.5, -0.5, -0.5),
    ]
    return 4, roots, simple


def roots_E8():
    """E8 root system: 248-dimensional algebra, 240 roots, rank 8."""
    roots = []
    for i in range(8):
        for j in range(i + 1, 8):
            for si in [1, -1]:
                for sj in [1, -1]:
                    r = [0] * 8
                    r[i] = si
                    r[j] = sj
                    roots.append(tuple(r))
    for bits in range(256):
        signs = [(bits >> k) & 1 for k in range(8)]
        if sum(signs) % 2 == 0:
            r = tuple((-1) ** s * 0.5 for s in signs)
            roots.append(r)
    assert len(roots) == 240, f"E8: expected 240 roots, got {len(roots)}"
    simple = [
        (0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5),
        (1, 1, 0, 0, 0, 0, 0, 0),
        (-1, 1, 0, 0, 0, 0, 0, 0),
        (0, -1, 1, 0, 0, 0, 0, 0),
        (0, 0, -1, 1, 0, 0, 0, 0),
        (0, 0, 0, -1, 1, 0, 0, 0),
        (0, 0, 0, 0, -1, 1, 0, 0),
        (0, 0, 0, 0, 0, -1, 1, 0),
    ]
    return 8, roots, simple


def roots_E7():
    """E7 root system: 133-dimensional algebra, 126 roots, rank 7."""
    _, e8_roots, _ = roots_E8()
    roots = [r for r in e8_roots if abs(r[6] + r[7]) < 1e-10]
    if len(roots) != 126:
        roots = [r for r in e8_roots if abs(r[6] - r[7]) < 1e-10]
    assert len(roots) == 126, f"E7: expected 126 roots, got {len(roots)}"
    simple = [
        (0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5),
        (1, 1, 0, 0, 0, 0, 0, 0),
        (-1, 1, 0, 0, 0, 0, 0, 0),
        (0, -1, 1, 0, 0, 0, 0, 0),
        (0, 0, -1, 1, 0, 0, 0, 0),
        (0, 0, 0, -1, 1, 0, 0, 0),
        (0, 0, 0, 0, -1, 1, 0, 0),
    ]
    return 7, roots, simple


def roots_E6():
    """E6 root system: 78-dimensional algebra, 72 roots, rank 6."""
    _, e8_roots, _ = roots_E8()

    # E6 roots: those E8 roots with x6=x7=x8 (last three components equal)
    roots = [r for r in e8_roots
             if abs(r[5] - r[6]) < 1e-10 and abs(r[6] - r[7]) < 1e-10]
    if len(roots) != 72:
        # E6 via constraint: perpendicular to (0,0,0,0,0,1,-1,0) and (0,0,0,0,0,0,1,-1)
        roots = [r for r in e8_roots
                 if abs(r[5] - r[6]) < 1e-10 and abs(r[6] - r[7]) < 1e-10]
    if len(roots) != 72:
        # Last resort: perp to e6-e7 and e6+e7+e8 (i.e. r5+r6+r7=0 and r5=r6)
        roots = [r for r in e8_roots
                 if abs(r[5] - r[6]) < 1e-10 and abs(r[5] + r[6] + r[7]) < 1e-10]
    assert len(roots) == 72, f"E6: expected 72 roots, got {len(roots)}"
    simple = [
        (0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5),
        (1, 1, 0, 0, 0, 0, 0, 0),
        (-1, 1, 0, 0, 0, 0, 0, 0),
        (0, -1, 1, 0, 0, 0, 0, 0),
        (0, 0, -1, 1, 0, 0, 0, 0),
        (0, 0, 0, -1, 1, 0, 0, 0),
    ]
    return 6, roots, simple


# ============================================================
# Adjacency matrix from root system (no matrices needed!)
# ============================================================

def adjacency_from_roots(rank, roots, simple_roots):
    """
    Build the adjacency matrix of the adjoint graph DIRECTLY from
    the root system. No matrix representations needed.

    Generators indexed as: H_0,...,H_{r-1}, E_{root_0},...,E_{root_{|Phi|-1}}
    dim = rank + |roots|

    Adjacency rules:
      H_i -- H_j  : 0  (Cartan subalgebra is abelian)
      H_i -- E_a  : 1  iff (a . simple_root_i) != 0
      E_a -- E_b  : 1  iff (a + b) in root_set  OR  b = -a
    """
    n_roots = len(roots)
    dim = rank + n_roots

    root_set = set()
    for r in roots:
        key = tuple(round(x, 6) for x in r)
        root_set.add(key)

    simple_arr = [np.array(s, dtype=np.float64) for s in simple_roots]
    roots_arr = [np.array(r, dtype=np.float64) for r in roots]

    adj = np.zeros((dim, dim), dtype=np.float64)

    # H_i -- E_alpha
    for i in range(rank):
        for j in range(n_roots):
            if abs(np.dot(roots_arr[j], simple_arr[i])) > 1e-10:
                adj[i, rank + j] = 1
                adj[rank + j, i] = 1

    # E_alpha -- E_beta
    for j in range(n_roots):
        for k in range(j + 1, n_roots):
            s = roots_arr[j] + roots_arr[k]
            if np.linalg.norm(s) < 1e-10:
                # beta = -alpha
                adj[rank + j, rank + k] = 1
                adj[rank + k, rank + j] = 1
            else:
                gamma = tuple(round(x, 6) for x in s)
                if gamma in root_set:
                    adj[rank + j, rank + k] = 1
                    adj[rank + k, rank + j] = 1

    return adj


# ============================================================
# Matrix-based adjacency (for cross-validation on small algebras)
# ============================================================

def su_generators(n):
    """su(n) generators: n^2-1 traceless hermitian matrices."""
    gens = []
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n), dtype=np.complex128)
            E[i, j] = 1; E[j, i] = 1
            gens.append(E / 2)
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n), dtype=np.complex128)
            E[i, j] = -1j; E[j, i] = 1j
            gens.append(E / 2)
    for k in range(1, n):
        D = np.zeros((n, n), dtype=np.complex128)
        for i in range(k):
            D[i, i] = 1
        D[k, k] = -k
        D /= np.sqrt(2 * k * (k + 1))
        gens.append(D)
    return gens


def adjacency_from_matrices(generators, threshold=1e-8):
    """Build adjacency from explicit matrix generators."""
    n = len(generators)
    adj = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            C = commutator(generators[i], generators[j])
            if np.linalg.norm(C) ** 2 > threshold:
                adj[i, j] = 1
                adj[j, i] = 1
    return adj


# ============================================================
# Spectral analysis
# ============================================================

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
    lambda_2 = eigenvalues[1] if n > 1 else 0.0
    ramanujan_bound = 2 * np.sqrt(max(q, 0))

    nontrivial = eigenvalues[1:]
    if n > 2 and is_regular and abs(eigenvalues[-1] + lambda_max) < 0.01:
        nontrivial = eigenvalues[1:-1]

    max_nontrivial = np.max(np.abs(nontrivial)) if len(nontrivial) > 0 else 0.0
    is_ramanujan = max_nontrivial <= ramanujan_bound + 1e-10

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
        'spectral_gap': lambda_max - lambda_2,
        'ratio': max_nontrivial / ramanujan_bound if ramanujan_bound > 0 else float('inf')
    }


# ============================================================
# Complete catalog
# ============================================================

def build_catalog(families=None, quick=False):
    """Build ordered catalog of all algebras to test."""
    if families is None:
        families = ['A', 'B', 'C', 'D', 'G', 'F', 'E']

    catalog = OrderedDict()

    if 'A' in families:
        n_max = 6 if quick else 10
        for n in range(2, n_max + 1):
            name = f"su({n})"
            catalog[name] = {
                'root_fn': lambda nn=n: roots_A(nn - 1),
                'expected_dim': n * n - 1,
                'family': 'A',
            }

    if 'B' in families:
        n_max = 4 if quick else 5
        for k in range(2, n_max + 1):
            name = f"so({2*k+1})"
            catalog[name] = {
                'root_fn': lambda kk=k: roots_B(kk),
                'expected_dim': k * (2 * k + 1),
                'family': 'B',
            }

    if 'C' in families:
        n_max = 3 if quick else 5
        for k in range(2, n_max + 1):
            name = f"sp({2*k})"
            catalog[name] = {
                'root_fn': lambda kk=k: roots_C(kk),
                'expected_dim': k * (2 * k + 1),
                'family': 'C',
            }

    if 'D' in families:
        n_max = 4 if quick else 6
        for k in range(3, n_max + 1):
            name = f"so({2*k})"
            catalog[name] = {
                'root_fn': lambda kk=k: roots_D(kk),
                'expected_dim': k * (2 * k - 1),
                'family': 'D',
            }

    if 'G' in families:
        catalog['G2'] = {
            'root_fn': roots_G2,
            'expected_dim': 14,
            'family': 'G',
        }

    if 'F' in families:
        catalog['F4'] = {
            'root_fn': roots_F4,
            'expected_dim': 52,
            'family': 'F',
        }

    if 'E' in families:
        catalog['E6'] = {
            'root_fn': roots_E6,
            'expected_dim': 78,
            'family': 'E',
        }
        if not quick:
            catalog['E7'] = {
                'root_fn': roots_E7,
                'expected_dim': 133,
                'family': 'E',
            }
            catalog['E8'] = {
                'root_fn': roots_E8,
                'expected_dim': 248,
                'family': 'E',
            }

    return catalog


# ============================================================
# Analysis
# ============================================================

def analyze_algebra(name, info):
    """Analyze a single algebra. Returns result dict or None on error."""
    t0 = time.time()

    try:
        rank, roots, simple = info['root_fn']()
    except Exception as e:
        print(f"  ERROR building root system for {name}: {e}")
        return None

    n_roots = len(roots)
    dim = rank + n_roots
    expected = info['expected_dim']

    if dim != expected:
        print(f"  WARNING: {name} has dim={dim} (expected {expected})")

    adj = adjacency_from_roots(rank, roots, simple)
    n_edges = int(adj.sum() // 2)

    spec = spectral_analysis(adj)
    t1 = time.time()

    ram_str = "YES" if spec['is_ramanujan'] else "NO"
    reg_str = "reg" if spec['is_regular'] else "non-reg"
    print(f"  {name:>8}: dim={dim:3d}, edges={n_edges:5d}, {reg_str}, "
          f"deg={spec['degree']:.1f}, ratio={spec['ratio']:.4f}, "
          f"Ramanujan: {ram_str}  [{t1-t0:.2f}s]")

    return {
        'dim': dim,
        'rank': rank,
        'n_roots': n_roots,
        'n_edges': n_edges,
        'spectral': spec,
        'family': info['family'],
        'time': t1 - t0,
    }


def cross_validate_small(results):
    """Cross-validate root-based adjacency with matrix-based for su(2..4)."""
    print("\n  Cross-validation (root system vs matrix representation):")
    for n in [2, 3, 4]:
        name = f"su({n})"
        if name not in results:
            continue
        gens = su_generators(n)
        adj_mat = adjacency_from_matrices(gens)
        rank_n, roots_n, simple_n = roots_A(n - 1)
        adj_root = adjacency_from_roots(rank_n, roots_n, simple_n)

        spec_mat = np.sort(np.linalg.eigvalsh(adj_mat))[::-1]
        spec_root = np.sort(np.linalg.eigvalsh(adj_root))[::-1]

        if np.allclose(spec_mat, spec_root, atol=1e-6):
            print(f"    {name}: spectra MATCH")
        else:
            print(f"    {name}: spectra DIFFER!")
            print(f"      Matrix:  {spec_mat[:6]}...")
            print(f"      Root:    {spec_root[:6]}...")


def run_analysis(families=None, quick=False):
    """Run the full extended analysis."""
    catalog = build_catalog(families, quick)

    print("=" * 80)
    print("  EXTENDED RAMANUJAN CHECK: ALL SIMPLE LIE ALGEBRAS")
    print("  Method: adjacency from root system (exact, no matrices)")
    print("=" * 80)

    results = OrderedDict()
    for name, info in catalog.items():
        res = analyze_algebra(name, info)
        if res is not None:
            results[name] = res

    return results


def print_summary(results):
    """Print formatted summary table."""
    print("\n" + "=" * 115)
    print("  SUMMARY: Ramanujan Property of Adjoint Commutation Graphs")
    print("=" * 115)
    hdr = (f"{'Algebra':>10} | {'Family':>6} | {'dim':>4} | {'Edges':>5} | "
           f"{'Regular':>7} | {'deg':>5} | {'lam_max':>7} | {'2sq(q)':>7} | "
           f"{'max|l_nt|':>9} | {'Ratio':>7} | {'Ramanujan':>10}")
    print(hdr)
    print("-" * 115)

    n_tested = 0
    n_ramanujan = 0
    max_dim_tested = 0

    for name, res in results.items():
        s = res['spectral']
        ram_str = 'YES' if s['is_ramanujan'] else '** NO **'
        reg_str = 'yes' if s['is_regular'] else 'no'
        n_tested += 1
        if s['is_ramanujan']:
            n_ramanujan += 1
        max_dim_tested = max(max_dim_tested, res['dim'])

        print(f"{name:>10} | {res['family']:>6} | {res['dim']:4d} | "
              f"{res['n_edges']:5d} | {reg_str:>7} | {s['degree']:5.1f} | "
              f"{s['lambda_max']:7.2f} | {s['ramanujan_bound']:7.2f} | "
              f"{s['max_nontrivial']:9.4f} | {s['ratio']:7.4f} | {ram_str:>10}")

    print("-" * 115)
    print(f"  TOTAL: {n_ramanujan}/{n_tested} Ramanujan   "
          f"(max dim tested: {max_dim_tested})")

    if n_ramanujan == n_tested:
        print("\n  *** ALL SIMPLE LIE ALGEBRAS TESTED SATISFY THE RAMANUJAN PROPERTY ***")
        if max_dim_tested >= 248:
            print("  *** INCLUDING E8 (dim=248) — CANNOT BE A SMALL-GRAPH COINCIDENCE ***")
    else:
        failures = [n for n, r in results.items() if not r['spectral']['is_ramanujan']]
        print(f"\n  FAILURES: {', '.join(failures)}")

    return n_tested, n_ramanujan


# ============================================================
# Visualization
# ============================================================

def plot_results(results, save=False):
    """Multi-panel visualization of extended results."""
    names = list(results.keys())
    dims = [results[n]['dim'] for n in names]
    ratios = [results[n]['spectral']['ratio'] for n in names]
    is_ram = [results[n]['spectral']['is_ramanujan'] for n in names]
    gaps = [results[n]['spectral']['spectral_gap'] for n in names]
    families_list = [results[n]['family'] for n in names]

    fig, axes = plt.subplots(2, 2, figsize=(18, 12))
    fig.suptitle('Extended Ramanujan Check: All Simple Lie Algebras\n'
                 '(Adjoint commutation graph from root system)',
                 fontsize=14, fontweight='bold')

    fam_colors = {'A': 'royalblue', 'B': 'crimson', 'C': 'forestgreen',
                  'D': 'darkorange', 'G': 'purple', 'F': 'brown', 'E': 'gold'}
    fam_markers = {'A': 'o', 'B': 's', 'C': '^', 'D': 'v',
                   'G': 'D', 'F': 'p', 'E': '*'}

    # Panel 1: Ratio vs dimension
    ax = axes[0, 0]
    for i, n in enumerate(names):
        f = families_list[i]
        ax.scatter(dims[i], ratios[i], c=fam_colors[f], marker=fam_markers[f],
                   s=120 if f in ('G', 'F', 'E') else 60, edgecolors='black',
                   linewidths=0.5, zorder=3)
    for i, n in enumerate(names):
        if dims[i] >= 50 or families_list[i] in ('G', 'F', 'E') or ratios[i] > 0.5:
            ax.annotate(n, (dims[i], ratios[i]), fontsize=7,
                        ha='left', va='bottom', xytext=(4, 4),
                        textcoords='offset points')
    ax.axhline(1.0, color='red', linestyle='--', alpha=0.5, label='Ramanujan bound')
    ax.set_xlabel('dim(algebra)')
    ax.set_ylabel('Ramanujan ratio  max|$\\lambda_{nt}$| / 2$\\sqrt{q}$')
    ax.set_title('Ramanujan Ratio vs Dimension')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0, top=max(1.2, max(ratios) * 1.1))

    # Panel 2: Ratio by family (connected)
    ax = axes[0, 1]
    fam_data = OrderedDict()
    for n, r in results.items():
        f = r['family']
        if f not in fam_data:
            fam_data[f] = ([], [])
        fam_data[f][0].append(r['dim'])
        fam_data[f][1].append(r['spectral']['ratio'])
    for f, (ds, rs) in fam_data.items():
        ax.plot(ds, rs, marker=fam_markers.get(f, 'o'), color=fam_colors.get(f, 'grey'),
                linestyle='-', alpha=0.8, label=f'Family {f}', markersize=8)
    ax.axhline(1.0, color='red', linestyle='--', alpha=0.3)
    ax.set_xlabel('dim(algebra)')
    ax.set_ylabel('Ramanujan ratio')
    ax.set_title('Scaling by Family')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)

    # Panel 3: Spectral gap vs dimension
    ax = axes[1, 0]
    for i, n in enumerate(names):
        f = families_list[i]
        ax.scatter(dims[i], gaps[i], c=fam_colors[f], marker=fam_markers[f],
                   s=120 if f in ('G', 'F', 'E') else 60, edgecolors='black',
                   linewidths=0.5)
    for i, n in enumerate(names):
        if dims[i] >= 50 or families_list[i] in ('G', 'F', 'E'):
            ax.annotate(n, (dims[i], gaps[i]), fontsize=7,
                        ha='left', va='bottom', xytext=(4, 4),
                        textcoords='offset points')
    ax.set_xlabel('dim(algebra)')
    ax.set_ylabel('Spectral gap ($\\lambda_1 - \\lambda_2$)')
    ax.set_title('Spectral Gap Growth')
    ax.grid(True, alpha=0.3)

    # Panel 4: Bar chart sorted by dimension
    ax = axes[1, 1]
    sorted_names = sorted(names, key=lambda n: results[n]['dim'])
    sorted_ratios = [results[n]['spectral']['ratio'] for n in sorted_names]
    sorted_colors = [fam_colors[results[n]['family']] for n in sorted_names]
    x_pos = range(len(sorted_names))
    ax.bar(x_pos, sorted_ratios, color=sorted_colors, alpha=0.8,
           edgecolor='black', linewidth=0.3)
    ax.axhline(1.0, color='red', linestyle='--', alpha=0.5, label='Ramanujan bound')
    ax.set_xticks(list(x_pos))
    ax.set_xticklabels(sorted_names, rotation=75, ha='right', fontsize=6)
    ax.set_ylabel('Ramanujan ratio')
    ax.set_title('All Algebras (sorted by dimension)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()

    if save:
        out_path = os.path.join(os.path.dirname(__file__), 'ramanujan_extended_results.png')
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        print(f"\n  Saved plot to {out_path}")
    else:
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
                        help='Quick mode: skip E7 and E8')
    parser.add_argument('--families', nargs='+', default=None,
                        help='Families to test: A B C D G F E')
    args = parser.parse_args()

    results = run_analysis(families=args.families, quick=args.quick)
    n_tested, n_ramanujan = print_summary(results)
    cross_validate_small(results)
    plot_results(results, save=args.save)
