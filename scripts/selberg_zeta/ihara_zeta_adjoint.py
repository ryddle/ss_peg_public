"""
Priority 3: Ihara Zeta Function of Adjoint Commutation Graphs
==============================================================

Computes the Ihara zeta function Z(u) for the adjoint commutation graph
of ALL simple Lie algebras, and verifies:

1. Sunada's equivalence: Ramanujan ⟺ Riemann Hypothesis for Z(u)
   (all nontrivial poles lie on |u| = 1/√q for regular graphs;
    generalised to the "Ramanujan circle" for irregular graphs)

2. Pole distribution: how the pole structure changes across families

3. N→∞ limit: how the pole density evolves with growing rank

THEORY:
  For a general (possibly irregular) graph G with adjacency matrix A and
  degree matrix D = diag(deg(v)):

    Z(u)^{-1} = det(I - uA + u²(D - I))

  (Bass's formula — no regularity needed.)

  The poles of Z(u) are the values u where Z(u)^{-1} = 0.
  This reduces to a polynomial eigenvalue problem:

    det(u²Q - uA + I) = 0,  Q = D - I

  which we linearise as a 2n×2n generalised eigenvalue problem.

  For REGULAR graphs (degree q+1):
    Q = qI, and the poles satisfy: 1 - uλ + qu² = 0
    where λ are eigenvalues of A.
    → u = (λ ± √(λ²-4q)) / (2q)

  RIEMANN HYPOTHESIS FOR IHARA ZETA:
    The nontrivial poles of Z(u) satisfy |u| = q^{-1/2}
    ⟺ the graph is Ramanujan.

  For IRREGULAR graphs, the generalised RH involves the spectral radius
  of the universal cover (Hashimoto operator). We use the empirical
  definition: all "bulk" poles lie on/inside a critical circle.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from ramanujan_extended import (
    roots_A, roots_B, roots_C, roots_D,
    roots_G2, roots_F4, roots_E6, roots_E7, roots_E8,
    adjacency_from_roots, spectral_analysis,
)


# ============================================================
# 1. Ihara zeta poles via Bass's determinant formula
# ============================================================

def ihara_poles_bass(adjacency):
    """
    Compute poles of the Ihara zeta Z(u) for a general graph via Bass's formula.

    Z(u)^{-1} = det(I - uA + u²Q)  where Q = D - I, D = degree matrix.

    Poles are solutions of: det(u²Q - uA + I) = 0
    Linearise:  M·v = u·N·v  where
      M = [[I, -A], [0, I]]
      N = [[0, -Q], [I, 0]]

    Returns (poles, pole_info) where pole_info has classification details.
    """
    n = adjacency.shape[0]
    D = np.diag(adjacency.sum(axis=1))
    Q = D - np.eye(n)
    I_n = np.eye(n)
    Z = np.zeros((n, n))

    M = np.block([[I_n, -adjacency], [Z, I_n]])
    N = np.block([[Z, -Q], [I_n, Z]])

    det_N = np.linalg.det(N)
    if abs(det_N) > 1e-12:
        # N is invertible: standard eigenvalue problem
        poles = np.linalg.eigvals(np.linalg.solve(N, M))
    else:
        # N singular: use generalised eigenvalue via scipy
        try:
            from scipy.linalg import eigvals as sp_eigvals
            poles = sp_eigvals(M, N)
            # Filter out infinities
            poles = poles[np.isfinite(poles)]
        except ImportError:
            # Fallback: eigenvalue approach using adjacency spectrum
            poles = _ihara_poles_from_spectrum(adjacency)

    return poles


def _ihara_poles_from_spectrum(adjacency):
    """Fallback: compute Ihara poles using adjacency eigenvalues (mean-field Q)."""
    degrees = adjacency.sum(axis=1)
    q_mean = np.mean(degrees) - 1
    eigenvalues = np.linalg.eigvals(adjacency)
    poles = []
    for lam in eigenvalues:
        disc = lam**2 - 4 * q_mean + 0j
        sqrt_d = np.sqrt(disc)
        if abs(q_mean) > 1e-10:
            poles.extend([(lam + sqrt_d) / (2 * q_mean),
                          (lam - sqrt_d) / (2 * q_mean)])
        elif abs(lam) > 1e-10:
            poles.append(1.0 / lam)
    return np.array(poles, dtype=complex)


def ihara_poles_hashimoto(adjacency):
    """
    Compute Ihara poles via the Hashimoto (edge adjacency) operator.

    The Hashimoto operator T_e on directed edges (i→j) is defined by:
      T_e[(i→j), (k→l)] = 1  iff j=k AND l≠i (non-backtracking step)

    Then: Z(u)^{-1} = det(I - u·T_e)
    and the poles are u = 1/μ where μ are eigenvalues of T_e.

    This is the EXACT formulation — valid for any graph, no regularity needed.
    """
    n = adjacency.shape[0]

    # Build directed edge list
    edges = []
    for i in range(n):
        for j in range(n):
            if adjacency[i, j] > 0.5:
                edges.append((i, j))

    m = len(edges)  # number of directed edges = 2|E| for undirected
    if m == 0:
        return np.array([]), edges

    # Build Hashimoto operator
    T_e = np.zeros((m, m), dtype=np.float64)
    edge_index = {e: idx for idx, e in enumerate(edges)}

    for idx, (i, j) in enumerate(edges):
        # Find all edges (j, k) where k != i
        for k in range(n):
            if k != i and adjacency[j, k] > 0.5:
                jk_idx = edge_index.get((j, k))
                if jk_idx is not None:
                    T_e[jk_idx, idx] = 1.0

    # Eigenvalues of T_e
    eigenvalues = np.linalg.eigvals(T_e)

    # Poles of Z(u) = 1/μ for nonzero μ
    nonzero = eigenvalues[np.abs(eigenvalues) > 1e-10]
    poles = 1.0 / nonzero

    return poles, eigenvalues


# ============================================================
# 2. Pole classification and Riemann Hypothesis check
# ============================================================

def classify_poles(poles, adjacency):
    """
    Classify Ihara poles into trivial and nontrivial, and verify the RH.

    For a graph with mean degree q+1:
      - Trivial poles: u = ±1/q (from λ = ±(q+1)) and u = ±1 (from Betti factor)
      - Nontrivial poles: all others
      - RH circle radius: 1/√q

    Ramanujan ⟺ all nontrivial poles satisfy |u| ≥ 1/√q
    """
    degrees = adjacency.sum(axis=1)
    d_mean = np.mean(degrees)
    q_mean = d_mean - 1
    d_max = np.max(degrees)
    d_min = np.min(degrees)

    # Reference circles
    ram_radius = 1.0 / np.sqrt(q_mean) if q_mean > 0 else float('inf')
    trivial_radius = 1.0 / q_mean if q_mean > 0 else float('inf')

    # Spectral info
    spec = spectral_analysis(adjacency)

    # Classify poles
    trivial_poles = []
    nontrivial_poles = []

    for p in poles:
        if not np.isfinite(p):
            continue
        ap = abs(p)
        # Trivial: near ±1/d_mean, ±1/(d_mean-2), or ±1
        if abs(ap - 1.0) < 0.01 or abs(ap - trivial_radius) < 0.01:
            trivial_poles.append(p)
        else:
            nontrivial_poles.append(p)

    nontrivial_poles = np.array(nontrivial_poles) if nontrivial_poles else np.array([], dtype=complex)
    trivial_poles = np.array(trivial_poles) if trivial_poles else np.array([], dtype=complex)

    # RH check: all nontrivial poles on or outside the Ramanujan circle
    if len(nontrivial_poles) > 0:
        min_nt_radius = np.min(np.abs(nontrivial_poles))
        max_nt_radius = np.max(np.abs(nontrivial_poles))
        rh_holds = min_nt_radius >= ram_radius - 1e-6
    else:
        min_nt_radius = float('inf')
        max_nt_radius = 0
        rh_holds = True

    return {
        'poles': poles,
        'trivial': trivial_poles,
        'nontrivial': nontrivial_poles,
        'n_trivial': len(trivial_poles),
        'n_nontrivial': len(nontrivial_poles),
        'ram_radius': ram_radius,
        'trivial_radius': trivial_radius,
        'min_nt_radius': min_nt_radius,
        'max_nt_radius': max_nt_radius,
        'rh_holds': rh_holds,
        'is_ramanujan': spec['is_ramanujan'],
        'd_mean': d_mean,
        'q_mean': q_mean,
        'd_max': d_max,
        'd_min': d_min,
        'spectral': spec,
    }


# ============================================================
# 3. Comprehensive analysis for a single algebra
# ============================================================

def analyze_algebra(name, rank, all_roots, simple):
    """Full Ihara analysis for one algebra."""
    adj = adjacency_from_roots(rank, all_roots, simple)
    dim = adj.shape[0]

    # Method 1: Bass determinant formula
    poles_bass = ihara_poles_bass(adj)

    # Method 2: Hashimoto (only for small graphs, memory: O(4|E|²))
    n_edges = int(np.sum(adj) / 2)
    n_directed = 2 * n_edges
    use_hashimoto = n_directed <= 4000  # Memory limit

    poles_hash = None
    hash_eigenvalues = None
    if use_hashimoto:
        poles_hash, hash_eigenvalues = ihara_poles_hashimoto(adj)

    # Classification
    info = classify_poles(poles_bass, adj)

    # Betti number
    betti = n_edges - dim + 1

    info.update({
        'name': name,
        'dim': dim,
        'n_edges': n_edges,
        'betti': betti,
        'poles_bass': poles_bass,
        'poles_hashimoto': poles_hash,
        'hash_eigenvalues': hash_eigenvalues,
        'adjacency': adj,
    })

    return info


# ============================================================
# 4. Build catalog and run analysis
# ============================================================

def build_ihara_catalog():
    """Build catalog of all algebras for Ihara analysis."""
    catalog = []

    # A_n = su(n+1)
    for n in range(1, 12):
        r, roots, simple = roots_A(n)
        name = f"su({n+1})"
        catalog.append((name, r, roots, simple))

    # B_n = so(2n+1)
    for n in range(2, 9):
        r, roots, simple = roots_B(n)
        name = f"so({2*n+1})"
        catalog.append((name, r, roots, simple))

    # C_n = sp(2n)
    for n in range(2, 9):
        r, roots, simple = roots_C(n)
        name = f"sp({2*n})"
        catalog.append((name, r, roots, simple))

    # D_n = so(2n)
    for n in range(3, 10):
        r, roots, simple = roots_D(n)
        name = f"so({2*n})"
        catalog.append((name, r, roots, simple))

    # Exceptionals
    for exc_name, fn, rnk in [('G2', roots_G2, 2), ('F4', roots_F4, 4),
                                ('E6', roots_E6, 6), ('E7', roots_E7, 7),
                                ('E8', roots_E8, 8)]:
        r, roots, simple = fn()
        catalog.append((exc_name, r, roots, simple))

    return catalog


def run_full_analysis():
    """Run Ihara zeta analysis for all algebras."""
    print("=" * 100)
    print("  PRIORITY 3: Ihara Zeta Function of Adjoint Commutation Graphs")
    print("=" * 100)
    print()

    catalog = build_ihara_catalog()
    results = []

    # Header
    hdr = (f"  {'Algebra':>10s}  {'dim':>4s}  {'|E|':>5s}  {'β₁':>4s}  "
           f"{'#poles':>6s}  {'#nt':>4s}  "
           f"{'1/√q':>7s}  {'min|u_nt|':>9s}  {'max|u_nt|':>9s}  "
           f"{'RH?':>4s}  {'Ram?':>4s}  {'≡':>2s}")
    print(hdr)
    print("  " + "-" * 96)

    for name, rank, all_roots, simple in catalog:
        info = analyze_algebra(name, rank, all_roots, simple)
        results.append(info)

        equiv = "✓" if info['rh_holds'] == info['is_ramanujan'] else "✗"
        print(f"  {name:>10s}  {info['dim']:4d}  {info['n_edges']:5d}  "
              f"{info['betti']:4d}  {len(info['poles_bass']):6d}  "
              f"{info['n_nontrivial']:4d}  "
              f"{info['ram_radius']:7.4f}  "
              f"{info['min_nt_radius']:9.4f}  "
              f"{info['max_nt_radius']:9.4f}  "
              f"{'YES' if info['rh_holds'] else 'NO':>4s}  "
              f"{'YES' if info['is_ramanujan'] else 'NO':>4s}  {equiv}")

    # Summary
    print()
    print("=" * 100)
    print("  SUNADA EQUIVALENCE VERIFICATION")
    print("=" * 100)
    n_total = len(results)
    n_equiv = sum(1 for r in results if r['rh_holds'] == r['is_ramanujan'])
    n_ram = sum(1 for r in results if r['is_ramanujan'])
    n_rh = sum(1 for r in results if r['rh_holds'])
    print(f"  Total algebras: {n_total}")
    print(f"  Ramanujan: {n_ram}/{n_total}")
    print(f"  RH holds: {n_rh}/{n_total}")
    print(f"  Equivalence (RH ⟺ Ram): {n_equiv}/{n_total}")
    if n_equiv == n_total:
        print("  → PERFECT: Sunada equivalence holds for ALL tested algebras")
    else:
        print("  → DISCREPANCIES found:")
        for r in results:
            if r['rh_holds'] != r['is_ramanujan']:
                print(f"    {r['name']}: RH={r['rh_holds']}, Ram={r['is_ramanujan']}")

    # Hashimoto cross-validation
    print()
    print("=" * 100)
    print("  HASHIMOTO CROSS-VALIDATION (where feasible)")
    print("=" * 100)
    for r in results:
        if r['poles_hashimoto'] is not None:
            p_bass = np.sort(np.abs(r['poles_bass']))
            p_hash = np.sort(np.abs(r['poles_hashimoto']))
            # Compare pole radii distributions
            if len(p_bass) > 0 and len(p_hash) > 0:
                # Find matching poles by radius
                min_bass = np.min(p_bass[p_bass > 1e-6]) if np.any(p_bass > 1e-6) else float('inf')
                min_hash = np.min(p_hash[p_hash > 1e-6]) if np.any(p_hash > 1e-6) else float('inf')
                print(f"  {r['name']:>10s}: Bass min|u|={min_bass:.6f}  "
                      f"Hash min|u|={min_hash:.6f}  "
                      f"#Bass={len(p_bass)}  #Hash={len(p_hash)}  "
                      f"{'MATCH' if abs(min_bass-min_hash)<0.01 else 'DIFFER'}")

    # Pole density analysis
    print()
    print("=" * 100)
    print("  POLE DISTRIBUTION ANALYSIS")
    print("=" * 100)
    print()
    print(f"  {'Algebra':>10s}  {'dim':>4s}  {'#poles':>6s}  {'#zero':>5s}  "
          f"{'1/√q':>7s}  {'med|u_nt|':>9s}  {'σ(|u|)':>8s}  "
          f"{'frac_on_circle':>14s}")
    print("  " + "-" * 80)

    for r in results:
        nt = r['nontrivial']
        if len(nt) == 0:
            continue
        radii = np.abs(nt)
        on_circle = np.sum(np.abs(radii - r['ram_radius']) < 0.01 * r['ram_radius'])
        frac = on_circle / len(nt) if len(nt) > 0 else 0
        n_zero = np.sum(np.abs(r['poles_bass']) < 1e-6)
        print(f"  {r['name']:>10s}  {r['dim']:4d}  {len(r['poles_bass']):6d}  "
              f"{n_zero:5d}  "
              f"{r['ram_radius']:7.4f}  {np.median(radii):9.4f}  "
              f"{np.std(radii):8.4f}  {frac:14.3f}")

    # N→∞ scaling of pole density
    print()
    print("=" * 100)
    print("  N→∞ SCALING: Pole distribution vs rank (A_n family)")
    print("=" * 100)
    print()
    a_results = [r for r in results if r['name'].startswith('su(')]
    for r in a_results:
        nt = r['nontrivial']
        if len(nt) == 0:
            print(f"  {r['name']:>10s}: no nontrivial poles")
            continue
        radii = np.abs(nt)
        radii_real = np.sort(radii)

        # Pole ring: what fraction lies in [1/√q - ε, 1/√q + ε]?
        eps = 0.05
        rc = r['ram_radius']
        in_ring = np.sum((radii >= rc - eps) & (radii <= rc + eps))
        frac_ring = in_ring / len(radii)

        # Innermost nontrivial pole
        inner = np.min(radii) if len(radii) > 0 else float('inf')

        print(f"  {r['name']:>10s}: dim={r['dim']:4d}  1/√q={rc:.4f}  "
              f"inner={inner:.4f}  inner/rc={inner/rc:.4f}  "
              f"frac_in_ring={frac_ring:.3f}  #nt={len(nt)}")

    # ========================================================
    # GENERALIZED RH via Hashimoto spectral radius
    # ========================================================
    print()
    print("=" * 100)
    print("  GENERALIZED RH: Hashimoto spectral radius ρ(T_e)")
    print("=" * 100)
    print()
    print("  For irregular graphs, the correct RH uses ρ = ρ(T_e), the spectral")
    print("  radius of the Hashimoto (non-backtracking) operator:")
    print("    Generalized Ramanujan ⟺ all nontrivial eigenvalues |μ| ≤ √ρ")
    print("    Equivalently: |u_pole| ≥ 1/√ρ for nontrivial poles")
    print()
    print(f"  {'Algebra':>10s}  {'dim':>4s}  {'ρ(T_e)':>8s}  {'1/√ρ':>7s}  "
          f"{'1/√q_mean':>9s}  {'min|u_nt|':>9s}  "
          f"{'gRH?':>5s}  {'Ram?':>4s}  {'≡':>2s}")
    print("  " + "-" * 80)

    n_gen_equiv = 0
    hash_results = {}
    for r in results:
        if r['hash_eigenvalues'] is not None:
            he = r['hash_eigenvalues']
            rho = np.max(np.abs(he))
            gen_ram_radius = 1.0 / np.sqrt(rho) if rho > 0 else float('inf')

            # Nontrivial: eigenvalues with |μ| < ρ - tiny tolerance
            nontrivial_eigs = he[np.abs(he) < rho - 1e-6]
            if len(nontrivial_eigs) > 0:
                max_nt_eig = np.max(np.abs(nontrivial_eigs))
                gen_rh = max_nt_eig <= np.sqrt(rho) + 1e-6
            else:
                max_nt_eig = 0
                gen_rh = True

            # Also check via poles
            nt_poles = r['nontrivial']
            if len(nt_poles) > 0:
                min_nt = np.min(np.abs(nt_poles))
                gen_rh_poles = min_nt >= gen_ram_radius - 1e-6
            else:
                min_nt = float('inf')
                gen_rh_poles = True

            equiv = "✓" if gen_rh == r['is_ramanujan'] else "✗"
            if gen_rh == r['is_ramanujan']:
                n_gen_equiv += 1

            hash_results[r['name']] = {
                'rho': rho, 'gen_ram_radius': gen_ram_radius,
                'gen_rh': gen_rh, 'max_nt_eig': max_nt_eig,
            }

            print(f"  {r['name']:>10s}  {r['dim']:4d}  {rho:8.3f}  "
                  f"{gen_ram_radius:7.4f}  {r['ram_radius']:9.4f}  "
                  f"{r['min_nt_radius']:9.4f}  "
                  f"{'YES' if gen_rh else 'NO':>5s}  "
                  f"{'YES' if r['is_ramanujan'] else 'NO':>4s}  {equiv}")

    n_hash_checked = sum(1 for r in results if r['hash_eigenvalues'] is not None)
    print()
    print(f"  Generalized RH equivalence: {n_gen_equiv}/{n_hash_checked} "
          f"(Hashimoto available for {n_hash_checked}/{len(results)} algebras)")

    # Store hash_results for plotting
    for r in results:
        r['hash_results'] = hash_results.get(r['name'])

    return results


# ============================================================
# 5. Visualization
# ============================================================

def plot_ihara_results(results):
    """Generate comprehensive visualization of Ihara zeta analysis."""
    fig = plt.figure(figsize=(24, 18))
    fig.suptitle("Ihara Zeta Function — Adjoint Commutation Graphs\n"
                 "Pole Distribution and Riemann Hypothesis Verification",
                 fontsize=16, fontweight='bold')

    gs = GridSpec(3, 4, hspace=0.35, wspace=0.30)

    # Panel 1-4: Pole plots for selected algebras
    selected = ['su(3)', 'su(7)', 'su(8)', 'E8']
    for idx, name in enumerate(selected):
        ax = fig.add_subplot(gs[0, idx])
        r = next((x for x in results if x['name'] == name), None)
        if r is None:
            continue
        poles = r['poles_bass']
        nt = r['nontrivial']
        triv = r['trivial']
        rc = r['ram_radius']

        # Ramanujan circle
        theta = np.linspace(0, 2*np.pi, 200)
        ax.plot(rc * np.cos(theta), rc * np.sin(theta), 'r-', alpha=0.5,
                linewidth=2, label=f'1/√q={rc:.3f}')
        # Unit circle
        ax.plot(np.cos(theta), np.sin(theta), 'k--', alpha=0.2, linewidth=1)

        # Plot poles
        if len(triv) > 0:
            ax.plot(triv.real, triv.imag, 'ks', markersize=4, alpha=0.3, label='trivial')
        if len(nt) > 0:
            ax.plot(nt.real, nt.imag, 'b.', markersize=3, alpha=0.5, label='nontrivial')

        ax.set_title(f"{name} (dim={r['dim']})\nRam={'✓' if r['is_ramanujan'] else '✗'}, "
                     f"RH={'✓' if r['rh_holds'] else '✗'}", fontsize=10)
        ax.set_aspect('equal')
        ax.legend(fontsize=7, loc='upper right')
        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-1.2, 1.2)
        ax.grid(True, alpha=0.2)

    # Panel 5: Ramanujan radius vs inner pole radius
    ax5 = fig.add_subplot(gs[1, 0])
    for r in results:
        color = 'green' if r['is_ramanujan'] else 'red'
        marker = 'o' if r['name'].startswith('su') else ('^' if r['name'].startswith('so') else 's')
        ax5.plot(r['ram_radius'], r['min_nt_radius'], marker, color=color,
                 markersize=5, alpha=0.7)
    lim = ax5.get_xlim()
    ax5.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='1/√q = min|u_nt|')
    ax5.set_xlabel('1/√q (Ramanujan radius)')
    ax5.set_ylabel('min|u_nt| (inner pole)')
    ax5.set_title('RH check: inner pole vs Ramanujan circle')
    ax5.legend(fontsize=8)
    ax5.grid(True, alpha=0.2)

    # Panel 6: Sunada equivalence summary
    ax6 = fig.add_subplot(gs[1, 1])
    names = [r['name'] for r in results]
    ram_vals = [1 if r['is_ramanujan'] else 0 for r in results]
    rh_vals = [1 if r['rh_holds'] else 0 for r in results]
    x = range(len(names))
    ax6.bar([i - 0.15 for i in x], ram_vals, 0.3, label='Ramanujan', color='blue', alpha=0.6)
    ax6.bar([i + 0.15 for i in x], rh_vals, 0.3, label='RH holds', color='orange', alpha=0.6)
    ax6.set_xticks(list(x)[::3])
    ax6.set_xticklabels([names[i] for i in range(0, len(names), 3)], rotation=45, fontsize=7)
    ax6.set_title('Sunada: Ramanujan ≡ RH for Ihara')
    ax6.legend(fontsize=8)
    ax6.set_ylabel('Yes/No')

    # Panel 7: Pole count scaling
    ax7 = fig.add_subplot(gs[1, 2])
    dims = [r['dim'] for r in results]
    n_poles = [len(r['poles_bass']) for r in results]
    n_nt = [r['n_nontrivial'] for r in results]
    ax7.plot(dims, n_poles, 'bo-', markersize=4, label='Total poles', alpha=0.6)
    ax7.plot(dims, n_nt, 'r^-', markersize=4, label='Nontrivial', alpha=0.6)
    ax7.set_xlabel('dim(algebra)')
    ax7.set_ylabel('Number of poles')
    ax7.set_title('Pole count scaling')
    ax7.legend(fontsize=8)
    ax7.grid(True, alpha=0.2)

    # Panel 8: Inner/outer pole ratio vs dimension
    ax8 = fig.add_subplot(gs[1, 3])
    for r in results:
        if r['min_nt_radius'] < 100:
            ratio_io = r['min_nt_radius'] / r['ram_radius'] if r['ram_radius'] > 0 else 0
            color = 'green' if r['is_ramanujan'] else 'red'
            ax8.plot(r['dim'], ratio_io, 'o', color=color, markersize=5, alpha=0.7)
    ax8.axhline(y=1.0, color='k', linestyle='--', alpha=0.3, label='RH boundary')
    ax8.set_xlabel('dim(algebra)')
    ax8.set_ylabel('min|u_nt| / (1/√q)')
    ax8.set_title('RH ratio vs dimension')
    ax8.legend(fontsize=8)
    ax8.grid(True, alpha=0.2)

    # Panels 9-12: Pole radius histograms for A_n family
    a_results = [r for r in results if r['name'].startswith('su(')]
    selected_a = [r for r in a_results if r['name'] in ['su(3)', 'su(6)', 'su(9)', 'su(12)']]
    for idx, r in enumerate(selected_a[:4]):
        ax = fig.add_subplot(gs[2, idx])
        nt = r['nontrivial']
        if len(nt) > 0:
            radii = np.abs(nt)
            ax.hist(radii, bins=40, density=True, alpha=0.7, color='steelblue',
                    edgecolor='navy', linewidth=0.5)
            ax.axvline(x=r['ram_radius'], color='red', linewidth=2,
                       label=f'1/√q={r["ram_radius"]:.3f}')
        ax.set_title(f"{r['name']} pole radius distribution", fontsize=10)
        ax.set_xlabel('|u|')
        ax.set_ylabel('density')
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.2)

    outpath = os.path.join(os.path.dirname(__file__), 'ihara_zeta_adjoint.png')
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\n  Saved visualization to {outpath}")
    plt.close(fig)


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    results = run_full_analysis()
    plot_ihara_results(results)
