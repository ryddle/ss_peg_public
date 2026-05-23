"""
Priority 4: Connection to West et al. Quantum Chaos Framework
==============================================================

West et al. (SciPost Phys. Core 8, 081 (2025), arXiv:2502.16404)
study the SAME graph object — the commutator graph of a DLA — and
link its graph-theoretic properties to:
  - OTOC (out-of-time-order correlator)
  - Frame potential
  - Krylov complexity
  - Graph complexity

This script computes the relevant graph-theoretic and spectral
quantities for all 37 simple Lie algebra adjoint graphs, then
correlates them with the Ramanujan property (P1-P2 results).

KEY CONNECTIONS:
  1. Average OTOC = 1 - 2·f_ac where f_ac = anticommutation fraction
     (West, Corollary 3)
  2. Average graph complexity at long times = average shortest path
     length from initial operator (West, Eq. 38)
  3. Krylov complexity ≥ graph complexity (West, Lemma 9)
  4. Spectral gap → mixing time → scrambling rate
  5. Ramanujan ⟹ optimal expander ⟹ fastest possible information
     spreading

THEORETICAL PREDICTIONS:
  - Ramanujan algebras should have smaller diameter/dim ratio
  - Ramanujan algebras should mix faster (lower mixing time bound)
  - The Ramanujan transition at h* ≈ 9.5 should be visible in
    ALL West et al. quantities
  - Standard Model algebras (h ≤ 6) → "optimally chaotic"

Usage:
    python west_connection.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from collections import OrderedDict
import os
import sys
import time

# Import root system and spectral tools from ramanujan_extended
sys.path.insert(0, os.path.dirname(__file__))
from ramanujan_extended import (
    roots_A, roots_B, roots_C, roots_D,
    roots_G2, roots_F4, roots_E6, roots_E7, roots_E8,
    adjacency_from_roots, spectral_analysis
)


# ============================================================
# 1. Graph-theoretic quantities (West et al. motivated)
# ============================================================

def shortest_path_matrix(adj):
    """
    All-pairs shortest paths via BFS on adjacency matrix.
    Returns distance matrix D where D[i,j] = shortest path length.
    Disconnected pairs get D[i,j] = -1.
    """
    n = adj.shape[0]
    dist = np.full((n, n), -1, dtype=np.int32)
    for source in range(n):
        # BFS from source
        visited = np.zeros(n, dtype=bool)
        visited[source] = True
        dist[source, source] = 0
        queue = [source]
        while queue:
            next_queue = []
            for u in queue:
                neighbors = np.where(adj[u] > 0.5)[0]
                for v in neighbors:
                    if not visited[v]:
                        visited[v] = True
                        dist[source, v] = dist[source, u] + 1
                        next_queue.append(v)
            queue = next_queue
    return dist


def graph_theoretic_analysis(adj, spec):
    """
    Compute all graph-theoretic quantities relevant to West et al.

    Returns dict with:
      - diameter: max shortest path
      - avg_shortest_path: mean over all pairs
      - avg_shortest_from_cartan: avg shortest path from Cartan generators
      - avg_shortest_from_roots: avg shortest path from root generators
      - edge_density: 2|E| / (n(n-1))
      - anticomm_fraction: fraction of pairs connected = edge_density
      - avg_otoc: average OTOC = 1 - 2·f_ac (West, Corollary 3)
      - mixing_time_bound: log(n-1) / log(λ_max / |λ_2|)
      - cheeger_lower: spectral_gap / (2 * d_max)
      - cheeger_upper: sqrt(2 * d_max * spectral_gap)
    """
    n = adj.shape[0]
    dist = shortest_path_matrix(adj)

    # Diameter
    diameter = np.max(dist[dist >= 0])

    # Average shortest path (excluding self-loops)
    mask = dist > 0
    avg_sp = np.mean(dist[mask]) if np.any(mask) else 0

    # Degree statistics
    degrees = adj.sum(axis=1)
    d_max = int(np.max(degrees))
    d_min = int(np.min(degrees))
    d_mean = np.mean(degrees)

    # Edge density and anticommutation fraction
    n_edges = int(np.sum(adj) / 2)
    n_pairs = n * (n - 1) / 2
    edge_density = n_edges / n_pairs if n_pairs > 0 else 0

    # Average OTOC (West, Corollary 3)
    # For simple Lie algebras, the adjoint graph is one connected component
    # (plus identity which we don't include).
    # avg_OTOC = 1 - 2 * f_anticommute
    # In root basis: f_anticommute = edge_density (fraction of pairs that don't commute)
    avg_otoc = 1 - 2 * edge_density

    # Mixing time bound from spectral gap
    lambda_max = spec['lambda_max']
    lambda_2 = spec['lambda_2']
    spectral_gap = spec['spectral_gap']

    if abs(lambda_2) > 1e-10 and lambda_max > abs(lambda_2):
        mixing_time = np.log(n - 1) / np.log(lambda_max / abs(lambda_2))
    else:
        mixing_time = float('inf')

    # Alternative mixing time: t_mix ~ 1/spectral_gap * log(n)
    # For lazy random walk on the graph: t_mix ≤ (d_max / spectral_gap) * log(n)
    mixing_time_rw = (d_max / spectral_gap * np.log(n)) if spectral_gap > 1e-10 else float('inf')

    # Cheeger bounds (from Cheeger inequality)
    # spectral_gap / (2 * d_max) ≤ h(G) ≤ sqrt(2 * d_max * spectral_gap)
    cheeger_lower = spectral_gap / (2 * d_max) if d_max > 0 else 0
    cheeger_upper = np.sqrt(2 * d_max * spectral_gap) if spectral_gap > 0 else 0

    # Average shortest path per node class
    rank = int(np.sum(degrees < d_mean - 0.5))  # Approx: Cartan elements have lower degree
    # More robust: Cartan generators are the first 'rank' nodes in our encoding
    # We know the structure: first rank nodes are H_i, rest are E_α

    # Average graph complexity at long times (West, Eq. 38)
    # For each node p: avg_G = (1/|C|) * sum_{q in C} ell(p,q)
    avg_graph_complexity = np.zeros(n)
    for i in range(n):
        component = dist[i] >= 0
        component[i] = False  # Exclude self
        if np.any(component):
            avg_graph_complexity[i] = np.mean(dist[i, component])

    mean_graph_complexity = np.mean(avg_graph_complexity)

    return {
        'diameter': diameter,
        'avg_shortest_path': avg_sp,
        'mean_graph_complexity': mean_graph_complexity,
        'edge_density': edge_density,
        'n_edges': n_edges,
        'anticomm_fraction': edge_density,
        'avg_otoc': avg_otoc,
        'mixing_time_spectral': mixing_time,
        'mixing_time_rw': mixing_time_rw,
        'cheeger_lower': cheeger_lower,
        'cheeger_upper': cheeger_upper,
        'd_max': d_max,
        'd_min': d_min,
        'd_mean': d_mean,
        'dist_matrix': dist,
        'avg_graph_complexity_per_node': avg_graph_complexity,
    }


# ============================================================
# 2. Coxeter numbers and algebra metadata
# ============================================================

COXETER_NUMBERS = {
    'su(2)': 2, 'su(3)': 3, 'su(4)': 4, 'su(5)': 5,
    'su(6)': 6, 'su(7)': 7, 'su(8)': 8, 'su(9)': 9,
    'su(10)': 10, 'su(11)': 11, 'su(12)': 12,
    'so(5)': 4, 'so(7)': 6, 'so(9)': 8, 'so(10)': 8,
    'so(11)': 10, 'so(12)': 10, 'so(13)': 12, 'so(14)': 12,
    'so(15)': 14, 'so(16)': 14, 'so(17)': 16, 'so(18)': 16,
    'sp(4)': 4, 'sp(6)': 6, 'sp(8)': 8, 'sp(10)': 10,
    'sp(12)': 12, 'sp(14)': 14, 'sp(16)': 16,
    'G2': 6, 'F4': 12, 'E6': 12, 'E7': 18, 'E8': 30
}


def get_family(name):
    if name.startswith('su'): return 'A'
    if name.startswith('so'): return 'B/D'
    if name.startswith('sp'): return 'C'
    return 'Exc'


# ============================================================
# 3. Build catalog
# ============================================================

def build_west_catalog():
    """Build catalog of all algebras for West connection analysis."""
    catalog = []

    # A_n: su(n+1) for n=1..11
    for n in range(1, 12):
        r, roots, simple = roots_A(n)
        catalog.append((f"su({n+1})", r, roots, simple))

    # B_n: so(2n+1) for n=2..8
    for n in range(2, 9):
        r, roots, simple = roots_B(n)
        catalog.append((f"so({2*n+1})", r, roots, simple))

    # C_n: sp(2n) for n=2..8
    for n in range(2, 9):
        r, roots, simple = roots_C(n)
        catalog.append((f"sp({2*n})", r, roots, simple))

    # D_n: so(2n) for n=3..9
    for n in range(3, 10):
        r, roots, simple = roots_D(n)
        catalog.append((f"so({2*n})", r, roots, simple))

    # Exceptionals
    for exc_name, fn in [('G2', roots_G2), ('F4', roots_F4),
                          ('E6', roots_E6), ('E7', roots_E7),
                          ('E8', roots_E8)]:
        r, roots, simple = fn()
        catalog.append((exc_name, r, roots, simple))

    return catalog


# ============================================================
# 4. Main analysis
# ============================================================

def run_full_analysis():
    """Run West connection analysis for all algebras."""
    print("=" * 110)
    print("  PRIORITY 4: Connection to West et al. Quantum Chaos Framework")
    print("  (SciPost Phys. Core 8, 081 (2025) — arXiv:2502.16404)")
    print("=" * 110)
    print()

    catalog = build_west_catalog()
    results = []

    # Header
    hdr = (f"  {'Algebra':>10s}  {'dim':>4s}  {'h':>3s}  "
           f"{'diam':>4s}  {'⟨ℓ⟩':>6s}  {'⟨G⟩':>6s}  "
           f"{'f_ac':>6s}  {'⟨OTOC⟩':>7s}  "
           f"{'λ_gap':>6s}  {'t_mix':>6s}  "
           f"{'h_low':>6s}  {'Ram':>4s}")
    print(hdr)
    print("  " + "-" * 106)

    for name, rank, all_roots, simple in catalog:
        t0 = time.time()

        # Build adjacency and get spectral properties
        adj = adjacency_from_roots(rank, all_roots, simple)
        dim = adj.shape[0]
        spec = spectral_analysis(adj)

        # Graph-theoretic analysis
        gt = graph_theoretic_analysis(adj, spec)

        h_cox = COXETER_NUMBERS.get(name, -1)

        result = {
            'name': name,
            'dim': dim,
            'rank': rank,
            'h_coxeter': h_cox,
            'family': get_family(name),
            'spectral': spec,
            'graph': gt,
            'time': time.time() - t0,
        }
        results.append(result)

        ram_str = "YES" if spec['is_ramanujan'] else "NO"
        t_mix = gt['mixing_time_rw']
        t_mix_str = f"{t_mix:6.1f}" if t_mix < 1e6 else "   inf"

        print(f"  {name:>10s}  {dim:4d}  {h_cox:3d}  "
              f"{gt['diameter']:4d}  {gt['avg_shortest_path']:6.3f}  "
              f"{gt['mean_graph_complexity']:6.3f}  "
              f"{gt['anticomm_fraction']:6.4f}  {gt['avg_otoc']:7.4f}  "
              f"{spec['spectral_gap']:6.3f}  {t_mix_str}  "
              f"{gt['cheeger_lower']:6.4f}  {ram_str:>4s}")

    # ================================================================
    # Analysis Section 1: Ramanujan vs Non-Ramanujan comparison
    # ================================================================
    print()
    print("=" * 110)
    print("  COMPARISON: Ramanujan vs Non-Ramanujan algebras")
    print("=" * 110)

    ram = [r for r in results if r['spectral']['is_ramanujan']]
    non_ram = [r for r in results if not r['spectral']['is_ramanujan']]

    def avg_key(lst, key_fn):
        vals = [key_fn(r) for r in lst]
        return np.mean(vals), np.std(vals)

    metrics = [
        ("Diameter / dim", lambda r: r['graph']['diameter'] / r['dim']),
        ("Avg shortest path / dim", lambda r: r['graph']['avg_shortest_path'] / r['dim']),
        ("Avg graph complexity / dim", lambda r: r['graph']['mean_graph_complexity'] / r['dim']),
        ("Edge density", lambda r: r['graph']['edge_density']),
        ("Average OTOC", lambda r: r['graph']['avg_otoc']),
        ("Spectral gap", lambda r: r['spectral']['spectral_gap']),
        ("Spectral gap / d_max", lambda r: r['spectral']['spectral_gap'] / r['graph']['d_max']),
        ("Mixing time (RW)", lambda r: min(r['graph']['mixing_time_rw'], 1e6)),
        ("Cheeger lower bound", lambda r: r['graph']['cheeger_lower']),
        ("Ratio max|λ_nt| / 2√q", lambda r: r['spectral']['ratio']),
    ]

    print(f"\n  {'Metric':<30s}  {'Ramanujan (mean±std)':>22s}  {'Non-Ramanujan':>22s}  {'Implication'}")
    print("  " + "-" * 106)

    for label, fn in metrics:
        if ram:
            m_r, s_r = avg_key(ram, fn)
            r_str = f"{m_r:.4f} ± {s_r:.4f}"
        else:
            r_str = "N/A"
        if non_ram:
            m_n, s_n = avg_key(non_ram, fn)
            n_str = f"{m_n:.4f} ± {s_n:.4f}"
        else:
            n_str = "N/A"

        # Determine implication
        if ram and non_ram:
            if "Diameter" in label or "shortest" in label or "complexity" in label:
                imp = "Ram → faster spreading" if m_r < m_n else "---"
            elif "OTOC" in label:
                imp = "Ram → more scrambling" if m_r < m_n else "---"
            elif "gap" in label or "Cheeger" in label:
                imp = "Ram → better expansion" if m_r > m_n else "---"
            elif "Mixing" in label:
                imp = "Ram → faster mixing" if m_r < m_n else "---"
            elif "Ratio" in label:
                imp = "Ram → tighter bound" if m_r < m_n else "---"
            else:
                imp = ""
        else:
            imp = ""

        print(f"  {label:<30s}  {r_str:>22s}  {n_str:>22s}  {imp}")

    # ================================================================
    # Analysis Section 2: Scaling with Coxeter number
    # ================================================================
    print()
    print("=" * 110)
    print("  SCALING WITH COXETER NUMBER h")
    print("=" * 110)

    h_vals = np.array([r['h_coxeter'] for r in results if r['h_coxeter'] > 0])
    dims = np.array([r['dim'] for r in results if r['h_coxeter'] > 0])
    diams = np.array([r['graph']['diameter'] for r in results if r['h_coxeter'] > 0])
    avg_sps = np.array([r['graph']['avg_shortest_path'] for r in results if r['h_coxeter'] > 0])
    gaps = np.array([r['spectral']['spectral_gap'] for r in results if r['h_coxeter'] > 0])
    otocs = np.array([r['graph']['avg_otoc'] for r in results if r['h_coxeter'] > 0])
    mix_rw = np.array([r['graph']['mixing_time_rw'] for r in results if r['h_coxeter'] > 0])
    is_ram = np.array([r['spectral']['is_ramanujan'] for r in results if r['h_coxeter'] > 0])

    # Diameter scaling
    print(f"\n  {'h':>3s}  {'dim':>4s}  {'diam':>4s}  {'diam/dim':>8s}  "
          f"{'⟨ℓ⟩/dim':>8s}  {'gap':>7s}  {'gap/d_max':>9s}  "
          f"{'⟨OTOC⟩':>7s}  {'t_mix':>7s}  {'Ram':>4s}")
    print("  " + "-" * 85)

    filtered = [r for r in results if r['h_coxeter'] > 0]
    filtered.sort(key=lambda r: r['h_coxeter'])
    for r in filtered:
        h = r['h_coxeter']
        d = r['dim']
        diam = r['graph']['diameter']
        sp = r['graph']['avg_shortest_path']
        g = r['spectral']['spectral_gap']
        dm = r['graph']['d_max']
        ot = r['graph']['avg_otoc']
        tm = r['graph']['mixing_time_rw']
        ram = "YES" if r['spectral']['is_ramanujan'] else "NO"
        tm_str = f"{tm:7.1f}" if tm < 1e6 else "    inf"
        print(f"  {h:3d}  {d:4d}  {diam:4d}  {diam/d:8.4f}  "
              f"{sp/d:8.4f}  {g:7.3f}  {g/dm:9.5f}  "
              f"{ot:7.4f}  {tm_str}  {ram:>4s}")

    # ================================================================
    # Analysis Section 3: Krylov complexity lower bound
    # ================================================================
    print()
    print("=" * 110)
    print("  KRYLOV COMPLEXITY LOWER BOUND (West, Lemma 9)")
    print("=" * 110)
    print()
    print("  K(p_t) ≥ G(p_t)  and  ⟨G⟩_∞ = (1/|C|) Σ_q ℓ(p,q)")
    print()

    print(f"  {'Algebra':>10s}  {'dim':>4s}  {'⟨G⟩_Cartan':>11s}  "
          f"{'⟨G⟩_root':>10s}  {'⟨G⟩_total':>10s}  "
          f"{'diam':>4s}  {'Ram':>4s}")
    print("  " + "-" * 70)

    for r in filtered:
        dim = r['dim']
        rank = r['rank']
        gc = r['graph']['avg_graph_complexity_per_node']
        # Cartan elements: first 'rank' nodes
        gc_cartan = np.mean(gc[:rank]) if rank > 0 else 0
        # Root generators: remaining nodes
        gc_root = np.mean(gc[rank:]) if len(gc) > rank else 0
        gc_total = np.mean(gc)
        diam = r['graph']['diameter']
        ram = "YES" if r['spectral']['is_ramanujan'] else "NO"
        print(f"  {r['name']:>10s}  {dim:4d}  {gc_cartan:11.4f}  "
              f"{gc_root:10.4f}  {gc_total:10.4f}  "
              f"{diam:4d}  {ram:>4s}")

    # ================================================================
    # Analysis Section 4: Critical Coxeter number and chaos transition
    # ================================================================
    print()
    print("=" * 110)
    print("  RAMANUJAN TRANSITION → CHAOS TRANSITION")
    print("=" * 110)
    print()
    print("  The Ramanujan transition at h* ≈ 9.5 marks a phase transition")
    print("  in the graph's expansion properties. Physical implications:")
    print()
    print("  h < h*: Ramanujan → optimal expander → fastest OTOC decay →")
    print("           maximally efficient scrambling")
    print("  h > h*: Non-Ramanujan → sub-optimal expansion → slower scrambling →")
    print("           reduced quantum chaos")
    print()

    # Standard Model algebras
    sm_algebras = ['su(2)', 'su(3)']
    gut_algebras = ['su(5)', 'so(10)', 'E6', 'E8']

    print("  PHYSICS-RELEVANT ALGEBRAS:")
    print(f"  {'Algebra':>10s}  {'Context':>20s}  {'h':>3s}  "
          f"{'Ram?':>4s}  {'⟨OTOC⟩':>7s}  {'gap':>6s}  {'t_mix':>6s}")
    print("  " + "-" * 70)

    physics_map = {
        'su(2)': 'Weak isospin',
        'su(3)': 'QCD color',
        'su(5)': 'Georgi-Glashow GUT',
        'so(10)': 'SO(10) GUT',
        'E6': 'E6 GUT',
        'E8': 'HE string / E8×E8',
    }

    for r in filtered:
        if r['name'] in physics_map:
            context = physics_map[r['name']]
            h = r['h_coxeter']
            ram = "YES" if r['spectral']['is_ramanujan'] else "NO"
            ot = r['graph']['avg_otoc']
            g = r['spectral']['spectral_gap']
            tm = r['graph']['mixing_time_rw']
            tm_str = f"{tm:6.1f}" if tm < 1e6 else "   inf"
            print(f"  {r['name']:>10s}  {context:>20s}  {h:3d}  "
                  f"{ram:>4s}  {ot:7.4f}  {g:6.3f}  {tm_str}")

    # ================================================================
    # Analysis Section 5: Quantitative correlations
    # ================================================================
    print()
    print("=" * 110)
    print("  QUANTITATIVE CORRELATIONS")
    print("=" * 110)

    valid = [r for r in results if r['h_coxeter'] > 0]
    x_h = np.array([r['h_coxeter'] for r in valid], dtype=float)
    y_gap_norm = np.array([r['spectral']['spectral_gap'] / r['graph']['d_max']
                           for r in valid])
    y_diam_norm = np.array([r['graph']['diameter'] / r['dim'] for r in valid])
    y_sp_norm = np.array([r['graph']['avg_shortest_path'] / r['dim'] for r in valid])
    y_otoc = np.array([r['graph']['avg_otoc'] for r in valid])
    y_ratio = np.array([r['spectral']['ratio'] for r in valid])

    # Correlation coefficients
    corrs = [
        ("h vs ratio (λ/2√q)", x_h, y_ratio),
        ("h vs gap/d_max", x_h, y_gap_norm),
        ("h vs diam/dim", x_h, y_diam_norm),
        ("h vs ⟨ℓ⟩/dim", x_h, y_sp_norm),
        ("h vs ⟨OTOC⟩", x_h, y_otoc),
        ("gap/d_max vs diam/dim", y_gap_norm, y_diam_norm),
        ("gap/d_max vs ⟨ℓ⟩/dim", y_gap_norm, y_sp_norm),
        ("ratio vs ⟨ℓ⟩/dim", y_ratio, y_sp_norm),
    ]

    print(f"\n  {'Correlation':>30s}  {'r (Pearson)':>12s}  {'Interpretation'}")
    print("  " + "-" * 80)

    for label, x, y in corrs:
        if len(x) > 2:
            r_val = np.corrcoef(x, y)[0, 1]
            if abs(r_val) > 0.7:
                strength = "STRONG"
            elif abs(r_val) > 0.4:
                strength = "moderate"
            else:
                strength = "weak"
            print(f"  {label:>30s}  {r_val:+12.4f}  {strength}")
        else:
            print(f"  {label:>30s}  {'N/A':>12s}")

    return results


# ============================================================
# 5. Visualization
# ============================================================

def create_visualization(results):
    """Create 9-panel visualization of West et al. connections."""
    fig = plt.figure(figsize=(22, 18))
    fig.suptitle("Priority 4: Connection to West et al. — Quantum Chaos ↔ Ramanujan Property",
                 fontsize=15, fontweight='bold', y=0.98)
    gs = GridSpec(3, 3, hspace=0.35, wspace=0.35)

    valid = [r for r in results if r['h_coxeter'] > 0]
    h_vals = np.array([r['h_coxeter'] for r in valid])
    dims = np.array([r['dim'] for r in valid])
    is_ram = np.array([r['spectral']['is_ramanujan'] for r in valid])
    colors = ['#2196F3' if ram else '#FF5722' for ram in is_ram]
    markers = ['o' if ram else 's' for ram in is_ram]

    # Panel 1: Ratio vs Coxeter number (context from P2)
    ax1 = fig.add_subplot(gs[0, 0])
    ratios = np.array([r['spectral']['ratio'] for r in valid])
    for i in range(len(valid)):
        ax1.scatter(h_vals[i], ratios[i], c=colors[i], marker=markers[i],
                    s=60, edgecolors='black', linewidth=0.5, zorder=3)
    ax1.axhline(y=1.0, color='green', linestyle='--', linewidth=1.5, label='Ramanujan boundary')
    ax1.axvline(x=9.5, color='gray', linestyle=':', linewidth=1.5, label='h* ≈ 9.5')
    ax1.set_xlabel('Coxeter number h')
    ax1.set_ylabel('max|λ_nt| / 2√q')
    ax1.set_title('(a) Ramanujan ratio vs h')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Panel 2: Average OTOC vs Coxeter number
    ax2 = fig.add_subplot(gs[0, 1])
    otocs = np.array([r['graph']['avg_otoc'] for r in valid])
    for i in range(len(valid)):
        ax2.scatter(h_vals[i], otocs[i], c=colors[i], marker=markers[i],
                    s=60, edgecolors='black', linewidth=0.5, zorder=3)
    ax2.axvline(x=9.5, color='gray', linestyle=':', linewidth=1.5)
    ax2.set_xlabel('Coxeter number h')
    ax2.set_ylabel('⟨OTOC⟩ = 1 - 2f_ac')
    ax2.set_title('(b) Average OTOC vs h (West, Cor. 3)')
    ax2.grid(True, alpha=0.3)

    # Panel 3: Spectral gap (normalized by d_max) vs h
    ax3 = fig.add_subplot(gs[0, 2])
    gap_norm = np.array([r['spectral']['spectral_gap'] / r['graph']['d_max'] for r in valid])
    for i in range(len(valid)):
        ax3.scatter(h_vals[i], gap_norm[i], c=colors[i], marker=markers[i],
                    s=60, edgecolors='black', linewidth=0.5, zorder=3)
    ax3.axvline(x=9.5, color='gray', linestyle=':', linewidth=1.5)
    ax3.set_xlabel('Coxeter number h')
    ax3.set_ylabel('Spectral gap / d_max')
    ax3.set_title('(c) Normalized spectral gap vs h')
    ax3.grid(True, alpha=0.3)

    # Panel 4: Diameter / dim vs h
    ax4 = fig.add_subplot(gs[1, 0])
    diam_norm = np.array([r['graph']['diameter'] / r['dim'] for r in valid])
    for i in range(len(valid)):
        ax4.scatter(h_vals[i], diam_norm[i], c=colors[i], marker=markers[i],
                    s=60, edgecolors='black', linewidth=0.5, zorder=3)
    ax4.axvline(x=9.5, color='gray', linestyle=':', linewidth=1.5)
    ax4.set_xlabel('Coxeter number h')
    ax4.set_ylabel('Diameter / dim')
    ax4.set_title('(d) Normalized diameter vs h')
    ax4.grid(True, alpha=0.3)

    # Panel 5: Average graph complexity / dim vs h
    ax5 = fig.add_subplot(gs[1, 1])
    gc_norm = np.array([r['graph']['mean_graph_complexity'] / r['dim'] for r in valid])
    for i in range(len(valid)):
        ax5.scatter(h_vals[i], gc_norm[i], c=colors[i], marker=markers[i],
                    s=60, edgecolors='black', linewidth=0.5, zorder=3)
    ax5.axvline(x=9.5, color='gray', linestyle=':', linewidth=1.5)
    ax5.set_xlabel('Coxeter number h')
    ax5.set_ylabel('⟨G⟩ / dim')
    ax5.set_title('(e) Avg graph complexity / dim vs h (Krylov bound)')
    ax5.grid(True, alpha=0.3)

    # Panel 6: Mixing time vs h
    ax6 = fig.add_subplot(gs[1, 2])
    mix_times = np.array([r['graph']['mixing_time_rw'] for r in valid])
    for i in range(len(valid)):
        ax6.scatter(h_vals[i], mix_times[i], c=colors[i], marker=markers[i],
                    s=60, edgecolors='black', linewidth=0.5, zorder=3)
    ax6.axvline(x=9.5, color='gray', linestyle=':', linewidth=1.5)
    ax6.set_xlabel('Coxeter number h')
    ax6.set_ylabel('t_mix bound (RW)')
    ax6.set_title('(f) Mixing time bound vs h')
    ax6.grid(True, alpha=0.3)
    ax6.set_yscale('log')

    # Panel 7: Spectral gap vs diameter (expansion ↔ spreading)
    ax7 = fig.add_subplot(gs[2, 0])
    gaps = np.array([r['spectral']['spectral_gap'] for r in valid])
    diams = np.array([r['graph']['diameter'] for r in valid])
    for i in range(len(valid)):
        ax7.scatter(gaps[i], diams[i], c=colors[i], marker=markers[i],
                    s=60, edgecolors='black', linewidth=0.5, zorder=3)
        ax7.annotate(valid[i]['name'], (gaps[i], diams[i]),
                     fontsize=5, alpha=0.7, ha='left')
    ax7.set_xlabel('Spectral gap')
    ax7.set_ylabel('Diameter')
    ax7.set_title('(g) Gap vs Diameter — expansion ↔ spreading')
    ax7.grid(True, alpha=0.3)

    # Panel 8: OTOC vs spectral gap / d_max
    ax8 = fig.add_subplot(gs[2, 1])
    for i in range(len(valid)):
        ax8.scatter(gap_norm[i], otocs[i], c=colors[i], marker=markers[i],
                    s=60, edgecolors='black', linewidth=0.5, zorder=3)
    ax8.set_xlabel('Spectral gap / d_max')
    ax8.set_ylabel('⟨OTOC⟩')
    ax8.set_title('(h) OTOC ↔ normalized spectral gap')
    ax8.grid(True, alpha=0.3)

    # Panel 9: Graph complexity vs Cheeger lower bound
    ax9 = fig.add_subplot(gs[2, 2])
    cheeger_l = np.array([r['graph']['cheeger_lower'] for r in valid])
    for i in range(len(valid)):
        ax9.scatter(cheeger_l[i], gc_norm[i], c=colors[i], marker=markers[i],
                    s=60, edgecolors='black', linewidth=0.5, zorder=3)
    ax9.set_xlabel('Cheeger lower bound')
    ax9.set_ylabel('⟨G⟩ / dim')
    ax9.set_title('(i) Krylov bound ↔ Cheeger constant')
    ax9.grid(True, alpha=0.3)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#2196F3',
               markersize=10, label='Ramanujan'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#FF5722',
               markersize=10, label='Non-Ramanujan'),
        Line2D([0], [0], color='gray', linestyle=':', linewidth=1.5,
               label='h* ≈ 9.5 (Ramanujan transition)'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3,
               fontsize=10, frameon=True, bbox_to_anchor=(0.5, 0.01))

    outpath = os.path.join(os.path.dirname(__file__), 'west_connection.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\n  [Saved figure: {outpath}]")
    plt.close()


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    results = run_full_analysis()
    create_visualization(results)
    print("\n  DONE. West et al. connection analysis complete.")
