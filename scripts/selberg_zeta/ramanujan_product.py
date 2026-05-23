#!/usr/bin/env python3
"""
P7: Ramanujan analysis of PRODUCT (direct-sum) Lie algebras.

All P1-P6 work tested simple Lie algebras only.  The Standard Model gauge
group is su(3) + su(2) + u(1) — a direct sum, NOT simple.  This script
answers:

  1. Is the SM gauge algebra "spectrally optimal" among physically motivated
     product algebras?
  2. Does the Ramanujan property compose under direct sums?
     (All factors Ramanujan => product Ramanujan?)
  3. How does disconnectedness of the adjoint graph affect the 6 definitions?
  4. How does the SM compare to GUT unifications (SU(5), SO(10), E6)?

The adjoint commutation graph of g1 + g2 is the DISJOINT UNION of the
individual adjoint graphs (because [g1, g2] = 0).  Thus:
  - The spectrum is the union of the component spectra
  - The graph is disconnected (one component per simple factor)
  - u(1) contributes an isolated vertex (degree 0)

Author: SS-PEG project
Date:   June 2025
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(__file__))
from ramanujan_extended import (
    roots_A, roots_B, roots_C, roots_D,
    roots_G2, roots_F4, roots_E6, roots_E7, roots_E8,
    adjacency_from_roots, spectral_analysis
)

HASH_DIM_LIMIT = 200


# ── Direct-sum adjacency ─────────────────────────────────────────────

def adjacency_direct_sum(*adj_list):
    """Block-diagonal adjacency for a direct sum of Lie algebras."""
    total = sum(a.shape[0] for a in adj_list)
    result = np.zeros((total, total), dtype=np.float64)
    offset = 0
    for a in adj_list:
        n = a.shape[0]
        result[offset:offset + n, offset:offset + n] = a
        offset += n
    return result


def adjacency_u1():
    """u(1) is abelian: single vertex, no edges."""
    return np.zeros((1, 1), dtype=np.float64)


# ── Simple-algebra adjacency by name ─────────────────────────────────

_ROOT_FNS = {
    'su(2)':  lambda: roots_A(1),
    'su(3)':  lambda: roots_A(2),
    'su(4)':  lambda: roots_A(3),
    'su(5)':  lambda: roots_A(4),
    'su(6)':  lambda: roots_A(5),
    'su(7)':  lambda: roots_A(6),
    'su(8)':  lambda: roots_A(7),
    'so(5)':  lambda: roots_B(2),
    'so(7)':  lambda: roots_B(3),
    'so(9)':  lambda: roots_B(4),
    'so(10)': lambda: roots_D(5),
    'so(6)':  lambda: roots_D(3),
    'sp(4)':  lambda: roots_C(2),
    'sp(6)':  lambda: roots_C(3),
    'G2':     roots_G2,
    'F4':     roots_F4,
    'E6':     roots_E6,
    'E7':     roots_E7,
    'E8':     roots_E8,
}


def get_factor_adj(name):
    """Return adjacency matrix for a simple Lie algebra or u(1)."""
    if name == 'u(1)':
        return adjacency_u1()
    fn = _ROOT_FNS.get(name)
    if fn is None:
        raise ValueError(f"Unknown algebra: {name}")
    rank, roots, simple = fn()
    return adjacency_from_roots(rank, roots, simple)


# ── Product algebra catalog ──────────────────────────────────────────

def build_product_catalog():
    """Physically motivated product algebras + simple GUT references."""
    cat = OrderedDict()

    # ── Standard Model and electroweak sub-groups ──
    cat['EW'] = {
        'desc': 'Electroweak SU(2)xU(1)',
        'factors': ['su(2)', 'u(1)'], 'category': 'SM'}
    cat['QCD+QED'] = {
        'desc': 'SU(3)xU(1)',
        'factors': ['su(3)', 'u(1)'], 'category': 'SM'}
    cat['SM_no_u1'] = {
        'desc': 'SM without U(1)',
        'factors': ['su(3)', 'su(2)'], 'category': 'SM'}
    cat['SM'] = {
        'desc': 'Standard Model SU(3)xSU(2)xU(1)',
        'factors': ['su(3)', 'su(2)', 'u(1)'], 'category': 'SM'}

    # ── GUT product groups ──
    cat['Pati-Salam'] = {
        'desc': 'SU(4)xSU(2)_LxSU(2)_R',
        'factors': ['su(4)', 'su(2)', 'su(2)'], 'category': 'GUT'}
    cat['Left-Right'] = {
        'desc': 'SU(3)xSU(2)_LxSU(2)_RxU(1)',
        'factors': ['su(3)', 'su(2)', 'su(2)', 'u(1)'], 'category': 'GUT'}
    cat['Trinification'] = {
        'desc': 'SU(3)^3',
        'factors': ['su(3)', 'su(3)', 'su(3)'], 'category': 'GUT'}
    cat['Flipped_SU5'] = {
        'desc': 'SU(5)xU(1)',
        'factors': ['su(5)', 'u(1)'], 'category': 'GUT'}

    # ── Simple GUT references (1 factor, for comparison) ──
    cat['SU(5)_GUT'] = {
        'desc': 'Georgi-Glashow SU(5)',
        'factors': ['su(5)'], 'category': 'GUT_simple'}
    cat['SO(10)_GUT'] = {
        'desc': 'SO(10)',
        'factors': ['so(10)'], 'category': 'GUT_simple'}
    cat['E6_GUT'] = {
        'desc': 'E6 GUT',
        'factors': ['E6'], 'category': 'GUT_simple'}

    # ── Alternative / "what-if" product groups ──
    cat['su2_x_su2'] = {
        'desc': 'SU(2)xSU(2) ~ Spin(4)',
        'factors': ['su(2)', 'su(2)'], 'category': 'Alt'}
    cat['su3_x_su3'] = {
        'desc': 'SU(3)xSU(3) chiral',
        'factors': ['su(3)', 'su(3)'], 'category': 'Alt'}
    cat['su4_x_su2_x_u1'] = {
        'desc': 'SU(4)xSU(2)xU(1)',
        'factors': ['su(4)', 'su(2)', 'u(1)'], 'category': 'Alt'}
    cat['su3_x_su3_x_u1'] = {
        'desc': 'SU(3)xSU(3)xU(1)',
        'factors': ['su(3)', 'su(3)', 'u(1)'], 'category': 'Alt'}
    cat['so5_x_su2_x_u1'] = {
        'desc': 'SO(5)xSU(2)xU(1)',
        'factors': ['so(5)', 'su(2)', 'u(1)'], 'category': 'Alt'}
    cat['G2_x_su2'] = {
        'desc': 'G2xSU(2)',
        'factors': ['G2', 'su(2)'], 'category': 'Alt'}
    cat['su6_x_su2'] = {
        'desc': 'SU(6)xSU(2)',
        'factors': ['su(6)', 'su(2)'], 'category': 'Alt'}
    cat['su8'] = {
        'desc': 'SU(8) simple (boundary)',
        'factors': ['su(8)'], 'category': 'Alt'}

    # Compute expected dimensions
    _dims = {
        'u(1)': 1, 'su(2)': 3, 'su(3)': 8, 'su(4)': 15, 'su(5)': 24,
        'su(6)': 35, 'su(7)': 48, 'su(8)': 63, 'so(5)': 10, 'so(6)': 15,
        'so(7)': 21, 'so(9)': 36, 'so(10)': 45, 'sp(4)': 10, 'sp(6)': 21,
        'G2': 14, 'F4': 52, 'E6': 78, 'E7': 133, 'E8': 248,
    }
    for info in cat.values():
        info['dim'] = sum(_dims[f] for f in info['factors'])

    return cat


# ── Connected-component decomposition ────────────────────────────────

def find_components(adj):
    """Return list of (indices, sub_adjacency) for each connected component."""
    n = adj.shape[0]
    visited = np.zeros(n, dtype=bool)
    components = []
    for start in range(n):
        if visited[start]:
            continue
        queue = [start]
        visited[start] = True
        comp = [start]
        while queue:
            node = queue.pop(0)
            for nb in range(n):
                if adj[node, nb] > 0.5 and not visited[nb]:
                    visited[nb] = True
                    queue.append(nb)
                    comp.append(nb)
        idx = sorted(comp)
        sub = adj[np.ix_(idx, idx)]
        components.append((idx, sub))
    return components


def per_component_analysis(adj):
    """Spectral analysis of each connected component separately."""
    components = find_components(adj)
    results = []
    for idx, sub in components:
        if len(idx) == 1:
            results.append({
                'size': 1, 'is_ramanujan': True, 'max_nontrivial': 0.0,
                'ramanujan_bound': 0.0, 'ratio': 0.0, 'degree': 0.0,
                'eigenvalues': np.array([0.0]), 'spectral_gap': 0.0,
                'q': 0.0, 'lambda_max': 0.0, 'is_regular': True,
            })
        else:
            sa = spectral_analysis(sub)
            sa['is_ramanujan'] = bool(sa['is_ramanujan'])
            sa['size'] = len(idx)
            results.append(sa)
    return results


# ── Component-aware nontrivial eigenvalues ───────────────────────────

def component_aware_max_nontrivial(adj):
    """
    For a disconnected graph, the "nontrivial" eigenvalues should exclude
    each component's trivial (largest) eigenvalue.  The naive whole-graph
    spectral_analysis only removes the global lambda_max, counting the
    lambda_max of smaller components as nontrivial.

    Returns (max_nt_aware, per_component_trivials, all_nontrivials).
    """
    components = find_components(adj)
    all_nt = []
    trivials = []
    for idx, sub in components:
        if len(idx) == 1:
            trivials.append(0.0)
            continue
        eigs = np.sort(np.linalg.eigvalsh(sub))[::-1]
        trivials.append(float(eigs[0]))
        # Nontrivial for this component
        nt = eigs[1:]
        degrees = sub.sum(axis=1)
        is_reg = np.allclose(degrees, degrees[0])
        if is_reg and len(eigs) > 2 and abs(eigs[-1] + eigs[0]) < 0.01:
            nt = eigs[1:-1]
        all_nt.extend(np.abs(nt).tolist())
    max_nt = max(all_nt) if all_nt else 0.0
    return max_nt, trivials, all_nt


# ── Whole-graph 6 definitions (disconnection-aware) ──────────────────

def whole_graph_6defs(adj):
    """
    Compute Ramanujan status under 6 definitions for the full product graph.
    Uses component-aware nontrivial extraction.
    """
    n = adj.shape[0]
    eigs = np.sort(np.linalg.eigvalsh(adj))[::-1]
    degrees = adj.sum(axis=1)

    # Component-aware max nontrivial
    max_nt_aware, trivials, _ = component_aware_max_nontrivial(adj)
    # Naive max nontrivial (for comparison)
    max_nt_naive = float(np.max(np.abs(eigs[1:]))) if n > 1 else 0.0

    # Degree statistics (only non-isolated vertices for harm/geom)
    d_mean = float(np.mean(degrees))
    d_max = float(np.max(degrees))
    nonzero = degrees[degrees > 0]
    if len(nonzero) > 0:
        d_harm = float(len(nonzero) / np.sum(1.0 / nonzero))
        d_geom_m1 = float(np.exp(np.mean(np.log(np.maximum(nonzero - 1, 1e-30)))))
    else:
        d_harm = 0.0
        d_geom_m1 = 0.0
    sorted_deg = np.sort(degrees)[::-1]
    d_2nd = float(sorted_deg[1]) if n > 1 else d_mean

    defs = {}
    for label, q_eff in [('D1_mean', d_mean - 1), ('D2_dmax', d_max - 1),
                          ('D3_geom', d_geom_m1), ('D5_harm', d_harm - 1),
                          ('D6_d2', d_2nd - 1)]:
        bound = 2 * np.sqrt(max(q_eff, 0))
        is_ram = bool(max_nt_aware <= bound + 1e-10)
        ratio = max_nt_aware / bound if bound > 0 else float('inf')
        defs[label] = {'q': q_eff, 'bound': bound, 'is_ram': is_ram, 'ratio': ratio}

    # D4 Hashimoto
    if n <= HASH_DIM_LIMIT:
        edges = []
        for i in range(n):
            for j in range(n):
                if adj[i, j] > 0.5:
                    edges.append((i, j))
        m = len(edges)
        if m > 0:
            T_e = np.zeros((m, m), dtype=np.float64)
            edge_index = {e: idx for idx, e in enumerate(edges)}
            for idx_e, (i, j) in enumerate(edges):
                for k in range(n):
                    if k != i and adj[j, k] > 0.5:
                        jk = edge_index.get((j, k))
                        if jk is not None:
                            T_e[jk, idx_e] = 1.0
            hash_eigs = np.linalg.eigvals(T_e)
            rho = float(np.max(np.abs(hash_eigs)))
            bound4 = 2 * np.sqrt(max(rho, 0))
            is_ram4 = bool(max_nt_aware <= bound4 + 1e-10)
            ratio4 = max_nt_aware / bound4 if bound4 > 0 else float('inf')
            defs['D4_hash'] = {'q': rho, 'bound': bound4, 'is_ram': is_ram4, 'ratio': ratio4}
        else:
            defs['D4_hash'] = {'q': 0, 'bound': 0, 'is_ram': True, 'ratio': 0}
    else:
        defs['D4_hash'] = None

    return {
        'max_nt_aware': max_nt_aware,
        'max_nt_naive': max_nt_naive,
        'lambda_max': float(eigs[0]),
        'trivials': trivials,
        'degrees': degrees,
        'd_mean': d_mean, 'd_max': d_max, 'd_harm': d_harm,
        'd_geom_m1': d_geom_m1, 'd_2nd': d_2nd,
        'eigs': eigs,
        'definitions': defs,
    }


# ── Composability ────────────────────────────────────────────────────

def test_composability(factor_specs, product_defs):
    """Test whether Ramanujan composes: all factors Ram => product Ram?"""
    all_factors_ram = all(bool(f['is_ramanujan']) for f in factor_specs)
    d1 = product_defs['definitions']['D1_mean']
    product_ram = d1['is_ram']
    return {
        'all_factors_ram': all_factors_ram,
        'product_ram': product_ram,
        'composes': all_factors_ram == product_ram,
        'fwd_fail': all_factors_ram and not product_ram,
        'bwd_fail': product_ram and not all_factors_ram,
    }


# ── Optimality score ─────────────────────────────────────────────────

def compute_scores(results):
    """Spectral optimality: lower = more Ramanujan."""
    scores = {}
    for name, res in results.items():
        defs = res['w6']['definitions']
        tested = [d for d in defs.values() if d is not None]
        n_ram = sum(1 for d in tested if d['is_ram'])
        ratios = [d['ratio'] for d in tested if d['ratio'] < float('inf')]
        avg_ratio = np.mean(ratios) if ratios else float('inf')
        comps = res['components']
        max_gap = max((c['spectral_gap'] for c in comps), default=0)
        scores[name] = {
            'n_ram': n_ram,
            'n_tested': len(tested),
            'avg_ratio': avg_ratio,
            'max_gap': max_gap,
            'dim': res['dim'],
            'n_comp': len(comps),
            'score': (len(tested) - n_ram) * 0.5 + avg_ratio,
        }
    return scores


# ── Visualization ────────────────────────────────────────────────────

def create_figure(results, scores, output_path):
    """6-panel visualization for P7."""
    fig, axes = plt.subplots(2, 3, figsize=(19, 12))
    fig.suptitle('P7: Ramanujan Analysis of Product Lie Algebras',
                 fontsize=14, fontweight='bold')

    names = list(results.keys())
    cats = [results[n]['category'] for n in names]
    cat_cmap = {'SM': '#2196F3', 'GUT': '#FF5722', 'GUT_simple': '#E64A19',
                'Alt': '#9C27B0'}

    def_keys = ['D1_mean', 'D2_dmax', 'D3_geom', 'D4_hash', 'D5_harm', 'D6_d2']
    def_labels = ['D1\nmean', 'D2\ndmax', 'D3\ngeom', 'D4\nHash', 'D5\nharm', 'D6\nd2']

    # ── Panel 1: Classification heatmap ──
    ax = axes[0, 0]
    mat = np.full((len(names), 6), 0.5)
    for i, n in enumerate(names):
        for j, dk in enumerate(def_keys):
            d = results[n]['w6']['definitions'].get(dk)
            if d is None:
                mat[i, j] = 0.5
            elif d['is_ram']:
                mat[i, j] = 1.0
            else:
                mat[i, j] = 0.0
    cmap3 = ListedColormap(['#EF5350', '#BDBDBD', '#66BB6A'])
    ax.imshow(mat, aspect='auto', cmap=cmap3, vmin=0, vmax=1, interpolation='nearest')
    ax.set_xticks(range(6))
    ax.set_xticklabels(['D1', 'D2', 'D3', 'D4', 'D5', 'D6'], fontsize=8)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_title('Classification Heatmap (6 defs)')
    for i, n in enumerate(names):
        if n == 'SM':
            ax.axhspan(i - 0.5, i + 0.5, color='gold', alpha=0.25)

    # ── Panel 2: Optimality ranking ──
    ax = axes[0, 1]
    srt = sorted(scores.keys(), key=lambda n: scores[n]['score'])
    vals = [scores[n]['score'] for n in srt]
    cols = [cat_cmap.get(results[n]['category'], '#607D8B') for n in srt]
    bars = ax.barh(range(len(srt)), vals, color=cols)
    ax.set_yticks(range(len(srt)))
    ax.set_yticklabels(srt, fontsize=7)
    ax.set_xlabel('Optimality score (lower = better)')
    ax.set_title('Spectral Optimality Ranking')
    for i, n in enumerate(srt):
        if n == 'SM':
            bars[i].set_edgecolor('gold')
            bars[i].set_linewidth(3)

    # ── Panel 3: D1 ratio bar chart ──
    ax = axes[0, 2]
    ratios_d1 = [results[n]['w6']['definitions']['D1_mean']['ratio'] for n in names]
    cols = [cat_cmap.get(results[n]['category'], '#607D8B') for n in names]
    ax.bar(range(len(names)), ratios_d1, color=cols)
    ax.axhline(1.0, color='red', ls='--', lw=1.5, label='Ramanujan boundary')
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=55, ha='right', fontsize=6)
    ax.set_ylabel('D1 ratio')
    ax.set_title('D1 (Mean-Degree) Ratio')
    ax.legend(fontsize=7)

    # ── Panel 4: Dimension vs ratio scatter ──
    ax = axes[1, 0]
    dims = [results[n]['dim'] for n in names]
    for c in cat_cmap:
        mask = [i for i, n in enumerate(names) if results[n]['category'] == c]
        if mask:
            ax.scatter([dims[i] for i in mask], [ratios_d1[i] for i in mask],
                       c=cat_cmap[c], label=c, s=80, edgecolors='black', zorder=3)
    ax.axhline(1.0, color='red', ls='--', lw=1.5)
    ax.set_xlabel('Dimension')
    ax.set_ylabel('D1 ratio')
    ax.set_title('Dimension vs Ramanujan Ratio')
    ax.legend(fontsize=7)
    for i, n in enumerate(names):
        ax.annotate(n, (dims[i], ratios_d1[i]), fontsize=5.5, alpha=0.7,
                    xytext=(3, 3), textcoords='offset points')

    # ── Panel 5: Composability test ──
    ax = axes[1, 1]
    prod_names = [n for n in names if len(results[n]['factors']) > 1]
    if prod_names:
        cdata = np.zeros((len(prod_names), 3))
        for i, n in enumerate(prod_names):
            c = results[n]['composability']
            cdata[i, 0] = 1 if c['all_factors_ram'] else 0
            cdata[i, 1] = 1 if c['product_ram'] else 0
            cdata[i, 2] = 1 if c['composes'] else 0
        cmap2 = ListedColormap(['#EF5350', '#66BB6A'])
        ax.imshow(cdata, aspect='auto', cmap=cmap2, vmin=0, vmax=1,
                  interpolation='nearest')
        ax.set_xticks(range(3))
        ax.set_xticklabels(['All factors\nRamanujan', 'Product\nRamanujan',
                            'Composes'], fontsize=8)
        ax.set_yticks(range(len(prod_names)))
        ax.set_yticklabels(prod_names, fontsize=7)
    ax.set_title('Composability: Factors -> Product')

    # ── Panel 6: Naive vs component-aware max nontrivial ──
    ax = axes[1, 2]
    naive = [results[n]['w6']['max_nt_naive'] for n in names]
    aware = [results[n]['w6']['max_nt_aware'] for n in names]
    x = np.arange(len(names))
    w = 0.35
    ax.bar(x - w / 2, naive, w, label='Naive (whole-graph)', color='#90CAF9')
    ax.bar(x + w / 2, aware, w, label='Component-aware', color='#1565C0')
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=55, ha='right', fontsize=6)
    ax.set_ylabel('max |lambda_nontrivial|')
    ax.set_title('Naive vs Component-Aware |lambda_nt|')
    ax.legend(fontsize=7)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved: {output_path}")


# ── Main analysis ────────────────────────────────────────────────────

def run_analysis():
    catalog = build_product_catalog()
    results = OrderedDict()

    print("=" * 90)
    print("P7: RAMANUJAN ANALYSIS OF PRODUCT LIE ALGEBRAS")
    print("=" * 90)

    # ── Phase 1: Factor spectral cache ──
    print("\n--- Phase 1: Individual Factor Analysis ---")
    factor_cache = {}
    for info in catalog.values():
        for f in info['factors']:
            if f not in factor_cache:
                a = get_factor_adj(f)
                sa = spectral_analysis(a)
                sa['is_ramanujan'] = bool(sa['is_ramanujan'])
                factor_cache[f] = {'adj': a, 'spec': sa}
                dim = a.shape[0]
                tag = "RAM" if sa['is_ramanujan'] else "non-Ram"
                print(f"  {f:10s}  dim={dim:3d}  d_mean={sa['degree']:.2f}  "
                      f"|lam_nt|={sa['max_nontrivial']:.4f}  "
                      f"bound={sa['ramanujan_bound']:.4f}  "
                      f"ratio={sa['ratio']:.4f}  {tag}")

    # ── Phase 2: Product graph construction & 6-definition analysis ──
    print("\n--- Phase 2: Product Analysis (component-aware) ---")
    hdr = (f"{'Name':20s} {'Dim':>4s} {'#C':>3s} {'d_mean':>6s} "
           f"{'|lam|_ca':>8s} {'|lam|_nv':>8s} "
           f"{'D1':>3s} {'D2':>3s} {'D3':>3s} {'D4':>3s} {'D5':>3s} {'D6':>3s}")
    print(hdr)
    print("-" * len(hdr))

    for name, info in catalog.items():
        adj_list = [factor_cache[f]['adj'] for f in info['factors']]
        prod_adj = adjacency_direct_sum(*adj_list) if len(adj_list) > 1 else adj_list[0]
        dim = prod_adj.shape[0]

        comps = per_component_analysis(prod_adj)
        w6 = whole_graph_6defs(prod_adj)
        fspecs = [factor_cache[f]['spec'] for f in info['factors']]
        comp_test = test_composability(fspecs, w6)

        results[name] = {
            'adj': prod_adj, 'dim': dim, 'factors': info['factors'],
            'category': info['category'], 'desc': info['desc'],
            'components': comps, 'w6': w6,
            'factor_specs': fspecs, 'composability': comp_test,
        }

        defs = w6['definitions']

        def sym(dk):
            d = defs.get(dk)
            return 'N/A' if d is None else ('Y' if d['is_ram'] else 'N')

        print(f"{name:20s} {dim:4d} {len(comps):3d} {w6['d_mean']:6.2f} "
              f"{w6['max_nt_aware']:8.4f} {w6['max_nt_naive']:8.4f} "
              f"{sym('D1_mean'):>3s} {sym('D2_dmax'):>3s} {sym('D3_geom'):>3s} "
              f"{sym('D4_hash'):>3s} {sym('D5_harm'):>3s} {sym('D6_d2'):>3s}")

    # ── Phase 3: Composability of Ramanujan property ──
    print("\n--- Phase 3: Composability ---")
    print(f"{'Product':20s} {'Factors':40s} {'AllFRam':>7s} {'ProdRam':>7s} {'=?':>3s}")
    print("-" * 82)
    compose_ok = compose_fail = 0
    for name, res in results.items():
        if len(res['factors']) > 1:
            c = res['composability']
            ftxt = ' x '.join(res['factors'])
            af = 'YES' if c['all_factors_ram'] else 'NO'
            pr = 'YES' if c['product_ram'] else 'NO'
            eq = 'YES' if c['composes'] else 'FAIL'
            print(f"{name:20s} {ftxt:40s} {af:>7s} {pr:>7s} {eq:>4s}")
            if c['composes']:
                compose_ok += 1
            else:
                compose_fail += 1

    print(f"\nComposability: {compose_ok} OK, {compose_fail} FAIL "
          f"out of {compose_ok + compose_fail} products")

    # ── Phase 4: Naive vs component-aware discrepancies ──
    print("\n--- Phase 4: Naive vs Component-Aware Nontrivial ---")
    print("(For disconnected graphs, naive max|lam_nt| can include the")
    print(" trivial eigenvalue of a smaller component.)\n")
    print(f"{'Name':20s} {'|lam|_naive':>11s} {'|lam|_aware':>11s} "
          f"{'Delta':>7s} {'Trivials of components':s}")
    print("-" * 80)
    discrepancies = 0
    for name, res in results.items():
        w = res['w6']
        delta = abs(w['max_nt_naive'] - w['max_nt_aware'])
        trv = ', '.join(f"{t:.2f}" for t in w['trivials'])
        flag = ' ***' if delta > 0.01 else ''
        print(f"{name:20s} {w['max_nt_naive']:11.4f} {w['max_nt_aware']:11.4f} "
              f"{delta:7.4f} [{trv}]{flag}")
        if delta > 0.01:
            discrepancies += 1
    print(f"\nDiscrepancies (Delta > 0.01): {discrepancies}")

    # ── Phase 5: SM vs GUT unification paths ──
    print("\n--- Phase 5: SM vs GUT Unification Paths ---")
    sm = results['SM']
    sm_d1 = sm['w6']['definitions']['D1_mean']
    sm_nram = sum(1 for d in sm['w6']['definitions'].values()
                  if d is not None and d['is_ram'])
    print(f"\n  SM:  dim={sm['dim']}  D1 ratio={sm_d1['ratio']:.4f}  "
          f"#Ramanujan defs={sm_nram}/6")
    print(f"       Components: {len(sm['components'])}  "
          f"Component-aware |lam_nt|={sm['w6']['max_nt_aware']:.4f}")

    gut_seq = ['Pati-Salam', 'Left-Right', 'Trinification', 'Flipped_SU5',
               'SU(5)_GUT', 'SO(10)_GUT', 'E6_GUT']
    for gn in gut_seq:
        if gn not in results:
            continue
        gr = results[gn]
        gd1 = gr['w6']['definitions']['D1_mean']
        gnr = sum(1 for d in gr['w6']['definitions'].values()
                  if d is not None and d['is_ram'])
        vs = "BETTER" if gd1['ratio'] < sm_d1['ratio'] else "WORSE"
        print(f"\n  {gn}:  dim={gr['dim']}  D1 ratio={gd1['ratio']:.4f}  "
              f"#Ram={gnr}/6  -> {vs} than SM")

    # ── Phase 6: Per-component structure ──
    print("\n--- Phase 6: Component Structure ---")
    for name, res in results.items():
        if len(res['factors']) <= 1:
            continue
        comps = res['components']
        wg = res['w6']
        print(f"\n{name} ({' x '.join(res['factors'])}): "
              f"{len(comps)} components, dim={res['dim']}")
        for i, c in enumerate(comps):
            if c['size'] == 1:
                print(f"  Comp {i + 1}: size=1 (isolated, u(1))")
            else:
                tag = "RAM" if c['is_ramanujan'] else "non-Ram"
                print(f"  Comp {i + 1}: size={c['size']}  d={c['degree']:.1f}  "
                      f"|lam_nt|={c['max_nontrivial']:.4f}  "
                      f"bound={c['ramanujan_bound']:.4f}  {tag}")
        # Compare
        all_cr = all(c['is_ramanujan'] for c in comps)
        w_ram = wg['definitions']['D1_mean']['is_ram']
        if all_cr != w_ram:
            print(f"  ** DISCREPANCY: per-component={all_cr}, "
                  f"whole-graph D1={w_ram}")
        else:
            print(f"  Consistent: per-component={all_cr}, whole-graph D1={w_ram}")

    # ── Phase 7: Optimality ranking ──
    scores = compute_scores(results)
    print("\n--- Phase 7: Spectral Optimality Ranking ---")
    srt = sorted(scores.keys(), key=lambda n: scores[n]['score'])
    print(f"{'Rank':>4s} {'Name':20s} {'Dim':>4s} {'#Ram':>5s} "
          f"{'AvgRatio':>9s} {'Score':>7s} {'Cat':>10s}")
    print("-" * 68)
    sm_rank = None
    for i, n in enumerate(srt):
        s = scores[n]
        mark = ' ***' if n == 'SM' else ''
        if n == 'SM':
            sm_rank = i + 1
        print(f"{i + 1:4d} {n:20s} {s['dim']:4d} {s['n_ram']:3d}/{s['n_tested']}"
              f" {s['avg_ratio']:9.4f} {s['score']:7.4f} {results[n]['category']:>10s}{mark}")

    # ── Summary ──
    print("\n" + "=" * 90)
    print("SUMMARY: P7 KEY FINDINGS")
    print("=" * 90)

    print(f"\n1. SM optimality rank: #{sm_rank}/{len(srt)}")
    print(f"   Score = {scores['SM']['score']:.4f}, "
          f"avg ratio = {scores['SM']['avg_ratio']:.4f}")

    # Count unanimous Ramanujan products
    n_all_ram = sum(
        1 for n in results
        if all(d['is_ram'] for d in results[n]['w6']['definitions'].values()
               if d is not None))
    print(f"\n2. Unanimous Ramanujan (all defs): {n_all_ram}/{len(results)}")

    print(f"\n3. Composability: {compose_ok} OK, {compose_fail} FAIL")

    print(f"\n4. Naive vs component-aware discrepancies: {discrepancies}")

    # SM vs simple GUTs
    for gn in ['SU(5)_GUT', 'SO(10)_GUT', 'E6_GUT']:
        if gn in results:
            gr = results[gn]['w6']['definitions']['D1_mean']['ratio']
            vs = "better" if gr < sm_d1['ratio'] else "worse"
            print(f"\n5. {gn}: ratio={gr:.4f} ({vs} than SM {sm_d1['ratio']:.4f})")

    # Alternatives better than SM
    better = [n for n in results
              if results[n]['category'] == 'Alt'
              and results[n]['w6']['definitions']['D1_mean']['ratio'] < sm_d1['ratio']]
    print(f"\n6. Alternative products with lower ratio than SM: {len(better)}")
    for b in better:
        r = results[b]['w6']['definitions']['D1_mean']['ratio']
        print(f"   - {b}: ratio={r:.4f}")

    # Figure
    out = os.path.join(os.path.dirname(__file__), 'ramanujan_product.png')
    create_figure(results, scores, out)


if __name__ == '__main__':
    run_analysis()
