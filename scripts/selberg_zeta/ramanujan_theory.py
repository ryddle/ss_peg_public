"""
Theoretical Analysis: Why the Ramanujan Ratio Grows with Rank
=============================================================

Investigates the relationship between root system invariants
and the spectral properties of the adjoint commutation graph.

Key invariants per algebra:
  - rank r
  - number of positive roots |Phi+|
  - dimension = r + 2|Phi+|
  - Coxeter number h
  - dual Coxeter number h_v
  - ratio of long/short root lengths
  - degree distribution {deg(H_i), deg(E_alpha)}
  - average degree, max degree, min degree
  - spectral radius, second eigenvalue, Ramanujan ratio

Goal: find an analytical formula ratio ~ f(root system invariants)
and understand the Ramanujan transition.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
from collections import OrderedDict

from ramanujan_extended import (
    roots_A, roots_B, roots_C, roots_D,
    roots_G2, roots_F4, roots_E6, roots_E7, roots_E8,
    adjacency_from_roots, spectral_analysis,
)


# ============================================================
# Root system invariants
# ============================================================

COXETER_NUMBERS = {
    'A': lambda n: n + 1,       # su(n+1): h = n+1
    'B': lambda n: 2 * n,       # so(2n+1): h = 2n
    'C': lambda n: 2 * n,       # sp(2n): h = 2n
    'D': lambda n: 2 * (n - 1), # so(2n): h = 2(n-1)
    'G2': 6,
    'F4': 12,
    'E6': 12,
    'E7': 18,
    'E8': 30,
}

DUAL_COXETER_NUMBERS = {
    'A': lambda n: n + 1,
    'B': lambda n: 2 * n - 1,
    'C': lambda n: n + 1,
    'D': lambda n: 2 * (n - 1),
    'G2': 4,
    'F4': 9,
    'E6': 12,
    'E7': 18,
    'E8': 30,
}

LACING_NUMBER = {
    'A': 1,   # simply-laced
    'B': 2,   # long/short ratio = sqrt(2)
    'C': 2,   # long/short ratio = sqrt(2)
    'D': 1,   # simply-laced
    'G2': 3,  # long/short ratio = sqrt(3)
    'F4': 2,  # long/short ratio = sqrt(2)
    'E6': 1,  # simply-laced
    'E7': 1,  # simply-laced
    'E8': 1,  # simply-laced
}


def get_coxeter(family, rank_param):
    """Get Coxeter number h."""
    key = family if family in ('G2', 'F4', 'E6', 'E7', 'E8') else family[0]
    val = COXETER_NUMBERS[key]
    return val if isinstance(val, int) else val(rank_param)


def get_dual_coxeter(family, rank_param):
    """Get dual Coxeter number h_v."""
    key = family if family in ('G2', 'F4', 'E6', 'E7', 'E8') else family[0]
    val = DUAL_COXETER_NUMBERS[key]
    return val if isinstance(val, int) else val(rank_param)


def get_lacing(family):
    """Get lacing number (1 = simply-laced)."""
    key = family if family in ('G2', 'F4', 'E6', 'E7', 'E8') else family[0]
    return LACING_NUMBER[key]


# ============================================================
# Degree distribution analysis
# ============================================================

def degree_analysis(rank, roots, simple_roots):
    """Detailed degree distribution from root system."""
    adj = adjacency_from_roots(rank, roots, simple_roots)
    degrees = adj.sum(axis=1).astype(int)

    # Cartan degrees (first `rank` vertices)
    cartan_degs = degrees[:rank]
    # Root degrees (remaining vertices)
    root_degs = degrees[rank:]

    root_arr = [np.array(r, dtype=np.float64) for r in roots]
    simple_arr = [np.array(s, dtype=np.float64) for s in simple_roots]

    # For each root, count:
    #   - n_cartan: connections to Cartan generators
    #   - n_sum: connections via alpha+beta in Phi
    #   - n_neg: connection to -alpha (always 1 if -alpha in roots)
    root_set = set(tuple(round(x, 6) for x in r) for r in roots)

    n_cartan_per_root = []
    n_sum_per_root = []
    n_neg_per_root = []

    for j, alpha in enumerate(root_arr):
        nc = sum(1 for s in simple_arr if abs(np.dot(alpha, s)) > 1e-10)
        n_cartan_per_root.append(nc)

        ns = 0
        nn = 0
        for k, beta in enumerate(root_arr):
            if k == j:
                continue
            s = alpha + beta
            if np.linalg.norm(s) < 1e-10:
                nn = 1
            else:
                gamma = tuple(round(x, 6) for x in s)
                if gamma in root_set:
                    ns += 1
        n_sum_per_root.append(ns)
        n_neg_per_root.append(nn)

    return {
        'cartan_degs': cartan_degs,
        'root_degs': root_degs,
        'n_cartan_per_root': np.array(n_cartan_per_root),
        'n_sum_per_root': np.array(n_sum_per_root),
        'n_neg_per_root': np.array(n_neg_per_root),
        'deg_min': degrees.min(),
        'deg_max': degrees.max(),
        'deg_mean': degrees.mean(),
        'deg_std': degrees.std(),
        'all_degrees': degrees,
    }


# ============================================================
# Root string analysis
# ============================================================

def root_string_lengths(roots):
    """For each pair (alpha, beta), compute the alpha-string through beta.
    The alpha-string through beta is: beta - p*alpha, ..., beta, ..., beta + q*alpha
    where p,q >= 0 and p - q = <beta, alpha_v> (= 2(beta.alpha)/(alpha.alpha)).
    The length of the string is p + q + 1.

    Returns the distribution of string lengths.
    """
    root_arr = [np.array(r, dtype=np.float64) for r in roots]
    root_set = set(tuple(round(x, 6) for x in r) for r in roots)
    string_lens = []

    for alpha in root_arr:
        a_norm_sq = np.dot(alpha, alpha)
        if a_norm_sq < 1e-10:
            continue
        for beta in root_arr:
            if np.allclose(alpha, beta) or np.allclose(alpha, -beta):
                continue
            # Count q: max k >= 0 s.t. beta + k*alpha in Phi
            q = 0
            while True:
                test = beta + (q + 1) * alpha
                key = tuple(round(x, 6) for x in test)
                if key in root_set:
                    q += 1
                else:
                    break
            # Count p: max k >= 0 s.t. beta - k*alpha in Phi
            p = 0
            while True:
                test = beta - (p + 1) * alpha
                key = tuple(round(x, 6) for x in test)
                if key in root_set:
                    p += 1
                else:
                    break
            string_lens.append(p + q + 1)

    return np.array(string_lens)


# ============================================================
# Build complete analysis table
# ============================================================

def build_analysis():
    """Complete invariant and spectral analysis for all algebras."""
    catalog = OrderedDict()

    # A family: su(n+1), rank param = n (su index = n+1, Lie rank = n)
    for n in range(1, 12):
        name = f"su({n+1})"
        rank, roots, simple = roots_A(n)
        catalog[name] = {
            'family': 'A', 'lie_rank': n, 'rank_param': n,
            'root_data': (rank, roots, simple),
        }

    # B family
    for k in range(2, 9):
        name = f"so({2*k+1})"
        rank, roots, simple = roots_B(k)
        catalog[name] = {
            'family': 'B', 'lie_rank': k, 'rank_param': k,
            'root_data': (rank, roots, simple),
        }

    # C family
    for k in range(2, 9):
        name = f"sp({2*k})"
        rank, roots, simple = roots_C(k)
        catalog[name] = {
            'family': 'C', 'lie_rank': k, 'rank_param': k,
            'root_data': (rank, roots, simple),
        }

    # D family
    for k in range(3, 10):
        name = f"so({2*k})"
        rank, roots, simple = roots_D(k)
        catalog[name] = {
            'family': 'D', 'lie_rank': k, 'rank_param': k,
            'root_data': (rank, roots, simple),
        }

    # Exceptionals
    for exc_name, exc_fn, exc_fam in [
        ('G2', roots_G2, 'G2'), ('F4', roots_F4, 'F4'),
        ('E6', roots_E6, 'E6'), ('E7', roots_E7, 'E7'), ('E8', roots_E8, 'E8')
    ]:
        rank, roots, simple = exc_fn()
        catalog[exc_name] = {
            'family': exc_fam, 'lie_rank': rank, 'rank_param': rank,
            'root_data': (rank, roots, simple),
        }

    results = OrderedDict()
    print("=" * 100)
    print("  THEORETICAL ANALYSIS: Root System Invariants vs Spectral Properties")
    print("=" * 100)

    for name, info in catalog.items():
        rank, roots, simple = info['root_data']
        n_roots = len(roots)
        dim = rank + n_roots

        h = get_coxeter(info['family'], info['rank_param'])
        hv = get_dual_coxeter(info['family'], info['rank_param'])
        lac = get_lacing(info['family'])

        adj = adjacency_from_roots(rank, roots, simple)
        spec = spectral_analysis(adj)
        deg_info = degree_analysis(rank, roots, simple)

        # Key derived quantities
        n_positive = n_roots // 2
        edge_density = int(adj.sum() // 2) / (dim * (dim - 1) / 2)

        results[name] = {
            'family': info['family'],
            'lie_rank': info['lie_rank'],
            'dim': dim,
            'n_roots': n_roots,
            'n_positive': n_positive,
            'h': h,
            'hv': hv,
            'lacing': lac,
            'spectral': spec,
            'degree': deg_info,
            'edge_density': edge_density,
        }

        ram_str = "YES" if spec['is_ramanujan'] else "NO"
        print(f"  {name:>8}: dim={dim:3d} r={rank} |Phi|={n_roots:3d} "
              f"h={h:2d} hv={hv:2d} lac={lac} "
              f"deg=[{deg_info['deg_min']},{deg_info['deg_max']}] "
              f"ratio={spec['ratio']:.4f} Ram={ram_str}")

    return results


# ============================================================
# Theoretical predictions
# ============================================================

def theoretical_analysis(results):
    """Derive scaling laws and theoretical bounds."""
    print("\n" + "=" * 100)
    print("  SCALING ANALYSIS")
    print("=" * 100)

    names = list(results.keys())
    dims = np.array([results[n]['dim'] for n in names])
    ratios = np.array([results[n]['spectral']['ratio'] for n in names])
    ranks = np.array([results[n]['lie_rank'] for n in names])
    h_vals = np.array([results[n]['h'] for n in names])
    hv_vals = np.array([results[n]['hv'] for n in names])
    n_roots = np.array([results[n]['n_roots'] for n in names])
    deg_max = np.array([results[n]['degree']['deg_max'] for n in names])
    deg_mean = np.array([results[n]['degree']['deg_mean'] for n in names])
    densities = np.array([results[n]['edge_density'] for n in names])
    lam_max = np.array([results[n]['spectral']['lambda_max'] for n in names])
    max_nt = np.array([results[n]['spectral']['max_nontrivial'] for n in names])
    ram_bound = np.array([results[n]['spectral']['ramanujan_bound'] for n in names])
    is_ram = np.array([results[n]['spectral']['is_ramanujan'] for n in names])
    lacing = np.array([results[n]['lacing'] for n in names])
    gaps = np.array([results[n]['spectral']['spectral_gap'] for n in names])

    # ---- Fit: ratio vs dim ----
    valid = dims > 3
    log_d = np.log(dims[valid])
    log_r = np.log(ratios[valid])
    # Linear fit in log-log: log(ratio) = a * log(dim) + b
    coeffs_dim = np.polyfit(log_d, log_r, 1)
    print(f"\n  log(ratio) = {coeffs_dim[0]:.4f} * log(dim) + {coeffs_dim[1]:.4f}")
    print(f"  => ratio ~ dim^{coeffs_dim[0]:.3f} * {np.exp(coeffs_dim[1]):.4f}")

    # ---- Fit: ratio vs Coxeter number ----
    log_h = np.log(h_vals[valid])
    coeffs_h = np.polyfit(log_h, log_r, 1)
    print(f"\n  log(ratio) = {coeffs_h[0]:.4f} * log(h) + {coeffs_h[1]:.4f}")
    print(f"  => ratio ~ h^{coeffs_h[0]:.3f} * {np.exp(coeffs_h[1]):.4f}")

    # ---- Fit: ratio vs rank ----
    log_rk = np.log(ranks[valid])
    coeffs_rk = np.polyfit(log_rk, log_r, 1)
    print(f"\n  log(ratio) = {coeffs_rk[0]:.4f} * log(rank) + {coeffs_rk[1]:.4f}")
    print(f"  => ratio ~ rank^{coeffs_rk[0]:.3f} * {np.exp(coeffs_rk[1]):.4f}")

    # ---- Key relationship: max|lambda_nt| vs deg_max ----
    print(f"\n  Degree analysis:")
    print(f"    max|lambda_nt| / deg_max: "
          f"min={min(max_nt/deg_max):.4f}, max={max(max_nt/deg_max):.4f}")
    print(f"    max|lambda_nt| / deg_mean: "
          f"min={min(max_nt/deg_mean):.4f}, max={max(max_nt/deg_mean):.4f}")

    # ---- Theoretical prediction: for A_n ----
    print(f"\n  Family A analysis (su(n)):")
    a_mask = np.array([results[n]['family'] == 'A' for n in names])
    if a_mask.sum() > 2:
        a_dims = dims[a_mask]
        a_ratios = ratios[a_mask]
        a_ranks = ranks[a_mask]
        a_h = h_vals[a_mask]
        a_max_nt = max_nt[a_mask]
        a_ram = ram_bound[a_mask]
        a_deg_mean = deg_mean[a_mask]

        for i, n in enumerate(np.array(names)[a_mask]):
            r = results[n]
            print(f"    {n:>8}: r={r['lie_rank']:2d} h={r['h']:2d} "
                  f"dim={r['dim']:3d} deg_mean={r['degree']['deg_mean']:.1f} "
                  f"max_nt={r['spectral']['max_nontrivial']:.3f} "
                  f"2sq(q)={r['spectral']['ramanujan_bound']:.3f} "
                  f"ratio={r['spectral']['ratio']:.4f}")

        # For A_n: |Phi| = n(n+1), dim = (n+1)^2 - 1, h = n+1
        # Average degree ~ 2*n (from data)
        # Ramanujan bound ~ 2*sqrt(2n - 1) ~ 2*sqrt(2)*sqrt(n)
        # max|lambda_nt| seems to grow ~ n
        # So ratio ~ n / (2*sqrt(2)*sqrt(n)) = sqrt(n) / (2*sqrt(2))
        # => ratio ~ sqrt(rank) / C
        a_log_r = np.log(a_ratios)
        a_log_rk = np.log(a_ranks)
        c_a = np.polyfit(a_log_rk, a_log_r, 1)
        print(f"\n    Family A fit: ratio ~ rank^{c_a[0]:.3f} * {np.exp(c_a[1]):.4f}")
        print(f"    Prediction: ratio = 1 at rank ~ "
              f"{np.exp(-c_a[1]/c_a[0]):.1f} (i.e. su({int(np.exp(-c_a[1]/c_a[0]))+1}))")

    # ---- Theoretical insight: degree vs dim scaling ----
    print(f"\n  Degree scaling:")
    print(f"    deg_mean / dim: {(deg_mean/dims).mean():.4f} "
          f"(range {(deg_mean/dims).min():.4f} - {(deg_mean/dims).max():.4f})")
    print(f"    deg_mean ~ O(dim^alpha):")
    c_deg = np.polyfit(np.log(dims[valid]), np.log(deg_mean[valid]), 1)
    print(f"    alpha = {c_deg[0]:.4f}")

    # ---- The key formula ----
    # Ramanujan bound = 2*sqrt(q) where q+1 = avg_degree
    # For large algebras: avg_degree ~ C * dim^alpha
    # => ram_bound ~ 2*sqrt(C * dim^alpha) ~ dim^(alpha/2)
    # max|lambda_nt| grows faster than dim^(alpha/2)
    print(f"\n  Ram bound scaling: 2*sqrt(q) ~ dim^?:")
    c_rb = np.polyfit(np.log(dims[valid]), np.log(ram_bound[valid]), 1)
    print(f"    ram_bound ~ dim^{c_rb[0]:.4f}")
    print(f"\n  max|lambda_nt| scaling: ~ dim^?:")
    c_nt = np.polyfit(np.log(dims[valid]), np.log(max_nt[valid]), 1)
    print(f"    max_nt ~ dim^{c_nt[0]:.4f}")
    print(f"\n  => ratio ~ dim^({c_nt[0]:.3f} - {c_rb[0]:.3f}) = dim^{c_nt[0]-c_rb[0]:.3f}")

    # ---- Edge density transition ----
    print(f"\n  Edge density analysis:")
    for n in names:
        r = results[n]
        if r['spectral']['is_ramanujan']:
            status = " "
        else:
            status = "*"
        print(f"  {status} {n:>8}: density={r['edge_density']:.4f} "
              f"ratio={r['spectral']['ratio']:.4f}")

    # ---- Simply-laced vs non-simply-laced ----
    print(f"\n  Lacing effect (at similar dim):")
    sl_mask = lacing == 1
    nsl_mask = lacing > 1
    if sl_mask.sum() > 0 and nsl_mask.sum() > 0:
        print(f"    Simply-laced: mean ratio = {ratios[sl_mask].mean():.4f}")
        print(f"    Non-simply-laced: mean ratio = {ratios[nsl_mask].mean():.4f}")

    return {
        'power_dim': coeffs_dim[0],
        'power_h': coeffs_h[0],
        'power_rank': coeffs_rk[0],
        'deg_scaling': c_deg[0],
        'ram_bound_scaling': c_rb[0],
        'max_nt_scaling': c_nt[0],
        'ratio_scaling': c_nt[0] - c_rb[0],
    }


# ============================================================
# Visualization
# ============================================================

def plot_theory(results, fits, save=False):
    """Multi-panel theoretical analysis visualization."""
    names = list(results.keys())
    dims = [results[n]['dim'] for n in names]
    ratios = [results[n]['spectral']['ratio'] for n in names]
    ranks = [results[n]['lie_rank'] for n in names]
    h_vals = [results[n]['h'] for n in names]
    max_nt = [results[n]['spectral']['max_nontrivial'] for n in names]
    ram_bound = [results[n]['spectral']['ramanujan_bound'] for n in names]
    deg_mean = [results[n]['degree']['deg_mean'] for n in names]
    deg_max = [results[n]['degree']['deg_max'] for n in names]
    deg_min = [results[n]['degree']['deg_min'] for n in names]
    densities = [results[n]['edge_density'] for n in names]
    families = [results[n]['family'] for n in names]
    is_ram = [results[n]['spectral']['is_ramanujan'] for n in names]

    fam_colors = {'A': 'royalblue', 'B': 'crimson', 'C': 'forestgreen',
                  'D': 'darkorange', 'G2': 'purple', 'F4': 'brown',
                  'E6': 'gold', 'E7': 'goldenrod', 'E8': 'darkgoldenrod'}
    fam_markers = {'A': 'o', 'B': 's', 'C': '^', 'D': 'v',
                   'G2': 'D', 'F4': 'p', 'E6': '*', 'E7': '*', 'E8': '*'}

    fig, axes = plt.subplots(2, 3, figsize=(22, 13))
    fig.suptitle('Theoretical Analysis: Root Invariants vs Spectral Properties\n'
                 'Adjoint Commutation Graph — Ramanujan Transition',
                 fontsize=14, fontweight='bold')

    def scatter(ax, xs, ys, annotate_big=True):
        for i, n in enumerate(names):
            f = families[i]
            c = fam_colors.get(f, 'grey')
            m = fam_markers.get(f, 'o')
            sz = 150 if f.startswith('E') or f in ('G2', 'F4') else 60
            ec = 'red' if not is_ram[i] else 'black'
            lw = 2 if not is_ram[i] else 0.5
            ax.scatter(xs[i], ys[i], c=c, marker=m, s=sz,
                       edgecolors=ec, linewidths=lw, zorder=3)
            if annotate_big and (dims[i] >= 50 or f in ('G2', 'F4', 'E6', 'E7', 'E8')):
                ax.annotate(n, (xs[i], ys[i]), fontsize=6,
                            ha='left', va='bottom', xytext=(4, 4),
                            textcoords='offset points')

    # Panel 1: log-log ratio vs dim with fit
    ax = axes[0, 0]
    scatter(ax, dims, ratios)
    d_fit = np.linspace(3, 260, 200)
    r_fit = np.exp(fits['power_dim'] * np.log(d_fit) +
                   np.polyfit(np.log([d for d in dims if d > 3]),
                              np.log([r for d, r in zip(dims, ratios) if d > 3]), 1)[1])
    ax.plot(d_fit, r_fit, 'k--', alpha=0.5,
            label=f'fit: dim$^{{{fits["power_dim"]:.2f}}}$')
    ax.axhline(1.0, color='red', linestyle='-', alpha=0.3)
    ax.set_xlabel('dim(algebra)')
    ax.set_ylabel('Ramanujan ratio')
    ax.set_title(f'Ratio ~ dim$^{{{fits["power_dim"]:.2f}}}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 2: ratio vs Coxeter number
    ax = axes[0, 1]
    scatter(ax, h_vals, ratios)
    ax.axhline(1.0, color='red', linestyle='-', alpha=0.3)
    ax.set_xlabel('Coxeter number h')
    ax.set_ylabel('Ramanujan ratio')
    ax.set_title(f'Ratio ~ h$^{{{fits["power_h"]:.2f}}}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)

    # Panel 3: max|lambda_nt| and 2sqrt(q) vs dim
    ax = axes[0, 2]
    ax.scatter(dims, max_nt, c='crimson', marker='x', s=40, label='max|$\\lambda_{nt}$|', zorder=3)
    ax.scatter(dims, ram_bound, c='blue', marker='+', s=40, label='2$\\sqrt{q}$', zorder=3)
    d_fit = np.linspace(3, 260, 200)
    ax.plot(d_fit, np.exp(fits['max_nt_scaling'] * np.log(d_fit) +
            np.polyfit(np.log([d for d in dims if d > 3]),
                       np.log([m for d, m in zip(dims, max_nt) if d > 3]), 1)[1]),
            'r--', alpha=0.5, label=f'dim$^{{{fits["max_nt_scaling"]:.2f}}}$')
    ax.plot(d_fit, np.exp(fits['ram_bound_scaling'] * np.log(d_fit) +
            np.polyfit(np.log([d for d in dims if d > 3]),
                       np.log([r for d, r in zip(dims, ram_bound) if d > 3]), 1)[1]),
            'b--', alpha=0.5, label=f'dim$^{{{fits["ram_bound_scaling"]:.2f}}}$')
    ax.set_xlabel('dim(algebra)')
    ax.set_ylabel('Spectral value')
    ax.set_title('Why the ratio grows: different exponents')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 4: Degree distribution (max, mean, min vs dim)
    ax = axes[1, 0]
    ax.scatter(dims, deg_max, c='red', marker='^', s=30, label='deg_max', alpha=0.7)
    ax.scatter(dims, deg_mean, c='blue', marker='o', s=30, label='deg_mean', alpha=0.7)
    ax.scatter(dims, deg_min, c='green', marker='v', s=30, label='deg_min', alpha=0.7)
    ax.set_xlabel('dim(algebra)')
    ax.set_ylabel('Degree')
    ax.set_title('Degree distribution vs dimension')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 5: Edge density vs dim
    ax = axes[1, 1]
    scatter(ax, dims, densities)
    ax.set_xlabel('dim(algebra)')
    ax.set_ylabel('Edge density |E| / (n choose 2)')
    ax.set_title('Sparsification with dimension')
    ax.grid(True, alpha=0.3)

    # Panel 6: Family-level ratio scaling with rank
    ax = axes[1, 2]
    fam_data = OrderedDict()
    for n, r in results.items():
        f = r['family']
        f_key = f if f in ('A', 'B', 'C', 'D') else f
        if f_key not in fam_data:
            fam_data[f_key] = ([], [])
        fam_data[f_key][0].append(r['lie_rank'])
        fam_data[f_key][1].append(r['spectral']['ratio'])
    for f, (rks, rs) in fam_data.items():
        c = fam_colors.get(f, 'grey')
        m = fam_markers.get(f, 'o')
        if len(rks) > 1:
            ax.plot(rks, rs, marker=m, color=c, linestyle='-', alpha=0.8,
                    label=f'Family {f}', markersize=8)
        else:
            ax.scatter(rks, rs, c=c, marker=m, s=100, label=f, zorder=3)
    ax.axhline(1.0, color='red', linestyle='--', alpha=0.3)
    ax.set_xlabel('Lie rank')
    ax.set_ylabel('Ramanujan ratio')
    ax.set_title('Scaling with rank per family')
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    if save:
        out = os.path.join(os.path.dirname(__file__), 'ramanujan_theory_analysis.png')
        fig.savefig(out, dpi=150, bbox_inches='tight')
        print(f"\n  Saved to {out}")
    else:
        plt.show()


# ============================================================
# Analytical bounds derivation
# ============================================================

def analytical_bounds(results):
    """Derive approximate analytical bounds for the Ramanujan ratio."""
    print("\n" + "=" * 100)
    print("  ANALYTICAL BOUNDS DERIVATION")
    print("=" * 100)

    # For A_n = su(n+1):
    #   rank = n, |Phi| = n(n+1), dim = n^2 + 2n = (n+1)^2 - 1
    #   Each root e_i - e_j connects to:
    #     - Cartan H_k if (e_i-e_j) . (e_k-e_{k+1}) != 0
    #       This is nonzero when k=i or k=i-1 or k=j or k=j-1
    #       For generic i,j: connects to min(4, ...) Cartan elements
    #     - Other roots e_k-e_l if (e_i-e_j)+(e_k-e_l) is a root
    #       i.e. if j=k (giving e_i-e_l, a root if i!=l)
    #       or i=l (giving e_k-e_j, a root if k!=j)
    #       Count: for fixed e_i-e_j, number of e_k-e_l with k=j: (n+1-1)=n choices for l!=j, minus l=i
    #       Plus l=i: n choices for k!=i, minus k=j
    #       Total root neighbors from sums: 2(n-1)
    #     - Plus the negative root -(e_i-e_j) = e_j-e_i: +1
    #     Total degree of E_{e_i-e_j} ~ 2(n-1) + 1 + (Cartan connections)
    #
    #   For Cartan H_k (= e_k - e_{k+1}):
    #     connects to E_alpha iff alpha . (e_k-e_{k+1}) != 0
    #     alpha = e_i - e_j: dot = delta_{ik} - delta_{jk} - delta_{i,k+1} + delta_{j,k+1}
    #     nonzero when i or j equals k or k+1
    #     Count: roots with i=k: n choices for j!=k, same for j=k, i=k+1, j=k+1
    #     ~ 4n - corrections

    print("\n  A_n = su(n+1) analytical:")
    print("  For E_alpha (alpha = e_i - e_j):")
    print("    Root-sum neighbors: (n-1) from j=k + (n-1) from l=i = 2(n-1)")
    print("    Negative root: 1")
    print("    Cartan neighbors: varies, but typically 2-4")
    print("    => deg(E_alpha) ~ 2n - 1 + Cartan_connections")
    print()
    print("  For H_k (= e_k - e_{k+1}):")
    print("    Connects to e_i-e_j when {i,j} intersects {k,k+1}")
    print("    Count ~ 4n - 4  (for interior k)")
    print()

    # Verify with data
    for n in range(1, 12):
        name = f"su({n+1})"
        if name not in results:
            continue
        r = results[name]
        d = r['degree']
        predicted_root_deg = 2 * n - 1 + 3  # approx Cartan contribution
        actual_mean_root = d['root_degs'].mean()
        print(f"    su({n+1:2d}): r={n:2d}, mean_root_deg={actual_mean_root:.1f}, "
              f"Cartan_degs={list(d['cartan_degs'])}, "
              f"n_sum_mean={d['n_sum_per_root'].mean():.1f}")

    # For non-regular graph: the issue is degree variance
    # avg_degree ~ q+1 is used in Ramanujan bound 2*sqrt(q)
    # But max nontrivial eigenvalue can be dominated by high-degree vertices
    # The Alon-Boppana bound for d-regular: lambda_2 >= 2*sqrt(d-1) - o(1)
    # For non-regular: Greenberg (1995) generalization

    print("\n  KEY THEORETICAL INSIGHT:")
    print("  The adjoint graph is NON-regular (Cartan vs root vertices have")
    print("  very different degrees). The Ramanujan bound 2*sqrt(q) uses the")
    print("  AVERAGE degree, but the actual max nontrivial eigenvalue is")
    print("  controlled by the MAXIMUM degree (Perron-Frobenius effect).")
    print()

    # Check: max_nt vs sqrt(deg_max)
    print("  Comparison: max|lambda_nt| vs sqrt(deg_max) and sqrt(deg_mean):")
    for n in names_sorted_by_dim(results):
        r = results[n]
        s = r['spectral']
        d = r['degree']
        print(f"    {n:>8}: max_nt={s['max_nontrivial']:7.3f}  "
              f"sqrt(dmax)={np.sqrt(d['deg_max']):6.3f}  "
              f"sqrt(dmean)={np.sqrt(d['deg_mean']):6.3f}  "
              f"2sqrt(q)={s['ramanujan_bound']:6.3f}  "
              f"dmax/dmean={d['deg_max']/d['deg_mean']:.3f}")


def names_sorted_by_dim(results):
    return sorted(results.keys(), key=lambda n: results[n]['dim'])


# ============================================================
# Entry point
# ============================================================

if __name__ == '__main__':
    results = build_analysis()
    fits = theoretical_analysis(results)
    analytical_bounds(results)
    plot_theory(results, fits, save=True)
