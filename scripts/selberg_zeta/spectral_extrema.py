"""
P10: Spectral Extrema Characterization
=======================================
P9 showed max|lambda_nt| = h-1 holds for only 4/36 classical algebras.
This script investigates: what algebraic invariant actually controls
the dominant nontrivial eigenvalue, and what distinguishes algebras
where the h-1 relation holds exactly?

Analyses:
  1. All 37 simple algebras + extended classical families: max|lambda_nt| data
  2. Regression against multiple candidate predictors
  3. Simply-laced vs non-simply-laced dichotomy
  4. Eigenvector analysis of the dominant mode
  5. Exact characterization of the h-1 locus
  6. Universal formula derivation
  7. 9-panel visualization
"""
import sys, os, time
import numpy as np
from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.dirname(__file__))
from ramanujan_extended import (
    roots_A, roots_B, roots_C, roots_D,
    roots_G2, roots_F4, roots_E6, roots_E7, roots_E8,
    adjacency_from_roots, spectral_analysis,
)

# ======================================================================
# Data: all simple Lie algebras with algebraic metadata
# ======================================================================

ALGEBRAS = OrderedDict()

# Classical families — extended ranges
for n in range(1, 14):  # A_n, n=1..13
    name = f"A_{n}"
    ALGEBRAS[name] = {
        'root_fn': roots_A, 'rank_arg': n,
        'family': 'A', 'simply_laced': True,
        'h': n + 1, 'h_dual': n + 1,
        'rank': n, 'n_roots': n * (n + 1),
    }

for n in range(2, 10):  # B_n, n=2..9
    name = f"B_{n}"
    ALGEBRAS[name] = {
        'root_fn': roots_B, 'rank_arg': n,
        'family': 'B', 'simply_laced': False,
        'h': 2 * n, 'h_dual': 2 * n - 1,
        'rank': n, 'n_roots': 2 * n * (2 * n - 1) // 1,  # will fix below
    }

for n in range(2, 10):  # C_n, n=2..9
    name = f"C_{n}"
    ALGEBRAS[name] = {
        'root_fn': roots_C, 'rank_arg': n,
        'family': 'C', 'simply_laced': False,
        'h': 2 * n, 'h_dual': n + 1,
        'rank': n, 'n_roots': 2 * n * (2 * n - 1) // 1,
    }

for n in range(3, 10):  # D_n, n=3..9
    name = f"D_{n}"
    ALGEBRAS[name] = {
        'root_fn': roots_D, 'rank_arg': n,
        'family': 'D', 'simply_laced': True,
        'h': 2 * (n - 1), 'h_dual': 2 * (n - 1),
        'rank': n, 'n_roots': 2 * n * (n - 1),
    }

# Exceptional algebras
ALGEBRAS['G_2'] = {
    'root_fn': roots_G2, 'rank_arg': None,
    'family': 'G', 'simply_laced': False,
    'h': 6, 'h_dual': 4, 'rank': 2, 'n_roots': 12,
}
ALGEBRAS['F_4'] = {
    'root_fn': roots_F4, 'rank_arg': None,
    'family': 'F', 'simply_laced': False,
    'h': 12, 'h_dual': 9, 'rank': 4, 'n_roots': 48,
}
ALGEBRAS['E_6'] = {
    'root_fn': roots_E6, 'rank_arg': None,
    'family': 'E', 'simply_laced': True,
    'h': 12, 'h_dual': 12, 'rank': 6, 'n_roots': 72,
}
ALGEBRAS['E_7'] = {
    'root_fn': roots_E7, 'rank_arg': None,
    'family': 'E', 'simply_laced': True,
    'h': 18, 'h_dual': 18, 'rank': 7, 'n_roots': 126,
}
ALGEBRAS['E_8'] = {
    'root_fn': roots_E8, 'rank_arg': None,
    'family': 'E', 'simply_laced': True,
    'h': 30, 'h_dual': 30, 'rank': 8, 'n_roots': 240,
}


def compute_algebra(name, info):
    """Compute full spectral data for one algebra."""
    if info['rank_arg'] is not None:
        rank, roots, simple = info['root_fn'](info['rank_arg'])
    else:
        rank, roots, simple = info['root_fn']()

    adj = adjacency_from_roots(rank, roots, simple)
    sa = spectral_analysis(adj)
    dim = adj.shape[0]

    # Degree statistics
    degs = np.sum(adj, axis=1).astype(float)
    d_mean = np.mean(degs)
    d_max = np.max(degs)
    d_min = np.min(degs)
    cv = np.std(degs) / d_mean if d_mean > 0 else 0

    # Eigenvector of dominant nontrivial eigenvalue
    eigenvalues, eigenvectors = np.linalg.eigh(adj.astype(float))
    idx_sorted = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx_sorted]
    eigenvectors = eigenvectors[:, idx_sorted]

    # lambda_max is eigenvalues[0] (trivial/Perron-Frobenius)
    # max_nontrivial is max(|eigenvalues[1]|, |eigenvalues[-1]|)
    max_nt = sa['max_nontrivial']
    h = info['h']
    h_dual = info['h_dual']

    # Which eigenvector corresponds to max_nontrivial?
    if abs(eigenvalues[1]) >= abs(eigenvalues[-1]):
        dom_idx = 1
        dom_eig = eigenvalues[1]
    else:
        # Find the most negative eigenvalue index
        dom_idx = len(eigenvalues) - 1
        dom_eig = eigenvalues[-1]

    dom_vec = eigenvectors[:, dom_idx]

    # Participation ratio of dominant eigenvector: PR = 1/(sum(v_i^4)) / dim
    # Measures how delocalized the eigenvector is
    pr = 1.0 / (dim * np.sum(dom_vec**4)) if np.sum(dom_vec**4) > 0 else 0

    # Cartan vs root components
    n_cartan = rank
    n_root = len(roots)
    cartan_weight = np.sum(dom_vec[:n_cartan]**2)
    root_weight = np.sum(dom_vec[n_cartan:]**2)

    return {
        'name': name,
        'family': info['family'],
        'simply_laced': info['simply_laced'],
        'dim': dim,
        'rank': rank,
        'n_roots': len(roots),
        'h': h,
        'h_dual': h_dual,
        'd_mean': d_mean,
        'd_max': int(d_max),
        'd_min': int(d_min),
        'cv_degree': cv,
        'ratio': sa['ratio'],
        'is_ramanujan': bool(sa['is_ramanujan']),
        'max_nontrivial': max_nt,
        'lambda_max': sa['lambda_max'],
        'spectral_gap': sa['spectral_gap'],
        'h_minus_1': h - 1,
        'deviation': max_nt - (h - 1),
        'rel_deviation': (max_nt - (h - 1)) / (h - 1) if h > 1 else 0,
        'max_nt_over_h1': max_nt / (h - 1) if h > 1 else float('inf'),
        'eigenvalues': eigenvalues,
        'dom_eigenvector': dom_vec,
        'participation_ratio': pr,
        'cartan_weight': cartan_weight,
        'root_weight': root_weight,
    }


def regression_analysis(results):
    """Test multiple predictors for max|lambda_nt|."""
    from scipy.stats import pearsonr

    y = np.array([r['max_nontrivial'] for r in results])

    predictors = {
        'h - 1': np.array([r['h'] - 1 for r in results]),
        'h_dual - 1': np.array([r['h_dual'] - 1 for r in results]),
        'rank': np.array([r['rank'] for r in results], dtype=float),
        'sqrt(dim)': np.sqrt([r['dim'] for r in results]),
        'd_mean': np.array([r['d_mean'] for r in results]),
        'sqrt(d_mean)': np.sqrt([r['d_mean'] for r in results]),
        'd_max': np.array([r['d_max'] for r in results], dtype=float),
        'n_roots/rank': np.array([r['n_roots'] / r['rank'] for r in results]),
        'dim/h': np.array([r['dim'] / r['h'] for r in results]),
        'rank*(h-1)': np.array([r['rank'] * (r['h'] - 1) for r in results], dtype=float),
    }

    print("\n=== Phase 2: Regression Analysis ===")
    print(f"{'Predictor':>20} {'Pearson r':>10} {'p-value':>12} {'MSE':>10}")
    print("-" * 55)

    best_r = 0
    best_name = ''
    for pname, x in sorted(predictors.items(), key=lambda kv: -abs(pearsonr(kv[1], y)[0])):
        r, p = pearsonr(x, y)
        # Linear fit
        coeffs = np.polyfit(x, y, 1)
        predicted = np.polyval(coeffs, x)
        mse = np.mean((y - predicted)**2)
        print(f"{pname:>20} {r:10.4f} {p:12.2e} {mse:10.4f}")
        if abs(r) > best_r:
            best_r = abs(r)
            best_name = pname

    print(f"\nBest linear predictor: {best_name} (r = {best_r:.4f})")

    # Power-law fits: max_nt = a * x^b
    print("\n=== Phase 2b: Power-law fits ===")
    print(f"{'Predictor':>20} {'a':>8} {'b':>8} {'R2':>8}")
    print("-" * 50)

    for pname in ['h - 1', 'h_dual - 1', 'd_mean', 'dim/h']:
        x = predictors[pname]
        mask = (x > 0) & (y > 0)
        if np.sum(mask) < 3:
            continue
        log_x = np.log(x[mask])
        log_y = np.log(y[mask])
        b, log_a = np.polyfit(log_x, log_y, 1)
        a = np.exp(log_a)
        predicted = a * x[mask]**b
        ss_res = np.sum((y[mask] - predicted)**2)
        ss_tot = np.sum((y[mask] - np.mean(y[mask]))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        print(f"{pname:>20} {a:8.4f} {b:8.4f} {r2:8.4f}")

    return predictors


def laced_analysis(results):
    """Simply-laced vs non-simply-laced comparison."""
    print("\n=== Phase 3: Simply-laced vs Non-simply-laced ===")

    sl = [r for r in results if r['simply_laced']]
    nsl = [r for r in results if not r['simply_laced']]

    print(f"\nSimply-laced ({len(sl)} algebras): A, D, E families")
    sl_devs = [r['rel_deviation'] for r in sl]
    print(f"  max|lambda_nt|/(h-1): mean={np.mean([r['max_nt_over_h1'] for r in sl if r['h'] > 1]):.3f}, "
          f"std={np.std([r['max_nt_over_h1'] for r in sl if r['h'] > 1]):.3f}")
    exact = [r['name'] for r in sl if abs(r['deviation']) < 0.01]
    print(f"  Exact h-1 matches: {exact}")
    near = [r['name'] for r in sl if abs(r['deviation']) < 0.5 and abs(r['deviation']) >= 0.01]
    print(f"  Near h-1 (|dev| < 0.5): {near}")

    print(f"\nNon-simply-laced ({len(nsl)} algebras): B, C, G, F families")
    print(f"  max|lambda_nt|/(h-1): mean={np.mean([r['max_nt_over_h1'] for r in nsl if r['h'] > 1]):.3f}, "
          f"std={np.std([r['max_nt_over_h1'] for r in nsl if r['h'] > 1]):.3f}")
    exact_nsl = [r['name'] for r in nsl if abs(r['deviation']) < 0.01]
    print(f"  Exact h-1 matches: {exact_nsl}")

    # Test: do simply-laced algebras follow max_nt = h-1 better?
    from scipy.stats import mannwhitneyu
    sl_ratios = [r['max_nt_over_h1'] for r in sl if r['h'] > 1]
    nsl_ratios = [r['max_nt_over_h1'] for r in nsl if r['h'] > 1]
    if len(sl_ratios) >= 3 and len(nsl_ratios) >= 3:
        U, p = mannwhitneyu(sl_ratios, nsl_ratios, alternative='two-sided')
        print(f"\n  Mann-Whitney U test (simply vs non-simply laced): p = {p:.4e}")


def eigenvector_analysis(results):
    """Analyze the dominant eigenvector structure."""
    print("\n=== Phase 4: Eigenvector Analysis ===")
    print(f"{'Name':>8} {'max_nt':>8} {'PR':>6} {'Cartan%':>8} {'Root%':>8} {'d_max':>6} {'exact':>6}")
    print("-" * 60)

    for r in sorted(results, key=lambda x: x['h']):
        exact = '***' if abs(r['deviation']) < 0.01 else ''
        print(f"{r['name']:>8} {r['max_nontrivial']:8.3f} {r['participation_ratio']:6.3f} "
              f"{r['cartan_weight']*100:7.1f}% {r['root_weight']*100:7.1f}% "
              f"{r['d_max']:6d} {exact:>6}")


def h1_locus_analysis(results):
    """Deep analysis of algebras where max|lambda_nt| = h-1 exactly."""
    print("\n=== Phase 5: The h-1 Locus ===")

    exact = [r for r in results if abs(r['deviation']) < 0.05]
    print(f"\nAlgebras with |max_nt - (h-1)| < 0.05 ({len(exact)}):")
    for r in exact:
        print(f"  {r['name']:>8}: h={r['h']}, max_nt={r['max_nontrivial']:.6f}, "
              f"h-1={r['h_minus_1']}, dev={r['deviation']:.6f}, "
              f"family={r['family']}, SL={r['simply_laced']}, "
              f"d_max={r['d_max']}, d_mean={r['d_mean']:.2f}")

    # What do exact-match algebras share?
    if exact:
        print("\n  Shared properties of h-1 locus:")
        all_sl = all(r['simply_laced'] for r in exact)
        print(f"    All simply-laced: {all_sl}")
        families = set(r['family'] for r in exact)
        print(f"    Families represented: {families}")
        h_vals = [r['h'] for r in exact]
        print(f"    Coxeter numbers: {h_vals}")

        # Check: does d_max = dim - 1 for any?
        for r in exact:
            print(f"    {r['name']}: d_max/dim = {r['d_max']}/{r['dim']} = {r['d_max']/r['dim']:.3f}, "
                  f"d_max/(h-1+rank) = {r['d_max']/(r['h']-1+r['rank']):.3f}")

    # Near-misses
    near = [r for r in results if 0.05 <= abs(r['deviation']) < 1.0]
    if near:
        print(f"\nNear-miss algebras (0.05 <= |dev| < 1.0): {len(near)}")
        for r in sorted(near, key=lambda x: abs(x['deviation'])):
            print(f"  {r['name']:>8}: dev={r['deviation']:+.3f}, family={r['family']}, "
                  f"SL={r['simply_laced']}")


def universal_formula(results):
    """Derive the best universal formula for max|lambda_nt|."""
    print("\n=== Phase 6: Universal Formula ===")
    from scipy.optimize import curve_fit

    # Filter: exclude A_1 (trivial, dim=3)
    data = [r for r in results if r['dim'] > 3]

    h = np.array([r['h'] for r in data], dtype=float)
    rank = np.array([r['rank'] for r in data], dtype=float)
    dim = np.array([r['dim'] for r in data], dtype=float)
    d_mean = np.array([r['d_mean'] for r in data])
    d_max = np.array([r['d_max'] for r in data], dtype=float)
    n_roots = np.array([r['n_roots'] for r in data], dtype=float)
    y = np.array([r['max_nontrivial'] for r in data])
    sl = np.array([1.0 if r['simply_laced'] else 0.0 for r in data])

    # Candidate formulas
    formulas = {}

    # 1. max_nt = a * h^b
    try:
        def f1(x, a, b): return a * x**b
        popt, _ = curve_fit(f1, h, y, p0=[1, 1])
        pred = f1(h, *popt)
        mse = np.mean((y - pred)**2)
        formulas['a*h^b'] = (popt, mse, f"max_nt = {popt[0]:.4f} * h^{popt[1]:.4f}")
    except Exception:
        pass

    # 2. max_nt = a * d_mean^b
    try:
        popt, _ = curve_fit(f1, d_mean, y, p0=[1, 0.5])
        pred = f1(d_mean, *popt)
        mse = np.mean((y - pred)**2)
        formulas['a*d_mean^b'] = (popt, mse, f"max_nt = {popt[0]:.4f} * d_mean^{popt[1]:.4f}")
    except Exception:
        pass

    # 3. max_nt = a * d_max^b
    try:
        popt, _ = curve_fit(f1, d_max, y, p0=[1, 0.5])
        pred = f1(d_max, *popt)
        mse = np.mean((y - pred)**2)
        formulas['a*d_max^b'] = (popt, mse, f"max_nt = {popt[0]:.4f} * d_max^{popt[1]:.4f}")
    except Exception:
        pass

    # 4. max_nt = a * (n_roots/rank)^b
    try:
        x4 = n_roots / rank
        popt, _ = curve_fit(f1, x4, y, p0=[1, 0.5])
        pred = f1(x4, *popt)
        mse = np.mean((y - pred)**2)
        formulas['a*(n_roots/rank)^b'] = (popt, mse,
            f"max_nt = {popt[0]:.4f} * (n_roots/rank)^{popt[1]:.4f}")
    except Exception:
        pass

    # 5. Multivariate: max_nt = a * h^b * rank^c
    try:
        def f5(X, a, b, c):
            return a * X[0]**b * X[1]**c
        X5 = np.vstack([h, rank])
        popt, _ = curve_fit(f5, X5, y, p0=[1, 1, 0])
        pred = f5(X5, *popt)
        mse = np.mean((y - pred)**2)
        formulas['a*h^b*rank^c'] = (popt, mse,
            f"max_nt = {popt[0]:.4f} * h^{popt[1]:.4f} * rank^{popt[2]:.4f}")
    except Exception:
        pass

    # 6. max_nt = a * sqrt(d_mean) + b
    try:
        def f6(x, a, b): return a * np.sqrt(x) + b
        popt, _ = curve_fit(f6, d_mean, y, p0=[1, 0])
        pred = f6(d_mean, *popt)
        mse = np.mean((y - pred)**2)
        formulas['a*sqrt(d_mean)+b'] = (popt, mse,
            f"max_nt = {popt[0]:.4f} * sqrt(d_mean) + {popt[1]:.4f}")
    except Exception:
        pass

    # 7. max_nt = a * d_max + b
    try:
        def f7(x, a, b): return a * x + b
        popt, _ = curve_fit(f7, d_max, y, p0=[0.5, 0])
        pred = f7(d_max, *popt)
        mse = np.mean((y - pred)**2)
        formulas['a*d_max+b'] = (popt, mse,
            f"max_nt = {popt[0]:.4f} * d_max + {popt[1]:.4f}")
    except Exception:
        pass

    print(f"\n{'Formula':>25} {'MSE':>10} {'Expression':>50}")
    print("-" * 90)
    for fname, (popt, mse, expr) in sorted(formulas.items(), key=lambda kv: kv[1][1]):
        print(f"{fname:>25} {mse:10.4f} {expr:>50}")

    return formulas, data


def create_figure(results, formulas, output_path):
    """9-panel visualization for P10."""
    fig = plt.figure(figsize=(20, 16))
    fig.suptitle("P10: Spectral Extrema Characterization", fontsize=16, fontweight='bold')
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.30)

    # Color map by family
    family_colors = {
        'A': '#1565C0', 'B': '#EF5350', 'C': '#66BB6A', 'D': '#AB47BC',
        'G': '#FF9800', 'F': '#795548', 'E': '#00BCD4',
    }
    family_markers = {
        'A': 'o', 'B': 's', 'C': '^', 'D': 'D',
        'G': 'p', 'F': 'h', 'E': '*',
    }

    def scatter_by_family(ax, x_fn, y_fn, **kwargs):
        for fam in ['A', 'B', 'C', 'D', 'G', 'F', 'E']:
            subset = [r for r in results if r['family'] == fam]
            if not subset:
                continue
            xx = [x_fn(r) for r in subset]
            yy = [y_fn(r) for r in subset]
            ax.scatter(xx, yy, c=family_colors[fam], marker=family_markers[fam],
                       s=60, label=fam, zorder=3, edgecolors='k', linewidth=0.3, **kwargs)

    # Panel 1: max|lambda_nt| vs h-1 (the conjecture)
    ax1 = fig.add_subplot(gs[0, 0])
    scatter_by_family(ax1, lambda r: r['h_minus_1'], lambda r: r['max_nontrivial'])
    hrange = np.linspace(0, 30, 100)
    ax1.plot(hrange, hrange, 'k--', alpha=0.5, label='y = h-1')
    ax1.set_xlabel('h - 1')
    ax1.set_ylabel(r'max$|\lambda_{nt}|$')
    ax1.set_title(r'max$|\lambda_{nt}|$ vs $h-1$')
    ax1.legend(fontsize=7, ncol=2)

    # Panel 2: max|lambda_nt| vs d_max
    ax2 = fig.add_subplot(gs[0, 1])
    scatter_by_family(ax2, lambda r: r['d_max'], lambda r: r['max_nontrivial'])
    drange = np.array([r['d_max'] for r in results])
    if 'a*d_max+b' in formulas:
        popt = formulas['a*d_max+b'][0]
        xfit = np.linspace(min(drange), max(drange), 100)
        ax2.plot(xfit, popt[0] * xfit + popt[1], 'r--', alpha=0.7,
                 label=f'{popt[0]:.3f}*d_max + {popt[1]:.2f}')
    ax2.set_xlabel(r'$d_{max}$')
    ax2.set_ylabel(r'max$|\lambda_{nt}|$')
    ax2.set_title(r'max$|\lambda_{nt}|$ vs $d_{max}$')
    ax2.legend(fontsize=7)

    # Panel 3: max|lambda_nt| vs d_mean
    ax3 = fig.add_subplot(gs[0, 2])
    scatter_by_family(ax3, lambda r: r['d_mean'], lambda r: r['max_nontrivial'])
    if 'a*d_mean^b' in formulas:
        popt = formulas['a*d_mean^b'][0]
        dfit = np.linspace(1, max(r['d_mean'] for r in results) * 1.05, 100)
        ax3.plot(dfit, popt[0] * dfit**popt[1], 'r--', alpha=0.7,
                 label=f'{popt[0]:.3f}$\\cdot d^{{{popt[1]:.3f}}}$')
    ax3.set_xlabel(r'$\langle d \rangle$')
    ax3.set_ylabel(r'max$|\lambda_{nt}|$')
    ax3.set_title(r'max$|\lambda_{nt}|$ vs mean degree')
    ax3.legend(fontsize=7)

    # Panel 4: Deviation from h-1 by family
    ax4 = fig.add_subplot(gs[1, 0])
    scatter_by_family(ax4, lambda r: r['h'], lambda r: r['deviation'])
    ax4.axhline(0, color='k', ls='--', alpha=0.4)
    ax4.set_xlabel('Coxeter number h')
    ax4.set_ylabel(r'max$|\lambda_{nt}|$ - (h-1)')
    ax4.set_title('Deviation from h-1 conjecture')

    # Panel 5: Ratio max_nt/(h-1) — simply-laced vs non-simply-laced
    ax5 = fig.add_subplot(gs[1, 1])
    sl_data = [r for r in results if r['simply_laced'] and r['h'] > 1]
    nsl_data = [r for r in results if not r['simply_laced'] and r['h'] > 1]
    ax5.scatter([r['h'] for r in sl_data],
                [r['max_nt_over_h1'] for r in sl_data],
                c='#1565C0', marker='o', s=60, label='Simply-laced', edgecolors='k', linewidth=0.3)
    ax5.scatter([r['h'] for r in nsl_data],
                [r['max_nt_over_h1'] for r in nsl_data],
                c='#EF5350', marker='s', s=60, label='Non-simply-laced', edgecolors='k', linewidth=0.3)
    ax5.axhline(1, color='k', ls='--', alpha=0.4, label='h-1 exact')
    ax5.set_xlabel('Coxeter number h')
    ax5.set_ylabel(r'max$|\lambda_{nt}|$ / (h-1)')
    ax5.set_title('Simply-laced vs non-simply-laced')
    ax5.legend(fontsize=7)

    # Panel 6: Participation ratio of dominant eigenvector
    ax6 = fig.add_subplot(gs[1, 2])
    scatter_by_family(ax6, lambda r: r['h'], lambda r: r['participation_ratio'])
    ax6.set_xlabel('Coxeter number h')
    ax6.set_ylabel('Participation ratio (PR)')
    ax6.set_title('Dominant eigenvector delocalization')

    # Panel 7: Cartan vs Root weight fraction
    ax7 = fig.add_subplot(gs[2, 0])
    scatter_by_family(ax7, lambda r: r['h'], lambda r: r['cartan_weight'] * 100)
    ax7.set_xlabel('Coxeter number h')
    ax7.set_ylabel('Cartan weight (%)')
    ax7.set_title('Eigenvector: Cartan component')

    # Panel 8: d_max vs dim — structural relationship
    ax8 = fig.add_subplot(gs[2, 1])
    scatter_by_family(ax8, lambda r: r['dim'], lambda r: r['d_max'])
    # Mark exact h-1 matches
    exact = [r for r in results if abs(r['deviation']) < 0.05]
    if exact:
        ax8.scatter([r['dim'] for r in exact], [r['d_max'] for r in exact],
                    c='none', marker='o', s=200, edgecolors='red', linewidth=2,
                    label='h-1 exact', zorder=5)
    ax8.set_xlabel('dim')
    ax8.set_ylabel(r'$d_{max}$')
    ax8.set_title(r'$d_{max}$ vs dim (red = h-1 exact)')
    ax8.legend(fontsize=7)

    # Panel 9: log-log best fit
    ax9 = fig.add_subplot(gs[2, 2])
    data_filt = [r for r in results if r['dim'] > 3]
    d_mean_arr = np.array([r['d_mean'] for r in data_filt])
    y_arr = np.array([r['max_nontrivial'] for r in data_filt])
    for fam in ['A', 'B', 'C', 'D', 'G', 'F', 'E']:
        subset = [r for r in data_filt if r['family'] == fam]
        if not subset:
            continue
        ax9.scatter(np.log([r['d_mean'] for r in subset]),
                    np.log([r['max_nontrivial'] for r in subset]),
                    c=family_colors[fam], marker=family_markers[fam],
                    s=60, label=fam, edgecolors='k', linewidth=0.3)
    if 'a*d_mean^b' in formulas:
        popt = formulas['a*d_mean^b'][0]
        xfit = np.linspace(np.log(min(d_mean_arr)), np.log(max(d_mean_arr)), 100)
        ax9.plot(xfit, np.log(popt[0]) + popt[1] * xfit, 'k--',
                 label=f'{popt[0]:.3f}*d^{{{popt[1]:.3f}}}')
    ax9.set_xlabel(r'$\ln(\langle d \rangle)$')
    ax9.set_ylabel(r'$\ln($max$|\lambda_{nt}|)$')
    ax9.set_title('Power-law fit (log-log)')
    ax9.legend(fontsize=7, ncol=2)

    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved: {output_path}")


def run_analysis():
    """Main P10 pipeline."""
    print("=" * 70)
    print("P10: Spectral Extrema Characterization")
    print("=" * 70)

    # Phase 1: Compute all algebras
    print("\n=== Phase 1: Computing all algebras ===")
    results = []
    t0 = time.time()
    for name, info in ALGEBRAS.items():
        try:
            r = compute_algebra(name, info)
            results.append(r)
        except Exception as e:
            print(f"  FAILED {name}: {e}")

    elapsed = time.time() - t0
    print(f"\nComputed {len(results)} algebras in {elapsed:.1f}s")

    # Print master table
    print(f"\n{'Name':>8} {'h':>4} {'dim':>5} {'d_mean':>7} {'d_max':>6} "
          f"{'max_nt':>8} {'h-1':>5} {'dev':>8} {'ratio':>7} {'SL':>3} {'Ram':>4}")
    print("-" * 80)
    for r in sorted(results, key=lambda x: x['h']):
        sl_str = 'Y' if r['simply_laced'] else 'N'
        ram_str = 'Y' if r['is_ramanujan'] else 'N'
        exact = ' ***' if abs(r['deviation']) < 0.05 else ''
        print(f"{r['name']:>8} {r['h']:4d} {r['dim']:5d} {r['d_mean']:7.2f} {r['d_max']:6d} "
              f"{r['max_nontrivial']:8.3f} {r['h_minus_1']:5d} {r['deviation']:+8.3f} "
              f"{r['ratio']:7.4f} {sl_str:>3} {ram_str:>4}{exact}")

    # Phase 2: Regression
    formulas_dict = regression_analysis(results)

    # Phase 3: Simply-laced analysis
    laced_analysis(results)

    # Phase 4: Eigenvector analysis
    eigenvector_analysis(results)

    # Phase 5: h-1 locus
    h1_locus_analysis(results)

    # Phase 6: Universal formula
    formulas, data = universal_formula(results)

    # Phase 7: Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    n_exact = sum(1 for r in results if abs(r['deviation']) < 0.05)
    n_near = sum(1 for r in results if abs(r['deviation']) < 0.5)
    print(f"Total algebras: {len(results)}")
    print(f"Exact h-1 matches (|dev| < 0.05): {n_exact}")
    print(f"Near h-1 (|dev| < 0.5): {n_near}")
    exact_names = [r['name'] for r in results if abs(r['deviation']) < 0.05]
    print(f"Exact locus: {exact_names}")

    sl_exact = [r['name'] for r in results if abs(r['deviation']) < 0.05 and r['simply_laced']]
    nsl_exact = [r['name'] for r in results if abs(r['deviation']) < 0.05 and not r['simply_laced']]
    print(f"  Simply-laced: {sl_exact}")
    print(f"  Non-simply-laced: {nsl_exact}")

    # Which formula wins?
    if formulas:
        best = min(formulas.items(), key=lambda kv: kv[1][1])
        print(f"\nBest formula: {best[0]} — MSE={best[1][1]:.4f}")
        print(f"  {best[1][2]}")

    # Figure
    output_path = os.path.join(os.path.dirname(__file__), 'spectral_extrema.png')
    create_figure(results, formulas, output_path)

    return results


if __name__ == '__main__':
    run_analysis()
