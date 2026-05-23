"""
P11: Order Parameter Phi vs Landscape Accessibility
=====================================================
Tests the hypothesis from 7.6.5: the composite order parameter Phi
should correlate with landscape basin accessibility.

Deep-ordered algebras (Phi < 0.2) → large basins → 100% success
Critical algebras (Phi ~ 0.4) → narrow basins → blind discovery fails
Disordered algebras (Phi > 0.5) → landscape-inaccessible

Data sources:
  - Phi: computed from Ihara zeta poles (P5 formula)
  - Spectral: adjacency graph eigenvalues (P1-P10)
  - Accessibility proxies: mixing time, spectral gap, expansion
  - Experimental ground truth: neural experiment success/failure

Analyses:
  1. Recompute Phi for all 37 algebras (same formula as P5)
  2. Compute accessibility proxies: mixing time, Cheeger bound, gap
  3. Correlate Phi with proxy measures
  4. Test against experimental ground truth
  5. Phase diagram overlay: Phi vs accessibility vs Ramanujan
  6. Predictive model: which algebras are landscape-accessible?
  7. 9-panel visualization
"""

import sys, os, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr, spearmanr, kendalltau
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(__file__))
from ramanujan_extended import (
    roots_A, roots_B, roots_C, roots_D,
    roots_G2, roots_F4, roots_E6, roots_E7, roots_E8,
    adjacency_from_roots, spectral_analysis,
)
from ihara_zeta_adjoint import ihara_poles_bass, classify_poles


# ============================================================
# Algebra catalog (same 37 as P5)
# ============================================================

COXETER = {
    'su(2)': 2, 'su(3)': 3, 'su(4)': 4, 'su(5)': 5, 'su(6)': 6,
    'su(7)': 7, 'su(8)': 8, 'su(9)': 9, 'su(10)': 10, 'su(11)': 11, 'su(12)': 12,
    'so(5)': 4, 'so(7)': 6, 'so(9)': 8, 'so(11)': 10, 'so(13)': 12,
    'so(15)': 14, 'so(17)': 16,
    'sp(4)': 4, 'sp(6)': 6, 'sp(8)': 8, 'sp(10)': 10, 'sp(12)': 12,
    'sp(14)': 14, 'sp(16)': 16,
    'so(6)': 4, 'so(8)': 6, 'so(10)': 8, 'so(12)': 10, 'so(14)': 12,
    'so(16)': 14, 'so(18)': 16,
    'G2': 6, 'F4': 12, 'E6': 12, 'E7': 18, 'E8': 30,
}

SIMPLY_LACED = {
    'su(2)': True, 'su(3)': True, 'su(4)': True, 'su(5)': True,
    'su(6)': True, 'su(7)': True, 'su(8)': True, 'su(9)': True,
    'su(10)': True, 'su(11)': True, 'su(12)': True,
    'so(5)': False, 'so(7)': False, 'so(9)': False, 'so(11)': False,
    'so(13)': False, 'so(15)': False, 'so(17)': False,
    'sp(4)': False, 'sp(6)': False, 'sp(8)': False, 'sp(10)': False,
    'sp(12)': False, 'sp(14)': False, 'sp(16)': False,
    'so(6)': True, 'so(8)': True, 'so(10)': True, 'so(12)': True,
    'so(14)': True, 'so(16)': True, 'so(18)': True,
    'G2': False, 'F4': False,
    'E6': True, 'E7': True, 'E8': True,
}

FAMILIES = {}
for n in range(1, 12):
    FAMILIES[f"su({n+1})"] = 'A'
for n in range(2, 9):
    FAMILIES[f"so({2*n+1})"] = 'B'
for n in range(2, 9):
    FAMILIES[f"sp({2*n})"] = 'C'
for n in range(3, 10):
    FAMILIES[f"so({2*n})"] = 'D'
FAMILIES['G2'] = 'G'
FAMILIES['F4'] = 'F'
FAMILIES['E6'] = 'E'
FAMILIES['E7'] = 'E'
FAMILIES['E8'] = 'E'

# Experimental ground truth from neural experiments
# Scale: 0 = inaccessible, 0.5 = hard, 1.0 = easy (100% success)
EXPERIMENTAL_DATA = {
    'su(2)': {'success': 1.0, 'label': '100% closure (SS-PEG)'},
    'su(3)': {'success': 1.0, 'label': '100% closure (SS-PEG)'},
    'su(4)': {'success': 1.0, 'label': '100% closure (SS-PEG)'},
    'su(5)': {'success': 1.0, 'label': '100% closure (SS-PEG)'},
    'G2':    {'success': 0.8, 'label': 'Discovered (11c, gamma>=5)'},
    'F4':    {'success': 0.3, 'label': '0/128 random; warm-start only (11d,11f)'},
    'E6':    {'success': 0.1, 'label': 'Failed blind (11e); warm-start partial'},
    'E8':    {'success': 0.0, 'label': 'Never attempted (presumed inaccessible)'},
}


def build_catalog():
    """Build 37-algebra catalog."""
    catalog = []
    for n in range(1, 12):
        r, roots, simple = roots_A(n)
        catalog.append((f"su({n+1})", r, roots, simple))
    for n in range(2, 9):
        r, roots, simple = roots_B(n)
        catalog.append((f"so({2*n+1})", r, roots, simple))
    for n in range(2, 9):
        r, roots, simple = roots_C(n)
        catalog.append((f"sp({2*n})", r, roots, simple))
    for n in range(3, 10):
        r, roots, simple = roots_D(n)
        catalog.append((f"so({2*n})", r, roots, simple))
    for exc_name, fn in [('G2', roots_G2), ('F4', roots_F4),
                          ('E6', roots_E6), ('E7', roots_E7), ('E8', roots_E8)]:
        r, roots, simple = fn()
        catalog.append((exc_name, r, roots, simple))
    return catalog


# ============================================================
# Core computation
# ============================================================

def compute_all():
    """Compute spectral + pole + Phi data for all 37 algebras."""
    catalog = build_catalog()
    results = []

    print("=" * 100)
    print("  P11: Order Parameter Phi vs Landscape Accessibility")
    print("=" * 100)
    print()

    header = (f"  {'Algebra':>10s}  {'dim':>4s}  {'h':>3s}  {'Ram':>4s}  "
              f"{'ratio':>7s}  {'inner/rc':>9s}  {'f_ring':>7s}  "
              f"{'sigma_r':>8s}  {'Phi':>7s}  {'Phase':>10s}  "
              f"{'gap':>7s}  {'t_mix':>7s}")
    print(header)
    print("  " + "-" * 110)

    for name, rank, roots, simple in catalog:
        adj = adjacency_from_roots(rank, roots, simple)
        dim = adj.shape[0]
        spec = spectral_analysis(adj)
        h = COXETER.get(name, -1)

        # --- Ihara zeta poles ---
        poles = ihara_poles_bass(adj)
        info = classify_poles(poles, adj)
        nt = info['nontrivial']
        radii = np.abs(nt) if len(nt) > 0 else np.array([])
        rc = info['ram_radius']

        eps = 0.05
        in_ring = (np.sum((radii >= rc*(1-eps)) & (radii <= rc*(1+eps)))
                   if len(radii) > 0 else 0)
        frac_ring = in_ring / len(radii) if len(radii) > 0 else 1.0
        inner_ratio = np.min(radii) / rc if len(radii) > 0 and rc > 0 else 1.0
        sigma_r = np.std(radii) if len(radii) > 0 else 0.0

        # --- Composite order parameter Phi (P5 formula) ---
        phi_inner = 1.0 - min(inner_ratio, 1.0)
        phi_frac = 1.0 - frac_ring
        phi_sigma = min(sigma_r / 0.05, 1.0)
        phi_ratio = min(max(spec['ratio'] - 1.0, 0.0) / 1.0, 1.0)
        Phi = 0.4*phi_inner + 0.3*phi_frac + 0.2*phi_sigma + 0.1*phi_ratio

        phase = "ordered" if Phi < 0.3 else ("critical" if Phi < 0.5 else "disordered")

        # --- Spectral accessibility proxies ---
        eigenvalues = np.sort(np.linalg.eigvalsh(adj.astype(float)))[::-1]
        lam1 = eigenvalues[0]   # Perron-Frobenius
        lam2 = eigenvalues[1]   # second largest
        lam_min = eigenvalues[-1]
        max_nt = spec['max_nontrivial']

        # Spectral gap (absolute)
        spectral_gap = lam1 - max_nt

        # Mixing time estimate: t_mix ~ log(dim) / log(lam1 / max_nt)
        if max_nt > 0 and lam1 > max_nt:
            t_mix = np.log(dim) / np.log(lam1 / max_nt)
        else:
            t_mix = float('inf')

        # Cheeger constant lower bound: h_C >= (lam1 - lam2) / 2
        # (Cheeger inequality for regular graphs)
        degs = np.sum(adj, axis=1).astype(float)
        d_mean = np.mean(degs)
        d_max = np.max(degs)
        is_regular = (np.max(degs) == np.min(degs))

        if is_regular:
            cheeger_lb = (lam1 - max_nt) / (2 * lam1)
        else:
            cheeger_lb = (lam1 - max_nt) / (2 * d_max)

        # Normalized spectral gap
        norm_gap = spectral_gap / lam1 if lam1 > 0 else 0

        # Ramanujan surplus: how far inside (positive) or outside (negative)
        # the Ramanujan bound is max_nt
        ram_bound = 2 * np.sqrt(lam1 - 1) if lam1 > 1 else 0
        ram_surplus = ram_bound - max_nt  # positive = Ramanujan, negative = not

        # Relative expansion rate
        expansion = spectral_gap / np.log(dim) if dim > 1 else 0

        r = {
            'name': name,
            'family': FAMILIES.get(name, '?'),
            'simply_laced': SIMPLY_LACED.get(name, True),
            'dim': dim,
            'rank': rank,
            'h': h,
            'd_mean': d_mean,
            'd_max': float(d_max),
            'is_regular': is_regular,
            'ratio': spec['ratio'],
            'is_ramanujan': bool(spec['is_ramanujan']),
            'max_nt': max_nt,
            'lam1': lam1,
            # Pole data
            'inner_ratio': inner_ratio,
            'frac_ring': frac_ring,
            'sigma_r': sigma_r,
            'rc': rc,
            'n_poles_nt': len(nt),
            # Phi
            'Phi': Phi,
            'phase': phase,
            # Accessibility proxies
            'spectral_gap': spectral_gap,
            't_mix': t_mix,
            'cheeger_lb': cheeger_lb,
            'norm_gap': norm_gap,
            'ram_surplus': ram_surplus,
            'expansion': expansion,
            # Experimental
            'has_exp': name in EXPERIMENTAL_DATA,
            'exp_success': EXPERIMENTAL_DATA.get(name, {}).get('success', None),
            'exp_label': EXPERIMENTAL_DATA.get(name, {}).get('label', ''),
        }
        results.append(r)

        is_ram_str = "YES" if r['is_ramanujan'] else "no"
        print(f"  {name:>10s}  {dim:4d}  {h:3d}  {is_ram_str:>4s}  "
              f"{spec['ratio']:7.4f}  {inner_ratio:9.4f}  {frac_ring:7.3f}  "
              f"{sigma_r:8.5f}  {Phi:7.4f}  {phase:>10s}  "
              f"{spectral_gap:7.3f}  {t_mix:7.2f}")

    return results


# ============================================================
# Analysis 1: Phi correlations
# ============================================================

def correlation_analysis(results):
    """Correlate Phi with all accessibility proxies."""
    print()
    print("=" * 90)
    print("  ANALYSIS 1: Phi Correlations with Accessibility Proxies")
    print("=" * 90)
    print()

    Phi = np.array([r['Phi'] for r in results])
    proxies = {
        'spectral_gap':  np.array([r['spectral_gap'] for r in results]),
        't_mix':         np.array([r['t_mix'] for r in results]),
        'norm_gap':      np.array([r['norm_gap'] for r in results]),
        'cheeger_lb':    np.array([r['cheeger_lb'] for r in results]),
        'ram_surplus':   np.array([r['ram_surplus'] for r in results]),
        'expansion':     np.array([r['expansion'] for r in results]),
        'ratio':         np.array([r['ratio'] for r in results]),
        'h':             np.array([r['h'] for r in results], dtype=float),
        'dim':           np.array([r['dim'] for r in results], dtype=float),
        'log_dim':       np.log(np.array([r['dim'] for r in results], dtype=float)),
    }

    print(f"  {'Proxy':>15s}  {'Pearson r':>10s}  {'p-value':>10s}  "
          f"{'Spearman rho':>13s}  {'p-value':>10s}  {'Direction':>10s}")
    print("  " + "-" * 80)

    corr_data = {}
    for pname, pvals in proxies.items():
        mask = np.isfinite(pvals) & np.isfinite(Phi)
        if mask.sum() < 5:
            continue
        pr, pp = pearsonr(Phi[mask], pvals[mask])
        sr, sp = spearmanr(Phi[mask], pvals[mask])
        direction = "Phi up => " + ("up" if pr > 0 else "down")
        print(f"  {pname:>15s}  {pr:10.4f}  {pp:10.2e}  {sr:13.4f}  {sp:10.2e}  {direction:>10s}")
        corr_data[pname] = {'pearson': pr, 'p_pearson': pp, 'spearman': sr, 'p_spearman': sp}

    return corr_data


# ============================================================
# Analysis 2: Phase boundary characterization
# ============================================================

def phase_analysis(results):
    """Characterize the three phases."""
    print()
    print("=" * 90)
    print("  ANALYSIS 2: Phase Boundary Characterization")
    print("=" * 90)
    print()

    phases = {'ordered': [], 'critical': [], 'disordered': []}
    for r in results:
        phases[r['phase']].append(r)

    for pname in ['ordered', 'critical', 'disordered']:
        group = phases[pname]
        if not group:
            continue
        gaps = [r['spectral_gap'] for r in group]
        tmix = [r['t_mix'] for r in group]
        phis = [r['Phi'] for r in group]
        rams = sum(1 for r in group if r['is_ramanujan'])
        print(f"  {pname.upper():>12s} (n={len(group):2d}):  "
              f"Phi = {np.mean(phis):.3f} +/- {np.std(phis):.3f},  "
              f"gap = {np.mean(gaps):.3f} +/- {np.std(gaps):.3f},  "
              f"t_mix = {np.mean(tmix):.2f} +/- {np.std(tmix):.2f},  "
              f"Ramanujan = {rams}/{len(group)}")

    return phases


# ============================================================
# Analysis 3: Experimental ground truth test
# ============================================================

def experimental_test(results):
    """Test Phi predictions against experimental data."""
    print()
    print("=" * 90)
    print("  ANALYSIS 3: Experimental Ground Truth Validation")
    print("=" * 90)
    print()

    exp_results = [r for r in results if r['has_exp']]
    exp_results.sort(key=lambda r: r['Phi'])

    print(f"  {'Algebra':>10s}  {'Phi':>7s}  {'Phase':>10s}  {'Ram':>4s}  "
          f"{'gap':>7s}  {'t_mix':>7s}  {'Exp success':>12s}  {'Experiment'}")
    print("  " + "-" * 105)
    for r in exp_results:
        ram_str = "YES" if r['is_ramanujan'] else "no"
        print(f"  {r['name']:>10s}  {r['Phi']:7.4f}  {r['phase']:>10s}  {ram_str:>4s}  "
              f"{r['spectral_gap']:7.3f}  {r['t_mix']:7.2f}  "
              f"{r['exp_success']:12.1f}  {r['exp_label']}")

    # Correlation: Phi vs experimental success (should be negative)
    phis = np.array([r['Phi'] for r in exp_results])
    successes = np.array([r['exp_success'] for r in exp_results])
    if len(phis) >= 4:
        pr, pp = pearsonr(phis, successes)
        sr, sp = spearmanr(phis, successes)
        kr, kp = kendalltau(phis, successes)
        print()
        print(f"  Phi vs experimental success:")
        print(f"    Pearson  r = {pr:.4f}, p = {pp:.4g}")
        print(f"    Spearman r = {sr:.4f}, p = {sp:.4g}")
        print(f"    Kendall  t = {kr:.4f}, p = {kp:.4g}")

    # Monotonicity check: does Phi strictly order success?
    print()
    print("  Monotonicity check (Phi increases => success decreases):")
    monotonic = True
    for i in range(len(exp_results) - 1):
        r1, r2 = exp_results[i], exp_results[i+1]
        ok = r1['exp_success'] >= r2['exp_success']
        marker = "OK" if ok else "VIOLATION"
        if not ok:
            monotonic = False
        print(f"    {r1['name']:>8s} (Phi={r1['Phi']:.3f}, s={r1['exp_success']:.1f}) -> "
              f"{r2['name']:>8s} (Phi={r2['Phi']:.3f}, s={r2['exp_success']:.1f})  [{marker}]")

    print(f"\n  Monotonic ordering: {'YES' if monotonic else 'NO (has violations)'}")

    return exp_results


# ============================================================
# Analysis 4: Accessibility predictor ranking
# ============================================================

def predictor_ranking(results):
    """Rank all candidate predictors of landscape accessibility."""
    print()
    print("=" * 90)
    print("  ANALYSIS 4: Accessibility Predictor Ranking")
    print("=" * 90)
    print()

    exp_results = [r for r in results if r['has_exp']]
    if len(exp_results) < 4:
        print("  Not enough experimental data points.")
        return {}

    successes = np.array([r['exp_success'] for r in exp_results])

    candidates = {
        'Phi':          np.array([r['Phi'] for r in exp_results]),
        'ratio':        np.array([r['ratio'] for r in exp_results]),
        'spectral_gap': np.array([r['spectral_gap'] for r in exp_results]),
        'norm_gap':     np.array([r['norm_gap'] for r in exp_results]),
        't_mix':        np.array([r['t_mix'] for r in exp_results]),
        'cheeger_lb':   np.array([r['cheeger_lb'] for r in exp_results]),
        'ram_surplus':  np.array([r['ram_surplus'] for r in exp_results]),
        'h':            np.array([r['h'] for r in exp_results], dtype=float),
        'dim':          np.array([r['dim'] for r in exp_results], dtype=float),
        'is_ramanujan': np.array([1.0 if r['is_ramanujan'] else 0.0 for r in exp_results]),
        'inner_ratio':  np.array([r['inner_ratio'] for r in exp_results]),
        'frac_ring':    np.array([r['frac_ring'] for r in exp_results]),
    }

    print(f"  {'Predictor':>15s}  {'Spearman rho':>13s}  {'p-value':>10s}  "
          f"{'|rho|':>6s}  {'Correct order':>14s}")
    print("  " + "-" * 70)

    ranking = []
    for pname, pvals in candidates.items():
        mask = np.isfinite(pvals) & np.isfinite(successes)
        if mask.sum() < 4:
            continue
        sr, sp = spearmanr(pvals[mask], successes[mask])
        # Count correct pairwise orderings
        n_pairs = 0
        n_correct = 0
        for i in range(mask.sum()):
            for j in range(i+1, mask.sum()):
                if successes[mask][i] != successes[mask][j]:
                    n_pairs += 1
                    pred_order = pvals[mask][i] < pvals[mask][j]  # lower pred = more accessible?
                    actual_order = successes[mask][i] > successes[mask][j]
                    # For Phi: higher Phi => less accessible => predict success inversely
                    # For gap: higher gap => more accessible => predict directly
                    if sr < 0:  # inverse predictor
                        if (pvals[mask][i] < pvals[mask][j]) == (successes[mask][i] > successes[mask][j]):
                            n_correct += 1
                    else:  # direct predictor
                        if (pvals[mask][i] > pvals[mask][j]) == (successes[mask][i] > successes[mask][j]):
                            n_correct += 1
        frac_correct = n_correct / n_pairs if n_pairs > 0 else 0
        print(f"  {pname:>15s}  {sr:13.4f}  {sp:10.4g}  {abs(sr):6.3f}  "
              f"{n_correct}/{n_pairs} ({frac_correct:.0%})")
        ranking.append((pname, abs(sr), sr, sp, frac_correct))

    ranking.sort(key=lambda x: x[1], reverse=True)
    print()
    print("  Top predictors (by |Spearman rho|):")
    for i, (pname, absr, sr, sp, frac) in enumerate(ranking[:5]):
        print(f"    {i+1}. {pname}: |rho| = {absr:.4f} (p = {sp:.4g})")

    return dict(ranking=ranking, candidates=candidates)


# ============================================================
# Analysis 5: Extrapolated predictions
# ============================================================

def extrapolated_predictions(results):
    """Predict accessibility for ALL algebras based on Phi."""
    print()
    print("=" * 90)
    print("  ANALYSIS 5: Extrapolated Accessibility Predictions")
    print("=" * 90)
    print()

    # Sort by Phi
    sorted_r = sorted(results, key=lambda r: r['Phi'])

    print(f"  {'Algebra':>10s}  {'Phi':>7s}  {'Phase':>10s}  {'Ram':>4s}  "
          f"{'gap':>7s}  {'t_mix':>7s}  {'Predicted':>10s}  {'Exp':>6s}")
    print("  " + "-" * 85)

    for r in sorted_r:
        ram_str = "YES" if r['is_ramanujan'] else "no"
        # Predicted accessibility based on Phi
        if r['Phi'] < 0.15:
            pred = "EASY"
        elif r['Phi'] < 0.30:
            pred = "accessible"
        elif r['Phi'] < 0.45:
            pred = "hard"
        elif r['Phi'] < 0.55:
            pred = "very hard"
        else:
            pred = "BLOCKED"

        exp_str = f"{r['exp_success']:.1f}" if r['has_exp'] else "---"
        print(f"  {r['name']:>10s}  {r['Phi']:7.4f}  {r['phase']:>10s}  {ram_str:>4s}  "
              f"{r['spectral_gap']:7.3f}  {r['t_mix']:7.2f}  {pred:>10s}  {exp_str:>6s}")

    # Count predictions per category
    print()
    cats = {'EASY': 0, 'accessible': 0, 'hard': 0, 'very hard': 0, 'BLOCKED': 0}
    for r in sorted_r:
        if r['Phi'] < 0.15:
            cats['EASY'] += 1
        elif r['Phi'] < 0.30:
            cats['accessible'] += 1
        elif r['Phi'] < 0.45:
            cats['hard'] += 1
        elif r['Phi'] < 0.55:
            cats['very hard'] += 1
        else:
            cats['BLOCKED'] += 1

    print("  Prediction summary:")
    for cat, count in cats.items():
        print(f"    {cat:>10s}: {count:2d}/37 ({count/37*100:.0f}%)")


# ============================================================
# Analysis 6: Phi decomposition
# ============================================================

def phi_decomposition(results):
    """Decompose Phi into its 4 components and analyze which drives the transition."""
    print()
    print("=" * 90)
    print("  ANALYSIS 6: Phi Component Decomposition")
    print("=" * 90)
    print()

    print(f"  {'Algebra':>10s}  {'Phi':>7s}  {'phi_inner':>10s}  {'phi_frac':>10s}  "
          f"{'phi_sigma':>10s}  {'phi_ratio':>10s}  {'Dominant':>10s}")
    print("  " + "-" * 85)

    components = {'inner': [], 'frac': [], 'sigma': [], 'ratio': []}
    for r in sorted(results, key=lambda x: x['Phi']):
        phi_inner = 1.0 - min(r['inner_ratio'], 1.0)
        phi_frac = 1.0 - r['frac_ring']
        phi_sigma = min(r['sigma_r'] / 0.05, 1.0)
        phi_ratio = min(max(r['ratio'] - 1.0, 0.0) / 1.0, 1.0)

        weighted = {
            'inner': 0.4 * phi_inner,
            'frac': 0.3 * phi_frac,
            'sigma': 0.2 * phi_sigma,
            'ratio': 0.1 * phi_ratio,
        }
        dominant = max(weighted, key=weighted.get)

        for k in components:
            components[k].append(weighted[k])

        print(f"  {r['name']:>10s}  {r['Phi']:7.4f}  "
              f"{0.4*phi_inner:10.4f}  {0.3*phi_frac:10.4f}  "
              f"{0.2*phi_sigma:10.4f}  {0.1*phi_ratio:10.4f}  {dominant:>10s}")

    # Which component contributes most variance?
    print()
    print("  Variance contribution by component:")
    for k, vals in components.items():
        v = np.var(vals)
        print(f"    {k:>10s}: var = {v:.6f}")

    return components


# ============================================================
# Visualization
# ============================================================

def create_figure(results, corr_data, phases, output_path):
    """9-panel visualization for P11."""
    fig = plt.figure(figsize=(20, 16))
    fig.suptitle(r"P11: Order Parameter $\Phi$ vs Landscape Accessibility",
                 fontsize=16, fontweight='bold')
    gs = GridSpec(3, 3, figure=fig, hspace=0.38, wspace=0.32)

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
                       s=60, label=fam, zorder=3, edgecolors='k', linewidth=0.3,
                       **kwargs)

    # --- Panel 1: Phi vs Coxeter number h ---
    ax1 = fig.add_subplot(gs[0, 0])
    scatter_by_family(ax1, lambda r: r['h'], lambda r: r['Phi'])
    ax1.axhline(0.3, color='orange', ls='--', alpha=0.5, label=r'$\Phi$ = 0.3')
    ax1.axhline(0.5, color='red', ls='--', alpha=0.5, label=r'$\Phi$ = 0.5')
    ax1.set_xlabel('Coxeter number $h$')
    ax1.set_ylabel(r'$\Phi$')
    ax1.set_title(r'$\Phi$ vs Coxeter number')
    ax1.legend(fontsize=7, ncol=2)

    # --- Panel 2: Phi vs spectral gap ---
    ax2 = fig.add_subplot(gs[0, 1])
    scatter_by_family(ax2, lambda r: r['Phi'], lambda r: r['spectral_gap'])
    ax2.axvline(0.3, color='orange', ls='--', alpha=0.4)
    ax2.axvline(0.5, color='red', ls='--', alpha=0.4)
    ax2.set_xlabel(r'$\Phi$')
    ax2.set_ylabel('Spectral gap')
    ax2.set_title(r'Spectral gap vs $\Phi$')
    ax2.legend(fontsize=7, ncol=2)

    # --- Panel 3: Phi vs mixing time ---
    ax3 = fig.add_subplot(gs[0, 2])
    scatter_by_family(ax3, lambda r: r['Phi'], lambda r: r['t_mix'])
    ax3.axvline(0.3, color='orange', ls='--', alpha=0.4)
    ax3.axvline(0.5, color='red', ls='--', alpha=0.4)
    ax3.set_xlabel(r'$\Phi$')
    ax3.set_ylabel(r'$t_{mix} \sim \ln(n) / \ln(\lambda_1 / \lambda_{nt})$')
    ax3.set_title(r'Mixing time vs $\Phi$')
    ax3.legend(fontsize=7, ncol=2)

    # --- Panel 4: Phi vs Ramanujan surplus ---
    ax4 = fig.add_subplot(gs[1, 0])
    scatter_by_family(ax4, lambda r: r['Phi'], lambda r: r['ram_surplus'])
    ax4.axhline(0, color='k', ls='--', alpha=0.4, label='Ramanujan boundary')
    ax4.axvline(0.3, color='orange', ls='--', alpha=0.4)
    ax4.axvline(0.5, color='red', ls='--', alpha=0.4)
    ax4.set_xlabel(r'$\Phi$')
    ax4.set_ylabel(r'Ram. surplus: $2\sqrt{q} - \max|\lambda_{nt}|$')
    ax4.set_title(r'Ramanujan surplus vs $\Phi$')
    ax4.legend(fontsize=7)

    # --- Panel 5: Phase diagram (Phi vs norm gap, colored by phase) ---
    ax5 = fig.add_subplot(gs[1, 1])
    phase_colors = {'ordered': '#2196F3', 'critical': '#FF9800', 'disordered': '#F44336'}
    for pname, group in phases.items():
        if not group:
            continue
        xx = [r['Phi'] for r in group]
        yy = [r['norm_gap'] for r in group]
        ax5.scatter(xx, yy, c=phase_colors[pname], s=70, label=pname,
                    edgecolors='k', linewidth=0.4, zorder=3)
    # Highlight experimental points
    for r in results:
        if r['has_exp']:
            ax5.annotate(r['name'], (r['Phi'], r['norm_gap']),
                        fontsize=6, ha='center', va='bottom',
                        fontweight='bold', color='black')
    ax5.set_xlabel(r'$\Phi$')
    ax5.set_ylabel(r'Normalized spectral gap $\Delta / \lambda_1$')
    ax5.set_title('Phase diagram')
    ax5.legend(fontsize=7)

    # --- Panel 6: Experimental validation ---
    ax6 = fig.add_subplot(gs[1, 2])
    exp_r = [r for r in results if r['has_exp']]
    exp_r.sort(key=lambda r: r['Phi'])
    phis = [r['Phi'] for r in exp_r]
    succs = [r['exp_success'] for r in exp_r]
    names_exp = [r['name'] for r in exp_r]
    colors_exp = [phase_colors.get(r['phase'], 'gray') for r in exp_r]

    ax6.bar(range(len(phis)), succs, color=colors_exp, edgecolor='k', linewidth=0.5)
    ax6.set_xticks(range(len(phis)))
    ax6.set_xticklabels(names_exp, rotation=45, ha='right', fontsize=7)
    ax6.set_ylabel('Experimental success')
    ax6.set_title(r'Experimental validation (ordered by $\Phi$)')
    # Add Phi values as secondary labels
    for i, (phi, s) in enumerate(zip(phis, succs)):
        ax6.text(i, s + 0.03, f'{phi:.2f}', ha='center', va='bottom', fontsize=6, color='gray')

    # --- Panel 7: Phi component decomposition (stacked bar) ---
    ax7 = fig.add_subplot(gs[2, 0])
    sorted_results = sorted(results, key=lambda r: r['Phi'])
    names_sorted = [r['name'] for r in sorted_results]
    inner_vals = []
    frac_vals = []
    sigma_vals = []
    ratio_vals = []
    for r in sorted_results:
        phi_inner = 0.4 * (1.0 - min(r['inner_ratio'], 1.0))
        phi_frac = 0.3 * (1.0 - r['frac_ring'])
        phi_sigma = 0.2 * min(r['sigma_r'] / 0.05, 1.0)
        phi_ratio = 0.1 * min(max(r['ratio'] - 1.0, 0.0) / 1.0, 1.0)
        inner_vals.append(phi_inner)
        frac_vals.append(phi_frac)
        sigma_vals.append(phi_sigma)
        ratio_vals.append(phi_ratio)

    x_pos = np.arange(len(sorted_results))
    b1 = ax7.bar(x_pos, inner_vals, color='#1565C0', label='0.4(1-inner/rc)')
    b2 = ax7.bar(x_pos, frac_vals, bottom=inner_vals, color='#66BB6A', label='0.3(1-f_ring)')
    b3_bottom = np.array(inner_vals) + np.array(frac_vals)
    b3 = ax7.bar(x_pos, sigma_vals, bottom=b3_bottom, color='#FF9800', label=r'0.2$\sigma_r$/0.05')
    b4_bottom = b3_bottom + np.array(sigma_vals)
    b4 = ax7.bar(x_pos, ratio_vals, bottom=b4_bottom, color='#F44336', label='0.1(R-1)+')
    ax7.set_xticks(x_pos[::3])
    ax7.set_xticklabels([names_sorted[i] for i in range(0, len(names_sorted), 3)],
                         rotation=45, ha='right', fontsize=5)
    ax7.set_ylabel(r'$\Phi$ components')
    ax7.set_title(r'$\Phi$ decomposition (stacked)')
    ax7.legend(fontsize=6, loc='upper left')

    # --- Panel 8: Expansion rate vs Phi ---
    ax8 = fig.add_subplot(gs[2, 1])
    scatter_by_family(ax8, lambda r: r['Phi'], lambda r: r['expansion'])
    ax8.axvline(0.3, color='orange', ls='--', alpha=0.4)
    ax8.axvline(0.5, color='red', ls='--', alpha=0.4)
    ax8.set_xlabel(r'$\Phi$')
    ax8.set_ylabel(r'Expansion rate $\Delta / \ln(n)$')
    ax8.set_title(r'Expansion rate vs $\Phi$')
    ax8.legend(fontsize=7, ncol=2)

    # --- Panel 9: Unified predictor comparison ---
    ax9 = fig.add_subplot(gs[2, 2])
    # For algebras with experimental data, plot normalized Phi vs success
    exp_r = [r for r in results if r['has_exp']]
    if exp_r:
        phis_e = np.array([r['Phi'] for r in exp_r])
        succs_e = np.array([r['exp_success'] for r in exp_r])
        ax9.scatter(phis_e, succs_e, c='#1565C0', s=100, zorder=5,
                    edgecolors='k', linewidth=0.5)
        for r in exp_r:
            ax9.annotate(r['name'], (r['Phi'], r['exp_success']),
                        fontsize=7, ha='left', va='bottom', xytext=(3, 3),
                        textcoords='offset points')

        # Fit sigmoid: success = 1 / (1 + exp(a*(Phi - b)))
        try:
            def sigmoid(x, a, b):
                return 1.0 / (1.0 + np.exp(a * (x - b)))
            popt, _ = curve_fit(sigmoid, phis_e, succs_e, p0=[10, 0.3], maxfev=5000)
            phi_fit = np.linspace(0, 0.7, 100)
            ax9.plot(phi_fit, sigmoid(phi_fit, *popt), 'r--', alpha=0.7,
                     label=rf'$\sigma(\Phi)$: $a$={popt[0]:.1f}, $\Phi_c$={popt[1]:.2f}')
        except Exception:
            pass

        # Add ALL algebras as gray background
        all_phis = [r['Phi'] for r in results]
        ax9.scatter(all_phis, [-0.05] * len(all_phis), c='gray', s=20, alpha=0.4,
                    marker='|', zorder=1)
        ax9.set_xlabel(r'$\Phi$')
        ax9.set_ylabel('Experimental success')
        ax9.set_title(r'$\Phi$ as accessibility predictor')
        ax9.set_ylim(-0.12, 1.15)
        ax9.legend(fontsize=7)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\n  Figure saved: {output_path}")
    plt.close(fig)


# ============================================================
# Main
# ============================================================

def main():
    t0 = time.time()

    results = compute_all()
    corr_data = correlation_analysis(results)
    phases = phase_analysis(results)
    exp_results = experimental_test(results)
    pred_data = predictor_ranking(results)
    extrapolated_predictions(results)
    components = phi_decomposition(results)

    outpath = os.path.join(os.path.dirname(__file__), "phi_basin_correlation.png")
    create_figure(results, corr_data, phases, outpath)

    dt = time.time() - t0
    print(f"\n  Total runtime: {dt:.1f} s")
    print(f"  Algebras processed: {len(results)}")


if __name__ == "__main__":
    main()
