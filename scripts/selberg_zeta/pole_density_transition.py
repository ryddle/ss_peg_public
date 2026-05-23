"""
Priority 5: Density Transition in Ihara Zeta Poles
====================================================

Tracks how the Ihara zeta pole distribution changes as the algebra
"crystallizes" through the Ramanujan transition at h* ≈ 9.5.

ANALYSES:
  1. Parametric pole flow across the A_n family (individual pole tracking)
  2. Radial density of states ρ(r) evolution
  3. Angular distribution uniformity (Rayleigh test)
  4. Nearest-neighbor spacing statistics (Poisson vs GUE level repulsion)
  5. Critical scaling of order parameters near h*
  6. Z(u) landscape on complex plane
  7. Cross-family comparison (A, B, C, D families)
  8. Unified phase diagram: all order parameters vs h

Builds on P3 (ihara_zeta_adjoint.py) infrastructure.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import ks_1samp, pearsonr
from scipy.optimize import curve_fit
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from ramanujan_extended import (
    roots_A, roots_B, roots_C, roots_D,
    roots_G2, roots_F4, roots_E6, roots_E7, roots_E8,
    adjacency_from_roots, spectral_analysis,
)
from ihara_zeta_adjoint import ihara_poles_bass, classify_poles

COXETER_NUMBERS = {
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


# ============================================================
# Build catalog
# ============================================================

def build_catalog():
    """Build catalog of all algebras."""
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


def compute_pole_data(name, rank, all_roots, simple):
    """Compute pole data for one algebra."""
    adj = adjacency_from_roots(rank, all_roots, simple)
    dim = adj.shape[0]
    poles = ihara_poles_bass(adj)
    info = classify_poles(poles, adj)
    spec = spectral_analysis(adj)
    h = COXETER_NUMBERS.get(name, -1)

    nt = info['nontrivial']
    radii = np.abs(nt) if len(nt) > 0 else np.array([])
    angles = np.angle(nt) if len(nt) > 0 else np.array([])

    rc = info['ram_radius']
    eps = 0.05
    in_ring = np.sum((radii >= rc * (1 - eps)) & (radii <= rc * (1 + eps))) if len(radii) > 0 else 0
    frac_ring = in_ring / len(radii) if len(radii) > 0 else 1.0
    inner_ratio = np.min(radii) / rc if len(radii) > 0 and rc > 0 else 1.0
    sigma_r = np.std(radii) if len(radii) > 0 else 0.0

    return {
        'name': name, 'dim': dim, 'h': h,
        'poles': poles, 'nontrivial': nt, 'radii': radii, 'angles': angles,
        'rc': rc, 'frac_ring': frac_ring, 'inner_ratio': inner_ratio,
        'sigma_r': sigma_r, 'is_ramanujan': spec['is_ramanujan'],
        'ratio': spec['ratio'], 'adjacency': adj,
        'n_edges': int(np.sum(adj) / 2),
    }


# ============================================================
# Analysis 1: Parametric pole flow (A_n family)
# ============================================================

def pole_flow_analysis(results):
    """Track how pole radii evolve across the A_n family."""
    print("=" * 90)
    print("  ANALYSIS 1: Parametric Pole Flow — A_n Family")
    print("=" * 90)
    print()

    a_results = [r for r in results if r['name'].startswith('su(')]
    a_results.sort(key=lambda x: x['dim'])

    print(f"  {'Algebra':>10s}  {'dim':>4s}  {'h':>3s}  {'1/√q':>7s}  "
          f"{'r_min':>7s}  {'r_max':>7s}  {'r_med':>7s}  "
          f"{'inner/rc':>9s}  {'frac_ring':>10s}  {'σ_r':>7s}  {'Ram':>4s}")
    print("  " + "-" * 95)

    for r in a_results:
        if len(r['radii']) == 0:
            print(f"  {r['name']:>10s}  {r['dim']:4d}  {r['h']:3d}  "
                  f"{r['rc']:7.4f}  (no nontrivial poles)")
            continue
        print(f"  {r['name']:>10s}  {r['dim']:4d}  {r['h']:3d}  {r['rc']:7.4f}  "
              f"{np.min(r['radii']):7.4f}  {np.max(r['radii']):7.4f}  "
              f"{np.median(r['radii']):7.4f}  "
              f"{r['inner_ratio']:9.4f}  {r['frac_ring']:10.3f}  "
              f"{r['sigma_r']:7.4f}  {'YES' if r['is_ramanujan'] else 'NO':>4s}")

    # Compute pole flow velocity (change in inner/rc per unit h)
    print()
    print("  Pole migration velocity (Δ(inner/rc) / Δh):")
    for i in range(1, len(a_results)):
        if a_results[i]['h'] > 0 and a_results[i-1]['h'] > 0:
            dh = a_results[i]['h'] - a_results[i-1]['h']
            if dh > 0:
                d_inner = a_results[i]['inner_ratio'] - a_results[i-1]['inner_ratio']
                print(f"    {a_results[i-1]['name']:>8s} → {a_results[i]['name']:>8s}: "
                      f"Δ(inner/rc)/Δh = {d_inner/dh:+.4f}")

    return a_results


# ============================================================
# Analysis 2: Radial density of states
# ============================================================

def radial_density_analysis(results):
    """Compute radial density profiles ρ(r) and track evolution."""
    print()
    print("=" * 90)
    print("  ANALYSIS 2: Radial Density of States ρ(r)")
    print("=" * 90)
    print()

    a_results = [r for r in results if r['name'].startswith('su(')]
    a_results.sort(key=lambda x: x['dim'])

    density_data = {}
    for r in a_results:
        if len(r['radii']) < 5:
            continue
        radii = r['radii']
        rc = r['rc']

        # Normalized radii: r/rc
        norm_radii = radii / rc

        # Kernel density estimation (simple histogram-based)
        bins = np.linspace(0, 3.0, 61)
        counts, edges = np.histogram(norm_radii, bins=bins, density=True)
        centers = 0.5 * (edges[1:] + edges[:-1])

        # Find peak position and width
        peak_idx = np.argmax(counts)
        peak_pos = centers[peak_idx]
        peak_height = counts[peak_idx]

        # FWHM estimate
        half_max = peak_height / 2
        above = counts >= half_max
        if np.any(above):
            left = centers[above][0]
            right = centers[above][-1]
            fwhm = right - left
        else:
            fwhm = 0.0

        # Bimodality: count distinct peaks
        from scipy.signal import find_peaks
        peaks, props = find_peaks(counts, height=peak_height * 0.2,
                                   distance=5)
        n_peaks = len(peaks)

        density_data[r['name']] = {
            'centers': centers, 'density': counts,
            'peak_pos': peak_pos, 'fwhm': fwhm, 'n_peaks': n_peaks,
        }

        print(f"  {r['name']:>10s}: peak at r/rc={peak_pos:.3f}  "
              f"FWHM={fwhm:.3f}  #peaks={n_peaks}  "
              f"Ram={'YES' if r['is_ramanujan'] else 'NO'}")

    return density_data


# ============================================================
# Analysis 3: Angular distribution uniformity
# ============================================================

def angular_analysis(results):
    """Test angular uniformity of poles via Rayleigh test."""
    print()
    print("=" * 90)
    print("  ANALYSIS 3: Angular Distribution of Poles")
    print("=" * 90)
    print()
    print("  Rayleigh test: R̄ = |Σ exp(iθ)| / N  (0 = uniform, 1 = clustered)")
    print()
    print(f"  {'Algebra':>10s}  {'#nt':>5s}  {'R̄':>7s}  {'p-value':>9s}  "
          f"{'Uniform?':>9s}  {'Ram':>4s}")
    print("  " + "-" * 60)

    angular_data = {}
    for r in results:
        if len(r['angles']) < 5:
            continue
        angles = r['angles']
        n = len(angles)

        # Rayleigh statistic
        R_bar = np.abs(np.sum(np.exp(1j * angles))) / n
        # Rayleigh test p-value (approximate)
        z = n * R_bar**2
        p_val = np.exp(-z) * (1 + (2*z - z**2) / (4*n) - (24*z - 132*z**2 + 76*z**3 - 9*z**4) / (288*n**2))
        p_val = max(p_val, 1e-15)

        uniform = "YES" if p_val > 0.05 else "NO"

        # Angular entropy (vs uniform)
        n_bins = min(36, n // 3)
        if n_bins >= 4:
            hist, _ = np.histogram(angles, bins=n_bins, range=(-np.pi, np.pi))
            hist_norm = hist / hist.sum()
            entropy = -np.sum(hist_norm[hist_norm > 0] * np.log(hist_norm[hist_norm > 0]))
            max_entropy = np.log(n_bins)
            rel_entropy = entropy / max_entropy if max_entropy > 0 else 1.0
        else:
            rel_entropy = 1.0

        angular_data[r['name']] = {
            'R_bar': R_bar, 'p_val': p_val, 'uniform': uniform,
            'rel_entropy': rel_entropy,
        }

        print(f"  {r['name']:>10s}  {n:5d}  {R_bar:7.4f}  {p_val:9.2e}  "
              f"{uniform:>9s}  {'YES' if r['is_ramanujan'] else 'NO':>4s}")

    return angular_data


# ============================================================
# Analysis 4: Nearest-neighbor spacing statistics
# ============================================================

def spacing_analysis(results):
    """Compute nearest-neighbor spacing distributions for pole radii."""
    print()
    print("=" * 90)
    print("  ANALYSIS 4: Nearest-Neighbor Spacing Statistics")
    print("=" * 90)
    print()
    print("  s = normalized spacing of sorted pole radii")
    print("  Poisson: P(s) = exp(-s)  →  ⟨r⟩ = 1")
    print("  GUE:     P(s) ~ s²exp(-4s²/π)  →  level repulsion")
    print()
    print(f"  {'Algebra':>10s}  {'#nt':>5s}  {'⟨s⟩':>7s}  {'var(s)':>8s}  "
          f"{'⟨s²⟩/⟨s⟩²':>10s}  {'Type':>10s}  {'Ram':>4s}")
    print("  " + "-" * 70)

    spacing_data = {}
    for r in results:
        if len(r['radii']) < 10:
            continue
        radii = np.sort(r['radii'])

        # Normalized spacings
        spacings = np.diff(radii)
        mean_s = np.mean(spacings)
        if mean_s > 1e-12:
            s_norm = spacings / mean_s
        else:
            continue

        mean_sn = np.mean(s_norm)
        var_sn = np.var(s_norm)
        ratio_sq = np.mean(s_norm**2) / mean_sn**2 if mean_sn > 0 else 0

        # Discriminant: Poisson → var/⟨s⟩² = 1, s²/s² = 2
        #               GUE → var/⟨s⟩² ≈ 0.29, s²/s² ≈ 1.29
        if ratio_sq < 1.5:
            spacing_type = "GUE-like"
        elif ratio_sq < 1.8:
            spacing_type = "intermediate"
        else:
            spacing_type = "Poisson"

        spacing_data[r['name']] = {
            'spacings': s_norm, 'mean': mean_sn, 'var': var_sn,
            'ratio_sq': ratio_sq, 'type': spacing_type,
        }

        print(f"  {r['name']:>10s}  {len(r['radii']):5d}  {mean_sn:7.4f}  "
              f"{var_sn:8.4f}  {ratio_sq:10.4f}  {spacing_type:>10s}  "
              f"{'YES' if r['is_ramanujan'] else 'NO':>4s}")

    return spacing_data


# ============================================================
# Analysis 5: Critical scaling near h*
# ============================================================

def critical_scaling(results):
    """Fit power-law scaling of order parameters near h* ≈ 9.5."""
    print()
    print("=" * 90)
    print("  ANALYSIS 5: Critical Scaling Near h* ≈ 9.5")
    print("=" * 90)
    print()

    h_star = 9.5

    # Collect order parameters
    hs = []
    inner_ratios = []
    frac_rings = []
    sigmas = []
    ratios = []

    for r in results:
        if r['h'] <= 0 or len(r['radii']) < 5:
            continue
        hs.append(r['h'])
        inner_ratios.append(r['inner_ratio'])
        frac_rings.append(r['frac_ring'])
        sigmas.append(r['sigma_r'])
        ratios.append(r['ratio'])

    hs = np.array(hs)
    inner_ratios = np.array(inner_ratios)
    frac_rings = np.array(frac_rings)
    sigmas = np.array(sigmas)
    ratios = np.array(ratios)

    # 1. Inner ratio vs h: power-law fit for h > h*
    mask_above = hs > h_star
    mask_below = hs <= h_star

    print(f"  Order parameter means: h < h* vs h > h*")
    print(f"  {'Quantity':>20s}  {'h < h*':>12s}  {'h > h*':>12s}  {'Ratio':>8s}")
    print("  " + "-" * 60)

    for qty_name, qty in [('inner/rc', inner_ratios), ('frac_in_ring', frac_rings),
                           ('σ(r)', sigmas), ('ratio λ/2√q', ratios)]:
        if np.any(mask_below) and np.any(mask_above):
            m_below = np.mean(qty[mask_below])
            m_above = np.mean(qty[mask_above])
            ratio_val = m_above / m_below if m_below > 1e-10 else float('inf')
            print(f"  {qty_name:>20s}  {m_below:12.4f}  {m_above:12.4f}  {ratio_val:8.3f}")

    # 2. Power-law fit: ratio ∝ h^α
    print()
    print("  Power-law fits: quantity = a · h^b")
    print("  " + "-" * 60)

    def power_law(x, a, b):
        return a * np.power(x, b)

    for qty_name, qty in [('ratio λ/2√q', ratios), ('σ(r)', sigmas)]:
        valid = (hs > 1) & np.isfinite(qty) & (qty > 0)
        if np.sum(valid) >= 3:
            try:
                popt, pcov = curve_fit(power_law, hs[valid], qty[valid],
                                        p0=[0.3, 0.5], maxfev=5000)
                a, b = popt
                residuals = qty[valid] - power_law(hs[valid], a, b)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((qty[valid] - np.mean(qty[valid]))**2)
                r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 0
                print(f"  {qty_name:>20s}: {a:.4f} · h^{b:.4f}  (R²={r_sq:.4f})")
            except Exception:
                print(f"  {qty_name:>20s}: fit failed")

    # 3. Critical exponent: |inner/rc - 1| ∝ |h - h*|^ν
    print()
    print("  Critical exponent: |inner/rc - c| ~ |h - h*|^ν")

    # Use only data near the transition
    near_mask = np.abs(hs - h_star) < 6
    if np.sum(near_mask) >= 4:
        eps_h = np.abs(hs[near_mask] - h_star) + 0.1  # regularize
        eps_order = np.abs(inner_ratios[near_mask] - np.mean(inner_ratios))
        valid_c = eps_order > 1e-6
        if np.sum(valid_c) >= 3:
            try:
                popt, _ = curve_fit(power_law,
                                     eps_h[valid_c], eps_order[valid_c],
                                     p0=[0.1, 0.5], maxfev=5000)
                print(f"    ν ≈ {popt[1]:.3f}  (amplitude = {popt[0]:.4f})")
            except Exception:
                print(f"    fit failed")

    return {'hs': hs, 'inner_ratios': inner_ratios, 'frac_rings': frac_rings,
            'sigmas': sigmas, 'ratios': ratios, 'h_star': h_star}


# ============================================================
# Analysis 6: Z(u) landscape on complex plane
# ============================================================

def zeta_landscape(results):
    """Compute log|Z(u)^{-1}| on a grid in the complex plane for selected algebras."""
    print()
    print("=" * 90)
    print("  ANALYSIS 6: Z(u) Landscape on Complex Plane")
    print("=" * 90)
    print()

    landscape_data = {}
    selected = ['su(3)', 'su(5)', 'su(8)', 'su(12)']

    for name in selected:
        r = next((x for x in results if x['name'] == name), None)
        if r is None:
            continue
        adj = r['adjacency']
        n = adj.shape[0]
        D = np.diag(adj.sum(axis=1))
        Q = D - np.eye(n)

        # Grid in complex plane
        rc = r['rc']
        grid_range = min(1.5, 3 * rc)
        n_grid = 100
        re = np.linspace(-grid_range, grid_range, n_grid)
        im = np.linspace(-grid_range, grid_range, n_grid)
        RE, IM = np.meshgrid(re, im)
        U = RE + 1j * IM

        # Compute log|Z^{-1}(u)| = log|det(I - uA + u²Q)|
        log_Z_inv = np.zeros_like(RE)
        for i in range(n_grid):
            for j in range(n_grid):
                u = U[i, j]
                M = np.eye(n) - u * adj + u**2 * Q
                sign, logdet = np.linalg.slogdet(M)
                log_Z_inv[i, j] = logdet

        landscape_data[name] = {
            'RE': RE, 'IM': IM, 'log_Z_inv': log_Z_inv,
            'rc': rc, 'poles': r['nontrivial'],
        }
        print(f"  {name}: computed Z(u) landscape on {n_grid}×{n_grid} grid, "
              f"range=[{-grid_range:.2f}, {grid_range:.2f}]")

    return landscape_data


# ============================================================
# Analysis 7: Cross-family comparison
# ============================================================

def cross_family_analysis(results):
    """Compare pole density transition across A, B, C, D families."""
    print()
    print("=" * 90)
    print("  ANALYSIS 7: Cross-Family Comparison")
    print("=" * 90)
    print()

    families = {
        'A_n (su)': [r for r in results if r['name'].startswith('su(')],
        'B_n (so_odd)': [r for r in results if r['name'].startswith('so(')
                          and int(r['name'][3:-1]) % 2 == 1
                          and int(r['name'][3:-1]) <= 17],
        'C_n (sp)': [r for r in results if r['name'].startswith('sp(')],
        'D_n (so_even)': [r for r in results if r['name'].startswith('so(')
                           and int(r['name'][3:-1]) % 2 == 0
                           and int(r['name'][3:-1]) >= 6],
    }

    cross_data = {}
    for family_name, fam_results in families.items():
        fam_results.sort(key=lambda x: x['dim'])
        print(f"  {family_name}:")
        print(f"    {'Algebra':>10s}  {'dim':>4s}  {'h':>3s}  {'inner/rc':>9s}  "
              f"{'frac_ring':>10s}  {'σ_r':>7s}  {'Ram':>4s}")

        fam_data = []
        for r in fam_results:
            if len(r['radii']) < 5:
                continue
            print(f"    {r['name']:>10s}  {r['dim']:4d}  {r['h']:3d}  "
                  f"{r['inner_ratio']:9.4f}  {r['frac_ring']:10.3f}  "
                  f"{r['sigma_r']:7.4f}  {'YES' if r['is_ramanujan'] else 'NO':>4s}")
            fam_data.append(r)

        cross_data[family_name] = fam_data

        # Find transition rank
        for i in range(1, len(fam_data)):
            if fam_data[i-1]['is_ramanujan'] and not fam_data[i]['is_ramanujan']:
                print(f"    → Transition between {fam_data[i-1]['name']} and {fam_data[i]['name']} "
                      f"(h: {fam_data[i-1]['h']} → {fam_data[i]['h']})")
                break
        else:
            if all(r['is_ramanujan'] for r in fam_data):
                print(f"    → All Ramanujan in range tested")
            elif not any(r['is_ramanujan'] for r in fam_data):
                print(f"    → None Ramanujan in range tested")
        print()

    return cross_data


# ============================================================
# Analysis 8: Unified phase diagram
# ============================================================

def phase_diagram(results, crit_data):
    """Construct unified phase diagram: all order parameters vs h."""
    print()
    print("=" * 90)
    print("  ANALYSIS 8: Unified Phase Diagram")
    print("=" * 90)
    print()

    h_star = crit_data['h_star']

    # Compute "generalized order parameter" combining inner/rc, frac_ring, and spectral ratio
    print(f"  {'Algebra':>10s}  {'h':>3s}  {'inner/rc':>9s}  {'frac_ring':>10s}  "
          f"{'σ_r':>7s}  {'ratio':>7s}  {'Φ':>7s}  {'Phase':>10s}")
    print("  " + "-" * 80)

    for r in results:
        if r['h'] <= 0 or len(r['radii']) < 5:
            continue

        # Composite order parameter Φ ∈ [0,1]:
        # Φ = 0 → deep Ramanujan (ordered), Φ = 1 → deep non-Ramanujan (disordered)
        phi_inner = 1.0 - min(r['inner_ratio'], 1.0)
        phi_frac = 1.0 - r['frac_ring']
        phi_sigma = min(r['sigma_r'] / 0.05, 1.0)
        phi_ratio = min(max(r['ratio'] - 1.0, 0.0) / 1.0, 1.0)

        Phi = 0.4 * phi_inner + 0.3 * phi_frac + 0.2 * phi_sigma + 0.1 * phi_ratio

        phase = "ordered" if Phi < 0.3 else ("critical" if Phi < 0.5 else "disordered")

        print(f"  {r['name']:>10s}  {r['h']:3d}  {r['inner_ratio']:9.4f}  "
              f"{r['frac_ring']:10.3f}  {r['sigma_r']:7.4f}  {r['ratio']:7.4f}  "
              f"{Phi:7.4f}  {phase:>10s}")


# ============================================================
# Visualization
# ============================================================

def create_visualization(results, density_data, angular_data, spacing_data,
                          crit_data, landscape_data, cross_data):
    """Create comprehensive 12-panel visualization."""
    fig = plt.figure(figsize=(28, 21))
    fig.suptitle("Priority 5: Density Transition in Ihara Zeta Poles\n"
                 "Phase Transition in Pole Distribution at h* ≈ 9.5",
                 fontsize=16, fontweight='bold')

    gs = GridSpec(3, 4, hspace=0.35, wspace=0.30)

    # ---- Panel (a): Pole flow — inner/rc vs h ----
    ax = fig.add_subplot(gs[0, 0])
    for r in results:
        if r['h'] <= 0 or len(r['radii']) < 5:
            continue
        color = '#2196F3' if r['is_ramanujan'] else '#F44336'
        marker = 'o'
        ax.plot(r['h'], r['inner_ratio'], marker, color=color,
                markersize=6, alpha=0.7)
    ax.axvline(x=9.5, color='gray', linestyle=':', alpha=0.5)
    ax.axhline(y=1.0, color='green', linestyle='--', alpha=0.3)
    ax.set_xlabel('Coxeter number h')
    ax.set_ylabel('inner/rc (min|u_nt| / (1/√q))')
    ax.set_title('(a) Inner pole migration vs h')
    ax.grid(True, alpha=0.2)

    # ---- Panel (b): frac_in_ring vs h ----
    ax = fig.add_subplot(gs[0, 1])
    for r in results:
        if r['h'] <= 0 or len(r['radii']) < 5:
            continue
        color = '#2196F3' if r['is_ramanujan'] else '#F44336'
        ax.plot(r['h'], r['frac_ring'], 'o', color=color,
                markersize=6, alpha=0.7)
    ax.axvline(x=9.5, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Coxeter number h')
    ax.set_ylabel('Fraction in ±5% ring')
    ax.set_title('(b) Pole concentration vs h')
    ax.grid(True, alpha=0.2)

    # ---- Panel (c): σ(r) vs h ----
    ax = fig.add_subplot(gs[0, 2])
    for r in results:
        if r['h'] <= 0 or len(r['radii']) < 5:
            continue
        color = '#2196F3' if r['is_ramanujan'] else '#F44336'
        ax.plot(r['h'], r['sigma_r'], 'o', color=color,
                markersize=6, alpha=0.7)
    ax.axvline(x=9.5, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Coxeter number h')
    ax.set_ylabel('σ(|u|) pole dispersion')
    ax.set_title('(c) Pole dispersion vs h')
    ax.grid(True, alpha=0.2)

    # ---- Panel (d): Radial density evolution (A_n family) ----
    ax = fig.add_subplot(gs[0, 3])
    cmap = plt.cm.viridis
    a_names_for_density = ['su(3)', 'su(5)', 'su(7)', 'su(9)', 'su(12)']
    for i, name in enumerate(a_names_for_density):
        if name in density_data:
            dd = density_data[name]
            c = cmap(i / len(a_names_for_density))
            ax.plot(dd['centers'], dd['density'], '-', color=c,
                    linewidth=1.5, label=name, alpha=0.8)
    ax.axvline(x=1.0, color='red', linestyle='--', alpha=0.3, label='r/rc=1')
    ax.set_xlabel('r / rc (normalized radius)')
    ax.set_ylabel('Density ρ(r/rc)')
    ax.set_title('(d) Radial density evolution (A_n)')
    ax.legend(fontsize=7, ncol=2)
    ax.set_xlim(0, 2.5)
    ax.grid(True, alpha=0.2)

    # ---- Panel (e): Z(u) landscape for su(3) ----
    if 'su(3)' in landscape_data:
        ax = fig.add_subplot(gs[1, 0])
        ld = landscape_data['su(3)']
        vmax = np.percentile(ld['log_Z_inv'][np.isfinite(ld['log_Z_inv'])], 95)
        vmin = np.percentile(ld['log_Z_inv'][np.isfinite(ld['log_Z_inv'])], 5)
        ax.contourf(ld['RE'], ld['IM'], ld['log_Z_inv'], levels=40,
                     cmap='RdYlBu', vmin=vmin, vmax=vmax)
        theta = np.linspace(0, 2*np.pi, 200)
        ax.plot(ld['rc']*np.cos(theta), ld['rc']*np.sin(theta),
                'k--', linewidth=1.5, label=f'1/√q')
        if len(ld['poles']) > 0:
            ax.plot(ld['poles'].real, ld['poles'].imag, 'k.', markersize=2, alpha=0.5)
        ax.set_aspect('equal')
        ax.set_title('(e) log|Z⁻¹(u)| — su(3) [Ram]')
        ax.set_xlabel('Re(u)')
        ax.set_ylabel('Im(u)')

    # ---- Panel (f): Z(u) landscape for su(8) (transition) ----
    if 'su(8)' in landscape_data:
        ax = fig.add_subplot(gs[1, 1])
        ld = landscape_data['su(8)']
        vmax = np.percentile(ld['log_Z_inv'][np.isfinite(ld['log_Z_inv'])], 95)
        vmin = np.percentile(ld['log_Z_inv'][np.isfinite(ld['log_Z_inv'])], 5)
        ax.contourf(ld['RE'], ld['IM'], ld['log_Z_inv'], levels=40,
                     cmap='RdYlBu', vmin=vmin, vmax=vmax)
        theta = np.linspace(0, 2*np.pi, 200)
        ax.plot(ld['rc']*np.cos(theta), ld['rc']*np.sin(theta),
                'k--', linewidth=1.5)
        if len(ld['poles']) > 0:
            ax.plot(ld['poles'].real, ld['poles'].imag, 'k.', markersize=2, alpha=0.5)
        ax.set_aspect('equal')
        ax.set_title('(f) log|Z⁻¹(u)| — su(8) [transition]')
        ax.set_xlabel('Re(u)')
        ax.set_ylabel('Im(u)')

    # ---- Panel (g): Z(u) landscape for su(12) ----
    if 'su(12)' in landscape_data:
        ax = fig.add_subplot(gs[1, 2])
        ld = landscape_data['su(12)']
        vmax = np.percentile(ld['log_Z_inv'][np.isfinite(ld['log_Z_inv'])], 95)
        vmin = np.percentile(ld['log_Z_inv'][np.isfinite(ld['log_Z_inv'])], 5)
        ax.contourf(ld['RE'], ld['IM'], ld['log_Z_inv'], levels=40,
                     cmap='RdYlBu', vmin=vmin, vmax=vmax)
        theta = np.linspace(0, 2*np.pi, 200)
        ax.plot(ld['rc']*np.cos(theta), ld['rc']*np.sin(theta),
                'k--', linewidth=1.5)
        if len(ld['poles']) > 0:
            ax.plot(ld['poles'].real, ld['poles'].imag, 'k.', markersize=2, alpha=0.5)
        ax.set_aspect('equal')
        ax.set_title('(g) log|Z⁻¹(u)| — su(12) [non-Ram]')
        ax.set_xlabel('Re(u)')
        ax.set_ylabel('Im(u)')

    # ---- Panel (h): Spacing statistics comparison ----
    ax = fig.add_subplot(gs[1, 3])
    for name in ['su(3)', 'su(7)', 'su(8)', 'su(12)']:
        if name in spacing_data:
            sd = spacing_data[name]
            bins = np.linspace(0, 4, 30)
            ax.hist(sd['spacings'], bins=bins, density=True, alpha=0.4, label=name)
    # Overlay Poisson
    s_th = np.linspace(0, 4, 100)
    ax.plot(s_th, np.exp(-s_th), 'k--', linewidth=1.5, label='Poisson', alpha=0.7)
    # Overlay GUE
    ax.plot(s_th, (32/np.pi**2) * s_th**2 * np.exp(-4*s_th**2/np.pi),
            'r-', linewidth=1.5, label='GUE', alpha=0.7)
    ax.set_xlabel('Normalized spacing s')
    ax.set_ylabel('P(s)')
    ax.set_title('(h) Nearest-neighbor spacing distribution')
    ax.legend(fontsize=7)
    ax.set_xlim(0, 4)
    ax.grid(True, alpha=0.2)

    # ---- Panel (i): Cross-family inner/rc comparison ----
    ax = fig.add_subplot(gs[2, 0])
    markers = {'A_n (su)': 'o', 'B_n (so_odd)': '^', 'C_n (sp)': 's', 'D_n (so_even)': 'D'}
    colors_fam = {'A_n (su)': '#2196F3', 'B_n (so_odd)': '#4CAF50',
                   'C_n (sp)': '#FF9800', 'D_n (so_even)': '#9C27B0'}
    for fam_name, fam_data in cross_data.items():
        hs_f = [r['h'] for r in fam_data if r['h'] > 0]
        ir_f = [r['inner_ratio'] for r in fam_data if r['h'] > 0]
        ax.plot(hs_f, ir_f, markers.get(fam_name, 'o') + '-',
                color=colors_fam.get(fam_name, 'gray'),
                markersize=5, label=fam_name, alpha=0.8)
    ax.axvline(x=9.5, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Coxeter number h')
    ax.set_ylabel('inner/rc')
    ax.set_title('(i) Cross-family pole migration')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.2)

    # ---- Panel (j): Cross-family frac_ring comparison ----
    ax = fig.add_subplot(gs[2, 1])
    for fam_name, fam_data in cross_data.items():
        hs_f = [r['h'] for r in fam_data if r['h'] > 0]
        fr_f = [r['frac_ring'] for r in fam_data if r['h'] > 0]
        ax.plot(hs_f, fr_f, markers.get(fam_name, 'o') + '-',
                color=colors_fam.get(fam_name, 'gray'),
                markersize=5, label=fam_name, alpha=0.8)
    ax.axvline(x=9.5, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Coxeter number h')
    ax.set_ylabel('Fraction in ±5% ring')
    ax.set_title('(j) Cross-family pole concentration')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.2)

    # ---- Panel (k): Angular uniformity R̄ vs h ----
    ax = fig.add_subplot(gs[2, 2])
    for r in results:
        if r['name'] in angular_data and r['h'] > 0:
            ad = angular_data[r['name']]
            color = '#2196F3' if r['is_ramanujan'] else '#F44336'
            ax.plot(r['h'], ad['R_bar'], 'o', color=color,
                    markersize=6, alpha=0.7)
    ax.axvline(x=9.5, color='gray', linestyle=':', alpha=0.5)
    ax.axhline(y=0.0, color='green', linestyle='--', alpha=0.3)
    ax.set_xlabel('Coxeter number h')
    ax.set_ylabel('Rayleigh R̄ (0=uniform)')
    ax.set_title('(k) Angular uniformity vs h')
    ax.grid(True, alpha=0.2)

    # ---- Panel (l): Unified phase diagram ----
    ax = fig.add_subplot(gs[2, 3])
    hs_all = crit_data['hs']
    ir_all = crit_data['inner_ratios']
    fr_all = crit_data['frac_rings']
    sg_all = crit_data['sigmas']

    ax.plot(hs_all, ir_all, 'o-', color='#2196F3', markersize=4,
            label='inner/rc', alpha=0.7)
    ax.plot(hs_all, fr_all, 's-', color='#4CAF50', markersize=4,
            label='frac_ring', alpha=0.7)
    ax.plot(hs_all, sg_all * 20, '^-', color='#FF9800', markersize=4,
            label='σ(r) × 20', alpha=0.7)
    ax.axvline(x=9.5, color='red', linestyle=':', linewidth=2,
               alpha=0.7, label='h* ≈ 9.5')
    ax.fill_betweenx([0, 1.5], 0, 9.5, alpha=0.05, color='blue')
    ax.fill_betweenx([0, 1.5], 9.5, 35, alpha=0.05, color='red')
    ax.text(5, 1.35, 'Ramanujan\n(ordered)', ha='center', fontsize=9,
            color='blue', alpha=0.6)
    ax.text(20, 1.35, 'Non-Ramanujan\n(disordered)', ha='center', fontsize=9,
            color='red', alpha=0.6)
    ax.set_xlabel('Coxeter number h')
    ax.set_ylabel('Order parameter')
    ax.set_title('(l) Unified phase diagram')
    ax.legend(fontsize=7, loc='center right')
    ax.set_ylim(0, 1.5)
    ax.grid(True, alpha=0.2)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#2196F3',
               markersize=8, label='Ramanujan'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#F44336',
               markersize=8, label='Non-Ramanujan'),
        Line2D([0], [0], color='gray', linestyle=':', label='h* ≈ 9.5'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3,
               fontsize=10, bbox_to_anchor=(0.5, -0.01))

    outpath = os.path.join(os.path.dirname(__file__), 'pole_density_transition.png')
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\n  [Saved figure: {outpath}]")
    plt.close(fig)


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 90)
    print("  PRIORITY 5: Density Transition in Ihara Zeta Poles")
    print("  Tracking crystallization of pole structure across h* ≈ 9.5")
    print("=" * 90)
    print()

    # Compute poles for all algebras
    catalog = build_catalog()
    results = []
    print("  Computing poles for all 37 algebras...")
    for i, (name, rank, all_roots, simple) in enumerate(catalog):
        data = compute_pole_data(name, rank, all_roots, simple)
        results.append(data)
        if (i + 1) % 10 == 0:
            print(f"    ... {i+1}/{len(catalog)} done")
    print(f"    ... {len(catalog)}/{len(catalog)} done")
    print()

    # Run all analyses
    a_results = pole_flow_analysis(results)
    density_data = radial_density_analysis(results)
    angular_data = angular_analysis(results)
    spacing_data = spacing_analysis(results)
    crit_data = critical_scaling(results)
    landscape_data = zeta_landscape(results)
    cross_data = cross_family_analysis(results)
    phase_diagram(results, crit_data)

    # Visualization
    create_visualization(results, density_data, angular_data, spacing_data,
                          crit_data, landscape_data, cross_data)

    print()
    print("  DONE. Density transition analysis complete.")


if __name__ == '__main__':
    main()
