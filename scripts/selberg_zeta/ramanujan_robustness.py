"""
P6: Robustness of the Ramanujan Property Under Alternative Definitions
======================================================================

The adjoint commutation graphs of Lie algebras are IRREGULAR. The standard
Ramanujan property (|λ_nt| ≤ 2√(d-1)) is defined for d-regular graphs.
For irregular graphs, multiple definitions exist:

  D1) MEAN-DEGREE (current):  |λ_nt| ≤ 2√(d̄ - 1)            [d̄ = mean degree]
  D2) MAX-DEGREE (Li–Solé):   |λ_nt| ≤ 2√(d_max - 1)         [most permissive]
  D3) GEOMETRIC-MEAN:         |λ_nt| ≤ 2·(Π(d_i-1))^(1/2n)   [degree product]
  D4) HASHIMOTO (ρ(T_e)):     |λ_nt| ≤ 2√ρ(T_e)              [non-backtracking]
  D5) ALON-BOPPANA:           |λ_nt| ≤ 2√(d̃ - 1)            [d̃ = harmonic mean]
  D6) FRIEDMAN (d₂):          |λ_nt| ≤ 2√(d₂ - 1)            [2nd largest degree]

The question: ARE THE P1–P5 CONCLUSIONS ROBUST across all definitions?
Specifically: Which algebras flip Ramanujan ↔ non-Ramanujan?

Usage:
    python ramanujan_robustness.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from ramanujan_extended import (
    roots_A, roots_B, roots_C, roots_D,
    roots_G2, roots_F4, roots_E6, roots_E7, roots_E8,
    adjacency_from_roots, spectral_analysis
)
from ihara_zeta_adjoint import ihara_poles_hashimoto


# ==============================================================
# Build catalog (same 37 algebras)
# ==============================================================

def build_catalog():
    catalog = []
    for n in range(1, 12):
        r, roots, simple = roots_A(n)
        catalog.append((f"su({n+1})", r, roots, simple, 'A'))
    for n in range(2, 9):
        r, roots, simple = roots_B(n)
        catalog.append((f"so({2*n+1})", r, roots, simple, 'B'))
    for n in range(2, 9):
        r, roots, simple = roots_C(n)
        catalog.append((f"sp({2*n})", r, roots, simple, 'C'))
    for n in range(3, 10):
        r, roots, simple = roots_D(n)
        catalog.append((f"so({2*n})", r, roots, simple, 'D'))
    for exc_name, fn, fam in [('G2', roots_G2, 'G'), ('F4', roots_F4, 'F'),
                               ('E6', roots_E6, 'E'), ('E7', roots_E7, 'E'),
                               ('E8', roots_E8, 'E')]:
        r, roots, simple = fn()
        catalog.append((exc_name, r, roots, simple, fam))
    return catalog


# ==============================================================
# Six Ramanujan definitions
# ==============================================================

def compute_all_definitions(adj, hash_eigs=None):
    """Compute Ramanujan bounds under 6 definitions for a given adjacency matrix."""
    n = adj.shape[0]
    eigs = np.linalg.eigvalsh(adj)
    eigs = np.sort(eigs)[::-1]

    degrees = adj.sum(axis=1)
    d_mean = np.mean(degrees)
    d_max = np.max(degrees)
    d_min = np.min(degrees)

    # Harmonic mean of degrees
    d_harm = n / np.sum(1.0 / degrees) if np.all(degrees > 0) else d_mean

    # Geometric mean of (d_i - 1)
    dm1 = degrees - 1
    dm1_pos = dm1[dm1 > 0]
    if len(dm1_pos) == n:
        d_geom_m1 = np.exp(np.mean(np.log(dm1_pos)))
    else:
        d_geom_m1 = d_mean - 1

    # Second-largest degree (Friedman)
    unique_degs = np.unique(degrees)[::-1]
    d_2nd = unique_degs[1] if len(unique_degs) > 1 else unique_degs[0]

    # Nontrivial eigenvalues
    is_regular = np.allclose(degrees, degrees[0])
    nontrivial = eigs[1:]
    if n > 2 and is_regular and abs(eigs[-1] + eigs[0]) < 0.01:
        nontrivial = eigs[1:-1]
    max_nt = np.max(np.abs(nontrivial)) if len(nontrivial) > 0 else 0.0

    # Hashimoto spectral radius
    rho_hash = np.max(np.abs(hash_eigs)) if hash_eigs is not None else None

    definitions = {}

    # D1: Mean-degree
    q1 = d_mean - 1
    bound1 = 2 * np.sqrt(max(q1, 0))
    definitions['D1_mean'] = {
        'q': q1, 'bound': bound1,
        'is_ram': bool(max_nt <= bound1 + 1e-10),
        'ratio': max_nt / bound1 if bound1 > 0 else float('inf'),
        'label': f'd̄={d_mean:.2f}'
    }

    # D2: Max-degree (Li-Solé) — most permissive
    q2 = d_max - 1
    bound2 = 2 * np.sqrt(max(q2, 0))
    definitions['D2_dmax'] = {
        'q': q2, 'bound': bound2,
        'is_ram': bool(max_nt <= bound2 + 1e-10),
        'ratio': max_nt / bound2 if bound2 > 0 else float('inf'),
        'label': f'd_max={d_max:.0f}'
    }

    # D3: Geometric mean of (d_i - 1)
    bound3 = 2 * np.sqrt(max(d_geom_m1, 0))
    definitions['D3_geom'] = {
        'q': d_geom_m1, 'bound': bound3,
        'is_ram': bool(max_nt <= bound3 + 1e-10),
        'ratio': max_nt / bound3 if bound3 > 0 else float('inf'),
        'label': f'q_geom={d_geom_m1:.2f}'
    }

    # D4: Hashimoto spectral radius — bound mapped to adjacency: 2√ρ(T_e)
    # For regular graphs ρ(T_e) = q = d-1, so 2√ρ = 2√(d-1) = standard bound
    if rho_hash is not None:
        bound4 = 2 * np.sqrt(rho_hash)
        definitions['D4_hash'] = {
            'q': rho_hash, 'bound': bound4,
            'is_ram': bool(max_nt <= bound4 + 1e-10),
            'ratio': max_nt / bound4 if bound4 > 0 else float('inf'),
            'label': f'ρ(T_e)={rho_hash:.2f}'
        }
    else:
        definitions['D4_hash'] = {
            'q': None, 'bound': None,
            'is_ram': None, 'ratio': None,
            'label': 'N/A (too large)'
        }

    # D5: Harmonic mean (Alon-Boppana style)
    q5 = d_harm - 1
    bound5 = 2 * np.sqrt(max(q5, 0))
    definitions['D5_harm'] = {
        'q': q5, 'bound': bound5,
        'is_ram': bool(max_nt <= bound5 + 1e-10),
        'ratio': max_nt / bound5 if bound5 > 0 else float('inf'),
        'label': f'd_harm={d_harm:.2f}'
    }

    # D6: Second-largest degree (Friedman-style)
    q6 = d_2nd - 1
    bound6 = 2 * np.sqrt(max(q6, 0))
    definitions['D6_d2'] = {
        'q': q6, 'bound': bound6,
        'is_ram': bool(max_nt <= bound6 + 1e-10),
        'ratio': max_nt / bound6 if bound6 > 0 else float('inf'),
        'label': f'd₂={d_2nd:.0f}'
    }

    return {
        'max_nt': max_nt,
        'd_mean': d_mean,
        'd_max': d_max,
        'd_min': d_min,
        'd_harm': d_harm,
        'd_geom_m1': d_geom_m1,
        'd_2nd': d_2nd,
        'is_regular': is_regular,
        'heterogeneity': d_max / d_min if d_min > 0 else float('inf'),
        'definitions': definitions
    }


# ==============================================================
# Main analysis
# ==============================================================

def run_analysis():
    print("=" * 120)
    print("  P6: ROBUSTNESS OF RAMANUJAN PROPERTY UNDER ALTERNATIVE DEFINITIONS")
    print("  FOR IRREGULAR ADJOINT COMMUTATION GRAPHS")
    print("=" * 120)
    print()

    catalog = build_catalog()
    results = []

    # Hashimoto dimension threshold (m = 2*|E| directed edges, T_e is m×m)
    HASH_DIM_LIMIT = 200  # adjacency dim limit for Hashimoto computation

    for name, rank, roots, simple, family in catalog:
        adj = adjacency_from_roots(rank, roots, simple)
        dim = adj.shape[0]

        # Hashimoto: only for algebras with dim ≤ limit
        hash_eigs = None
        if dim <= HASH_DIM_LIMIT:
            try:
                _, h_eigs = ihara_poles_hashimoto(adj)
                hash_eigs = h_eigs
            except Exception:
                pass

        defs = compute_all_definitions(adj, hash_eigs)
        coxeter = _coxeter_number(name)
        results.append({
            'name': name, 'dim': dim, 'family': family,
            'coxeter': coxeter,
            **defs
        })

    # ============================
    # TABLE 1: Full comparison
    # ============================
    print("─" * 120)
    print(f"  {'Algebra':>10s}  {'dim':>4s}  {'h':>3s}  "
          f"{'d̄':>6s}  {'d_max':>5s}  {'d_min':>5s}  {'het':>5s}  "
          f"{'|λ_nt|':>7s}  "
          f"{'D1':>4s}  {'D2':>4s}  {'D3':>4s}  {'D4':>4s}  {'D5':>4s}  {'D6':>4s}  "
          f"{'flip?':>5s}")
    print("─" * 120)

    flip_count = 0
    all_agree = 0
    for r in results:
        d = r['definitions']
        statuses = []
        for key in ['D1_mean', 'D2_dmax', 'D3_geom', 'D4_hash', 'D5_harm', 'D6_d2']:
            val = d[key]['is_ram']
            if val is None:
                statuses.append('--')
            elif val:
                statuses.append('✓')
            else:
                statuses.append('✗')

        # Check for flips: any definition disagrees with D1
        d1_val = d['D1_mean']['is_ram']
        has_flip = False
        for key in ['D2_dmax', 'D3_geom', 'D4_hash', 'D5_harm', 'D6_d2']:
            v = d[key]['is_ram']
            if v is not None and v != d1_val:
                has_flip = True
                break

        if has_flip:
            flip_count += 1
            flip_str = "⚠ FLIP"
        else:
            all_agree += 1
            flip_str = ""

        print(f"  {r['name']:>10s}  {r['dim']:4d}  {r['coxeter']:3d}  "
              f"{r['d_mean']:6.1f}  {r['d_max']:5.0f}  {r['d_min']:5.0f}  "
              f"{r['heterogeneity']:5.2f}  "
              f"{r['max_nt']:7.3f}  "
              f"{'  '.join(statuses)}  "
              f"{flip_str}")

    total = len(results)
    print("─" * 120)
    print(f"\n  SUMMARY: {all_agree}/{total} algebras have unanimous classification, "
          f"{flip_count}/{total} have at least one flip")
    print()

    # ============================
    # TABLE 2: Ratio comparison
    # ============================
    print("=" * 120)
    print("  RATIOS: |λ_nt| / bound  (< 1 = Ramanujan)")
    print("=" * 120)
    print(f"  {'Algebra':>10s}  {'dim':>4s}  {'h':>3s}  "
          f"{'R(D1)':>7s}  {'R(D2)':>7s}  {'R(D3)':>7s}  "
          f"{'R(D4)':>7s}  {'R(D5)':>7s}  {'R(D6)':>7s}  "
          f"{'spread':>7s}")
    print("─" * 120)

    for r in results:
        d = r['definitions']
        ratios = []
        ratio_strs = []
        for key in ['D1_mean', 'D2_dmax', 'D3_geom', 'D4_hash', 'D5_harm', 'D6_d2']:
            v = d[key]['ratio']
            if v is not None:
                ratios.append(v)
                ratio_strs.append(f"{v:7.4f}")
            else:
                ratio_strs.append(f"{'--':>7s}")

        spread = max(ratios) - min(ratios) if len(ratios) > 1 else 0.0
        print(f"  {r['name']:>10s}  {r['dim']:4d}  {r['coxeter']:3d}  "
              f"{'  '.join(ratio_strs)}  "
              f"{spread:7.4f}")

    # ============================
    # TABLE 3: Classification by definition
    # ============================
    print()
    print("=" * 120)
    print("  CLASSIFICATION COUNTS")
    print("=" * 120)
    def_names = ['D1_mean', 'D2_dmax', 'D3_geom', 'D4_hash', 'D5_harm', 'D6_d2']
    labels = ['D1 (d̄)', 'D2 (d_max)', 'D3 (geom)', 'D4 (Hash)', 'D5 (harm)', 'D6 (d₂)']

    for dname, label in zip(def_names, labels):
        n_ram = sum(1 for r in results if r['definitions'][dname]['is_ram'] is True)
        n_nram = sum(1 for r in results if r['definitions'][dname]['is_ram'] is False)
        n_na = sum(1 for r in results if r['definitions'][dname]['is_ram'] is None)
        na_str = f" ({n_na} N/A)" if n_na > 0 else ""
        print(f"  {label:>15s}: {n_ram:2d} Ramanujan, {n_nram:2d} non-Ramanujan{na_str}")

    # ============================
    # ANALYSIS 1: Flip algebras
    # ============================
    print()
    print("=" * 120)
    print("  FLIP ANALYSIS: Algebras whose classification changes between definitions")
    print("=" * 120)
    print()

    flip_algebras = []
    for r in results:
        d = r['definitions']
        d1_val = d['D1_mean']['is_ram']
        statuses = {}
        for key in def_names:
            v = d[key]['is_ram']
            if v is not None:
                statuses[key] = v

        unique_vals = set(statuses.values())
        if len(unique_vals) > 1:
            flip_algebras.append(r)
            print(f"  {r['name']} (dim={r['dim']}, h={r['coxeter']}):")
            for key, label in zip(def_names, labels):
                v = d[key]
                if v['is_ram'] is not None:
                    sym = "✓ RAM" if v['is_ram'] else "✗ non"
                    print(f"    {label:>15s}: {sym}  ratio={v['ratio']:.4f}  "
                          f"bound={v['bound']:.4f}  ({v['label']})")
            print()

    if not flip_algebras:
        print("  NO FLIPS — all definitions agree for all 37 algebras!")
    else:
        print(f"  Total flips: {len(flip_algebras)} algebras")

    # ============================
    # ANALYSIS 2: Degree heterogeneity structure
    # ============================
    print()
    print("=" * 120)
    print("  DEGREE HETEROGENEITY OF ADJOINT GRAPHS")
    print("=" * 120)
    print()
    print(f"  {'Algebra':>10s}  {'dim':>4s}  {'d̄':>7s}  {'d_max':>5s}  "
          f"{'d_min':>5s}  {'d_harm':>7s}  {'q_geom':>7s}  "
          f"{'het(max/min)':>11s}  {'gap(max-mean)':>13s}")
    print("─" * 100)

    for r in results:
        gap = r['d_max'] - r['d_mean']
        print(f"  {r['name']:>10s}  {r['dim']:4d}  {r['d_mean']:7.2f}  "
              f"{r['d_max']:5.0f}  {r['d_min']:5.0f}  {r['d_harm']:7.2f}  "
              f"{r['d_geom_m1']:7.2f}  "
              f"{r['heterogeneity']:11.3f}  {gap:13.2f}")

    # ============================
    # ANALYSIS 3: Ordering of bounds
    # ============================
    print()
    print("=" * 120)
    print("  BOUND ORDERING: Which definition is strictest/loosest?")
    print("=" * 120)
    print()
    print("  Mathematical ordering (for irregular graphs):")
    print("    d̄ ≤ d_max  →  2√(d̄-1) ≤ 2√(d_max-1)  →  D1 is STRICTER than D2")
    print("    d_harm ≤ d̄  →  D5 is STRICTER than D1")
    print("    Generally: D5(harm) ≤ D1(mean) ≤ D3(geom) ≤ D2(dmax)")
    print("    D4(Hash) and D6(d₂) depend on graph structure")
    print()

    # Verify empirically
    n_d5_strictest = 0
    n_d2_loosest = 0
    for r in results:
        d = r['definitions']
        bounds = {}
        for key in def_names:
            if d[key]['bound'] is not None:
                bounds[key] = d[key]['bound']
        if len(bounds) < 3:
            continue
        min_key = min(bounds, key=bounds.get)
        max_key = max(bounds, key=bounds.get)
        if min_key == 'D5_harm':
            n_d5_strictest += 1
        if max_key == 'D2_dmax':
            n_d2_loosest += 1

    print(f"  Empirical: D5(harm) is strictest for {n_d5_strictest}/{total} algebras")
    print(f"  Empirical: D2(dmax) is loosest for {n_d2_loosest}/{total} algebras")

    # ============================
    # ANALYSIS 4: Critical transition h* under each definition
    # ============================
    print()
    print("=" * 120)
    print("  CRITICAL COXETER NUMBER h* UNDER EACH DEFINITION")
    print("=" * 120)
    print()

    # For A_n family: find where Ramanujan → non-Ramanujan under each def
    a_results = [r for r in results if r['family'] == 'A']
    a_results.sort(key=lambda r: r['coxeter'])

    print(f"  A_n family transition (su(n+1)):")
    for dname, label in zip(def_names, labels):
        last_ram = None
        first_nram = None
        for r in a_results:
            v = r['definitions'][dname]['is_ram']
            if v is True:
                last_ram = r
            elif v is False and first_nram is None:
                first_nram = r
        if last_ram and first_nram:
            print(f"    {label:>15s}: last Ramanujan = {last_ram['name']} (h={last_ram['coxeter']}), "
                  f"first non-Ram = {first_nram['name']} (h={first_nram['coxeter']})  →  "
                  f"h* ∈ [{last_ram['coxeter']}, {first_nram['coxeter']}]")
        elif first_nram is None:
            print(f"    {label:>15s}: ALL Ramanujan in A_n range (h up to {a_results[-1]['coxeter']})")
        elif last_ram is None:
            print(f"    {label:>15s}: NONE Ramanujan in A_n range")

    # For all families
    print()
    print(f"  Cross-family transition:")
    for dname, label in zip(def_names, labels):
        ram_h = [r['coxeter'] for r in results if r['definitions'][dname]['is_ram'] is True]
        nram_h = [r['coxeter'] for r in results if r['definitions'][dname]['is_ram'] is False]
        if ram_h and nram_h:
            max_ram_h = max(ram_h)
            min_nram_h = min(nram_h)
            print(f"    {label:>15s}: max h(Ram)={max_ram_h:3d}, "
                  f"min h(non-Ram)={min_nram_h:3d}, "
                  f"overlap={'YES' if min_nram_h <= max_ram_h else 'NO'}")

    # ============================
    # ANALYSIS 5: P1-P5 conclusions robustness check
    # ============================
    print()
    print("=" * 120)
    print("  P1-P5 CONCLUSIONS ROBUSTNESS CHECK")
    print("=" * 120)
    print()

    sm_algebras = ['su(2)', 'su(3)']
    gut_algebras = ['su(5)', 'so(10)']
    problem_algebras = ['E6', 'E7', 'E8']

    for group_name, group in [("Standard Model", sm_algebras),
                               ("GUT candidates", gut_algebras),
                               ("Problem cases", problem_algebras)]:
        print(f"  {group_name}:")
        for aname in group:
            r = next((x for x in results if x['name'] == aname), None)
            if r is None:
                continue
            d = r['definitions']
            all_ram = all(d[k]['is_ram'] is True for k in def_names if d[k]['is_ram'] is not None)
            all_nram = all(d[k]['is_ram'] is False for k in def_names if d[k]['is_ram'] is not None)
            if all_ram:
                verdict = "ROBUST Ramanujan (all definitions)"
            elif all_nram:
                verdict = "ROBUST non-Ramanujan (all definitions)"
            else:
                ram_defs = [l for k, l in zip(def_names, labels) if d[k]['is_ram'] is True]
                nram_defs = [l for k, l in zip(def_names, labels) if d[k]['is_ram'] is False]
                verdict = (f"MIXED: Ram by {','.join(ram_defs)}; "
                          f"non-Ram by {','.join(nram_defs)}")
            print(f"    {aname:>10s}: {verdict}")
        print()

    # Key conclusions from P1-P5
    checks = [
        ("SM algebras are Ramanujan",
         all(all(results[i]['definitions'][k]['is_ram'] is True
                 for k in def_names if results[i]['definitions'][k]['is_ram'] is not None)
             for i, r in enumerate(results) if r['name'] in sm_algebras)),
        ("E₈ is NOT Ramanujan",
         all(results[i]['definitions'][k]['is_ram'] is False
             for i, r in enumerate(results) if r['name'] == 'E8'
             for k in def_names if results[i]['definitions'][k]['is_ram'] is not None)),
        ("Transition exists in A_n family",
         any(results[i]['definitions']['D1_mean']['is_ram'] is True for i, r in enumerate(results)
             if r['family'] == 'A') and
         any(results[i]['definitions']['D1_mean']['is_ram'] is False for i, r in enumerate(results)
             if r['family'] == 'A')),
    ]

    print("  KEY CONCLUSION CHECKS:")
    for desc, holds in checks:
        print(f"    {'✓' if holds else '✗'} {desc}")

    # ============================
    # ANALYSIS 6: Sensitivity band — algebras near the boundary
    # ============================
    print()
    print("=" * 120)
    print("  SENSITIVITY BAND: Algebras with ratio ∈ [0.85, 1.15] under any definition")
    print("=" * 120)
    print()

    for r in results:
        d = r['definitions']
        near_boundary = []
        for key, label in zip(def_names, labels):
            v = d[key]['ratio']
            if v is not None and 0.85 <= v <= 1.15:
                near_boundary.append((label, v))
        if near_boundary:
            print(f"  {r['name']:>10s} (h={r['coxeter']:2d}): "
                  + ", ".join(f"{l}: {v:.4f}" for l, v in near_boundary))

    # ============================
    # Visualization
    # ============================
    create_visualization(results, def_names, labels)

    return results


def _coxeter_number(name):
    """Return Coxeter number for a simple Lie algebra."""
    coxeter_map = {
        'su(2)': 2, 'su(3)': 3, 'su(4)': 4, 'su(5)': 5, 'su(6)': 6,
        'su(7)': 7, 'su(8)': 8, 'su(9)': 9, 'su(10)': 10, 'su(11)': 11,
        'su(12)': 12,
        'so(5)': 4, 'so(7)': 6, 'so(9)': 8, 'so(11)': 10, 'so(13)': 12,
        'so(15)': 14, 'so(17)': 16,
        'sp(4)': 4, 'sp(6)': 6, 'sp(8)': 8, 'sp(10)': 10, 'sp(12)': 12,
        'sp(14)': 14, 'sp(16)': 16,
        'so(6)': 4, 'so(8)': 6, 'so(10)': 8, 'so(12)': 10, 'so(14)': 12,
        'so(16)': 14, 'so(18)': 16,
        'G2': 6, 'F4': 12, 'E6': 12, 'E7': 18, 'E8': 30
    }
    return coxeter_map.get(name, 0)


def create_visualization(results, def_names, labels):
    """Create 9-panel visualization of robustness analysis."""
    fig = plt.figure(figsize=(24, 20))
    fig.suptitle("P6: Robustness of Ramanujan Property — 6 Definitions for Irregular Graphs",
                 fontsize=16, fontweight='bold')
    gs = GridSpec(3, 3, hspace=0.35, wspace=0.30)

    colors_def = {
        'D1_mean': '#2196F3',
        'D2_dmax': '#4CAF50',
        'D3_geom': '#FF9800',
        'D4_hash': '#9C27B0',
        'D5_harm': '#F44336',
        'D6_d2':   '#795548',
    }
    short_labels = {
        'D1_mean': 'D1(d̄)',
        'D2_dmax': 'D2(d_max)',
        'D3_geom': 'D3(geom)',
        'D4_hash': 'D4(Hash)',
        'D5_harm': 'D5(harm)',
        'D6_d2':   'D6(d₂)',
    }

    # Panel 1: Ratio vs Coxeter h for all definitions
    ax1 = fig.add_subplot(gs[0, 0])
    for dname in def_names:
        h_vals = []
        r_vals = []
        for r in results:
            v = r['definitions'][dname]['ratio']
            if v is not None:
                h_vals.append(r['coxeter'])
                r_vals.append(v)
        ax1.scatter(h_vals, r_vals, c=colors_def[dname], s=20, alpha=0.7,
                   label=short_labels[dname])
    ax1.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, linewidth=2)
    ax1.set_xlabel('Coxeter number h')
    ax1.set_ylabel('Ratio |λ_nt|/bound')
    ax1.set_title('Ratio vs h (all definitions)')
    ax1.legend(fontsize=7, loc='upper left')
    ax1.grid(True, alpha=0.3)

    # Panel 2: Classification agreement heatmap
    ax2 = fig.add_subplot(gs[0, 1])
    alg_names = [r['name'] for r in results]
    data = np.zeros((len(results), len(def_names)))
    for i, r in enumerate(results):
        for j, dname in enumerate(def_names):
            v = r['definitions'][dname]['is_ram']
            if v is True:
                data[i, j] = 1.0
            elif v is False:
                data[i, j] = 0.0
            else:
                data[i, j] = 0.5  # N/A

    im = ax2.imshow(data, aspect='auto', cmap='RdYlGn', vmin=0, vmax=1,
                    interpolation='nearest')
    ax2.set_xticks(range(len(def_names)))
    ax2.set_xticklabels([short_labels[d] for d in def_names], rotation=45, fontsize=7)
    ax2.set_yticks(range(len(results)))
    ax2.set_yticklabels(alg_names, fontsize=5)
    ax2.set_title('Classification heatmap\n(green=Ram, red=non-Ram, yellow=N/A)')
    plt.colorbar(im, ax=ax2, shrink=0.5)

    # Panel 3: Bound values comparison for selected algebras
    ax3 = fig.add_subplot(gs[0, 2])
    selected = ['su(3)', 'su(5)', 'su(8)', 'F4', 'E6', 'E8']
    x_pos = np.arange(len(selected))
    width = 0.12
    for j, dname in enumerate(def_names):
        bounds = []
        for sname in selected:
            r = next((x for x in results if x['name'] == sname), None)
            if r and r['definitions'][dname]['bound'] is not None:
                bounds.append(r['definitions'][dname]['bound'])
            else:
                bounds.append(0)
        ax3.bar(x_pos + j * width, bounds, width, label=short_labels[dname],
               color=colors_def[dname], alpha=0.8)
    # Add |λ_nt| as horizontal markers
    for i, sname in enumerate(selected):
        r = next((x for x in results if x['name'] == sname), None)
        if r:
            ax3.plot([x_pos[i] - 0.05, x_pos[i] + 6 * width + 0.05],
                    [r['max_nt'], r['max_nt']], 'k--', linewidth=1.5)
    ax3.set_xticks(x_pos + 2.5 * width)
    ax3.set_xticklabels(selected, rotation=45)
    ax3.set_ylabel('Bound value / |λ_nt|')
    ax3.set_title('Bounds vs |λ_nt| (dashed)')
    ax3.legend(fontsize=6, loc='upper left')

    # Panel 4: Heterogeneity vs ratio (D1)
    ax4 = fig.add_subplot(gs[1, 0])
    het = [r['heterogeneity'] for r in results]
    ratio_d1 = [r['definitions']['D1_mean']['ratio'] for r in results]
    colors_ram = ['blue' if r['definitions']['D1_mean']['is_ram'] else 'red' for r in results]
    ax4.scatter(het, ratio_d1, c=colors_ram, s=30, alpha=0.7, edgecolors='black', linewidth=0.3)
    for r in results:
        if r['name'] in ['su(8)', 'F4', 'E6', 'E8', 'su(3)']:
            ax4.annotate(r['name'], (r['heterogeneity'], r['definitions']['D1_mean']['ratio']),
                        fontsize=6, ha='left')
    ax4.axhline(y=1.0, color='red', linestyle='--', alpha=0.5)
    ax4.set_xlabel('Degree heterogeneity (d_max / d_min)')
    ax4.set_ylabel('Ratio (D1, mean-degree)')
    ax4.set_title('Heterogeneity vs Ramanujan ratio')
    ax4.grid(True, alpha=0.3)

    # Panel 5: Ratio spread (max - min across definitions) vs h
    ax5 = fig.add_subplot(gs[1, 1])
    h_vals = [r['coxeter'] for r in results]
    spreads = []
    for r in results:
        ratios = [r['definitions'][k]['ratio'] for k in def_names
                  if r['definitions'][k]['ratio'] is not None]
        spreads.append(max(ratios) - min(ratios) if len(ratios) > 1 else 0)
    ax5.scatter(h_vals, spreads, c=colors_ram, s=30, alpha=0.7, edgecolors='black', linewidth=0.3)
    for r, sp in zip(results, spreads):
        if sp > 0.15 or r['name'] in ['F4', 'E8']:
            ax5.annotate(r['name'], (r['coxeter'], sp), fontsize=6, ha='left')
    ax5.set_xlabel('Coxeter number h')
    ax5.set_ylabel('Ratio spread (max - min)')
    ax5.set_title('Definition sensitivity vs h')
    ax5.grid(True, alpha=0.3)

    # Panel 6: A_n family transition under all definitions
    ax6 = fig.add_subplot(gs[1, 2])
    a_results = sorted([r for r in results if r['family'] == 'A'], key=lambda r: r['coxeter'])
    for dname in def_names:
        h_a = [r['coxeter'] for r in a_results]
        r_a = [r['definitions'][dname]['ratio'] for r in a_results
               if r['definitions'][dname]['ratio'] is not None]
        h_a_valid = [r['coxeter'] for r in a_results
                     if r['definitions'][dname]['ratio'] is not None]
        ax6.plot(h_a_valid, r_a, 'o-', color=colors_def[dname], markersize=4,
                label=short_labels[dname], alpha=0.8)
    ax6.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, linewidth=2)
    ax6.set_xlabel('Coxeter number h')
    ax6.set_ylabel('Ratio')
    ax6.set_title('A_n family: ratio vs h (all definitions)')
    ax6.legend(fontsize=7)
    ax6.grid(True, alpha=0.3)

    # Panel 7: Number of Ramanujan algebras by definition
    ax7 = fig.add_subplot(gs[2, 0])
    counts_ram = []
    for dname in def_names:
        n_ram = sum(1 for r in results if r['definitions'][dname]['is_ram'] is True)
        counts_ram.append(n_ram)
    bars = ax7.bar(range(len(def_names)), counts_ram,
                   color=[colors_def[d] for d in def_names])
    ax7.set_xticks(range(len(def_names)))
    ax7.set_xticklabels([short_labels[d] for d in def_names], rotation=45, fontsize=8)
    ax7.set_ylabel('# Ramanujan algebras')
    ax7.set_title(f'Ramanujan count / {len(results)} algebras')
    ax7.axhline(y=19, color='gray', linestyle=':', alpha=0.5, label='D1 baseline (19)')
    for b, c in zip(bars, counts_ram):
        ax7.text(b.get_x() + b.get_width()/2, b.get_height() + 0.3,
                str(c), ha='center', fontsize=10, fontweight='bold')
    ax7.legend(fontsize=7)

    # Panel 8: Degree distribution for key algebras
    ax8 = fig.add_subplot(gs[2, 1])
    for sname, color in [('su(3)', 'blue'), ('su(8)', 'orange'),
                          ('F4', 'green'), ('E8', 'red')]:
        entry = next((c for c in build_catalog() if c[0] == sname), None)
        if entry:
            adj = adjacency_from_roots(entry[1], entry[2], entry[3])
            degs = adj.sum(axis=1)
            unique_d, counts_d = np.unique(degs, return_counts=True)
            ax8.bar(unique_d + ({'su(3)': -0.3, 'su(8)': -0.1,
                                  'F4': 0.1, 'E8': 0.3}.get(sname, 0)),
                    counts_d / len(degs), width=0.2, alpha=0.8, label=sname,
                    color=color)
    ax8.set_xlabel('Degree')
    ax8.set_ylabel('Fraction of vertices')
    ax8.set_title('Degree distributions')
    ax8.legend(fontsize=8)
    ax8.grid(True, alpha=0.3)

    # Panel 9: SM/GUT/Exceptional comparison radar
    ax9 = fig.add_subplot(gs[2, 2])
    groups = {
        'SM: su(2),su(3)': ['su(2)', 'su(3)'],
        'GUT: su(5),so(10)': ['su(5)', 'so(10)'],
        'Exceptional: E₆,E₇,E₈': ['E6', 'E7', 'E8'],
    }
    x_grp = np.arange(len(groups))
    width = 0.12
    for j, dname in enumerate(def_names):
        means = []
        for gname, members in groups.items():
            ratios_g = []
            for m in members:
                r = next((x for x in results if x['name'] == m), None)
                if r and r['definitions'][dname]['ratio'] is not None:
                    ratios_g.append(r['definitions'][dname]['ratio'])
            means.append(np.mean(ratios_g) if ratios_g else 0)
        ax9.bar(x_grp + j * width, means, width, label=short_labels[dname],
               color=colors_def[dname], alpha=0.8)
    ax9.axhline(y=1.0, color='red', linestyle='--', alpha=0.5)
    ax9.set_xticks(x_grp + 2.5 * width)
    ax9.set_xticklabels(list(groups.keys()), fontsize=8)
    ax9.set_ylabel('Mean ratio')
    ax9.set_title('SM / GUT / Exceptional comparison')
    ax9.legend(fontsize=6, loc='upper left')

    plt.savefig(os.path.join(os.path.dirname(__file__), 'ramanujan_robustness.png'),
                dpi=150, bbox_inches='tight')
    print(f"\n  Figure saved: ramanujan_robustness.png")


if __name__ == '__main__':
    run_analysis()
