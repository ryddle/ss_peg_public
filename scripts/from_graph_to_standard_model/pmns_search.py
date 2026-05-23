#!/usr/bin/env python3
"""
pmns_search.py — Systematic search for PMNS matrix parameters from graph invariants
=====================================================================================

The PMNS (Pontecorvo-Maki-Nakagawa-Sakata) matrix describes neutrino flavor mixing.
It has 4 independent parameters in the standard Dirac parametrization:
    - 3 mixing angles: theta_12 (solar), theta_23 (atmospheric), theta_13 (reactor)
    - 1 CP-violating phase: delta_CP

Unlike the CKM matrix (small angles, strongly hierarchical), the PMNS matrix has
LARGE mixing angles (theta_12 ~ 33°, theta_23 ~ 49°, theta_13 ~ 8.5°).

We search for graph formulas matching these parameters and derived quantities
(|U_alpha_i| matrix elements, mass-squared ratios, Jarlskog invariant).

Building blocks: the same 5 spectral invariants from the SS-PEG program.

Author: SS-PEG program
Date: 2026
"""

import numpy as np
from itertools import product as cart_product
import time

np.random.seed(42)

# ============================================================
# §0. GRAPH BUILDING BLOCKS
# ============================================================

print("=" * 72)
print("PMNS MATRIX PARAMETER SEARCH FROM GRAPH INVARIANTS")
print("=" * 72)

# Core spectral invariants
R2   = 1 / 2
R3   = (np.sqrt(57) - 3) / (2 * np.sqrt(17))  # 0.55175
R_EE = 1 / np.sqrt(2)                           # 0.70711
hv3  = 3
hv2  = 2

blocks = {
    "R2":   R2,
    "R3":   R3,
    "REE":  R_EE,
    "hv3":  hv3,
    "hv2":  hv2,
}

print("\nBuilding blocks:")
for name, val in blocks.items():
    print(f"  {name:6s} = {val:.6f}  (ln = {np.log(val):.6f})")

# ============================================================
# §1. PMNS EXPERIMENTAL VALUES (NuFIT 5.2, Normal Ordering)
# ============================================================

print("\n" + "=" * 72)
print("§1. PMNS EXPERIMENTAL VALUES (NuFIT 5.2, Normal Ordering)")
print("=" * 72)

# Mixing angles (degrees -> radians)
# Source: NuFIT 5.2 (November 2022), Normal Ordering
# Reference: Esteban et al., JHEP 09 (2020) 178 [updated Nov 2022]
theta_12_deg = 33.41   # solar angle,      +0.75 / -0.72
theta_23_deg = 49.1    # atmospheric angle, +1.0  / -1.3
theta_13_deg = 8.54    # reactor angle,     +0.11 / -0.12
delta_CP_deg = 197.0   # CP phase,          +42   / -25

theta_12 = np.radians(theta_12_deg)
theta_23 = np.radians(theta_23_deg)
theta_13 = np.radians(theta_13_deg)
delta_CP = np.radians(delta_CP_deg)

sin_t12 = np.sin(theta_12)
sin_t23 = np.sin(theta_23)
sin_t13 = np.sin(theta_13)
cos_t12 = np.cos(theta_12)
cos_t23 = np.cos(theta_23)
cos_t13 = np.cos(theta_13)

sin2_t12 = sin_t12**2   # ~0.303
sin2_t23 = sin_t23**2   # ~0.573
sin2_t13 = sin_t13**2   # ~0.022

# PMNS matrix elements (magnitudes)
U_e1 = cos_t12 * cos_t13
U_e2 = sin_t12 * cos_t13
U_e3 = sin_t13
U_mu1 = abs(-sin_t12*cos_t23 - cos_t12*sin_t23*sin_t13*np.exp(1j*delta_CP))
U_mu2 = abs( cos_t12*cos_t23 - sin_t12*sin_t23*sin_t13*np.exp(1j*delta_CP))
U_mu3 = sin_t23 * cos_t13
U_tau1 = abs( sin_t12*sin_t23 - cos_t12*cos_t23*sin_t13*np.exp(1j*delta_CP))
U_tau2 = abs(-cos_t12*sin_t23 - sin_t12*cos_t23*sin_t13*np.exp(1j*delta_CP))
U_tau3 = cos_t23 * cos_t13

# Neutrino mass-squared differences (NuFIT 5.2, NO)
Dm2_21 = 7.42e-5   # eV^2 (solar)
Dm2_31 = 2.515e-3  # eV^2 (atmospheric, NO)
Dm2_32 = Dm2_31 - Dm2_21

# Dimensionless mass-squared ratio
r_mass = Dm2_21 / Dm2_31   # ~0.0295

# Leptonic Jarlskog invariant
J_PMNS = cos_t12 * sin_t12 * cos_t23 * sin_t23 * cos_t13**2 * sin_t13 * np.sin(delta_CP)
J_PMNS_abs = abs(J_PMNS)

print(f"\n  Mixing angles (NuFIT 5.2, Normal Ordering):")
print(f"    theta_12 = {theta_12_deg}°  (solar)")
print(f"    theta_23 = {theta_23_deg}°   (atmospheric)")
print(f"    theta_13 = {theta_13_deg}°   (reactor)")
print(f"    delta_CP = {delta_CP_deg}°   (Dirac CP phase)")

print(f"\n  sin^2 theta values:")
print(f"    sin^2 theta_12 = {sin2_t12:.5f}")
print(f"    sin^2 theta_23 = {sin2_t23:.5f}")
print(f"    sin^2 theta_13 = {sin2_t13:.5f}")

print(f"\n  PMNS magnitudes:")
print(f"    |U_e1|  = {U_e1:.5f}   |U_e2|  = {U_e2:.5f}   |U_e3|  = {U_e3:.5f}")
print(f"    |U_mu1| = {U_mu1:.5f}   |U_mu2| = {U_mu2:.5f}   |U_mu3| = {U_mu3:.5f}")
print(f"    |U_tau1|= {U_tau1:.5f}   |U_tau2|= {U_tau2:.5f}   |U_tau3|= {U_tau3:.5f}")

print(f"\n  Mass-squared differences:")
print(f"    Dm^2_21 = {Dm2_21:.2e} eV^2  (solar)")
print(f"    Dm^2_31 = {Dm2_31:.4e} eV^2  (atmospheric)")
print(f"    r = Dm^2_21/Dm^2_31 = {r_mass:.5f}")

print(f"\n  Leptonic Jarlskog invariant:")
print(f"    |J_PMNS| = {J_PMNS_abs:.5f}")
print(f"    J_PMNS   = {J_PMNS:.5f}  (negative: sin(197°) < 0)")

# ============================================================
# §2. TARGET QUANTITIES FOR SEARCH
# ============================================================

print("\n" + "=" * 72)
print("§2. TARGET QUANTITIES")
print("=" * 72)

targets = {
    # --- Mixing angles (sin^2, sin, and angle/pi) ---
    "sin2_t12":         sin2_t12,                     # 0.3034 (solar)
    "sin2_t23":         sin2_t23,                     # 0.5730 (atmospheric)
    "sin2_t13":         sin2_t13,                     # 0.0220 (reactor)
    "sin_t12":          sin_t12,                      # 0.5508
    "sin_t23":          sin_t23,                      # 0.7570
    "sin_t13":          sin_t13,                      # 0.1485
    "theta12/pi":       theta_12/np.pi,               # 0.1856
    "theta23/pi":       theta_23/np.pi,               # 0.2728
    "theta13/pi":       theta_13/np.pi,               # 0.04744

    # --- PMNS matrix elements (magnitudes) ---
    "|U_e1|":           U_e1,
    "|U_e2|":           U_e2,
    "|U_e3|":           U_e3,
    "|U_mu1|":          U_mu1,
    "|U_mu2|":          U_mu2,
    "|U_mu3|":          U_mu3,
    "|U_tau1|":         U_tau1,
    "|U_tau2|":         U_tau2,
    "|U_tau3|":         U_tau3,

    # --- tan^2 theta (useful ratios) ---
    "tan2_t12":         sin2_t12/(1-sin2_t12),        # ~0.436
    "tan2_t23":         sin2_t23/(1-sin2_t23),        # ~1.342

    # --- CP violation ---
    "|J_PMNS|":         J_PMNS_abs,                   # ~0.030
    "delta/pi":         delta_CP_deg/180,              # ~1.094
    "|sin_delta|":      abs(np.sin(delta_CP)),         # ~0.292

    # --- Mass-squared ratio (dimensionless) ---
    "r_Dm2":            r_mass,                        # ~0.0295

    # --- Quark-lepton complementarity ---
    #     theta_12_PMNS + theta_C_CKM ~ pi/4 (Raidal, 2004)
    "t12+tC":           theta_12 + np.radians(13.04),  # ~0.811 rad

    # --- Notable combinations ---
    "sin2_12*sin2_23":  sin2_t12 * sin2_t23,           # product
    "sin2_13/sin2_12":  sin2_t13 / sin2_t12,           # hierarchy
    "sin2_13/sin2_23":  sin2_t13 / sin2_t23,           # hierarchy
}

print(f"\n  {'Target':22s} {'Value':>14s} {'log10':>10s}")
print(f"  {'-'*50}")
for name, val in sorted(targets.items(), key=lambda x: abs(x[1])):
    if val > 0:
        print(f"  {name:22s} {val:14.6f} {np.log10(val):10.4f}")
    else:
        print(f"  {name:22s} {val:14.6f}  (negative)")

# ============================================================
# §3. EXHAUSTIVE SEARCH — CORE BLOCKS [-4, +4]
# ============================================================

print("\n" + "=" * 72)
print("§3. EXHAUSTIVE SEARCH — CORE BLOCKS [-4, +4]")
print("=" * 72)

def exhaustive_search(blocks_dict, exp_range, targets_dict, tolerance=0.05, max_display=5):
    """Search all formulas prod(b_i^n_i) matching each target."""
    names = list(blocks_dict.keys())
    values = np.array(list(blocks_dict.values()))
    n_blocks = len(values)
    log_vals = np.log(values)

    exp_list = list(range(exp_range[0], exp_range[1]+1))
    all_exps = np.array(list(cart_product(exp_list, repeat=n_blocks)))
    n_formulas = len(all_exps)

    all_log_results = all_exps @ log_vals

    results = {}
    for target_name, target_val in targets_dict.items():
        if target_val <= 0:
            results[target_name] = []
            continue
        log_target = np.log(target_val)

        diffs = np.abs(all_log_results - log_target)
        mask = diffs < np.log(1 + tolerance)
        indices = np.where(mask)[0]

        matches = []
        for idx in indices:
            exps = all_exps[idx]
            val = np.exp(all_log_results[idx])
            err = abs(val - target_val) / target_val * 100
            complexity = int(np.sum(np.abs(exps)))

            parts = []
            for bname, exp in zip(names, exps):
                if exp == 0:
                    continue
                elif exp == 1:
                    parts.append(bname)
                elif exp == -1:
                    parts.append(f"1/{bname}")
                elif exp > 0:
                    parts.append(f"{bname}^{exp}")
                else:
                    parts.append(f"{bname}^({exp})")
            formula = " * ".join(parts) if parts else "1"
            matches.append((formula, val, err, complexity, list(exps)))

        matches.sort(key=lambda x: (x[2], x[3]))
        results[target_name] = matches[:max_display]

    return results, n_formulas

print(f"\n  Blocks: {{R2, R3, REE, hv3, hv2}}")
print(f"  Exponents: [-4, +4], tolerance 5%")
print(f"  Search space: 9^5 = {9**5:,} formulas\n")

t0 = time.time()
results_core, n_core = exhaustive_search(blocks, (-4, 4), targets, tolerance=0.05)
dt = time.time() - t0
print(f"  Searched {n_core:,} formulas in {dt:.1f}s\n")

n_found = 0
n_sub1 = 0
summary = []

for target_name in targets:
    target_val = targets[target_name]
    matches = results_core.get(target_name, [])

    print(f"  --- {target_name} = {target_val:.6f} ---")
    if not matches:
        print(f"    (no matches within 5%)")
    else:
        for formula, val, err, complexity, exps in matches:
            star = " ***" if err < 0.5 else " **" if err < 1 else " *" if err < 2 else ""
            print(f"    {formula:50s} = {val:.6f}  err {err:6.3f}%  (C={complexity}){star}")
        best_err = matches[0][2]
        if best_err < 5:
            n_found += 1
        if best_err < 1:
            n_sub1 += 1
        summary.append((target_name, target_val, matches[0]))
    print()

# ============================================================
# §4. EXTENDED SEARCH — WIDER EXPONENTS [-6, +6]
# ============================================================

print("=" * 72)
print("§4. EXTENDED SEARCH — WIDER EXPONENTS [-6, +6]")
print("=" * 72)

small_targets = {}
for name, val in targets.items():
    if val > 0 and val < 0.05:
        small_targets[name] = val

if small_targets:
    print(f"\n  Searching {len(small_targets)} small targets with exponents [-6, +6]")
    print(f"  Search space: 13^5 = {13**5:,} formulas\n")

    t0 = time.time()
    results_wide, n_wide = exhaustive_search(blocks, (-6, 6), small_targets, tolerance=0.05)
    dt = time.time() - t0
    print(f"  Searched {n_wide:,} formulas in {dt:.1f}s\n")

    for target_name in small_targets:
        target_val = small_targets[target_name]
        matches = results_wide.get(target_name, [])

        print(f"  --- {target_name} = {target_val:.6e} ---")
        if not matches:
            print(f"    (no matches within 5%)")
        else:
            for formula, val, err, complexity, exps in matches:
                star = " ***" if err < 0.5 else " **" if err < 1 else " *" if err < 2 else ""
                print(f"    {formula:55s} = {val:.6e}  err {err:6.3f}%  (C={complexity}){star}")
            # Update summary if better
            for i, (sn, sv, sm) in enumerate(summary):
                if sn == target_name and matches[0][2] < sm[2]:
                    summary[i] = (sn, sv, matches[0])
        print()

# ============================================================
# §5. SEARCH INCLUDING pi FACTORS
# ============================================================

print("=" * 72)
print("§5. SEARCH INCLUDING pi FACTORS")
print("=" * 72)

print("\n  Testing formulas of the form: prefactor * prod(blocks^n)\n")

prefactors = {
    "1":         1.0,
    "1/(4pi)":   1/(4*np.pi),
    "1/(4pi)^2": 1/(4*np.pi)**2,
    "1/(2pi)":   1/(2*np.pi),
    "1/pi":      1/np.pi,
    "pi":        np.pi,
    "pi/4":      np.pi/4,
    "sqrt(pi)":  np.sqrt(np.pi),
}

pmns_key = {
    "sin2_t12":   sin2_t12,
    "sin2_t23":   sin2_t23,
    "sin2_t13":   sin2_t13,
    "sin_t12":    sin_t12,
    "sin_t23":    sin_t23,
    "sin_t13":    sin_t13,
    "|J_PMNS|":   J_PMNS_abs,
    "r_Dm2":      r_mass,
}

for pf_name, pf_val in prefactors.items():
    modified_targets = {k: v/pf_val for k, v in pmns_key.items() if v/pf_val > 0}
    results_pf, _ = exhaustive_search(blocks, (-4, 4), modified_targets, tolerance=0.02, max_display=3)

    any_found = False
    for target_name, matches in results_pf.items():
        if matches and matches[0][2] < 1.0:
            if not any_found:
                print(f"  Prefactor: {pf_name}")
                any_found = True
            orig_val = pmns_key[target_name]
            formula = matches[0][0]
            err = matches[0][2]
            star = " ***" if err < 0.5 else " **" if err < 1 else ""
            print(f"    {target_name:16s} = {pf_name} * {formula:40s}  err {err:6.3f}%{star}")

    if any_found:
        print()

# ============================================================
# §6. SEARCH WITH RATIONAL PREFACTORS
# ============================================================

print("=" * 72)
print("§6. SEARCH WITH RATIONAL PREFACTORS")
print("=" * 72)

print("""
  Testing: (p/q) * prod(blocks^n) for small p, q.
  Particularly interesting for PMNS since tribimaximal mixing uses
  simple fractions: sin^2(t12) = 1/3, sin^2(t23) = 1/2, sin^2(t13) = 0.
""")

rational_prefactors = set()
for p in range(1, 28):
    for q in range(1, 65):
        rational_prefactors.add((p, q))

pmns_primary = {
    "sin2_t12":  sin2_t12,
    "sin2_t23":  sin2_t23,
    "sin2_t13":  sin2_t13,
    "sin_t12":   sin_t12,
    "sin_t23":   sin_t23,
    "sin_t13":   sin_t13,
    "|J_PMNS|":  J_PMNS_abs,
    "r_Dm2":     r_mass,
}

best_per_target = {k: (None, None, 100.0) for k in pmns_primary}

block_names = list(blocks.keys())
block_vals = np.array(list(blocks.values()))
log_block_vals = np.log(block_vals)

exp_list = list(range(-4, 5))
all_exps = np.array(list(cart_product(exp_list, repeat=5)))
all_log_products = all_exps @ log_block_vals

print(f"  Testing {len(rational_prefactors)} rational prefactors x {len(all_exps)} formulas")
print(f"  = {len(rational_prefactors) * len(all_exps):,} total evaluations\n")

t0 = time.time()
for p, q in rational_prefactors:
    pq = p / q
    log_pq = np.log(pq)
    all_log_with_pq = all_log_products + log_pq

    for target_name, target_val in pmns_primary.items():
        if target_val <= 0:
            continue
        log_target = np.log(target_val)
        diffs = np.abs(all_log_with_pq - log_target)
        best_idx = np.argmin(diffs)
        best_err_log = diffs[best_idx]

        if best_err_log < np.log(1.005):  # < 0.5% error
            val = pq * np.exp(all_log_products[best_idx])
            err = abs(val - target_val) / target_val * 100

            if err < best_per_target[target_name][2]:
                exps = all_exps[best_idx]
                parts = []
                for bname, exp in zip(block_names, exps):
                    if exp == 0:
                        continue
                    elif exp == 1:
                        parts.append(bname)
                    elif exp == -1:
                        parts.append(f"1/{bname}")
                    elif exp > 0:
                        parts.append(f"{bname}^{exp}")
                    else:
                        parts.append(f"{bname}^({exp})")
                block_str = " * ".join(parts) if parts else "1"
                formula = f"({p}/{q}) * {block_str}" if q != 1 else f"{p} * {block_str}"
                best_per_target[target_name] = (formula, val, err)

dt = time.time() - t0
print(f"  Completed in {dt:.1f}s\n")

print("  Best results with rational prefactors:")
print(f"  {'Target':16s} {'Formula':55s} {'Value':>12s} {'Error':>8s}")
print(f"  {'-'*95}")
for target_name, (formula, val, err) in best_per_target.items():
    if formula is not None:
        star = " ***" if err < 0.5 else " **" if err < 1 else ""
        print(f"  {target_name:16s} {formula:55s} {val:12.6f} {err:7.3f}%{star}")
    else:
        print(f"  {target_name:16s} {'(no match < 0.5%)':55s}")

# ============================================================
# §7. TRIBIMAXIMAL / SPECIAL PATTERN COMPARISON
# ============================================================

print("\n" + "=" * 72)
print("§7. TRIBIMAXIMAL / SPECIAL PATTERN COMPARISON")
print("=" * 72)

print("""
  Tribimaximal mixing (Harrison-Perkins-Scott, 2002):
    sin^2 t12 = 1/3,  sin^2 t23 = 1/2,  sin^2 t13 = 0

  Actual deviations from tribimaximal could be graph-computable.
""")

# Deviation parameters
d12 = sin2_t12 - 1/3
d23 = sin2_t23 - 1/2
d13 = sin2_t13 - 0

print(f"  Deviations from tribimaximal:")
print(f"    sin^2 t12 - 1/3 = {d12:+.5f}  ({d12/(1/3)*100:+.1f}% of 1/3)")
print(f"    sin^2 t23 - 1/2 = {d23:+.5f}  ({d23/(1/2)*100:+.1f}% of 1/2)")
print(f"    sin^2 t13 - 0   = {d13:+.5f}  (= reactor angle)")

# Search for deviations
dev_targets = {}
if abs(d12) > 0:
    dev_targets["d_12 = s2_12 - 1/3"] = abs(d12)
if abs(d23) > 0:
    dev_targets["d_23 = s2_23 - 1/2"] = abs(d23)

if dev_targets:
    print(f"\n  Searching graph formulas for TB deviations...")
    results_tb, _ = exhaustive_search(blocks, (-4, 4), dev_targets, tolerance=0.05, max_display=5)

    for target_name, matches in results_tb.items():
        print(f"\n  --- {target_name} = {dev_targets[target_name]:.6f} ---")
        if matches:
            for formula, val, err, complexity, exps in matches:
                star = " ***" if err < 0.5 else " **" if err < 1 else " *" if err < 2 else ""
                print(f"    {formula:50s} = {val:.6f}  err {err:6.3f}%  (C={complexity}){star}")
        else:
            print(f"    (no matches within 5%)")

# Check: is sin^2(t12) = 1/3 * (some block formula)?
print(f"\n  Testing: sin^2 t12 = (1/3) * F(blocks)")
test_target = {"F = 3*sin2_t12": 3 * sin2_t12}
results_f, _ = exhaustive_search(blocks, (-4, 4), test_target, tolerance=0.02, max_display=5)
for matches in results_f.values():
    if matches:
        for formula, val, err, complexity, exps in matches:
            star = " ***" if err < 0.5 else " **" if err < 1 else ""
            print(f"    sin^2 t12 = (1/3) * {formula:40s}  err {err:6.3f}%{star}")

print(f"\n  Testing: sin^2 t23 = (1/2) * F(blocks)")
test_target2 = {"F = 2*sin2_t23": 2 * sin2_t23}
results_f2, _ = exhaustive_search(blocks, (-4, 4), test_target2, tolerance=0.02, max_display=5)
for matches in results_f2.values():
    if matches:
        for formula, val, err, complexity, exps in matches:
            star = " ***" if err < 0.5 else " **" if err < 1 else ""
            print(f"    sin^2 t23 = (1/2) * {formula:40s}  err {err:6.3f}%{star}")

# ============================================================
# §8. SUMMARY
# ============================================================

print("\n" + "=" * 72)
print("§8. SUMMARY")
print("=" * 72)

print("\n  Best graph formulas for PMNS parameters:")
print(f"  {'Target':22s} {'Experimental':>12s} {'Best formula':>50s} {'Err':>8s}")
print(f"  {'-'*97}")

all_best = {}
for target_name, target_val in targets.items():
    matches = results_core.get(target_name, [])
    if matches:
        all_best[target_name] = (matches[0][0], matches[0][1], matches[0][2], matches[0][4])

for target_name in sorted(all_best.keys(), key=lambda x: all_best[x][2]):
    formula, val, err, exps = all_best[target_name]
    target_val = targets[target_name]
    star = " ***" if err < 0.5 else " **" if err < 1 else " *" if err < 2 else ""
    print(f"  {target_name:22s} {target_val:12.6f} {formula:>50s} {err:7.3f}%{star}")

n_sub05 = sum(1 for v in all_best.values() if v[2] < 0.5)
n_sub1 = sum(1 for v in all_best.values() if v[2] < 1)
n_sub2 = sum(1 for v in all_best.values() if v[2] < 2)
n_sub5 = sum(1 for v in all_best.values() if v[2] < 5)

print(f"\n  Found with < 0.5% error: {n_sub05}")
print(f"  Found with < 1% error:   {n_sub1}")
print(f"  Found with < 2% error:   {n_sub2}")
print(f"  Found with < 5% error:   {n_sub5}")
print(f"  Total targets: {len(targets)}")

# Compare with CKM statistics
print(f"\n  CKM comparison: 23/24 sub-1% (quark sector)")
print(f"  PMNS:           {n_sub1}/{len(targets)} sub-1% (lepton sector)")

# ============================================================
# §9. INDEPENDENCE CHECK
# ============================================================

print("\n" + "=" * 72)
print("§9. INDEPENDENCE CHECK")
print("=" * 72)

print("\n  Exponent vectors of sub-1% formulas:")
print(f"  {'Target':22s} {'R2':>4s} {'R3':>4s} {'REE':>4s} {'hv3':>4s} {'hv2':>4s}")
print(f"  {'-'*50}")

sub1_exps = []
sub1_names = []
for target_name in sorted(all_best.keys(), key=lambda x: all_best[x][2]):
    formula, val, err, exps = all_best[target_name]
    if err < 1.0:
        sub1_exps.append(exps)
        sub1_names.append(target_name)
        print(f"  {target_name:22s} {exps[0]:4d} {exps[1]:4d} {exps[2]:4d} {exps[3]:4d} {exps[4]:4d}")

if sub1_exps:
    M = np.array(sub1_exps, dtype=float)
    rank = np.linalg.matrix_rank(M)
    print(f"\n  Exponent matrix rank: {rank} (from {len(sub1_exps)} formulas)")
    print(f"  New independent directions: {rank}")

    # SVD for insight
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    print(f"\n  Singular values: {', '.join(f'{s:.3f}' for s in S)}")

    print(f"\n  If PMNS exponent vectors span directions NOT covered by")
    print(f"  CKM + mass ratio formulas, they add independent sigma.")
else:
    print("\n  No sub-1% formulas found in core search.")
    print("  Check rational prefactor results (§6) for sub-1% matches.")

print("\n" + "=" * 72)
print("PMNS SEARCH COMPLETE")
print("=" * 72)
