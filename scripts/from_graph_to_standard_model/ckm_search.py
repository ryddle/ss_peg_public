#!/usr/bin/env python3
"""
ckm_search.py — Systematic search for CKM matrix parameters from graph invariants
===================================================================================

The CKM (Cabibbo-Kobayashi-Maskawa) matrix describes quark flavor mixing.
It has 4 independent parameters in the standard parametrization:
    - 3 mixing angles: theta_12, theta_23, theta_13
    - 1 CP-violating phase: delta_CP

We search for graph formulas matching these parameters and derived quantities
(individual |V_ij| matrix elements, Wolfenstein parameters, Jarlskog invariant).

Building blocks: the same 5 spectral invariants that produced 21/21 predictions.

Author: SS-PEG program
Date: April 2026
"""

import numpy as np
from itertools import product as cart_product
import time

np.random.seed(42)

# ============================================================
# $0. GRAPH BUILDING BLOCKS
# ============================================================

print("=" * 72)
print("CKM MATRIX PARAMETER SEARCH FROM GRAPH INVARIANTS")
print("=" * 72)

# Core spectral invariants
R2  = 1 / 2
R3  = (np.sqrt(57) - 3) / (2 * np.sqrt(17))  # 0.55175
R_EE = 1 / np.sqrt(2)                          # 0.70711
hv3 = 3
hv2 = 2

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
# $1. CKM EXPERIMENTAL VALUES (PDG 2024)
# ============================================================

print("\n" + "=" * 72)
print("$1. CKM EXPERIMENTAL VALUES (PDG 2024)")
print("=" * 72)

# Standard parametrization angles
sin_t12 = 0.22500   # sin(theta_12) = lambda (Cabibbo)
sin_t23 = 0.04182   # sin(theta_23)
sin_t13 = 0.003660  # sin(theta_13)
delta_CP = 1.144     # CP phase in radians (~65.5 deg)

cos_t12 = np.sqrt(1 - sin_t12**2)
cos_t23 = np.sqrt(1 - sin_t23**2)
cos_t13 = np.sqrt(1 - sin_t13**2)

# CKM matrix elements (magnitudes)
V_ud = cos_t12 * cos_t13                                    # 0.97435
V_us = sin_t12 * cos_t13                                    # 0.22500
V_ub = sin_t13                                               # 0.00366
V_cd_mag = abs(-sin_t12*cos_t23 - cos_t12*sin_t23*sin_t13*np.exp(1j*delta_CP))
V_cs_mag = abs(cos_t12*cos_t23 - sin_t12*sin_t23*sin_t13*np.exp(1j*delta_CP))
V_cb = sin_t23 * cos_t13                                    # 0.04182
V_td_mag = abs(sin_t12*sin_t23 - cos_t12*cos_t23*sin_t13*np.exp(1j*delta_CP))
V_ts_mag = abs(-cos_t12*sin_t23 - sin_t12*cos_t23*sin_t13*np.exp(1j*delta_CP))
V_tb = cos_t23 * cos_t13                                    # 0.99913

# Wolfenstein parameters
lam_W  = sin_t12                      # lambda ~ 0.225
A_W    = sin_t23 / lam_W**2           # A ~ 0.826
rho_bar = 0.159                        # rho-bar (from global fit)
eta_bar = 0.348                        # eta-bar (from global fit)

# Jarlskog invariant (CP violation measure)
J_CP = cos_t12 * cos_t23 * cos_t13**2 * sin_t12 * sin_t23 * sin_t13 * np.sin(delta_CP)
# J ~ 3.08e-5

# Cabibbo angle
theta_C = np.arcsin(sin_t12)  # ~13 degrees
sin_tC = sin_t12
cos_tC = cos_t12

print("\n  Standard parametrization:")
print(f"    sin theta_12 (Cabibbo) = {sin_t12:.5f}")
print(f"    sin theta_23           = {sin_t23:.5f}")
print(f"    sin theta_13           = {sin_t13:.6f}")
print(f"    delta_CP               = {delta_CP:.4f} rad = {np.degrees(delta_CP):.1f} deg")
print(f"\n  CKM magnitudes:")
print(f"    |V_ud| = {V_ud:.5f}   |V_us| = {V_us:.5f}   |V_ub| = {V_ub:.6f}")
print(f"    |V_cd| = {V_cd_mag:.5f}   |V_cs| = {V_cs_mag:.5f}   |V_cb| = {V_cb:.5f}")
print(f"    |V_td| = {V_td_mag:.6f}  |V_ts| = {V_ts_mag:.5f}   |V_tb| = {V_tb:.5f}")
print(f"\n  Wolfenstein:")
print(f"    lambda = {lam_W:.5f},  A = {A_W:.4f},  rho_bar = {rho_bar:.4f},  eta_bar = {eta_bar:.4f}")
print(f"\n  Jarlskog invariant:")
print(f"    J_CP = {J_CP:.4e}")

# ============================================================
# $2. TARGET QUANTITIES FOR SEARCH
# ============================================================

print("\n" + "=" * 72)
print("$2. TARGET QUANTITIES")
print("=" * 72)

# We search for multiple representations of the CKM parameters.
# Some may be more natural than others from the graph perspective.

targets = {
    # --- Mixing angles (sin, sin^2, and angle/pi) ---
    "sin_t12":        sin_t12,          # 0.22500 (Cabibbo)
    "sin^2_t12":      sin_t12**2,       # 0.05063
    "sin_t23":        sin_t23,          # 0.04182
    "sin^2_t23":      sin_t23**2,       # 0.001749
    "sin_t13":        sin_t13,          # 0.003660
    "sin^2_t13":      sin_t13**2,       # 1.340e-5
    
    # --- CKM matrix elements ---
    "|V_us|":         V_us,             # 0.22500 (= sin_t12)
    "|V_cb|":         V_cb,             # 0.04182
    "|V_ub|":         V_ub,             # 0.003660
    "|V_td|":         V_td_mag,         # ~0.00854
    "|V_ts|":         V_ts_mag,         # ~0.04110
    
    # --- Wolfenstein parameters ---
    "lambda_W":       lam_W,            # 0.22500
    "A_Wolf":         A_W,              # 0.826
    
    # --- Ratios of CKM elements (scale-free) ---
    "|V_ub/V_cb|":    V_ub/V_cb,        # 0.0875
    "|V_td/V_ts|":    V_td_mag/V_ts_mag,# ~0.208
    "|V_us/V_ud|":    V_us/V_ud,        # 0.2310 (= tan theta_12)
    "|V_cb/V_us|":    V_cb/V_us,        # 0.1859
    "|V_ub/V_us|":    V_ub/V_us,        # 0.01627
    
    # --- CP violation ---
    "J_CP":           J_CP,             # 3.08e-5
    "delta/pi":       delta_CP/np.pi,   # 0.3641
    "sin_delta":      np.sin(delta_CP), # 0.9096
    
    # --- Hierarchical structure ---
    "lam^2":          lam_W**2,         # 0.05063  (controls V_cb ~ A*lam^2)
    "lam^3":          lam_W**3,         # 0.01139  (controls V_ub ~ lam^3)
    "lam^4":          lam_W**4,         # 0.002563 (controls V_td)
}

print(f"\n  {'Target':16s} {'Value':>14s} {'log10':>10s}")
print(f"  {'-'*44}")
for name, val in sorted(targets.items(), key=lambda x: abs(x[1])):
    if val > 0:
        print(f"  {name:16s} {val:14.6f} {np.log10(val):10.4f}")
    else:
        print(f"  {name:16s} {val:14.6f}  (negative)")

# ============================================================
# $3. EXHAUSTIVE SEARCH — CORE BLOCKS [-4, +4]
# ============================================================

print("\n" + "=" * 72)
print("$3. EXHAUSTIVE SEARCH — CORE BLOCKS")
print("=" * 72)

def exhaustive_search(blocks_dict, exp_range, targets_dict, tolerance=0.05, max_display=5):
    """
    Search all formulas prod(b_i^n_i) matching each target.
    """
    names = list(blocks_dict.keys())
    values = np.array(list(blocks_dict.values()))
    n_blocks = len(values)
    log_vals = np.log(values)
    
    exp_list = list(range(exp_range[0], exp_range[1]+1))
    all_exps = np.array(list(cart_product(exp_list, repeat=n_blocks)))
    n_formulas = len(all_exps)
    
    # Compute all formula values in log space
    all_log_results = all_exps @ log_vals
    
    results = {}
    for target_name, target_val in targets_dict.items():
        if target_val <= 0:
            results[target_name] = []
            continue
        log_target = np.log(target_val)
        
        # Find matches within tolerance
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

# Run search
print(f"\n  Blocks: {{R2, R3, REE, hv3, hv2}}")
print(f"  Exponents: [-4, +4], tolerance 5%")
print(f"  Search space: 9^5 = {9**5:,} formulas\n")

t0 = time.time()
results_core, n_core = exhaustive_search(blocks, (-4, 4), targets, tolerance=0.05)
dt = time.time() - t0
print(f"  Searched {n_core:,} formulas in {dt:.1f}s\n")

# Display results
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
            print(f"    {formula:45s} = {val:.6f}  err {err:6.3f}%  (C={complexity}){star}")
        best_err = matches[0][2]
        if best_err < 5:
            n_found += 1
        if best_err < 1:
            n_sub1 += 1
        summary.append((target_name, target_val, matches[0]))
    print()

# ============================================================
# $4. EXTENDED SEARCH — WIDER EXPONENTS [-6, +6]
# ============================================================

print("=" * 72)
print("$4. EXTENDED SEARCH — WIDER EXPONENTS [-6, +6]")
print("=" * 72)

# For the very small quantities (sin_t13, J_CP), wider exponents may help
small_targets = {}
for name, val in targets.items():
    if val > 0 and val < 0.01:
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
                print(f"    {formula:50s} = {val:.6e}  err {err:6.3f}%  (C={complexity}){star}")
            # Update summary if better
            for i, (sn, sv, sm) in enumerate(summary):
                if sn == target_name and matches[0][2] < sm[2]:
                    summary[i] = (sn, sv, matches[0])
        print()

# ============================================================
# $5. SEARCH WITH 4*pi PREFACTOR
# ============================================================

print("=" * 72)
print("$5. SEARCH INCLUDING 4*pi FACTORS")
print("=" * 72)

# The master formula uses 1/(4*pi) as a universal prefactor.
# CKM elements might use similar factors.

print("\n  Testing formulas of the form: f(blocks) / (4*pi)^k for k=0,1,2\n")

prefactors = {
    "1":        1.0,
    "1/(4pi)":  1/(4*np.pi),
    "1/(4pi)^2": 1/(4*np.pi)**2,
    "1/(2pi)":  1/(2*np.pi),
    "1/pi":     1/np.pi,
    "pi":       np.pi,
}

# Key CKM targets
ckm_key = {
    "sin_t12":   sin_t12,
    "sin_t23":   sin_t23,
    "sin_t13":   sin_t13,
    "|V_ub/V_cb|": V_ub/V_cb,
    "|V_td/V_ts|": V_td_mag/V_ts_mag,
    "J_CP":      J_CP,
    "lambda_W":  lam_W,
    "A_Wolf":    A_W,
}

for pf_name, pf_val in prefactors.items():
    modified_targets = {k: v/pf_val for k, v in ckm_key.items() if v/pf_val > 0}
    results_pf, _ = exhaustive_search(blocks, (-4, 4), modified_targets, tolerance=0.02, max_display=3)
    
    any_found = False
    for target_name, matches in results_pf.items():
        if matches and matches[0][2] < 1.0:
            if not any_found:
                print(f"  Prefactor: {pf_name}")
                any_found = True
            orig_val = ckm_key[target_name]
            formula = matches[0][0]
            err = matches[0][2]
            star = " ***" if err < 0.5 else " **" if err < 1 else ""
            print(f"    {target_name:16s} = {pf_name} * {formula:35s}  err {err:6.3f}%{star}")
    
    if any_found:
        print()

# ============================================================
# $6. SEARCH OVER RATIONALS: a/b * prod(blocks^n)
# ============================================================

print("=" * 72)
print("$6. SEARCH WITH RATIONAL PREFACTORS")
print("=" * 72)

print("""
  The Weinberg angle formula sin^2(theta_W) = 27*R3/64 uses a rational
  prefactor 27/64. CKM angles may need similar rational coefficients.
  
  Testing: (p/q) * prod(blocks^n) for small p, q.
""")

# Generate rational prefactors p/q with p in [1,27], q in [1,64]
rational_prefactors = set()
for p in range(1, 28):
    for q in range(1, 65):
        rational_prefactors.add((p, q))

ckm_primary = {
    "sin_t12":   sin_t12,
    "sin_t23":   sin_t23,
    "sin_t13":   sin_t13,
    "J_CP":      J_CP,
}

best_per_target = {k: (None, None, 100.0) for k in ckm_primary}

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
    
    for target_name, target_val in ckm_primary.items():
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
# $7. SUMMARY AND SIGNIFICANCE
# ============================================================

print("\n" + "=" * 72)
print("$7. SUMMARY")
print("=" * 72)

print("\n  Best graph formulas for CKM parameters:")
print(f"  {'Target':16s} {'Experimental':>12s} {'Best formula':>45s} {'Err':>8s}")
print(f"  {'-'*85}")

# Collect all best results
all_best = {}
for target_name, target_val in targets.items():
    matches = results_core.get(target_name, [])
    if matches:
        all_best[target_name] = (matches[0][0], matches[0][1], matches[0][2], matches[0][4])

for target_name in sorted(all_best.keys(), key=lambda x: all_best[x][2]):
    formula, val, err, exps = all_best[target_name]
    target_val = targets[target_name]
    star = " ***" if err < 0.5 else " **" if err < 1 else " *" if err < 2 else ""
    print(f"  {target_name:16s} {target_val:12.6f} {formula:>45s} {err:7.3f}%{star}")

print(f"\n  Found with < 1% error: {sum(1 for v in all_best.values() if v[2] < 1)}")
print(f"  Found with < 2% error: {sum(1 for v in all_best.values() if v[2] < 2)}")
print(f"  Found with < 5% error: {sum(1 for v in all_best.values() if v[2] < 5)}")
print(f"  Total targets: {len(targets)}")

# ============================================================
# $8. INDEPENDENCE CHECK
# ============================================================

print("\n" + "=" * 72)
print("$8. INDEPENDENCE CHECK")
print("=" * 72)

# Check which formulas are algebraically independent of the mass ratio formulas
print("\n  Exponent vectors of sub-1% formulas:")
print(f"  {'Target':16s} {'R2':>4s} {'R3':>4s} {'REE':>4s} {'hv3':>4s} {'hv2':>4s}")
print(f"  {'-'*44}")

sub1_exps = []
sub1_names = []
for target_name in sorted(all_best.keys(), key=lambda x: all_best[x][2]):
    formula, val, err, exps = all_best[target_name]
    if err < 1.0:
        sub1_exps.append(exps)
        sub1_names.append(target_name)
        print(f"  {target_name:16s} {exps[0]:4d} {exps[1]:4d} {exps[2]:4d} {exps[3]:4d} {exps[4]:4d}")

if sub1_exps:
    M = np.array(sub1_exps, dtype=float)
    rank = np.linalg.matrix_rank(M)
    print(f"\n  Exponent matrix rank: {rank} (from {len(sub1_exps)} formulas)")
    print(f"  New independent directions: {rank} (potential new sigma contribution)")
else:
    print("\n  No sub-1% formulas found.")

print("\n" + "=" * 72)
print("SEARCH COMPLETE")
print("=" * 72)
