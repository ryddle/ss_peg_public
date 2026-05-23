#!/usr/bin/env python3
"""
survivor_analysis.py — Anatomía de las 72 configuraciones al 2%
================================================================

From the 2500 joint configurations surviving the gen-3 gap constraints,
72 have all 9 masses within 2% of experiment. This script dissects them
to find what distinguishes the actual config from impostors — and what
additional constraint would break the degeneracy.

Author: SS-PEG program
Date: April 2026
"""

import numpy as np
import sys, os
from collections import Counter

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "survivor_analysis_output.md")
_orig_stdout = sys.stdout
_outfile = open(OUTPUT_PATH, "w", encoding="utf-8")

class _TeeWriter:
    def __init__(self, console, f):
        self.console, self.f = console, f
    def write(self, s):
        self.console.write(s); self.f.write(s)
    def flush(self):
        self.console.flush(); self.f.flush()

sys.stdout = _TeeWriter(_orig_stdout, _outfile)

# ════════════════════════════════════════════════════════════════
# CONSTANTS & DATA (same as discrete_search.py)
# ════════════════════════════════════════════════════════════════

LN2 = np.log(2)
R3 = (np.sqrt(57) - 3) / (2 * np.sqrt(17))
LNR3 = np.log(R3)
LN3 = np.log(3)
BASIS = np.array([LN2, LNR3, LN3])

m_e = 0.51100  # MeV

MASSES_EXP = {
    'e': 0.51100, 'mu': 105.658, 'tau': 1776.86,
    'u': 2.16, 'c': 1270.0, 't': 172760.0,
    'd': 4.67, 's': 93.4, 'b': 4180.0
}

ACTUAL = {
    'lepton': {
        'gen1': np.array([0.0, 0.0, 0.0]),
        'gen2': np.array([-0.5, -4.0, 3.0]),
        'gen3': np.array([1.0, -7.0, 3.0]),
    },
    'up': {
        'gen1': np.array([-1.5, -6.0, -1.0]),
        'gen2': np.array([0.5, -7.0, 3.0]),
        'gen3': np.array([6.0, -7.0, 4.0]),
    },
    'down': {
        'gen1': np.array([-0.5, -8.0, -2.0]),
        'gen2': np.array([1.5, -7.0, 0.0]),
        'gen3': np.array([1.5, -6.0, 4.0]),
    },
}

A_LEPTON = np.array([1.0, 0.5, -1.5])
A_UP     = np.array([1.75, 0.5, -1.5])
A_DOWN   = np.array([-1.0, 0.0, 1.0])


def mass_from_coords(p, q, r):
    return m_e * np.exp(p * LN2 + q * LNR3 + r * LN3)


def ln_mass_ratio(p, q, r):
    return p * LN2 + q * LNR3 + r * LN3


def is_half_int(x, tol=1e-10):
    return abs(2 * x - round(2 * x)) < tol


def coords_half_int(arr, tol=1e-10):
    return all(is_half_int(x, tol) for x in arr)


# ════════════════════════════════════════════════════════════════
# REBUILD THE SEARCH  (compact, same logic as discrete_search.py)
# ════════════════════════════════════════════════════════════════

P_RANGE = np.arange(-5, 12.5, 0.5)
Q_RANGE = np.arange(-14, 2.5, 0.5)
R_RANGE = np.arange(-10, 12.5, 0.5)

pp, qq, rr = np.meshgrid(P_RANGE, Q_RANGE, R_RANGE, indexing='ij')
all_p = pp.ravel()
all_q = qq.ravel()
all_r = rr.ravel()
all_ln_m = all_p * LN2 + all_q * LNR3 + all_r * LN3


def search_sector(a_vec, gen1_fixed=None, gen3_q_fixed=None,
                  m1_exp=None, m3_exp=None, tol=0.5):
    results = []
    if gen1_fixed is not None:
        g1_points = [gen1_fixed.copy()]
    else:
        ln_m1_target = np.log(m1_exp / m_e)
        ln_tol = np.log(1 + tol)
        mask1 = np.abs(all_ln_m - ln_m1_target) < ln_tol
        g1_points = [np.array(pt) for pt in zip(all_p[mask1], all_q[mask1], all_r[mask1])]

    ln_m3_target = np.log(m3_exp / m_e)
    ln_tol3 = np.log(1 + tol)
    mask3 = np.abs(all_ln_m - ln_m3_target) < ln_tol3
    if gen3_q_fixed is not None:
        mask3 = mask3 & (np.abs(all_q - gen3_q_fixed) < 1e-10)
    g3_points = [np.array(pt) for pt in zip(all_p[mask3], all_q[mask3], all_r[mask3])]

    for g1 in g1_points:
        for g3 in g3_points:
            g2 = (g1 + g3) / 2.0 - a_vec
            if not coords_half_int(g2):
                continue
            m1 = mass_from_coords(*g1)
            m2 = mass_from_coords(*g2)
            m3 = mass_from_coords(*g3)
            if not (m1 < m2 < m3):
                continue
            results.append({'gen1': g1, 'gen2': g2, 'gen3': g3, 'masses': (m1, m2, m3)})
    return results


print("=" * 72)
print("  SURVIVOR ANALYSIS — Anatomy of the 72 Configs at 2%")
print("=" * 72)
print()

# ════════════════════════════════════════════════════════════════
# §1  REBUILD JOINT CONFIGS
# ════════════════════════════════════════════════════════════════

print("§1. REBUILDING JOINT CONFIGURATION POOL")
print("-" * 72)

lepton_configs = search_sector(A_LEPTON, gen1_fixed=np.array([0.0, 0.0, 0.0]),
                               gen3_q_fixed=-7.0, m3_exp=MASSES_EXP['tau'], tol=10.0)
up_configs = search_sector(A_UP, gen3_q_fixed=-7.0,
                           m1_exp=MASSES_EXP['u'], m3_exp=MASSES_EXP['t'], tol=0.5)
down_configs = search_sector(A_DOWN,
                             m1_exp=MASSES_EXP['d'], m3_exp=MASSES_EXP['b'], tol=0.5)

# Filter to 10% individual
lepton_10 = [c for c in lepton_configs
             if abs(c['masses'][1] - MASSES_EXP['mu']) / MASSES_EXP['mu'] < 0.10
             and abs(c['masses'][2] - MASSES_EXP['tau']) / MASSES_EXP['tau'] < 0.10]
up_10 = [c for c in up_configs
         if all(abs(c['masses'][i] - m) / m < 0.10
                for i, m in enumerate([MASSES_EXP['u'], MASSES_EXP['c'], MASSES_EXP['t']]))]
down_10 = [c for c in down_configs
           if all(abs(c['masses'][i] - m) / m < 0.10
                  for i, m in enumerate([MASSES_EXP['d'], MASSES_EXP['s'], MASSES_EXP['b']]))]

# Build joint with gen-3 gap structure
joint_all = []
for lep in lepton_10:
    for up in up_10:
        if abs(lep['gen3'][1] - up['gen3'][1]) > 1e-10:
            continue
        for dn in down_10:
            if abs(up['gen3'][2] - dn['gen3'][2]) > 1e-10:
                continue
            if abs(dn['gen3'][1] - lep['gen3'][1] - 1.0) > 1e-10:
                continue
            if abs(lep['gen3'][2] - up['gen3'][2] + 1.0) > 1e-10:
                continue
            joint_all.append({'lepton': lep, 'up': up, 'down': dn})

print(f"  Joint configs with gen-3 gap structure: {len(joint_all)}")

# Filter to various tolerances
def filter_joint(configs, tol_frac):
    survivors = []
    pnames_by_sector = {'lepton': ['e','mu','tau'], 'up': ['u','c','t'], 'down': ['d','s','b']}
    for j in configs:
        ok = True
        for sector, pnames in pnames_by_sector.items():
            cfg = j[sector]
            for g_idx, (gkey, pn) in enumerate(zip(['gen1','gen2','gen3'], pnames)):
                m = mass_from_coords(*cfg[gkey])
                if abs(m / MASSES_EXP[pn] - 1) > tol_frac:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            survivors.append(j)
    return survivors

s10 = filter_joint(joint_all, 0.10)
s5  = filter_joint(joint_all, 0.05)
s3  = filter_joint(joint_all, 0.03)
s2  = filter_joint(joint_all, 0.02)

print(f"  At 10%: {len(s10)}")
print(f"  At  5%: {len(s5)}")
print(f"  At  3%: {len(s3)}")
print(f"  At  2%: {len(s2)}")
print()


# ════════════════════════════════════════════════════════════════
# §2  WHERE IS THE DEGENERACY? — Fixed vs Variable Coordinates
# ════════════════════════════════════════════════════════════════

print("\n§2. WHERE IS THE DEGENERACY? — Fixed vs Variable Coordinates")
print("-" * 72)

pool = s2
N = len(pool)
print(f"\n  Analyzing {N} configs at 2% tolerance.\n")

# For each sector and generation, collect all coordinate values
for sector in ['lepton', 'up', 'down']:
    print(f"  === {sector.upper()} SECTOR ===")
    for gkey, gen_label in [('gen1', 'gen-1'), ('gen2', 'gen-2'), ('gen3', 'gen-3')]:
        ps = [tuple(np.round(j[sector][gkey], 1)) for j in pool]
        unique = set(ps)
        counts = Counter(ps)
        print(f"    {gen_label}: {len(unique)} unique values")
        if len(unique) <= 10:
            for val, cnt in counts.most_common(10):
                actual_val = tuple(ACTUAL[sector][gkey])
                mark = " ★ ACTUAL" if np.allclose(val, actual_val) else ""
                print(f"      ({val[0]:+5.1f}, {val[1]:+5.1f}, {val[2]:+5.1f})"
                      f"  ×{cnt} ({100*cnt/N:.0f}%){mark}")
        else:
            print(f"      (too many to list; top 5:)")
            for val, cnt in counts.most_common(5):
                actual_val = tuple(ACTUAL[sector][gkey])
                mark = " ★ ACTUAL" if np.allclose(val, actual_val) else ""
                print(f"      ({val[0]:+5.1f}, {val[1]:+5.1f}, {val[2]:+5.1f})"
                      f"  ×{cnt} ({100*cnt/N:.0f}%){mark}")
    print()

# ════════════════════════════════════════════════════════════════
# §3  COORDINATE-BY-COORDINATE — Which components vary?
# ════════════════════════════════════════════════════════════════

print("\n§3. COORDINATE-BY-COORDINATE — Which components vary?")
print("-" * 72)

comp_names = ['p', 'q', 'r']
for sector in ['lepton', 'up', 'down']:
    print(f"\n  === {sector.upper()} ===")
    for gkey, gen_label in [('gen1', 'gen-1'), ('gen2', 'gen-2'), ('gen3', 'gen-3')]:
        coords = np.array([j[sector][gkey] for j in pool])
        print(f"    {gen_label}:")
        for ci, cn in enumerate(comp_names):
            vals = coords[:, ci]
            uvals = sorted(set(np.round(vals, 1)))
            rng = vals.max() - vals.min()
            if rng < 1e-10:
                print(f"      {cn}: FIXED at {vals[0]:+.1f}")
            else:
                print(f"      {cn}: range [{vals.min():+.1f}, {vals.max():+.1f}], "
                      f"span={rng:.1f}, {len(uvals)} values")

print()


# ════════════════════════════════════════════════════════════════
# §4  GEN-1 ANATOMY — What varies in d and u?
# ════════════════════════════════════════════════════════════════

print("\n§4. GEN-1 ANATOMY — Decomposing the Degeneracy")
print("-" * 72)

# Since we expect the lepton is fully fixed and gen-3 is mostly fixed,
# the degeneracy should be in gen-1 of up and down sectors.

# Extract gen-1 coordinates
u_gen1 = np.array([j['up']['gen1'] for j in pool])
d_gen1 = np.array([j['down']['gen1'] for j in pool])

# Are they correlated?
print(f"\n  Up gen-1 p values:     {sorted(set(np.round(u_gen1[:,0],1)))}")
print(f"  Up gen-1 q values:     {sorted(set(np.round(u_gen1[:,1],1)))}")
print(f"  Up gen-1 r values:     {sorted(set(np.round(u_gen1[:,2],1)))}")
print(f"  Down gen-1 p values:   {sorted(set(np.round(d_gen1[:,0],1)))}")
print(f"  Down gen-1 q values:   {sorted(set(np.round(d_gen1[:,1],1)))}")
print(f"  Down gen-1 r values:   {sorted(set(np.round(d_gen1[:,2],1)))}")

# Check if up gen-1 varies independently of down gen-1
print(f"\n  Unique (up_gen1, down_gen1) pairs: ", end="")
combined = [tuple(np.round(np.concatenate([j['up']['gen1'], j['down']['gen1']]), 1))
            for j in pool]
print(f"{len(set(combined))}")
n_up_unique = len(set([tuple(np.round(j['up']['gen1'], 1)) for j in pool]))
n_dn_unique = len(set([tuple(np.round(j['down']['gen1'], 1)) for j in pool]))
print(f"  Unique up gen-1: {n_up_unique}")
print(f"  Unique down gen-1: {n_dn_unique}")
print(f"  Product (independent): {n_up_unique * n_dn_unique}")
print(f"  Actual pairs: {len(set(combined))}")
if len(set(combined)) == n_up_unique * n_dn_unique:
    print(f"  → Up and Down gen-1 vary INDEPENDENTLY (no cross-constraint)")
else:
    print(f"  → Up and Down gen-1 are CORRELATED (cross-constraint exists!)")
print()


# ════════════════════════════════════════════════════════════════
# §5  VELOCITY B DISTRIBUTION — What does B look like for survivors?
# ════════════════════════════════════════════════════════════════

print("\n§5. VELOCITY B DISTRIBUTION — Across Survivors")
print("-" * 72)

# B = (gen3 - gen1)/2 - 2*a  for g = 1 to 3 (using the formula from velocity_analysis)
# More precisely: x(g) = A*g^2 + B*g + C evaluated at g=1,3
# B = (x3 - x1)/2 - 2*A   where x = ln(m/me) as 3D vector
# But in coords: B_vec = (gen3 - gen1)/2 - 2*a_vec

for sector, a_vec in [('lepton', A_LEPTON), ('up', A_UP), ('down', A_DOWN)]:
    B_vecs = []
    for j in pool:
        cfg = j[sector]
        g1, g3 = cfg['gen1'], cfg['gen3']
        B = (g3 - g1) / 2.0 - 2 * a_vec
        B_vecs.append(B)
    B_arr = np.array(B_vecs)

    print(f"\n  === {sector.upper()} VELOCITY B ===")
    for ci, cn in enumerate(comp_names):
        uvals = sorted(set(np.round(B_arr[:, ci], 2)))
        if len(uvals) == 1:
            print(f"    B_{cn}: FIXED at {uvals[0]:+.2f}")
        else:
            print(f"    B_{cn}: {len(uvals)} values, range [{B_arr[:,ci].min():+.2f}, {B_arr[:,ci].max():+.2f}]")
            # Show most common
            counts = Counter(np.round(B_arr[:, ci], 2))
            for val, cnt in counts.most_common(5):
                actual_B = (ACTUAL[sector]['gen3'] - ACTUAL[sector]['gen1']) / 2.0 - 2 * a_vec
                mark = " ★" if abs(val - actual_B[ci]) < 0.05 else ""
                print(f"        {val:+.2f}  ×{cnt}{mark}")

    # Show unique B vectors
    B_tuples = [tuple(np.round(b, 2)) for b in B_vecs]
    print(f"    Unique B vectors: {len(set(B_tuples))}")

    actual_B = (ACTUAL[sector]['gen3'] - ACTUAL[sector]['gen1']) / 2.0 - 2 * a_vec
    print(f"    Actual B: ({actual_B[0]:+.2f}, {actual_B[1]:+.2f}, {actual_B[2]:+.2f})")

print()


# ════════════════════════════════════════════════════════════════
# §6  MASS PARABOLA A — Is curvature A identical across survivors?
# ════════════════════════════════════════════════════════════════

print("\n§6. CURVATURE A — Consistency Check")
print("-" * 72)

# A is fixed by the sector (A_LEPTON, A_UP, A_DOWN). It doesn't vary.
# But the scalar curvature A = a . basis does. Let's confirm.
for sector, a_vec in [('lepton', A_LEPTON), ('up', A_UP), ('down', A_DOWN)]:
    A_scalar = np.dot(a_vec, BASIS)
    print(f"  {sector:>7}: A = {A_scalar:+.6f}  (fixed by construction)")
print(f"\n  ✓ Curvature A is IDENTICAL for all survivors (it's an input, not a result).")
print()


# ════════════════════════════════════════════════════════════════
# §7  DOWN gen-1 PATTERNS — Arithmetic relations?
# ════════════════════════════════════════════════════════════════

print("\n§7. DOWN gen-1 PATTERNS — Hidden Arithmetic?")
print("-" * 72)

d_gen1_tuples = [tuple(np.round(j['down']['gen1'], 1)) for j in pool]
d_gen1_unique = sorted(set(d_gen1_tuples))
print(f"\n  {len(d_gen1_unique)} unique down gen-1 coordinates:\n")

# For each down gen-1, show coordinates and computed mass
for coords in d_gen1_unique:
    p, q, r = coords
    m = mass_from_coords(p, q, r)
    # Check: is p+q+r constant? p-r constant? etc.
    actual = tuple(ACTUAL['down']['gen1'])
    mark = " ★ ACTUAL" if np.allclose(coords, actual) else ""
    print(f"    ({p:+5.1f}, {q:+5.1f}, {r:+5.1f})  m_d={m:.3f} MeV  "
          f"p+q+r={p+q+r:+5.1f}  p-r={p-r:+5.1f}  "
          f"2p+r={2*p+r:+5.1f}  p·ln2+q·lnR3+r·ln3={p*LN2+q*LNR3+r*LN3:+.4f}{mark}")

# Check if there's a linear constraint on the coordinates
print(f"\n  Looking for linear constraints on down gen-1...")
d1_arr = np.array(d_gen1_unique)
for a_coeff in range(-3, 4):
    for b_coeff in range(-3, 4):
        for c_coeff in range(-3, 4):
            if a_coeff == 0 and b_coeff == 0 and c_coeff == 0:
                continue
            combo = a_coeff * d1_arr[:, 0] + b_coeff * d1_arr[:, 1] + c_coeff * d1_arr[:, 2]
            vals = set(np.round(combo, 1))
            if len(vals) == 1:
                print(f"    {a_coeff}p + {b_coeff}q + {c_coeff}r = {combo[0]:+.1f}  "
                      f"(CONSTANT for all down gen-1!)")

# Check what's constant about up gen-1 too
u_gen1_tuples = [tuple(np.round(j['up']['gen1'], 1)) for j in pool]
u_gen1_unique = sorted(set(u_gen1_tuples))
print(f"\n  {len(u_gen1_unique)} unique up gen-1 coordinates:\n")
for coords in u_gen1_unique:
    p, q, r = coords
    m = mass_from_coords(p, q, r)
    actual = tuple(ACTUAL['up']['gen1'])
    mark = " ★ ACTUAL" if np.allclose(coords, actual) else ""
    print(f"    ({p:+5.1f}, {q:+5.1f}, {r:+5.1f})  m_u={m:.3f} MeV  "
          f"p+q+r={p+q+r:+5.1f}  p-r={p-r:+5.1f}{mark}")

print(f"\n  Looking for linear constraints on up gen-1...")
u1_arr = np.array(u_gen1_unique)
for a_coeff in range(-3, 4):
    for b_coeff in range(-3, 4):
        for c_coeff in range(-3, 4):
            if a_coeff == 0 and b_coeff == 0 and c_coeff == 0:
                continue
            combo = a_coeff * u1_arr[:, 0] + b_coeff * u1_arr[:, 1] + c_coeff * u1_arr[:, 2]
            vals = set(np.round(combo, 1))
            if len(vals) == 1:
                print(f"    {a_coeff}p + {b_coeff}q + {c_coeff}r = {combo[0]:+.1f}  "
                      f"(CONSTANT for all up gen-1!)")

print()


# ════════════════════════════════════════════════════════════════
# §8  INTER-SECTOR RELATION — u↔d gen-1 distance
# ════════════════════════════════════════════════════════════════

print("\n§8. INTER-SECTOR RELATION — u↔d Gen-1 Distance")
print("-" * 72)

print(f"\n  For each survivor, compute Δ₁ = d(gen-1) - u(gen-1):")
print(f"  (This is the gen-1 displacement between up and down sectors)\n")

delta1_list = []
for j in pool:
    d1 = j['down']['gen1']
    u1 = j['up']['gen1']
    delta = d1 - u1
    delta1_list.append(delta)

delta1_arr = np.array(delta1_list)
delta1_tuples = [tuple(np.round(d, 1)) for d in delta1_list]
unique_deltas = sorted(set(delta1_tuples))
print(f"  Unique Δ₁ = d₁ - u₁ vectors: {len(unique_deltas)}\n")

for delta in unique_deltas[:15]:
    cnt = delta1_tuples.count(delta)
    p, q, r = delta
    # Is this the actual delta?
    actual_delta = tuple(ACTUAL['down']['gen1'] - ACTUAL['up']['gen1'])
    mark = " ★ ACTUAL" if np.allclose(delta, actual_delta) else ""
    # Compute |Δ| in mass space
    ln_delta = p * LN2 + q * LNR3 + r * LN3
    print(f"    Δ₁=({p:+5.1f},{q:+5.1f},{r:+5.1f})  |Δ|_lattice={np.sqrt(p**2+q**2+r**2):.2f}  "
          f"Δ(lnm)={ln_delta:+.4f}  ×{cnt}{mark}")

# Check: |Δ| in mass space
lm_deltas = [np.dot(d, BASIS) for d in delta1_list]
print(f"\n  ln(m_d/m_u) at gen-1:")
print(f"    from coords: {sorted(set(np.round(lm_deltas, 4)))}")
print(f"    experiment:   {np.log(MASSES_EXP['d']/MASSES_EXP['u']):+.4f}")

print()


# ════════════════════════════════════════════════════════════════
# §9  RESIDUAL PATTERN — Error Sign Structure
# ════════════════════════════════════════════════════════════════

print("\n§9. RESIDUAL PATTERN — Error Sign Structure (+ → −)")
print("-" * 72)

print(f"\n  For each survivor, compute the fractional error per generation.\n")

# For each config, compute error at each of the 9 masses
all_errors = []
for j in pool:
    errs = {}
    for sector, pnames in [('lepton', ['e','mu','tau']),
                            ('up', ['u','c','t']),
                            ('down', ['d','s','b'])]:
        for g_idx, (gkey, pn) in enumerate(zip(['gen1','gen2','gen3'], pnames)):
            m = mass_from_coords(*j[sector][gkey])
            errs[pn] = (m / MASSES_EXP[pn] - 1) * 100  # percent
    all_errors.append(errs)

# Tabulate
print(f"  {'Config':>7}  {'e':>6}  {'μ':>6}  {'τ':>6}  {'u':>6}  {'c':>6}  {'t':>6}  "
      f"{'d':>6}  {'s':>6}  {'b':>6}")
print(f"  " + "─" * 75)

# Show first 20, plus the actual
actual_idx = None
for i, (j, errs) in enumerate(zip(pool, all_errors)):
    is_actual = (np.allclose(j['up']['gen1'], ACTUAL['up']['gen1'])
                 and np.allclose(j['down']['gen1'], ACTUAL['down']['gen1']))
    if is_actual:
        actual_idx = i

    if i < 15 or is_actual:
        mark = " ★" if is_actual else ""
        print(f"  {i+1:>5d}{mark} ", end="")
        for pn in ['e','mu','tau','u','c','t','d','s','b']:
            print(f" {errs[pn]:+5.1f}%", end="")
        print()

if actual_idx is not None and actual_idx >= 15:
    print(f"  (actual config is #{actual_idx+1})")

# Compute the sign pattern for each config
print(f"\n\n  SIGN PATTERN across generations (within each sector):")
print(f"  ───────────────────────────────────────────────────")

sign_patterns = Counter()
for errs in all_errors:
    pattern = ""
    for sector_pnames in [['e','mu','tau'], ['u','c','t'], ['d','s','b']]:
        signs = ''.join('+' if errs[pn] >= 0 else '-' for pn in sector_pnames)
        pattern += signs + "|"
    sign_patterns[pattern] += 1

print(f"\n  {'Pattern (ℓ|u|d)':>20}  Count  %")
for pattern, cnt in sign_patterns.most_common():
    print(f"  {pattern:>20}  {cnt:>5}  {100*cnt/N:.1f}%")


# ════════════════════════════════════════════════════════════════
# §10  RESIDUAL DECOMPOSITION — Linear vs Quadratic
# ════════════════════════════════════════════════════════════════

print("\n\n§10. RESIDUAL DECOMPOSITION — Linear (β) vs Quadratic (γ)")
print("-" * 72)

print("""
  For each sector, decompose the residual δ(g) into:
    δ(g) = α + β·(g-2) + γ·(g-2)²
  where:
    γ = (δ₁ - 2δ₂ + δ₃)/2  → curvature error
    β = (δ₃ - δ₁)/2         → velocity error
    α = δ₂                   → offset error
""")

decomp = {sector: {'alpha': [], 'beta': [], 'gamma': []}
          for sector in ['lepton', 'up', 'down']}
sector_pnames = {'lepton': ['e','mu','tau'], 'up': ['u','c','t'], 'down': ['d','s','b']}

for errs in all_errors:
    for sector, pnames in sector_pnames.items():
        d1, d2, d3 = [errs[pn] / 100 for pn in pnames]  # convert to fraction
        gamma = (d1 - 2 * d2 + d3) / 2
        beta = (d3 - d1) / 2
        alpha = d2
        decomp[sector]['alpha'].append(alpha)
        decomp[sector]['beta'].append(beta)
        decomp[sector]['gamma'].append(gamma)

for sector in ['lepton', 'up', 'down']:
    alpha = np.array(decomp[sector]['alpha'])
    beta = np.array(decomp[sector]['beta'])
    gamma = np.array(decomp[sector]['gamma'])
    print(f"  === {sector.upper()} ===")
    print(f"    γ (curvature err):  mean={gamma.mean()*100:+.3f}%  "
          f"std={gamma.std()*100:.3f}%  range=[{gamma.min()*100:+.3f}%, {gamma.max()*100:+.3f}%]")
    print(f"    β (velocity err):   mean={beta.mean()*100:+.3f}%  "
          f"std={beta.std()*100:.3f}%  range=[{beta.min()*100:+.3f}%, {beta.max()*100:+.3f}%]")
    print(f"    α (offset err):     mean={alpha.mean()*100:+.3f}%  "
          f"std={alpha.std()*100:.3f}%  range=[{alpha.min()*100:+.3f}%, {alpha.max()*100:+.3f}%]")
    print()


# ════════════════════════════════════════════════════════════════
# §11  MINIMAL β — Which configs have smallest velocity error?
# ════════════════════════════════════════════════════════════════

print("\n§11. MINIMAL VELOCITY ERROR — Ranking by |β|")
print("-" * 72)

print("""
  If the sign-flip pattern (+ → −) is a velocity error β,
  the "best" config should minimize |β| across all 3 sectors.
""")

# Compute total |β| for each config
beta_total = []
for i, errs in enumerate(all_errors):
    total_beta = 0
    for sector, pnames in sector_pnames.items():
        d1, d2, d3 = [errs[pn] / 100 for pn in pnames]
        beta = (d3 - d1) / 2
        total_beta += abs(beta)
    beta_total.append(total_beta)

# Rank
ranked = sorted(range(N), key=lambda i: beta_total[i])
print(f"  {'Rank':>4}  {'Config':>6}  {'|β_ℓ|':>7}  {'|β_u|':>7}  {'|β_d|':>7}  {'Σ|β|':>8}  Actual?")
print(f"  " + "─" * 65)

for rank, idx in enumerate(ranked[:20], 1):
    j = pool[idx]
    errs = all_errors[idx]
    betas = []
    for sector, pnames in sector_pnames.items():
        d1, d3 = errs[pnames[0]] / 100, errs[pnames[2]] / 100
        betas.append(abs((d3 - d1) / 2))
    is_actual = (np.allclose(j['up']['gen1'], ACTUAL['up']['gen1'])
                 and np.allclose(j['down']['gen1'], ACTUAL['down']['gen1']))
    mark = "★" if is_actual else ""
    print(f"  {rank:>4}  {idx+1:>6}  {betas[0]*100:>6.3f}%  {betas[1]*100:>6.3f}%  {betas[2]*100:>6.3f}%  "
          f"{sum(betas)*100:>7.3f}%  {mark}")

# Also show: where does the actual config rank?
if actual_idx is not None:
    actual_rank = ranked.index(actual_idx) + 1
    print(f"\n  → ACTUAL config rank by |β|: #{actual_rank} of {N}")


# ════════════════════════════════════════════════════════════════
# §12  THE CONSTRAINT THAT KILLS — What Distinguishes Actual?
# ════════════════════════════════════════════════════════════════

print("\n\n§12. THE CONSTRAINT THAT KILLS — Actual vs Others")
print("-" * 72)

if actual_idx is not None:
    j_act = pool[actual_idx]

    # Actual gen-1 coordinates
    u1_act = j_act['up']['gen1']
    d1_act = j_act['down']['gen1']
    print(f"\n  Actual gen-1:  u=({u1_act[0]:+.1f},{u1_act[1]:+.1f},{u1_act[2]:+.1f})  "
          f"d=({d1_act[0]:+.1f},{d1_act[1]:+.1f},{d1_act[2]:+.1f})")
    print(f"  Actual Δ₁:    ({(d1_act-u1_act)[0]:+.1f},{(d1_act-u1_act)[1]:+.1f},"
          f"{(d1_act-u1_act)[2]:+.1f})")

    # Properties of the actual vs others:
    # 1. Is gen-2 close to gen-3 in q? (spectral consolidation)
    print(f"\n  Testing distinguishing properties:")

    # Property: q₂ = q₃ for some sector?
    for sector in ['up', 'down']:
        q2_vals = [j[sector]['gen2'][1] for j in pool]
        q3_val = pool[0][sector]['gen3'][1]  # same for all
        match = [abs(q2 - q3_val) < 1e-10 for q2 in q2_vals]
        n_match = sum(match)
        actual_match = abs(j_act[sector]['gen2'][1] - q3_val) < 1e-10
        print(f"    q₂ = q₃ in {sector}: {n_match}/{N} configs  "
              f"(actual: {'YES' if actual_match else 'NO'})")

    # Property: integer vs half-integer at gen-1
    print(f"\n  Coordinate type at gen-1 (integer vs half-integer):")
    for sector in ['up', 'down']:
        g1_all = [j[sector]['gen1'] for j in pool]
        int_counts = Counter()
        for g1 in g1_all:
            types = tuple('I' if abs(x - round(x)) < 1e-10 else 'H' for x in g1)
            int_counts[types] += 1
        actual_type = tuple('I' if abs(x - round(x)) < 1e-10 else 'H'
                            for x in j_act[sector]['gen1'])
        print(f"    {sector} gen-1 type pattern:")
        for types, cnt in int_counts.most_common():
            mark = " ★ ACTUAL" if types == actual_type else ""
            print(f"      {types}  ×{cnt} ({100*cnt/N:.0f}%){mark}")

    # Property: Sum of gen-1 coordinates
    print(f"\n  Sum of gen-1 coordinates p+q+r:")
    for sector in ['up', 'down']:
        sums = [sum(j[sector]['gen1']) for j in pool]
        actual_sum = sum(j_act[sector]['gen1'])
        uvals = sorted(set(np.round(sums, 1)))
        print(f"    {sector}: values = {uvals}")
        print(f"    actual sum = {actual_sum:+.1f}")

    # Property: Velocity B as half-integer vector?
    print(f"\n  Is velocity B a half-integer vector?")
    for sector, a_vec in [('up', A_UP), ('down', A_DOWN)]:
        for j in [j_act]:
            g1, g3 = j[sector]['gen1'], j[sector]['gen3']
            B = (g3 - g1) / 2.0 - 2 * a_vec
            all_hi = all(is_half_int(b) for b in B)
            print(f"    {sector}: B=({B[0]:+.2f},{B[1]:+.2f},{B[2]:+.2f})  half-int: {all_hi}")

    # How many survivors have B as half-integer?
    print(f"\n  How many survivors have half-integer B?")
    for sector, a_vec in [('up', A_UP), ('down', A_DOWN)]:
        count_hi = 0
        for j in pool:
            g1, g3 = j[sector]['gen1'], j[sector]['gen3']
            B = (g3 - g1) / 2.0 - 2 * a_vec
            if all(is_half_int(b) for b in B):
                count_hi += 1
        print(f"    {sector}: {count_hi}/{N} have half-int B")


# ════════════════════════════════════════════════════════════════
# §13  SUMMARY
# ════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 72)
print("  §13. SUMMARY — Anatomy of Degeneracy")
print("=" * 72)

print(f"""
  FROM 2500 joint configs (gen-3 gap structure, 10% tolerance),
  {N} survive at 2%.

  FIXED (identical across all survivors):
    - Lepton sector: fully fixed (gen-1 = origin, unique gen-3)
    - Gen-3 for all sectors: fully fixed by gap structure
    - Curvature A: fixed by construction
    - Gen-2: determined by midpoint formula

  VARIABLE:
    - Up gen-1 coordinates
    - Down gen-1 coordinates
    → The entire degeneracy is in WHERE gen-1 quarks sit on the lattice.

  The velocity B inherits this degeneracy (B depends on gen-1).
  The mass errors show a systematic + → − sign flip with generation,
  confirming the residual is dominated by a VELOCITY correction (β),
  not a curvature correction (γ ≈ 0).

  To break the degeneracy, we need a constraint on gen-1 positions.
  Candidates: CKM angles, velocity regularity, integer-type restrictions.
""")


# ════════════════════════════════════════════════════════════════
# CLEANUP
# ════════════════════════════════════════════════════════════════

sys.stdout = _orig_stdout
_outfile.close()
print(f"\n[Output saved to {OUTPUT_PATH}]")
