#!/usr/bin/env python3
"""
ckm_geometry.py — CKM from curve geometry + q₂=q₃ constraint
=============================================================

1. Impose q₂ = q₃ in the up sector → reduce 72 → 9 configs
2. Study inter-sector geometry (tangent angles, distances)
3. Test whether CKM angles (Cabibbo = √2/(2π)) select the actual config
4. Look for additional constraints that break the down-sector degeneracy

Author: SS-PEG program
Date: April 2026
"""

import numpy as np
import sys, os
from collections import Counter

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "ckm_geometry_output.md")
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
# CONSTANTS
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

CABIBBO_PRED = np.sqrt(2) / (2 * np.pi)  # ≈ 0.22508

# CKM experimental values (PDG 2024, magnitudes)
CKM_EXP = {
    'V_ud': 0.97373, 'V_us': 0.2243, 'V_ub': 0.00382,
    'V_cd': 0.221,   'V_cs': 0.975,  'V_cb': 0.0408,
    'V_td': 0.0080,  'V_ts': 0.0388, 'V_tb': 0.9991,
}

def mass_from_coords(p, q, r):
    return m_e * np.exp(p * LN2 + q * LNR3 + r * LN3)

def is_half_int(x, tol=1e-10):
    return abs(2 * x - round(2 * x)) < tol

def coords_half_int(arr, tol=1e-10):
    return all(is_half_int(x, tol) for x in arr)

def is_integer(x, tol=1e-10):
    return abs(x - round(x)) < tol

def coords_integer(arr, tol=1e-10):
    return all(is_integer(x, tol) for x in arr)

# ════════════════════════════════════════════════════════════════
# REBUILD 72 SURVIVORS (same logic as survivor_analysis.py)
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


# Build pool
lepton_configs = search_sector(A_LEPTON, gen1_fixed=np.array([0.0, 0.0, 0.0]),
                               gen3_q_fixed=-7.0, m3_exp=MASSES_EXP['tau'], tol=10.0)
up_configs = search_sector(A_UP, gen3_q_fixed=-7.0,
                           m1_exp=MASSES_EXP['u'], m3_exp=MASSES_EXP['t'], tol=0.5)
down_configs = search_sector(A_DOWN,
                             m1_exp=MASSES_EXP['d'], m3_exp=MASSES_EXP['b'], tol=0.5)

lepton_10 = [c for c in lepton_configs
             if abs(c['masses'][1] - MASSES_EXP['mu']) / MASSES_EXP['mu'] < 0.10
             and abs(c['masses'][2] - MASSES_EXP['tau']) / MASSES_EXP['tau'] < 0.10]
up_10 = [c for c in up_configs
         if all(abs(c['masses'][i] - m) / m < 0.10
                for i, m in enumerate([MASSES_EXP['u'], MASSES_EXP['c'], MASSES_EXP['t']]))]
down_10 = [c for c in down_configs
           if all(abs(c['masses'][i] - m) / m < 0.10
                  for i, m in enumerate([MASSES_EXP['d'], MASSES_EXP['s'], MASSES_EXP['b']]))]

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
                    ok = False; break
            if not ok: break
        if ok: survivors.append(j)
    return survivors

s2 = filter_joint(joint_all, 0.02)

# ════════════════════════════════════════════════════════════════
print("=" * 72)
print("  CKM GEOMETRY — From 72 to the Actual Config")
print("=" * 72)
print()
print(f"  Starting pool: {len(s2)} configs at 2% tolerance")
print()


# ════════════════════════════════════════════════════════════════
# §1  IMPOSE q₂ = q₃ IN UP SECTOR
# ════════════════════════════════════════════════════════════════

print("§1. CONSTRAINT: q₂ = q₃ in UP sector (q(c) = q(t) = -7)")
print("-" * 72)

# For up sector, gen-3 q is fixed at -7.0. Require gen-2 q = -7.0 too.
s2_q2eq3 = [j for j in s2 if abs(j['up']['gen2'][1] - j['up']['gen3'][1]) < 1e-10]

print(f"\n  Before: {len(s2)} configs")
print(f"  After q₂ = q₃ (up): {len(s2_q2eq3)} configs")
print(f"  Reduction: {len(s2)} → {len(s2_q2eq3)} ({100*len(s2_q2eq3)/len(s2):.1f}%)")

# Verify: all share the SAME up sector
up_g1_strs = set()
for j in s2_q2eq3:
    up_g1_strs.add(tuple(np.round(j['up']['gen1'], 1)))
print(f"  Unique up gen-1 in survivors: {len(up_g1_strs)}")
if len(up_g1_strs) == 1:
    print(f"  → UP SECTOR IS FULLY FIXED: gen-1 = {list(up_g1_strs)[0]}")
    actual_up_g1 = tuple(ACTUAL['up']['gen1'])
    if np.allclose(list(up_g1_strs)[0], actual_up_g1):
        print(f"    ✓ Matches the actual up gen-1!")
    else:
        print(f"    ✗ Does NOT match actual up gen-1 = {actual_up_g1}")

# Check actual is present
actual_present = False
for j in s2_q2eq3:
    if (np.allclose(j['up']['gen1'], ACTUAL['up']['gen1']) and
        np.allclose(j['down']['gen1'], ACTUAL['down']['gen1'])):
        actual_present = True
        break
print(f"\n  Actual config present: {actual_present}")
print()


# ════════════════════════════════════════════════════════════════
# §2  THE REMAINING 9 — Down gen-1 Manifold
# ════════════════════════════════════════════════════════════════

print("\n§2. THE 9 SURVIVORS — Down gen-1 Manifold")
print("-" * 72)

pool9 = s2_q2eq3
N9 = len(pool9)

print(f"\n  All {N9} configs share identical lepton and up sectors.")
print(f"  The ONLY variable is the down gen-1 coordinate.\n")

# List all 9 with coordinates, masses, errors, velocity B
print(f"  {'#':>3}  {'d gen-1 (p,q,r)':^24}  {'m_d MeV':>8}  {'m_s MeV':>8}  "
      f"{'m_b MeV':>8}  {'err_d':>6}  {'err_s':>6}  {'err_b':>6}  "
      f"{'B_d (p,q,r)':^24}  {'B all-int?':>10}")
print("  " + "─" * 150)

for idx, j in enumerate(pool9):
    d1 = j['down']['gen1']
    d2 = j['down']['gen2']
    d3 = j['down']['gen3']

    m_d = mass_from_coords(*d1)
    m_s = mass_from_coords(*d2)
    m_b = mass_from_coords(*d3)

    err_d = (m_d / MASSES_EXP['d'] - 1) * 100
    err_s = (m_s / MASSES_EXP['s'] - 1) * 100
    err_b = (m_b / MASSES_EXP['b'] - 1) * 100

    B = (d3 - d1) / 2.0 - 2 * A_DOWN
    all_int = coords_integer(B)

    actual_d = np.allclose(d1, ACTUAL['down']['gen1'])
    mark = "★" if actual_d else " "

    print(f"  {idx+1:>2}{mark} ({d1[0]:+5.1f},{d1[1]:+5.1f},{d1[2]:+5.1f})  "
          f"{m_d:8.3f}  {m_s:8.2f}  {m_b:8.1f}  "
          f"{err_d:+5.1f}%  {err_s:+5.1f}%  {err_b:+5.1f}%  "
          f"({B[0]:+5.2f},{B[1]:+5.2f},{B[2]:+5.2f})  "
          f"{'YES' if all_int else 'no':>10}")

print()

# Count how many have all-integer B
n_int_B = sum(1 for j in pool9 if coords_integer(
    (j['down']['gen3'] - j['down']['gen1']) / 2.0 - 2 * A_DOWN))
print(f"  Down B all-integer: {n_int_B}/{N9}")

# Check: is "all-integer B" unique?
if n_int_B == 1:
    print(f"  → ALL-INTEGER B IS UNIQUE — selects exactly the actual config!")
elif n_int_B == 0:
    print(f"  → No config has all-integer B")
else:
    print(f"  → {n_int_B} configs have all-integer B (not unique)")
print()


# ════════════════════════════════════════════════════════════════
# §3  VELOCITY STRUCTURE — B Properties
# ════════════════════════════════════════════════════════════════

print("\n§3. VELOCITY B STRUCTURE — Detailed Properties")
print("-" * 72)

# For each surviving config, compute down B
B_list = []
for j in pool9:
    d1, d3 = j['down']['gen1'], j['down']['gen3']
    B = (d3 - d1) / 2.0 - 2 * A_DOWN
    B_list.append(B)

B_arr = np.array(B_list)

# Systematic: what integer/half-int properties?
print(f"\n  For each down gen-1, the velocity B = (d₃ - d₁)/2 - 2A_d:")
print(f"  (A_d = {tuple(A_DOWN)})\n")

for idx, (j, B) in enumerate(zip(pool9, B_list)):
    d1 = j['down']['gen1']
    actual_d = np.allclose(d1, ACTUAL['down']['gen1'])
    mark = "★" if actual_d else " "

    # Component types
    types = []
    for v in B:
        if is_integer(v):
            types.append('Z')  # integer
        elif is_half_int(v):
            types.append('H')  # half-integer
        else:
            types.append(f'{v:.3f}')
    type_str = ','.join(types)

    # |B| in mass space
    B_scalar = np.dot(B, BASIS)
    B_norm_latt = np.linalg.norm(B)

    print(f"  {idx+1:>2}{mark} B=({B[0]:+5.2f},{B[1]:+5.2f},{B[2]:+5.2f})  "
          f"type=({type_str})  |B|_latt={B_norm_latt:.3f}  B·basis={B_scalar:+.4f}  "
          f"sum={sum(B):+.2f}")

print()


# ════════════════════════════════════════════════════════════════
# §4  TANGENT VECTORS AT EACH GENERATION
# ════════════════════════════════════════════════════════════════

print("\n§4. TANGENT VECTORS — Parabolic Curves x(g) = A·g² + B·g + C")
print("-" * 72)

# Tangent: v(g) = dx/dg = 2A·g + B
# At gen g: v(g) = 2·a_vec·g + B
# We care about up and down tangent vectors at each generation

# Up is now fixed (all 9 have same up):
up_cfg = pool9[0]['up']
up_B = (up_cfg['gen3'] - up_cfg['gen1']) / 2.0 - 2 * A_UP

print(f"\n  UP SECTOR (fixed for all 9):")
print(f"    A_up = ({A_UP[0]:+.2f}, {A_UP[1]:+.2f}, {A_UP[2]:+.2f})")
print(f"    B_up = ({up_B[0]:+.2f}, {up_B[1]:+.2f}, {up_B[2]:+.2f})")
for g in [1, 2, 3]:
    v = 2 * A_UP * g + up_B
    print(f"    v_up(g={g}) = ({v[0]:+.2f}, {v[1]:+.2f}, {v[2]:+.2f})  "
          f"|v|={np.linalg.norm(v):.3f}")

print()

# For each down config, compute tangent vectors and angle with up
print(f"  DOWN SECTOR (varies across 9):")
print(f"    A_down = ({A_DOWN[0]:+.2f}, {A_DOWN[1]:+.2f}, {A_DOWN[2]:+.2f})")
print()

for g in [1, 2, 3]:
    gen_label = ['', '1st', '2nd', '3rd'][g]
    v_up = 2 * A_UP * g + up_B
    v_up_norm = np.linalg.norm(v_up)

    print(f"  === Generation {g} ({gen_label}) ===")
    print(f"  {'#':>3}  {'v_down(g)':^24}  {'|v_d|':>7}  {'cos θ':>8}  {'θ (deg)':>8}  {'sin θ':>8}")
    print("  " + "─" * 80)

    for idx, (j, B_d) in enumerate(zip(pool9, B_list)):
        v_down = 2 * A_DOWN * g + B_d
        v_d_norm = np.linalg.norm(v_down)
        cos_theta = np.dot(v_up, v_down) / (v_up_norm * v_d_norm)
        cos_theta = np.clip(cos_theta, -1, 1)
        theta = np.arccos(cos_theta)
        sin_theta = np.sin(theta)

        actual_d = np.allclose(j['down']['gen1'], ACTUAL['down']['gen1'])
        mark = "★" if actual_d else " "

        print(f"  {idx+1:>2}{mark} ({v_down[0]:+5.2f},{v_down[1]:+5.2f},{v_down[2]:+5.2f})  "
              f"{v_d_norm:7.3f}  {cos_theta:+8.5f}  {np.degrees(theta):8.3f}°  "
              f"{sin_theta:8.5f}")

    # Does sin θ ≈ Cabibbo for any config?
    print(f"\n  Cabibbo prediction: sin θ_C = √2/(2π) ≈ {CABIBBO_PRED:.5f}")
    print(f"  → θ_C ≈ {np.degrees(np.arcsin(CABIBBO_PRED)):.3f}°\n")


# ════════════════════════════════════════════════════════════════
# §5  CKM FROM ROTATION BETWEEN CURVE FRAMES
# ════════════════════════════════════════════════════════════════

print("\n§5. CKM FROM ROTATION BETWEEN UP AND DOWN FRAMES")
print("-" * 72)

# At each generation g, we have a tangent vector v(g) and an acceleration
# a_vec (constant). Together with v×a they form a Frenet-like frame.
# The CKM matrix could be the rotation between up and down frames.
#
# But our tangents are 3D (p,q,r coordinates), and the CKM is a 3×3
# rotation in generation space. Let's try a different approach:
#
# Define the "mass-space tangent" at each generation:
#   v_scalar(g) = v(g) · basis = (2A·g + B) · basis
#
# The mixing angle between up and down at generation g in mass space:
#   tan θ(g) = |v_u × v_d| / (v_u · v_d)  (but in 3D this is the exterior angle)

print(f"\n  APPROACH: Compare tangent vectors in MASS space (scalar projection)\n")

# Mass-space tangent: v_mass(g) = 2*A_scalar*g + B_scalar
A_up_s = np.dot(A_UP, BASIS)
B_up_s = np.dot(up_B, BASIS)
A_dn_s = np.dot(A_DOWN, BASIS)

print(f"  Up:   A_s = {A_up_s:+.6f},  B_s = {B_up_s:+.6f}")
print(f"  Down: A_s = {A_dn_s:+.6f}")
print()

for idx, (j, B_d) in enumerate(zip(pool9, B_list)):
    B_dn_s = np.dot(B_d, BASIS)
    actual_d = np.allclose(j['down']['gen1'], ACTUAL['down']['gen1'])
    mark = "★" if actual_d else " "

    # Mass-space tangent ratio at each generation
    print(f"  {idx+1:>2}{mark} B_d·basis = {B_dn_s:+.6f}")

    for g in [1, 2, 3]:
        v_u_s = 2 * A_up_s * g + B_up_s
        v_d_s = 2 * A_dn_s * g + B_dn_s
        ratio = v_d_s / v_u_s if abs(v_u_s) > 1e-15 else float('inf')
        diff = v_d_s - v_u_s
        print(f"      g={g}: v_u={v_u_s:+.4f}, v_d={v_d_s:+.4f}, "
              f"ratio={ratio:+.4f}, diff={diff:+.4f}")
    print()

print()


# ════════════════════════════════════════════════════════════════
# §6  CKM FROM GENERATION-SPACE ANGLES
# ════════════════════════════════════════════════════════════════

print("\n§6. CKM FROM GENERATION-SPACE ROTATION")
print("-" * 72)

# Alternative: The CKM matrix is the rotation between the up and down
# mass eigenstates in generation space. In our framework, the curve
# x(g) gives mass eigenvalues at g=1,2,3. The "generation-space"
# direction is the (normalized) vector of ln masses in gen space.
#
# Define for each sector the mass vector M = (ln m₁, ln m₂, ln m₃)
# The CKM rotation R satisfies: M_down = R · M_up (approximately)
# → R = M_down · M_up† (pseudoinverse in 3D)
#
# But more physically: the mixing angle θ₁₂ is the angle between
# up and down curves projected onto the gen-1↔gen-2 plane.

print(f"\n  For each config, compute the mass vectors:\n")

# Up masses (fixed)
up_lnm = np.array([np.log(mass_from_coords(*up_cfg[gk]) / m_e) for gk in ['gen1', 'gen2', 'gen3']])
print(f"  Up ln(m/mₑ):  ({up_lnm[0]:+.4f}, {up_lnm[1]:+.4f}, {up_lnm[2]:+.4f})")

# Lepton masses (fixed)
lep_cfg = pool9[0]['lepton']
lep_lnm = np.array([np.log(mass_from_coords(*lep_cfg[gk]) / m_e) for gk in ['gen1', 'gen2', 'gen3']])
print(f"  Lep ln(m/mₑ): ({lep_lnm[0]:+.4f}, {lep_lnm[1]:+.4f}, {lep_lnm[2]:+.4f})")

print()

# For down sector: build mass-ratio hierarchies and compare
print(f"  {'#':>3}  {'ln(m_d/mₑ)':>10}  {'ln(m_s/mₑ)':>10}  {'ln(m_b/mₑ)':>10}  "
      f"{'θ₁₂ (°)':>8}  {'sin θ₁₂':>9}  {'θ₂₃ (°)':>8}  {'sin θ₂₃':>9}  "
      f"{'θ₁₃ (°)':>8}  {'sin θ₁₃':>9}")
print("  " + "─" * 120)

for idx, j in enumerate(pool9):
    d_cfg = j['down']
    d_lnm = np.array([np.log(mass_from_coords(*d_cfg[gk]) / m_e) for gk in ['gen1', 'gen2', 'gen3']])

    # The CKM angle θ₁₂ relates to the rotation needed to go from
    # up mass hierarchy to down mass hierarchy in gen-1↔gen-2 subspace.
    #
    # Use the parametric approach: if up and down curves differ by a
    # rotation θ in generation space, then:
    #   sin θ₁₂ ≈ |Δm₁₂_d · Δm₁₂_u_perp| / |Δm₁₂_u|
    #
    # Simpler: use the mass ratios directly.
    # θ₁₂ = arctan(√(m_d/m_s) / √(m_u/m_c))  (Gatto relation variant)

    m_d_val = mass_from_coords(*d_cfg['gen1'])
    m_s_val = mass_from_coords(*d_cfg['gen2'])
    m_b_val = mass_from_coords(*d_cfg['gen3'])
    m_u_val = mass_from_coords(*up_cfg['gen1'])
    m_c_val = mass_from_coords(*up_cfg['gen2'])
    m_t_val = mass_from_coords(*up_cfg['gen3'])

    # Gatto relation: sin θ_C ≈ √(m_d/m_s)
    sin_12_gatto = np.sqrt(m_d_val / m_s_val)
    theta_12_gatto = np.arcsin(min(sin_12_gatto, 1.0))

    # Wolfenstein-like: sin θ₂₃ ≈ √(m_s/m_b)
    sin_23 = np.sqrt(m_s_val / m_b_val)
    theta_23 = np.arcsin(min(sin_23, 1.0))

    # sin θ₁₃ ≈ √(m_d/m_b)
    sin_13 = np.sqrt(m_d_val / m_b_val)
    theta_13 = np.arcsin(min(sin_13, 1.0))

    actual_d = np.allclose(d_cfg['gen1'], ACTUAL['down']['gen1'])
    mark = "★" if actual_d else " "

    print(f"  {idx+1:>2}{mark} {d_lnm[0]:+10.4f}  {d_lnm[1]:+10.4f}  {d_lnm[2]:+10.4f}  "
          f"{np.degrees(theta_12_gatto):8.3f}  {sin_12_gatto:9.5f}  "
          f"{np.degrees(theta_23):8.3f}  {sin_23:9.5f}  "
          f"{np.degrees(theta_13):8.3f}  {sin_13:9.5f}")

print(f"\n  Experimental CKM:")
print(f"    |V_us| = {CKM_EXP['V_us']:.4f}  (Cabibbo, sin θ₁₂)")
print(f"    |V_cb| = {CKM_EXP['V_cb']:.4f}  (sin θ₂₃)")
print(f"    |V_ub| = {CKM_EXP['V_ub']:.5f} (sin θ₁₃)")
print(f"    Cabibbo pred = √2/(2π) = {CABIBBO_PRED:.5f}")
print()


# ════════════════════════════════════════════════════════════════
# §7  CABIBBO CONSTRAINT — Does sin θ₁₂ ≈ √2/(2π) select one?
# ════════════════════════════════════════════════════════════════

print("\n§7. CABIBBO CONSTRAINT — Testing sin θ₁₂ ≈ Cabibbo")
print("-" * 72)

print(f"\n  Using Gatto relation: sin θ₁₂ = √(m_d/m_s)")
print(f"  Target: √2/(2π) = {CABIBBO_PRED:.6f}")
print(f"  Experiment: |V_us| = {CKM_EXP['V_us']:.4f}\n")

# Compute for each config
cabibbo_diffs = []
for idx, j in enumerate(pool9):
    d_cfg = j['down']
    m_d = mass_from_coords(*d_cfg['gen1'])
    m_s = mass_from_coords(*d_cfg['gen2'])
    sin_12 = np.sqrt(m_d / m_s)

    diff_pred = abs(sin_12 - CABIBBO_PRED)
    diff_exp = abs(sin_12 - CKM_EXP['V_us'])

    actual_d = np.allclose(d_cfg['gen1'], ACTUAL['down']['gen1'])
    mark = "★" if actual_d else " "

    cabibbo_diffs.append((idx, sin_12, diff_pred, diff_exp, actual_d))

    print(f"  {idx+1:>2}{mark} sin θ₁₂ = {sin_12:.6f}  "
          f"|Δ from √2/(2π)| = {diff_pred:.6f}  "
          f"|Δ from exp| = {diff_exp:.6f}")

# Rank by closeness to Cabibbo prediction
cabibbo_diffs.sort(key=lambda x: x[2])
print(f"\n  Ranked by |sin θ₁₂ - √2/(2π)|:")
for rank, (idx, sin12, dp, de, is_actual) in enumerate(cabibbo_diffs, 1):
    mark = "★" if is_actual else " "
    print(f"    #{rank}: config {idx+1}{mark}  sin θ₁₂ = {sin12:.6f}  "
          f"|Δ_pred| = {dp:.6f}  |Δ_exp| = {de:.6f}")

print()


# ════════════════════════════════════════════════════════════════
# §8  LATTICE DISTANCE — Gen-1 Separation u ↔ d
# ════════════════════════════════════════════════════════════════

print("\n§8. GEN-1 LATTICE DISPLACEMENT — Δ₁ = d₁ - u₁")
print("-" * 72)

# Up gen-1 is now fixed at (-1.5, -6.0, -1.0) for all 9
up_g1 = pool9[0]['up']['gen1']
print(f"\n  Up gen-1 (fixed): ({up_g1[0]:+.1f}, {up_g1[1]:+.1f}, {up_g1[2]:+.1f})")
print()

print(f"  {'#':>3}  {'Δ₁ = d₁-u₁':^24}  {'|Δ|_latt':>8}  {'Δ·basis':>8}  "
      f"{'p-r family':>10}  {'Δ sum':>6}")
print("  " + "─" * 85)

for idx, j in enumerate(pool9):
    d1 = j['down']['gen1']
    delta = d1 - up_g1
    d_latt = np.linalg.norm(delta)
    d_mass = np.dot(delta, BASIS)
    pr_fam = d1[0] - d1[2]

    actual_d = np.allclose(d1, ACTUAL['down']['gen1'])
    mark = "★" if actual_d else " "

    print(f"  {idx+1:>2}{mark} ({delta[0]:+5.1f},{delta[1]:+5.1f},{delta[2]:+5.1f})  "
          f"{d_latt:8.3f}  {d_mass:+8.4f}  {pr_fam:+10.1f}  {sum(delta):+6.1f}")

print()


# ════════════════════════════════════════════════════════════════
# §9  COMBINED CONSTRAINTS — Intersection of Conditions
# ════════════════════════════════════════════════════════════════

print("\n§9. COMBINED CONSTRAINTS — Hierarchy of Conditions")
print("-" * 72)

# Test multiple constraints and how they intersect
constraints = {}

for idx, j in enumerate(pool9):
    d_cfg = j['down']
    d1 = d_cfg['gen1']
    d3 = d_cfg['gen3']
    B_d = (d3 - d1) / 2.0 - 2 * A_DOWN
    m_d = mass_from_coords(*d1)
    m_s = mass_from_coords(*d_cfg['gen2'])
    sin_12 = np.sqrt(m_d / m_s)
    is_actual = np.allclose(d1, ACTUAL['down']['gen1'])

    # C1: B is all-integer
    c1 = coords_integer(B_d)
    # C2: Cabibbo within 5%
    c2 = abs(sin_12 / CABIBBO_PRED - 1) < 0.05
    # C3: Cabibbo within 10%
    c3 = abs(sin_12 / CABIBBO_PRED - 1) < 0.10
    # C4: d gen-1 in p-r family A (p-r = +1.5)
    c4 = abs((d1[0] - d1[2]) - 1.5) < 0.01
    # C5: B components have same sign
    c5 = all(b > 0 for b in B_d) or all(b < 0 for b in B_d)
    # C6: sum(d1) = -10.5
    c6 = abs(sum(d1) - (-10.5)) < 0.01
    # C7: |B|² is integer
    B_sq = np.dot(B_d, B_d)
    c7 = is_integer(B_sq)
    # C8: B · A_down is special (integer or zero)
    BA = np.dot(B_d, A_DOWN)
    c8 = is_integer(BA)

    constraints[idx] = {
        'c1_B_int': c1, 'c2_cab_5pct': c2, 'c3_cab_10pct': c3,
        'c4_family_A': c4, 'c5_B_same_sign': c5, 'c6_sum_d1': c6,
        'c7_B_sq_int': c7, 'c8_BA_int': c8,
        'is_actual': is_actual,
        'sin12': sin_12, 'B': B_d, 'B_sq': B_sq, 'BA': BA
    }

# Print constraint matrix
print(f"\n  {'#':>3}  {'B_int':>5}  {'Cab5%':>5}  {'Cab10%':>6}  {'FamA':>5}  "
      f"{'B>0':>4}  {'Σ=-10.5':>7}  {'B²∈Z':>5}  {'B·A∈Z':>6}  "
      f"{'sin θ₁₂':>8}  {'|B|²':>7}  {'B·A':>7}")
print("  " + "─" * 110)

for idx in range(N9):
    c = constraints[idx]
    mark = "★" if c['is_actual'] else " "
    print(f"  {idx+1:>2}{mark} "
          f"{'✓' if c['c1_B_int'] else '·':>5}  "
          f"{'✓' if c['c2_cab_5pct'] else '·':>5}  "
          f"{'✓' if c['c3_cab_10pct'] else '·':>6}  "
          f"{'✓' if c['c4_family_A'] else '·':>5}  "
          f"{'✓' if c['c5_B_same_sign'] else '·':>4}  "
          f"{'✓' if c['c6_sum_d1'] else '·':>7}  "
          f"{'✓' if c['c7_B_sq_int'] else '·':>5}  "
          f"{'✓' if c['c8_BA_int'] else '·':>6}  "
          f"{c['sin12']:8.5f}  {c['B_sq']:7.2f}  {c['BA']:+7.3f}")

# Count intersections
print(f"\n  Constraint                      Pass   Intersection with previous")
print("  " + "─" * 60)

labels = [
    ('q₂=q₃ (up)', lambda _: True),  # already applied
    ('B_d all-integer', lambda c: c['c1_B_int']),
    ('Cabibbo 10%', lambda c: c['c3_cab_10pct']),
    ('Cabibbo 5%', lambda c: c['c2_cab_5pct']),
    ('Family A (p-r=+1.5)', lambda c: c['c4_family_A']),
    ('B > 0 (all positive)', lambda c: c['c5_B_same_sign']),
    ('B·A ∈ Z', lambda c: c['c8_BA_int']),
]

# Cumulative intersection
current_set = set(range(N9))
for lbl, test in labels:
    pass_set = set(i for i in range(N9) if test(constraints[i]))
    current_set = current_set & pass_set
    actual_in = any(constraints[i]['is_actual'] for i in current_set)
    print(f"  {lbl:<30}  {len(pass_set):>3}/{N9}  "
          f"cumulative: {len(current_set):>2}  actual: {'✓' if actual_in else '✗'}")

print()


# ════════════════════════════════════════════════════════════════
# §10  DOWN FAMILY STRUCTURE — Two Families Revisited
# ════════════════════════════════════════════════════════════════

print("\n§10. DOWN FAMILY STRUCTURE — p-r Discriminator Revisited")
print("-" * 72)

# We know from survivor_analysis: 2 families by p-r
# Family A: p-r = +1.5 (5 members), Family B: p-r = +7.5 (4 members)
# Within each, coords advance as (+1,+3,+1)

fam_A = []
fam_B = []
for idx, j in enumerate(pool9):
    d1 = j['down']['gen1']
    pr = d1[0] - d1[2]
    if abs(pr - 1.5) < 0.1:
        fam_A.append(idx)
    elif abs(pr - 7.5) < 0.1:
        fam_B.append(idx)
    else:
        print(f"  WARNING: config {idx+1} has p-r = {pr:.1f} (neither family!)")

print(f"\n  Family A (p-r = +1.5): {len(fam_A)} configs — indices {[i+1 for i in fam_A]}")
print(f"  Family B (p-r = +7.5): {len(fam_B)} configs — indices {[i+1 for i in fam_B]}")
print()

for fam_label, fam_indices in [("A (p-r=+1.5)", fam_A), ("B (p-r=+7.5)", fam_B)]:
    print(f"  Family {fam_label}:")
    for idx in fam_indices:
        c = constraints[idx]
        d1 = pool9[idx]['down']['gen1']
        mark = "★" if c['is_actual'] else " "
        B = c['B']
        print(f"    {idx+1:>2}{mark} d₁=({d1[0]:+5.1f},{d1[1]:+5.1f},{d1[2]:+5.1f})  "
              f"B=({B[0]:+5.2f},{B[1]:+5.2f},{B[2]:+5.2f})  "
              f"B_int={'Y' if c['c1_B_int'] else 'N'}  "
              f"sin θ₁₂={c['sin12']:.5f}  Cab5%={'Y' if c['c2_cab_5pct'] else 'N'}")
    print()


# ════════════════════════════════════════════════════════════════
# §11  INTER-SECTOR SYMMETRY — u↔d at each generation
# ════════════════════════════════════════════════════════════════

print("\n§11. INTER-SECTOR SYMMETRY — u↔d Coordinate Differences per Generation")
print("-" * 72)

up_cfg = pool9[0]['up']
lep_cfg = pool9[0]['lepton']

for idx, j in enumerate(pool9):
    d_cfg = j['down']
    actual_d = np.allclose(d_cfg['gen1'], ACTUAL['down']['gen1'])
    mark = "★" if actual_d else " "

    print(f"\n  Config {idx+1}{mark}:")
    for g, gk in enumerate(['gen1', 'gen2', 'gen3'], 1):
        du = d_cfg[gk] - up_cfg[gk]
        dl = d_cfg[gk] - lep_cfg[gk]
        ul = up_cfg[gk] - lep_cfg[gk]
        print(f"    g={g}: d-u=({du[0]:+5.1f},{du[1]:+5.1f},{du[2]:+5.1f})  "
              f"d-ℓ=({dl[0]:+5.1f},{dl[1]:+5.1f},{dl[2]:+5.1f})  "
              f"u-ℓ=({ul[0]:+5.1f},{ul[1]:+5.1f},{ul[2]:+5.1f})")

print()


# ════════════════════════════════════════════════════════════════
# §12  COORDINATE INTEGRALITY — Type pattern of gen-2
# ════════════════════════════════════════════════════════════════

print("\n§12. COORDINATE INTEGRALITY — Type Pattern at gen-2 (down)")
print("-" * 72)

print(f"\n  For each config, classify each component as Integer (Z) or Half-int (H):\n")

for idx, j in enumerate(pool9):
    d2 = j['down']['gen2']
    types = []
    for v in d2:
        if is_integer(v):
            types.append('Z')
        elif is_half_int(v):
            types.append('H')
        else:
            types.append('?')

    actual_d = np.allclose(j['down']['gen1'], ACTUAL['down']['gen1'])
    mark = "★" if actual_d else " "

    print(f"  {idx+1:>2}{mark} gen-2=({d2[0]:+5.1f},{d2[1]:+5.1f},{d2[2]:+5.1f})  "
          f"type=({','.join(types)})")

# Also check: is gen-2 all-integer? (special?)
n_all_int_g2 = sum(1 for j in pool9
                   if coords_integer(j['down']['gen2']))
n_all_hint = sum(1 for j in pool9
                 if all(is_half_int(v) and not is_integer(v) for v in j['down']['gen2']))
print(f"\n  Down gen-2 all-integer: {n_all_int_g2}/{N9}")
print()


# ════════════════════════════════════════════════════════════════
# §13  THE ACTUAL CONFIG — why is it special?
# ════════════════════════════════════════════════════════════════

print("\n§13. THE ACTUAL CONFIG — What Makes It Special?")
print("-" * 72)

actual_idx = None
for idx, j in enumerate(pool9):
    if np.allclose(j['down']['gen1'], ACTUAL['down']['gen1']):
        actual_idx = idx
        break

if actual_idx is not None:
    c = constraints[actual_idx]
    j = pool9[actual_idx]

    print(f"\n  Actual config is #{actual_idx+1} of {N9}")
    print(f"  Down gen-1: ({j['down']['gen1'][0]:+.1f}, {j['down']['gen1'][1]:+.1f}, {j['down']['gen1'][2]:+.1f})")
    print(f"  Down gen-2: ({j['down']['gen2'][0]:+.1f}, {j['down']['gen2'][1]:+.1f}, {j['down']['gen2'][2]:+.1f})")
    print(f"  Down gen-3: ({j['down']['gen3'][0]:+.1f}, {j['down']['gen3'][1]:+.1f}, {j['down']['gen3'][2]:+.1f})")
    print(f"  Down B:     ({c['B'][0]:+.2f}, {c['B'][1]:+.2f}, {c['B'][2]:+.2f})")
    print()
    print(f"  Properties that MIGHT distinguish it:")
    print(f"    B all-integer:        {c['c1_B_int']}  (unique? {n_int_B == 1 if n_int_B is not None else '?'})")
    print(f"    Cabibbo within 5%:    {c['c2_cab_5pct']}  sin θ₁₂ = {c['sin12']:.6f}")
    print(f"    Cabibbo within 10%:   {c['c3_cab_10pct']}")
    print(f"    Family A (p-r=+1.5):  {c['c4_family_A']}")
    print(f"    B all positive:       {c['c5_B_same_sign']}")
    print(f"    B·A ∈ Z:              {c['c8_BA_int']}  (B·A = {c['BA']:+.3f})")
    print(f"    |B|² ∈ Z:             {c['c7_B_sq_int']}  (|B|² = {c['B_sq']:.2f})")

    # Find the UNIQUE combination that selects only the actual
    print(f"\n  Searching for MINIMAL constraint set that selects ONLY actual...")

    constraint_list = [
        ('B_int', lambda i: constraints[i]['c1_B_int']),
        ('Cab10%', lambda i: constraints[i]['c3_cab_10pct']),
        ('Cab5%', lambda i: constraints[i]['c2_cab_5pct']),
        ('FamA', lambda i: constraints[i]['c4_family_A']),
        ('B>0', lambda i: constraints[i]['c5_B_same_sign']),
        ('B·A∈Z', lambda i: constraints[i]['c8_BA_int']),
        ('B²∈Z', lambda i: constraints[i]['c7_B_sq_int']),
    ]

    # Try all pairs
    from itertools import combinations
    print()
    for k in range(1, len(constraint_list) + 1):
        for combo in combinations(constraint_list, k):
            survivors = [i for i in range(N9) if all(test(i) for _, test in combo)]
            if len(survivors) == 1 and survivors[0] == actual_idx:
                names = " ∧ ".join(n for n, _ in combo)
                print(f"    SELECT with {k} constraint{'s' if k>1 else ''}: {names}")
        # Stop after finding minimal
        if any(
            len([i for i in range(N9) if all(test(i) for _, test in combo)]) == 1
            and [i for i in range(N9) if all(test(i) for _, test in combo)][0] == actual_idx
            for combo in combinations(constraint_list, k)
        ):
            break

print()


# ════════════════════════════════════════════════════════════════
# §14  CROSS-CHECK: DOES B_int + q₂=q₃ SELECT FROM 72?
# ════════════════════════════════════════════════════════════════

print("\n§14. CROSS-CHECK — q₂=q₃(up) + B_int(down) from full 72")
print("-" * 72)

# Apply both constraints to the full 72
final_survivors = []
for j in s2:
    # Constraint 1: q₂ = q₃ in up
    if abs(j['up']['gen2'][1] - j['up']['gen3'][1]) > 1e-10:
        continue
    # Constraint 2: down B is all-integer
    B_d = (j['down']['gen3'] - j['down']['gen1']) / 2.0 - 2 * A_DOWN
    if not coords_integer(B_d):
        continue
    final_survivors.append(j)

print(f"\n  72 → q₂=q₃(up) → 9 → B_int(down) → {len(final_survivors)}")

if len(final_survivors) == 1:
    fs = final_survivors[0]
    print(f"\n  ✓ UNIQUE CONFIG SELECTED!")
    for sector in ['lepton', 'up', 'down']:
        print(f"\n  {sector.upper()}:")
        for gk, gl in [('gen1','gen-1'), ('gen2','gen-2'), ('gen3','gen-3')]:
            c = fs[sector][gk]
            m = mass_from_coords(*c)
            print(f"    {gl}: ({c[0]:+5.1f},{c[1]:+5.1f},{c[2]:+5.1f})  m = {m:.4f} MeV")
    print()
    # Verify it's the actual
    matches_actual = (
        np.allclose(fs['up']['gen1'], ACTUAL['up']['gen1']) and
        np.allclose(fs['down']['gen1'], ACTUAL['down']['gen1'])
    )
    print(f"  Matches physical reality: {'✓ YES' if matches_actual else '✗ NO'}")
else:
    print(f"\n  {len(final_survivors)} configs survive (not unique)")
    for i, fs in enumerate(final_survivors):
        d1 = fs['down']['gen1']
        u1 = fs['up']['gen1']
        actual = (np.allclose(u1, ACTUAL['up']['gen1']) and
                  np.allclose(d1, ACTUAL['down']['gen1']))
        mark = "★" if actual else " "
        print(f"    {i+1}{mark} u₁=({u1[0]:+.1f},{u1[1]:+.1f},{u1[2]:+.1f})  "
              f"d₁=({d1[0]:+.1f},{d1[1]:+.1f},{d1[2]:+.1f})")

print()


# ════════════════════════════════════════════════════════════════
# §15  SUMMARY
# ════════════════════════════════════════════════════════════════

print("=" * 72)
print("  §15. SUMMARY")
print("=" * 72)

print(f"""
  DEGENERACY BREAKING CHAIN:
  ──────────────────────────
  2500 configs (gen-3 gaps, 10%)
    │
    ├── 2% tolerance → 72
    │
    ├── q₂ = q₃ in up sector → 9   (up fully fixed)
    │
    ├── B_d all-integer → ???       (see above)
    │
    └── ... → 1 (= physical reality)

  KEY FINDINGS:
    1. q₂ = q₃ in the up sector fixes ALL up-sector coords to the actual values
    2. The remaining 9 configs differ ONLY in down gen-1 position
    3. Down gen-1 forms two families: A (p-r=+1.5, 5 members) and B (p-r=+7.5, 4 members)
    4. The actual config is in Family A

  The question is: WHAT additional constraint selects the actual
  down gen-1 = (-0.5, -8.0, -2.0) from the remaining 9?
""")


# ════════════════════════════════════════════════════════════════
# CLEANUP
# ════════════════════════════════════════════════════════════════

sys.stdout = _orig_stdout
_outfile.close()
print(f"Output written to {OUTPUT_PATH}")
