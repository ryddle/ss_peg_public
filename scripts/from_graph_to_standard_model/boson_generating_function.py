#!/usr/bin/env python3
"""
boson_generating_function.py — Boson Sector in the Generating Function Framework
==================================================================================

The generating function F(2I₃, Nc, g) reproduces all 36 fermion lattice
coordinates and predicts neutrino Dirac masses (§14 of 00_DERIVATION.md).
Prior work (boson_extension_output.md) showed that bosons do NOT sit on
any fermion parabola, nor can any (2I₃, Nc) pair reproduce boson coordinates.

This script asks the NEXT questions:
  - Do massive bosons (W, Z, H) form their OWN parabola?
  - Which ordering satisfies the universal flatness condition b = -5a?
  - How do the bosonic parabola coefficients compare to the fermion manifold?
  - What is the kinetic action of the boson sector?
  - Can we extend the GF with spin S as a third quantum number?
  - What structural relationships link bosons to fermions?

Structure:
  §1. Boson lattice coordinates — verification
  §2. Ordering analysis — all 6 permutations, flatness test
  §3. The flat ordering: parabola coefficients and properties
  §4. Fermion manifold projection — effective quantum numbers
  §5. Spin-extended generating function
  §6. Kinetic action — all 5 sectors compared
  §7. Boson-fermion structural analysis
  §8. Building-block formulas for boson mass ratios
  §9. Summary

Author: SS-PEG program
Date: 2026
"""

import numpy as np
from fractions import Fraction
from itertools import permutations, product as cart_product
import sys, os

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "boson_generating_function_output.md")
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
# BUILDING BLOCKS
# ════════════════════════════════════════════════════════════════

R2   = Fraction(1, 2)
R3   = (np.sqrt(57) - 3) / (2 * np.sqrt(17))  # 0.55228...
R_EE = 1 / np.sqrt(2)
hv3  = 3.0
hv2  = 2.0

LN2  = np.log(2)
LNR3 = np.log(R3)
LN3  = np.log(3)

BLOCKS = np.array([0.5, R3, R_EE, hv3, hv2])
BNAMES = ["R₂", "R₃", "R_EE", "h∨₃", "h∨₂"]
LOG_BLOCKS = np.log(BLOCKS)

m_e = 0.51099895  # MeV

# PDG 2024 masses (MeV)
m_W = 80379.0
m_Z = 91187.6
m_H = 125100.0
m_t = 172500.0

BOSON_MASSES = {"W": m_W, "Z": m_Z, "H": m_H}

# Known boson 3D coordinates (from bosonic_significance.py spanning tree)
BOSON_COORDS = {
    "W": (Fraction(11, 2), Fraction(-10), Fraction(2)),
    "Z": (Fraction(8),     Fraction(-11), Fraction(0)),
    "H": (Fraction(7),     Fraction(-9),  Fraction(2)),
}


def mass_from_coords(p, q, r):
    """Compute mass in MeV from lattice coordinates."""
    ln_ratio = float(p) * LN2 + float(q) * LNR3 + float(r) * LN3
    return m_e * np.exp(ln_ratio)


def formula_str(exps):
    parts = []
    for name, e in zip(BNAMES, exps):
        if e == 0: continue
        elif e == 1: parts.append(name)
        elif e == -1: parts.append(f"{name}⁻¹")
        else:
            sup = str(e).replace("-", "⁻")
            parts.append(f"{name}^{sup}")
    return " · ".join(parts) if parts else "1"


def make_grid(exp_range=(-6, 6)):
    exps_list = list(range(exp_range[0], exp_range[1] + 1))
    grid = np.array(list(cart_product(exps_list, repeat=5)))
    log_vals = grid @ LOG_BLOCKS
    return grid, log_vals

GRID, LOG_GRID = make_grid((-6, 6))


def best_formula(target_val, n=5):
    log_target = np.log(abs(target_val))
    errors = np.abs(LOG_GRID - log_target)
    top_idx = np.argsort(errors)[:n]
    results = []
    for idx in top_idx:
        exps = tuple(int(x) for x in GRID[idx])
        val = np.exp(LOG_GRID[idx])
        err = abs(val - abs(target_val)) / abs(target_val) * 100
        compl = sum(abs(e) for e in exps)
        results.append((exps, val, err, compl))
    return results


# ════════════════════════════════════════════════════════════════
# GENERATING FUNCTION (from §10.5 of 00_DERIVATION.md)
# ════════════════════════════════════════════════════════════════

GF_COEFFS = {
    "a_p": (Fraction(11, 8),  Fraction(-1),    Fraction(27, 8)),
    "a_q": (Fraction(1, 4),   Fraction(-1, 4), Fraction(1)),
    "a_r": (Fraction(-5, 4),  Fraction(5, 4),  Fraction(-4)),
    "b_p": (Fraction(-33, 8), Fraction(17, 4), Fraction(-95, 8)),
    "b_q": (Fraction(-7, 4),  Fraction(13, 4), Fraction(-21, 2)),
    "b_r": (Fraction(19, 4),  Fraction(-17, 4),Fraction(33, 2)),
    "c_p": (Fraction(9, 4),   Fraction(-7, 2), Fraction(33, 4)),
    "c_q": (Fraction(5, 2),   Fraction(-7),    Fraction(29, 2)),
    "c_r": (Fraction(-3),     Fraction(2),     Fraction(-11)),
}


def eval_gf_coeff(name, I3_2, Nc):
    alpha, beta, gamma = GF_COEFFS[name]
    return alpha * I3_2 + beta * Nc + gamma


def get_sector_coeffs(I3_2, Nc):
    result = {}
    for coord in ["p", "q", "r"]:
        a = eval_gf_coeff(f"a_{coord}", I3_2, Nc)
        b = eval_gf_coeff(f"b_{coord}", I3_2, Nc)
        c = eval_gf_coeff(f"c_{coord}", I3_2, Nc)
        result[coord] = (a, b, c)
    return result


# ════════════════════════════════════════════════════════════════
# FERMION SECTOR COEFFICIENTS (for comparison)
# ════════════════════════════════════════════════════════════════

FERMION_SECTORS = {
    "lepton":   {"2I3": Fraction(-1), "Nc": Fraction(1)},
    "up":       {"2I3": Fraction(+1), "Nc": Fraction(3)},
    "down":     {"2I3": Fraction(-1), "Nc": Fraction(3)},
    "neutrino": {"2I3": Fraction(+1), "Nc": Fraction(1)},
}

FERMION_COEFFS = {}
for name, qn in FERMION_SECTORS.items():
    FERMION_COEFFS[name] = get_sector_coeffs(qn["2I3"], qn["Nc"])


# ════════════════════════════════════════════════════════════════
print("=" * 78)
print("  BOSON SECTOR IN THE GENERATING FUNCTION FRAMEWORK")
print("  SS-PEG Research Program — 2026")
print("=" * 78)

# ════════════════════════════════════════════════════════════════
# §1. BOSON LATTICE COORDINATES — VERIFICATION
# ════════════════════════════════════════════════════════════════

print("\n" + "=" * 78)
print("§1. BOSON LATTICE COORDINATES — VERIFICATION")
print("=" * 78)

print(f"\n  {'Boson':>6} {'p':>8} {'q':>8} {'r':>8} │ {'m_pred (MeV)':>14} {'m_exp (MeV)':>12} {'Err%':>8}")
print(f"  {'─' * 68}")

for name in ["W", "Z", "H"]:
    p, q, r = BOSON_COORDS[name]
    m_pred = mass_from_coords(p, q, r)
    m_exp = BOSON_MASSES[name]
    err = abs(m_pred - m_exp) / m_exp * 100
    print(f"  {name:>6} {str(p):>8} {str(q):>8} {str(r):>8} │ {m_pred:>14.1f} {m_exp:>12.1f} {err:>7.3f}%")

print(f"\n  Note: All coordinates are multiples of 1/2 (half-integer lattice).")
print(f"  This is the same quantization as fermion coordinates.")


# ════════════════════════════════════════════════════════════════
# §2. ORDERING ANALYSIS — ALL 6 PERMUTATIONS
# ════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 78)
print("§2. ORDERING ANALYSIS — FLATNESS TEST FOR ALL PERMUTATIONS")
print("=" * 78)

print(f"""
  For fermions, each sector has a FLAT coordinate where b/a = -5
  (universal turning point g* = 5/2). We test all 6 orderings of
  (W, Z, H) as (g=1, g=2, g=3) to find which one(s) reproduce
  this flatness condition.
""")

boson_names = ["W", "Z", "H"]
boson_coords_list = [BOSON_COORDS[b] for b in boson_names]


def fit_parabola_exact(pts):
    """Fit exact rational parabola through 3 points at g=1,2,3."""
    y1, y2, y3 = pts
    # y = a*g² + b*g + c with g = 1, 2, 3
    # y1 = a + b + c
    # y2 = 4a + 2b + c
    # y3 = 9a + 3b + c
    a = (y1 - 2*y2 + y3) / 2
    b = (-5*y1 + 8*y2 - 3*y3) / 2
    c = 3*y1 - 3*y2 + y3
    return a, b, c


ordering_results = []
print(f"  {'Ordering':>15} │ {'Coord':>5} {'a':>8} {'b':>8} {'c':>8} │ {'b/a':>8} {'Flat?':>6}")
print(f"  {'─' * 70}")

for perm in permutations(range(3)):
    ordering = [boson_names[i] for i in perm]
    label = "→".join(ordering)
    pts_p = [BOSON_COORDS[ordering[g]][0] for g in range(3)]
    pts_q = [BOSON_COORDS[ordering[g]][1] for g in range(3)]
    pts_r = [BOSON_COORDS[ordering[g]][2] for g in range(3)]

    flat_count = 0
    flat_coords = []

    for coord_name, pts in [("p", pts_p), ("q", pts_q), ("r", pts_r)]:
        a, b, c = fit_parabola_exact(pts)
        if a != 0:
            ratio = b / a
            is_flat = (ratio == -5)
        else:
            ratio = None
            is_flat = False

        if is_flat:
            flat_count += 1
            flat_coords.append(coord_name)

        ratio_str = f"{float(ratio):>8.3f}" if ratio is not None else " linear"
        flat_str = "✓ FLAT" if is_flat else ""
        print(f"  {label:>15} │ {coord_name:>5} {str(a):>8} {str(b):>8} {str(c):>8} │ {ratio_str} {flat_str:>6}")

    ordering_results.append((label, perm, flat_count, flat_coords))
    print()

print(f"\n  SUMMARY:")
for label, perm, fc, fcoords in ordering_results:
    if fc > 0:
        print(f"    ★ {label}: {fc} flat coordinate(s): {', '.join(fcoords)}")
    else:
        print(f"      {label}: no flatness")


# Find the flat ordering(s)
flat_orderings = [(label, perm, fc, fcoords) for label, perm, fc, fcoords in ordering_results if fc > 0]

if flat_orderings:
    best_ordering = max(flat_orderings, key=lambda x: x[2])
    FLAT_ORDER_LABEL = best_ordering[0]
    FLAT_ORDER_PERM = best_ordering[1]
    FLAT_ORDER_NAMES = [boson_names[i] for i in FLAT_ORDER_PERM]
    print(f"\n  ★ Selected ordering: {FLAT_ORDER_LABEL} (matches fermion flatness condition)")
else:
    # Fallback: use Z→H→W (known from prior analysis)
    FLAT_ORDER_NAMES = ["Z", "H", "W"]
    FLAT_ORDER_LABEL = "Z→H→W"
    print(f"\n  (No flat ordering found; using Z→H→W from prior analysis)")


# ════════════════════════════════════════════════════════════════
# §3. THE FLAT ORDERING: PARABOLA COEFFICIENTS AND PROPERTIES
# ════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 78)
print(f"§3. THE {FLAT_ORDER_LABEL} PARABOLA: COEFFICIENTS AND PROPERTIES")
print("=" * 78)

# Compute exact parabola coefficients for the flat ordering
boson_parabola = {}
for coord_idx, coord_name in enumerate(["p", "q", "r"]):
    pts = [BOSON_COORDS[FLAT_ORDER_NAMES[g]][coord_idx] for g in range(3)]
    a, b, c = fit_parabola_exact(pts)
    boson_parabola[coord_name] = (a, b, c)

print(f"\n  Assignment: g=1 → {FLAT_ORDER_NAMES[0]}, g=2 → {FLAT_ORDER_NAMES[1]}, g=3 → {FLAT_ORDER_NAMES[2]}")
print(f"\n  {'Coord':>5} │ {'a':>8} {'b':>8} {'c':>8} │ {'b/a':>8} {'g*':>8} │ {'Notes'}")
print(f"  {'─' * 72}")

for coord_name in ["p", "q", "r"]:
    a, b, c = boson_parabola[coord_name]
    if a != 0:
        ratio = b / a
        g_star = -b / (2 * a)
        ratio_str = f"{float(ratio):>8.3f}"
        g_star_str = f"{float(g_star):>8.4f}"
    else:
        ratio_str = "   ∞"
        g_star_str = "   ∞"

    is_flat = (a != 0 and b / a == -5)
    notes = "FLAT (b = -5a, g* = 5/2)" if is_flat else ""
    print(f"  {coord_name:>5} │ {str(a):>8} {str(b):>8} {str(c):>8} │ {ratio_str} {g_star_str} │ {notes}")

# Rationality check
print(f"\n  Rationality check:")
all_denoms = set()
for coord_name in ["p", "q", "r"]:
    a, b, c = boson_parabola[coord_name]
    for coeff_name, val in [("a", a), ("b", b), ("c", c)]:
        d = val.denominator
        all_denoms.add(d)
        print(f"    {coord_name}_{coeff_name} = {val} (denominator {d})")

max_denom = max(all_denoms)
print(f"\n  All denominators: {sorted(all_denoms)}")
print(f"  Maximum denominator: {max_denom}")
ferm_str = "SAME as fermion coefficients (multiples of 1/4)" if max_denom <= 4 else "DIFFERENT from fermion coefficients"
print(f"  → {ferm_str}")

# Electron origin condition: c = -a - b?
print(f"\n  Electron origin condition (c = -a - b):")
for coord_name in ["p", "q", "r"]:
    a, b, c = boson_parabola[coord_name]
    expected = -a - b
    holds = (c == expected)
    print(f"    {coord_name}: c = {c}, -a-b = {expected} → {'✓ HOLDS' if holds else '✗ FAILS'}")


# ════════════════════════════════════════════════════════════════
# §4. FERMION MANIFOLD PROJECTION
# ════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 78)
print("§4. FERMION MANIFOLD PROJECTION — EFFECTIVE QUANTUM NUMBERS")
print("=" * 78)

print(f"""
  The fermion GF has X(s) = α·(2I₃) + β·Nc + γ for each of 9 coefficients.
  For the boson parabola, we ask: what (2I₃_eff, Nc_eff) would produce
  the observed coefficients? This is an OVERDETERMINED system (9 equations,
  2 unknowns), so we measure the residual.
""")

# Build the system: for each coefficient, α·x + β·y + γ = observed
# where x = 2I₃_eff, y = Nc_eff
coeff_names = ["a_p", "a_q", "a_r", "b_p", "b_q", "b_r", "c_p", "c_q", "c_r"]
boson_coeff_vals = []
for cname in coeff_names:
    level = cname[0]  # 'a', 'b', or 'c'
    coord = cname[2]  # 'p', 'q', or 'r'
    idx = {"a": 0, "b": 1, "c": 2}[level]
    boson_coeff_vals.append(float(boson_parabola[coord][idx]))

# A matrix: [α_i, β_i] for each coefficient
# b vector: observed - γ_i
A_mat = np.zeros((9, 2))
b_vec = np.zeros(9)
for i, cname in enumerate(coeff_names):
    alpha, beta, gamma = GF_COEFFS[cname]
    A_mat[i, 0] = float(alpha)
    A_mat[i, 1] = float(beta)
    b_vec[i] = boson_coeff_vals[i] - float(gamma)

# Least-squares solution
x_eff, residuals, rank, sv = np.linalg.lstsq(A_mat, b_vec, rcond=None)
I3_2_eff, Nc_eff = x_eff

print(f"  Least-squares effective quantum numbers:")
print(f"    2I₃_eff = {I3_2_eff:.6f}")
print(f"    Nc_eff  = {Nc_eff:.6f}")

# Predicted vs observed
print(f"\n  {'Coefficient':>12} │ {'Observed':>10} {'GF Predicted':>12} {'Residual':>10} {'Rel.Err%':>10}")
print(f"  {'─' * 62}")

predicted = A_mat @ x_eff + np.array([float(GF_COEFFS[cn][2]) for cn in coeff_names])
total_sq_res = 0
for i, cname in enumerate(coeff_names):
    obs = boson_coeff_vals[i]
    pred = predicted[i]
    res = obs - pred
    total_sq_res += res**2
    rel_err = abs(res / obs) * 100 if obs != 0 else float('inf')
    print(f"  {cname:>12} │ {obs:>10.4f} {pred:>12.4f} {res:>10.4f} {rel_err:>9.2f}%")

rms_res = np.sqrt(total_sq_res / 9)
print(f"\n  RMS residual: {rms_res:.6f}")
print(f"  Total sum of squared residuals: {total_sq_res:.6f}")

# Per-equation solution (each pair of equations → different effective QN)
print(f"\n  Per-coefficient-pair effective (2I₃, Nc):")
print(f"  {'Pair':>18} │ {'2I₃_eff':>10} {'Nc_eff':>10} │ {'Notes'}")
print(f"  {'─' * 60}")

for i in range(0, 9, 1):
    for j in range(i+1, 9):
        A_sub = A_mat[[i, j], :]
        b_sub = b_vec[[i, j]]
        det = np.linalg.det(A_sub)
        if abs(det) > 1e-10:
            x_sub = np.linalg.solve(A_sub, b_sub)
            label = f"({coeff_names[i]}, {coeff_names[j]})"
            # Check if this solution is close to any fermion sector
            close_to = ""
            for sname, sqn in FERMION_SECTORS.items():
                d = np.sqrt((x_sub[0] - float(sqn["2I3"]))**2 + (x_sub[1] - float(sqn["Nc"]))**2)
                if d < 2.0:
                    close_to = f"near {sname}"
                    break
            # Only print a few representative ones (first 15)
            if i * 9 + j < 45:
                pass  # Will print below

# Instead, find the range of solutions
all_solutions = []
for i in range(9):
    for j in range(i+1, 9):
        A_sub = A_mat[[i, j], :]
        b_sub = b_vec[[i, j]]
        det = np.linalg.det(A_sub)
        if abs(det) > 1e-10:
            x_sub = np.linalg.solve(A_sub, b_sub)
            all_solutions.append((coeff_names[i], coeff_names[j], x_sub[0], x_sub[1]))

if all_solutions:
    I3_vals = [s[2] for s in all_solutions]
    Nc_vals = [s[3] for s in all_solutions]
    print(f"\n  Range of pairwise solutions ({len(all_solutions)} solvable pairs):")
    print(f"    2I₃_eff ∈ [{min(I3_vals):.4f}, {max(I3_vals):.4f}], mean = {np.mean(I3_vals):.4f}")
    print(f"    Nc_eff  ∈ [{min(Nc_vals):.4f}, {max(Nc_vals):.4f}], mean = {np.mean(Nc_vals):.4f}")
    print(f"    Spread(2I₃) = {max(I3_vals) - min(I3_vals):.4f}")
    print(f"    Spread(Nc)  = {max(Nc_vals) - min(Nc_vals):.4f}")

    # Find the pair with minimum residual over ALL 9 equations
    best_pair_res = float('inf')
    best_pair = None
    for c1, c2, i3, nc in all_solutions:
        pred_all = np.array([float(GF_COEFFS[cn][0]) * i3 + float(GF_COEFFS[cn][1]) * nc + float(GF_COEFFS[cn][2]) for cn in coeff_names])
        res_all = np.sum((np.array(boson_coeff_vals) - pred_all)**2)
        if res_all < best_pair_res:
            best_pair_res = res_all
            best_pair = (c1, c2, i3, nc)

    if best_pair:
        print(f"\n    Best pairwise fit: ({best_pair[0]}, {best_pair[1]})")
        print(f"    → 2I₃_eff = {best_pair[2]:.6f}, Nc_eff = {best_pair[3]:.6f}")
        print(f"    → Total squared residual = {best_pair_res:.6f}")

print(f"\n  CONCLUSION: The large spread confirms bosons are NOT on the fermion")
print(f"  linear manifold. No single (2I₃, Nc) reproduces all 9 boson coefficients.")


# ════════════════════════════════════════════════════════════════
# §5. SPIN-EXTENDED GENERATING FUNCTION
# ════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 78)
print("§5. SPIN-EXTENDED GENERATING FUNCTION")
print("=" * 78)

print(f"""
  All fermion sectors have spin S = 1/2 (2S = 1). The boson sector
  has spin S = 1 for W,Z and S = 0 for H, but they form a SINGLE
  parabola. Can we extend: X = α·(2I₃) + β·Nc + δ·(2S) + γ?

  Problem: All 4 fermion sectors have 2S = 1, so δ is degenerate
  with γ in the fermion fit. We need the boson sector to break the
  degeneracy.

  Approach: Assign effective 2S_eff to the boson sector and solve
  the 5-sector system. But the boson sector must be assigned
  effective (2I₃_B, Nc_B) as well.
""")

# Try: assign bosons their "natural" quantum numbers
# The boson parabola treats W,Z,H as a single sector.
# What are the "average" quantum numbers?
# W: I=1, I₃=+1, Nc=1, 2S=2
# Z: I₃=0, Nc=1, 2S=2
# H: I₃=-1/2, Nc=1, 2S=0
# These are not uniform — fundamental problem.

# Alternative: what if we use I_total (isospin of the representation)?
# W: adjoint of SU(2), I=1 → 2I=2
# Z: mixture, effective I=0 → 2I=0
# H: doublet, I=1/2 → 2I=1 (but I₃=-1/2)

# The cleanest approach: solve for (2I₃_eff, Nc_eff, 2S_eff) from
# the overdetermined system, now with 3 unknowns instead of 2.

A_mat3 = np.zeros((9, 3))
b_vec3 = np.zeros(9)
for i, cname in enumerate(coeff_names):
    alpha, beta, gamma = GF_COEFFS[cname]
    A_mat3[i, 0] = float(alpha)   # 2I₃ coefficient
    A_mat3[i, 1] = float(beta)    # Nc coefficient
    # For fermions, δ·(2S) was absorbed into γ since 2S=1 for all.
    # So the fermion γ = γ_true + δ·1. We don't know δ separately.
    # We need to re-fit the fermion GF with spin included.
    b_vec3[i] = boson_coeff_vals[i] - float(gamma)
    # The 3rd column is the spin contribution, but we need the true γ.
    # Since γ_observed = γ_true + δ, and fermion 2S=1:
    # For bosons: X = α·(2I₃) + β·Nc + δ·(2S_B) + γ_true
    #           = α·(2I₃) + β·Nc + δ·(2S_B) + (γ_obs - δ)
    # So: X_boson = α·(2I₃_B) + β·(Nc_B) + γ_obs + δ·(2S_B - 1)
    A_mat3[i, 2] = 1.0  # coefficient for δ: multiply by (2S_B - 1)
    # But we don't know 2S_B! It could be 0, 1, or 2.

# Try different values of 2S_B
print(f"  Testing effective spin values for the boson sector:")
print(f"  {'2S_B':>6} │ {'2I₃_eff':>10} {'Nc_eff':>10} {'δ':>10} │ {'RMS res':>10}")
print(f"  {'─' * 58}")

for S2_B_test in [0, 1, 2, 3, 4]:
    # X_boson = α·x + β·y + γ_obs + δ·(2S_B - 1)
    # Rearrange: α·x + β·y + δ·(2S_B - 1) = X_boson - γ_obs
    # So: A_mat3[:, 2] = (2S_B - 1) for all rows
    A3 = A_mat3.copy()
    A3[:, 2] = float(S2_B_test - 1)
    b3 = b_vec3.copy()  # = X_boson - γ_obs

    x_eff3, _, _, _ = np.linalg.lstsq(A3, b3, rcond=None)
    pred3 = A3 @ x_eff3 + np.array([float(GF_COEFFS[cn][2]) for cn in coeff_names])
    # Wait, that's wrong. Let me redo:
    # X_boson = α·(2I₃_B) + β·(Nc_B) + δ·(2S_B - 1) + γ_obs
    # So: X_boson - γ_obs = α·(2I₃_B) + β·(Nc_B) + δ·(2S_B - 1)
    # Which is: b_vec3 = A3 @ [2I₃_B, Nc_B, δ]
    pred3_vals = []
    for idx, cn in enumerate(coeff_names):
        alpha, beta, gamma = GF_COEFFS[cn]
        pred = float(alpha) * x_eff3[0] + float(beta) * x_eff3[1] + x_eff3[2] * (S2_B_test - 1) + float(gamma)
        pred3_vals.append(pred)
    res3 = np.array(boson_coeff_vals) - np.array(pred3_vals)
    rms3 = np.sqrt(np.mean(res3**2))
    print(f"  {S2_B_test:>6} │ {x_eff3[0]:>10.4f} {x_eff3[1]:>10.4f} {x_eff3[2]:>10.4f} │ {rms3:>10.6f}")

print(f"\n  Note: The RMS residual measures how well the extended GF fits.")
print(f"  A significant reduction from the 2-parameter fit (§4) would indicate")
print(f"  that spin is a relevant quantum number for the generating function.")
print(f"  2-parameter RMS: {rms_res:.6f}")


# ════════════════════════════════════════════════════════════════
# §6. KINETIC ACTION — ALL SECTORS COMPARED
# ════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 78)
print("§6. KINETIC ACTION — ALL SECTORS COMPARED")
print("=" * 78)

print(f"""
  The kinetic action S_K = ∫₁³ ½·|ẋ|² dg measures the "cost" of
  traversing a sector's parabola from g=1 to g=3.
  S_K^(k) = (52/3)·a² + 8·a·b + b² per coordinate.
  (Using ∫₁³ g² dg = 26/3, ∫₁³ g dg = 4, ∫₁³ dg = 2)
""")

def kinetic_action_coord(a, b):
    """Kinetic action for one coordinate: ∫₁³ ½(2ag+b)² dg."""
    a, b = Fraction(a), Fraction(b)
    return Fraction(52, 3) * a**2 + 8 * a * b + b**2


def kinetic_action_sector(coeffs_dict):
    """Total kinetic action for a sector (sum over p, q, r)."""
    total = Fraction(0)
    per_coord = {}
    for coord in ["p", "q", "r"]:
        a, b, c = coeffs_dict[coord]
        sk = kinetic_action_coord(a, b)
        per_coord[coord] = sk
        total += sk
    return total, per_coord

# Fermion sectors
print(f"  {'Sector':>12} │ {'S_K(p)':>12} {'S_K(q)':>12} {'S_K(r)':>12} │ {'S_K(total)':>14} {'Decimal':>10}")
print(f"  {'─' * 80}")

all_SK = {}
for sname, coeffs in FERMION_COEFFS.items():
    total, per_coord = kinetic_action_sector(coeffs)
    all_SK[sname] = total
    print(f"  {sname:>12} │ {str(per_coord['p']):>12} {str(per_coord['q']):>12} {str(per_coord['r']):>12} │ {str(total):>14} {float(total):>10.4f}")

# Boson sector
total_B, per_coord_B = kinetic_action_sector(boson_parabola)
all_SK["boson"] = total_B
print(f"  {'─' * 80}")
print(f"  {'★ boson':>12} │ {str(per_coord_B['p']):>12} {str(per_coord_B['q']):>12} {str(per_coord_B['r']):>12} │ {str(total_B):>14} {float(total_B):>10.4f}")

# Grand total
print(f"\n  Sector kinetic actions (ordered):")
for sname in sorted(all_SK, key=lambda k: float(all_SK[k])):
    val = all_SK[sname]
    print(f"    {sname:>12}: S_K = {str(val):>12} ≈ {float(val):>10.4f}")

# Check if boson S_K has a nice rational value
print(f"\n  Boson S_K = {total_B} = {float(total_B):.6f}")
print(f"  Denominator: {total_B.denominator}")

# Ratios between sectors
print(f"\n  Kinetic action ratios:")
for sname in ["lepton", "up", "down", "neutrino"]:
    ratio = all_SK[sname] / total_B
    print(f"    S_K({sname}) / S_K(boson) = {ratio} ≈ {float(ratio):.6f}")

# Sum of all fermion sectors
SK_fermion_total = sum(all_SK[s] for s in ["lepton", "up", "down"])
SK_all_charged = SK_fermion_total
print(f"\n  S_K(charged fermions) = {SK_fermion_total} = {float(SK_fermion_total):.4f}")
print(f"  S_K(all 5 sectors) = {SK_fermion_total + all_SK['neutrino'] + total_B} = {float(SK_fermion_total + all_SK['neutrino'] + total_B):.4f}")

# Building-block formula for S_K(boson)
print(f"\n  Building-block formula search for S_K(boson) = {float(total_B):.6f}:")
results = best_formula(float(total_B), n=5)
for exps, val, err, compl in results:
    print(f"    {formula_str(exps):>40} = {val:.6f} (err {err:.4f}%, complexity {compl})")


# ════════════════════════════════════════════════════════════════
# §7. BOSON-FERMION STRUCTURAL ANALYSIS
# ════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 78)
print("§7. BOSON-FERMION STRUCTURAL ANALYSIS")
print("=" * 78)

# 7a. Coefficient comparison table
print(f"\n  7a. Parabola coefficient comparison (all 5 sectors):")
print(f"\n  {'Coeff':>8} │ {'Lepton':>8} {'Up':>8} {'Down':>8} {'Neutrino':>8} │ {'Boson':>8}")
print(f"  {'─' * 60}")

for level_name, level_idx in [("a", 0), ("b", 1), ("c", 2)]:
    for coord in ["p", "q", "r"]:
        label = f"{level_name}_{coord}"
        vals = []
        for sname in ["lepton", "up", "down", "neutrino"]:
            v = FERMION_COEFFS[sname][coord][level_idx]
            vals.append(str(v))
        bval = str(boson_parabola[coord][level_idx])
        print(f"  {label:>8} │ {vals[0]:>8} {vals[1]:>8} {vals[2]:>8} {vals[3]:>8} │ {bval:>8}")
    print()

# 7b. Curvature vector comparison
print(f"\n  7b. Curvature vectors a = (a_p, a_q, a_r):")
print(f"  {'Sector':>12} │ {'a_p':>8} {'a_q':>8} {'a_r':>8} │ {'|a|²':>10}")
print(f"  {'─' * 52}")

for sname in ["lepton", "up", "down", "neutrino", "boson"]:
    if sname == "boson":
        coeffs = boson_parabola
    else:
        coeffs = FERMION_COEFFS[sname]
    ap = coeffs["p"][0]
    aq = coeffs["q"][0]
    ar = coeffs["r"][0]
    norm_sq = ap**2 + aq**2 + ar**2
    print(f"  {sname:>12} │ {str(ap):>8} {str(aq):>8} {str(ar):>8} │ {str(norm_sq):>10}")

# 7c. Displacement from boson g=1 to each fermion g=3
print(f"\n  7c. Displacement vectors: boson(g=1) to fermion(g=3):")

boson_g1 = BOSON_COORDS[FLAT_ORDER_NAMES[0]]  # g=1 particle
print(f"  Boson g=1 ({FLAT_ORDER_NAMES[0]}): ({boson_g1[0]}, {boson_g1[1]}, {boson_g1[2]})")

# Fermion g=3 coordinates
ferm_g3 = {}
for sname in ["lepton", "up", "down", "neutrino"]:
    coeffs = FERMION_COEFFS[sname]
    p = coeffs["p"][0] * 9 + coeffs["p"][1] * 3 + coeffs["p"][2]
    q = coeffs["q"][0] * 9 + coeffs["q"][1] * 3 + coeffs["q"][2]
    r = coeffs["r"][0] * 9 + coeffs["r"][1] * 3 + coeffs["r"][2]
    ferm_g3[sname] = (p, q, r)
    dp = p - boson_g1[0]
    dq = q - boson_g1[1]
    dr = r - boson_g1[2]
    print(f"    → {sname:>10}(g=3) = ({p}, {q}, {r}), Δ = ({dp}, {dq}, {dr})")

# 7d. Is the boson sector the "spin completion" of any fermion sector?
print(f"\n  7d. Fermion sector closest to boson parabola (coefficient distance):")

for sname in ["lepton", "up", "down", "neutrino"]:
    dist_sq = Fraction(0)
    for coord in ["p", "q", "r"]:
        for idx in range(3):
            diff = boson_parabola[coord][idx] - FERMION_COEFFS[sname][coord][idx]
            dist_sq += diff**2
    print(f"    d²(boson, {sname:>10}) = {dist_sq} ≈ {float(dist_sq):.4f}")

# 7e. Flatness comparison
print(f"\n  7e. Flatness summary (all sectors):")
print(f"  {'Sector':>12} │ {'Flat coord':>10} {'a_flat':>8} {'b_flat':>8} {'b/a':>8} │ {'g*':>8}")
print(f"  {'─' * 62}")

flat_map_fermion = {"lepton": "r", "up": "q", "down": "p", "neutrino": None}
for sname in ["lepton", "up", "down", "neutrino"]:
    flat_coord = flat_map_fermion[sname]
    if flat_coord is None:
        # Neutrino has NO flat coordinate
        print(f"  {sname:>12} │ {'NONE':>10} {'—':>8} {'—':>8} {'—':>8} │ {'—':>8}")
        continue
    a_f = FERMION_COEFFS[sname][flat_coord][0]
    b_f = FERMION_COEFFS[sname][flat_coord][1]
    ratio = b_f / a_f if a_f != 0 else None
    g_star = -b_f / (2 * a_f) if a_f != 0 else None
    print(f"  {sname:>12} │ {flat_coord:>10} {str(a_f):>8} {str(b_f):>8} {str(ratio):>8} │ {str(g_star):>8}")

# Boson flat coordinate
for coord in ["p", "q", "r"]:
    a, b, c = boson_parabola[coord]
    if a != 0 and b / a == -5:
        g_star = -b / (2 * a)
        print(f"  {'★ boson':>12} │ {coord:>10} {str(a):>8} {str(b):>8} {str(b/a):>8} │ {str(g_star):>8}")
        break


# ════════════════════════════════════════════════════════════════
# §8. BUILDING-BLOCK FORMULAS FOR BOSON MASS RATIOS
# ════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 78)
print("§8. BUILDING-BLOCK FORMULAS FOR BOSON MASS RATIOS")
print("=" * 78)

print(f"""
  Even though bosons don't sit on the fermion GF, their mass ratios
  (both inter-boson and boson-fermion) may be expressible as products
  of the 5 building blocks R₂, R₃, R_EE, h∨₃, h∨₂.
""")

# Known ratios from bosonic_significance.py
known_ratios = {
    "m_t/m_W": ([0, 3, -1, 2, 0], m_t / m_W),
    "m_H/m_W": ([0, 1, -1, 0, 1], m_H / m_W),
    "m_H/m_Z": ([1, 2, 0, 2, 0], m_H / m_Z),
    "m_H/m_t": ([-1, -2, 0, -2, 0], m_H / m_t),  # Note: H < t
}

print(f"  Known boson mass ratios (from spanning tree):")
print(f"  {'Ratio':>12} │ {'Experimental':>12} {'Formula':>30} {'Predicted':>12} {'Err%':>8}")
print(f"  {'─' * 80}")

for name, (exps, exp_ratio) in known_ratios.items():
    pred_val = np.prod(BLOCKS ** np.array(exps, dtype=float))
    err = abs(pred_val - exp_ratio) / exp_ratio * 100
    fstr = formula_str(exps)
    print(f"  {name:>12} │ {exp_ratio:>12.6f} {fstr:>30} {pred_val:>12.6f} {err:>7.3f}%")

# Additional boson ratios
print(f"\n  Additional boson mass ratios:")
additional = {
    "m_Z/m_W": m_Z / m_W,
    "m_W/m_e": m_W / m_e,
    "m_Z/m_e": m_Z / m_e,
    "m_H/m_e": m_H / m_e,
}

for name, ratio in additional.items():
    print(f"\n  {name} = {ratio:.6f}:")
    results = best_formula(ratio, n=3)
    for exps, val, err, compl in results:
        tag = " ★" if err < 1.0 else ""
        print(f"    {formula_str(exps):>40} = {val:.6f} (err {err:.4f}%, |e|={compl}){tag}")

# Inter-boson ratios as lattice displacements
print(f"\n  Inter-boson coordinate displacements (in the {FLAT_ORDER_LABEL} ordering):")
for g1 in range(3):
    for g2 in range(g1+1, 3):
        n1 = FLAT_ORDER_NAMES[g1]
        n2 = FLAT_ORDER_NAMES[g2]
        dp = BOSON_COORDS[n2][0] - BOSON_COORDS[n1][0]
        dq = BOSON_COORDS[n2][1] - BOSON_COORDS[n1][1]
        dr = BOSON_COORDS[n2][2] - BOSON_COORDS[n1][2]
        mass_ratio = BOSON_MASSES[n2] / BOSON_MASSES[n1]
        # This ratio = 2^dp * R3^dq * 3^dr
        pred_ratio = 2**float(dp) * R3**float(dq) * 3**float(dr)
        print(f"    {n1}→{n2}: Δ(p,q,r) = ({dp}, {dq}, {dr}), m_{n2}/m_{n1} = {mass_ratio:.6f}, pred = {pred_ratio:.6f}")


# ════════════════════════════════════════════════════════════════
# §9. TRANSITION MATRIX AND 2-ADIC STRUCTURE
# ════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 78)
print("§9. BOSON TRANSITION MATRIX")
print("=" * 78)

print(f"""
  The fermion transition matrix N has entries n_sk = 2(x(2) - x(1))
  with Smith normal form diag(1,2,4). Does the boson sector share
  this structure?
""")

# Boson displacement: x(2) - x(1)
boson_delta = {}
for coord in ["p", "q", "r"]:
    a, b, c = boson_parabola[coord]
    x1 = a * 1 + b * 1 + c  # g=1
    x2 = a * 4 + b * 2 + c  # g=2
    delta = x2 - x1
    n = 2 * delta
    boson_delta[coord] = (delta, n)
    print(f"  {coord}: x(1) = {x1}, x(2) = {x2}, Δ₁₂ = {delta}, n = 2Δ = {n}")

print(f"\n  Boson N-row: ({boson_delta['p'][1]}, {boson_delta['q'][1]}, {boson_delta['r'][1]})")

# Check 2-adic properties
print(f"\n  2-adic valuations:")
for coord in ["p", "q", "r"]:
    n_val = boson_delta[coord][1]
    v2 = 0
    temp = n_val
    while temp != 0 and temp.denominator == 1 and int(temp) % 2 == 0:
        v2 += 1
        temp = temp / 2
    print(f"    n_{coord} = {n_val}, v₂ = {v2}")

# Compare with fermion N-rows
print(f"\n  Comparison with fermion N-rows:")
print(f"  {'Sector':>12} │ {'n_p':>6} {'n_q':>6} {'n_r':>6}")
print(f"  {'─' * 36}")

for sname in ["lepton", "up", "down"]:
    coeffs = FERMION_COEFFS[sname]
    n_row = []
    for coord in ["p", "q", "r"]:
        a, b, c = coeffs[coord]
        x1 = a + b + c
        x2 = 4*a + 2*b + c
        n = 2 * (x2 - x1)
        n_row.append(n)
    print(f"  {sname:>12} │ {str(n_row[0]):>6} {str(n_row[1]):>6} {str(n_row[2]):>6}")

print(f"  {'─' * 36}")
print(f"  {'★ boson':>12} │ {str(boson_delta['p'][1]):>6} {str(boson_delta['q'][1]):>6} {str(boson_delta['r'][1]):>6}")


# ════════════════════════════════════════════════════════════════
# §10. BOSON CURVATURE MATRIX AND SPECTRAL ANALYSIS
# ════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 78)
print("§10. BOSON CURVATURE MATRIX AND SPECTRAL ANALYSIS")
print("=" * 78)

# The curvature vector for the boson sector
a_boson = np.array([float(boson_parabola[k][0]) for k in ["p", "q", "r"]])
b_boson = np.array([float(boson_parabola[k][1]) for k in ["p", "q", "r"]])
c_boson = np.array([float(boson_parabola[k][2]) for k in ["p", "q", "r"]])

print(f"\n  Boson curvature vector: a = ({a_boson[0]}, {a_boson[1]}, {a_boson[2]})")
print(f"  Boson velocity vector:  b = ({b_boson[0]}, {b_boson[1]}, {b_boson[2]})")
print(f"  Boson offset vector:    c = ({c_boson[0]}, {c_boson[1]}, {c_boson[2]})")

# Inner products
print(f"\n  Inner products:")
print(f"    a·b = {np.dot(a_boson, b_boson):.6f}")
print(f"    a·c = {np.dot(a_boson, c_boson):.6f}")
print(f"    b·c = {np.dot(b_boson, c_boson):.6f}")
print(f"    |a|² = {np.dot(a_boson, a_boson):.6f}")
print(f"    |b|² = {np.dot(b_boson, b_boson):.6f}")
print(f"    |c|² = {np.dot(c_boson, c_boson):.6f}")

# Cross products with fermion curvature vectors
print(f"\n  Curvature alignment (cos θ between a_boson and a_fermion):")
for sname in ["lepton", "up", "down", "neutrino"]:
    a_f = np.array([float(FERMION_COEFFS[sname][k][0]) for k in ["p", "q", "r"]])
    cos_theta = np.dot(a_boson, a_f) / (np.linalg.norm(a_boson) * np.linalg.norm(a_f))
    theta_deg = np.degrees(np.arccos(np.clip(cos_theta, -1, 1)))
    print(f"    cos θ(boson, {sname:>10}) = {cos_theta:>+.6f}  (θ = {theta_deg:.1f}°)")

# Transfer matrix: b = a · M^T for boson sector (1×3 system)
# If b = M · a for some matrix, then M = b ⊗ a / |a|² (rank-1 approximation)
# More precisely, we need b_k = Σ_j M_kj a_j
# With one sector, M is underdetermined. But we can check if b/a ratio is constant.
print(f"\n  b/a ratios per coordinate:")
for i, coord in enumerate(["p", "q", "r"]):
    if a_boson[i] != 0:
        ratio = b_boson[i] / a_boson[i]
        print(f"    b_{coord}/a_{coord} = {ratio:.6f}")
    else:
        print(f"    b_{coord}/a_{coord} = ∞ (a=0)")


# ════════════════════════════════════════════════════════════════
# §11. SUMMARY
# ════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 78)
print("§11. SUMMARY AND CONCLUSIONS")
print("=" * 78)

print(f"""
  KEY FINDINGS:

  1. FLATNESS CONDITION: Among all 6 orderings of (W, Z, H), exactly""")

flat_count_total = sum(1 for _, _, fc, _ in ordering_results if fc > 0)
print(f"     {flat_count_total} ordering(s) satisfy the universal flatness condition b = -5a.")
if flat_orderings:
    for label, _, fc, fcoords in flat_orderings:
        print(f"     → {label}: flat in {', '.join(fcoords)} (same turning point g* = 5/2)")

print(f"""
  2. RATIONAL ARITHMETIC: All 9 boson parabola coefficients are exact
     rationals with denominators in {sorted(all_denoms)}, matching the
     fermion arithmetic (multiples of 1/4).

  3. FERMION MANIFOLD: The boson coefficients are NOT on the fermion
     linear manifold X = α·(2I₃) + β·Nc + γ.
     → Least-squares effective QN: 2I₃_eff = {I3_2_eff:.4f}, Nc_eff = {Nc_eff:.4f}
     → RMS residual = {rms_res:.6f} (significant deviation)

  4. KINETIC ACTION: S_K(boson) = {total_B} ≈ {float(total_B):.4f}""")

# Ordering
ordered = sorted(all_SK.items(), key=lambda x: float(x[1]))
for i, (sname, val) in enumerate(ordered):
    if sname == "boson":
        print(f"     → Rank {i+1}/{len(all_SK)} ({'lowest' if i == 0 else 'highest' if i == len(all_SK)-1 else 'intermediate'})")
        break

print(f"""
  5. BOSON-FERMION RELATIONSHIP: Bosons form a SEPARATE sector with
     its own parabolic structure. The flat coordinate ({[c for c in ['p','q','r'] if boson_parabola[c][0] != 0 and boson_parabola[c][1]/boson_parabola[c][0] == -5][0] if any(boson_parabola[c][0] != 0 and boson_parabola[c][1]/boson_parabola[c][0] == -5 for c in ['p','q','r']) else 'none'})
     shares the universal g* = 5/2 turning point with all 3 charged-fermion
     sectors (leptons→r, up→q, down→p).

  6. INTER-SECTOR PATTERN: With the boson sector, we now have 5 sectors:
     ℓ(r flat), u(q flat), d(p flat), ν(NO flat), B(r flat).
     The boson sector shares the r-flatness with leptons — both are
     color singlets (Nc = 1).

  7. MASS RATIOS: All boson-fermion mass ratios (m_t/m_W, m_H/m_W,
     m_H/m_Z, m_H/m_t) are exact building-block formulas at sub-1%
     accuracy, confirming that bosons participate in the same algebraic
     structure as fermions even though they sit on a separate parabola.

  PHYSICAL INTERPRETATION:
  The generating function F(2I₃, Nc, g) governs FERMION parabolas.
  Bosons require a separate treatment — they are not "generations" of
  the same particle but distinct particles with different quantum numbers.
  Nevertheless, the universal features (rational coefficients, flatness
  at g* = 5/2, building-block mass ratios) persist, suggesting a deeper
  unified structure that encompasses both statistics.
""")

# Close output file
sys.stdout = _orig_stdout
_outfile.close()
print(f"\nOutput saved to: {OUTPUT_PATH}")
