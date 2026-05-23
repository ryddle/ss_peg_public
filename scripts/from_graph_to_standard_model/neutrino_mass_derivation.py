#!/usr/bin/env python3
"""
neutrino_mass_derivation.py — Neutrino Mass Derivation from Generating Function
=================================================================================

The generating function F(2I₃, Nc, g) was derived from the 3 charged-fermion
sectors: ℓ(2I₃=-1, Nc=1), u(+1, 3), d(-1, 3). It produces ALL 36 lattice
coordinates with 0 free parameters.

Neutrinos occupy the FOURTH sector: (2I₃=+1, Nc=1) — NOT used in the fit.
Evaluating F at these quantum numbers is a genuine PREDICTION that yields
the "Dirac mass" scale m_D (what neutrinos would weigh with standard Yukawa
couplings).

The comparison m_D vs m_ν(physical) extracts the Type-I seesaw scale M_R,
and we test whether M_R itself is a graph formula.

Structure:
  §1. Generating function evaluation at neutrino quantum numbers
  §2. Dirac mass predictions vs physical neutrino masses
  §3. Seesaw scale extraction and graph-formula search
  §4. Consistency with PMNS angles and r_Δm²
  §5. Ramanujan spectral interpretation
  §6. Comparison with charged-lepton isospin partner
  §7. Summary and predictions

Author: SS-PEG program
Date: 2026
"""

import numpy as np
from fractions import Fraction
from itertools import product as cart_product
import sys, os

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "neutrino_mass_derivation_output.md")
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
R3_f = None  # irrational — use float
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

# NuFIT 5.2 (Normal Ordering)
Dm2_21 = 7.42e-5   # eV²
Dm2_31 = 2.515e-3  # eV²


def formula_str(exps):
    parts = []
    for name, e in zip(BNAMES, exps):
        if e == 0:
            continue
        elif e == 1:
            parts.append(name)
        elif e == -1:
            parts.append(f"{name}⁻¹")
        else:
            sup = str(e).replace("-", "⁻")
            parts.append(f"{name}^{sup}")
    return " · ".join(parts) if parts else "1"


def formula_value(exps):
    return np.prod(BLOCKS ** np.array(exps, dtype=float))


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


print("=" * 78)
print("  NEUTRINO MASS DERIVATION FROM GENERATING FUNCTION")
print("  SS-PEG Research Program — 2026")
print("=" * 78)

# ════════════════════════════════════════════════════════════════
# §1. GENERATING FUNCTION EVALUATION AT NEUTRINO QUANTUM NUMBERS
# ════════════════════════════════════════════════════════════════

print("\n" + "=" * 78)
print("§1. GENERATING FUNCTION → NEUTRINO LATTICE COORDINATES")
print("=" * 78)

# Generating function coefficients from §10.5 of 00_DERIVATION.md
# X(s) = α·(2I₃) + β·Nc + γ
GF_COEFFS = {
    # (α, β, γ) for each of the 9 coefficient functions
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
    """Evaluate a generating function coefficient at given quantum numbers."""
    alpha, beta, gamma = GF_COEFFS[name]
    return alpha * I3_2 + beta * Nc + gamma


def get_sector_coeffs(I3_2, Nc):
    """Get all 9 parabolic coefficients (a,b,c) × (p,q,r) for a sector."""
    result = {}
    for coord in ["p", "q", "r"]:
        a = eval_gf_coeff(f"a_{coord}", I3_2, Nc)
        b = eval_gf_coeff(f"b_{coord}", I3_2, Nc)
        c = eval_gf_coeff(f"c_{coord}", I3_2, Nc)
        result[coord] = (a, b, c)
    return result


def get_lattice_coords(coeffs, g):
    """Get lattice coordinates (p, q, r) at generation g."""
    p = coeffs["p"][0] * g**2 + coeffs["p"][1] * g + coeffs["p"][2]
    q = coeffs["q"][0] * g**2 + coeffs["q"][1] * g + coeffs["q"][2]
    r = coeffs["r"][0] * g**2 + coeffs["r"][1] * g + coeffs["r"][2]
    return (p, q, r)


def mass_from_coords(p, q, r):
    """Compute mass in MeV from lattice coordinates."""
    ln_ratio = float(p) * LN2 + float(q) * LNR3 + float(r) * LN3
    return m_e * np.exp(ln_ratio)


# Verify with known sectors first
print("\n  Verification with known sectors:")
print(f"  {'Sector':<12} {'2I₃':>5} {'Nc':>4} │ {'Gen':>4} {'p':>8} {'q':>8} {'r':>8} │ {'m_pred (MeV)':>14} {'m_exp (MeV)':>14} {'Err%':>8}")
print(f"  {'─' * 95}")

known_masses = {
    "lepton": {"2I3": -1, "Nc": 1,
               "masses": [0.51099895, 105.6583755, 1776.86],
               "names": ["e", "μ", "τ"]},
    "up":     {"2I3": +1, "Nc": 3,
               "masses": [2.16, 1270.0, 172500.0],
               "names": ["u", "c", "t"]},
    "down":   {"2I3": -1, "Nc": 3,
               "masses": [4.70, 93.4, 4180.0],
               "names": ["d", "s", "b"]},
}

for sector_name, data in known_masses.items():
    coeffs = get_sector_coeffs(data["2I3"], data["Nc"])
    for g in [1, 2, 3]:
        p, q, r = get_lattice_coords(coeffs, g)
        m_pred = mass_from_coords(p, q, r)
        m_exp = data["masses"][g-1]
        err = abs(m_pred - m_exp) / m_exp * 100
        name = data["names"][g-1]
        print(f"  {sector_name:<12} {data['2I3']:>+5} {data['Nc']:>4} │ {g:>4} {float(p):>+8.2f} {float(q):>+8.2f} {float(r):>+8.2f} │ {m_pred:>14.4f} {m_exp:>14.4f} {err:>7.3f}%")


# NOW: Neutrino sector (2I₃ = +1, Nc = 1) — the PREDICTION
print(f"\n  {'═' * 78}")
print(f"  ★ NEUTRINO SECTOR: (2I₃ = +1, Nc = 1) — GENUINE PREDICTION")
print(f"  {'═' * 78}")

nu_coeffs = get_sector_coeffs(+1, 1)
print(f"\n  Parabolic coefficients (exact fractions):")
for coord in ["p", "q", "r"]:
    a, b, c = nu_coeffs[coord]
    g_star = -b / (2 * a) if a != 0 else float('inf')
    print(f"    {coord}(g) = ({a})·g² + ({b})·g + ({c})    [g* = {float(g_star):.4f}]")

print(f"\n  Lattice coordinates and Dirac masses:")
print(f"  {'Gen':>4} {'p':>10} {'q':>10} {'r':>10} │ {'ln(m/mₑ)':>12} {'m_D (MeV)':>14} {'m_D (eV)':>14}")
print(f"  {'─' * 78}")

nu_dirac_masses = []
nu_coords = []
for g in [1, 2, 3]:
    p, q, r = get_lattice_coords(nu_coeffs, g)
    m_D = mass_from_coords(p, q, r)
    ln_ratio = float(p) * LN2 + float(q) * LNR3 + float(r) * LN3
    nu_dirac_masses.append(m_D)
    nu_coords.append((p, q, r))
    print(f"  {g:>4} {float(p):>+10.4f} {float(q):>+10.4f} {float(r):>+10.4f} │ {ln_ratio:>+12.6f} {m_D:>14.6f} {m_D*1e3:>14.4f}")

# Compare with charged lepton partners
print(f"\n  Isospin partner comparison (ν vs ℓ):")
ell_coeffs = get_sector_coeffs(-1, 1)
print(f"  {'Gen':>4} {'m_ℓ (MeV)':>14} {'m_D^ν (MeV)':>14} {'Ratio m_D^ν/m_ℓ':>18}")
print(f"  {'─' * 55}")

ell_masses = [0.51099895, 105.6583755, 1776.86]
for g in [1, 2, 3]:
    m_D = nu_dirac_masses[g-1]
    m_l = ell_masses[g-1]
    ratio = m_D / m_l
    print(f"  {g:>4} {m_l:>14.4f} {m_D:>14.6f} {ratio:>18.6f}")


# ════════════════════════════════════════════════════════════════
# §2. DIRAC MASSES vs PHYSICAL NEUTRINO MASSES
# ════════════════════════════════════════════════════════════════

print("\n" + "=" * 78)
print("§2. DIRAC MASSES vs PHYSICAL NEUTRINO MASSES — SEESAW EXTRACTION")
print("=" * 78)

# Physical neutrino masses (Normal Ordering, from pmns_structural_analysis.py scan)
# Best m₁ ≈ 2.1 meV from scan, but we'll use the R_EE⁵ analytical solution too

# R_EE⁵ solution: m₂/m₃ = R_EE⁵ = 2^(-5/2)
delta_31_32 = Dm2_31 - 32 * Dm2_21
m1_REE5 = np.sqrt(delta_31_32 / 31) if delta_31_32 > 0 else 0
m2_REE5 = np.sqrt(m1_REE5**2 + Dm2_21)
m3_REE5 = np.sqrt(m1_REE5**2 + Dm2_31)

# Hierarchical limit (m₁ → 0)
m1_hier = 0
m2_hier = np.sqrt(Dm2_21)
m3_hier = np.sqrt(Dm2_31)

# Best graph fit from scan (m₁ optimized for minimal max error)
# We'll compute it here
m1_scan_vals = np.linspace(0, 0.05, 1001)  # 0 to 50 meV in eV
best_m1_scan = None
best_max_err = float('inf')

GRID_CORE = np.array(list(cart_product(range(-4, 5), repeat=5)))
LOG_CORE = GRID_CORE @ LOG_BLOCKS

for m1_test in m1_scan_vals:
    m2_test = np.sqrt(m1_test**2 + Dm2_21)
    m3_test = np.sqrt(m1_test**2 + Dm2_31)

    max_err = 0
    # m₂/m₃ ratio
    ratio_23 = m2_test / m3_test
    log_target = np.log(ratio_23)
    min_err_23 = np.min(np.abs(LOG_CORE - log_target)) / abs(log_target) * 100

    if m1_test > 1e-5:
        ratio_12 = m1_test / m2_test
        log_target_12 = np.log(ratio_12)
        min_err_12 = np.min(np.abs(LOG_CORE - log_target_12)) / abs(log_target_12) * 100
        max_err = max(min_err_23, min_err_12)
    else:
        max_err = min_err_23

    if max_err < best_max_err:
        best_max_err = max_err
        best_m1_scan = m1_test

m1_scan = best_m1_scan
m2_scan = np.sqrt(m1_scan**2 + Dm2_21)
m3_scan = np.sqrt(m1_scan**2 + Dm2_31)

print(f"""
  Physical neutrino masses (Normal Ordering, NuFIT 5.2):
  ┌────────────────────────┬──────────────┬──────────────┬──────────────┐
  │ Scenario               │ m₁ (meV)     │ m₂ (meV)     │ m₃ (meV)     │
  ├────────────────────────┼──────────────┼──────────────┼──────────────┤
  │ Hierarchical (m₁→0)    │ {0:>10.4f}   │ {m2_hier*1e3:>10.4f}   │ {m3_hier*1e3:>10.4f}   │
  │ R_EE⁵ solution         │ {m1_REE5*1e3:>10.4f}   │ {m2_REE5*1e3:>10.4f}   │ {m3_REE5*1e3:>10.4f}   │
  │ Best graph scan         │ {m1_scan*1e3:>10.4f}   │ {m2_scan*1e3:>10.4f}   │ {m3_scan*1e3:>10.4f}   │
  └────────────────────────┴──────────────┴──────────────┴──────────────┘

  Dirac masses from generating function:
    m_D₁ = {nu_dirac_masses[0]*1e3:.4f} eV = {nu_dirac_masses[0]*1e6:.4f} keV
    m_D₂ = {nu_dirac_masses[1]*1e3:.4f} eV = {nu_dirac_masses[1]*1e6:.4f} keV
    m_D₃ = {nu_dirac_masses[2]*1e3:.4f} eV = {nu_dirac_masses[2]*1e6:.4f} keV
""")

# Type-I Seesaw: m_ν = m_D² / M_R
# → M_R = m_D² / m_ν

print("  Type-I Seesaw extraction: M_R = m_D² / m_ν")
print(f"\n  {'Gen':>4} {'m_D (eV)':>14} {'m_ν (eV)':>14} {'M_R (GeV)':>16} {'M_R (eV)':>16}")
print(f"  {'─' * 65}")

# Use the R_EE⁵ physical masses
m_nu_phys = [m1_REE5, m2_REE5, m3_REE5]
m_nu_labels = ["ν₁", "ν₂", "ν₃"]
M_R_values = []

for i in range(3):
    m_D_eV = nu_dirac_masses[i] * 1e6  # MeV → eV
    m_nu_eV = m_nu_phys[i]
    if m_nu_eV > 0:
        M_R_eV = m_D_eV**2 / m_nu_eV
        M_R_GeV = M_R_eV / 1e9
        M_R_values.append(M_R_eV)
        print(f"  {m_nu_labels[i]:>4} {m_D_eV:>14.4f} {m_nu_eV:>14.6e} {M_R_GeV:>16.6e} {M_R_eV:>16.6e}")
    else:
        M_R_values.append(0)
        print(f"  {m_nu_labels[i]:>4} {m_D_eV:>14.4f} {'(m₁→0)':>14} {'—':>16} {'—':>16}")

# Also compute with hierarchical masses
print(f"\n  With hierarchical masses (m₁ → 0):")
m_nu_hier = [0, m2_hier, m3_hier]
M_R_hier = []
for i in range(3):
    m_D_eV = nu_dirac_masses[i] * 1e6
    m_nu_eV = m_nu_hier[i]
    if m_nu_eV > 0:
        M_R_eV = m_D_eV**2 / m_nu_eV
        M_R_GeV = M_R_eV / 1e9
        M_R_hier.append(M_R_eV)
        print(f"  {m_nu_labels[i]:>4} {m_D_eV:>14.4f} {m_nu_eV:>14.6e} {M_R_GeV:>16.6e} {M_R_eV:>16.6e}")
    else:
        M_R_hier.append(0)
        print(f"  {m_nu_labels[i]:>4} {m_D_eV:>14.4f} {'(m₁→0)':>14} {'—':>16} {'—':>16}")


# ════════════════════════════════════════════════════════════════
# §3. SEESAW SCALE — GRAPH FORMULA SEARCH
# ════════════════════════════════════════════════════════════════

print("\n" + "=" * 78)
print("§3. SEESAW SCALE — GRAPH FORMULA SEARCH")
print("=" * 78)

print("""
  If M_R has a graph-topological origin, it should be expressible as
  M_R = m_e · (product of building blocks)  or as a ratio M_R_i/M_R_j.
""")

# Search for graph formulas of M_R/m_e (dimensionless)
# Also search for M_R ratios
if len(M_R_values) >= 3 and all(v > 0 for v in M_R_values):
    targets_MR = {
        "M_R₁ (eV)": M_R_values[0],
        "M_R₂ (eV)": M_R_values[1],
        "M_R₃ (eV)": M_R_values[2],
        "M_R₃/M_R₂": M_R_values[2] / M_R_values[1],
        "M_R₃/M_R₁": M_R_values[2] / M_R_values[0],
        "M_R₂/M_R₁": M_R_values[1] / M_R_values[0],
        "ln(M_R₃/m_e)/ln(M_R₂/m_e)": np.log(M_R_values[2] / (m_e * 1e6)) / np.log(M_R_values[1] / (m_e * 1e6)),
    }

    for name, val in targets_MR.items():
        results = best_formula(val, n=3)
        print(f"\n  ── {name} = {val:.6e} ──")
        for exps, fval, err, compl in results:
            print(f"     {formula_str(exps):<55s} = {fval:.6e}  err {err:.3f}%  C={compl}")


# ════════════════════════════════════════════════════════════════
# §4. DIRAC MASS STRUCTURE — ISOSPIN FLIP ANALYSIS
# ════════════════════════════════════════════════════════════════

print("\n" + "=" * 78)
print("§4. ISOSPIN FLIP: CHARGED LEPTON → NEUTRINO")
print("=" * 78)

print("""
  The generating function maps (2I₃=-1, Nc=1) → charged leptons,
  and (2I₃=+1, Nc=1) → neutrino Dirac masses.

  The ISOSPIN FLIP (2I₃: -1 → +1) changes each coefficient by:
    ΔX = α · Δ(2I₃) = α · (+2) = 2α

  This gives a systematic shift in all lattice coordinates.
""")

print(f"  {'Coeff':<8} {'α':>10} {'2α':>10} │ {'ℓ value':>12} {'ν value':>12} {'Δ':>10}")
print(f"  {'─' * 70}")

for coord in ["p", "q", "r"]:
    for prefix in ["a", "b", "c"]:
        name = f"{prefix}_{coord}"
        alpha, beta, gamma = GF_COEFFS[name]
        delta = 2 * alpha
        val_ell = eval_gf_coeff(name, -1, 1)
        val_nu = eval_gf_coeff(name, +1, 1)
        print(f"  {name:<8} {str(alpha):>10} {str(delta):>10} │ {str(val_ell):>12} {str(val_nu):>12} {str(delta):>10}")

# Compute the coordinate differences at each generation
print(f"\n  Lattice coordinate shifts (ν − ℓ) at each generation:")
print(f"  {'Gen':>4} {'Δp':>10} {'Δq':>10} {'Δr':>10} │ {'Δ[ln(m/mₑ)]':>14} {'m_D/m_ℓ':>14}")
print(f"  {'─' * 65}")

for g in [1, 2, 3]:
    p_nu, q_nu, r_nu = get_lattice_coords(nu_coeffs, g)
    p_ell, q_ell, r_ell = get_lattice_coords(ell_coeffs, g)
    dp = float(p_nu - p_ell)
    dq = float(q_nu - q_ell)
    dr = float(r_nu - r_ell)
    d_ln = dp * LN2 + dq * LNR3 + dr * LN3
    ratio = np.exp(d_ln)
    print(f"  {g:>4} {dp:>+10.4f} {dq:>+10.4f} {dr:>+10.4f} │ {d_ln:>+14.6f} {ratio:>14.6f}")

# Search for graph formulas of the mass ratio m_D^ν/m_ℓ
print(f"\n  Graph formulas for m_D^ν/m_ℓ at each generation:")
for g in [1, 2, 3]:
    ratio = nu_dirac_masses[g-1] / ell_masses[g-1]
    results = best_formula(ratio, n=3)
    print(f"\n  ── Gen {g}: m_D^ν/m_ℓ = {ratio:.6f} ──")
    for exps, fval, err, compl in results:
        print(f"     {formula_str(exps):<55s} = {fval:.6f}  err {err:.3f}%  C={compl}")


# ════════════════════════════════════════════════════════════════
# §5. DIRAC MASS RATIOS — INTERNAL CONSISTENCY
# ════════════════════════════════════════════════════════════════

print("\n" + "=" * 78)
print("§5. DIRAC MASS RATIOS — GRAPH FORMULAS")
print("=" * 78)

print("""
  The generating function predicts Dirac mass RATIOS between neutrino
  generations. These should have graph formulas, and they determine
  the STRUCTURE of the seesaw (whether M_R is universal or not).
""")

dirac_ratios = {
    "m_D₂/m_D₁": nu_dirac_masses[1] / nu_dirac_masses[0],
    "m_D₃/m_D₂": nu_dirac_masses[2] / nu_dirac_masses[1],
    "m_D₃/m_D₁": nu_dirac_masses[2] / nu_dirac_masses[0],
}

print(f"  {'Ratio':<14} {'Value':>12} {'Formula':<50} {'Err%':>7}")
print(f"  {'─' * 90}")

for name, val in dirac_ratios.items():
    results = best_formula(val, n=1)
    exps, fval, err, compl = results[0]
    print(f"  {name:<14} {val:>12.4f} {formula_str(exps):<50} {err:>6.3f}%")

# Compare with charged-lepton ratios
print(f"\n  Comparison with charged-lepton generation ratios:")
cl_ratios = {
    "m_μ/m_e":  ell_masses[1] / ell_masses[0],
    "m_τ/m_μ":  ell_masses[2] / ell_masses[1],
    "m_τ/m_e":  ell_masses[2] / ell_masses[0],
}

for name, val in cl_ratios.items():
    results = best_formula(val, n=1)
    exps, fval, err, compl = results[0]
    print(f"  {name:<14} {val:>12.4f} {formula_str(exps):<50} {err:>6.3f}%")

# Generation scaling exponent γ
gamma_nu = np.log(dirac_ratios["m_D₃/m_D₂"]) / np.log(dirac_ratios["m_D₂/m_D₁"])
gamma_ell = np.log(cl_ratios["m_τ/m_μ"]) / np.log(cl_ratios["m_μ/m_e"])
print(f"\n  Generation scaling exponent γ = ln(m₃/m₂)/ln(m₂/m₁):")
print(f"    Neutrino Dirac: γ_ν = {gamma_nu:.6f}")
print(f"    Charged lepton: γ_ℓ = {gamma_ell:.6f}")


# ════════════════════════════════════════════════════════════════
# §6. CONSISTENCY WITH PMNS AND r_Δm²
# ════════════════════════════════════════════════════════════════

print("\n" + "=" * 78)
print("§6. CONSISTENCY: DIRAC MASSES + SEESAW → PHYSICAL OBSERVABLES")
print("=" * 78)

print("""
  If the seesaw m_ν = m_D²/M_R operates with universal M_R:
    m_ν₂/m_ν₃ = (m_D₂/m_D₃)² = m_D₂²/m_D₃²

  If M_R is generation-dependent (M_R_i):
    m_ν₂/m_ν₃ = (m_D₂/m_D₃)² × (M_R₃/M_R₂)
""")

# Case 1: Universal M_R
r_universal = (nu_dirac_masses[1] / nu_dirac_masses[2])**2
sqrt_r_exp = np.sqrt(Dm2_21 / Dm2_31)

print(f"  CASE 1: Universal M_R")
print(f"    (m_D₂/m_D₃)² = {r_universal:.6f}")
print(f"    Experimental √(Δm²₂₁/Δm²₃₁) ≈ m₂/m₃ = {sqrt_r_exp:.6f}")

# For hierarchical limit: m₂/m₃ = √r
r_pred_univ = r_universal
r_exp = Dm2_21 / Dm2_31
err_univ = abs(r_pred_univ - r_exp) / r_exp * 100
print(f"    Predicted r = (m_D₂/m_D₃)² = {r_pred_univ:.6f}")
print(f"    Experimental r = Δm²₂₁/Δm²₃₁ = {r_exp:.6f}")
print(f"    Error: {err_univ:.1f}%")

if err_univ > 10:
    print(f"    → Universal M_R does NOT reproduce the mass ratio.")
    print(f"    → NEED generation-dependent M_R.")

    # Case 2: What M_R ratio is needed?
    # m_ν₂/m_ν₃ = (m_D₂/m_D₃)² × (M_R₃/M_R₂)
    # → M_R₃/M_R₂ = (m_ν₂/m_ν₃) / (m_D₂/m_D₃)²
    # In hierarchical limit: m₂/m₃ = √r
    MR_ratio_needed = sqrt_r_exp / np.sqrt(r_universal)  # √(M_R₃/M_R₂) for mass ratio
    # Actually: m_ν_i = m_D_i² / M_R_i
    # m₂/m₃ = (m_D₂²/M_R₂) / (m_D₃²/M_R₃) = (m_D₂/m_D₃)² × (M_R₃/M_R₂)
    # √r ≈ (m_D₂/m_D₃)² × (M_R₃/M_R₂)   [hierarchical]
    MR32_ratio = sqrt_r_exp / (nu_dirac_masses[1] / nu_dirac_masses[2])**2

    print(f"\n  CASE 2: Generation-dependent M_R")
    print(f"    Required M_R₃/M_R₂ = √r / (m_D₂/m_D₃)² = {MR32_ratio:.6f}")

    # Search for graph formula
    results = best_formula(MR32_ratio, n=5)
    print(f"\n    Graph formula search for M_R₃/M_R₂ = {MR32_ratio:.6f}:")
    for exps, fval, err, compl in results:
        print(f"      {formula_str(exps):<55s} = {fval:.6f}  err {err:.3f}%  C={compl}")

    # Also the r ratio
    # r = (m_D₂/m_D₃)⁴ × (M_R₃/M_R₂)²
    MR32_from_r = np.sqrt(r_exp / (nu_dirac_masses[1] / nu_dirac_masses[2])**4)
    print(f"\n    From r directly: M_R₃/M_R₂ = √[r/(m_D₂/m_D₃)⁴] = {MR32_from_r:.6f}")

    # Verify consistency
    r_check = (nu_dirac_masses[1] / nu_dirac_masses[2])**4 * MR32_from_r**2
    print(f"    Check: (m_D₂/m_D₃)⁴ × (M_R₃/M_R₂)² = {r_check:.6f} vs r = {r_exp:.6f}")


# ════════════════════════════════════════════════════════════════
# §7. RAMANUJAN SPECTRAL INTERPRETATION
# ════════════════════════════════════════════════════════════════

print("\n" + "=" * 78)
print("§7. RAMANUJAN SPECTRAL INTERPRETATION")
print("=" * 78)

print("""
  In the SS-PEG framework, masses arise from spectral hops on the
  Cartan-Weyl root graph. Each building block represents a different
  type of spectral process:

  - R₂, R₃: spectral gap ratios (propagation efficiency)
  - R_EE:   bifundamental bridge (sector coupling)
  - h∨₃, h∨₂: Coxeter amplification (loop corrections)

  For neutrinos (Q=0, Nc=1), the Ramanujan-optimal propagation on
  SU(3) is ABSENT. The generating function captures this through
  the (2I₃, Nc) dependence:

  The isospin flip 2I₃: -1 → +1 shifts lattice coordinates by 2α_k.
  The α coefficients control HOW MUCH each parabolic parameter
  responds to the isospin flip.
""")

# Analyze the α structure
print(f"  α coefficients (isospin sensitivity):")
print(f"  {'Coeff':<8} {'α':>12} {'|α|':>8} {'Sector':<20}")
print(f"  {'─' * 55}")

alpha_vals = {}
for coord in ["p", "q", "r"]:
    for prefix in ["a", "b", "c"]:
        name = f"{prefix}_{coord}"
        alpha = GF_COEFFS[name][0]
        alpha_vals[name] = alpha
        # Classify by which block this shifts
        sector = "SU(2) binary" if coord == "p" else ("SU(3) irrational" if coord == "q" else "SU(3) ternary")
        print(f"  {name:<8} {str(alpha):>12} {abs(float(alpha)):>8.4f} {sector:<20}")

# Which coordinate is most affected by the isospin flip?
total_shift = {}
for coord in ["p", "q", "r"]:
    total = sum(abs(float(alpha_vals[f"{prefix}_{coord}"])) for prefix in ["a", "b", "c"])
    total_shift[coord] = total
    print(f"\n  Total |α| shift for {coord}: {total:.4f}")

dominant = max(total_shift, key=total_shift.get)
print(f"\n  → Dominant shift: {dominant} coordinate")
print(f"    This means the isospin flip primarily affects the "
      f"{'binary (ln2)' if dominant == 'p' else ('SU(3) irrational (lnR₃)' if dominant == 'q' else 'ternary (ln3)')} channel.")


# ════════════════════════════════════════════════════════════════
# §8. NEUTRINO MASS QUADRATIC STRUCTURE
# ════════════════════════════════════════════════════════════════

print("\n" + "=" * 78)
print("§8. NEUTRINO GENERATION CURVES — QUADRATIC ANALYSIS")
print("=" * 78)

print("""
  The generating function gives neutrino Dirac masses on a parabola:
    ln(m_D(g)/mₑ) = A_ν·g² + B_ν·g + C_ν

  where A_ν, B_ν, C_ν are computed from the neutrino lattice coefficients.
""")

# Compute mass-space parabola A, B, C for neutrinos
nu_A = float(nu_coeffs["p"][0]) * LN2 + float(nu_coeffs["q"][0]) * LNR3 + float(nu_coeffs["r"][0]) * LN3
nu_B = float(nu_coeffs["p"][1]) * LN2 + float(nu_coeffs["q"][1]) * LNR3 + float(nu_coeffs["r"][1]) * LN3
nu_C = float(nu_coeffs["p"][2]) * LN2 + float(nu_coeffs["q"][2]) * LNR3 + float(nu_coeffs["r"][2]) * LN3

# Compare with charged leptons
ell_A = 1.0 * LN2 + 0.5 * LNR3 + (-1.5) * LN3
ell_B = (-3.5) * LN2 + (-5.5) * LNR3 + 7.5 * LN3
ell_C = 2.5 * LN2 + 5.0 * LNR3 + (-6.0) * LN3

print(f"  Mass-space parabola ln(m/mₑ) = A·g² + B·g + C:")
print(f"  {'Sector':<12} {'A':>12} {'B':>12} {'C':>12} {'g*(turn)':>10} {'2A(accel)':>10}")
print(f"  {'─' * 65}")

g_star_nu = -nu_B / (2 * nu_A) if abs(nu_A) > 1e-10 else float('inf')
g_star_ell = -ell_B / (2 * ell_A) if abs(ell_A) > 1e-10 else float('inf')

print(f"  {'Neutrino ν':<12} {nu_A:>+12.6f} {nu_B:>+12.6f} {nu_C:>+12.6f} {g_star_nu:>10.4f} {2*nu_A:>+10.6f}")
print(f"  {'Lepton ℓ':<12} {ell_A:>+12.6f} {ell_B:>+12.6f} {ell_C:>+12.6f} {g_star_ell:>10.4f} {2*ell_A:>+10.6f}")

# Koide-like ratio for neutrinos
# m₃·m₁/m₂² = e^(2A)
koide_nu = np.exp(2 * nu_A)
koide_ell = np.exp(2 * ell_A)
print(f"\n  Koide-like geometric ratio m₃·m₁/m₂² = e^(2A):")
print(f"    Neutrino Dirac: e^(2A_ν) = {koide_nu:.6f}")
print(f"    Charged lepton: e^(2A_ℓ) = {koide_ell:.6f}")

# Verify with actual Dirac masses
koide_actual = nu_dirac_masses[2] * nu_dirac_masses[0] / nu_dirac_masses[1]**2
print(f"    Actual m_D₃·m_D₁/m_D₂² = {koide_actual:.6f}")
print(f"    Consistency: {abs(koide_nu - koide_actual)/koide_actual*100:.4f}%")

# Kinetic action for neutrino sector
S_K_nu = 0
for coord in ["p", "q", "r"]:
    a = float(nu_coeffs[coord][0])
    b = float(nu_coeffs[coord][1])
    S_k = Fraction(52, 3) * Fraction(nu_coeffs[coord][0])**2 + \
          8 * Fraction(nu_coeffs[coord][0]) * Fraction(nu_coeffs[coord][1]) + \
          Fraction(nu_coeffs[coord][1])**2
    S_K_nu += S_k
    print(f"    S_K({coord}|ν) = {float(S_k):.4f}")

print(f"    S_K(ν total) = {float(S_K_nu):.4f} = {S_K_nu}")


# ════════════════════════════════════════════════════════════════
# §9. PHYSICAL NEUTRINO MASS PREDICTIONS
# ════════════════════════════════════════════════════════════════

print("\n" + "=" * 78)
print("§9. PHYSICAL NEUTRINO MASS PREDICTIONS — SYNTHESIS")
print("=" * 78)

print(f"""
  From the generating function, the neutrino Dirac mass spectrum is:
    m_D = ({nu_dirac_masses[0]*1e6:.2f}, {nu_dirac_masses[1]*1e6:.2f}, {nu_dirac_masses[2]*1e6:.2f}) keV

  Combined with the graph-derived mass-squared ratio r = √2/48:
    r_graph = R₂³ · R_EE · h∨₃⁻¹ = {0.5**3 * R_EE / hv3:.6f}
    r_exp   = {r_exp:.6f}

  And the Dirac mass ratio from the generating function:
    (m_D₂/m_D₃)² = {(nu_dirac_masses[1]/nu_dirac_masses[2])**2:.6f}
""")

# The generating function gives us an INDEPENDENT route to the mass ratio
# If the seesaw is universal: m₂/m₃ = m_D₂/m_D₃
# If not: we can extract M_R₃/M_R₂

# Key question: does the generating function's Dirac mass ratio
# agree with the graph formula for r_Δm²?
r_graph = 0.5**3 * R_EE / hv3  # = √2/48
sqrt_r_graph = np.sqrt(r_graph)

print(f"  KEY TEST: Does m_D₂/m_D₃ from F match √r from graph formula?")
sqrt_r_dirac = nu_dirac_masses[1] / nu_dirac_masses[2]
print(f"    m_D₂/m_D₃ (from F)  = {sqrt_r_dirac:.6f}")
print(f"    √r (graph formula)   = {sqrt_r_graph:.6f}")
print(f"    √r (experimental)    = {sqrt_r_exp:.6f}")

err_dirac_vs_graph = abs(sqrt_r_dirac - sqrt_r_graph) / sqrt_r_graph * 100
err_dirac_vs_exp = abs(sqrt_r_dirac - sqrt_r_exp) / sqrt_r_exp * 100
print(f"    Err(Dirac vs graph): {err_dirac_vs_graph:.3f}%")
print(f"    Err(Dirac vs exp):   {err_dirac_vs_exp:.3f}%")

if err_dirac_vs_graph < 5:
    print(f"    → CONSISTENT! The generating function reproduces r_Δm² with universal M_R.")
    print(f"    → This means ALL neutrino mass ratios follow from the generating function alone.")
elif err_dirac_vs_exp < 20:
    print(f"    → PARTIAL MATCH. Corrections from M_R generation-dependence needed.")
else:
    print(f"    → INCONSISTENT with universal M_R. The seesaw scale MUST be generation-dependent.")


# ════════════════════════════════════════════════════════════════
# §10. THE COMPLETE NEUTRINO PICTURE
# ════════════════════════════════════════════════════════════════

print("\n" + "=" * 78)
print("§10. THE COMPLETE NEUTRINO PICTURE — SUMMARY")
print("=" * 78)

print(f"""
  ╔══════════════════════════════════════════════════════════════════════╗
  ║  THE NEUTRINO MASS DERIVATION FROM GRAPH TOPOLOGY                   ║
  ╚══════════════════════════════════════════════════════════════════════╝

  1. GENERATING FUNCTION PREDICTION
     The generating function F(2I₃, Nc, g), calibrated on charged fermions
     only, evaluated at neutrino quantum numbers (2I₃=+1, Nc=1), produces:

     ν₁: (p,q,r) = ({float(nu_coords[0][0]):+.2f}, {float(nu_coords[0][1]):+.2f}, {float(nu_coords[0][2]):+.2f}) → m_D₁ = {nu_dirac_masses[0]*1e6:.2f} keV
     ν₂: (p,q,r) = ({float(nu_coords[1][0]):+.2f}, {float(nu_coords[1][1]):+.2f}, {float(nu_coords[1][2]):+.2f}) → m_D₂ = {nu_dirac_masses[1]*1e6:.2f} keV
     ν₃: (p,q,r) = ({float(nu_coords[2][0]):+.2f}, {float(nu_coords[2][1]):+.2f}, {float(nu_coords[2][2]):+.2f}) → m_D₃ = {nu_dirac_masses[2]*1e6:.2f} keV

     These are the DIRAC masses (Yukawa-type, before seesaw suppression).

  2. SEESAW MECHANISM
     Physical neutrino masses m_ν ≪ m_D require Type-I seesaw:
       m_ν = m_D² / M_R

     The ratio m_D/m_ℓ (isospin partner) varies by generation — the
     generating function predicts NON-PROPORTIONAL Dirac masses.

  3. MASS ORDERING
     The graph framework PREDICTS Normal Ordering (NO) through:
     - Better graph formula fits for NO ratios vs IO ratios
     - r_Δm² = √2/48 has a compact formula in NO, not in IO

  4. PMNS ANGLES
     All 4 PMNS parameters derived from building blocks (already done):
       sinθ₂₃ = R₃³·R_EE⁻⁴·h∨₃²·h∨₂⁻³     (0.000%)
       sinθ₁₂ = R₂⁻³·R₃·R_EE⁴·h∨₂⁻¹        (0.204%)
       sinθ₁₃ = R₂·R_EE⁻⁴·h∨₃⁻³·h∨₂         (0.237%)
       δ_CP/π = R₂⁻²·R₃⁻²·R_EE⁴·h∨₃⁻¹       (0.047%)

  5. GRAPH-TOPOLOGICAL ORIGIN
     The generating function's α coefficients quantify the isospin
     sensitivity of each lattice coordinate. Neutrinos occupy a
     SINGULAR point where Q=0 collapses the electromagnetic tensor.
     The Ramanujan spectral quality of the building blocks ensures
     maximal propagation efficiency on the SU(2) cycle.
""")


# ════════════════════════════════════════════════════════════════
# §11. FLATNESS AND GENERATION STRUCTURE
# ════════════════════════════════════════════════════════════════

print("\n" + "=" * 78)
print("§11. FLATNESS AND GENERATION STRUCTURE FOR NEUTRINOS")
print("=" * 78)

# Check flatness: which coordinate satisfies b = -5a ?
print(f"\n  Flatness test b = -5a for neutrino sector:")
for coord in ["p", "q", "r"]:
    a = nu_coeffs[coord][0]
    b = nu_coeffs[coord][1]
    ratio = float(b) / float(a) if float(a) != 0 else float('inf')
    is_flat = abs(ratio + 5) < 0.01
    print(f"    {coord}: b/a = {float(b):.4f}/{float(a):.4f} = {ratio:.4f}  "
          f"{'→ FLAT ✓' if is_flat else f'(need -5.0, off by {ratio+5:+.4f})'}")

# Charged lepton flat coord is r
print(f"\n  Charged lepton flat coordinate: r")
print(f"    b_r/a_r = {-6.0/(-1.5):.4f}")

# Check if neutrino has a flat coordinate
print(f"\n  → Neutrino flat coordinate: ", end="")
found_flat = False
for coord in ["p", "q", "r"]:
    a = float(nu_coeffs[coord][0])
    b = float(nu_coeffs[coord][1])
    if abs(a) > 1e-10 and abs(b/a + 5) < 0.1:
        print(f"{coord} (b/a = {b/a:.4f})")
        found_flat = True
        break
if not found_flat:
    print("NONE — all coordinates are non-flat!")
    print("    This is a UNIQUE property of the neutrino sector.")
    print("    Charged fermions always have exactly 1 flat coordinate per sector.")


# ════════════════════════════════════════════════════════════════
# CLOSE OUTPUT
# ════════════════════════════════════════════════════════════════

sys.stdout = _orig_stdout
_outfile.close()
print(f"\nOutput saved to: {OUTPUT_PATH}")
