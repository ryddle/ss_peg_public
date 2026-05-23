#!/usr/bin/env python3
"""
velocity_analysis.py — Integration & Boundary Conditions
=========================================================

Given that ln(m(g)/mₑ) = A·g² + B·g + C per sector,
and A (acceleration/curvature) is constrained by:
  - Isocurvature: a_q(ℓ)=a_q(u), a_r(ℓ)=a_r(u)
  - Conservation: |a_q + a_r| = 1
  - Down linearity: a_q(d) = 0

What fixes B (velocity) and C (offset)?

This script performs the formal integration analysis:
  1. What boundary conditions determine B and C?
  2. What is the MINIMAL input to derive all 9 masses?
  3. Is there a regularity principle that fixes B?

Author: SS-PEG program
Date: April 2026
"""

import numpy as np
import matplotlib.pyplot as plt
import sys, os

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "velocity_analysis_output.md")
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
# DATA
# ════════════════════════════════════════════════════════════════

LN2 = np.log(2)
R3 = (np.sqrt(57) - 3) / (2 * np.sqrt(17))
LNR3 = np.log(R3)
LN3 = np.log(3)

# Lattice coordinates from consistent_lattice.py
COORDS = {
    'e':  ( 0.0,  0.0,  0.0),
    'mu': (-0.5, -4.0,  3.0),
    'tau':( 1.0, -7.0,  3.0),
    'u':  (-1.5, -6.0, -1.0),
    'c':  ( 0.5, -7.0,  3.0),
    't':  ( 6.0, -7.0,  4.0),
    'd':  (-0.5, -8.0, -2.0),
    's':  ( 1.5, -7.0,  0.0),
    'b':  ( 1.5, -6.0,  4.0),
}

SECTORS = {
    'lepton': ['e', 'mu', 'tau'],
    'up':     ['u', 'c', 't'],
    'down':   ['d', 's', 'b'],
}

# Quadratic coefficients from transition_curves.py
# coefficients[sector][component] = (a, b, c) where comp(g) = a·g² + b·g + c
coefficients = {
    'lepton': [( 1.0, -3.5,  2.5),   # p
               ( 0.5, -5.5,  5.0),   # q
               (-1.5,  7.5, -6.0)],  # r
    'up':     [( 1.75, -3.25, 0.0),
               ( 0.5,  -2.5, -4.0),
               (-1.5,   8.5, -8.0)],
    'down':   [(-1.0,   5.0, -4.5),
               ( 0.0,   1.0, -9.0),
               ( 1.0,  -1.0, -2.0)],
}

# Experimental masses (MeV)
MASSES = {
    'e': 0.51100, 'mu': 105.658, 'tau': 1776.86,
    'u': 2.16, 'c': 1270.0, 't': 172760.0,
    'd': 4.67, 's': 93.4, 'b': 4180.0
}
m_e = MASSES['e']

# ════════════════════════════════════════════════════════════════
# §1  FORMAL INTEGRATION STRUCTURE
# ════════════════════════════════════════════════════════════════

print("=" * 72)
print("  VELOCITY ANALYSIS — Integration & Boundary Conditions")
print("=" * 72)

print("\n§1. FORMAL INTEGRATION — From Acceleration to Mass")
print("-" * 72)
print("""
  The mass parabola: ln(m(g)/mₑ) = A·g² + B·g + C

  Reading this as kinematics in "generation time" g:
    Position:     x(g) = A·g² + B·g + C     [= ln(m/mₑ)]
    Velocity:     v(g) = 2A·g + B            [= d/dg ln(m/mₑ)]
    Acceleration: a(g) = 2A = constant       [= d²/dg² ln(m/mₑ)]

  The acceleration 2A is CONSTANT (uniform) — each sector has a
  fixed "generational acceleration" in mass-space.

  INTEGRATION PROCESS:
  ────────────────────
  Given: acceleration 2A (from graph topology)
  
  First integration:
    v(g) = ∫ 2A dg = 2A·g + B
    → B is the FIRST integration constant = initial velocity v(0)
    
  Second integration:
    x(g) = ∫ v(g) dg = A·g² + B·g + C
    → C is the SECOND integration constant = initial position x(0)

  To fix B and C we need TWO boundary conditions per sector.
""")

# Compute mass-space A, B, C for each sector
mass_ABC = {}
for sector_name in ["lepton", "up", "down"]:
    ap, bp, cp = coefficients[sector_name][0]
    aq, bq, cq = coefficients[sector_name][1]
    ar, br, cr = coefficients[sector_name][2]

    A = ap * LN2 + aq * LNR3 + ar * LN3
    B = bp * LN2 + bq * LNR3 + br * LN3
    C = cp * LN2 + cq * LNR3 + cr * LN3
    mass_ABC[sector_name] = (A, B, C)


# ════════════════════════════════════════════════════════════════
# §2  BOUNDARY CONDITIONS — What fixes B and C?
# ════════════════════════════════════════════════════════════════

print("\n§2. BOUNDARY CONDITIONS — What Fixes the Integration Constants?")
print("-" * 72)
print("""
  Each sector needs 2 boundary conditions. Natural choices:

  OPTION A: "Fix gen-1 and gen-2" (know the 2 lightest masses)
    x(1) = ln(m₁/mₑ)  →  A + B + C = known
    x(2) = ln(m₂/mₑ)  →  4A + 2B + C = known

    Solving: B = ln(m₂/mₑ) - 3A - ln(m₁/mₑ) = ln(m₂/m₁) - 3A
             C = ln(m₁/mₑ) - A - B

    → gen-3 mass is PREDICTED:
      ln(m₃/mₑ) = 9A + 3B + C = 2A + 2·ln(m₂/mₑ) - ln(m₁/mₑ)

    Equivalently:  ln(m₃/m₁) = 2A + 2·ln(m₂/m₁)
    Or:            m₃/m₁ = e^(2A) · (m₂/m₁)²
""")

print("  ★ PREDICTION FORMULA: m₃/m₁ = e^(2A) · (m₂/m₁)²")
print()
print(f"  {'Sector':<10} {'A':>10} {'e^(2A)':>10} {'(m₂/m₁)²':>12} "
      f"{'Predicted':>12} {'Actual m₃/m₁':>14} {'Error':>8}")
print(f"  {'─'*70}")

for sector_name, particles in [("lepton", ["e","mu","tau"]),
                                 ("up", ["u","c","t"]),
                                 ("down", ["d","s","b"])]:
    A = mass_ABC[sector_name][0]
    m1, m2, m3 = [MASSES[p] for p in particles]
    predicted = np.exp(2*A) * (m2/m1)**2
    actual = m3/m1
    err = (predicted - actual) / actual * 100
    print(f"  {sector_name:<10} {A:>+10.6f} {np.exp(2*A):>10.6f} {(m2/m1)**2:>12.2f} "
          f"{predicted:>12.2f} {actual:>14.2f} {err:>+8.4f}%")

print(f"""
  The formula m₃/m₁ = e^(2A) · (m₂/m₁)² is EXACT by construction
  (3 points determine a unique quadratic). The errors are from
  the lattice discretization (~0.3%).

  BUT: this formula tells us something profound.
  Rearranging: m₃ · m₁ / m₂² = e^(2A)

  ★ The "Koide-like" ratio m₃·m₁/m₂² is FIXED by the curvature A!
    It equals e^(2A) — a pure function of graph topology.
""")

print(f"  Koide-like geometric ratio m₃·m₁/m₂²:")
for sector_name, particles in [("lepton", ["e","mu","tau"]),
                                 ("up", ["u","c","t"]),
                                 ("down", ["d","s","b"])]:
    A = mass_ABC[sector_name][0]
    m1, m2, m3 = [MASSES[p] for p in particles]
    ratio = m3 * m1 / m2**2
    pred = np.exp(2*A)
    print(f"  {sector_name:<10}: m₃·m₁/m₂² = {ratio:.6f},  e^(2A) = {pred:.6f}")

print()


# ════════════════════════════════════════════════════════════════
# §3  WHAT DETERMINES A? — Curvature from Graph
# ════════════════════════════════════════════════════════════════

print("\n§3. CURVATURE A FROM GRAPH — Decomposition")
print("-" * 72)
print("""
  A = aₚ·ln2 + a_q·lnR₃ + a_r·ln3

  Constraints determine a_q, a_r:
    ℓ/u: a_q = +0.5, a_r = -1.5  (isocurvature + conservation)
    d:   a_q = 0,    a_r = +1.0  (linearity + conservation)

  Only aₚ differs per sector:
    aₚ(ℓ) = 1.0,  aₚ(u) = 1.75,  aₚ(d) = -1.0

  So: A = aₚ·ln2 + FIXED_PART
  where FIXED_PART depends only on the sector type (normal/inverted).
""")

# Fixed parts
FIXED_lu = 0.5 * LNR3 + (-1.5) * LN3     # for ℓ and u
FIXED_d  = 0.0 * LNR3 + 1.0 * LN3        # for down

print(f"  Fixed part (ℓ/u): 0.5·lnR₃ − 1.5·ln3 = {FIXED_lu:.6f}")
print(f"  Fixed part (d):   0.0·lnR₃ + 1.0·ln3 = {FIXED_d:.6f}")
print()

for sector_name in ["lepton", "up", "down"]:
    ap = coefficients[sector_name][0][0]
    A = mass_ABC[sector_name][0]
    fixed = FIXED_lu if sector_name != "down" else FIXED_d
    print(f"  A_{sector_name[:1]} = {ap:+.2f}·ln2 + {fixed:.6f} = {ap*LN2:.6f} + {fixed:.6f} = {A:.6f}")

print(f"""
  ★ The ONLY free parameter in the curvature A is aₚ (per sector).
    If we can derive aₚ from quantum numbers, A is FULLY determined.

  Can aₚ = f(I₃, Nₒ)?
    ℓ: I₃ = −1/2, Nₒ = 1  → aₚ = 1.00
    u: I₃ = +1/2, Nₒ = 3  → aₚ = 1.75
    d: I₃ = −1/2, Nₒ = 3  → aₚ = −1.00
""")

# Fit aₚ = α·(2I₃) + β·Nₒ + γ
I3_vals = [-0.5, 0.5, -0.5]
Nc_vals = [1, 3, 3]
ap_vals = [1.0, 1.75, -1.0]

# Solve 3×3 system
M = np.array([[2*I3_vals[i], Nc_vals[i], 1] for i in range(3)])
rhs = np.array(ap_vals)
alpha, beta, gamma = np.linalg.solve(M, rhs)

print(f"  Linear fit: aₚ = {alpha:.4f}·(2I₃) + {beta:.4f}·Nₒ + {gamma:.4f}")
print()

# Check with fractions
from fractions import Fraction
alpha_f = Fraction(alpha).limit_denominator(16)
beta_f = Fraction(beta).limit_denominator(16)
gamma_f = Fraction(gamma).limit_denominator(16)
print(f"  As fractions: aₚ = ({alpha_f})·(2I₃) + ({beta_f})·Nₒ + ({gamma_f})")
print()

for name, I3, Nc, ap_true in [("ℓ", -0.5, 1, 1.0),
                                 ("u", 0.5, 3, 1.75),
                                 ("d", -0.5, 3, -1.0)]:
    ap_pred = alpha * (2*I3) + beta * Nc + gamma
    print(f"  {name}: aₚ = {alpha_f}·({2*I3:+.0f}) + ({beta_f})·{Nc} + ({gamma_f}) "
          f"= {ap_pred:.4f}  [actual: {ap_true:.4f}]  ✓")

print(f"""
  This works but has 3 parameters for 3 data points (trivially exact).
  
  HOWEVER: the coefficients {alpha_f}, {beta_f}, {gamma_f} are rational.
  Let's check if they relate to our building blocks:
    {alpha_f} = ?
    {beta_f} = ?  
    {gamma_f} = ?
""")

# Alternative parametrizations
# Try aₚ as function of Q (electric charge)
Q_vals = [-1, 2/3, -1/3]  # e, u, d
print("  Alternative: aₚ vs electric charge Q:")
for name, Q, ap in zip(["ℓ","u","d"], Q_vals, ap_vals):
    print(f"    {name}: Q = {Q:+.2f}, aₚ = {ap:+.4f}")

# Fit linear in Q
M_Q = np.array([[Q_vals[i], 1] for i in range(3)])
# Overdetermined, use least squares
res_Q = np.linalg.lstsq(M_Q, ap_vals, rcond=None)
print(f"  Linear fit: aₚ = {res_Q[0][0]:.4f}·Q + {res_Q[0][1]:.4f}")
print(f"  Residual norm: {np.sqrt(res_Q[1][0]) if len(res_Q[1]) > 0 else 'N/A'}")

# Try: aₚ = α·(2I₃)·Nₒ + β  (interaction term)
# ℓ: α·(-1)·1 + β = 1   → -α + β = 1
# u: α·(1)·3 + β = 1.75  → 3α + β = 1.75
# From these: 4α = 0.75 → α = 3/16... check d:
# d: α·(-1)·3 + β = ?  →  -3α + β = -3·3/16 + (1 + 3/16) = -9/16 + 19/16 = 10/16 = 5/8
# But aₚ(d) = -1. So 5/8 ≠ -1. Doesn't work.

print(f"""
  Full decomposition with (2I₃) and Nₒ:

  aₚ = {alpha_f}·(2I₃) + ({beta_f})·Nₒ + {gamma_f}

  The (2I₃) term = {alpha_f} could be {alpha_f}
  The Nₒ term = {beta_f} (negative: more color → more negative → lighter gen-2!)
  The constant = {gamma_f}
""")


# ════════════════════════════════════════════════════════════════
# §4  VELOCITY B — The Integration Constant
# ════════════════════════════════════════════════════════════════

print("\n§4. VELOCITY B — What Fixes the First Integration Constant?")
print("-" * 72)

print("""
  B = bₚ·ln2 + b_q·lnR₃ + b_r·ln3

  Velocity in mass-space (= initial log-mass growth rate):
""")

print(f"  {'Sector':<10} {'B (mass v.)':>12} {'v(g=1)':>12} {'v(g=2)':>12} {'v(g=3)':>12}")
print(f"  {'─'*55}")

for sector_name in ["lepton", "up", "down"]:
    A, B, C = mass_ABC[sector_name]
    v1 = 2*A*1 + B
    v2 = 2*A*2 + B
    v3 = 2*A*3 + B
    print(f"  {sector_name:<10} {B:>+12.6f} {v1:>+12.4f} {v2:>+12.4f} {v3:>+12.4f}")

print(f"""
  v(g) = 2A·g + B = generational mass velocity at generation g.
  
  v(g) > 0 means mass INCREASES with generation (always true for g=1,2,3).
  v decreases for ℓ and u (negative A → decelerating!).
  v increases for down (positive A → accelerating!).
""")

# Lattice-level velocity vectors
print("  Lattice velocity vectors b = (bₚ, b_q, b_r):")
print(f"  {'Sector':<10} {'bₚ':>8} {'b_q':>8} {'b_r':>8} {'|b|':>10}")
print(f"  {'─'*40}")

b_vectors = {}
for sector_name in ["lepton", "up", "down"]:
    bp = coefficients[sector_name][0][1]
    bq = coefficients[sector_name][1][1]
    br = coefficients[sector_name][2][1]
    norm_b = np.sqrt(bp**2 + bq**2 + br**2)
    b_vectors[sector_name] = np.array([bp, bq, br])
    print(f"  {sector_name:<10} {bp:>+8.2f} {bq:>+8.2f} {br:>+8.2f} {norm_b:>10.4f}")

print()

# Are the velocity vectors related by a rotation/reflection?
print("  Angles between velocity vectors:")
for s1, s2 in [("lepton","up"), ("lepton","down"), ("up","down")]:
    v1, v2 = b_vectors[s1], b_vectors[s2]
    cos_th = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.degrees(np.arccos(np.clip(cos_th, -1, 1)))
    print(f"    {s1}↔{s2}: {angle:.1f}°")

print()


# ════════════════════════════════════════════════════════════════
# §5  THE MINIMAL INPUT PROBLEM
# ════════════════════════════════════════════════════════════════

print("\n§5. MINIMAL INPUT — How Many Masses Must We Know?")
print("-" * 72)

print(f"""
  SCENARIO 1: A from graph, know gen-1 and gen-2 masses → predict gen-3
  ─────────────────────────────────────────────────────────────────────
  Input:  3 curvatures (A) + 3 gen-1 masses + 3 gen-2 masses = 9
  Output: 3 gen-3 masses (predicted)
  Status: TRIVIALLY EXACT (3 points → unique quadratic)

  SCENARIO 2: A from graph, know gen-1 masses → predict gen-2 AND gen-3
  ─────────────────────────────────────────────────────────────────────
  Input:  3 curvatures + 3 gen-1 masses = 6 parameters
  Output: 6 masses (3 gen-2 + 3 gen-3)
  Status: NEED 3 velocity values B → 3 additional constraints
  
  SCENARIO 3: A from graph, B from regularity → only gen-1 masses needed
  ─────────────────────────────────────────────────────────────────────
  Input:  3 curvatures + 3 velocities + 3 gen-1 masses = 9 (but A,B from graph!)
  Output: 6 masses from gen-1 anchors alone
  Status: IDEAL — but what regularity principle fixes B?

  SCENARIO 4: Everything from graph (no experimental input)
  ─────────────────────────────────────────────────────────
  Input:  A from graph + B from regularity + C from quantum numbers
  Output: ALL 9 masses (relative to mₑ)
  Status: ULTIMATE GOAL
""")

# Show what scenario 2 requires
print("  SCENARIO 2 in practice — B as function of gen-1 mass:")
print(f"  {'Sector':<10} {'A':>10} {'ln(m₁/mₑ)':>12} {'B':>10} {'B + 3A':>10}")
print(f"  {'─'*55}")

for sector_name, particles in [("lepton", ["e","mu","tau"]),
                                 ("up", ["u","c","t"]),
                                 ("down", ["d","s","b"])]:
    A, B, C = mass_ABC[sector_name]
    m1 = MASSES[particles[0]]
    ln_m1 = np.log(m1 / m_e)
    # B = (x(2) - x(1)) - 3A = ln(m₂/m₁) - 3A
    print(f"  {sector_name:<10} {A:>+10.6f} {ln_m1:>+12.6f} {B:>+10.6f} {B+3*A:>+10.6f}")

print(f"""
  B + 3A = v(1) + A = velocity at g=1 plus half the acceleration.
  This equals ln(m₂/m₁) = logarithmic mass jump from gen-1 to gen-2.

  If we could derive B from A (or from graph invariants alone),
  knowing the gen-1 mass ratio would give ALL masses.
""")


# ════════════════════════════════════════════════════════════════
# §6  SEARCHING FOR VELOCITY REGULARITY
# ════════════════════════════════════════════════════════════════

print("\n§6. VELOCITY REGULARITY — Testing Hypotheses for B")
print("-" * 72)

print("""
  Hypothesis 1: Is B proportional to A?
  (Constant ratio → same "damping" across sectors)
""")

for sector_name in ["lepton", "up", "down"]:
    A, B, C = mass_ABC[sector_name]
    if abs(A) > 1e-10:
        ratio = B / A
        print(f"  {sector_name:<10}: B/A = {B:.6f}/{A:.6f} = {ratio:.4f}")
    else:
        print(f"  {sector_name:<10}: A ≈ 0, ratio undefined")

print()

print("""
  Hypothesis 2: Is the arc-length per generation constant?
  Arc length = ∫|γ'(g)|dg in lattice space
""")

for sector_name in ["lepton", "up", "down"]:
    names = SECTORS[sector_name]
    total_arc = 0
    arcs = []
    for g in range(1, 3):
        p1 = np.array(COORDS[names[g-1]], dtype=float)
        p2 = np.array(COORDS[names[g]], dtype=float)
        arc = np.linalg.norm(p2 - p1)
        arcs.append(arc)
        total_arc += arc
    ratio_arc = arcs[1] / arcs[0] if arcs[0] != 0 else float('inf')
    print(f"  {sector_name:<10}: arc₁₂ = {arcs[0]:.4f}, arc₂₃ = {arcs[1]:.4f}, "
          f"ratio = {ratio_arc:.4f}")

print()

print("""
  Hypothesis 3: Minimal curvature energy (∫κ²ds)?
  For a quadratic Bézier, curvature is related to |a|.
""")

for sector_name in ["lepton", "up", "down"]:
    ap = coefficients[sector_name][0][0]
    aq = coefficients[sector_name][1][0]
    ar = coefficients[sector_name][2][0]
    a_vec = np.array([ap, aq, ar])
    curv_energy = np.sum(a_vec**2)
    print(f"  {sector_name:<10}: |a|² = {curv_energy:.4f}, |a| = {np.sqrt(curv_energy):.4f}")

print()

print("""
  Hypothesis 4: Does B live on the lattice?
  Check if B-components are expressible in the basis {ln2, lnR₃, ln3}.
""")

# B = bₚ·ln2 + b_q·lnR₃ + b_r·ln3
# We already know the b components are rational (half-integers)
print("  b-coefficients (should be rational if lattice-native):")
print(f"  {'Sector':<10} {'bₚ':>10} {'b_q':>10} {'b_r':>10}")
print(f"  {'─'*42}")
for sector_name in ["lepton", "up", "down"]:
    bp = coefficients[sector_name][0][1]
    bq = coefficients[sector_name][1][1]
    br = coefficients[sector_name][2][1]
    print(f"  {sector_name:<10} {bp:>+10.2f} {bq:>+10.2f} {br:>+10.2f}")

# Check: are all b values multiples of 0.25?
print()
all_b = []
for sector_name in ["lepton", "up", "down"]:
    for idx in range(3):
        all_b.append(coefficients[sector_name][idx][1])

min_step = min(abs(b) for b in all_b if abs(b) > 0)
print(f"  Smallest |b| = {min_step}")
print(f"  All b as multiples of 0.25:")
for sector_name in ["lepton", "up", "down"]:
    vals = [coefficients[sector_name][idx][1] for idx in range(3)]
    mults = [v / 0.25 for v in vals]
    print(f"  {sector_name:<10}: {mults[0]:>+6.0f}/4, {mults[1]:>+6.0f}/4, {mults[2]:>+6.0f}/4")


# ════════════════════════════════════════════════════════════════
# §7  VELOCITY FROM CURVATURE — Structural Relations
# ════════════════════════════════════════════════════════════════

print("\n\n§7. VELOCITY FROM CURVATURE — Are b and a Related?")
print("-" * 72)

print("""
  The curvature (a) and velocity (b) vectors per sector:
  If the lattice imposes a relationship b = f(a), we can derive B from A.
""")

print(f"  {'Sector':<10} {'aₚ':>6} {'a_q':>6} {'a_r':>6} │ {'bₚ':>8} {'b_q':>8} {'b_r':>8} │ {'b/a ratio':>10}")
print(f"  {'─'*75}")

for sector_name in ["lepton", "up", "down"]:
    a_vec = np.array([coefficients[sector_name][i][0] for i in range(3)])
    b_vec = np.array([coefficients[sector_name][i][1] for i in range(3)])
    ratios = []
    for i in range(3):
        if abs(a_vec[i]) > 1e-10:
            ratios.append(f"{b_vec[i]/a_vec[i]:+.2f}")
        else:
            ratios.append("  ∞")
    print(f"  {sector_name:<10} {a_vec[0]:>+6.2f} {a_vec[1]:>+6.2f} {a_vec[2]:>+6.2f} │ "
          f"{b_vec[0]:>+8.2f} {b_vec[1]:>+8.2f} {b_vec[2]:>+8.2f} │ "
          f"{ratios[0]:>5} {ratios[1]:>5} {ratios[2]:>5}")

print()

# Check: b + k·a = something simple?
print("  Testing b + k·a for small integers k:")
for k in range(-6, 7):
    results = {}
    for sector_name in ["lepton", "up", "down"]:
        a_vec = np.array([coefficients[sector_name][i][0] for i in range(3)])
        b_vec = np.array([coefficients[sector_name][i][1] for i in range(3)])
        combo = b_vec + k * a_vec
        results[sector_name] = combo
    # Check if any combination is shared or simple
    print(f"  k={k:+2d}: ℓ=({results['lepton'][0]:+6.2f},{results['lepton'][1]:+6.2f},{results['lepton'][2]:+6.2f}) "
          f" u=({results['up'][0]:+6.2f},{results['up'][1]:+6.2f},{results['up'][2]:+6.2f}) "
          f" d=({results['down'][0]:+6.2f},{results['down'][1]:+6.2f},{results['down'][2]:+6.2f})")

print()

# Look for the tangent at the "crossing point" g=3
# The tangent v(3) = 2a·3 + b = 6a + b
print("  Tangent at g=3 (unification point): v(3) = 6a + b")
for sector_name in ["lepton", "up", "down"]:
    a_vec = np.array([coefficients[sector_name][i][0] for i in range(3)])
    b_vec = np.array([coefficients[sector_name][i][1] for i in range(3)])
    v3 = 6 * a_vec + b_vec
    print(f"  {sector_name:<10}: v(3) = ({v3[0]:+6.2f}, {v3[1]:+6.2f}, {v3[2]:+6.2f})")

print()

# Tangent at g=0 (extrapolated "origin")
print("  Tangent at g=0 (extrapolated pre-generation): v(0) = b")
for sector_name in ["lepton", "up", "down"]:
    b_vec = np.array([coefficients[sector_name][i][1] for i in range(3)])
    print(f"  {sector_name:<10}: v(0) = ({b_vec[0]:+6.2f}, {b_vec[1]:+6.2f}, {b_vec[2]:+6.2f})")

print()

# Position at g=0 (extrapolated)
print("  Position at g=0 (extrapolated): γ(0) = c")
for sector_name in ["lepton", "up", "down"]:
    c_vec = np.array([coefficients[sector_name][i][2] for i in range(3)])
    print(f"  {sector_name:<10}: γ(0) = ({c_vec[0]:+6.2f}, {c_vec[1]:+6.2f}, {c_vec[2]:+6.2f})")


# ════════════════════════════════════════════════════════════════
# §8  VELOCITY CONSTRAINT FROM UNIFICATION
# ════════════════════════════════════════════════════════════════

print("\n\n§8. VELOCITY CONSTRAINT — Unification at g=3 Restricts B")
print("-" * 72)

print("""
  We know q_ℓ(3) = q_u(3) = −7 (exact unification at gen 3).
  
  This constrains the q-velocity:
    q(g) = a_q·g² + b_q·g + c_q
    q(3) = 9·a_q + 3·b_q + c_q = −7
    
  For ℓ: 9(0.5) + 3(-5.5) + 5.0 = 4.5 - 16.5 + 5.0 = -7.0 ✓
  For u: 9(0.5) + 3(-2.5) +(-4.0) = 4.5 - 7.5 - 4.0 = -7.0 ✓
  
  These share a_q = 0.5 but have DIFFERENT b_q: -5.5 vs -2.5.
  The difference b_q(u) − b_q(ℓ) = 3.0 is compensated by 
  the difference c_q(u) − c_q(ℓ) = -9.0.
  
  Constraint: 3·Δb_q + Δc_q = 0 → Δc_q = −3·Δb_q
  (The offset absorbs 3× the velocity difference to meet at g=3.)
""")

db_q = coefficients['up'][1][1] - coefficients['lepton'][1][1]
dc_q = coefficients['up'][1][2] - coefficients['lepton'][1][2]
print(f"  Δb_q = {db_q:+.1f}, Δc_q = {dc_q:+.1f}, −3·Δb_q = {-3*db_q:+.1f}  ✓")

print(f"""
  ★ The q-unification at g=3 couples velocity and offset:
    "faster start" (larger b_q) requires "deeper offset" (more negative c_q)
    to arrive at the same point at generation 3.
    
  This is a BOUNDARY VALUE problem, not an initial value problem!
  The sectors are like beams anchored at both ends (g=1 and g=3),
  not free projectiles.
""")


# ════════════════════════════════════════════════════════════════
# §9  RETHINKING: ENDPOINT CONSTRAINTS
# ════════════════════════════════════════════════════════════════

print("\n§9. RETHINKING — From IVP to BVP (Boundary Value Problem)")
print("-" * 72)

print("""
  We've been thinking "initial value": know acceleration + initial conditions.
  But the LATTICE suggests "boundary value": know endpoints + curvature.

  For each sector, the curve passes through 3 lattice points.
  These are NOT arbitrary — they are constrained by:
  
  1. INTEGER/HALF-INTEGER coordinates (lattice quantization)
  2. Curvature conservation |a_q + a_r| = 1
  3. q-unification at g=3 (at least ℓ↔u)
  4. Electron at origin (gen-1 lepton)
  
  So the CORRECT formulation is:
  
  ┌─────────────────────────────────────────────────────────┐
  │  Given: lattice {ln2, lnR₃, ln3}                       │
  │         curvature law |a_q + a_r| = 1                   │
  │         quantization: coordinates ∈ ½ℤ                  │
  │         electron = (0,0,0)                              │
  │         unification q_ℓ(3) = q_u(3)                     │
  │                                                         │
  │  Find: all lattice positions that satisfy these AND      │
  │        produce physical mass ordering m₁ < m₂ < m₃      │
  └─────────────────────────────────────────────────────────┘
  
  This is a DISCRETE optimization: find integer points on
  quadratic Bézier curves with fixed curvature constraints!
""")

# How many candidate lattice points are there?
# gen-1 positions: constrained by quantum numbers (roughly)
# gen-2 positions: constrained by curvature + gen-1 + gen-3
# gen-3 positions: constrained by q-unification at g=3

print("  LATTICE QUANTIZATION CHECK:")
print(f"  Are all coordinates in ½ℤ?")
print()
all_half_int = True
for pname, (p, q, r) in sorted(COORDS.items()):
    is_half = all(abs(2*x - round(2*x)) < 1e-10 for x in [p, q, r])
    mark = "✓" if is_half else "✗"
    all_half_int = all_half_int and is_half
    print(f"  {pname:<5}: ({p:+5.1f}, {q:+5.1f}, {r:+5.1f})  {mark}")

print(f"\n  All ∈ ½ℤ: {'YES' if all_half_int else 'NO'}")

# Actually check: all coords are integers or k/2
print("\n  Finer check — coordinate types:")
types = set()
for pname, (p, q, r) in COORDS.items():
    for x in [p, q, r]:
        if abs(x - round(x)) < 1e-10:
            types.add("integer")
        elif abs(2*x - round(2*x)) < 1e-10:
            types.add("half-integer")
        elif abs(4*x - round(4*x)) < 1e-10:
            types.add("quarter-integer")
        else:
            types.add(f"other({x})")
print(f"  Types found: {types}")


# ════════════════════════════════════════════════════════════════
# §10  DISCRETE SEARCH — How Many Solutions Exist?
# ════════════════════════════════════════════════════════════════

print("\n\n§10. DISCRETE SEARCH — Counting Compatible Lattice Configurations")
print("-" * 72)

print("""
  Setup: Fix the curvatures (a_p, a_q, a_r) per sector.
  Then for each sector, the gen-2 particle is determined by gen-1 and gen-3:
    γ(1) + γ(3) − 2·γ(2) = 2a  ⇒  coord_2 = (coord_1 + coord_3)/2 − a
    
  So the FREE choices are the gen-1 and gen-3 positions (gen-2 follows).
  
  Constraints:
    - All coordinates in ½ℤ (after applying curvature)
    - Mass ordering: m₁ < m₂ < m₃ within each sector
    - Lepton gen-1 = (0,0,0)
    - q_ℓ(3) = q_u(3) (unification)
    
  Search space: gen-1 and gen-3 positions for each sector.
  Let's enumerate all valid configurations in a bounded region
  and see how constrained the system really is.
""")

# For generation 3 (tau, t, b), we need q = −7 for ℓ and u (unification)
# Lepton gen-1 fixed at (0,0,0)
# Lepton gen-3 = (p₃, −7, r₃) for integer/half-int p₃, r₃

# The mass at generation g is: m = mₑ · exp(p·ln2 + q·lnR₃ + r·ln3)
def mass_from_coords(p, q, r):
    return m_e * np.exp(p * LN2 + q * LNR3 + r * LN3)

# For leptons: gen-1 = (0,0,0), gen-3 fixed by search
# gen-2 = (gen1 + gen3)/2 + a/2
a_lepton = np.array([1.0, 0.5, -1.5])
a_up = np.array([1.75, 0.5, -1.5])
a_down = np.array([-1.0, 0.0, 1.0])

# Enumerate lepton solutions
P_RANGE = np.arange(-3, 10.5, 0.5)  # half-integer range for p
R_RANGE = np.arange(-3, 10.5, 0.5)  # half-integer range for r

lepton_solutions = []
for p3 in P_RANGE:
    for r3 in R_RANGE:
        gen1 = np.array([0, 0, 0], dtype=float)
        gen3 = np.array([p3, -7, r3], dtype=float)
        gen2 = (gen1 + gen3) / 2 - a_lepton  # γ(2) = (γ(1)+γ(3))/2 - a

        # Check gen2 is half-integer
        if not all(abs(2*x - round(2*x)) < 1e-10 for x in gen2):
            continue

        # Check mass ordering
        m1 = mass_from_coords(*gen1)
        m2 = mass_from_coords(*gen2)
        m3 = mass_from_coords(*gen3)

        if m1 < m2 < m3:
            lepton_solutions.append({
                'gen1': tuple(gen1), 'gen2': tuple(gen2), 'gen3': tuple(gen3),
                'masses': (m1, m2, m3),
                'p3': p3, 'r3': r3
            })

print(f"  LEPTON sector: gen-1 = (0,0,0), gen-3 q = −7")
print(f"  Curvature: a = ({a_lepton[0]:+.1f}, {a_lepton[1]:+.1f}, {a_lepton[2]:+.1f})")
print(f"  Search range: p₃ ∈ [{P_RANGE[0]}, {P_RANGE[-1]}], r₃ ∈ [{R_RANGE[0]}, {R_RANGE[-1]}]")
print(f"  Valid configurations (m₁ < m₂ < m₃, half-integer): {len(lepton_solutions)}")
print()

# Show the actual solution
actual_lepton = None
for sol in lepton_solutions:
    if abs(sol['p3'] - 1.0) < 0.01 and abs(sol['r3'] - 3.0) < 0.01:
        actual_lepton = sol
        break

if actual_lepton:
    print(f"  ★ ACTUAL lepton configuration:")
    print(f"    gen-1: {actual_lepton['gen1']} → m = {actual_lepton['masses'][0]:.3f} MeV")
    print(f"    gen-2: {actual_lepton['gen2']} → m = {actual_lepton['masses'][1]:.3f} MeV")
    print(f"    gen-3: {actual_lepton['gen3']} → m = {actual_lepton['masses'][2]:.3f} MeV")

print()

# How many match τ mass within 10%?
tau_exp = MASSES['tau']
close_to_tau = [s for s in lepton_solutions
                if abs(s['masses'][2] - tau_exp) / tau_exp < 0.1]
print(f"  Solutions with m₃ within 10% of m_τ ({tau_exp:.1f} MeV): {len(close_to_tau)}")

# And matching μ too?
mu_exp = MASSES['mu']
close_to_both = [s for s in lepton_solutions
                 if abs(s['masses'][1] - mu_exp) / mu_exp < 0.1
                 and abs(s['masses'][2] - tau_exp) / tau_exp < 0.1]
print(f"  Solutions with BOTH m₂ ≈ m_μ AND m₃ ≈ m_τ (10%): {len(close_to_both)}")

if close_to_both:
    print("\n  These are:")
    for sol in close_to_both:
        print(f"    gen-3 = ({sol['p3']:+.1f}, −7, {sol['r3']:+.1f}), "
              f"m₂ = {sol['masses'][1]:.1f} MeV, m₃ = {sol['masses'][2]:.1f} MeV")

# Total count analysis
print(f"""
  ★ CONSTRAINT POWER:
    Before constraints: {len(P_RANGE)} × {len(R_RANGE)} = {len(P_RANGE)*len(R_RANGE)} lattice points
    After mass ordering: {len(lepton_solutions)} valid
    After mass matching (10%): {len(close_to_both)} valid
    
    The lattice + curvature law reduces the solution space dramatically.
""")


# ════════════════════════════════════════════════════════════════
# §11  VELOCITY DECOMPOSITION — b in terms of graph invariants
# ════════════════════════════════════════════════════════════════

print("\n§11. VELOCITY DECOMPOSITION — Can We Express b?")
print("-" * 72)

print("""
  Like aₚ = f(I₃, Nc), can we express bₚ, b_q, b_r in terms
  of quantum numbers and/or graph invariants?
  
  The b-coefficients:
    Lepton: (−3.5, −5.5, +7.5)
    Up:     (−3.25, −2.5, +8.5)
    Down:   (+5.0, +1.0, −1.0)
    
  Differences (relative to lepton):
    Δb(u−ℓ): (+0.25, +3.0, +1.0)
    Δb(d−ℓ): (+8.5, +6.5, −8.5)
""")

# Is Δb(u−ℓ) related to Δa(u−ℓ)?
da_ul = np.array([1.75-1.0, 0.5-0.5, -1.5-(-1.5)])  # = (0.75, 0, 0)
db_ul = np.array([-3.25-(-3.5), -2.5-(-5.5), 8.5-7.5])  # = (0.25, 3.0, 1.0)

print(f"  Δa(u−ℓ) = ({da_ul[0]:+.2f}, {da_ul[1]:+.2f}, {da_ul[2]:+.2f})")
print(f"  Δb(u−ℓ) = ({db_ul[0]:+.2f}, {db_ul[1]:+.2f}, {db_ul[2]:+.2f})")
print()

# The q-unification constraint: 3·Δb_q + Δc_q = 0
# Also for r: if r_ℓ(3) ≈ r_u(3), then 3·Δb_r + Δc_r = 0 as well?
dr_at_3_lu = (coefficients['up'][2][0]*9 + coefficients['up'][2][1]*3 + coefficients['up'][2][2]) - \
             (coefficients['lepton'][2][0]*9 + coefficients['lepton'][2][1]*3 + coefficients['lepton'][2][2])
print(f"  r_u(3) − r_ℓ(3) = {dr_at_3_lu:+.1f}")
print(f"  (ℓ and u are close at gen 3 in r too — gap of only 1)")
print()

# What if we impose EXACT r-unification at g=3?
# Then a_r shared → Δa_r = 0, so:
# r_u(3) = r_ℓ(3) requires 3·Δb_r + Δc_r = 0
# Currently: 3(1.0) + (-8-(-6)) = 3 - 2 = 1 ≠ 0
# So r-unification is APPROXIMATE (gap = 1), not exact like q.
print("  EXACT unification test at g=3:")
for idx, comp in enumerate(['p', 'q', 'r']):
    vals = {}
    for sector_name in ["lepton", "up", "down"]:
        a, b, c = coefficients[sector_name][idx]
        vals[sector_name] = a*9 + b*3 + c
    print(f"  {comp}: ℓ={vals['lepton']:+.1f}, u={vals['up']:+.1f}, d={vals['down']:+.1f}, "
          f"Δ(u−ℓ)={vals['up']-vals['lepton']:+.1f}, "
          f"Δ(d−ℓ)={vals['down']-vals['lepton']:+.1f}")

print(f"""
  At g = 3: q-unification is EXACT (Δ=0), r-gap = 1, p-gap very large.
  
  The q-unification is a HARD constraint: 3·Δb_q + Δc_q = 0 (ℓ↔u).
  The r near-unification is SOFT: gap = 1.
  
  If we ALSO imposed r-unification, Δb_r would be further constrained.
  Currently Δb_r(u−ℓ) = 1.0, Δc_r(u−ℓ) = −2.0, giving gap = 3(1)+(-2) = 1.
  For exact: need Δc_r = −3·Δb_r → Δb_r = Δb_r, Δc_r = −3·Δb_r
  With current Δb_r = 1: need Δc_r = −3 (currently −2).
  → Missing 1 unit of offset.

  ★ The r-gap of exactly 1 at gen 3 may be the TOPOLOGICAL UNIT
    that prevents full C¹ continuity — the "quantum" of spectral gap.
""")


# ════════════════════════════════════════════════════════════════
# §12  SUMMARY — The Velocity Problem
# ════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("  §12. SUMMARY — The Status of Mass Derivation")
print("=" * 72)

print(f"""
  WHAT WE HAVE (from graph topology alone):
  ──────────────────────────────────────────
  ✓ 3D lattice basis: {{ln2, lnR₃, ln3}}
  ✓ Curvature law: |a_q + a_r| = 1
  ✓ Isocurvature: a_q(ℓ) = a_q(u), a_r(ℓ) = a_r(u)
  ✓ Down linearity: a_q(d) = 0
  ✓ q-unification: q_ℓ(3) = q_u(3) = −7
  ✓ Electron = origin
  ✓ Quantization: all coords ∈ ½ℤ

  WHAT THIS DETERMINES:
  ─────────────────────
  ✓ Curvature vector a per sector (given aₚ = f(I₃,Nc))
  ✓ Relationship: gen-2 = (gen-1 + gen-3)/2 − a  (Bézier midpoint)
  ✓ BVP coupling: Δc_q = −3·Δb_q between ℓ and u

  WHAT REMAINS FREE:
  ──────────────────
  ? aₚ per sector: expressible as {alpha_f}·(2I₃)+({beta_f})·Nc+{gamma_f}
    but not yet derived from first principles
  ? Velocity b per sector per component (9 params, reduced by constraints)
  ? Gen-1 positions for up/down sectors (4 coords, after e=origin)
  ? Gen-3 positions (constrained by unification for q)

  THE CORRECT PICTURE:
  ────────────────────
  This is NOT an initial value problem (know start → predict end).
  This is a BOUNDARY VALUE problem on a discrete lattice:
  
  Find half-integer coordinates that satisfy:
    1. Curvature constraints (from gauge structure)
    2. Mass ordering m₁ < m₂ < m₃ (from stability)
    3. q-unification at gen 3 (from spectral convergence)
    4. Bézier midpoint relation (quadratic structure)
    
  The solution space is HIGHLY constrained:
  for leptons, only {len(close_to_both)} configuration(s) match
  the observed masses within 10%.
  
  NEXT STEP: Run the discrete search for ALL THREE sectors
  simultaneously, with inter-sector constraints from CKM mixing
  (the C⁰ crossing points). This should further reduce solutions.
""")


# Close output
_outfile.close()
sys.stdout = _orig_stdout
print(f"\nOutput saved to: {OUTPUT_PATH}")
