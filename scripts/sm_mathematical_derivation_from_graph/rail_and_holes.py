#!/usr/bin/env python3
"""
rail_and_holes.py — The marble-on-a-rail picture
==================================================
The parabola x(g) = ag² + bg + c traces a continuous curve ("rail") in
3D lattice space (ln2, lnR₃, ln3). Physical particles sit at g=1,2,3.

Questions:
  1. How close does the rail pass to lattice points at NON-integer g?
  2. Does the lattice create "potential wells" at g=1,2,3?
  3. What is the "depth" of each hole?
"""

import numpy as np
from fractions import Fraction as F
from collections import OrderedDict

# ── Data ──
PARTICLES = OrderedDict([
    ("e",   (F(0),    F(0),    F(0))),
    ("mu",  (F(-1,2), F(-4),   F(3))),
    ("tau", (F(1),    F(-7),   F(3))),
    ("u",   (F(-3,2), F(-6),   F(-1))),
    ("c",   (F(1,2),  F(-7),   F(3))),
    ("t",   (F(6),    F(-7),   F(4))),
    ("d",   (F(-1,2), F(-8),   F(-2))),
    ("s",   (F(3,2),  F(-7),   F(0))),
    ("b",   (F(3,2),  F(-6),   F(4))),
])

SECTORS = OrderedDict([
    ("lepton", ["e", "mu", "tau"]),
    ("up",     ["u", "c", "t"]),
    ("down",   ["d", "s", "b"]),
])

def fit_parabola(x1, x2, x3):
    x1, x2, x3 = F(x1), F(x2), F(x3)
    a = (x1 - 2*x2 + x3) / 2
    b = (-5*x1 + 8*x2 - 3*x3) / 2
    c = 3*x1 - 3*x2 + x3
    return a, b, c

# Build coefficients
coeffs = {}  # sector -> [(a_p,b_p,c_p), (a_q,b_q,c_q), (a_r,b_r,c_r)]
for sname, names in SECTORS.items():
    coords = [PARTICLES[n] for n in names]
    sector_coeffs = []
    for k in range(3):
        a, b, c = fit_parabola(coords[0][k], coords[1][k], coords[2][k])
        sector_coeffs.append((float(a), float(b), float(c)))
    coeffs[sname] = sector_coeffs

def x_vec(sector, g):
    """Return 3D position on the rail at continuous parameter g."""
    return np.array([a*g**2 + b*g + c for a, b, c in coeffs[sector]])

def v_vec(sector, g):
    """Return 3D velocity on the rail at continuous parameter g."""
    return np.array([2*a*g + b for a, b, c in coeffs[sector]])

# ════════════════════════════════════════════════════════════════
# §1. THE RAIL IN 3D: continuous trajectory
# ════════════════════════════════════════════════════════════════
print("=" * 75)
print("  §1. THE RAIL: continuous parabolic trajectory in lattice space")
print("=" * 75)

for sname in SECTORS:
    names = SECTORS[sname]
    print(f"\n  {sname.upper()} sector:")
    print(f"  {'g':>5}  {'p':>8}  {'q':>8}  {'r':>8}  {'|x|':>8}  {'|v|':>8}")
    print("  " + "-" * 55)
    for g_val in np.arange(0.5, 3.51, 0.25):
        x = x_vec(sname, g_val)
        v = v_vec(sname, g_val)
        marker = " <-- GEN" if g_val in [1.0, 2.0, 3.0] else ""
        print(f"  {g_val:5.2f}  {x[0]:8.3f}  {x[1]:8.3f}  {x[2]:8.3f}  "
              f"{np.linalg.norm(x):8.3f}  {np.linalg.norm(v):8.3f}{marker}")

# ════════════════════════════════════════════════════════════════
# §2. DISTANCE TO NEAREST LATTICE POINT along the rail
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 75)
print("  §2. LATTICE PROXIMITY: distance to nearest half-integer point")
print("=" * 75)

def dist_to_lattice(x):
    """Distance from x to nearest point with half-integer coordinates."""
    # Our lattice has coords in multiples of 1/2 (quarter-integer denominators)
    # but actual particle coords are multiples of 1/2
    nearest = np.round(2 * x) / 2  # snap to nearest half-integer
    return np.linalg.norm(x - nearest), nearest

print("""
  The lattice has half-integer spacing (all particle coords are in Z/2).
  We measure how close the continuous rail passes to lattice nodes.
""")

for sname in SECTORS:
    print(f"\n  {sname.upper()} sector:")
    print(f"  {'g':>5}  {'d_lattice':>10}  {'nearest (p,q,r)':>25}  {'actual?':>8}")
    print("  " + "-" * 60)

    for g_val in np.arange(0.5, 3.51, 0.1):
        x = x_vec(sname, g_val)
        d, nearest = dist_to_lattice(x)
        is_gen = ""
        if abs(g_val - round(g_val)) < 0.01 and round(g_val) in [1, 2, 3]:
            is_gen = f"GEN {round(g_val)}"
        if d < 0.5 or is_gen:  # only show close approaches
            n_str = f"({nearest[0]:5.1f},{nearest[1]:5.1f},{nearest[2]:5.1f})"
            print(f"  {g_val:5.2f}  {d:10.4f}  {n_str:>25}  {is_gen:>8}")

# ════════════════════════════════════════════════════════════════
# §3. THE "POTENTIAL WELL" — lattice distance as function of g
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 75)
print("  §3. POTENTIAL WELL: V(g) = d²(x(g), lattice)")
print("=" * 75)

print("""
  If V(g) = ||x(g) - nearest_lattice(x(g))||², then minima of V(g) are
  the "holes" where the marble settles. Do g=1,2,3 minimize V(g)?
""")

for sname in SECTORS:
    print(f"\n  {sname.upper()} sector:")
    g_fine = np.arange(0.0, 4.01, 0.01)
    V_vals = []
    for g_val in g_fine:
        x = x_vec(sname, g_val)
        d, _ = dist_to_lattice(x)
        V_vals.append(d**2)

    V_vals = np.array(V_vals)

    # Find all local minima
    minima = []
    for i in range(1, len(V_vals) - 1):
        if V_vals[i] < V_vals[i-1] and V_vals[i] < V_vals[i+1]:
            minima.append((g_fine[i], V_vals[i]))

    print(f"  Local minima of V(g) = d²:")
    for g_m, v_m in minima:
        is_gen = ""
        for gi in [1, 2, 3]:
            if abs(g_m - gi) < 0.05:
                is_gen = f"  <-- GEN {gi}"
        print(f"    g = {g_m:.2f},  V = {v_m:.6f}{is_gen}")

    # Values at exact g=1,2,3
    print(f"  Exact values at generations:")
    for gi in [1, 2, 3]:
        x = x_vec(sname, gi)
        d, nearest = dist_to_lattice(x)
        print(f"    g={gi}: V = {d**2:.6f}, d = {d:.6f}, "
              f"x = ({x[0]:.2f},{x[1]:.2f},{x[2]:.2f})")

# ════════════════════════════════════════════════════════════════
# §4. FRACTIONAL PART ANALYSIS — why g=integer is special
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 75)
print("  §4. FRACTIONAL PARTS: integrality at g=1,2,3")
print("=" * 75)

print("""
  For x(g) = ag² + bg + c with a,b,c in Z/4, when is x(g) in Z/2?
  
  x(g) = ag² + bg + c.  If g is integer:
    ag² in Z/4 (since a in Z/4 and g² in Z)
    bg  in Z/4 (since b in Z/4 and g in Z)
    c   in Z/4
    Sum in Z/4... but actual coords are in Z/2, so we need tighter.
""")

# Check: for which rational g in [0,4] is x(g) in (Z/2)^3 for ALL sectors?
print("  Scanning g in Q: for which g is x(g) in (Z/2)^3 for all sectors?")
print()

from fractions import Fraction

# Exact coefficients
exact_coeffs = {}
for sname, names in SECTORS.items():
    coords = [PARTICLES[n] for n in names]
    sector_ex = []
    for k in range(3):
        a, b, c = fit_parabola(coords[0][k], coords[1][k], coords[2][k])
        sector_ex.append((a, b, c))
    exact_coeffs[sname] = sector_ex

def is_half_integer(x):
    """Check if x is in Z/2."""
    return (2 * x).denominator == 1

count = 0
valid_g = []
# Scan g = n/d for d up to 12
for denom in range(1, 13):
    for numer in range(0, 4 * denom + 1):
        g = Fraction(numer, denom)
        if g > 4:
            break
        all_half_int = True
        for sname in SECTORS:
            for a, b, c in exact_coeffs[sname]:
                val = a * g**2 + b * g + c
                if not is_half_integer(val):
                    all_half_int = False
                    break
            if not all_half_int:
                break
        if all_half_int:
            valid_g.append(g)

print(f"  g values in [0,4] with denominator <= 12 where ALL 9 coords are in Z/2:")
print()
for g in sorted(set(valid_g)):
    coords_str = []
    for sname in SECTORS:
        x = [a * g**2 + b * g + c for a, b, c in exact_coeffs[sname]]
        coords_str.append(f"({x[0]},{x[1]},{x[2]})")
    is_phys = " <-- PHYSICAL" if g in [1, 2, 3] else ""
    print(f"    g = {str(g):>5} = {float(g):5.3f}  "
          f"l={coords_str[0]:>20}  u={coords_str[1]:>20}  d={coords_str[2]:>20}{is_phys}")

print(f"\n  Total: {len(set(valid_g))} valid g values found.")

# ════════════════════════════════════════════════════════════════
# §5. THE HOLES: what makes g=1,2,3 special vs other valid g?
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 75)
print("  §5. WHY EXACTLY 3 HOLES? Integer vs half-integer g")
print("=" * 75)

print("""
  Among the valid g values, g=1,2,3 are the INTEGER ones in the physical
  range. But are there special properties of integers?
""")

# For integer g, check which coords land on Z (not just Z/2)
print("  Integrality level at each valid g:")
print(f"  {'g':>5}  {'#Z coords':>10}  {'#Z/2 only':>10}  {'detail'}")
print("  " + "-" * 60)

for g in sorted(set(valid_g)):
    if g < 0 or g > 4:
        continue
    n_int = 0
    n_half = 0
    for sname in SECTORS:
        for a, b, c in exact_coeffs[sname]:
            val = a * g**2 + b * g + c
            if val.denominator == 1:
                n_int += 1
            else:
                n_half += 1
    tag = ""
    if g in [Fraction(1), Fraction(2), Fraction(3)]:
        tag = " <-- PHYSICAL"
    if n_int == 9:
        tag += " (ALL integer)"
    print(f"  {str(g):>5}  {n_int:>10}  {n_half:>10}  {tag}")

# ════════════════════════════════════════════════════════════════
# §6. SPEED ON THE RAIL at each generation
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 75)
print("  §6. MARBLE SPEED: |v(g)| at each generation")
print("=" * 75)

print("""
  Speed = |dx/dg| tells us how fast the marble moves through each hole.
  Slow = deep hole (long dwell time), Fast = shallow hole.
""")

for sname in SECTORS:
    names = SECTORS[sname]
    print(f"\n  {sname.upper()} ({', '.join(names)}):")
    speeds = []
    for gi in [1, 2, 3]:
        v = v_vec(sname, gi)
        speed = np.linalg.norm(v)
        speeds.append(speed)
        print(f"    g={gi} ({names[gi-1]:>3}): |v| = {speed:.4f}, "
              f"v = ({v[0]:+.2f}, {v[1]:+.2f}, {v[2]:+.2f})")
    print(f"    Speed ratio gen1/gen3 = {speeds[0]/speeds[2]:.3f}")
    print(f"    Slowest: gen {np.argmin(speeds)+1} ({names[np.argmin(speeds)]})")

# ════════════════════════════════════════════════════════════════
# §7. SUMMARY
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 75)
print("  §7. SUMMARY: The Marble-on-a-Rail Picture")
print("=" * 75)
print("""
  RAIL: Each sector has a parabolic curve x(g) in 3D lattice space.
        The shape is mostly fixed by quantum numbers (curvature a)
        and flatness (velocity b for 3 coords).

  HOLES: The rail passes through lattice nodes ONLY at specific g values.
         g = 1, 2, 3 are the physical particles (generations).
         The lattice acts as a PERIODIC POTENTIAL that quantizes g.

  MARBLE: Rolls along the rail with speed |v(g)|.
          - Fast through gen 1 (shallow hole, light particle)
          - Slow through gen 2-3 for flat coords (deep hole, heavy particles degenerate)
          - The 6 non-flat delta_b values control the tilt of the rail
            and thus the relative depth of the holes.

  MISSING: The equation of motion for the marble.
           What we know:  S = (8/3)sum(a²) + sum(delta_b²)
           What we need:  V(g) that creates the holes at g = integer
           Candidate:     V(g) ~ d²(x(g), lattice) = sum of fractional parts²

  THE KEY INSIGHT: The lattice IS the potential.
  The marble doesn't need an external force — the discrete structure
  of the Cartan-Weyl graph creates the wells naturally.
""")
