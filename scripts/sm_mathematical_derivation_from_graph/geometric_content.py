#!/usr/bin/env python3
"""
geometric_content.py — What do b=-4a and b=-5a tell us?
========================================================
b = -4a  =>  g* = 2    =>  x(1) = x(3)  =>  gen 1-3 palindrome, zero mean velocity
b = -5a  =>  g* = 5/2  =>  x(2) = x(3)  =>  gen 2-3 degeneracy
"""

from fractions import Fraction as F

# Exact coefficients (a, b) from parameter_reduction.py
data = {
    'l_p': (F(1),      F(-7, 2)),
    'l_q': (F(1, 2),   F(-11, 2)),
    'l_r': (F(-3, 2),  F(15, 2)),    # FLAT
    'u_p': (F(7, 4),   F(-13, 4)),
    'u_q': (F(1, 2),   F(-5, 2)),    # FLAT
    'u_r': (F(-3, 2),  F(17, 2)),
    'd_p': (F(-1),     F(5)),         # FLAT
    'd_q': (F(0),      F(1)),
    'd_r': (F(1),      F(-1)),
}

flat_coords = {'l_r', 'u_q', 'd_p'}

print("=" * 80)
print("  GEOMETRIC CONTENT OF b/a RATIOS")
print("=" * 80)

# ── Table ──
hdr = f"{'coord':>6}  {'a':>7}  {'b':>7}  {'ratio':>7}  {'g*':>6}  {'x1=x3':>7}  {'x2=x3':>7}  {'flat'}"
print(f"\n{hdr}")
print("-" * 80)

for name, (a, b) in data.items():
    c = F(0)
    x1 = a + b + c
    x2 = 4*a + 2*b + c
    x3 = 9*a + 3*b + c

    if a != 0:
        ratio_s = f"{float(b/a):.2f}"
        gstar_s = f"{float(-b/(2*a)):.2f}"
    else:
        ratio_s = "inf"
        gstar_s = "inf"

    eq13 = "YES" if x1 == x3 else "no"
    eq23 = "YES" if x2 == x3 else "no"
    fl = "FLAT" if name in flat_coords else ""

    print(f"{name:>6}  {str(a):>7}  {str(b):>7}  {ratio_s:>7}  {gstar_s:>6}  {eq13:>7}  {eq23:>7}  {fl}")

# ── Flatness verification ──
print("\n── Flatness verification (b + 5a = 0?) ──")
for name in sorted(flat_coords):
    a, b = data[name]
    print(f"  {name}: b + 5a = {b + 5*a}")

# ── Deviation from kinetic optimum ──
print("\n── Deviation delta_b = b + 4a ──")
total_db2 = F(0)
for name, (a, b) in data.items():
    db = b + 4*a
    total_db2 += db**2
    fl = " (flat: delta_b = -a)" if name in flat_coords else ""
    print(f"  {name}: delta_b = {db:>7} = {float(db):>+8.4f}{fl}")
print(f"\n  Sum(delta_b^2) = {total_db2} = {float(total_db2):.4f}")

# ── Mean velocity ──
print("\n── Mean velocity: v_bar = (1/2) int_1^3 x_dot dg = 4a + b ──")
for name, (a, b) in data.items():
    vbar = 4*a + b
    fl = " *FLAT*" if name in flat_coords else ""
    eq = " => zero (kinetic opt)" if vbar == 0 else ""
    print(f"  {name}: v_bar = {str(vbar):>7} = {float(vbar):>+8.4f}{fl}{eq}")

# ── What x(1)=x(3) and x(2)=x(3) mean physically ──
print("\n" + "=" * 80)
print("  PHYSICAL INTERPRETATION")
print("=" * 80)
print("""
  In log-mass space, x_k(g) is coordinate k for generation g.

  b = -4a  =>  x(1) = x(3)
    Generation 1 and 3 have identical coordinate value.
    The curve is a palindrome around g=2 (gen 2 is the extremum).
    Equivalent to ZERO MEAN VELOCITY over [1,3].

  b = -5a  =>  x(2) = x(3)
    Generation 2 and 3 have identical coordinate value.
    The heavy generations are DEGENERATE in this coordinate.
    The turning point g*=5/2 lies between gen 2 and gen 3.
    Mean velocity = -a (non-zero drift).

  TENSION: The kinetic optimum (minimum action) wants b=-4a, but flatness
  forces b=-5a for 3 coordinates. The deviation delta_b = b+4a measures
  how far each coordinate is from the minimum-action configuration.

  For flat coords: delta_b = -a  (fixed by flatness)
  For non-flat:    delta_b is UNEXPLAINED (6 free parameters)
""")

# ── The integers -4 and -5 ──
print("── Origin of the integers -4 and -5 ──")
print(f"""
  -4 comes from: <g>_[1,3] = (1/2) int_1^3 g dg = 4/2 = 2
    => b* = -2 * <g> * a = -4a
    => midpoint (mean) of the generation range

  -5 comes from: x(2) = x(3) => a(4-9) + b(2-3) = 0
    => b = -5a = -(3^2 - 2^2)/(3 - 2) * a
    => encodes the QUADRATIC gap between gen 2 and gen 3

  Ratio: -5/-4 = 5/4  (the "flatness decentering factor")
  The flatness condition shifts g* from 2 to 5/2, i.e., by +1/2.
""")

# ── Verify numerically with actual particle coordinates ──
print("── Verification with actual particle coordinates ──")

# c offsets
c_data = {
    'l_p': F(5, 2), 'l_q': F(5), 'l_r': F(-6),
    'u_p': F(0),    'u_q': F(-4), 'u_r': F(-8),
    'd_p': F(-9,2), 'd_q': F(-9), 'd_r': F(-2),
}

print(f"\n{'coord':>6}  {'x(1)':>10}  {'x(2)':>10}  {'x(3)':>10}  {'x1=x3?':>7}  {'x2=x3?':>7}  {'note'}")
print("-" * 80)

for name, (a, b) in data.items():
    c = c_data[name]
    x1 = a*1 + b*1 + c
    x2 = a*4 + b*2 + c
    x3 = a*9 + b*3 + c

    eq13 = "YES" if x1 == x3 else "no"
    eq23 = "YES" if x2 == x3 else "no"

    fl = "FLAT" if name in flat_coords else ""
    note = ""
    if eq13 == "YES":
        note = "palindrome"
    elif eq23 == "YES":
        note = "gen2-3 degen"

    print(f"{name:>6}  {str(x1):>10}  {str(x2):>10}  {str(x3):>10}  {eq13:>7}  {eq23:>7}  {fl} {note}")
