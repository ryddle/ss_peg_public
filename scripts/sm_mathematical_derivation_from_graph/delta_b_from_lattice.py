#!/usr/bin/env python3
"""
delta_b_from_lattice.py
=======================
Can the 6 free δb values be DERIVED from the lattice + curvature alone?

The lattice constrains δb ∈ a + Z/2  (so that Δx_{12} = -a + δb ∈ Z/2).
This leaves infinitely many candidates spaced by 1/2.

We test several selection principles to see which one (if any)
picks the *physical* δb uniquely.

Candidates:
  P1. Minimize |δb|            (closest to kinetic optimum b=-4a)
  P2. Minimize |c|             (proximity of g=0 to origin)
  P3. Minimize barrier height  (max V between consecutive holes)
  P4. Minimize integrated V    (∫V dg over [0,4])
  P5. Minimize |x(1)|          (first-generation proximity to origin)
  P6. Some cross-sector constraint (sum rules, etc.)
"""

from fractions import Fraction as F
import numpy as np
import os, pathlib

# ── data ──────────────────────────────────────────────────────────────
#  a_k and b_k for each (sector, coord)
#  Sector order: lepton, up, down.  Coord order: p, q, r.
a_data = {
    ('lepton','p'): F(1),    ('lepton','q'): F(1,2),  ('lepton','r'): F(-3,2),
    ('up','p'):     F(7,4),  ('up','q'):     F(1,2),   ('up','r'):     F(-3,2),
    ('down','p'):   F(-1),   ('down','q'):   F(0),     ('down','r'):   F(1),
}

b_data = {
    ('lepton','p'): F(-7,2), ('lepton','q'): F(-11,2), ('lepton','r'): F(15,2),
    ('up','p'):     F(-13,4),('up','q'):     F(-5,2),  ('up','r'):     F(17,2),
    ('down','p'):   F(5),    ('down','q'):   F(1),     ('down','r'):   F(-1),
}

c_data = {
    ('lepton','p'): F(5,2),  ('lepton','q'): F(5),     ('lepton','r'): F(-6),
    ('up','p'):     F(0),    ('up','q'):     F(-4),     ('up','r'):     F(-8),
    ('down','p'):   F(-9,2), ('down','q'):   F(-9),     ('down','r'):   F(-2),
}

# Flat coordinates (δb forced by flatness, not free)
flat = {('lepton','r'), ('up','q'), ('down','p')}

# Physical δb = b + 4a
db_physical = {}
for key in a_data:
    db_physical[key] = b_data[key] + 4 * a_data[key]

print("="*70)
print("  LATTICE CONSTRAINT ON δb")
print("="*70)
print()
print("  For x(g)=ag²+bg+c to land on Z/2 at g=1,2,3,")
print("  we need δb ≡ a (mod 1/2),  i.e. δb ∈ a + Z/2.")
print()

for key in sorted(a_data.keys()):
    s, k = key
    a = a_data[key]
    db = db_physical[key]
    tag = " [FLAT]" if key in flat else ""
    # Show the coset
    a_mod = a - F(int(2*a), 2)  # fractional part mod 1/2
    print(f"  {s:6s}/{k}: a={float(a):+6.2f}  δb_phys={float(db):+6.2f}"
          f"  coset: δb ∈ {float(a):.2f} + Z/2"
          f"  check: (δb-a) mod 1/2 = {float((db - a) % F(1,2)):.4f}{tag}")

# ── scan δb candidates for the 6 free coordinates ────────────────────
print()
print("="*70)
print("  SELECTION PRINCIPLES — scanning δb ∈ a + Z/2, range [-8, +8]")
print("="*70)

free_keys = sorted([k for k in a_data if k not in flat])

def x_of_g(a, b, c, g):
    """Parabola value (float)."""
    return float(a)*g**2 + float(b)*g + float(c)

def v_of_g(a, b, g):
    """Velocity = dx/dg (float)."""
    return 2*float(a)*g + float(b)

def dist_to_lattice(val):
    """Distance to nearest element of Z/2."""
    return abs(val - round(2*val)/2)

def potential_V(a_vec, b_vec, c_vec, g, coords=3):
    """V(g) = Σ d²(x_k(g), Z/2) for one sector's 3 coords."""
    V = 0.0
    for i in range(coords):
        x = x_of_g(a_vec[i], b_vec[i], c_vec[i], g)
        V += dist_to_lattice(x)**2
    return V

def barrier_height(a_vec, b_vec, c_vec, g1, g2, npts=500):
    """Max V(g) between g1 and g2."""
    gs = np.linspace(g1, g2, npts)
    return max(potential_V(a_vec, b_vec, c_vec, g, 3) for g in gs)

def integrated_V(a_vec, b_vec, c_vec, g1, g2, npts=1000):
    """∫V dg from g1 to g2 (trapezoidal)."""
    gs = np.linspace(g1, g2, npts)
    Vs = [potential_V(a_vec, b_vec, c_vec, g, 3) for g in gs]
    return np.trapezoid(Vs, gs)

# For each free coordinate, scan δb and evaluate metrics
print()
for key in free_keys:
    s, k = key
    a = a_data[key]
    db_phys = db_physical[key]
    c_phys = c_data[key]

    # Generate candidate δb in a + Z/2
    candidates = []
    n_range = 20  # ±10 half-integers
    for n in range(-n_range, n_range+1):
        db_cand = a + F(n, 2)
        candidates.append(db_cand)

    print(f"  {s}/{k}:  a = {float(a):.4f},  δb_phys = {float(db_phys):.4f}")
    print(f"  {'δb':>8s}  {'|δb|':>8s}  {'b':>8s}  {'c(x1∈Z/2)':>10s}  "
          f"{'|c|':>8s}  {'|x(1)|':>8s}  {'|x(0)|':>8s}  {'Σ|x|²':>8s}")

    results = []
    for db_cand in candidates:
        b_cand = db_cand - 4*a
        # c must satisfy x(1) = a + b + c ∈ Z/2
        # Fix c so x(1) matches the nearest half-integer to the physical x(1)
        # Actually, c is free (modulo lattice). Let's fix x(1) = x(1)_physical
        # to isolate the effect of δb on the OTHER generations.
        #
        # Better: for given (a, δb), c is determined by choosing x(1).
        # x(1) = a + b + c = a + (δb-4a) + c = -3a + δb + c
        # So c = x(1) + 3a - δb.
        # x(1) must be in Z/2. The simplest: pick x(1) to minimize |x(1)|.
        # x(1) candidate: nearest Z/2 to 0, given that c must be "reasonable".
        # Actually, let's just derive c from the constraint and scan x(1).

        # For ANY choice of x(1) ∈ Z/2:
        #   c = x(1) + 3a - δb
        #   x(2) = x(1) - a + δb
        #   x(3) = x(1) - 2a + 2δb  (= x(1) + 2(-a + δb))

        # Let's pick x(1) = the PHYSICAL value and see what x(2), x(3) become
        x1_phys = a + (db_phys - 4*a) + c_phys  # = physical x(1)
        # If we change δb but keep x(1) fixed:
        c_cand = x1_phys + 3*a - db_cand
        x2_cand = x1_phys - a + db_cand
        x3_cand = x1_phys - 2*a + 2*db_cand

        # Check lattice membership
        ok2 = (x2_cand.denominator <= 2)
        ok3 = (x3_cand.denominator <= 2)

        if ok2 and ok3:
            abs_db = abs(float(db_cand))
            abs_c = abs(float(c_cand))
            abs_x1 = abs(float(x1_phys))
            abs_x0 = abs(float(x_of_g(a, b_cand, c_cand, 0)))
            sum_x2 = float(x1_phys)**2 + float(x2_cand)**2 + float(x3_cand)**2
            is_phys = "  ◄ PHYSICAL" if db_cand == db_phys else ""
            results.append((db_cand, abs_db, abs_c, abs_x1, abs_x0, sum_x2, is_phys))

    # Sort by |δb|
    results.sort(key=lambda r: abs(float(r[0])))
    # Show only nearby candidates
    for r in results:
        if abs(float(r[0])) <= 6:
            db, adb, ac, ax1, ax0, sx2, tag = r
            print(f"  {float(db):+8.2f}  {adb:8.2f}  "
                  f"{float(db - 4*a_data[key]):+8.2f}  "
                  f"{'✓':>10s}  {ac:8.2f}  {ax1:8.2f}  {ax0:8.2f}  {sx2:8.2f}{tag}")
    print()

# ══════════════════════════════════════════════════════════════════════
# P1: Does |δb| minimization pick the physical value?
# ══════════════════════════════════════════════════════════════════════
print("="*70)
print("  P1: MINIMUM |δb| (closest to kinetic optimum)")
print("="*70)
print()
for key in free_keys:
    a = a_data[key]
    db_phys = db_physical[key]
    # Closest to zero in a + Z/2
    # Offset from zero: a mod 1/2 gives the fractional part
    a_f = float(a)
    # Nearest element of a + Z/2 to zero: round(-a / 0.5) * 0.5 + a... 
    # or equivalently, round(2*0) closest such that db = a + n/2 minimizes |db|
    best_n = round(-2*a_f)  # n such that a + n/2 ≈ 0
    db_min = a + F(int(best_n), 2)
    match = "✓ MATCH" if db_min == db_phys else f"✗ predicts {float(db_min):+.2f}"
    s, k = key
    print(f"  {s:6s}/{k}: δb_phys={float(db_phys):+6.2f}  "
          f"min|δb|={float(db_min):+6.2f}  {match}")

# ══════════════════════════════════════════════════════════════════════
# Now the REAL test: Cross-sector, using the full V(g) = Σ_k d²(xk, Z/2)
# For each sector, vary the free δb values and compute integrated V
# The flat δb is fixed. We scan the 2 free δb per sector.
# ══════════════════════════════════════════════════════════════════════
print()
print("="*70)
print("  P4: MINIMUM INTEGRATED POTENTIAL per sector")
print("      V(g) = Σ_k d²(x_k(g), Z/2), integrate over [0.5, 3.5]")
print("="*70)
print()

sectors = ['lepton', 'up', 'down']
coords = ['p', 'q', 'r']

for sector in sectors:
    # Identify flat and free coords for this sector
    flat_coord = None
    free_coords = []
    for k in coords:
        if (sector, k) in flat:
            flat_coord = k
        else:
            free_coords.append(k)

    print(f"  SECTOR: {sector}  (flat: {flat_coord}, free: {free_coords})")

    # Get fixed (flat) coordinate data
    a_all = [a_data[(sector, k)] for k in coords]
    b_phys_all = [b_data[(sector, k)] for k in coords]
    c_phys_all = [c_data[(sector, k)] for k in coords]
    db_phys_all = [db_physical[(sector, k)] for k in coords]

    # x(1) for each coord — these are the physical lattice points
    x1_all = [a_all[i] + b_phys_all[i] + c_phys_all[i] for i in range(3)]

    # For scanning: fix x(1) to physical value, vary δb for free coords
    # c = x(1) + 3a - δb, so b = δb - 4a, c = x(1) + 3a - δb
    # This ensures x(1) stays on the lattice.

    # Scan the 2 free δb simultaneously
    i_free = [coords.index(k) for k in free_coords]
    a_free = [a_all[i] for i in i_free]
    db_phys_free = [db_phys_all[i] for i in i_free]

    scan_range = range(-12, 13)  # ±6 half-integers around physical
    best_intV = float('inf')
    best_db = None
    phys_intV = None

    results_sector = []

    for n0 in scan_range:
        for n1 in scan_range:
            db0 = a_free[0] + F(int(round(float(db_phys_free[0] - a_free[0]) * 2)) + n0, 2)
            db1 = a_free[1] + F(int(round(float(db_phys_free[1] - a_free[1]) * 2)) + n1, 2)

            # Build full coefficient vectors
            a_vec = [float(a_all[i]) for i in range(3)]
            b_vec = [float(b_phys_all[i]) for i in range(3)]
            c_vec = [float(c_phys_all[i]) for i in range(3)]

            # Override free coords
            free_dbs = [db0, db1]
            for idx_f, idx_c in enumerate(i_free):
                new_b = free_dbs[idx_f] - 4*a_all[idx_c]
                new_c = x1_all[idx_c] + 3*a_all[idx_c] - free_dbs[idx_f]
                b_vec[idx_c] = float(new_b)
                c_vec[idx_c] = float(new_c)

            # Compute integrated V over [0.5, 3.5]
            intV = integrated_V(a_vec, b_vec, c_vec, 0.5, 3.5, npts=500)

            is_phys = (n0 == 0 and n1 == 0)
            if is_phys:
                phys_intV = intV

            results_sector.append((float(db0), float(db1), intV, is_phys))

            if intV < best_intV:
                best_intV = intV
                best_db = (float(db0), float(db1))

    # Show results
    phys_db = (float(db_phys_free[0]), float(db_phys_free[1]))
    print(f"    Physical δb: ({phys_db[0]:+.2f}, {phys_db[1]:+.2f})  "
          f"∫V = {phys_intV:.6f}")
    print(f"    Minimum ∫V:  ({best_db[0]:+.2f}, {best_db[1]:+.2f})  "
          f"∫V = {best_intV:.6f}")
    match = "✓ MATCH!" if (abs(best_db[0] - phys_db[0]) < 0.01 and
                            abs(best_db[1] - phys_db[1]) < 0.01) else "✗ different"
    print(f"    {match}")

    # Show top 5 lowest ∫V
    results_sector.sort(key=lambda r: r[2])
    print(f"    Top 5 lowest ∫V:")
    for i, (d0, d1, iv, ip) in enumerate(results_sector[:5]):
        tag = " ◄ PHYSICAL" if ip else ""
        print(f"      ({d0:+6.2f}, {d1:+6.2f})  ∫V = {iv:.6f}{tag}")
    print()

# ══════════════════════════════════════════════════════════════════════
# P5: Minimum BARRIER height between generations
# ══════════════════════════════════════════════════════════════════════
print("="*70)
print("  P5: MINIMUM MAX-BARRIER between generations")
print("="*70)
print()

for sector in sectors:
    flat_coord = None
    free_coords = []
    for k in coords:
        if (sector, k) in flat:
            flat_coord = k
        else:
            free_coords.append(k)

    a_all = [a_data[(sector, k)] for k in coords]
    b_phys_all = [b_data[(sector, k)] for k in coords]
    c_phys_all = [c_data[(sector, k)] for k in coords]
    db_phys_all = [db_physical[(sector, k)] for k in coords]
    x1_all = [a_all[i] + b_phys_all[i] + c_phys_all[i] for i in range(3)]

    i_free = [coords.index(k) for k in free_coords]
    a_free = [a_all[i] for i in i_free]
    db_phys_free = [db_phys_all[i] for i in i_free]

    scan_range = range(-12, 13)
    best_barr = float('inf')
    best_db = None
    phys_barr = None

    results_sector = []

    for n0 in scan_range:
        for n1 in scan_range:
            db0 = a_free[0] + F(int(round(float(db_phys_free[0] - a_free[0]) * 2)) + n0, 2)
            db1 = a_free[1] + F(int(round(float(db_phys_free[1] - a_free[1]) * 2)) + n1, 2)

            a_vec = [float(a_all[i]) for i in range(3)]
            b_vec = [float(b_phys_all[i]) for i in range(3)]
            c_vec = [float(c_phys_all[i]) for i in range(3)]

            free_dbs = [db0, db1]
            for idx_f, idx_c in enumerate(i_free):
                new_b = free_dbs[idx_f] - 4*a_all[idx_c]
                new_c = x1_all[idx_c] + 3*a_all[idx_c] - free_dbs[idx_f]
                b_vec[idx_c] = float(new_b)
                c_vec[idx_c] = float(new_c)

            # Max barrier across the 3 inter-generation gaps
            barr = max(
                barrier_height(a_vec, b_vec, c_vec, 0.5, 1.5, 200),
                barrier_height(a_vec, b_vec, c_vec, 1.5, 2.5, 200),
                barrier_height(a_vec, b_vec, c_vec, 2.5, 3.5, 200),
            )

            is_phys = (n0 == 0 and n1 == 0)
            if is_phys:
                phys_barr = barr

            results_sector.append((float(db0), float(db1), barr, is_phys))

            if barr < best_barr:
                best_barr = barr
                best_db = (float(db0), float(db1))

    phys_db = (float(db_phys_free[0]), float(db_phys_free[1]))
    print(f"  {sector}: Physical ({phys_db[0]:+.2f}, {phys_db[1]:+.2f})  "
          f"max_barrier = {phys_barr:.6f}")
    print(f"  {sector}: Minimum  ({best_db[0]:+.2f}, {best_db[1]:+.2f})  "
          f"max_barrier = {best_barr:.6f}")
    match = "✓ MATCH!" if (abs(best_db[0] - phys_db[0]) < 0.01 and
                            abs(best_db[1] - phys_db[1]) < 0.01) else "✗ different"
    print(f"  {match}")

    results_sector.sort(key=lambda r: r[2])
    print(f"  Top 5:")
    for i, (d0, d1, br, ip) in enumerate(results_sector[:5]):
        tag = " ◄ PHYSICAL" if ip else ""
        print(f"    ({d0:+6.2f}, {d1:+6.2f})  barrier = {br:.6f}{tag}")
    print()

# ══════════════════════════════════════════════════════════════════════
# P6: Vary x(1) too — full scan of (x1, δb) combinations
#     This asks: can the lattice + curvature alone (without knowing x(1))
#     determine the trajectory?
# ══════════════════════════════════════════════════════════════════════
print("="*70)
print("  P6: FULL SCAN — vary x(1) AND δb simultaneously")
print("      Minimize ∫V over [0.1, 3.9]")
print("="*70)
print()

for sector in sectors:
    a_all = [a_data[(sector, k)] for k in coords]
    b_phys_all = [b_data[(sector, k)] for k in coords]
    c_phys_all = [c_data[(sector, k)] for k in coords]
    db_phys_all = [db_physical[(sector, k)] for k in coords]
    x1_phys = [float(a_all[i] + b_phys_all[i] + c_phys_all[i]) for i in range(3)]

    # For each coord, generate candidates:
    # x(1) ∈ Z/2, δb ∈ a + Z/2, c = x(1) + 3a - δb
    # Scan x(1) in [-10, 16] and δb in [-8, 8]

    best_intV = float('inf')
    best_config = None
    phys_intV = None

    # Build candidate lists per coordinate
    cand_per_coord = []
    for i in range(3):
        a_f = float(a_all[i])
        cands = []
        for x1_n in range(-20, 33):  # x(1) ∈ Z/2 from -10 to 16
            x1 = x1_n / 2.0
            for db_n in range(-16, 17):  # δb in a + Z/2
                db = a_f + db_n / 2.0
                b = db - 4*a_f
                c = x1 + 3*a_f - db
                cands.append((x1, db, b, c))
        cand_per_coord.append(cands)

    # Full combinatorial scan is too large (53^3 × 33^3 per coord).
    # Instead, optimize per coordinate independently (since V = Σ d²_k,
    # and candidates per coord are independent if we minimize ∫V).

    print(f"  SECTOR: {sector}")
    print(f"  Optimizing each coordinate independently (V = Σ V_k):")
    for i, k in enumerate(coords):
        a_f = float(a_all[i])
        tag = " [FLAT]" if (sector, k) in flat else ""

        best_intVk = float('inf')
        best_conf = None
        phys_conf = None

        for x1_n in range(-20, 33):
            x1 = x1_n / 2.0
            for db_n in range(-16, 17):
                db = a_f + db_n / 2.0
                b = db - 4*a_f
                c = x1 + 3*a_f - db

                # Integrated V_k = ∫ d²(x_k(g), Z/2) dg over [0.1, 3.9]
                gs = np.linspace(0.1, 3.9, 300)
                xs = a_f * gs**2 + b * gs + c
                ds = np.minimum(xs % 0.5, 0.5 - (xs % 0.5))
                iv = np.trapezoid(ds**2, gs)

                is_phys = (abs(x1 - x1_phys[i]) < 0.01 and
                           abs(db - float(db_phys_all[i])) < 0.01)
                if is_phys:
                    phys_conf = (x1, db, iv)

                if iv < best_intVk:
                    best_intVk = iv
                    best_conf = (x1, db, iv)

        phys_tag = ""
        if phys_conf:
            x1p, dbp, ivp = phys_conf
            x1b, dbb, ivb = best_conf
            if abs(x1p - x1b) < 0.01 and abs(dbp - dbb) < 0.01:
                phys_tag = " ✓ MATCH"
            else:
                phys_tag = " ✗"
            print(f"    {k}: physical x(1)={x1p:+6.1f} δb={dbp:+6.2f} ∫V={ivp:.6f}")
            print(f"    {k}: best     x(1)={x1b:+6.1f} δb={dbb:+6.2f} ∫V={ivb:.6f}{phys_tag}{tag}")
        else:
            x1b, dbb, ivb = best_conf
            print(f"    {k}: best     x(1)={x1b:+6.1f} δb={dbb:+6.2f} ∫V={ivb:.6f}{tag}")
    print()

# ══════════════════════════════════════════════════════════════════════
# Summary table
# ══════════════════════════════════════════════════════════════════════
print("="*70)
print("  SUMMARY: which principle recovers the physical δb?")
print("="*70)
print()
print("  Principle              | Recovers δb?")
print("  -----------------------|-------------")
print("  P1: min |δb|           | see above ")
print("  P4: min ∫V (fixed x1)  | see above ")
print("  P5: min max-barrier    | see above ")
print("  P6: min ∫V (free x1)   | see above ")
print()
print("  If P6 works: the lattice potential + curvature ALONE")
print("  determine ALL parameters (δb AND x(1)).")
print("  The marble finds the smoothest path through the lattice.")
