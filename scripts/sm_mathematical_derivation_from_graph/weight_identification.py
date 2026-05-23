#!/usr/bin/env python3
"""
weight_identification.py — Lie Algebra Weight Identification
=============================================================

Tests whether the 12-particle coordinates (p, q, r) can be identified
as weights of a representation of a rank-3 Lie algebra.

Sections:
  §1  Root system reconstruction from inter-generational differences
  §2  Dynkin embedding test: project onto simple roots of A₃, B₃, C₃
  §3  Shifted lattice analysis: test m + ρ for Weyl symmetry
  §4  Sector-wise representation decomposition
  §5  Extension to 12 points (bosons)

Author: SS-PEG program
Date: April 2026
"""

import numpy as np
from itertools import combinations, product as cart_product
from collections import Counter
import sys, os

np.random.seed(42)

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "weight_identification_output.md")
_orig_stdout = sys.stdout
_outfile = open(OUTPUT_PATH, "w", encoding="utf-8")

class _TeeWriter:
    def __init__(self, console, f):
        self.console = console
        self.file = f
    def write(self, text):
        self.console.write(text)
        self.file.write(text)
    def flush(self):
        self.console.flush()
        self.file.flush()

sys.stdout = _TeeWriter(_orig_stdout, _outfile)

# ════════════════════════════════════════════════════════════════
# PARTICLE DATA
# ════════════════════════════════════════════════════════════════

PARTICLES = {
    "e":   {"pqr": ( 0.0,   0, 0), "gen": 1, "I3": -0.5, "Y": -1.0, "Q": -1,   "Nc": 1},
    "mu":  {"pqr": (-0.5,  -4, 3), "gen": 2, "I3": -0.5, "Y": -1.0, "Q": -1,   "Nc": 1},
    "tau": {"pqr": ( 1.0,  -7, 3), "gen": 3, "I3": -0.5, "Y": -1.0, "Q": -1,   "Nc": 1},
    "u":   {"pqr": (-1.5,  -6,-1), "gen": 1, "I3":  0.5, "Y":  1/3, "Q":  2/3, "Nc": 3},
    "c":   {"pqr": ( 0.5,  -7, 3), "gen": 2, "I3":  0.5, "Y":  1/3, "Q":  2/3, "Nc": 3},
    "t":   {"pqr": ( 6.0,  -7, 4), "gen": 3, "I3":  0.5, "Y":  1/3, "Q":  2/3, "Nc": 3},
    "d":   {"pqr": (-0.5,  -8,-2), "gen": 1, "I3": -0.5, "Y":  1/3, "Q": -1/3, "Nc": 3},
    "s":   {"pqr": ( 1.5,  -7, 0), "gen": 2, "I3": -0.5, "Y":  1/3, "Q": -1/3, "Nc": 3},
    "b":   {"pqr": ( 1.5,  -6, 4), "gen": 3, "I3": -0.5, "Y":  1/3, "Q": -1/3, "Nc": 3},
    "W":   {"pqr": ( 5.5, -10, 2), "gen": 0, "I3":  0.0, "Y":  0.0, "Q":  0.0, "Nc": 0},
    "Z":   {"pqr": ( 8.0, -11, 0), "gen": 0, "I3":  0.0, "Y":  0.0, "Q":  0.0, "Nc": 0},
    "H":   {"pqr": ( 7.0,  -9, 2), "gen": 0, "I3":  0.0, "Y":  0.0, "Q":  0.0, "Nc": 0},
}

FERMIONS = {k: v for k, v in PARTICLES.items() if v["gen"] > 0}
BOSONS   = {k: v for k, v in PARTICLES.items() if v["gen"] == 0}

SECTORS = {
    "lepton": ["e", "mu", "tau"],
    "up":     ["u", "c", "t"],
    "down":   ["d", "s", "b"],
}

# ════════════════════════════════════════════════════════════════
# §1. ROOT SYSTEM FROM INTER-GENERATIONAL DIFFERENCES
# ════════════════════════════════════════════════════════════════

print("=" * 72)
print("  §1. ROOT SYSTEM RECONSTRUCTION")
print("=" * 72)

print("""
  Strategy: In a Lie algebra representation, the difference between
  consecutive weights along a root chain is a root α. We extract
  all gen→gen+1 differences as candidate roots.
""")

# Extract generational roots per sector
gen_roots = {}
for sname, names in SECTORS.items():
    c = [np.array(PARTICLES[n]["pqr"]) for n in names]
    d12 = c[1] - c[0]
    d23 = c[2] - c[1]
    gen_roots[sname] = (d12, d23)
    print(f"  {sname}:")
    print(f"    α₁₂ = {names[0]}→{names[1]} = ({d12[0]:+.1f}, {d12[1]:+.0f}, {d12[2]:+.0f})")
    print(f"    α₂₃ = {names[1]}→{names[2]} = ({d23[0]:+.1f}, {d23[1]:+.0f}, {d23[2]:+.0f})")

# Collect all 6 generational roots
all_roots = []
for sname in SECTORS:
    all_roots.extend(gen_roots[sname])

# Check relationships between roots
print(f"\n  Root relationships:")
for i in range(len(all_roots)):
    for j in range(i+1, len(all_roots)):
        ri, rj = all_roots[i], all_roots[j]
        dot = np.dot(ri, rj)
        ni, nj = np.linalg.norm(ri), np.linalg.norm(rj)
        cos_angle = dot / (ni * nj) if ni > 0 and nj > 0 else 0
        angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
        if abs(ni) > 1e-10 and abs(nj) > 1e-10:
            ratio = ni / nj
            # Cartan integer: 2 * dot(αi, αj) / dot(αj, αj)
            cartan_ij = 2 * dot / np.dot(rj, rj) if np.dot(rj, rj) > 0 else 0
            cartan_ji = 2 * dot / np.dot(ri, ri) if np.dot(ri, ri) > 0 else 0

# Cross-sector roots: connecting different sectors at same generation
print(f"\n  Cross-sector roots (same generation):")
for g in [1, 2, 3]:
    particles_g = [k for k, v in FERMIONS.items() if v["gen"] == g]
    for i, n1 in enumerate(particles_g):
        for n2 in particles_g[i+1:]:
            delta = np.array(PARTICLES[n2]["pqr"]) - np.array(PARTICLES[n1]["pqr"])
            print(f"    {n1}→{n2} (gen {g}): ({delta[0]:+.1f}, {delta[1]:+.0f}, {delta[2]:+.0f})")

# ════════════════════════════════════════════════════════════════
# §2. DYNKIN EMBEDDING TEST
# ════════════════════════════════════════════════════════════════

print()
print("=" * 72)
print("  §2. DYNKIN EMBEDDING TEST")
print("=" * 72)

print("""
  Testing: can we find a linear map L: R³→R³ such that L(pqr)
  are Dynkin labels (integers) of some rank-3 Lie algebra?

  For A₃: weights are (a₁, a₂, a₃) with aᵢ ∈ Z
  For B₃: similar but with different norms for short/long roots
""")

# Use fermion coordinates as columns of the weight matrix
ferm_names_ordered = ["e", "mu", "tau", "u", "c", "t", "d", "s", "b"]
W = np.array([list(PARTICLES[n]["pqr"]) for n in ferm_names_ordered])  # 9×3

# A₃ simple roots in standard basis (e₁-e₂, e₂-e₃, e₃-e₄)
# We need to find if our roots embed into A₃
# A₃ Cartan matrix:
A3_cartan = np.array([[ 2, -1,  0],
                       [-1,  2, -1],
                       [ 0, -1,  2]])

B3_cartan = np.array([[ 2, -1,  0],
                       [-1,  2, -1],
                       [ 0, -1,  2]])  # B3 differs in root length ratio

C3_cartan = np.array([[ 2, -1,  0],
                       [-1,  2, -2],
                       [ 0, -1,  2]])

print(f"\n  Sector-wise Cartan matrix construction:")
print(f"  Using α₁₂ and α₂₃ as candidate simple roots per sector")

for sname in SECTORS:
    alpha1, alpha2 = gen_roots[sname]
    n1, n2 = np.linalg.norm(alpha1), np.linalg.norm(alpha2)

    if n1 > 1e-10 and n2 > 1e-10:
        a12 = 2 * np.dot(alpha1, alpha2) / np.dot(alpha2, alpha2)
        a21 = 2 * np.dot(alpha1, alpha2) / np.dot(alpha1, alpha1)
        ratio = n1 / n2

        print(f"\n  {sname}:")
        print(f"    α₁ = {alpha1}, |α₁| = {n1:.4f}")
        print(f"    α₂ = {alpha2}, |α₂| = {n2:.4f}")
        print(f"    |α₁|/|α₂| = {ratio:.4f}")
        print(f"    ⟨α₁,α₂⟩/⟨α₂,α₂⟩ × 2 = {a12:.4f}")
        print(f"    ⟨α₁,α₂⟩/⟨α₁,α₁⟩ × 2 = {a21:.4f}")
        print(f"    Cartan-like: [[2, {a12:.2f}], [{a21:.2f}, 2]]")

        # Check if Cartan integers are integer (necessary for Lie algebra)
        a12_int = abs(a12 - round(a12)) < 0.1
        a21_int = abs(a21 - round(a21)) < 0.1
        print(f"    Integer Cartan? {a12_int and a21_int}")

# ════════════════════════════════════════════════════════════════
# §3. SHIFTED LATTICE & WEYL SYMMETRY
# ════════════════════════════════════════════════════════════════

print()
print("=" * 72)
print("  §3. SHIFTED LATTICE & WEYL SYMMETRY")
print("=" * 72)

print("""
  Weyl character formula uses shifted weights: λ + ρ where ρ = half
  sum of positive roots. Test symmetry of shifted weights.
""")

# Compute center of mass of fermion coordinates
ferm_coords = np.array([list(PARTICLES[n]["pqr"]) for n in ferm_names_ordered])
centroid = np.mean(ferm_coords, axis=0)
print(f"  Fermion centroid: ({centroid[0]:.4f}, {centroid[1]:.4f}, {centroid[2]:.4f})")

# Shift all weights by -centroid
shifted = ferm_coords - centroid
print(f"\n  Shifted coordinates (w - centroid):")
print(f"  {'Particle':<8} {'Δp':>8} {'Δq':>8} {'Δr':>8} {'|Δ|':>8}")
print(f"  {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*8}")
for i, name in enumerate(ferm_names_ordered):
    s = shifted[i]
    print(f"  {name:<8} {s[0]:>+8.4f} {s[1]:>+8.4f} {s[2]:>+8.4f} {np.linalg.norm(s):>8.4f}")

# Check for reflection symmetries
print(f"\n  Reflection symmetry test:")
print(f"  For each shifted weight w, check if -w is also present (up to tolerance)")
for i, name in enumerate(ferm_names_ordered):
    neg = -shifted[i]
    dists = np.linalg.norm(shifted - neg, axis=1)
    min_dist = np.min(dists)
    match_idx = np.argmin(dists)
    match_name = ferm_names_ordered[match_idx]
    is_match = min_dist < 0.5
    print(f"    {name} ↔ -{name} closest to {match_name} (dist={min_dist:.4f}) {'✓' if is_match else '✗'}")

# Per-sector shifted analysis
print(f"\n  Per-sector shifted coordinates (relative to sector centroid):")
for sname, names in SECTORS.items():
    coords = np.array([list(PARTICLES[n]["pqr"]) for n in names])
    sc = np.mean(coords, axis=0)
    s = coords - sc
    print(f"\n  {sname} (centroid: {sc}):")
    for i, name in enumerate(names):
        print(f"    {name}: ({s[i][0]:+.4f}, {s[i][1]:+.4f}, {s[i][2]:+.4f})  |w|={np.linalg.norm(s[i]):.4f}")
    # Check if gen1 + gen3 = 2 * gen2 (linear embedding)
    midpoint_err = np.linalg.norm(s[0] + s[2])  # should be 0 if gen2 = centroid
    print(f"    gen1 + gen3 - 2·gen2 deviation: {midpoint_err:.6f}")

# ════════════════════════════════════════════════════════════════
# §4. SECTOR-WISE REPRESENTATION DECOMPOSITION
# ════════════════════════════════════════════════════════════════

print()
print("=" * 72)
print("  §4. SECTOR-WISE REPRESENTATION DECOMPOSITION")
print("=" * 72)

print("""
  Each sector has 3 particles (generations 1, 2, 3). In SU(2)
  language, 3 states = spin-1 (triplet). Test if the generational
  structure is consistent with SU(2)_generation.
""")

# For SU(2) triplet (j=1): weights are +1, 0, -1
# Check if generational coordinates are equally spaced along some axis
for sname, names in SECTORS.items():
    coords = np.array([list(PARTICLES[n]["pqr"]) for n in names])
    # Direction from gen1 to gen3
    direction = coords[2] - coords[0]
    norm = np.linalg.norm(direction)
    direction_hat = direction / norm if norm > 0 else direction

    # Project all 3 onto this direction
    projs = [np.dot(c - coords[0], direction_hat) for c in coords]
    expected = [0, norm/2, norm]  # equal spacing

    print(f"\n  {sname}: gen1→gen3 direction = {direction}, |d| = {norm:.4f}")
    print(f"    Projections: {[f'{p:.4f}' for p in projs]}")
    print(f"    Expected (equal): {[f'{e:.4f}' for e in expected]}")
    deviation = abs(projs[1] - norm/2)
    print(f"    Gen2 midpoint deviation: {deviation:.6f}")
    print(f"    → {'SU(2) triplet consistent ✓' if deviation < 0.01 else 'NOT an SU(2) triplet ✗'}")

    # Perpendicular component of gen2
    perp = coords[1] - coords[0] - projs[1] * direction_hat
    print(f"    Gen2 perpendicular component: {perp} (|⊥| = {np.linalg.norm(perp):.6f})")

# ════════════════════════════════════════════════════════════════
# §5. EXTENSION TO 12 POINTS (INCLUDING BOSONS)
# ════════════════════════════════════════════════════════════════

print()
print("=" * 72)
print("  §5. BOSONIC EXTENSION")
print("=" * 72)

print("""
  The 3 bosons (W, Z, H) are connected to the fermion lattice via
  the spanning tree. Test their position relative to the fermion
  weight structure.
""")

all_names = ["e", "mu", "tau", "u", "c", "t", "d", "s", "b", "W", "Z", "H"]
all_coords = np.array([list(PARTICLES[n]["pqr"]) for n in all_names])

# Boson coordinates relative to fermion centroid
print(f"  Boson positions relative to fermion centroid {centroid}:")
for name in ["W", "Z", "H"]:
    c = np.array(PARTICLES[name]["pqr"])
    shifted_c = c - centroid
    print(f"    {name}: pqr={c}, shifted=({shifted_c[0]:+.4f}, {shifted_c[1]:+.4f}, {shifted_c[2]:+.4f})")

# Are bosons in the convex hull of fermions?
print(f"\n  Convex hull test (are bosons inside the fermion polytope?):")
from scipy.spatial import ConvexHull
try:
    hull = ConvexHull(ferm_coords)
    print(f"    Fermion convex hull volume: {hull.volume:.4f}")
    for name in ["W", "Z", "H"]:
        c = np.array(PARTICLES[name]["pqr"])
        # Test if point is inside convex hull
        # Add point, check if hull volume increases
        extended = np.vstack([ferm_coords, c.reshape(1, -1)])
        hull2 = ConvexHull(extended)
        inside = abs(hull2.volume - hull.volume) < 1e-6
        print(f"    {name}: {'inside ✓' if inside else 'outside ✗'} (ΔV = {hull2.volume - hull.volume:.6f})")
except Exception as e:
    print(f"    ConvexHull not available: {e}")

# Boson-fermion distance analysis
print(f"\n  Nearest fermion to each boson:")
for bname in ["W", "Z", "H"]:
    bc = np.array(PARTICLES[bname]["pqr"])
    min_dist = float('inf')
    near_name = ""
    for fname in ferm_names_ordered:
        fc = np.array(PARTICLES[fname]["pqr"])
        d = np.linalg.norm(bc - fc)
        if d < min_dist:
            min_dist = d
            near_name = fname
    delta = np.array(PARTICLES[near_name]["pqr"]) - bc
    print(f"    {bname} ↔ {near_name}: dist={min_dist:.4f}, Δ=({delta[0]:+.1f}, {delta[1]:+.0f}, {delta[2]:+.0f})")

# Test: are boson differences from t/τ simple root-like vectors?
print(f"\n  Boson chains (from spanning tree):")
spanning_chain = [("t", "W"), ("W", "H"), ("H", "Z")]
for n1, n2 in spanning_chain:
    c1, c2 = np.array(PARTICLES[n1]["pqr"]), np.array(PARTICLES[n2]["pqr"])
    delta = c2 - c1
    print(f"    {n1}→{n2}: ({delta[0]:+.1f}, {delta[1]:+.0f}, {delta[2]:+.0f})  |Δ|={np.linalg.norm(delta):.4f}")

# 12-point Gram matrix
print(f"\n  Full 12-point coordinate matrix analysis:")
print(f"    Rank: {np.linalg.matrix_rank(all_coords, tol=1e-6)}")
U, S, Vt = np.linalg.svd(all_coords, full_matrices=False)
print(f"    Singular values: {S}")
print(f"    Principal directions:")
for i, v in enumerate(Vt):
    print(f"      PC{i+1}: ({v[0]:+.4f}, {v[1]:+.4f}, {v[2]:+.4f})  explains {100*S[i]**2/np.sum(S**2):.1f}%")

# ════════════════════════════════════════════════════════════════
# SUMMARY
# ════════════════════════════════════════════════════════════════

print()
print("=" * 72)
print("  SUMMARY")
print("=" * 72)

# Collect key results
for sname, names in SECTORS.items():
    coords = np.array([list(PARTICLES[n]["pqr"]) for n in names])
    direction = coords[2] - coords[0]
    norm = np.linalg.norm(direction)
    direction_hat = direction / norm if norm > 0 else direction
    projs = [np.dot(c - coords[0], direction_hat) for c in coords]
    deviation = abs(projs[1] - norm/2)
    equal_spaced = deviation < 0.01

print(f"""
  §1: Root system reconstruction
      6 generational roots extracted (2 per sector)
      Cross-sector roots also tabulated

  §2: Dynkin embedding
      Cartan-like matrix computed per sector
      Integer Cartan condition tested

  §3: Weyl symmetry
      Fermion centroid: {centroid}
      Reflection matches evaluated

  §4: SU(2)_gen decomposition
      Tested equally-spaced triplet hypothesis per sector

  §5: Bosonic extension
      Spanning tree chain t→W→H→Z analyzed
      Boson positions relative to fermion convex hull determined
      12-point coordinate matrix rank: {np.linalg.matrix_rank(all_coords, tol=1e-6)}
""")

# Close tee
sys.stdout = _orig_stdout
_outfile.close()
print(f"Output saved to: {OUTPUT_PATH}")
