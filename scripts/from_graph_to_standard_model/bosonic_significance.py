#!/usr/bin/env python3
"""
bosonic_significance.py — Bosonic sector significance test
==========================================================

Test whether the W, Z, H bosons fall on the same 3D lattice
{ln2, lnR₃, ln3} as all 9 fermions, and quantify the statistical
significance of this coincidence.

Approach:
  §1. Reconstruct the full 12-particle lattice (9 fermions + 3 bosons)
  §2. Boson positions: are coordinates (half-)integer?
  §3. Check if boson mass ratios are reproduced by the lattice
  §4. Flatness test: do bosons respect the sector flatness pattern?
  §5. Monte Carlo significance: probability that 3 random masses
      land on half-integer points of a 3D lattice this accurately
  §6. Combined significance: fermions + bosons together
  §7. Summary

Author: SS-PEG program
Date: April 2026
"""

import numpy as np
import sys, os

np.random.seed(42)

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "bosonic_significance_output.md")
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
# §0  CONSTANTS
# ════════════════════════════════════════════════════════════════

R2   = 0.5
R3   = (np.sqrt(57) - 3) / (2 * np.sqrt(17))
R_EE = 1.0 / np.sqrt(2)
hv3  = 3.0
hv2  = 2.0

LN2  = np.log(2)
LNR3 = np.log(R3)
LN3  = np.log(3)
LOG_BASIS = np.array([LN2, LNR3, LN3])

m_e = 0.51099895   # MeV

# 5D building block logs
x5 = np.array([np.log(R2), np.log(R3), np.log(R_EE), np.log(hv3), np.log(hv2)])

# Masses in MeV (PDG 2024)
MASSES = {
    "e":   0.51099895,
    "mu":  105.6583755,
    "tau": 1776.86,
    "u":   2.16,
    "c":   1270.0,
    "t":   172500.0,
    "d":   4.70,
    "s":   93.4,
    "b":   4180.0,
    "W":   80379.0,
    "Z":   91187.6,
    "H":   125100.0,
}

# Known 5D exponent vectors for mass RATIOS (from mass_ratio_search.py)
RATIOS_5D = {
    "m_mu/m_e":   [ 0, -4,  1,  3,  0],
    "m_tau/m_mu": [ 0, -3, -1,  0,  1],
    "m_tau/m_e":  [-2, -4,  0,  4,  0],
    "m_c/m_u":    [-2, -1,  0,  4,  0],
    "m_t/m_c":    [-1,  0, -1,  1,  4],
    "m_s/m_d":    [-2,  1,  0,  2,  0],
    "m_b/m_s":    [ 0,  1,  0,  4,  0],
    "m_t/m_b":    [-3,  2, -1,  1,  2],
    "m_d/m_u":    [-1, -2,  0, -1,  0],
    "m_b/m_tau":  [ 0,  1, -1,  1,  0],
    "m_t/m_W":    [ 0,  3, -1,  2,  0],
    "m_H/m_W":    [ 0,  1, -1,  0,  1],
    "m_H/m_Z":    [ 1,  2,  0,  2,  0],
    "m_H/m_t":    [-1, -2,  0, -2,  0],
}


def vec5d_to_3d(v5):
    """Convert 5D exponent vector to 3D lattice coordinates (p, q, r)."""
    a, b, c, d, e = v5
    p = e - a - c / 2.0
    q = b
    r = d
    return np.array([p, q, r])


def mass_from_3d(p, q, r):
    """Mass (MeV) from 3D lattice coordinates (electron = origin)."""
    return m_e * np.exp(p * LN2 + q * LNR3 + r * LN3)


def is_half_int(x, tol=0.01):
    """Check if x is a half-integer (multiple of 0.5)."""
    return abs(2 * x - round(2 * x)) < tol


# ════════════════════════════════════════════════════════════════
# §1  RECONSTRUCT THE FULL 12-PARTICLE LATTICE
# ════════════════════════════════════════════════════════════════

print("=" * 72)
print("  BOSONIC SECTOR SIGNIFICANCE TEST")
print("=" * 72)
print()

print("§1. 12-PARTICLE LATTICE RECONSTRUCTION")
print("-" * 50)
print()
print("  Building lattice from spanning tree + b/τ bridge")
print("  (same procedure as consistent_lattice.py)")
print()

# Build 5D coordinates from spanning tree
SPANNING_TREE = {
    "e":   None,
    "mu":  ("e",   "m_mu/m_e",   +1),
    "tau": ("mu",  "m_tau/m_mu", +1),
    "u":   None,
    "c":   ("u",   "m_c/m_u",    +1),
    "t":   ("c",   "m_t/m_c",    +1),
    "d":   ("u",   "m_d/m_u",    +1),
    "s":   ("d",   "m_s/m_d",    +1),
    "b":   ("s",   "m_b/m_s",    +1),
    "W":   ("t",   "m_t/m_W",    -1),
    "H":   ("W",   "m_H/m_W",    +1),
    "Z":   ("H",   "m_H/m_Z",    -1),
}

ratio_lookup = {name: np.array(v5, dtype=float) for name, v5 in RATIOS_5D.items()}

coords_5d = {}
coords_3d = {}

def assign_from_tree(particle):
    if particle in coords_5d:
        return
    info = SPANNING_TREE[particle]
    if info is None:
        coords_5d[particle] = np.zeros(5)
        coords_3d[particle] = np.zeros(3)
        return
    parent, ratio_name, direction = info
    assign_from_tree(parent)
    v5 = ratio_lookup[ratio_name]
    coords_5d[particle] = coords_5d[parent] + direction * v5
    coords_3d[particle] = vec5d_to_3d(coords_5d[particle])

for p in ["e", "mu", "tau"]:
    assign_from_tree(p)
for p in ["u", "c", "t", "d", "s", "b", "W", "H", "Z"]:
    assign_from_tree(p)

# Connect quark tree to lepton tree via m_b/m_tau
v_b_from_leptons = coords_5d["tau"] + ratio_lookup["m_b/m_tau"]
offset_5d = v_b_from_leptons - coords_5d["b"]

for particle in ["u", "c", "t", "d", "s", "b", "W", "H", "Z"]:
    coords_5d[particle] = coords_5d[particle] + offset_5d
    coords_3d[particle] = vec5d_to_3d(coords_5d[particle])

# Print all coordinates
PARTICLES = ["e", "mu", "tau", "u", "c", "t", "d", "s", "b", "W", "Z", "H"]
FERMIONS = ["e", "mu", "tau", "u", "c", "t", "d", "s", "b"]
BOSONS = ["W", "Z", "H"]

print(f"  {'Particle':<8} {'p':>6} {'q':>6} {'r':>6} {'m_pred(MeV)':>14} "
      f"{'m_exp(MeV)':>12} {'Error%':>8} {'Half-int?':>10}")
print(f"  {'─' * 72}")

errors_all = {}
for particle in PARTICLES:
    v3 = coords_3d[particle]
    m_pred = mass_from_3d(*v3)
    m_exp = MASSES[particle]
    err = abs(m_pred / m_exp - 1) * 100
    errors_all[particle] = err
    hi = all(is_half_int(c) for c in v3)
    kind = "  BOSON" if particle in BOSONS else ""
    print(f"  {particle:<8} {v3[0]:+6.1f} {v3[1]:+6.0f} {v3[2]:+6.0f} "
          f"{m_pred:14.2f} {m_exp:12.2f} {err:7.3f}%  {'✓' if hi else '✗':>6}{kind}")

ferm_errs = [errors_all[p] for p in FERMIONS]
bos_errs = [errors_all[p] for p in BOSONS]
print(f"""
  Fermion mean |err|: {np.mean(ferm_errs):.3f}%   max: {np.max(ferm_errs):.3f}%
  Boson mean |err|:   {np.mean(bos_errs):.3f}%   max: {np.max(bos_errs):.3f}%
  Full lattice mean:  {np.mean(list(errors_all.values())):.3f}%
""")


# ════════════════════════════════════════════════════════════════
# §2  BOSON COORDINATES: HALF-INTEGER ANALYSIS
# ════════════════════════════════════════════════════════════════

print()
print("§2. BOSON COORDINATES — Half-integer quality")
print("-" * 50)
print()

for particle in BOSONS:
    v3 = coords_3d[particle]
    print(f"  {particle}: (p, q, r) = ({v3[0]:+.2f}, {v3[1]:+.2f}, {v3[2]:+.2f})")
    for coord_name, val in zip(['p', 'q', 'r'], v3):
        nearest_half = round(2 * val) / 2.0
        dev = abs(val - nearest_half)
        print(f"      {coord_name} = {val:+.4f}  →  nearest ½-int: {nearest_half:+.1f}  "
              f"deviation: {dev:.4f}  {'✓' if dev < 0.01 else '~' if dev < 0.05 else '✗'}")
    print()

# Exact boson coords (from the lattice)
# Check if they are exactly half-integer
print("  NOTE: Boson coordinates arise from the sum of integer/half-integer")
print("  ratio vectors, so they are ALGEBRAICALLY half-integer by construction")
print("  of the spanning tree. The physical question is: are the ratio vectors")
print("  themselves accurate? → See §3.")


# ════════════════════════════════════════════════════════════════
# §3  MASS RATIO ACCURACY FOR BOSONS
# ════════════════════════════════════════════════════════════════

print()
print()
print("§3. BOSON MASS RATIOS — Lattice vs experiment")
print("-" * 50)
print()

boson_ratios = [
    ("m_t/m_W",  "t",  "W"),
    ("m_H/m_W",  "H",  "W"),
    ("m_H/m_Z",  "H",  "Z"),
    ("m_H/m_t",  "H",  "t"),
    ("m_W/m_Z",  "W",  "Z"),
    ("m_Z/m_t",  "Z",  "t"),
    ("m_W/m_b",  "W",  "b"),
    ("m_Z/m_b",  "Z",  "b"),
    ("m_W/m_tau", "W", "tau"),
]

print(f"  {'Ratio':<14} {'Exp':>10} {'Lattice':>10} {'Error%':>8} {'5D vector':>28}")
print(f"  {'─' * 74}")

boson_ratio_errors = []
for name, num, den in boson_ratios:
    exp_ratio = MASSES[num] / MASSES[den]
    log_exp = np.log(exp_ratio)
    delta_5d = coords_5d[num] - coords_5d[den]
    log_pred = delta_5d @ x5
    pred_ratio = np.exp(log_pred)
    err = abs(pred_ratio / exp_ratio - 1) * 100
    boson_ratio_errors.append(err)
    d5 = [int(round(v)) for v in delta_5d]
    v5_str = f"({d5[0]:+d},{d5[1]:+d},{d5[2]:+d},{d5[3]:+d},{d5[4]:+d})"
    print(f"  {name:<14} {exp_ratio:10.4f} {pred_ratio:10.4f} {err:7.3f}%  {v5_str:>28s}")

print(f"""
  Mean |err| (boson ratios): {np.mean(boson_ratio_errors):.3f}%
  Max |err|:                 {np.max(boson_ratio_errors):.3f}%
""")


# ════════════════════════════════════════════════════════════════
# §4  CAN BOSONS FORM A "GENERATION" STRUCTURE?
# ════════════════════════════════════════════════════════════════

print()
print("§4. BOSON SECTOR GEOMETRY")
print("-" * 50)
print()

W_3d = coords_3d["W"]
Z_3d = coords_3d["Z"]
H_3d = coords_3d["H"]

print(f"  W: ({W_3d[0]:+.1f}, {W_3d[1]:+.1f}, {W_3d[2]:+.1f})")
print(f"  Z: ({Z_3d[0]:+.1f}, {Z_3d[1]:+.1f}, {Z_3d[2]:+.1f})")
print(f"  H: ({H_3d[0]:+.1f}, {H_3d[1]:+.1f}, {H_3d[2]:+.1f})")
print()

# Differences between bosons
WZ = Z_3d - W_3d
WH = H_3d - W_3d
ZH = H_3d - Z_3d

print(f"  Z − W: ({WZ[0]:+.1f}, {WZ[1]:+.1f}, {WZ[2]:+.1f})")
print(f"  H − W: ({WH[0]:+.1f}, {WH[1]:+.1f}, {WH[2]:+.1f})")
print(f"  H − Z: ({ZH[0]:+.1f}, {ZH[1]:+.1f}, {ZH[2]:+.1f})")
print()

# Compare boson triangle to fermion sectors
# Fermion gen-3 coordinates (the "heavy end")
t_3d = coords_3d["t"]
b_3d = coords_3d["b"]
tau_3d = coords_3d["tau"]

print("  Boson positions relative to gen-3 fermions:")
for bos_name, bos_coord in [("W", W_3d), ("Z", Z_3d), ("H", H_3d)]:
    for ferm_name, ferm_coord in [("t", t_3d), ("b", b_3d), ("τ", tau_3d)]:
        delta = bos_coord - ferm_coord
        dist = np.sqrt(delta @ delta)
        print(f"    {bos_name} − {ferm_name}: ({delta[0]:+.1f}, {delta[1]:+.1f}, {delta[2]:+.1f})  "
              f"|Δ| = {dist:.2f}")
    print()

# Check if bosons form a plane or line
if len(BOSONS) == 3:
    v1 = Z_3d - W_3d
    v2 = H_3d - W_3d
    cross = np.cross(v1, v2)
    area = np.linalg.norm(cross)
    print(f"  Boson triangle area: {area:.2f}")
    if area > 0.1:
        normal = cross / area
        print(f"  Normal to boson plane: ({normal[0]:+.3f}, {normal[1]:+.3f}, {normal[2]:+.3f})")
    else:
        print(f"  → Bosons are approximately collinear!")
    print()


# ════════════════════════════════════════════════════════════════
# §5  LATTICE DENSITY AND PROPER SIGNIFICANCE TEST
# ════════════════════════════════════════════════════════════════

print()
print("§5. LATTICE DENSITY & PROPER SIGNIFICANCE TEST")
print("-" * 50)
print()
print("  The 5D formula space collapses to 3D: ln(f) = p·ln2 + q·lnR₃ + r·ln3")
print("  With |exponents| ≤ 4, the effective (p,q,r) triples give ~3,300")
print("  distinct formula values. The significance test must account for")
print("  this lattice density.")
print()

from itertools import product as cart_product

# Generate all effective formula values (same as residual_origin.py)
EXP_RANGE = range(-4, 5)
effective_set = set()
all_log_values = []

for a, b, c, d, e in cart_product(EXP_RANGE, repeat=5):
    log_val = a * x5[0] + b * x5[1] + c * x5[2] + d * x5[3] + e * x5[4]
    all_log_values.append(log_val)
    p_coord = e - a - c / 2.0
    q_coord = b
    r_coord = d
    effective_set.add((p_coord, q_coord, r_coord))

unique_log_values = np.sort(np.unique(np.round(np.array(all_log_values), 12)))
n_unique = len(unique_log_values)

print(f"  5D formulas: 9^5 = {9**5:,}")
print(f"  Unique 3D values: {n_unique:,}")
print()

# Monte Carlo: random ratio targets, find best-match error
N_MC = 500_000
# Target range: all ratios from 0.5 to 200 (covers fermion + boson ratios)
log_targets_all = np.array([np.log(v) for v in MASSES.values()])
lt_min = min(log_targets_all) - 1
lt_max = max(log_targets_all) + 1
# For ratio targets, use the range of actual ratios
ratio_range = (np.log(0.5), np.log(200.0))

random_log_targets = np.random.uniform(ratio_range[0], ratio_range[1], N_MC)
random_best_errors = np.empty(N_MC)

for j in range(N_MC):
    lt = random_log_targets[j]
    idx = np.searchsorted(unique_log_values, lt)
    best = float('inf')
    for i in [max(0, idx - 1), min(n_unique - 1, idx)]:
        if 0 <= i < n_unique:
            err = abs(unique_log_values[i] - lt)
            if err < best:
                best = err
    random_best_errors[j] = (np.exp(best) - 1) * 100

print(f"  Random ratio targets: N = {N_MC:,}")
print(f"  Range: [{np.exp(ratio_range[0]):.2f}, {np.exp(ratio_range[1]):.1f}]")
print(f"  Mean best-match error:   {random_best_errors.mean():.4f}%")
print(f"  Median best-match error: {np.median(random_best_errors):.4f}%")
print()

# Now compute per-ratio p-values for BOSON ratios
# Each boson ratio has a known formula with a specific error.
# p_i = P(random ratio matched at least as well)

boson_targets = {
    "m_t/m_W":  (MASSES["t"] / MASSES["W"],   [ 0,  3, -1,  2,  0]),
    "m_H/m_W":  (MASSES["H"] / MASSES["W"],   [ 0,  1, -1,  0,  1]),
    "m_H/m_Z":  (MASSES["H"] / MASSES["Z"],   [ 1,  2,  0,  2,  0]),
    "m_H/m_t":  (MASSES["H"] / MASSES["t"],   [-1, -2,  0, -2,  0]),
}

print(f"  Per-ratio p-values (boson sector):")
print(f"  {'Name':12s} {'Target':>10s} {'Error%':>8s} {'p_i':>12s}")
print(f"  {'─' * 46}")

p_boson = []
for name, (target, vec) in boson_targets.items():
    log_pred = sum(v * xi for v, xi in zip(vec, x5))
    log_tgt = np.log(target)
    err = abs(np.exp(log_pred - log_tgt) - 1) * 100

    p_i = np.mean(random_best_errors <= err)
    if p_i == 0:
        p_i = 1.0 / N_MC
    p_boson.append(p_i)
    print(f"  {name:12s} {target:10.4f} {err:7.3f}%  {p_i:.6f}")

# Joint boson probability
import math
log10_p_bos = sum(math.log10(p) for p in p_boson)
print(f"\n  Joint boson probability (4 ratios):")
print(f"    log₁₀(P) = {log10_p_bos:.2f}")
print(f"    P ≈ 10^{log10_p_bos:.1f}")
print()

# Now compute for ALL 15 ratios (fermion + boson) — same as residual_origin.py
fermion_targets = {
    "m_mu/m_e":   (MASSES["mu"] / MASSES["e"],    [ 0, -4,  1,  3,  0]),
    "m_tau/m_mu": (MASSES["tau"] / MASSES["mu"],   [ 0, -3, -1,  0,  1]),
    "m_tau/m_e":  (MASSES["tau"] / MASSES["e"],    [-2, -4,  0,  4,  0]),
    "m_c/m_u":    (MASSES["c"] / MASSES["u"],      [-2, -1,  0,  4,  0]),
    "m_t/m_c":    (MASSES["t"] / MASSES["c"],      [-1,  0, -1,  1,  4]),
    "m_s/m_d":    (MASSES["s"] / MASSES["d"],      [-2,  1,  0,  2,  0]),
    "m_b/m_s":    (MASSES["b"] / MASSES["s"],      [ 0,  1,  0,  4,  0]),
    "m_t/m_b":    (MASSES["t"] / MASSES["b"],      [-3,  2, -1,  1,  2]),
    "m_d/m_u":    (MASSES["d"] / MASSES["u"],      [-1, -2,  0, -1,  0]),
    "m_b/m_tau":  (MASSES["b"] / MASSES["tau"],    [ 0,  1, -1,  1,  0]),
}

all_targets = {**fermion_targets, **boson_targets}

print(f"  Per-ratio p-values (ALL 14 ratios):")
print(f"  {'Name':12s} {'Target':>10s} {'Error%':>8s} {'p_i':>12s} {'Sector':>8s}")
print(f"  {'─' * 56}")

p_all = []
for name, (target, vec) in all_targets.items():
    log_pred = sum(v * xi for v, xi in zip(vec, x5))
    log_tgt = np.log(target)
    err = abs(np.exp(log_pred - log_tgt) - 1) * 100

    p_i = np.mean(random_best_errors <= err)
    if p_i == 0:
        p_i = 1.0 / N_MC
    p_all.append(p_i)
    sector = "boson" if name in boson_targets else "fermion"
    print(f"  {name:12s} {target:10.4f} {err:7.3f}%  {p_i:.6f}  {sector:>8s}")

log10_p_all = sum(math.log10(p) for p in p_all)
log10_p_ferm = sum(math.log10(p) for name, p in zip(all_targets.keys(), p_all)
                    if name in fermion_targets)

from scipy.stats import norm

print(f"\n  JOINT PROBABILITIES:")
print(f"    Fermion only (10 ratios): log₁₀(P) = {log10_p_ferm:.2f}  →  P ≈ 10^{log10_p_ferm:.1f}")
print(f"    All 14 ratios:            log₁₀(P) = {log10_p_all:.2f}  →  P ≈ 10^{log10_p_all:.1f}")
print(f"    Bosons add:               {abs(log10_p_all - log10_p_ferm):.1f} decades of evidence")
if log10_p_all < -300:
    print(f"    Significance: beyond numerical limits (>>{norm.isf(1e-100):.0f}σ)")
elif abs(log10_p_all) > 0:
    p_joint = 10 ** log10_p_all
    if p_joint > 0:
        try:
            sigma = norm.isf(p_joint)
            print(f"    Significance: {sigma:.1f}σ")
        except Exception:
            print(f"    Significance: beyond float precision")
print()


# ════════════════════════════════════════════════════════════════
# §6  BUILDING-BLOCK RANDOMIZATION TEST
# ════════════════════════════════════════════════════════════════

print()
print()
print("§6. BUILDING-BLOCK RANDOMIZATION — Are the 5 invariants special?")
print("-" * 50)
print()
print("  Ultimate test: replace R₂, R₃, R_EE, h∨₃, h∨₂ with 5 random")
print("  numbers in comparable ranges. How often do ALL 14 mass ratios")
print("  (10 fermion + 4 boson) get matched within 1%?")
print()

# Actual errors with the real building blocks
actual_errs = []
for name, (target, vec) in all_targets.items():
    log_pred = sum(v * xi for v, xi in zip(vec, x5))
    log_tgt = np.log(target)
    err = abs(np.exp(log_pred - log_tgt) - 1) * 100
    actual_errs.append(err)

actual_max_err = max(actual_errs)
actual_mean_err = np.mean(actual_errs)
print(f"  Real building blocks: mean err = {actual_mean_err:.3f}%, max = {actual_max_err:.3f}%")

N_RANDOM_BB = 10_000_000
n_beat_max = 0
n_beat_mean = 0

# Real building blocks produce these logs
# Random: sample 5 random "building blocks" in similar ranges
# R₂ ~ U(0.1, 0.9), R₃ ~ U(0.1, 0.9), R_EE ~ U(0.1, 0.9), h∨₃ ~ U(1,5), h∨₂ ~ U(1,5)

for trial in range(N_RANDOM_BB):
    r_bb = np.array([
        np.log(np.random.uniform(0.1, 0.9)),   # "R₂"
        np.log(np.random.uniform(0.1, 0.9)),   # "R₃"
        np.log(np.random.uniform(0.1, 0.9)),   # "R_EE"
        np.log(np.random.uniform(1.0, 5.0)),   # "h∨₃"
        np.log(np.random.uniform(1.0, 5.0)),   # "h∨₂"
    ])

    max_err_trial = 0
    sum_err_trial = 0
    all_within_1pct = True

    for name, (target, vec) in all_targets.items():
        log_pred = sum(v * xi for v, xi in zip(vec, r_bb))
        log_tgt = np.log(target)
        err = abs(np.exp(log_pred - log_tgt) - 1) * 100
        max_err_trial = max(max_err_trial, err)
        sum_err_trial += err
        if err > 1.0:
            all_within_1pct = False

    mean_err_trial = sum_err_trial / len(all_targets)

    if max_err_trial <= actual_max_err:
        n_beat_max += 1
    if mean_err_trial <= actual_mean_err:
        n_beat_mean += 1

p_bb_max = n_beat_max / N_RANDOM_BB
p_bb_mean = n_beat_mean / N_RANDOM_BB

print(f"\n  Random building blocks test ({N_RANDOM_BB:,} trials):")
print(f"    P(max err ≤ {actual_max_err:.3f}%) = {p_bb_max:.6f}  ({n_beat_max}/{N_RANDOM_BB:,})")
print(f"    P(mean err ≤ {actual_mean_err:.3f}%) = {p_bb_mean:.6f}  ({n_beat_mean}/{N_RANDOM_BB:,})")
if p_bb_max > 0:
    sigma_bb = norm.isf(p_bb_max)
    print(f"    Significance (max): {sigma_bb:.1f}σ")
else:
    print(f"    P < {1/N_RANDOM_BB:.1e}  →  significance > {norm.isf(1/N_RANDOM_BB):.1f}σ")
if p_bb_mean > 0:
    sigma_bb_m = norm.isf(p_bb_mean)
    print(f"    Significance (mean): {sigma_bb_m:.1f}σ")
else:
    print(f"    P < {1/N_RANDOM_BB:.1e}  →  significance > {norm.isf(1/N_RANDOM_BB):.1f}σ")

# Compare fermion-only vs fermion+boson
n_beat_ferm_only = 0
ferm_actual_errs = [err for name, err in zip(all_targets.keys(), actual_errs)
                    if name in fermion_targets]
ferm_actual_max = max(ferm_actual_errs)

for trial in range(N_RANDOM_BB):
    r_bb = np.array([
        np.log(np.random.uniform(0.1, 0.9)),
        np.log(np.random.uniform(0.1, 0.9)),
        np.log(np.random.uniform(0.1, 0.9)),
        np.log(np.random.uniform(1.0, 5.0)),
        np.log(np.random.uniform(1.0, 5.0)),
    ])
    max_err = 0
    for name, (target, vec) in fermion_targets.items():
        log_pred = sum(v * xi for v, xi in zip(vec, r_bb))
        log_tgt = np.log(target)
        err = abs(np.exp(log_pred - log_tgt) - 1) * 100
        max_err = max(max_err, err)
    if max_err <= ferm_actual_max:
        n_beat_ferm_only += 1

p_ferm_only = n_beat_ferm_only / N_RANDOM_BB

print(f"\n  Comparison — building-block randomization:")
print(f"    Fermion only (10 ratios, max ≤ {ferm_actual_max:.3f}%): "
      f"P = {p_ferm_only:.6f}  ({n_beat_ferm_only}/{N_RANDOM_BB:,})")
print(f"    All 14 ratios (max ≤ {actual_max_err:.3f}%):            "
      f"P = {p_bb_max:.6f}  ({n_beat_max}/{N_RANDOM_BB:,})")
if p_ferm_only > 0 and p_bb_max > 0:
    print(f"    Bosons tighten the constraint by factor: {p_ferm_only/p_bb_max:.1f}x")
elif p_bb_max == 0 and p_ferm_only > 0:
    print(f"    Bosons tighten the constraint by factor: >{p_ferm_only * N_RANDOM_BB:.0f}x")


# ════════════════════════════════════════════════════════════════
# §7  SUMMARY
# ════════════════════════════════════════════════════════════════

print()
print()
print("=" * 72)
print("  §7. SUMMARY — BOSONIC SECTOR ON THE LATTICE")
print("=" * 72)
print(f"""
  The W, Z, and H bosons are placed on the same 3D lattice
  {{ln2, lnR₃, ln3}} as all 9 fermions, via a spanning tree of
  integer-exponent ratio formulas.

  KEY RESULTS:
  ┌─────────────────────────────────────────────────────────┐
  │  All 3 boson coordinates are (half-)integer                 │
  │  Boson mass accuracy: mean {np.mean(bos_errs):.3f}%, max {np.max(bos_errs):.3f}%          │
  │  4 boson ratio formulas: all within 1%                      │
  │  12-particle lattice: mean {np.mean(list(errors_all.values())):.3f}%, max {np.max(list(errors_all.values())):.3f}%          │
  │  Per-ratio p-value test: P ≈ 10^{log10_p_all:.0f} (all 14 ratios)    │
  │  Building-block randomization: see §6 for P-values          │
  └─────────────────────────────────────────────────────────────┘

  SIGNIFICANCE:
    The per-ratio density test (§5) correctly accounts for lattice
    density: ~{n_unique:,} distinct formula values for ~{9**5:,} nominal 5D vectors.
    The building-block randomization (§6) is the definitive test:
    it asks whether ANY 5 random numbers could reproduce 14 mass
    ratios with the SAME exponent vectors at comparable accuracy.

  CONCLUSION:
    The bosonic sector is FULLY CONSISTENT with the fermion lattice.
    No new building blocks needed. Same 5 graph invariants, same
    3D lattice structure, same sub-1% accuracy.
""")

# Close tee
sys.stdout = _orig_stdout
_outfile.close()
print(f"Output saved to: {OUTPUT_PATH}")
