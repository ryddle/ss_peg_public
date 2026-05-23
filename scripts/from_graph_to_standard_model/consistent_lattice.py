#!/usr/bin/env python3
"""
consistent_lattice.py — Consistent Mass Lattice from Transitivity
==================================================================

Instead of searching for INDIVIDUAL mass ratio formulas, assign a
"spectral address" v_i = (p_i, q_i, r_i) to EACH particle in the
3D effective lattice {ln 2, ln R₃, ln 3}, and find the GLOBAL optimum
that minimizes total error across ALL ratios simultaneously.

Key idea:
  ln(m_i/m_j)_predicted = (p_i-p_j)·ln2 + (q_i-q_j)·lnR₃ + (r_i-r_j)·ln3

This FORCES transitivity: m_i/m_k = (m_i/m_j)·(m_j/m_k) is automatic.

Phase 1: Build consistent lattice via integer optimization
Phase 2: Compare errors with individual-search results
Phase 3: Test if coordinates follow f(gen, I₃, Y, N_c)
Phase 4: Predict UNMEASURED ratios from the lattice

Author: SS-PEG program
Date: April 2026
"""

import numpy as np
from itertools import product as cart_product
import sys, os

np.random.seed(42)

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "consistent_lattice_output.md")
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
# §0  CONSTANTS AND DATA
# ════════════════════════════════════════════════════════════════

R2   = 0.5
R3   = (np.sqrt(57) - 3) / (2 * np.sqrt(17))
R_EE = 1 / np.sqrt(2)
hv3  = 3.0
hv2  = 2.0

LN2  = np.log(2)
LNR3 = np.log(R3)
LN3  = np.log(3)
LOG_BASIS = np.array([LN2, LNR3, LN3])

# 5D building block logs (for back-converting)
x5 = np.array([np.log(R2), np.log(R3), np.log(R_EE), np.log(hv3), np.log(hv2)])
BNAMES = ["R₂", "R₃", "R_EE", "h∨₃", "h∨₂"]

# Masses in MeV (PDG 2024)
MASSES = {
    "e":   0.51099895,
    "mu":  105.6583755,
    "tau": 1776.86,
    "u":   2.16,
    "c":   1270.0,
    "b":   4180.0,
    "d":   4.70,
    "s":   93.4,
    "t":   172500.0,
    "W":   80379.0,
    "Z":   91187.6,
    "H":   125100.0,
}

# All measured ratios and their known (searched) exponent vectors
# Format: (name, numerator, denominator, searched_5d_vector)
RATIOS = [
    ("m_mu/m_e",   "mu",  "e",   [ 0, -4,  1,  3,  0]),
    ("m_tau/m_mu", "tau", "mu",  [ 0, -3, -1,  0,  1]),
    ("m_tau/m_e",  "tau", "e",   [-2, -4,  0,  4,  0]),
    ("m_c/m_u",    "c",   "u",   [-2, -1,  0,  4,  0]),
    ("m_t/m_c",    "t",   "c",   [-1,  0, -1,  1,  4]),
    ("m_s/m_d",    "s",   "d",   [-2,  1,  0,  2,  0]),
    ("m_b/m_s",    "b",   "s",   [ 0,  1,  0,  4,  0]),
    ("m_t/m_b",    "t",   "b",   [-3,  2, -1,  1,  2]),
    ("m_d/m_u",    "d",   "u",   [-1, -2,  0, -1,  0]),
    ("m_b/m_tau",  "b",   "tau", [ 0,  1, -1,  1,  0]),
    ("m_t/m_W",    "t",   "W",   [ 0,  3, -1,  2,  0]),
    ("m_H/m_W",    "H",   "W",   [ 0,  1, -1,  0,  1]),
    ("m_H/m_Z",    "H",   "Z",   [ 1,  2,  0,  2,  0]),
    ("m_H/m_t",    "H",   "t",   [-1, -2,  0, -2,  0]),
]


def vec5d_to_3d(v5):
    """Convert 5D exponent vector to 3D effective lattice coordinates (p, q, r)."""
    a, b, c, d, e = v5
    p = e - a - c / 2.0  # binary coordinate (half-integer if c odd)
    q = b                 # SU(3) spectral coordinate
    r = d                 # Coxeter coordinate
    return np.array([p, q, r])


def log_ratio_from_3d(coord_i, coord_j):
    """Predicted ln(m_i/m_j) from 3D lattice coordinates."""
    delta = coord_i - coord_j
    return delta @ LOG_BASIS


def percent_error(predicted_log, experimental_log):
    """Percent error in the ratio itself."""
    return abs(np.exp(predicted_log - experimental_log) - 1) * 100


print("=" * 72)
print("  CONSISTENT MASS LATTICE FROM TRANSITIVITY CONSTRAINTS")
print("=" * 72)

# ════════════════════════════════════════════════════════════════
# §1  INDIVIDUAL SEARCH RESULTS (baseline)
# ════════════════════════════════════════════════════════════════

print("\n§1. BASELINE — Individual search results (from mass_ratio_search.py)")
print("-" * 72)

print(f"\n  {'Ratio':<14} {'Exp.value':>10} {'Formula':>10} {'Error%':>8} {'5D vector':>24} {'3D (p,q,r)':>16}")
print(f"  {'-'*86}")

individual_errors = []
for name, num, den, v5 in RATIOS:
    exp_ratio = MASSES[num] / MASSES[den]
    log_exp = np.log(exp_ratio)
    log_pred = np.array(v5, dtype=float) @ x5
    err = percent_error(log_pred, log_exp)
    individual_errors.append(err)
    v3 = vec5d_to_3d(v5)
    v5_str = f"({v5[0]:+d},{v5[1]:+d},{v5[2]:+d},{v5[3]:+d},{v5[4]:+d})"
    v3_str = f"({v3[0]:+.1f},{v3[1]:+.0f},{v3[2]:+.0f})"
    print(f"  {name:<14} {exp_ratio:10.4f} {np.exp(log_pred):10.4f} {err:7.3f}%  {v5_str:>24s}  {v3_str:>16s}")

print(f"""
  Mean |error|:  {np.mean(individual_errors):.4f}%
  Max  |error|:  {np.max(individual_errors):.4f}%
  RMS  error:    {np.sqrt(np.mean(np.array(individual_errors)**2)):.4f}%
""")


# ════════════════════════════════════════════════════════════════
# §2  TRANSITIVITY CHECK ON SEARCHED VECTORS
# ════════════════════════════════════════════════════════════════

print("\n§2. TRANSITIVITY CHECK — Are searched vectors internally consistent?")
print("-" * 72)

# Build particle coordinates from the searched vectors
# If searches were consistent, v_i - v_j = searched_vector for ALL ratios,
# and v_i = sum of vectors along any path from anchor to i.

# Use electron as anchor
# Build a spanning tree: e → mu → tau, u → c → t, d → s → b
# Plus cross-sector: u → d, tau → b, t → W, W → H, Z → H

# First: assign coordinates from the FIRST PATH found (spanning tree)
SPANNING_TREE = {
    # particle: (parent, ratio_name, direction)
    # direction = +1 if ratio is particle/parent, -1 if parent/particle
    "e":   None,
    "mu":  ("e",   "m_mu/m_e",   +1),
    "tau": ("mu",  "m_tau/m_mu", +1),
    "u":   None,  # need cross-sector to connect to leptons
    "c":   ("u",   "m_c/m_u",    +1),
    "t":   ("c",   "m_t/m_c",    +1),
    "d":   ("u",   "m_d/m_u",    +1),
    "s":   ("d",   "m_s/m_d",    +1),
    "b":   ("s",   "m_b/m_s",    +1),
    "W":   ("t",   "m_t/m_W",    -1),  # m_t/m_W → W = t - v(m_t/m_W)
    "H":   ("W",   "m_H/m_W",    +1),
    "Z":   ("H",   "m_H/m_Z",    -1),  # m_H/m_Z → Z = H - v(m_H/m_Z)
}

# Lookup ratio vectors
ratio_lookup = {name: np.array(v5, dtype=float) for name, _, _, v5 in RATIOS}

# Build 5D and 3D coordinates from spanning tree
coords_5d = {}
coords_3d = {}

def assign_from_tree(particle):
    if particle in coords_5d:
        return
    info = SPANNING_TREE[particle]
    if info is None:
        # Disconnected root — needs a cross-sector link
        # For quarks, we connect u to the lepton sector via m_b/m_tau
        # Chain: b → tau (via m_b/m_tau), and b is reachable from u via d→s→b
        # But we need u connected to e somehow. There's no direct m_u/m_e ratio.
        # So quarks form a separate connected component in the spanning tree.
        # We'll handle this after building both trees.
        coords_5d[particle] = np.zeros(5)
        coords_3d[particle] = np.zeros(3)
        return

    parent, ratio_name, direction = info
    assign_from_tree(parent)

    v5_ratio = ratio_lookup[ratio_name]
    if direction == +1:
        # particle/parent → particle = parent + v5_ratio (in log space)
        coords_5d[particle] = coords_5d[parent] + v5_ratio
    else:
        # ratio is particle/parent inverted → particle = parent - v5_ratio
        coords_5d[particle] = coords_5d[parent] - v5_ratio

    coords_3d[particle] = vec5d_to_3d(coords_5d[particle])

# Assign lepton tree (rooted at e)
for p in ["e", "mu", "tau"]:
    assign_from_tree(p)

# Assign quark+boson tree (rooted at u)
for p in ["u", "c", "t", "d", "s", "b", "W", "H", "Z"]:
    assign_from_tree(p)

# Now connect the quark tree to the lepton tree via m_b/m_tau
# m_b/m_tau → v(b)_lepton_frame = v(tau) + v(m_b/m_tau)
v_b_from_leptons = coords_5d["tau"] + ratio_lookup["m_b/m_tau"]
v_b_from_quarks  = coords_5d["b"]

# The OFFSET between the two trees
offset_5d = v_b_from_leptons - v_b_from_quarks

print(f"""
  Two disconnected trees: leptons (rooted at e) and quarks+bosons (rooted at u).
  Connected via m_b/m_τ:

  v(b) from lepton tree: v(τ) + v(m_b/m_τ) = {v_b_from_leptons}
  v(b) from quark tree:  v(u) + v(d/u) + v(s/d) + v(b/s) = {v_b_from_quarks}

  Offset (difference) = {offset_5d}
  This offset IS the "spectral address" of u relative to e:
    v(u) - v(e) = {offset_5d}
    In 3D: {vec5d_to_3d(offset_5d)}
""")

# Apply offset: shift all quark+boson coordinates to lepton frame
for particle in ["u", "c", "t", "d", "s", "b", "W", "H", "Z"]:
    coords_5d[particle] = coords_5d[particle] + offset_5d
    coords_3d[particle] = vec5d_to_3d(coords_5d[particle])


# Print all particle coordinates
print("  Particle coordinates (5D, from spanning tree + b/τ bridge):")
print(f"  {'Particle':<8} {'5D vector':>28} {'3D (p,q,r)':>18} {'m_pred (MeV)':>14} {'m_exp (MeV)':>12} {'Error%':>8}")
print(f"  {'-'*92}")

m_anchor = MASSES["e"]
for particle in ["e", "mu", "tau", "u", "d", "c", "s", "t", "b", "W", "Z", "H"]:
    v5 = coords_5d[particle]
    v3 = coords_3d[particle]
    log_mass_pred = np.log(m_anchor) + v5 @ x5
    m_pred = np.exp(log_mass_pred)
    m_exp = MASSES[particle]
    err = abs(m_pred / m_exp - 1) * 100
    v5i = [int(round(vi)) for vi in v5]
    v5_str = f"({v5i[0]:+d},{v5i[1]:+d},{v5i[2]:+d},{v5i[3]:+d},{v5i[4]:+d})"
    v3_str = f"({v3[0]:+5.1f},{v3[1]:+4.0f},{v3[2]:+4.0f})"
    print(f"  {particle:<8} {v5_str:>28s} {v3_str:>18s} {m_pred:14.2f} {m_exp:12.2f} {err:7.2f}%")


# ════════════════════════════════════════════════════════════════
# §3  TRANSITIVITY ERRORS
# ════════════════════════════════════════════════════════════════

print("\n\n§3. TRANSITIVITY ERRORS — Predicted vs experimental for ALL ratios")
print("-" * 72)

print(f"""
  Now ALL ratios are DERIVED from particle coordinates (not searched individually).
  For redundant ratios (not in spanning tree), errors reveal transitivity violations.
""")

print(f"  {'Ratio':<14} {'Exp.value':>10} {'Lattice':>10} {'Error%':>8} {'Searched':>10} {'SrchErr%':>9} {'In tree?':>9}")
print(f"  {'-'*74}")

lattice_errors = []
for name, num, den, v5_searched in RATIOS:
    exp_ratio = MASSES[num] / MASSES[den]
    log_exp = np.log(exp_ratio)

    # From lattice coordinates
    delta_5d = coords_5d[num] - coords_5d[den]
    log_lattice = delta_5d @ x5
    err_lattice = percent_error(log_lattice, log_exp)
    lattice_errors.append(err_lattice)

    # From individual search
    log_searched = np.array(v5_searched, dtype=float) @ x5
    err_searched = percent_error(log_searched, log_exp)

    # Is this ratio in the spanning tree?
    in_tree = False
    for p, info in SPANNING_TREE.items():
        if info and info[1] == name:
            in_tree = True
            break
    # m_b/m_tau used as bridge
    if name == "m_b/m_tau":
        in_tree = True

    tree_str = "YES" if in_tree else "no"
    better = " ◄" if err_lattice < err_searched and not in_tree else ""

    print(f"  {name:<14} {exp_ratio:10.4f} {np.exp(log_lattice):10.4f} {err_lattice:7.3f}%  "
          f"{np.exp(log_searched):10.4f} {err_searched:7.3f}%  {tree_str:>8s}{better}")

print(f"""
  LATTICE (transitivity-forced):
    Mean |error|:  {np.mean(lattice_errors):.4f}%
    Max  |error|:  {np.max(lattice_errors):.4f}%
    RMS:           {np.sqrt(np.mean(np.array(lattice_errors)**2)):.4f}%

  INDIVIDUAL SEARCH (no transitivity):
    Mean |error|:  {np.mean(individual_errors):.4f}%
    Max  |error|:  {np.max(individual_errors):.4f}%
    RMS:           {np.sqrt(np.mean(np.array(individual_errors)**2)):.4f}%
""")

# Key comparison
ratio_mean = np.mean(lattice_errors) / np.mean(individual_errors)
print(f"  RATIO (lattice/individual): {ratio_mean:.2f}×")
if ratio_mean < 1.5:
    print(f"  → Transitivity-forced lattice is COMPARABLE to individual search.")
    print(f"    The lattice structure is real — consistency does NOT degrade quality.")
elif ratio_mean < 3:
    print(f"  → Transitivity causes moderate degradation ({ratio_mean:.1f}×).")
    print(f"    Some formulas may need refinement, but structure is plausible.")
else:
    print(f"  → Transitivity causes severe degradation ({ratio_mean:.1f}×).")
    print(f"    Individual formulas are largely independent — no global lattice.")


# ════════════════════════════════════════════════════════════════
# §4  GLOBAL OPTIMIZATION — Can we do better?
# ════════════════════════════════════════════════════════════════

print("\n\n§4. GLOBAL OPTIMIZATION — Best consistent lattice")
print("-" * 72)

print(f"""
  The spanning tree approach assigns coordinates from one specific path.
  Different spanning trees (different bridges) give different results.
  A GLOBAL optimization finds the particle coordinates that minimize
  total squared error across ALL ratios simultaneously.

  Variables: 3D coordinates (p,q,r) for 11 particles (e fixed at origin)
  Constraints: q ∈ Z, r ∈ Z, p ∈ Z/2
  Objective: minimize Σ (ln(predicted_ratio) - ln(experimental_ratio))²
""")

# List of particles (e is anchor at origin)
PARTICLES = ["e", "mu", "tau", "u", "c", "t", "d", "s", "b", "W", "Z", "H"]
N_PARTICLES = len(PARTICLES)
N_FREE = N_PARTICLES - 1  # 11
pidx = {p: i for i, p in enumerate(PARTICLES)}

# Build the observation equations
# For each ratio m_num/m_den:
#   (p_num - p_den)·ln2 + (q_num - q_den)·lnR₃ + (r_num - r_den)·ln3 = ln(exp_ratio)
# In matrix form: A·v = b  where v = [p_1,q_1,r_1, p_2,q_2,r_2, ..., p_11,q_11,r_11]

N_RATIOS = len(RATIOS)
A = np.zeros((N_RATIOS, 3 * N_FREE))
b = np.zeros(N_RATIOS)

for row, (name, num, den, _) in enumerate(RATIOS):
    exp_ratio = MASSES[num] / MASSES[den]
    b[row] = np.log(exp_ratio)

    # Numerator contribution
    idx_num = pidx[num]
    if idx_num > 0:  # not the anchor
        col_base = 3 * (idx_num - 1)
        A[row, col_base + 0] = LN2    # p coefficient
        A[row, col_base + 1] = LNR3   # q coefficient
        A[row, col_base + 2] = LN3    # r coefficient

    # Denominator contribution (subtract)
    idx_den = pidx[den]
    if idx_den > 0:
        col_base = 3 * (idx_den - 1)
        A[row, col_base + 0] -= LN2
        A[row, col_base + 1] -= LNR3
        A[row, col_base + 2] -= LN3


# Step 1: Solve the RELAXED (continuous) problem
v_real, residuals, rank, sv = np.linalg.lstsq(A, b, rcond=None)

print(f"  System: {N_RATIOS} equations, {3*N_FREE} unknowns, rank = {rank}")
print(f"  Relaxed solution residual: {np.sqrt(np.sum((A @ v_real - b)**2)):.6f}")

# Extract coordinates
coords_relaxed = {"e": np.array([0.0, 0.0, 0.0])}
for i, particle in enumerate(PARTICLES[1:]):
    coords_relaxed[particle] = v_real[3*i:3*i+3]

print(f"\n  Relaxed (continuous) particle coordinates:")
print(f"  {'Particle':<8} {'p':>8} {'q':>8} {'r':>8}")
print(f"  {'-'*36}")
for particle in PARTICLES:
    p, q, r = coords_relaxed[particle]
    print(f"  {particle:<8} {p:8.3f} {q:8.3f} {r:8.3f}")


def total_error(coords, ratios, masses):
    """Compute total squared log-error for given particle coordinates."""
    total_sq = 0
    errors = []
    for name, num, den, _ in ratios:
        log_exp = np.log(masses[num] / masses[den])
        delta = coords[num] - coords[den]
        log_pred = delta @ LOG_BASIS
        err_pct = percent_error(log_pred, log_exp)
        total_sq += (log_pred - log_exp) ** 2
        errors.append(err_pct)
    return total_sq, errors


# Step 2: Use SPANNING TREE coordinates as starting point (much better than rounding)
# The spanning tree already gives 1.04× — start optimization from there
print(f"\n  Starting optimization from SPANNING TREE coordinates (not from rounded relaxed).")
print(f"  Rationale: rounding the relaxed solution lands far from the correct lattice point.")
print(f"  The spanning tree already achieves {np.mean(lattice_errors):.4f}% mean error.\n")

# Convert spanning tree 5D coords to 3D for the optimizer
seed_coords = {}
for particle in PARTICLES:
    seed_coords[particle] = vec5d_to_3d(coords_5d[particle])

print(f"  Seed coordinates (from spanning tree):")
print(f"  {'Particle':<8} {'p':>8} {'q':>8} {'r':>8}")
print(f"  {'-'*32}")
for particle in PARTICLES:
    p, q, r = seed_coords[particle]
    print(f"  {particle:<8} {p:8.1f} {q:8.0f} {r:8.0f}")

best_coords = {k: v.copy() for k, v in seed_coords.items()}
best_total, best_errors = total_error(best_coords, RATIOS, MASSES)

improved = True
iteration = 0
while improved:
    improved = False
    iteration += 1
    for particle in PARTICLES[1:]:  # skip anchor
        current = best_coords[particle].copy()
        # Try perturbations
        for dp in [-0.5, 0, 0.5]:
            for dq in [-1, 0, 1]:
                for dr in [-1, 0, 1]:
                    if dp == 0 and dq == 0 and dr == 0:
                        continue
                    trial = current + np.array([dp, dq, dr])
                    best_coords[particle] = trial
                    trial_total, _ = total_error(best_coords, RATIOS, MASSES)
                    if trial_total < best_total - 1e-12:
                        best_total = trial_total
                        current = trial.copy()
                        improved = True
                    else:
                        best_coords[particle] = current.copy()

print(f"  Converged after {iteration} iterations.")

# Evaluate optimized lattice
_, optimized_errors = total_error(best_coords, RATIOS, MASSES)

print(f"\n  Optimized lattice particle coordinates:")
print(f"  {'Particle':<8} {'p':>8} {'q':>8} {'r':>8}")
print(f"  {'-'*32}")
for particle in PARTICLES:
    p, q, r = best_coords[particle]
    print(f"  {particle:<8} {p:8.1f} {q:8.0f} {r:8.0f}")

print(f"\n  Optimized lattice errors:")
print(f"  {'Ratio':<14} {'Exp.value':>10} {'Lattice':>10} {'Error%':>8} {'Searched':>10} {'SrchErr%':>9} {'Better?':>8}")
print(f"  {'-'*80}")

n_better = 0
for i, (name, num, den, v5_searched) in enumerate(RATIOS):
    exp_ratio = MASSES[num] / MASSES[den]
    delta = best_coords[num] - best_coords[den]
    pred_ratio = np.exp(delta @ LOG_BASIS)
    log_searched = np.array(v5_searched, dtype=float) @ x5
    srch_ratio = np.exp(log_searched)
    srch_err = percent_error(log_searched, np.log(exp_ratio))
    better = "YES" if optimized_errors[i] < srch_err else ""
    if optimized_errors[i] < srch_err:
        n_better += 1
    print(f"  {name:<14} {exp_ratio:10.4f} {pred_ratio:10.4f} {optimized_errors[i]:7.3f}%  "
          f"{srch_ratio:10.4f} {srch_err:7.3f}%  {better:>7s}")

print(f"""
  OPTIMIZED LATTICE (consistent):
    Mean |error|:  {np.mean(optimized_errors):.4f}%
    Max  |error|:  {np.max(optimized_errors):.4f}%
    RMS:           {np.sqrt(np.mean(np.array(optimized_errors)**2)):.4f}%
    #{n_better}/{len(RATIOS)} ratios are BETTER than individual search

  INDIVIDUAL SEARCH:
    Mean |error|:  {np.mean(individual_errors):.4f}%

  RATIO: {np.mean(optimized_errors)/np.mean(individual_errors):.2f}×
""")


# ════════════════════════════════════════════════════════════════
# §5  CONVERT BACK TO 5D VECTORS
# ════════════════════════════════════════════════════════════════

print("\n§5. CONSISTENT 5D EXPONENT VECTORS")
print("-" * 72)

print(f"""
  The optimized 3D coordinates (p,q,r) define the LATTICE consistently.
  For each ratio, the 3D difference translates to a FAMILY of 5D vectors:
    p = e - a - c/2, q = b, r = d
  
  Multiple (a,c,e) triples give the same p, but with different physical
  meaning (number of bridge crossings c vs weak bounces a vs weak height e).
  
  We choose the 5D vector with minimum total |exponent| (Occam).
""")

def best_5d_from_3d(delta_3d):
    """Find the simplest 5D vector consistent with given 3D difference."""
    p_target, q_target, r_target = delta_3d
    b = int(round(q_target))
    d = int(round(r_target))

    # p = e - a - c/2 → search over small c values
    best_v5 = None
    best_norm = 999
    for c in range(-6, 7):
        # p_target = e - a - c/2
        # Need e - a = p_target + c/2
        remainder = p_target + c / 2.0
        if abs(remainder - round(remainder)) > 0.01:
            continue
        # e - a = remainder (integer)
        # To minimize |a| + |e|, choose a and e close to 0
        ea_diff = int(round(remainder))
        # Try a range of a values
        for a in range(-6, 7):
            e = a + ea_diff
            if abs(e) > 6:
                continue
            v5 = [a, b, c, d, e]
            norm = sum(abs(x) for x in v5)
            if norm < best_norm:
                best_norm = norm
                best_v5 = v5

    return best_v5


print(f"  {'Ratio':<14} {'Searched 5D':>28} {'Consistent 5D':>28} {'Same?':>6}")
print(f"  {'-'*80}")

n_same = 0
consistent_vectors = {}
for name, num, den, v5_searched in RATIOS:
    delta_3d = best_coords[num] - best_coords[den]
    v5_consistent = best_5d_from_3d(delta_3d)
    consistent_vectors[name] = v5_consistent

    v5s_str = f"({v5_searched[0]:+d},{v5_searched[1]:+d},{v5_searched[2]:+d},{v5_searched[3]:+d},{v5_searched[4]:+d})"
    if v5_consistent:
        v5c_str = f"({v5_consistent[0]:+d},{v5_consistent[1]:+d},{v5_consistent[2]:+d},{v5_consistent[3]:+d},{v5_consistent[4]:+d})"
        same = "YES" if v5_consistent == list(v5_searched) else "differ"
        if same == "YES":
            n_same += 1
    else:
        v5c_str = "(no solution)"
        same = "N/A"

    print(f"  {name:<14} {v5s_str:>28s} {v5c_str:>28s} {same:>6s}")

print(f"\n  {n_same}/{len(RATIOS)} vectors are IDENTICAL to individual search.")
print(f"  {len(RATIOS)-n_same} vectors DIFFER — these are the transitivity corrections.")


# ════════════════════════════════════════════════════════════════
# §6  QUANTUM NUMBER REGRESSION
# ════════════════════════════════════════════════════════════════

print("\n\n§6. QUANTUM NUMBER REGRESSION — Can we derive coordinates from (gen, I₃, Y, Nc)?")
print("-" * 72)

# Particle quantum numbers: (generation, I₃, Y, Q, Nc)
QN = {
    "e":   (1, -0.5, -1.0, -1,    1),
    "mu":  (2, -0.5, -1.0, -1,    1),
    "tau": (3, -0.5, -1.0, -1,    1),
    "u":   (1,  0.5,  1/3,  2/3,  3),
    "c":   (2,  0.5,  1/3,  2/3,  3),
    "t":   (3,  0.5,  1/3,  2/3,  3),
    "d":   (1, -0.5,  1/3, -1/3,  3),
    "s":   (2, -0.5,  1/3, -1/3,  3),
    "b":   (3, -0.5,  1/3, -1/3,  3),
}

# Bosons don't have generation or standard fermion quantum numbers, skip them
FERMIONS = ["e", "mu", "tau", "u", "c", "t", "d", "s", "b"]

# Build regression: v_i = α + β·gen + γ·I₃ + δ·Y + ε·Nc + ...
# Also try quadratic: + ζ·gen²
print(f"""
  For fermions only, test the ansatz:
    v_i(p,q,r) = α₀ + α₁·gen + α₂·gen² + α₃·I₃ + α₄·Y + α₅·Nc

  If this fits well, the lattice coordinates are DERIVED from quantum numbers,
  not searched. Each coordinate p,q,r would be a known function of the particle.
""")

# Feature matrix
X = np.zeros((len(FERMIONS), 6))
Y_targets = np.zeros((len(FERMIONS), 3))

for i, pname in enumerate(FERMIONS):
    gen, I3, Y, Q, Nc = QN[pname]
    X[i] = [1, gen, gen**2, I3, Y, Nc]
    Y_targets[i] = best_coords[pname]

feature_names = ["const", "gen", "gen²", "I₃", "Y", "Nc"]

print(f"  {'Coord':<6} ", end="")
for fn in feature_names:
    print(f" {fn:>8s}", end="")
print(f"  {'R²':>7s}  {'RMSE':>8s}")
print(f"  {'-'*76}")

regression_coeffs = {}
for k, coord_name in enumerate(["p", "q", "r"]):
    y = Y_targets[:, k]
    beta, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    y_pred = X @ beta
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 1e-12 else float('nan')
    rmse = np.sqrt(ss_res / len(y))
    regression_coeffs[coord_name] = beta

    print(f"  {coord_name:<6} ", end="")
    for coeff in beta:
        print(f" {coeff:>+8.3f}", end="")
    print(f"  {R2:7.4f}  {rmse:8.4f}")


# Print predicted vs actual coordinates
print(f"\n  Predicted vs Actual coordinates:")
print(f"  {'Particle':<8} {'p_act':>8} {'p_pred':>8} {'q_act':>8} {'q_pred':>8} {'r_act':>8} {'r_pred':>8}")
print(f"  {'-'*56}")

for i, pname in enumerate(FERMIONS):
    v_actual = best_coords[pname]
    v_pred = X[i] @ np.column_stack([regression_coeffs[c] for c in ["p", "q", "r"]])
    print(f"  {pname:<8} {v_actual[0]:8.1f} {v_pred[0]:8.2f}  {v_actual[1]:8.0f} {v_pred[1]:8.2f}  {v_actual[2]:8.0f} {v_pred[2]:8.2f}")


# ════════════════════════════════════════════════════════════════
# §7  SECTOR ANALYSIS — What distinguishes the three fermion sectors?
# ════════════════════════════════════════════════════════════════

print("\n\n§7. SECTOR ANALYSIS — Generation stepping in the consistent lattice")
print("-" * 72)

print(f"""
  For each sector, compute the lattice STEP between consecutive generations:
    Δv = v(gen+1) - v(gen)

  If the lattice is real, Δv should have a systematic sector-dependent pattern.
""")

sector_particles = {
    "Charged leptons (I₃=-½, Nc=1)": ["e", "mu", "tau"],
    "Up quarks (I₃=+½, Nc=3)":       ["u", "c", "t"],
    "Down quarks (I₃=-½, Nc=3)":      ["d", "s", "b"],
}

for sector_name, particles in sector_particles.items():
    print(f"\n  {sector_name}:")
    for i in range(len(particles) - 1):
        p1, p2 = particles[i], particles[i+1]
        v1, v2 = best_coords[p1], best_coords[p2]
        delta = v2 - v1
        log_ratio = delta @ LOG_BASIS
        mass_ratio = np.exp(log_ratio)
        print(f"    {p2}/{p1}: Δv = ({delta[0]:+.1f}, {delta[1]:+.0f}, {delta[2]:+.0f})  "
              f"→ ratio = {mass_ratio:.3f}")

    # Generation acceleration
    if len(particles) == 3:
        v1 = best_coords[particles[0]]
        v2 = best_coords[particles[1]]
        v3 = best_coords[particles[2]]
        delta_12 = v2 - v1
        delta_23 = v3 - v2
        gamma_p = delta_23[0] / delta_12[0] if abs(delta_12[0]) > 0.01 else float('inf')
        gamma_q = delta_23[1] / delta_12[1] if abs(delta_12[1]) > 0.01 else float('inf')
        gamma_r = delta_23[2] / delta_12[2] if abs(delta_12[2]) > 0.01 else float('inf')
        print(f"    Acceleration γ = Δ(2→3)/Δ(1→2): p={gamma_p:.3f}, q={gamma_q:.3f}, r={gamma_r:.3f}")


# ════════════════════════════════════════════════════════════════
# §8  SYNTHESIS
# ════════════════════════════════════════════════════════════════

print("\n\n§8. SYNTHESIS")
print("=" * 72)

lattice_mean = np.mean(optimized_errors)
search_mean = np.mean(individual_errors)
ratio = lattice_mean / search_mean

print(f"""
  1. TRANSITIVITY TEST:
     Consistent lattice mean error: {lattice_mean:.4f}%
     Individual search mean error:  {search_mean:.4f}%
     Ratio: {ratio:.2f}×
""")

if ratio < 1.5:
    print(f"     → The mass lattice IS consistent. Forcing transitivity does NOT")
    print(f"       degrade accuracy significantly. The 15 formulas ARE connected")
    print(f"       through a common set of particle coordinates.\n")
else:
    print(f"     → Transitivity costs {ratio:.1f}× accuracy. The lattice has")
    print(f"       imperfections, but may still reveal structure.\n")

# Count how many lattice vectors differ from searched
print(f"  2. VECTOR CHANGES:")
print(f"     {n_same}/{len(RATIOS)} vectors unchanged by transitivity")
print(f"     {len(RATIOS)-n_same}/{len(RATIOS)} vectors modified to enforce consistency")
print(f"     Modified vectors are the CORRECT ones (if the lattice is real).\n")

# Regression quality
print(f"  3. QUANTUM NUMBER DERIVABILITY:")
for k, coord_name in enumerate(["p", "q", "r"]):
    y = Y_targets[:, k]
    y_pred = X @ regression_coeffs[coord_name]
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 1e-12 else float('nan')
    status = "excellent" if R2 > 0.95 else "good" if R2 > 0.8 else "moderate" if R2 > 0.5 else "poor"
    print(f"     {coord_name}: R² = {R2:.4f} ({status})")

print(f"""
  If ALL R² > 0.95, the lattice coordinates ARE derivable from quantum
  numbers, and the mass formulas can be DERIVED (not searched).
  If not, additional structure (non-linear terms) is needed.
""")


# ════════════════════════════════════════════════════════════════
# §9  UNMEASURED PREDICTIONS
# ════════════════════════════════════════════════════════════════

print("\n§9. PREDICTIONS — Ratios NOT used in the fit")
print("-" * 72)

# Additional ratios we can predict from the lattice
EXTRA_RATIOS = [
    ("m_t/m_e",   "t",   "e"),
    ("m_t/m_u",   "t",   "u"),
    ("m_b/m_d",   "b",   "d"),
    ("m_c/m_s",   "c",   "s"),
    ("m_c/m_d",   "c",   "d"),
    ("m_b/m_u",   "b",   "u"),
    ("m_tau/m_d", "tau", "d"),
    ("m_mu/m_u",  "mu",  "u"),
    ("m_W/m_Z",   "W",   "Z"),
    ("m_Z/m_t",   "Z",   "t"),
]

print(f"  {'Ratio':<14} {'Exp.value':>10} {'Lattice':>10} {'Error%':>8}")
print(f"  {'-'*46}")

for name, num, den in EXTRA_RATIOS:
    exp_ratio = MASSES[num] / MASSES[den]
    delta = best_coords[num] - best_coords[den]
    pred_ratio = np.exp(delta @ LOG_BASIS)
    err = abs(pred_ratio / exp_ratio - 1) * 100
    print(f"  {name:<14} {exp_ratio:10.4f} {pred_ratio:10.4f} {err:7.3f}%")


# Close output file
sys.stdout = _orig_stdout
_outfile.close()
print(f"\nOutput saved to {OUTPUT_PATH}")
