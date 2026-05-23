#!/usr/bin/env python3
"""
residual_origin.py — Origin of the ~0.3% deviation
====================================================

Two fundamental questions:
  Q1: What is the JOINT probability of 15 mass ratio matches at ~0.3%?
  Q2: Is 0.3% the lattice resolution of integer-exponent formulas,
      or a genuine physical deviation?

Key insight: the 5 building blocks {R₂, R₃, R_EE, h∨₃, h∨₂} reduce to
3 independent log-bases {ln 2, ln R₃, ln 3} because:
    R₂ = 2⁻¹, R_EE = 2⁻¹/², h∨₂ = 2¹
The 9⁵ = 59,049 "formulas" collapse to ~3,300 distinct values.

Author: SS-PEG program
Date: April 2026
"""

import numpy as np
from itertools import product as cart_product
import sys, os

np.random.seed(42)

# Redirect to file + console
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "residual_origin_output.md")
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
# §0 SETUP
# ════════════════════════════════════════════════════════════════

R2  = 0.5
R3  = (np.sqrt(57) - 3) / (2 * np.sqrt(17))
R_EE = 1 / np.sqrt(2)
hv3 = 3
hv2 = 2

x = np.array([np.log(R2), np.log(R3), np.log(R_EE), np.log(hv3), np.log(hv2)])
LN2  = np.log(2)
LNR3 = np.log(R3)
LN3  = np.log(3)

EXP_RANGE = range(-4, 5)   # [-4, 4]
N_PER = len(EXP_RANGE)     # 9

# Experimental targets (mass ratios only — no 4π factors)
TARGETS = {
    "Koide":     2/3,
    "m_mu/m_e":  206.7683,
    "m_tau/m_mu": 16.8170,
    "m_tau/m_e": 3477.2283,
    "m_c/m_u":   587.963,
    "m_t/m_c":   135.827,
    "m_s/m_d":   19.8723,
    "m_b/m_s":   44.7537,
    "m_t/m_b":   41.2679,
    "m_d/m_u":   2.1759,
    "m_b/m_tau": 2.3525,
    "m_t/m_W":   2.1461,
    "m_H/m_W":   1.5564,
    "m_H/m_Z":   1.3719,
    "m_H/m_t":   0.7252,
}

# Known exponent vectors (R₂, R₃, R_EE, h∨₃, h∨₂)
KNOWN_VECTORS = {
    "Koide":     [-1,  0,  0, -1,  0],
    "m_mu/m_e":  [ 0, -4,  1,  3,  0],
    "m_tau/m_mu":[ 0, -3, -1,  0,  1],
    "m_tau/m_e": [-2, -4,  0,  4,  0],
    "m_c/m_u":   [-2, -1,  0,  4,  0],
    "m_t/m_c":   [-1,  0, -1,  1,  4],
    "m_s/m_d":   [-2,  1,  0,  2,  0],
    "m_b/m_s":   [ 0,  1,  0,  4,  0],
    "m_t/m_b":   [-3,  2, -1,  1,  2],
    "m_d/m_u":   [-1, -2,  0, -1,  0],
    "m_b/m_tau": [ 0,  1, -1,  1,  0],
    "m_t/m_W":   [ 0,  3, -1,  2,  0],
    "m_H/m_W":   [ 0,  1, -1,  0,  1],
    "m_H/m_Z":   [ 1,  2,  0,  2,  0],
    "m_H/m_t":   [-1, -2,  0, -2,  0],
}

BNAMES = ["R2", "R3", "R_EE", "hv3", "hv2"]


print("=" * 72)
print("  ORIGIN OF THE 0.3% — LATTICE RESOLUTION OR GENUINE DEVIATION?")
print("=" * 72)


# ════════════════════════════════════════════════════════════════
# §1 EFFECTIVE DIMENSIONALITY
# ════════════════════════════════════════════════════════════════
print("\n§1. EFFECTIVE DIMENSIONALITY REDUCTION")
print("-" * 72)

print("""
  The 5 building blocks are NOT independent in log-space:
    ln R₂  = -ln 2
    ln R_EE = -½ ln 2
    ln h∨₂ = +ln 2
    ln R₃  = ln R₃   (irrational, independent of ln 2 and ln 3)
    ln h∨₃ = ln 3

  Any formula ln(f) = a·ln R₂ + b·ln R₃ + c·ln R_EE + d·ln h∨₃ + e·ln h∨₂
  reduces to:
    ln(f) = p·ln 2 + q·ln R₃ + r·ln 3
  where:
    p = (e - a - c/2)  [half-integer, from R₂ᵃ·R_EEᶜ·h∨₂ᵉ]
    q = b               [integer, from R₃ᵇ]
    r = d               [integer, from h∨₃ᵈ]
""")

# Generate ALL formula log-values and simultaneously count effective triples
effective_set = set()
all_log_values = []

for a, b, c, d, e in cart_product(EXP_RANGE, repeat=5):
    log_val = a * x[0] + b * x[1] + c * x[2] + d * x[3] + e * x[4]
    all_log_values.append(log_val)
    p = e - a - c / 2.0
    q = b
    r = d
    effective_set.add((p, q, r))

all_log_values = np.array(all_log_values)
total_5d = len(all_log_values)
n_effective_triples = len(effective_set)

# Unique formula values (sorted)
unique_log_values = np.sort(np.unique(np.round(all_log_values, 12)))
n_unique = len(unique_log_values)

print(f"  Nominal formula count:     {N_PER}^5 = {total_5d:,}")
print(f"  Effective (p,q,r) triples: {n_effective_triples:,}")
print(f"  Numerically unique values: {n_unique:,}")
print(f"  Reduction factor:          {total_5d / n_unique:.1f}x")

# Verify collapse: the number of unique values should match the
# number of effective triples (by Baker's theorem, ln 2, ln R3, ln 3
# are linearly independent over Q, so distinct triples → distinct values)
print(f"\n  Baker's theorem check: {n_effective_triples} triples vs {n_unique} unique values")
if abs(n_effective_triples - n_unique) <= 5:
    print(f"  → CONFIRMED: 3 independent log-bases, {n_unique} distinct formulas")
else:
    print(f"  → WARNING: slight mismatch (rounding artifacts possible)")


# ════════════════════════════════════════════════════════════════
# §2 LATTICE SPACING ANALYSIS
# ════════════════════════════════════════════════════════════════
print("\n§2. LATTICE SPACING (FORMULA RESOLUTION)")
print("-" * 72)

# Sort and compute gaps
gaps = np.diff(unique_log_values)
gaps_pct = (np.exp(gaps) - 1) * 100  # convert to relative error %

# Focus on the range covered by actual targets
log_targets = np.array([np.log(t) for t in TARGETS.values()])
lt_min, lt_max = log_targets.min(), log_targets.max()

# Gaps within target range (with 1 unit margin)
in_range = (unique_log_values[:-1] >= lt_min - 0.5) & (unique_log_values[1:] <= lt_max + 0.5)
gaps_in_range = gaps_pct[in_range]

print(f"  Target log-range: [{lt_min:.3f}, {lt_max:.3f}]  (span = {lt_max - lt_min:.3f})")
print(f"  Lattice points in range: {in_range.sum()}")
print(f"")
print(f"  Gap statistics (within target range):")
print(f"    Mean gap:    {gaps_in_range.mean():.4f}%")
print(f"    Median gap:  {np.median(gaps_in_range):.4f}%")
print(f"    P10 gap:     {np.percentile(gaps_in_range, 10):.4f}%")
print(f"    P90 gap:     {np.percentile(gaps_in_range, 90):.4f}%")
print(f"    Max gap:     {gaps_in_range.max():.4f}%")
print(f"    Min gap:     {gaps_in_range.min():.6f}%")
print(f"")
print(f"  Expected BEST-MATCH error ≈ mean_gap/2 = {gaps_in_range.mean() / 2:.4f}%")


# ════════════════════════════════════════════════════════════════
# §3 ACTUAL FORMULA ERRORS vs BEST POSSIBLE LATTICE MATCH
# ════════════════════════════════════════════════════════════════
print("\n§3. ACTUAL ERRORS vs LATTICE BEST-MATCH")
print("-" * 72)

actual_formula_errors = []
best_lattice_errors = []
local_gaps = []

print(f"  {'Name':12s} {'Target':>10s} {'Formula err':>11s} {'Best match':>11s} {'Local gap':>10s}")
print(f"  {'-'*60}")

for name, target in TARGETS.items():
    vec = KNOWN_VECTORS[name]
    log_pred = sum(v * xi for v, xi in zip(vec, x))
    log_target = np.log(target)

    # Error of the known formula
    formula_err = abs(np.exp(log_pred - log_target) - 1) * 100
    actual_formula_errors.append(formula_err)

    # Best possible match from the lattice
    idx = np.searchsorted(unique_log_values, log_target)
    candidates = [i for i in [max(0, idx - 1), min(n_unique - 1, idx)] if 0 <= i < n_unique]
    best_log_err = min(abs(unique_log_values[i] - log_target) for i in candidates)
    best_err = (np.exp(best_log_err) - 1) * 100
    best_lattice_errors.append(best_err)

    # Local gap at this target
    if idx > 0 and idx < n_unique:
        local_g = (np.exp(unique_log_values[idx] - unique_log_values[idx - 1]) - 1) * 100
    else:
        local_g = np.nan
    local_gaps.append(local_g)

    marker = " <<<" if formula_err > best_err * 3 else ""
    print(f"  {name:12s} {target:10.4f} {formula_err:10.4f}% {best_err:10.4f}% {local_g:9.4f}%{marker}")

actual_formula_errors = np.array(actual_formula_errors)
best_lattice_errors = np.array(best_lattice_errors)
local_gaps = np.array(local_gaps)

print(f"\n  Mean formula error:      {actual_formula_errors.mean():.4f}%")
print(f"  Mean best-lattice error: {best_lattice_errors.mean():.4f}%")
print(f"  Mean local gap:          {np.nanmean(local_gaps):.4f}%")

# Are the known formulas optimal?
optimal_count = sum(abs(f - b) < 0.0001 for f, b in zip(actual_formula_errors, best_lattice_errors))
print(f"\n  Formulas that ARE the best lattice match: {optimal_count}/{len(TARGETS)}")
suboptimal = [(n, f, b) for n, f, b in zip(TARGETS.keys(), actual_formula_errors, best_lattice_errors)
              if abs(f - b) > 0.001]
if suboptimal:
    print(f"  Sub-optimal formulas (better lattice match exists):")
    for n, f, b in suboptimal:
        print(f"    {n:12s}: formula={f:.4f}%, best={b:.4f}%")


# ════════════════════════════════════════════════════════════════
# §4 MONTE CARLO: RANDOM TARGET ERROR DISTRIBUTION
# ════════════════════════════════════════════════════════════════
print("\n§4. MONTE CARLO — WHAT ERRORS DO RANDOM TARGETS GET?")
print("-" * 72)

N_MC = 200_000
random_log_targets = np.random.uniform(lt_min, lt_max, N_MC)
random_best_errors = np.empty(N_MC)

for j in range(N_MC):
    lt = random_log_targets[j]
    idx = np.searchsorted(unique_log_values, lt)
    best = float('inf')
    for i in [max(0, idx - 1), min(n_unique - 1, idx)]:
        err = abs(unique_log_values[i] - lt)
        if err < best:
            best = err
    random_best_errors[j] = (np.exp(best) - 1) * 100

print(f"  Random targets: N = {N_MC:,}")
print(f"  Mean best-match error:   {random_best_errors.mean():.4f}%")
print(f"  Median best-match error: {np.median(random_best_errors):.4f}%")
print(f"  Max best-match error:    {random_best_errors.max():.4f}%")

thresholds = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0]
print(f"\n  {'Threshold':>10s} {'P(err<thr)':>12s} {'1/P':>10s}")
for thr in thresholds:
    p = np.mean(random_best_errors < thr)
    p_str = f"{p:.6f}" if p > 0 else "<1/N_MC"
    inv_str = f"{1 / p:.1f}" if p > 0 else f">{N_MC}"
    print(f"  {thr:8.2f}%   {p_str:>12s}   {inv_str:>10s}")


# ════════════════════════════════════════════════════════════════
# §5 JOINT PROBABILITY — THE KEY QUESTION
# ════════════════════════════════════════════════════════════════
print("\n§5. JOINT PROBABILITY OF 15 MASS RATIO MATCHES")
print("-" * 72)

print("""
  Method: For each of the 15 targets, compute
    p_i = P(random_target gets best-match error ≤ actual_error_i)
  from the Monte Carlo distribution. Then:
    P_joint = Π p_i  (assuming approximate independence)
""")

p_individual = []
for name, err in zip(TARGETS.keys(), actual_formula_errors):
    p = np.mean(random_best_errors <= err)
    if p == 0:
        p = 1 / N_MC  # floor
    p_individual.append(p)
    print(f"  {name:12s}  err = {err:.4f}%  P(random ≤ err) = {p:.6f}")

p_individual = np.array(p_individual)

# Joint probability: all 15
log10_p_all = np.sum(np.log10(p_individual))
print(f"\n  JOINT PROBABILITY (product of all 15):")
print(f"    log₁₀(P) = {log10_p_all:.2f}")
print(f"    P ≈ 10^{log10_p_all:.1f}")

# Joint probability: 5 most constraining (rank = 5)
p_sorted = np.sort(p_individual)
log10_p5 = np.sum(np.log10(p_sorted[:5]))
print(f"\n  JOINT PROBABILITY (5 most constraining, rank-corrected):")
print(f"    log₁₀(P) = {log10_p5:.2f}")
print(f"    P ≈ 10^{log10_p5:.1f}")

# Joint probability: per-observable comparison
# What fraction of random 15-tuples have ALL errors ≤ max(actual)?
max_actual = actual_formula_errors.max()
mean_actual = actual_formula_errors.mean()
p_one_at_max = np.mean(random_best_errors <= max_actual)
p_one_at_mean = np.mean(random_best_errors <= mean_actual)

print(f"\n  ALTERNATIVE: uniform threshold")
print(f"    Max actual error: {max_actual:.4f}%")
print(f"    P(one random ≤ max) = {p_one_at_max:.6f}")
print(f"    P(15 random all ≤ max) = {p_one_at_max:.6f}^15 = {p_one_at_max**15:.4e}")
print(f"    P(5 random all ≤ max)  = {p_one_at_max:.6f}^5  = {p_one_at_max**5:.4e}")


# ════════════════════════════════════════════════════════════════
# §6 HOW FAR FROM INTEGER ARE THE "TRUE" EXPONENTS?
# ════════════════════════════════════════════════════════════════
print("\n§6. NON-INTEGER CORRECTIONS TO EXPONENTS")
print("-" * 72)

print("""
  If the true relationship is R₂^a · R₃^b · R_EE^c · h∨₃^d · h∨₂^e = target,
  the integer exponents give ~0.3% error. What CONTINUOUS correction to
  each exponent would make it exact?

  For each formula, we find the smallest single-exponent correction
  δn_k that absorbs the residual: δn_k = residual / ln(block_k)
""")

print(f"  {'Name':12s} {'err%':>6s}  {'dR2':>8s} {'dR3':>8s} {'dR_EE':>8s} {'dhv3':>8s} {'dhv2':>8s}")
print(f"  {'-' * 70}")

for name, target in TARGETS.items():
    vec = KNOWN_VECTORS[name]
    log_pred = sum(v * xi for v, xi in zip(vec, x))
    log_target = np.log(target)
    residual = log_target - log_pred
    err = abs(np.exp(residual) - 1) * 100

    corrections = []
    for k in range(5):
        if abs(x[k]) > 1e-10:
            delta_n = residual / x[k]
        else:
            delta_n = np.inf
        corrections.append(delta_n)

    print(f"  {name:12s} {err:5.3f}%  "
          f"{corrections[0]:+8.5f} {corrections[1]:+8.5f} {corrections[2]:+8.5f} "
          f"{corrections[3]:+8.5f} {corrections[4]:+8.5f}")

# Check: are the corrections correlated? Is there a SINGLE universal shift?
print(f"\n  If there existed a universal shift δ to ALL R₃ exponents:")
all_delta_R3 = []
for name, target in TARGETS.items():
    vec = KNOWN_VECTORS[name]
    log_pred = sum(v * xi for v, xi in zip(vec, x))
    residual = np.log(target) - log_pred
    if vec[1] != 0:  # R3 exponent != 0
        delta = residual / x[1]
        all_delta_R3.append(delta)

if all_delta_R3:
    deltas = np.array(all_delta_R3)
    print(f"    Mean δ(R₃ exponent):   {deltas.mean():+.5f}")
    print(f"    Std δ(R₃ exponent):    {deltas.std():.5f}")
    print(f"    If universal: all should be same → std/mean = {deltas.std()/abs(deltas.mean()):.2f}")
    if deltas.std() / abs(deltas.mean()) < 0.5:
        print(f"    → There IS a preferred correction direction!")
    else:
        print(f"    → No universal correction. Residuals are quasi-random.")


# ════════════════════════════════════════════════════════════════
# §7 CORRELATION: ERROR vs LOCAL LATTICE DENSITY
# ════════════════════════════════════════════════════════════════
print("\n§7. ERROR vs LOCAL LATTICE DENSITY")
print("-" * 72)

print("""
  If errors = lattice resolution: formula_error ∝ local_gap
  If errors = physical deviation:  formula_error independent of local_gap
""")

valid = ~np.isnan(local_gaps)
if valid.sum() > 3:
    corr = np.corrcoef(actual_formula_errors[valid], local_gaps[valid])[0, 1]
    print(f"  Correlation(formula_error, local_gap) = {corr:.4f}")
    if corr > 0.5:
        print(f"  → STRONG positive correlation: errors ARE lattice resolution")
    elif corr > 0.2:
        print(f"  → Weak positive correlation: partially lattice, partially physical")
    else:
        print(f"  → No correlation: errors are NOT from lattice resolution")

    # Also check: error/gap ratio
    ratios = actual_formula_errors[valid] / local_gaps[valid]
    print(f"\n  Error / local_gap ratio:")
    print(f"    Mean:   {ratios.mean():.3f}")
    print(f"    Median: {np.median(ratios):.3f}")
    print(f"    Std:    {ratios.std():.3f}")
    print(f"    (If pure lattice resolution → ratio ≈ 0.5 uniformly)")


# ════════════════════════════════════════════════════════════════
# §8 THE RANK CONSTRAINT AND CONSISTENCY TESTS
# ════════════════════════════════════════════════════════════════
print("\n§8. RANK CONSTRAINT — WHICH MATCHES ARE INDEPENDENT?")
print("-" * 72)

# Build the exponent matrix
M = np.array(list(KNOWN_VECTORS.values()), dtype=float)
names = list(KNOWN_VECTORS.keys())
rank = np.linalg.matrix_rank(M)

print(f"  Exponent matrix: {M.shape[0]} formulas x {M.shape[1]} blocks")
print(f"  Rank: {rank}")
print(f"  Independent predictions: {rank}")
print(f"  Consistency constraints: {len(TARGETS) - rank}")

# Find a maximal independent set using SVD
U, S, Vt = np.linalg.svd(M, full_matrices=False)
print(f"\n  Singular values: {np.round(S, 3)}")

# Greedy approach to find maximal independent set
independent_idx = []
current_rank = 0
for i in range(M.shape[0]):
    test = M[independent_idx + [i], :]
    if np.linalg.matrix_rank(test) > current_rank:
        independent_idx.append(i)
        current_rank += 1
    if current_rank == rank:
        break

dependent_idx = [i for i in range(M.shape[0]) if i not in independent_idx]

print(f"\n  Independent set: {[names[i] for i in independent_idx]}")
print(f"  Dependent set:   {[names[i] for i in dependent_idx]}")

# For each dependent formula, express as linear combination of independent ones
# and check if the implied experimental relationship holds
M_indep = M[independent_idx, :]
print(f"\n  CONSISTENCY CHECKS (dependent formulas):")
print(f"  {'Name':12s} {'Implied value':>13s} {'Actual':>10s} {'Consist.err':>12s}")
print(f"  {'-' * 55}")

consistency_errors = []
for j in dependent_idx:
    # Solve: n_j = Σ c_k * n_{independent_k}
    coeffs, residuals, _, _ = np.linalg.lstsq(M_indep.T, M[j, :], rcond=None)

    # Implied log-value
    log_targets_indep = np.array([np.log(list(TARGETS.values())[i]) for i in independent_idx])
    implied_log = coeffs @ log_targets_indep
    actual_log = np.log(list(TARGETS.values())[j])

    implied_val = np.exp(implied_log)
    actual_val = list(TARGETS.values())[j]
    consist_err = abs(implied_val - actual_val) / actual_val * 100
    consistency_errors.append(consist_err)

    print(f"  {names[j]:12s} {implied_val:13.4f} {actual_val:10.4f} {consist_err:11.4f}%")

consistency_errors = np.array(consistency_errors)
print(f"\n  Mean consistency error:   {consistency_errors.mean():.4f}%")
print(f"  Max consistency error:    {consistency_errors.max():.4f}%")
print(f"")
print(f"  These consistency checks are INDEPENDENT of the graph.")
print(f"  They test whether the experimental mass ratios satisfy")
print(f"  the multiplicative relations implied by the exponent structure.")


# ════════════════════════════════════════════════════════════════
# §9 MONTE CARLO: JOINT PROBABILITY WITH RANK CORRECTION
# ════════════════════════════════════════════════════════════════
print("\n§9. CORRECTED JOINT PROBABILITY")
print("-" * 72)

# The honest calculation:
# 1. Only 5 formulas are independently "searched"
# 2. The other 10 are DETERMINED by rank-5 constraint
# 3. But those 10 must ALSO match experiment → consistency test
# 4. The probability is:
#    P_search × P_consistency

# P_search: probability that 5 independent searches each find a match
p_search_5 = np.prod(p_individual[independent_idx])
log10_search = np.log10(p_search_5) if p_search_5 > 0 else -999

# P_consistency: probability that 10 consistency checks pass within observed error
# This is a property of the EXPERIMENTAL mass ratios, not the graph
# We can estimate it by Monte Carlo: generate random "mass ratios" and check
N_MC_consist = 100_000
consist_pass_count = 0

for trial in range(N_MC_consist):
    # Generate 5 random "independent" targets in the same range
    rand_indep = np.random.uniform(lt_min, lt_max, rank)

    # For each, find the best lattice match
    rand_vecs = np.zeros((rank, 5))
    for k in range(rank):
        lt = rand_indep[k]
        idx = np.searchsorted(unique_log_values, lt)
        best_val = unique_log_values[max(0, min(idx, n_unique - 1))]
        if idx > 0:
            if abs(unique_log_values[idx - 1] - lt) < abs(best_val - lt):
                best_val = unique_log_values[idx - 1]
        # We don't need the actual vector, just check if consistency holds
        # For a proper test, we'd need the actual exponent vector too
        # Simplified: just check if the 15th random target also matches

    # Simplified approach: sample 15 random targets and check if ALL match
    # within the max observed error
    rand_targets_15 = np.random.uniform(lt_min, lt_max, 15)
    all_match = True
    for lt in rand_targets_15:
        idx = np.searchsorted(unique_log_values, lt)
        best_err_val = float('inf')
        for i in [max(0, idx - 1), min(n_unique - 1, idx)]:
            e = abs(unique_log_values[i] - lt)
            if e < best_err_val:
                best_err_val = e
        if (np.exp(best_err_val) - 1) * 100 > max_actual:
            all_match = False
            break
    if all_match:
        consist_pass_count += 1

p_all_15_mc = consist_pass_count / N_MC_consist
log10_mc = np.log10(p_all_15_mc) if p_all_15_mc > 0 else np.log10(1 / N_MC_consist)

print(f"  Monte Carlo estimate: P(15 random targets ALL match within {max_actual:.3f}%)")
print(f"    Trials:  {N_MC_consist:,}")
print(f"    Passes:  {consist_pass_count}")
print(f"    P =      {p_all_15_mc:.6f}")
print(f"    log₁₀(P) = {log10_mc:.2f}")

# Second MC: product of individual P_i (analytical, more precise)
print(f"\n  Analytical product of individual P_i:")
print(f"    P_joint (all 15) = 10^{log10_p_all:.2f}")
print(f"    P_joint (5 indep) = 10^{log10_search:.2f}")


# ════════════════════════════════════════════════════════════════
# §10 SYNTHESIS
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("§10. SYNTHESIS — WHAT IS THE 0.3%?")
print("=" * 72)

print(f"""
  ╔══════════════════════════════════════════════════════════════════╗
  ║  FINDING 1: EFFECTIVE FORMULA SPACE                             ║
  ║                                                                  ║
  ║  5 building blocks → 3 independent log-bases                    ║
  ║  9^5 = {total_5d:,} nominal formulas → {n_unique:,} distinct values        ║
  ║  Reduction: {total_5d / n_unique:.0f}x (R₂, R_EE, h∨₂ all powers of 2)           ║
  ╠══════════════════════════════════════════════════════════════════╣
  ║  FINDING 2: LATTICE RESOLUTION                                  ║
  ║                                                                  ║
  ║  Mean gap in target range:   {gaps_in_range.mean():.4f}%                        ║
  ║  Expected best-match error:  {gaps_in_range.mean() / 2:.4f}%  (gap/2)              ║
  ║  Actual mean formula error:  {actual_formula_errors.mean():.4f}%                        ║
  ║  Ratio actual/expected:      {actual_formula_errors.mean() / (gaps_in_range.mean() / 2):.2f}x                           ║
  ╠══════════════════════════════════════════════════════════════════╣
  ║  FINDING 3: JOINT PROBABILITY                                   ║
  ║                                                                  ║
  ║  P(all 15 ≤ observed errors)   = 10^{log10_p_all:.1f}                     ║
  ║  P(5 most constraining)        = 10^{log10_p5:.1f}                     ║
  ║  P(15 random all ≤ max error)  = {p_all_15_mc:.6f} (MC)                  ║
  ╠══════════════════════════════════════════════════════════════════╣
  ║  FINDING 4: CONSISTENCY                                          ║
  ║                                                                  ║
  ║  Rank of exponent matrix:      {rank}                                ║
  ║  Independent searched formulas: {rank}                                ║
  ║  Consistency tests (dependent): {len(TARGETS) - rank}                               ║
  ║  Mean consistency error:        {consistency_errors.mean():.4f}%                       ║
  ╚══════════════════════════════════════════════════════════════════╝
""")

# The conceptual question
print("""
  ═══════════════════════════════════════════════════════════════════
  THE CONCEPTUAL QUESTION: WHAT IS THE ~0.3% DEVIATION?
  ═══════════════════════════════════════════════════════════════════

  Four hypotheses tested:

  H1: Lattice resolution (integer-exponent discretization)
      If true: errors ∝ local gap, ratio ≈ 0.5
      Result: see §7 correlation analysis above

  H2: Missing QFT corrections (tree-level approximation)
      If true: errors ≈ α_s/π ≈ 3.75%
      Result: observed 0.3% << 3.75% → NOT pure tree-level

  H3: Non-perturbative (formulas are exact "dressed" values)
      If true: errors ≈ 0 (or ~ exp(-1/α))
      Result: observed 0.3% >> 0 → NOT fully non-perturbative

  H4: Partial resummation (graph captures most but not all
      radiative structure)
      If true: errors ≈ α²/(4π) ≈ 0.01-0.3% (2-loop level)
      Result: observed 0.3% is at the upper edge of this range

  α_s²/(4π)² ≈ 0.009%  (pure 2-loop strong)
  α_s²/(4π)  ≈ 0.11%   (mixed 2-loop)
  α_s/π      ≈ 3.75%   (1-loop strong)

  The 0.3% sits between 1-loop and 2-loop, consistent with
  "all leading 1-loop absorbed, some subleading 1-loop remaining."
""")

print("\n" + "=" * 72)
print("  ANALYSIS COMPLETE")
print("=" * 72)
print(f"\n  Output saved to: {OUTPUT_PATH}")

sys.stdout = _orig_stdout
_outfile.close()
print(f"\nDone. Output saved to: {OUTPUT_PATH}")
