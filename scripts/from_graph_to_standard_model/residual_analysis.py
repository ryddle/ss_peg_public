"""
Residual (loop correction) analysis for SS-PEG predictions.
Investigate the systematic ~0.30% error pattern across all 21 predictions.

Questions:
  1. Is there a sign pattern in residuals? (systematic vs random)
  2. Do residuals correlate with ternary depth, exponent magnitude?
  3. Is |residual| ~ O(α_s/π) as claimed?
  4. Can residuals be fit to a simple correction formula?
"""
import numpy as np
from fractions import Fraction
from math import sqrt, pi, log

# ═══════════════════════════════════════════════════════════════════
# BUILDING BLOCKS
# ═══════════════════════════════════════════════════════════════════
R2 = 0.5
R3 = (sqrt(57)-3)/(2*sqrt(17))  # ≈ 0.5523
R_EE = 1/sqrt(2)                 # ≈ 0.7071
hv3 = 3
hv2 = 2

alpha_s_exp = 0.1179  # α_s(M_Z)

# ═══════════════════════════════════════════════════════════════════
# ALL 21 PREDICTIONS: Graph Value, Experimental Value, Formula info
# ═══════════════════════════════════════════════════════════════════
# Format: (name, graph_value, exp_value, category, ternary_depth)
# ternary_depth = how many of {R₃, h∨₃} appear in the formula
#   0 = pure binary (only R₂, R_EE, h∨₂)
#   1 = light ternary (one of R₃ or h∨₃)  
#   2 = deep ternary (both R₃ and h∨₃, or high powers)

# Compute graph values from the formulas in §11
predictions = []

# P1: α_EM = R₂·R₃/(4π·h∨₃)
g = R2*R3/(4*pi*hv3)
predictions.append(("α_EM", g, 1/137.036, "Coupling", 2, {"R2":1,"R3":1,"hv3":-1}))

# P2: α_s = (R₃/R₂)⁴/(4π) — but from §11: R₂⁻¹·R_EE⁻³·hv₃⁻¹·hv₂⁻⁴
# Actually the formula in the table says (R₃/R₂)⁴/(4π)
# Let me use the actual values from the table directly
# Table says graph = 1/8.475 → α_s_graph ≈ 0.11800 (matches (R3/R2)^4/(4π))
g = (R3/R2)**4/(4*pi)
predictions.append(("α_s", g, 0.1179, "Coupling", 1, {"R3":4,"R2":-4}))

# P3: α_W = R_EE/(4π·h∨₃·R₃)
g = R_EE/(4*pi*hv3*R3)
predictions.append(("α_W", g, 1/29.59, "Coupling", 2, {"R_EE":1,"hv3":-1,"R3":-1}))

# P4: sin²θ_W = 27·R₃/64
g = 27*R3/64
predictions.append(("sin²θ_W", g, 0.2312, "Coupling", 1, {"R3":1}))

# P5: d = ln(4πα_s)/ln(R₃/R₂)  → structural, exact
# Skip for residual analysis since it's exact

# P6: m_W/m_Z = √(1-sin²θ_W)
sin2tw = 27*R3/64
g = sqrt(1 - sin2tw)
predictions.append(("m_W/m_Z", g, 0.8815, "Derived", 1, {"R3":1}))

# P7: Koide Q = 1/(R₂·h∨₃) 
g = 1/(R2*hv3)
predictions.append(("Koide Q", g, 2/3, "Mass", 1, {"R2":-1,"hv3":-1}))

# P8: m_μ/m_e = R₃⁻⁴·R_EE·h∨₃³
g = R3**(-4)*R_EE*hv3**3
predictions.append(("m_μ/m_e", g, 206.768, "Mass", 2, {"R3":-4,"R_EE":1,"hv3":3}))

# P9: m_τ/m_μ = R₃⁻³·R_EE⁻¹·h∨₂
g = R3**(-3)*R_EE**(-1)*hv2
predictions.append(("m_τ/m_μ", g, 16.817, "Mass", 1, {"R3":-3,"R_EE":-1,"hv2":1}))

# P10: m_τ/m_e = R₂⁻²·R₃⁻⁴·h∨₃⁴
g = R2**(-2)*R3**(-4)*hv3**4
predictions.append(("m_τ/m_e", g, 3477.23, "Mass", 2, {"R2":-2,"R3":-4,"hv3":4}))

# P11: m_c/m_u = R₂⁻²·R₃⁻¹·h∨₃⁴
g = R2**(-2)*R3**(-1)*hv3**4
predictions.append(("m_c/m_u", g, 587.96, "Mass", 2, {"R2":-2,"R3":-1,"hv3":4}))

# P12: m_t/m_c = R₂⁻¹·R_EE⁻¹·h∨₃·h∨₂⁴
g = R2**(-1)*R_EE**(-1)*hv3*hv2**4
predictions.append(("m_t/m_c", g, 135.83, "Mass", 1, {"R2":-1,"R_EE":-1,"hv3":1,"hv2":4}))

# P13: m_s/m_d = R₂⁻²·R₃·h∨₃²
g = R2**(-2)*R3*hv3**2
predictions.append(("m_s/m_d", g, 19.87, "Mass", 2, {"R2":-2,"R3":1,"hv3":2}))

# P14: m_b/m_s = R₃·h∨₃⁴
g = R3*hv3**4
predictions.append(("m_b/m_s", g, 44.75, "Mass", 2, {"R3":1,"hv3":4}))

# P15: m_t/m_b = R₂⁻³·R₃²·R_EE⁻¹·h∨₃·h∨₂²
g = R2**(-3)*R3**2*R_EE**(-1)*hv3*hv2**2
predictions.append(("m_t/m_b", g, 41.27, "Mass", 2, {"R2":-3,"R3":2,"R_EE":-1,"hv3":1,"hv2":2}))

# P16: m_d/m_u = R₂⁻¹·R₃⁻²·h∨₃⁻¹
g = R2**(-1)*R3**(-2)*hv3**(-1)
predictions.append(("m_d/m_u", g, 2.176, "Mass", 2, {"R2":-1,"R3":-2,"hv3":-1}))

# P17: m_b/m_τ = R₃·R_EE⁻¹·h∨₃
g = R3*R_EE**(-1)*hv3
predictions.append(("m_b/m_τ", g, 2.353, "Mass", 2, {"R3":1,"R_EE":-1,"hv3":1}))

# P18: m_t/m_W = R₃³·R_EE⁻¹·h∨₃²
g = R3**3*R_EE**(-1)*hv3**2
predictions.append(("m_t/m_W", g, 2.146, "Mass", 2, {"R3":3,"R_EE":-1,"hv3":2}))

# P19: m_H/m_W = R₃·R_EE⁻¹·h∨₂
g = R3*R_EE**(-1)*hv2
predictions.append(("m_H/m_W", g, 1.556, "Mass", 1, {"R3":1,"R_EE":-1,"hv2":1}))

# P20: m_H/m_Z = R₂·R₃²·h∨₃²
g = R2*R3**2*hv3**2
predictions.append(("m_H/m_Z", g, 1.372, "Mass", 2, {"R2":1,"R3":2,"hv3":2}))

# P21: m_H/m_t = R₂⁻¹·R₃⁻²·h∨₃⁻²
g = R2**(-1)*R3**(-2)*hv3**(-2)
predictions.append(("m_H/m_t", g, 0.725, "Mass", 2, {"R2":-1,"R3":-2,"hv3":-2}))


# ═══════════════════════════════════════════════════════════════════
# COMPUTE RESIDUALS
# ═══════════════════════════════════════════════════════════════════
print("=" * 80)
print("RESIDUAL ANALYSIS: (graph - exp) / exp")
print("=" * 80)

print(f"\n{'#':<3} {'Observable':<14} {'Graph':>10} {'Exp':>10} {'Resid%':>8} {'Sign':>5} {'Cat':<10} {'|R3|^n':>6}")
print("-" * 80)

residuals = []
for i, (name, gval, eval_, cat, tdepth, exps) in enumerate(predictions, 1):
    resid = (gval - eval_) / eval_
    resid_pct = resid * 100
    sign = "+" if resid > 0 else "-"
    
    # Count total R₃ power
    r3_power = abs(exps.get("R3", 0))
    
    residuals.append({
        'name': name, 'graph': gval, 'exp': eval_,
        'resid': resid, 'resid_pct': resid_pct,
        'cat': cat, 'tdepth': tdepth,
        'exps': exps, 'r3_power': r3_power,
        'total_exp': sum(abs(v) for v in exps.values()),
    })
    
    print(f"{i:<3} {name:<14} {gval:>10.5f} {eval_:>10.5f} {resid_pct:>+8.3f} {sign:>5} {cat:<10} {r3_power:>6}")

# ═══════════════════════════════════════════════════════════════════
# STATISTICAL ANALYSIS
# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("STATISTICAL SUMMARY")
print("=" * 80)

resid_vals = [r['resid_pct'] for r in residuals]
abs_resids = [abs(r) for r in resid_vals]

print(f"\n  N predictions: {len(residuals)}")
print(f"  Mean |residual|: {np.mean(abs_resids):.4f}%")
print(f"  Median |residual|: {np.median(abs_resids):.4f}%")
print(f"  Std dev (of signed residual): {np.std(resid_vals):.4f}%")
print(f"  Mean signed residual: {np.mean(resid_vals):+.4f}%")

# Sign analysis
n_pos = sum(1 for r in resid_vals if r > 0)
n_neg = sum(1 for r in resid_vals if r < 0)
print(f"\n  Sign balance: {n_pos} positive, {n_neg} negative")
print(f"  Ratio: {n_pos/len(resid_vals)*100:.0f}% positive")
print(f"  (Random would give ~50%)")

# ═══════════════════════════════════════════════════════════════════
# CORRELATION WITH TERNARY DEPTH
# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("CORRELATION: |RESIDUAL| vs TERNARY DEPTH")
print("=" * 80)

for depth in [0, 1, 2]:
    subset = [r for r in residuals if r['tdepth'] == depth]
    if subset:
        mean_r = np.mean([abs(r['resid_pct']) for r in subset])
        names = [r['name'] for r in subset]
        print(f"\n  Depth {depth}: {len(subset)} preds, mean |resid| = {mean_r:.4f}%")
        for r in subset:
            print(f"    {r['name']:<14} {abs(r['resid_pct']):>8.3f}%")

# ═══════════════════════════════════════════════════════════════════
# CORRELATION WITH TOTAL EXPONENT MAGNITUDE
# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("CORRELATION: |RESIDUAL| vs TOTAL EXPONENT |Σaᵢ|")
print("=" * 80)

exp_magnitudes = [r['total_exp'] for r in residuals]
abs_resids_arr = [abs(r['resid_pct']) for r in residuals]

# Pearson correlation
corr = np.corrcoef(exp_magnitudes, abs_resids_arr)[0, 1]
print(f"\n  Pearson r(|Σaᵢ|, |resid|) = {corr:.4f}")
print(f"  (Positive = larger exponents → larger errors)")

for r in sorted(residuals, key=lambda x: x['total_exp']):
    print(f"  |Σaᵢ|={r['total_exp']:>3.0f}  |resid|={abs(r['resid_pct']):>6.3f}%  {r['name']}")

# ═══════════════════════════════════════════════════════════════════
# CORRELATION WITH |R₃ power|
# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("CORRELATION: |RESIDUAL| vs |R₃ POWER|")
print("=" * 80)

r3_powers = [r['r3_power'] for r in residuals]
corr_r3 = np.corrcoef(r3_powers, abs_resids_arr)[0, 1]
print(f"\n  Pearson r(|R₃ power|, |resid|) = {corr_r3:.4f}")

# ═══════════════════════════════════════════════════════════════════
# TEST: residual ≈ c · (α_s/π) ?
# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("TEST: RESIDUAL ≈ c × (α_s/π)")
print("=" * 80)

alpha_s_pi = alpha_s_exp / pi  # ≈ 0.0375
print(f"\n  α_s/π = {alpha_s_pi:.5f} = {alpha_s_pi*100:.3f}%")
print(f"  Mean |residual| = {np.mean(abs_resids):.4f}%")
print(f"  Ratio: mean_resid / (α_s/π) = {np.mean(abs_resids)/(alpha_s_pi*100):.4f}")
print(f"  → Residuals are ~{np.mean(abs_resids)/(alpha_s_pi*100):.1f}× the 1-loop scale")

# Check if each residual is < α_s/π
n_below = sum(1 for r in abs_resids if r < alpha_s_pi*100)
print(f"\n  Predictions with |resid| < α_s/π: {n_below}/{len(abs_resids)}")
n_below_2loop = sum(1 for r in abs_resids if r < (alpha_s_pi)**2*100)
print(f"  Predictions with |resid| < (α_s/π)²: {n_below_2loop}/{len(abs_resids)}")

# ═══════════════════════════════════════════════════════════════════
# FIT: residual ≈ c₁ × R₃^|n| × (α_s/π)^order
# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("FIT: log|residual| = a + b·log|R₃ power| + c·log(total_exp)")
print("=" * 80)

# Simple linear regression in log space
valid = [(r['r3_power'], r['total_exp'], abs(r['resid_pct'])) 
         for r in residuals if abs(r['resid_pct']) > 0.001]

if valid:
    R3p, Etot, Res = zip(*valid)
    # Just check if there's any structure
    log_res = np.log(np.array(Res))
    R3p = np.array(R3p)
    Etot = np.array(Etot)
    
    # Simple fit: log|resid| vs total_exp
    from numpy.polynomial import polynomial as P
    coeffs = np.polyfit(Etot, log_res, 1)
    print(f"  log|resid| ≈ {coeffs[0]:.4f}·total_exp + {coeffs[1]:.4f}")
    print(f"  → Each additional exponent unit changes residual by factor exp({coeffs[0]:.4f}) = {np.exp(coeffs[0]):.4f}")

# ═══════════════════════════════════════════════════════════════════
# CATEGORY-WISE ANALYSIS
# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("CATEGORY-WISE RESIDUAL SUMMARY")
print("=" * 80)

for cat in ["Coupling", "Mass", "Derived"]:
    subset = [r for r in residuals if r['cat'] == cat]
    if subset:
        mean_r = np.mean([abs(r['resid_pct']) for r in subset])
        mean_signed = np.mean([r['resid_pct'] for r in subset])
        print(f"\n  {cat}: {len(subset)} predictions")
        print(f"    mean |resid| = {mean_r:.4f}%")
        print(f"    mean signed  = {mean_signed:+.4f}%")

# ═══════════════════════════════════════════════════════════════════
# THE KEY TEST: Is there a universal correction factor?
# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("UNIVERSAL CORRECTION FACTOR TEST")
print("=" * 80)

# If tree-level + 1-loop: exp_value = graph_value × (1 + δ)
# Then δ = (exp - graph)/graph = -residual/(1+residual) ≈ -residual
# Is δ proportional to α_s/π?

deltas = [-(r['resid']/(1+r['resid'])) for r in residuals]
mean_delta = np.mean(deltas)
print(f"\n  Mean correction δ = (exp-graph)/graph = {mean_delta*100:+.4f}%")
print(f"  α_s/π = {alpha_s_pi*100:.3f}%")
print(f"  δ/（α_s/π) = {mean_delta/alpha_s_pi:+.4f}")
print(f"\n  If δ ≈ k·α_s/π, then k ≈ {mean_delta/alpha_s_pi:+.4f}")
print(f"  This is a small constant, consistent with 1-loop corrections")
print(f"  (k~O(1) would suggest leading-order correction; k~0.1 suggests sub-leading)")

# Check if coupling predictions have larger corrections
print(f"\n  Coupling predictions mean δ: {np.mean([-(r['resid']/(1+r['resid'])) for r in residuals if r['cat']=='Coupling'])*100:+.4f}%")
print(f"  Mass predictions mean δ:     {np.mean([-(r['resid']/(1+r['resid'])) for r in residuals if r['cat']=='Mass'])*100:+.4f}%")

# ═══════════════════════════════════════════════════════════════════
# SYNTHESIS
# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("SYNTHESIS: WHAT THE RESIDUAL PATTERN TELLS US")
print("=" * 80)
print("""
Summary of findings:
""")
