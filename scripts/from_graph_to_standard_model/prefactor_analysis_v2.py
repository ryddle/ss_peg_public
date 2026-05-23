"""
DEEPER prefactor analysis: Explore the emergent structure.

Key observation from first pass:
  - 3/9 prefactors are PURE h∨₃/h∨₂:
    sin²θ_W = 27/64 = h∨₃³ · h∨₂⁻⁶
    J_CP     = 9/64  = h∨₃² · h∨₂⁻⁶
    sin θ₂₃(PMNS) = 1/36 = h∨₃⁻² · h∨₂⁻²
  - The other 6 contain primes: 5, 7, 11, 17, 19, 23, 47, 61

Hypothesis: 5 = h₂+h₃, 11 = dim(su₂)+dim(su₃), so maybe extra primes
come from SUMS of group invariants rather than products.
"""
import numpy as np
from fractions import Fraction
from math import sqrt

# ──────────────────────────────────────────────────────────────────
# Graph / group constants
# ──────────────────────────────────────────────────────────────────
h2 = 2;  hv2 = 2   # Coxeter / dual Coxeter SU(2)
h3 = 3;  hv3 = 3   # Coxeter / dual Coxeter SU(3)
dim_su2 = 3         # 2^2 - 1
dim_su3 = 8         # 3^2 - 1
rank_su2 = 1
rank_su3 = 2
pos_roots_su2 = 1   # N(N-1)/2 for SU(N): SU(2)=1
pos_roots_su3 = 3   # SU(3)=3
fund_dim_su2 = 2
fund_dim_su3 = 3
adj_dim_su2 = 3
adj_dim_su3 = 8

# Composite numbers from SU(2)×SU(3)
h_sum = h2 + h3                 # 5
dim_sum = dim_su2 + dim_su3     # 11
rank_sum = rank_su2 + rank_su3  # 3
total_positive_roots = pos_roots_su2 + pos_roots_su3  # 4
dim_product = dim_su2 * dim_su3  # 24
weyl_product = 2 * 6            # 12
fund_product = fund_dim_su2 * fund_dim_su3  # 6

print("=" * 80)
print("NATURAL NUMBERS FROM SU(2)×SU(3):")
print("=" * 80)
numbers = {
    'h₂': h2, 'h₃': h3,
    'h₂+h₃': h_sum,
    'dim(su₂)': dim_su2, 'dim(su₃)': dim_su3,
    'dim(su₂)+dim(su₃)': dim_sum,
    'dim(su₂)×dim(su₃)': dim_product,
    'rank_sum': rank_sum,
    'total_pos_roots': total_positive_roots,
    '|W₂|×|W₃|': weyl_product,
    'fund₂×fund₃': fund_product,
    'h₂×h₃': h2*h3,
    'h_sum × dim_sum': h_sum * dim_sum,  # 5×11 = 55
    'dim_sum × h_sum': dim_sum * h_sum,
}
for name, val in numbers.items():
    print(f"  {name} = {val}")

# ──────────────────────────────────────────────────────────────────
# Check: are the "exotic" primes expressible as sums/products of 
# low-level group invariants?
# ──────────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("EXOTIC PRIMES AND THEIR ORIGINS:")
print("=" * 80)

exotic_primes = [5, 7, 11, 17, 19, 23, 47, 61]
# Candidate expressions from group theory
expressions = {
    'h₂+h₃': 2+3,
    'h₃+total_pos_roots': 3+4,
    'dim(su₂)+dim(su₃)': 3+8,
    'fund₂×fund₃+dim(su₃)+fund₃': 6+8+3,
    'h_sum²−h₃': 25-3,
    'dim_product−h_sum': 24-5,
    'dim_sum+h₂×h₃': 11+6,
    'dim_sum×h₂−h₃': 11*2-3,
    'h_sum²+h₂': 25+2,
    'dim_sum²−dim_product×h_sum': 121-120,
    'dim(su₂)²−rank_su₂': 9-1,
    'dim(su₃)+pos_roots_su3²': 8+9,
    'dim(su₃)²−h_sum²−dim_su2': 64-25-3,
    'dim_product−h_sum': 24-5,
    'fund₃!²−h₂': 36-2,
}

# Actually, let me be more systematic: try a+b, a-b, a*b, a*b+c, a*b-c
# for small group invariants
basic = {'h₂': 2, 'h₃': 3, 'dim₂': 3, 'dim₃': 8, 'r₂': 1, 'r₃': 2, 
         'Δ₂': 1, 'Δ₃': 3, 'f₂': 2, 'f₃': 3}

print("\nSystematic search: which sums/differences/products of 2 basic invariants give exotic primes?")
for target in exotic_primes:
    found_exprs = []
    names = list(basic.keys())
    vals = list(basic.values())
    for i in range(len(names)):
        for j in range(len(names)):
            a, b = vals[i], vals[j]
            na, nb = names[i], names[j]
            if a + b == target:
                found_exprs.append(f"{na}+{nb} = {a}+{b}")
            if a - b == target and a > b:
                found_exprs.append(f"{na}-{nb} = {a}-{b}")
            if a * b == target:
                found_exprs.append(f"{na}×{nb} = {a}×{b}")
            if a * b + 1 == target:
                found_exprs.append(f"{na}×{nb}+1")
            if a * b - 1 == target:
                found_exprs.append(f"{na}×{nb}-1")
            # depth 3: a*b+c, a*b-c
            for k in range(len(names)):
                c = vals[k]
                nc = names[k]
                if a * b + c == target:
                    found_exprs.append(f"{na}×{nb}+{nc}")
                if a * b - c == target and a*b > c:
                    found_exprs.append(f"{na}×{nb}-{nc}")
    
    # Remove duplicates
    unique = list(dict.fromkeys(found_exprs))
    print(f"\n  {target}: {'  |  '.join(unique[:8])}")

# ──────────────────────────────────────────────────────────────────
# KEY OBSERVATION: Relations between prefactors
# ──────────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("INTER-PREFACTOR RELATIONS:")
print("=" * 80)

prefactors = {
    'sin2_tW':      Fraction(27, 64),
    'sin_t12_CKM':  Fraction(1, 5),
    'sin_t23_CKM':  Fraction(17, 55),
    'sin_t13_CKM':  Fraction(7, 47),
    'J_CP':         Fraction(9, 64),
    'sin2_t12_PMNS': Fraction(11, 3),
    'sin_t23_PMNS': Fraction(1, 36),
    'sin2_t13_PMNS': Fraction(23, 55),
    'J_PMNS':       Fraction(19, 61),
}

# Check J_CP = sin2_tW / h∨₃ ?
print(f"\nJ_CP / sin²θ_W = {prefactors['J_CP'] / prefactors['sin2_tW']} = {float(prefactors['J_CP'] / prefactors['sin2_tW']):.6f}")
print(f"  1/h∨₃      = {Fraction(1,3)} = {1/3:.6f}")
print(f"  -> J_CP = sin²θ_W / h∨₃ ? {prefactors['J_CP'] == prefactors['sin2_tW'] / 3}")
print(f"  -> J_CP = sin²θ_W × (h∨₃/h∨₂⁶) ? sin²θ_W = h∨₃³/h∨₂⁶, J_CP = h∨₃²/h∨₂⁶")
print(f"  -> CONFIRMED: J_CP_prefactor = sin²θ_W_prefactor / h∨₃")

# Check all pairwise ratios
print("\nAll pairwise ratios (simplified):")
names = list(prefactors.keys())
interesting_ratios = []
for i in range(len(names)):
    for j in range(i+1, len(names)):
        ratio = prefactors[names[i]] / prefactors[names[j]]
        rn, rd = ratio.numerator, ratio.denominator
        # Check if it's simple (small numerator and denominator)
        if rn <= 50 and rd <= 50:
            interesting_ratios.append((names[i], names[j], ratio))
            
for n1, n2, r in interesting_ratios:
    print(f"  {n1} / {n2} = {r}")

# ──────────────────────────────────────────────────────────────────
# The 55 = 5×11 pattern
# ──────────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("THE 55 = (h₂+h₃) × dim(su₂⊕su₃) PATTERN:")
print("=" * 80)
print(f"\n55 = 5 × 11 = (h₂+h₃) × (dim su₂ + dim su₃)")
print(f"Appears in: 17/55 (CKM θ₂₃) and 23/55 (PMNS θ₁₃)")
print(f"  17/55 = 17 / [(h₂+h₃)(dim su₂ + dim su₃)]")
print(f"  23/55 = 23 / [(h₂+h₃)(dim su₂ + dim su₃)]")
print()

# Check if 17 and 23 also have meaning
print("Numerator analysis for the /55 family:")
print(f"  17 = dim₃×h₃ - h_sum×h₂ = {8*3 - 5*2}? NO, = {8*3-5*2}")
# Try more
print(f"  17 = dim₃×h₂ + r₃ - 1 = {8*2+2-1}? = {8*2+2-1}")
print(f"  17 = dim₃×h₃ - dim₃ + h₂ = {8*3-8+2}? = {8*3-8+2}")
print(f"  17 = dim₃ + dim₃ + r₃ - 1 = {8+8+2-1}? = {8+8+2-1}")
print(f"  17 = f₂×dim₃ + rank₃ - 1 = {2*8+2-1}")
print(f"  17 = h₃×h_sum + h₂ = {3*5+2}")
print(f"  23 = h₃×dim₃ - 1 = {3*8-1}")
print(f"  23 = dim₃×h₃ - h₂ + h₃ - 1 = {8*3-2+3-1}? = {8*3-2+3-1}")

# ──────────────────────────────────────────────────────────────────
# Connection to R₃ = (√57-3)/(2√17)
# ──────────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("CONNECTION TO R₃ = (√57−3)/(2√17):")
print("=" * 80)
print(f"  57 = 3 × 19    (h₃ × 19)")
print(f"  17 = prime      (appears in sin_t23 CKM numerator!)")
print(f"  R₃² = (57-6√57+9)/(4×17) = (66-6√57)/(68)")
R3 = (sqrt(57)-3)/(2*sqrt(17))
print(f"  R₃ = {R3:.10f}")
print(f"  R₃² = {R3**2:.10f}")

# Check: does 17/55 involve R₃ indirectly?
print(f"\n  17/55 × 55/17 = 1 (trivial)")
print(f"  17 is already in R₃ = (√57-3)/(2√17)")
print(f"  19 is in 57 = 3×19, and also in J_PMNS = 19/61")
print(f"  47 has no obvious group-theory origin")
print(f"  61 has no obvious group-theory origin")

# ──────────────────────────────────────────────────────────────────
# CLASSIFICATION: Which sector is each prefactor?
# ──────────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("CLASSIFICATION OF PREFACTORS BY STRUCTURE:")
print("=" * 80)

print("""
Tier 1 — PURE COXETER (only h∨₂, h∨₃):
  sin²θ_W   = 27/64 = h∨₃³/h∨₂⁶           "Weinberg = (h₃/h₂²)³"
  J_CP(CKM) = 9/64  = h∨₃²/h∨₂⁶           "J_CP = (h₃/h₂²)² × h₂⁻²"  [= sin²θ_W / h∨₃]
  sinθ₂₃(PMNS) = 1/36 = 1/(h∨₃²·h∨₂²)    "1 over area h₃²×h₂²"

Tier 2 — SU(2)×SU(3) ADDITIVE INVARIANTS (involve h₂+h₃=5 or dim(g)=11):
  sinθ₁₂(CKM)  = 1/5  = 1/(h₂+h₃)        "inverse of total Coxeter number"
  sinθ₂₃(CKM)  = 17/55                     denominator = (h₂+h₃)×dim(g)
  sin²θ₁₂(PMNS) = 11/3 = dim(g)/h₃        "dimension of full algebra / Coxeter h₃"
  sin²θ₁₃(PMNS) = 23/55                    denominator = (h₂+h₃)×dim(g)

Tier 3 — INVOLVE SPECTRAL GRAPH INVARIANT (primes from R₃):
  sinθ₁₃(CKM)  = 7/47                      47 prime (no obvious group origin)
  J_PMNS        = 19/61                      19∈57=3×19 (in R₃); 61 prime
""")

# ──────────────────────────────────────────────────────────────────
# VERIFY Tier 2 interpretations
# ──────────────────────────────────────────────────────────────────
print("=" * 80)
print("VERIFICATION OF TIER 2 INTERPRETATIONS:")
print("=" * 80)
print(f"\n1/5 = 1/(h₂+h₃): {Fraction(1,5) == Fraction(1, h2+h3)}")
print(f"11/3 = dim(g)/h₃: {Fraction(11,3) == Fraction(dim_sum, h3)}")
print(f"17/55 denominator = (h₂+h₃)×dim(g) = {h_sum}×{dim_sum} = {h_sum*dim_sum}: {55 == h_sum*dim_sum}")
print(f"23/55 denominator = (h₂+h₃)×dim(g) = {h_sum}×{dim_sum} = {h_sum*dim_sum}: {55 == h_sum*dim_sum}")

# ──────────────────────────────────────────────────────────────────
# Are 17 and 23 (numerators of the /55 fractions) also expressible?
# ──────────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("NUMERATORS 17 AND 23 — DEEPER SEARCH:")
print("=" * 80)

# Extended set of group numbers
ext = {
    'h₂': 2, 'h₃': 3, 'dim₂': 3, 'dim₃': 8, 'r₂': 1, 'r₃': 2,
    'Δ₂': 1, 'Δ₃': 3, 'f₂': 2, 'f₃': 3, 'dim_g': 11, 'h_sum': 5,
    'r_sum': 3, 'Δ_tot': 4, 'f_prod': 6, '|W|': 12, 'dim_prod': 24,
}

for target in [17, 23, 47, 61, 7, 19]:
    found = []
    ns = list(ext.keys())
    vs = list(ext.values())
    # a*b ± c
    for i in range(len(ns)):
        for j in range(len(ns)):
            for k in range(len(ns)):
                if vs[i]*vs[j] + vs[k] == target:
                    found.append(f"{ns[i]}×{ns[j]}+{ns[k]}")
                if vs[i]*vs[j] - vs[k] == target and vs[i]*vs[j] > vs[k]:
                    found.append(f"{ns[i]}×{ns[j]}-{ns[k]}")
    # a² ± b
    for i in range(len(ns)):
        for j in range(len(ns)):
            if vs[i]**2 + vs[j] == target:
                found.append(f"{ns[i]}²+{ns[j]}")
            if vs[i]**2 - vs[j] == target:
                found.append(f"{ns[i]}²-{ns[j]}")
    # a+b+c
    for i in range(len(ns)):
        for j in range(i, len(ns)):
            for k in range(j, len(ns)):
                if vs[i]+vs[j]+vs[k] == target:
                    found.append(f"{ns[i]}+{ns[j]}+{ns[k]}")
    
    unique = list(dict.fromkeys(found))
    # Filter for most natural / simplest
    print(f"\n  {target}:")
    for u in unique[:10]:
        print(f"    = {u}")

# ──────────────────────────────────────────────────────────────────
# SUMMARY TABLE
# ──────────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("SUMMARY: STRUCTURAL INTERPRETATION OF ALL 9 PREFACTORS")
print("=" * 80)
print(f"""
{'Observable':<20} {'Prefactor':<10} {'Structural form':<40} {'Tier':<6}
{'─'*80}
sin²θ_W              27/64      h∨₃³ · h∨₂⁻⁶                          1
J_CP (CKM)           9/64       h∨₃² · h∨₂⁻⁶ = sin²θ_W/h∨₃           1
sinθ₂₃ (PMNS)        1/36       (h∨₃·h∨₂)⁻²                           1
sinθ₁₂ (CKM)         1/5        1/(h₂+h₃)                             2
sin²θ₁₂ (PMNS)       11/3       dim(su₂⊕su₃)/h₃                      2
sinθ₂₃ (CKM)         17/55      17/[(h₂+h₃)·dim(g)]                   2*
sin²θ₁₃ (PMNS)       23/55      23/[(h₂+h₃)·dim(g)]                   2*
sinθ₁₃ (CKM)         7/47       ?/?                                    3
|J_PMNS|              19/61      ?/?                                    3

* Denominator fully explained; numerators (17, 23) need deeper analysis.
  17 = h₃×h_sum+h₂ = h₃(h₂+h₃)+h₂ = h₃²+h₃h₂+h₂
  23 = dim₃×h₃-Δ₂ = 8×3-1 = dim₃×h∨₃ - pos_roots(SU₂)
""")
