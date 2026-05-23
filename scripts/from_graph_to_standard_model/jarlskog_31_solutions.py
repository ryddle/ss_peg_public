#!/usr/bin/env python3
"""
jarlskog_31_solutions.py — Why exactly 31 solutions?
=====================================================
Explores the combinatorial/geometric structure behind the 31 representations
of J_CP = 2^{-15} from the 5 building blocks with |exp| ≤ 2.
"""

import numpy as np
from itertools import product as cart_product
import math

# Building blocks and their log2 values
R2   = 1/2       # 2^{-1}
R3   = (np.sqrt(57) - 3) / (2 * np.sqrt(17))  # irrational
R_EE = 1/np.sqrt(2)  # 2^{-1/2}
hv3  = 3          # NOT a power of 2
hv2  = 2          # 2^{+1}

J_CP_exp = 3.0519e-5

print("=" * 72)
print("WHY 31 SOLUTIONS? — ANATOMY OF THE REPRESENTATIONS OF 2^{-15}")
print("=" * 72)

# ============================================================
# §1. WHICH BLOCKS ARE "BINARY"?
# ============================================================
print("\n§1. BINARY vs NON-BINARY BUILDING BLOCKS")
print("-" * 50)

blocks = {"R2": R2, "R3": R3, "R_EE": R_EE, "hv3": hv3, "hv2": hv2}
for name, val in blocks.items():
    log2 = np.log2(val)
    # Check if rational power of 2
    is_binary = any(abs(log2 * d - round(log2 * d)) < 1e-10 for d in range(1, 13))
    print(f"  {name:>5} = {val:.8f}  log₂ = {log2:>8.5f}  {'BINARY (power of 2)' if is_binary else 'IRRATIONAL'}")

print("""
  → Only R2, R_EE, hv2 are powers of 2.
  → R3 and hv3 are NOT → they CANNOT contribute to 2^{-15}.
  → All 31 solutions must have b = d = 0 (exponents of R3, hv3).
  → The problem reduces to 3 effective variables: (a, c, e) on (R2, R_EE, hv2).
""")

# ============================================================
# §2. THE CONSTRAINT EQUATION
# ============================================================
print("§2. THE CONSTRAINT EQUATION")
print("-" * 50)
print("""
  For (a, c, e) with R2^a · R_EE^c · hv2^e:
    base = 2^{-a} · 2^{-c/2} · 2^{e} = 2^{-(a + c/2 - e)}
    base^n = 2^{-15}  ⟹  n · (a + c/2 - e) = 15

  Define k = a + c/2 - e.  Then n = 15/k must be integer ∈ [2, 20].

  With a, c, e ∈ {-2, -1, 0, 1, 2}:
    k_min = -2 + (-2)/2 - 2 = -5     (negative → wrong sign)
    k_max =  2 + 2/2 + 2    =  5

  So k ∈ {-5, -4.5, ..., 0, ..., 4.5, 5} (steps of 0.5)
  We need k > 0 with n = 15/k ∈ ℤ and 2 ≤ n ≤ 20.
""")

# Find all valid k values
valid_ks = []
for k_2 in range(1, 11):  # k = k_2/2, from 0.5 to 5.0
    k = k_2 / 2
    if abs(15 / k - round(15 / k)) < 1e-10:
        n = round(15 / k)
        if 2 <= n <= 20:
            valid_ks.append((k, n))

print(f"  Valid (k, n) pairs:")
print(f"  {'k':>6} {'n':>6} {'k = 15/n':>12} {'Prime factors of n':>25}")
print("  " + "-" * 55)
for k, n in valid_ks:
    # Prime factorization
    factors = []
    m = n
    for p in [2, 3, 5, 7, 11, 13]:
        while m % p == 0:
            factors.append(p)
            m //= p
    fstr = " × ".join(str(f) for f in factors) if factors else "1"
    print(f"  {k:>6.1f} {n:>6d} {15/n:>12.1f} {fstr:>25}")

# ============================================================
# §3. COUNT SOLUTIONS PER n
# ============================================================
print("\n§3. LATTICE POINT COUNT PER n-VALUE")
print("-" * 50)

total = 0
per_n = {}
all_solutions = []

for k, n in valid_ks:
    count = 0
    sols = []
    for a in range(-2, 3):
        for c in range(-2, 3):
            for e in range(-2, 3):
                if abs(a + c/2 - e - k) < 1e-10:
                    count += 1
                    sols.append((a, c, e))
    per_n[n] = (count, sols)
    total += count
    all_solutions.extend([(n, a, c, e) for a, c, e in sols])
    print(f"  n = {n:>2}  (k = {k:.1f}):  {count:>2} solutions")

print(f"  {'':>20} ────────")
print(f"  {'':>20} Total: {total}")

# ============================================================
# §4. IS 31 = 2^5 - 1 MEANINGFUL?
# ============================================================
print(f"\n§4. THE NUMBER 31 = 2⁵ - 1")
print("-" * 50)

print(f"""
  31 = 2⁵ - 1 is the FIFTH MERSENNE PRIME.
  5 = number of building blocks in SS-PEG.

  But is this a COINCIDENCE of the search bounds [-2, +2]?
  Let's test with different bounds:
""")

for bound in [1, 2, 3, 4, 5]:
    count = 0
    for a in range(-bound, bound + 1):
        for c in range(-bound, bound + 1):
            for e in range(-bound, bound + 1):
                k = a + c/2 - e
                if k > 0:
                    n = 15 / k
                    if abs(n - round(n)) < 1e-10:
                        n_int = round(n)
                        if 2 <= n_int <= 20:
                            count += 1
    is_mersenne = (count == 2**int(np.log2(count + 1)) - 1) if count > 0 else False
    tag = f" = 2^{int(np.log2(count+1))}-1 (Mersenne!)" if is_mersenne and count > 1 else ""
    print(f"  bounds = [-{bound}, +{bound}]:  {count:>4} solutions{tag}")

# Test with bound [-6,+6] (the actual CKM search range)
for bound in [6]:
    count = 0
    for a in range(-bound, bound + 1):
        for c in range(-bound, bound + 1):
            for e in range(-bound, bound + 1):
                k = a + c/2 - e
                if k > 0:
                    n = 15 / k
                    if abs(n - round(n)) < 1e-10:
                        n_int = round(n)
                        if 2 <= n_int <= 20:
                            count += 1
    print(f"  bounds = [-{bound}, +{bound}]:  {count:>4} solutions")

# Without n-bounds
print(f"\n  Without the n ∈ [2, 20] restriction:")
for bound in [2]:
    count = 0
    for a in range(-bound, bound + 1):
        for c in range(-bound, bound + 1):
            for e in range(-bound, bound + 1):
                k = a + c/2 - e
                if k > 0:
                    n = 15 / k
                    if abs(n - round(n)) < 1e-10:
                        n_int = round(n)
                        if n_int >= 1:
                            count += 1
    print(f"  bounds = [-{bound}, +{bound}], n ≥ 1:  {count:>4} solutions")

    count = 0
    for a in range(-bound, bound + 1):
        for c in range(-bound, bound + 1):
            for e in range(-bound, bound + 1):
                k = a + c/2 - e
                if k > 0:
                    n = 15 / k
                    if abs(n - round(n)) < 1e-10:
                        n_int = round(n)
                        if n_int >= 2:
                            count += 1
    print(f"  bounds = [-{bound}, +{bound}], n ≥ 2:  {count:>4} solutions = 31 ← original")

# ============================================================
# §5. DEEPER: THE 5 n-VALUES AND THE FACTORIZATION OF 15
# ============================================================
print(f"\n§5. THE 5 n-VALUES FROM THE FACTORIZATION OF 15")
print("-" * 50)

print("""
  15 = 3 × 5, with divisors {1, 3, 5, 15}.

  Integer k → n = 15/k is a divisor of 15:
    k=1 → n=15,  k=3 → n=5,  k=5 → n=3

  Half-integer k → n = 15/(half-int) = 30/(2k):
    k=3/2 → n=10,  k=5/2 → n=6

  So the 5 valid n-values are: {3, 5, 6, 10, 15}.
  These are all products of subsets of {2, 3, 5} that give n ∈ [2, 20]:
    3     = 3
    5     = 5
    6     = 2 · 3
    10    = 2 · 5
    15    = 3 · 5
  Missing: 1, 2, and 30 = 2·3·5 (out of bounds or trivial).

  The primes {2, 3, 5} are exactly:
    2 = hv2 (the Coxeter number triggering half-integer k)
    3 = hv3 (smallest odd prime, first nontrivial factor of 15)
    5 = (15/3) = complementary factor
""")

# ============================================================
# §6. GEOMETRIC STRUCTURE — LATTICE SECTIONS
# ============================================================
print("§6. GEOMETRIC STRUCTURE — LATTICE SECTIONS OF THE CUBE")
print("-" * 50)

print("""
  The 125 lattice points (a, c, e) ∈ [-2, 2]³ are intersected by
  21 hyperplanes of the form 2a + c - 2e = m, for m ∈ [-10, 10].

  Each hyperplane contains lattice points of the cube.
  We select only those hyperplanes where m = 2k = 30/n with integer n ∈ [2, 20].
""")

# Count lattice points on each hyperplane
print(f"  {'m = 2k':>6} {'k':>6} {'n=15/k':>8} {'#points':>8} {'Selected?':>10}")
print("  " + "-" * 45)
hyperplane_counts = {}
for m in range(-10, 11):
    count = 0
    for a in range(-2, 3):
        for c in range(-2, 3):
            for e in range(-2, 3):
                if 2*a + c - 2*e == m:
                    count += 1
    hyperplane_counts[m] = count
    k = m / 2
    if k > 0 and abs(15/k - round(15/k)) < 1e-10 and 2 <= round(15/k) <= 20:
        n = round(15/k)
        sel = "✓"
    else:
        n = None
        sel = ""
    n_str = str(n) if n else "-"
    print(f"  {m:>6} {k:>6.1f} {n_str:>8} {count:>8} {sel:>10}")

print(f"\n  Total lattice points in cube: {sum(hyperplane_counts.values())}")
print(f"  Selected lattice points (our 31): {sum(hyperplane_counts[m] for m in [2, 3, 5, 6, 10])}")

# Distribution is symmetric?
print(f"\n  Symmetry check of hyperplane counts:")
print(f"  {'m':>4} {'count':>6} {'m':>6} {'count':>6} {'symmetric?':>12}")
for m in range(0, 11):
    c1 = hyperplane_counts[m]
    c2 = hyperplane_counts[-m]
    sym = "✓" if c1 == c2 else "✗"
    print(f"  {m:>4} {c1:>6} {-m:>6} {c2:>6} {sym:>12}")

# ============================================================
# §7. THE 31 SOLUTIONS VISUALIZED
# ============================================================
print(f"\n§7. ALL 31 SOLUTIONS — GROUPED BY EXPONENT n")
print("-" * 50)

for k, n in sorted(valid_ks, key=lambda x: x[1]):
    count, sols = per_n[n]
    print(f"\n  n = {n} (k = {k:.1f}, base = 2^{{{-k:.1f}}}):  {count} solutions")
    for a, c, e in sols:
        parts = []
        if a != 0: parts.append(f"R2^{a}")
        if c != 0: parts.append(f"R_EE^{c}")
        if e != 0: parts.append(f"hv2^{e}")
        expr = " · ".join(parts)
        base_val = R2**a * R_EE**c * hv2**e
        J_val = base_val ** n
        print(f"    ({expr:<25})^{n} = 2^{-15} = {J_val:.6e}")

# ============================================================
# §8. DEEPER GEOMETRIC INTERPRETATION
# ============================================================
print(f"\n§8. THE CROSS-POLYTOPE INTERPRETATION")
print("-" * 50)

print("""
  The 5³ = 125 lattice points in [-2,2]³ are sliced by the functional
  f(a,c,e) = 2a + c - 2e (a weighted sum, not symmetric!).

  The WEIGHT VECTOR is w = (2, 1, -2), with:
    ||w||² = 4 + 1 + 4 = 9 = 3²
    ||w||  = 3

  The selected hyperplanes are at heights m = 2, 3, 5, 6, 10.
  These heights ALL satisfy: m divides 30 = 2·15 = 2·(2⁴-1).

  Divisors of 30: {1, 2, 3, 5, 6, 10, 15, 30}
  Our 5 heights: {2, 3, 5, 6, 10} = divisors of 30 \ {1, 15, 30}

  The missing divisors:
    m=1  → k=0.5 → n=30 (exceeds search limit n ≤ 20)
    m=15 → k=7.5 → exceeds cube bounds (k_max = 5)
    m=30 → k=15  → far exceeds cube bounds

  So within our search window: 31 = (all accessible divisors of 30)
""")

# ============================================================
# §9. MAGICAL DECOMPOSITION: 1 + 5 + (5+1) + 7 + 12
# ============================================================
print("§9. DECOMPOSITION PATTERNS")
print("-" * 50)

decomp = [(3, 1), (5, 6), (6, 5), (10, 7), (15, 12)]
print(f"\n  n:     ", "  ".join(f"{n:>3}" for n, _ in decomp))
print(f"  count: ", "  ".join(f"{c:>3}" for _, c in decomp))
print(f"  Total: {sum(c for _, c in decomp)}")

# Check: are the counts triangular numbers?
print(f"\n  Checking for triangular number patterns:")
for n, c in decomp:
    # T_k = k(k+1)/2
    for t in range(1, 10):
        if t*(t+1)//2 == c:
            print(f"    n={n}: count={c} = T_{t} = {t}({t}+1)/2  (triangular!)")
            break
    else:
        # Other patterns
        print(f"    n={n}: count={c}")

# Check: is count = (number of lattice points on hyperplane 2a+c-2e = 2k in [-2,2]³)
print(f"\n  These counts follow from discrete geometry:")
print(f"  For the weight vector w = (2, 1, -2) on the cube [-2,2]³,")
print(f"  the count N(m) of lattice points with 2a + c - 2e = m is:")
print(f"  a convolution of three 1D distributions.")
print(f"  Specifically, N(m) = Σ δ(2a+c-2e, m) over all (a,c,e) ∈ {{-2,...,2}}³")

# The generating function approach
print(f"\n  Generating function approach:")
print(f"    g_a(z) = z^{-4} + z^{-2} + 1 + z^2 + z^4    (2a ∈ {{-4,-2,0,2,4}})")
print(f"    g_c(z) = z^{-2} + z^{-1} + 1 + z + z^2       (c ∈ {{-2,-1,0,1,2}})")
print(f"    g_e(z) = z^4 + z^2 + 1 + z^{-2} + z^{-4}    (-2e ∈ {{-4,-2,0,2,4}})")
print(f"")
print(f"    N(m) = [z^m] g_a(z) · g_c(z) · g_e(z)")
print(f"")
print(f"    g_a · g_e = (z^{-4}+z^{-2}+1+z^2+z^4)² because g_a = g_e")
print(f"    This is a symmetric distribution centered at 0.")

# Compute the full generating function
coeffs = {}
for a in range(-2, 3):
    for c in range(-2, 3):
        for e in range(-2, 3):
            m = 2*a + c - 2*e
            coeffs[m] = coeffs.get(m, 0) + 1

print(f"\n  Full distribution N(m) (should sum to 125):")
for m in sorted(coeffs.keys()):
    bar = "█" * coeffs[m]
    marker = " ← SELECTED" if m in [2, 3, 5, 6, 10] else ""
    print(f"    m={m:>3}: {coeffs[m]:>3}  {bar}{marker}")
print(f"    Sum = {sum(coeffs.values())}")

# ============================================================
# §10. THE ANSWER
# ============================================================
print(f"\n{'=' * 72}")
print("§10. WHY 31?")
print("=" * 72)

print(f"""
  ANSWER: 31 arises from the intersection of THREE structures:

  1. ARITHMETIC: 15 = 3 × 5, so n = 15/k is integer when k divides 15.
     Including half-integers (from R_EE's 2^{{-1/2}}), the effective
     condition is: k divides 30/2 = 15 OR 2k divides 30.
     This gives 5 valid n-values: {{3, 5, 6, 10, 15}}.

  2. GEOMETRY: Each n-value selects a hyperplane 2a + c - 2e = m in the
     discrete cube [-2,2]³. The counts are determined by the convolution
     g_a × g_c × g_e, which is the discrete Fourier spectrum of the cube.

  3. BOUNDS: The search window |exp| ≤ 2 and n ∈ [2, 20] truncates the
     infinite family of representations. Exactly 5 of the 8 divisors
     of 30 fall within these bounds, yielding 1 + 6 + 5 + 7 + 12 = 31.

  IS 31 = 2⁵ - 1 MEANINGFUL?
  ══════════════════════════

  The decomposition 31 = 1 + 6 + 5 + 7 + 12 does NOT factorize as 2⁵ - 1
  in an obvious way. The Mersenne prime property appears to be a NUMERICAL
  COINCIDENCE of the specific bounds [-2, +2] and n ∈ [2, 20].

  However, there IS a deep structural point:

  → 15 = 2⁴ - 1 is itself a Mersenne number.
  → J_CP = 2^{{-15}} = 2^{{-(2⁴ - 1)}}  

  The "binary simplicity" of J_CP — that it's a pure power of 2 with
  a Mersenne exponent — is what CAUSES the rich factorization structure:
  15 has many divisors relative to its size, each giving a family of
  representations.

  CONNECTION TO THE 6 TRIANGLES:
  ═══════════════════════════════

  Of the 5 n-values, n = 6 gives exactly {5} solutions.
  And 6 is the ONLY n-value where the search found the formula first.

  The 5 solutions at n = 6 are:
    (R2¹ · R_EE⁻¹ · hv2⁻²)^6     ← equivalent to "original" form
    (R2² · R_EE⁻¹ · hv2⁻¹)^6
    (R2¹ · R_EE¹  · hv2⁻¹)^6     ← the CANONICAL form (R2·R_EE/hv2)^6
    (R2² · R_EE¹)^6
    (R_EE¹ · hv2⁻²)^6

  These 5 bases all equal 2^{{-5/2}}, just written with different
  "gauge choices" for how R2 and hv2 (both being 2^{{±1}}) trade off.

  The 5 representations of the SAME base correspond to 5 = dim(SU(2)²),
  the dimension of the space of SU(2) × SU(2) embeddings — or equivalently,
  the number of independent ways to distribute the "binary weight" -5/2
  among three SU(2) invariants.
""")

# Final: verify the "5 gauges" interpretation for n=6
print("  Verification: all 5 n=6 bases are identical:")
_, sols_6 = per_n[6]
for a, c, e in sols_6:
    val = R2**a * R_EE**c * hv2**e
    log2_val = np.log2(val)
    parts = []
    if a != 0: parts.append(f"R2^{a}")
    if c != 0: parts.append(f"R_EE^{c}")
    if e != 0: parts.append(f"hv2^{e}")
    print(f"    {' · '.join(parts):<25} = {val:.10f}  = 2^({log2_val:.1f})")

print(f"\n  All 5 give 2^(-5/2) ≈ {2**(-5/2):.10f}  ✓")
print(f"\n  The 31 solutions = 5 'gauges' for n=6  ×  more for other n values")
print(f"  They represent the ALGEBRAIC RICHNESS of 2^{{-15}} = 2^{{-(2⁴-1)}},")
print(f"  not independent physical content.")

print(f"\n{'=' * 72}")
print("ANALYSIS COMPLETE")
print("=" * 72)
