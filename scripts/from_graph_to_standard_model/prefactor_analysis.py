"""
Analyze rational prefactors from SS-PEG tower to determine
if they decompose into products of Coxeter numbers h∨₃=3, h∨₂=2
or involve deeper combinatorial structure.
"""
import numpy as np
from fractions import Fraction
from itertools import product as cprod
from math import comb, factorial

# Building blocks
hv3 = 3  # dual Coxeter SU(3)
hv2 = 2  # dual Coxeter SU(2)
R2  = Fraction(1, 2)
R3_approx = 0.5523  # (sqrt(57)-3)/(2*sqrt(17))

# The rational prefactors from COMPREHENSIVE_RESULTS.md
prefactors = {
    'sin2_tW (Weinberg)':  Fraction(27, 64),
    'sin_t12 (CKM)':      Fraction(1, 5),     # 3/15 reduced
    'sin_t23 (CKM)':      Fraction(17, 55),
    'sin_t13 (CKM)':      Fraction(7, 47),
    'J_CP (CKM)':         Fraction(9, 64),
    'sin2_t12 (PMNS)':    Fraction(11, 3),
    'sin_t23 (PMNS)':     Fraction(1, 36),
    'sin2_t13 (PMNS)':    Fraction(23, 55),
    'J_PMNS':             Fraction(19, 61),
}

def factorize(n):
    """Factorize into primes."""
    if n == 0:
        return {0: 1}
    if n == 1:
        return {}
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = 1
    return factors


print("=" * 80)
print("SECTION 1: PRIME FACTORIZATION OF EACH PREFACTOR")
print("=" * 80)
for name, frac in prefactors.items():
    p, q = frac.numerator, frac.denominator
    pf = factorize(p)
    qf = factorize(q)
    pstr = ' * '.join(f'{k}^{v}' if v > 1 else str(k) for k, v in sorted(pf.items())) or '1'
    qstr = ' * '.join(f'{k}^{v}' if v > 1 else str(k) for k, v in sorted(qf.items())) or '1'
    
    only_23 = all(k in [2, 3] for k in pf) and all(k in [2, 3] for k in qf)
    
    if only_23:
        a3 = pf.get(3, 0) - qf.get(3, 0)
        a2 = pf.get(2, 0) - qf.get(2, 0)
        tag = f"PURE: 3^{a3} * 2^{a2}"
    else:
        extra = set(k for k in list(pf.keys()) + list(qf.keys()) if k not in [2, 3])
        tag = f"Extra primes: {extra}"
    
    print(f"\n{name}:")
    print(f"  {frac} = {pstr} / {qstr}")
    print(f"  -> {tag}")

print("\n" + "=" * 80)
print("SECTION 2: ATTEMPT h∨₃^a * h∨₂^b DECOMPOSITION")
print("=" * 80)
for name, frac in prefactors.items():
    val = float(frac)
    best_err = 1e10
    best_ab = None
    for a in range(-8, 9):
        for b in range(-8, 9):
            cand_val = 3.0**a * 2.0**b
            err = abs(cand_val - val) / max(abs(val), 1e-15)
            if err < best_err:
                best_err = err
                best_ab = (a, b)
    
    cand_val = 3.0**best_ab[0] * 2.0**best_ab[1]
    match = best_err < 1e-10
    status = "EXACT MATCH" if match else f"nearest, err={best_err*100:.2f}%"
    print(f"  {name:<25} {str(frac):<10} -> h3^{best_ab[0]}*h2^{best_ab[1]} = {cand_val:.6f}  [{status}]")

print("\n" + "=" * 80)
print("SECTION 3: REPRESENTATION DIMENSION RATIOS")
print("  dim(fund SU(2)) = 2, dim(adj SU(2)) = 3")
print("  dim(fund SU(3)) = 3, dim(adj SU(3)) = 8")
print("  Casimir SU(2) fund = 3/4, adj = 2")
print("  Casimir SU(3) fund = 4/3, adj = 3")
print("=" * 80)

# Representation dimensions and their ratios
dims = {
    'fund_SU2': 2, 'adj_SU2': 3,
    'fund_SU3': 3, 'adj_SU3': 8,
    'fund_SU2_x_SU3': 6, 'adj_SU2_x_SU3': 24,
}

# Also try Casimir values
casimirs = {
    'C2_fund_SU2': Fraction(3, 4),
    'C2_adj_SU2':  Fraction(2, 1),
    'C2_fund_SU3': Fraction(4, 3),
    'C2_adj_SU3':  Fraction(3, 1),
}

# Try all dim_A / dim_B ratios
print("\nChecking if any prefactor = dim(R_i)/dim(R_j):")
for name, frac in prefactors.items():
    matches = []
    for n1, d1 in dims.items():
        for n2, d2 in dims.items():
            if n1 != n2 and Fraction(d1, d2) == frac:
                matches.append(f"{n1}/{n2} = {d1}/{d2}")
    if matches:
        print(f"  {name}: {frac} = {matches}")
    
    # Also check Casimir ratios
    for n1, c1 in casimirs.items():
        for n2, c2 in casimirs.items():
            if n1 != n2 and c1/c2 == frac:
                matches.append(f"C: {n1}/{n2} = {c1}/{c2}")
    if not matches:
        pass  # Will check below

print("\n" + "=" * 80)
print("SECTION 4: COMBINATORIAL DECOMPOSITION")
print("  Try p/q = C(n,k) / C(m,l) * 3^a * 2^b for small n,k,m,l,a,b")
print("=" * 80)

# Generate all small binomial coefficients
binom_vals = {}
for n in range(1, 13):
    for k in range(0, n+1):
        c = comb(n, k)
        if c not in binom_vals:
            binom_vals[c] = []
        binom_vals[c].append((n, k))

for name, frac in prefactors.items():
    found = []
    # p/q = C(n1,k1)/C(n2,k2) * 3^a * 2^b
    for a in range(-4, 5):
        for b in range(-4, 5):
            target = frac / Fraction(3**max(a,0), 3**max(-a,0)) / Fraction(2**max(b,0), 2**max(-b,0))
            if target.numerator in binom_vals and target.denominator == 1:
                for nk in binom_vals[target.numerator]:
                    found.append(f"C({nk[0]},{nk[1]}) * 3^{a} * 2^{b}")
            elif target.denominator in binom_vals and target.numerator == 1:
                for nk in binom_vals[target.denominator]:
                    found.append(f"3^{a} * 2^{b} / C({nk[0]},{nk[1]})")
            elif target.numerator in binom_vals and target.denominator in binom_vals:
                for nk1 in binom_vals[target.numerator]:
                    for nk2 in binom_vals[target.denominator]:
                        if nk1[0] <= 8 and nk2[0] <= 8:
                            found.append(f"C({nk1[0]},{nk1[1]})/C({nk2[0]},{nk2[1]}) * 3^{a} * 2^{b}")
    
    # Keep only simple ones (short description)
    simple = [f for f in found if len(f) < 40][:5]
    if simple:
        print(f"\n  {name} = {frac}:")
        for s in simple:
            print(f"    {s}")
    else:
        print(f"\n  {name} = {frac}: No simple combinatorial decomposition found")

print("\n" + "=" * 80)
print("SECTION 5: DEEPER PATTERNS - GROUP THEORY NUMBERS")
print("=" * 80)

# SU(N) group theory numbers
# Weyl group order: |W(SU(N))| = N!
# Number of positive roots: N(N-1)/2
# Dimension of Lie algebra: N^2-1
# Index of fundamental rep: 1/2
# Dynkin index: T(fund)=1/2, T(adj)=N
# Center: Z_N

group_numbers = {
    '|W(SU(2))|': 2,      # 2!
    '|W(SU(3))|': 6,      # 3!
    '|W(SU(2)xSU(3))|': 12,  # 2! * 3!
    'pos_roots_SU2': 1,
    'pos_roots_SU3': 3,
    'dim_su2': 3,          # 2^2-1
    'dim_su3': 8,          # 3^2-1
    'rank_SU2': 1,
    'rank_SU3': 2,
    'center_SU2': 2,
    'center_SU3': 3,
    'h_SU2(Coxeter)': 2,
    'h_SU3(Coxeter)': 3,
    'hv_SU2(dual_Cox)': 2,
    'hv_SU3(dual_Cox)': 3,
}

# Try each prefactor as ratio of two group numbers * 3^a * 2^b
print("\nSearching: p/q = (group_num_1 / group_num_2) * 3^a * 2^b")
for name, frac in prefactors.items():
    found = []
    for n1, v1 in group_numbers.items():
        for n2, v2 in group_numbers.items():
            ratio = Fraction(v1, v2)
            rem = frac / ratio
            # Check if remainder is 3^a * 2^b
            if rem == 0:
                continue
            num, den = rem.numerator, rem.denominator
            nf = factorize(num)
            df = factorize(den)
            if all(k in [2, 3] for k in nf) and all(k in [2, 3] for k in df):
                a3 = nf.get(3, 0) - df.get(3, 0)
                a2 = nf.get(2, 0) - df.get(2, 0)
                if abs(a3) <= 4 and abs(a2) <= 4:
                    found.append(f"({n1}/{n2}) * 3^{a3} * 2^{a2}")
    
    if found:
        print(f"\n  {name} = {frac}:")
        for f in found[:5]:
            print(f"    {f}")

print("\n" + "=" * 80)
print("SECTION 6: NUMEROLOGY - WHICH PREFACTORS INVOLVE WHICH PRIMES?")
print("=" * 80)
print("\nPrime content summary of ALL prefactors:")
all_primes = set()
prime_map = {}
for name, frac in prefactors.items():
    p, q = frac.numerator, frac.denominator
    primes = set(factorize(p).keys()) | set(factorize(q).keys())
    prime_map[name] = primes
    all_primes |= primes

print(f"All primes appearing: {sorted(all_primes)}")
print()
for name, primes in prime_map.items():
    beyond = primes - {2, 3}
    status = "PURE {2,3}" if not beyond else f"contains {beyond}"
    print(f"  {name:<25} primes: {sorted(primes):<20}  {status}")

# Count occurrences of each prime
print("\nPrime frequency across all prefactors:")
for p in sorted(all_primes):
    count = sum(1 for primes in prime_map.values() if p in primes)
    print(f"  prime {p}: appears in {count}/{len(prefactors)} prefactors")
