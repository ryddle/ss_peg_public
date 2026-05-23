#!/usr/bin/env python3
"""
two_adic_structure.py
=====================
Deep 2-adic analysis of the N matrix (n = 2·Δx₁₂).

2-adic integers Z₂: every integer has a unique expansion Σ aₖ 2ᵏ (aₖ∈{0,1}).
The 2-adic valuation v₂(n) = max power of 2 dividing n.
The 2-adic norm |n|₂ = 2^{−v₂(n)}.

We analyze:
  1. v₂ matrix and its structure
  2. Smith normal form of N over Z (2-adic invariant factors)
  3. N modulo 2ᵏ for k = 1,2,3
  4. 2-adic expansion of each entry
  5. The sign-valuation decomposition: n = (-1)^ε · 2^v · u, u odd
  6. Whether the v₂ matrix encodes quantum numbers
"""
from fractions import Fraction as F
import numpy as np
from math import gcd
from itertools import product as iprod

# ── data ──────────────────────────────────────────────────────────────
sectors = ['lepton', 'up', 'down']
coords = ['p', 'q', 'r']

# n matrix (integers)
n_vals = {
    ('lepton','p'): -1, ('lepton','q'): -8, ('lepton','r'):  6,
    ('up','p'):      4, ('up','q'):     -2, ('up','r'):      8,
    ('down','p'):    4, ('down','q'):    2, ('down','r'):     4,
}
flat = {('lepton','r'), ('up','q'), ('down','p')}

N = np.array([[n_vals[(s,k)] for k in coords] for s in sectors])

# ── 2-adic valuation ─────────────────────────────────────────────────
def v2(n):
    """2-adic valuation of integer n."""
    if n == 0: return float('inf')
    n = abs(n)
    v = 0
    while n % 2 == 0:
        v += 1
        n //= 2
    return v

def odd_part(n):
    """Odd part: n / 2^v₂(n)."""
    if n == 0: return 0
    sign = 1 if n > 0 else -1
    n = abs(n)
    while n % 2 == 0:
        n //= 2
    return sign * n

print("=" * 70)
print("  2-ADIC ANALYSIS OF THE N MATRIX")
print("=" * 70)

# ── 1. Basic decomposition: n = (-1)^ε · 2^v · u ────────────────────
print("\n§1. SIGN-VALUATION-UNIT DECOMPOSITION: n = (-1)^ε · 2^v · u")
print("-" * 70)

V2 = np.zeros((3,3), dtype=int)
eps = np.zeros((3,3), dtype=int)  # sign bit
unit = np.zeros((3,3), dtype=int)  # odd part (always ±1 except one)

print(f"\n  {'key':>12s}  {'n':>4s}  {'ε':>2s}  {'v₂':>3s}  {'u':>3s}  decomposition")
for i, s in enumerate(sectors):
    for j, k in enumerate(coords):
        n = n_vals[(s,k)]
        e = 0 if n > 0 else 1
        v = v2(n)
        u = odd_part(n)
        V2[i,j] = v
        eps[i,j] = e
        unit[i,j] = u
        tag = " [F]" if (s,k) in flat else ""
        print(f"  {s:>6s}/{k}:  {n:>+4d}   {e}   {v:>2d}  {u:>+3d}  "
              f"(-1)^{e}·2^{v}·{abs(u)}{tag}")

# ── 2. The v₂ matrix ─────────────────────────────────────────────────
print(f"\n§2. 2-ADIC VALUATION MATRIX v₂(N)")
print("-" * 70)
print(f"\n          p    q    r")
for i, s in enumerate(sectors):
    row = "".join(f"   {V2[i,j]:>2d}{'*' if (s,coords[j]) in flat else ' '}" for j in range(3))
    print(f"  {s:>6s}  {row}")

print(f"\n  Row sums:    {[sum(V2[i,:]) for i in range(3)]}")
print(f"  Col sums:    {[sum(V2[:,j]) for j in range(3)]}")
print(f"  Total:       {V2.sum()}")
print(f"  Trace:       {np.trace(V2)}")
print(f"  Anti-diag:   {V2[0,2] + V2[1,1] + V2[2,0]}")

# ── 3. Quantum number connections ────────────────────────────────────
print(f"\n§3. v₂ AND QUANTUM NUMBERS")
print("-" * 70)

# Quantum numbers
I3 = {'lepton': 0, 'up': F(1,2), 'down': F(-1,2)}
Q = {'lepton': {'p': -1, 'q': 0, 'r': 0},
     'up':     {'p': F(2,3), 'q': F(2,3), 'r': F(2,3)},
     'down':   {'p': F(-1,3), 'q': F(-1,3), 'r': F(-1,3)}}
Nc = {'lepton': 1, 'up': 3, 'down': 3}

print("\n  Testing: v₂ = f(I₃, Nc, coord)?")
print(f"  {'key':>12s}  {'v₂':>3s}  {'2I₃':>4s}  {'Nc':>3s}  {'Q':>5s}")
for i, s in enumerate(sectors):
    for j, k in enumerate(coords):
        tag = "*" if (s,k) in flat else " "
        print(f"  {s:>6s}/{k}{tag} {V2[i,j]:>3d}  "
              f"{float(2*I3[s]):>+4.0f}  {Nc[s]:>3d}  {float(Q[s][k]):>+5.2f}")

# ── 4. Sign matrix ───────────────────────────────────────────────────
print(f"\n§4. SIGN MATRIX ε (n < 0 → ε = 1)")
print("-" * 70)
print(f"\n          p    q    r")
for i, s in enumerate(sectors):
    row = "".join(f"   {eps[i,j]:>2d}{'*' if (s,coords[j]) in flat else ' '}" 
                  for j in range(3))
    print(f"  {s:>6s}  {row}")
print(f"\n  Negative entries: ", end="")
negs = [(s,k) for s in sectors for k in coords if n_vals[(s,k)] < 0]
print(", ".join(f"{s}/{k}" for s,k in negs))
print(f"  = {{l/p, l/q, u/q}} = lepton row (minus flat) + up/q (flat)")
print(f"  # negative: {sum(eps.flat)} out of 9")

# ── 5. Unit (odd part) matrix ────────────────────────────────────────
print(f"\n§5. ODD PART (UNIT) MATRIX u = n/2^v₂(n)")
print("-" * 70)
print(f"\n          p    q    r")
for i, s in enumerate(sectors):
    row = "".join(f"  {unit[i,j]:>+3d}{'*' if (s,coords[j]) in flat else ' '}" 
                  for j in range(3))
    print(f"  {s:>6s}  {row}")
print(f"\n  ALL units are ±1 EXCEPT lepton/r = +3 (flat)")
print(f"  det(unit matrix) = {int(np.linalg.det(unit.astype(float)))}")

# ── 6. Smith normal form over Z ──────────────────────────────────────
print(f"\n§6. SMITH NORMAL FORM OF N OVER Z")
print("-" * 70)

def smith_normal_form(M):
    """Compute Smith normal form of integer matrix M.
    Returns (D, L, R) such that L M R = D diagonal."""
    M = M.copy().tolist()
    n = len(M)
    L = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    R = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    
    for col in range(n):
        # Find pivot (smallest nonzero in submatrix)
        for iteration in range(100):
            # Find min abs nonzero in remaining submatrix
            min_val = None
            min_pos = None
            for i in range(col, n):
                for j in range(col, n):
                    if M[i][j] != 0:
                        if min_val is None or abs(M[i][j]) < min_val:
                            min_val = abs(M[i][j])
                            min_pos = (i, j)
            
            if min_pos is None:
                break
            
            pi, pj = min_pos
            # Swap to pivot position
            if pi != col:
                M[pi], M[col] = M[col], M[pi]
                L[pi], L[col] = L[col], L[pi]
            if pj != col:
                for i in range(n):
                    M[i][pj], M[i][col] = M[i][col], M[i][pj]
                    R[i][pj], R[i][col] = R[i][col], R[i][pj]
            
            # Eliminate column
            changed = False
            for i in range(col+1, n):
                if M[i][col] != 0:
                    q = M[i][col] // M[col][col]
                    for j in range(n):
                        M[i][j] -= q * M[col][j]
                        L[i][j] -= q * L[col][j]
                    if M[i][col] != 0:
                        changed = True
            
            # Eliminate row
            for j in range(col+1, n):
                if M[col][j] != 0:
                    q = M[col][j] // M[col][col]
                    for i in range(n):
                        M[i][j] -= q * M[i][col]
                        R[i][j] -= q * R[i][col]
                    if M[col][j] != 0:
                        changed = True
            
            if not changed:
                # Check if pivot divides all remaining
                divides_all = True
                for i in range(col+1, n):
                    for j in range(col+1, n):
                        if M[i][j] % M[col][col] != 0:
                            divides_all = False
                            # Add row i to row col
                            for jj in range(n):
                                M[col][jj] += M[i][jj]
                                L[col][jj] += L[i][jj]
                            break
                    if not divides_all:
                        break
                if divides_all:
                    break
        
        # Make diagonal positive
        if col < n and M[col][col] < 0:
            for j in range(n):
                M[col][j] = -M[col][j]
                L[col][j] = -L[col][j]
    
    return np.array(M), np.array(L), np.array(R)

D, L, R = smith_normal_form(N.copy())
diag = [D[i,i] for i in range(3)]
print(f"\n  N = {N.tolist()}")
print(f"  Smith normal form D = diag({diag})")
print(f"  Invariant factors: {diag}")
print(f"  v₂ of invariant factors: {[v2(d) for d in diag]}")
print(f"  Product = det = {np.prod(diag)} (expected {int(np.linalg.det(N))})")

# Verify
LNR = L @ N @ R
print(f"  Verification L·N·R = {[LNR[i,i] for i in range(3)]} (diag)")

# ── 7. N mod 2^k ─────────────────────────────────────────────────────
print(f"\n§7. N MODULO POWERS OF 2")
print("-" * 70)

for k in range(1, 5):
    m = 2**k
    Nmod = N % m
    print(f"\n  N mod {m} (= 2^{k}):")
    print(f"          p    q    r")
    for i, s in enumerate(sectors):
        row = "".join(f"  {Nmod[i,j]:>4d}" for j in range(3))
        print(f"  {s:>6s}  {row}")
    rnk = np.linalg.matrix_rank(Nmod.astype(float) % m)
    det_mod = int(round(np.linalg.det(Nmod.astype(float)))) % m
    print(f"  det(N) mod {m} = {det_mod}")

# ── 8. 2-adic expansion ──────────────────────────────────────────────
print(f"\n§8. 2-ADIC EXPANSION OF EACH ENTRY")
print("-" * 70)
print(f"\n  n = Σ aₖ·2ᵏ  (for positive; negate for sign)")

def two_adic_expansion(n, bits=4):
    """Return 2-adic digits of |n|."""
    m = abs(n)
    digits = []
    for _ in range(bits):
        digits.append(m % 2)
        m //= 2
    return digits

print(f"\n  {'key':>12s}  {'n':>4s}  2⁰ 2¹ 2² 2³  binary")
for i, s in enumerate(sectors):
    for j, k in enumerate(coords):
        n = n_vals[(s,k)]
        d = two_adic_expansion(n, 4)
        sign = "+" if n > 0 else "-"
        tag = "*" if (s,k) in flat else " "
        binary = f"{sign}{''.join(str(b) for b in reversed(d))}"
        print(f"  {s:>6s}/{k}{tag} {n:>+4d}   "
              f"{d[0]}  {d[1]}  {d[2]}  {d[3]}  {binary}")

# ── 9. The v₂ matrix as a LATIN SQUARE? ──────────────────────────────
print(f"\n§9. STRUCTURE OF THE v₂ MATRIX")
print("-" * 70)

print(f"\n  v₂ matrix:")
print(f"          p    q    r")
for i, s in enumerate(sectors):
    row = "".join(f"   {V2[i,j]:>2d}" for j in range(3))
    print(f"  {s:>6s}  {row}")

# Check: is it a permutation matrix + constant?
print(f"\n  Is v₂ a permutation of {0,1,2,3} per row/col?")
for i, s in enumerate(sectors):
    vals = sorted(V2[i,:])
    print(f"    {s} row: {list(V2[i,:])} sorted={vals}")
for j, k in enumerate(coords):
    vals = sorted(V2[:,j])
    print(f"    {k} col: {list(V2[:,j])} sorted={vals}")

# The v₂ matrix is:
#   0  3  1
#   2  1  3
#   2  1  2
# Row multisets: {0,1,3}, {1,2,3}, {1,2,2}
# Col multisets: {0,2,2}, {1,1,3}, {1,2,3}

# Check: replace flat entries — what are the FREE v₂ values?
print(f"\n  FREE entries of v₂ (non-flat):")
print(f"          p    q    r")
for i, s in enumerate(sectors):
    row = ""
    for j, k in enumerate(coords):
        if (s,k) in flat:
            row += f"   {'--':>2s}"
        else:
            row += f"   {V2[i,j]:>2d}"
    print(f"  {s:>6s}  {row}")

free_v2 = []
for s in sectors:
    for k in coords:
        if (s,k) not in flat:
            free_v2.append(V2[sectors.index(s), coords.index(k)])
free_v2_sorted = sorted(free_v2)
print(f"\n  Free v₂ values: {free_v2} → sorted: {free_v2_sorted}")
print(f"  = {{0, 1, 2, 2, 3, 3}}")
print(f"  Sum: {sum(free_v2)}")
print(f"  Mean: {sum(free_v2)/6:.2f}")

# ── 10. The ANTI-DIAGONAL structure ──────────────────────────────────
print(f"\n§10. DIAGONAL AND ANTI-DIAGONAL")
print("-" * 70)

diag_v2 = [V2[i,i] for i in range(3)]
anti_v2 = [V2[i,2-i] for i in range(3)]
print(f"  Main diagonal of v₂: {diag_v2} = {{0, 1, 2}}")
print(f"  Anti-diagonal of v₂: {anti_v2} = {{1, 1, 2}}")

diag_n = [N[i,i] for i in range(3)]
anti_n = [N[i,2-i] for i in range(3)]
print(f"  Main diagonal of N:  {diag_n} (product = {np.prod(diag_n)})")
print(f"  Anti-diagonal of N:  {anti_n} (product = {np.prod(anti_n)})")

# ── 11. The v₂ matrix vs flat position ───────────────────────────────
print(f"\n§11. FLAT COORD POSITION AND v₂")
print("-" * 70)
print("\n  Flat positions: (l,r)=(0,2), (u,q)=(1,1), (d,p)=(2,0)")
print("  = the ANTI-DIAGONAL!")
print(f"  v₂ on anti-diagonal (= flat): {anti_v2}")
print(f"  v₂ on main diagonal (= free): {diag_v2}")
print()
print("  The flat coords sit on the anti-diagonal of the 3×3 matrix.")
print("  This is the permutation σ = (p→r, q→q, r→p) applied to sectors.")

# ── 12. N over F₂ (the field with 2 elements) ───────────────────────
print(f"\n§12. N OVER F₂ (mod 2)")
print("-" * 70)
N2 = N % 2
print(f"\n  N mod 2:")
print(f"          p    q    r")
for i, s in enumerate(sectors):
    row = "".join(f"  {N2[i,j]:>4d}" for j in range(3))
    print(f"  {s:>6s}  {row}")

print(f"\n  Only entry ≡ 1 mod 2: lepton/p (n = -1)")
print(f"  rank(N mod 2) = {np.linalg.matrix_rank(N2 % 2)}")
print(f"  N mod 2 has rank 1 over F₂ — nearly singular!")

# ── 13. 2-adic filtration ────────────────────────────────────────────
print(f"\n§13. 2-ADIC FILTRATION")
print("-" * 70)
print("\n  Entries by 2-adic valuation level:")
for v in range(4):
    entries = [(s,k,n_vals[(s,k)]) for s in sectors for k in coords if v2(n_vals[(s,k)]) == v]
    if entries:
        names = ", ".join(f"{s}/{k}={n}" for s,k,n in entries)
        print(f"  v₂ = {v}: {names}")

print("\n  Layer decomposition N = N₀ + 2·N₁ + 4·N₂ + 8·N₃:")
# N₀ = entries with v₂ = 0 (only lepton/p = -1)
# The rest contribute at higher 2-powers

layers = {}
for v_level in range(4):
    layer = np.zeros((3,3), dtype=int)
    for i, s in enumerate(sectors):
        for j, k in enumerate(coords):
            n = n_vals[(s,k)]
            d = two_adic_expansion(n, 4)
            sgn = 1 if n > 0 else -1
            # The contribution at level v_level
            if v2(n) == v_level:
                layer[i,j] = odd_part(n)
    layers[v_level] = layer
    print(f"\n  Layer v₂={v_level} (entries with exactly this valuation, shown as odd part):")
    print(f"          p    q    r")
    for i, s in enumerate(sectors):
        row = "".join(f"  {layer[i,j]:>+4d}" for j in range(3))
        print(f"  {s:>6s}  {row}")

# ── 14. 2-adic distance between rows ─────────────────────────────────
print(f"\n§14. 2-ADIC DISTANCES BETWEEN SECTORS")
print("-" * 70)

def padic_norm(n, p=2):
    if n == 0: return 0.0
    v = 0
    m = abs(n)
    while m % p == 0:
        v += 1
        m //= p
    return p**(-v)

print("\n  |n₁ − n₂|₂ = 2-adic norm of difference:")
for i in range(3):
    for j in range(i+1, 3):
        diffs = N[i] - N[j]
        norms = [padic_norm(int(d)) for d in diffs]
        max_norm = max(norms)
        print(f"  d₂({sectors[i]}, {sectors[j]}):")
        print(f"    Δn = {diffs.tolist()}")
        print(f"    |Δn|₂ per coord = {norms}")
        print(f"    max |Δn|₂ = {max_norm} = 2^{-v2(int(np.min(np.abs(diffs[diffs!=0]))))}")

# ── 15. Characteristic polynomial ────────────────────────────────────
print(f"\n§15. CHARACTERISTIC POLYNOMIAL OF N")
print("-" * 70)

# det(N - λI) = -λ³ + tr·λ² - (sum of 2×2 minors)·λ + det
tr = int(np.trace(N))
det_N = int(round(np.linalg.det(N)))

# 2×2 minors (cofactors of diagonal)
M11 = N[1,1]*N[2,2] - N[1,2]*N[2,1]  # minor (0,0)
M22 = N[0,0]*N[2,2] - N[0,2]*N[2,0]  # minor (1,1)
M33 = N[0,0]*N[1,1] - N[0,1]*N[1,0]  # minor (2,2)
sum_minors = int(M11 + M22 + M33)

print(f"\n  N = {N.tolist()}")
print(f"  tr(N) = {tr}")
print(f"  Σ(2×2 minors on diag) = {sum_minors}")
print(f"  det(N) = {det_N}")
print(f"\n  χ(λ) = λ³ − {tr}λ² + ({sum_minors})λ − ({det_N})")
print(f"       = λ³ − λ² + ({sum_minors})λ + {-det_N}")
print(f"\n  2-adic valuations of coefficients:")
print(f"    v₂(1) = 0   [λ³]")
print(f"    v₂({tr}) = {v2(tr)}   [λ²]")
print(f"    v₂({sum_minors}) = {v2(sum_minors)}   [λ¹]")
print(f"    v₂({det_N}) = {v2(det_N)}   [λ⁰]")

# Newton polygon (2-adic)
print(f"\n  Newton polygon for p=2:")
points = [(0, v2(det_N)), (1, v2(sum_minors)), (2, v2(tr)), (3, 0)]
print(f"    Points (degree, v₂): {points}")
# Lower convex hull
print(f"    Slopes: ", end="")
for i in range(len(points)-1):
    slope = (points[i+1][1] - points[i][1]) / (points[i+1][0] - points[i][0])
    print(f"({points[i][0]}→{points[i+1][0]}): {slope:.2f}  ", end="")
print()
print(f"\n  Newton polygon slopes give 2-adic valuations of eigenvalues.")

# ── 16. Key summary ──────────────────────────────────────────────────
print(f"\n{'='*70}")
print(f"  SUMMARY: 2-ADIC ANATOMY OF N")
print(f"{'='*70}")
print(f"""
  N = [-1  -8   6]     v₂(N) = [0  3  1]     sign = [-  -  +]
      [ 4  -2   8]             [2  1  3]            [+  -  +]
      [ 4   2   4]             [2  1  2]            [+  +  +]

  1. det(N) = -8 = -2³  (pure power of 2)
     tr(N)  = 1          (2-adic unit)
     
  2. Smith normal form: diag({', '.join(str(d) for d in diag)})
     → Z³/N·Z³ ≅ Z/{diag[0]} × Z/{diag[1]} × Z/{diag[2]}
     
  3. Over F₂: rank 1.  N is "almost zero mod 2".
     Only lepton/p = -1 is odd.
     
  4. Flat entries on ANTI-DIAGONAL: (l,r), (u,q), (d,p)
     Their v₂: {anti_v2} → come from -4a (curvatures)
     
  5. Free entries off anti-diagonal:
     v₂ multiset = {{0, 1, 2, 2, 3, 3}}, sum = 11
     ALL have odd part ±1 (pure 2-powers with sign)
     
  6. The ONLY entry with odd part ≠ ±1 is lepton/r = 6 = 2·3 [FLAT]
     The factor 3 comes from a_r(lepton) = -3/2.
     
  7. Negative entries: l/p, l/q, u/q (= 2 free + 1 flat)
     All involve either lepton-sector or q-coordinate.
""")
