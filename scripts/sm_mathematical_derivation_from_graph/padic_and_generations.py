#!/usr/bin/env python3
"""
padic_and_generations.py
========================
1. Does the Smith tower {1,2,4} encode generations?
2. Full 3-adic analysis: does 3 really cancel, or is it hiding?
3. Combined p-adic portrait for p=2 and p=3.
"""
from fractions import Fraction as F
import numpy as np
from math import gcd
from itertools import product as iprod

# ── data ──────────────────────────────────────────────────────────────
sectors = ['lepton', 'up', 'down']
coords  = ['p', 'q', 'r']
flat = {('lepton','r'), ('up','q'), ('down','p')}

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

def x_at(s, k, g):
    return a_data[(s,k)]*g**2 + b_data[(s,k)]*g + c_data[(s,k)]

# n = 2·Δx₁₂
n_vals = {}
for s in sectors:
    for k in coords:
        n_vals[(s,k)] = int(2*(x_at(s,k,2) - x_at(s,k,1)))

N = np.array([[n_vals[(s,k)] for k in coords] for s in sectors])

# ── p-adic utilities ──────────────────────────────────────────────────
def vp(n, p):
    """p-adic valuation."""
    if n == 0: return float('inf')
    n = abs(n)
    v = 0
    while n % p == 0:
        v += 1
        n //= p
    return v

def unit_p(n, p):
    """p-free part: n / p^vp(n)."""
    if n == 0: return 0
    sign = 1 if n > 0 else -1
    m = abs(n)
    while m % p == 0:
        m //= p
    return sign * m

def smith_normal_form_int(M):
    """Smith normal form over Z. Returns diagonal entries."""
    M = [list(row) for row in M]
    n = len(M)
    
    for col in range(n):
        for _ in range(200):
            # find smallest nonzero in submatrix
            best = None
            for i in range(col, n):
                for j in range(col, n):
                    if M[i][j] != 0 and (best is None or abs(M[i][j]) < abs(best[0])):
                        best = (M[i][j], i, j)
            if best is None:
                break
            _, pi, pj = best
            # swap to pivot
            if pi != col:
                M[pi], M[col] = M[col], M[pi]
            if pj != col:
                for i in range(n):
                    M[i][pj], M[i][col] = M[i][col], M[i][pj]
            
            changed = False
            for i in range(col+1, n):
                if M[i][col] != 0:
                    q = M[i][col] // M[col][col]
                    for j in range(n):
                        M[i][j] -= q * M[col][j]
                    if M[i][col] != 0:
                        changed = True
            for j in range(col+1, n):
                if M[col][j] != 0:
                    q = M[col][j] // M[col][col]
                    for i in range(n):
                        M[i][j] -= q * M[i][col]
                    if M[col][j] != 0:
                        changed = True
            
            if not changed:
                ok = True
                for i in range(col+1, n):
                    for j in range(col+1, n):
                        if M[i][j] % M[col][col] != 0:
                            for jj in range(n):
                                M[col][jj] += M[i][jj]
                            ok = False
                            break
                    if not ok:
                        break
                if ok:
                    break
        
        if col < n and M[col][col] < 0:
            for j in range(n):
                M[col][j] = -M[col][j]
    
    return [M[i][i] for i in range(n)]

# ══════════════════════════════════════════════════════════════════════
# PART 1: SMITH TOWER AND GENERATIONS
# ══════════════════════════════════════════════════════════════════════
print("=" * 70)
print("  PART 1: SMITH TOWER {1, 2, 4} AND GENERATIONS")
print("=" * 70)

# The Smith form diag(1,2,4) means:
# There exist unimodular L,R (det ±1) such that L·N·R = diag(1,2,4)
# Equivalently: N defines a map Z³ → Z³ whose cokernel is Z/1 × Z/2 × Z/4

print("\n§1a. Generation matrices: what does N do to each generation?")
print("-" * 70)

# x(g) is a 3×3 matrix (sector × coord) for each g
# Δx₁₂ = N/2
# Δx₂₃ = x(3) - x(2)

# Build the three generation vectors in the 9-dim space
# But more relevantly: N acts on the COORDINATE space {p,q,r}
# For each sector, N gives the jump Δ(2x)₁₂

# Actually let's look at ALL transition matrices:
# T(g→g+1) = 2·(x(g+1) - x(g)) per sector
print("\n  Transition matrix T₁₂ = 2·Δx₁₂ = N:")
print(f"  {N}")

T23 = np.array([[int(2*(x_at(s,k,3) - x_at(s,k,2))) for k in coords] for s in sectors])
print(f"\n  Transition matrix T₂₃ = 2·Δx₂₃:")
print(f"  {T23}")

T13 = N + T23  # 2·(x(3) - x(1))
print(f"\n  Full jump T₁₃ = 2·(x(3) - x(1)) = N + T₂₃:")
print(f"  {T13}")

print(f"\n  Key: T₂₃ - T₁₂ = 4·a (curvature):")
diff = T23 - N
a_check = np.array([[int(4*a_data[(s,k)]) for k in coords] for s in sectors])
print(f"  T₂₃ - N = {diff.tolist()}")
print(f"  4a       = {a_check.tolist()}")
print(f"  Match: {np.array_equal(diff, a_check)}")

# Smith forms of all three
print(f"\n§1b. Smith normal forms:")
print("-" * 70)
for label, mat in [("T₁₂ = N", N), ("T₂₃", T23), ("T₁₃", T13)]:
    d = smith_normal_form_int(mat)
    det_val = int(round(np.linalg.det(mat.astype(float))))
    print(f"  {label}: Smith = diag({d}), det = {det_val}")
    print(f"    v₂ of factors: {[vp(x,2) for x in d]}")
    print(f"    v₃ of factors: {[vp(x,3) for x in d]}")

# ── Does the Smith form change with generation? ──────────────────────
print(f"\n§1c. Curvature matrix 4A and its Smith form:")
print("-" * 70)
A4 = a_check
d_A = smith_normal_form_int(A4)
print(f"  4A = {A4.tolist()}")
print(f"  Smith(4A) = diag({d_A})")
print(f"  det(4A) = {int(round(np.linalg.det(A4.astype(float))))}")

# ── The key relation: T(g→g+1) = N + (g-1)·4A ──────────────────────
# T₁₂ = N, T₂₃ = N + 4A, T₃₄ = N + 8A (if g=4 existed)
print(f"\n§1d. T(g→g+1) = N + (g−1)·4A")
print("-" * 70)
print(f"  T₁₂ = N + 0·4A = N")
print(f"  T₂₃ = N + 1·4A")
print(f"  T₃₄ = N + 2·4A  (hypothetical g=4)")

T34 = N + 2*A4
d34 = smith_normal_form_int(T34)
det34 = int(round(np.linalg.det(T34.astype(float))))
print(f"\n  T₃₄ = {T34.tolist()}")
print(f"  Smith(T₃₄) = diag({d34}), det = {det34}")

print(f"\n  Summary of Smith forms by generation transition:")
print(f"  {'transition':>12s}  {'Smith':>20s}  {'det':>6s}  {'v₂(det)':>8s}  {'v₃(det)':>8s}")
for label, mat in [("g=1→2 (N)", N), ("g=2→3", T23), ("g=3→4*", T34)]:
    d = smith_normal_form_int(mat)
    det_val = int(round(np.linalg.det(mat.astype(float))))
    print(f"  {label:>12s}  diag({d})      {det_val:>+6d}  "
          f"  {vp(abs(det_val),2):>5d}     {vp(abs(det_val),3):>5d}")

# ── What generation does the lattice select? ──────────────────────────
print(f"\n§1e. Which integer g makes T(g→g+1) special?")
print("-" * 70)
print(f"  det(N + (g−1)·4A) as function of g:")
for g in range(0, 7):
    T = N + (g-1)*A4
    det_val = int(round(np.linalg.det(T.astype(float))))
    d = smith_normal_form_int(T)
    v2_det = vp(abs(det_val), 2) if det_val != 0 else "∞"
    v3_det = vp(abs(det_val), 3) if det_val != 0 else "∞"
    marker = " ← physical" if g in [1,2,3] else ""
    print(f"  g={g}: det = {det_val:>8d}  v₂={str(v2_det):>3s}  v₃={str(v3_det):>3s}  "
          f"Smith=diag({d}){marker}")

# ══════════════════════════════════════════════════════════════════════
# PART 2: 3-ADIC ANALYSIS
# ══════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 2: 3-ADIC ANALYSIS — WHERE IS THE 3?")
print(f"{'='*70}")

print(f"\n§2a. 3-adic valuation of N entries")
print("-" * 70)

V3 = np.zeros((3,3), dtype=int)
U3 = np.zeros((3,3), dtype=int)
print(f"\n  {'key':>12s}  {'n':>4s}  {'v₃':>3s}  {'3-free':>7s}")
for i, s in enumerate(sectors):
    for j, k in enumerate(coords):
        n = n_vals[(s,k)]
        V3[i,j] = vp(n, 3)
        U3[i,j] = unit_p(n, 3)
        tag = " [F]" if (s,k) in flat else ""
        print(f"  {s:>6s}/{k}:  {n:>+4d}   {V3[i,j]:>2d}  {U3[i,j]:>+6d}{tag}")

print(f"\n  v₃ matrix:")
print(f"          p    q    r")
for i, s in enumerate(sectors):
    row = "".join(f"   {V3[i,j]:>2d}{'*' if (s,coords[j]) in flat else ' '}" for j in range(3))
    print(f"  {s:>6s}  {row}")

print(f"\n  ALL v₃ = 0 EXCEPT lepton/r = 6 (v₃ = 1) ← the SAME anomalous entry!")
print(f"  det(N) = -8: v₃(det) = {vp(8,3)} → 3 is INVISIBLE in the determinant")

print(f"\n§2b. 3-free (unit) part matrix:")
print(f"          p    q    r")
for i, s in enumerate(sectors):
    row = "".join(f"  {U3[i,j]:>+4d}{'*' if (s,coords[j]) in flat else ' '}" for j in range(3))
    print(f"  {s:>6s}  {row}")

# ── 3 in the COORDINATES themselves (not just Δx₁₂) ─────────────────
print(f"\n§2c. Where does 3 ACTUALLY live? In 2x(g)!")
print("-" * 70)
print(f"  Remember: the lattice basis is {{ln2, lnR₃, ln3}}.")
print(f"  The r-coordinate IS the ln3 direction.")
print(f"  Let's trace v₃ through all matrices 2x(g):")

for g in [1, 2, 3]:
    print(f"\n  v₃(2x(g={g})):")
    print(f"          p    q    r")
    for i, s in enumerate(sectors):
        row = ""
        for j, k in enumerate(coords):
            val = int(2*x_at(s, k, g))
            v3_val = vp(val, 3) if val != 0 else "-"
            row += f"  {str(v3_val):>4s}"
        print(f"  {s:>6s}  {row}")

print(f"\n  And v₃ of the 2x values themselves:")
for g in [1, 2, 3]:
    names = {1: ['e','u','d'], 2: ['mu','c','s'], 3: ['tau','t','b']}
    print(f"\n  2x(g={g}) [{'/'.join(names[g])}]:")
    print(f"  {'':>8s}     p       q       r     | v₃(p)  v₃(q)  v₃(r)")
    for i, s in enumerate(sectors):
        vals = [int(2*x_at(s,k,g)) for k in coords]
        v3s = [vp(v,3) if v != 0 else "-" for v in vals]
        row_v = "".join(f"  {str(v):>5s} " for v in v3s)
        row_n = "".join(f"  {v:>+5d} " for v in vals)
        print(f"  {s:>8s}  {row_n}   |{row_v}")

# ── 3 in the curvature ───────────────────────────────────────────────
print(f"\n§2d. 3-adic content of curvatures (4a):")
print("-" * 70)
print(f"          p    q    r")
for i, s in enumerate(sectors):
    row = ""
    for j, k in enumerate(coords):
        v = int(4*a_data[(s,k)])
        row += f"  {v:>+4d}(v₃={vp(v,3) if v!=0 else '-'}) "
    print(f"  {s:>6s}  {row}")

# ── Where 3 enters: the actual masses ln(m/me) have ln3 component ────
print(f"\n§2e. The r-coordinate extracts the 3-ADIC content of mass ratios")
print("-" * 70)
print(f"  ln(m/mₑ) = p·ln2 + q·lnR₃ + r·ln3")
print(f"  The POWER OF 3 in the mass ratio is determined by r alone:")
print(f"  m/mₑ = 2^p · R₃^q · 3^r   →   v₃(m/mₑ) = r  (when r is integer)")
print()

for g in [1, 2, 3]:
    names = {1: ['e','u','d'], 2: ['μ','c','s'], 3: ['τ','t','b']}
    for i, s in enumerate(sectors):
        r_val = x_at(s, 'r', g)
        name = names[g][i]
        is_int = (r_val.denominator == 1)
        print(f"  {name:>3s}: r = {float(r_val):>+5.1f}  "
              f"{'→ 3^' + str(int(r_val)) if is_int else '(half-integer)'}")

# ── 3-adic N vs 2-adic N ─────────────────────────────────────────────
print(f"\n§2f. N over F₃ (mod 3):")
print("-" * 70)
N3 = N % 3
# Handle negatives properly
N3_proper = np.array([[n_vals[(s,k)] % 3 for k in coords] for s in sectors])
print(f"          p    q    r")
for i, s in enumerate(sectors):
    row = "".join(f"  {N3_proper[i,j]:>4d}" for j in range(3))
    print(f"  {s:>6s}  {row}")
det3 = int(round(np.linalg.det(N3_proper.astype(float)))) % 3
rank3 = np.linalg.matrix_rank(N3_proper.astype(float))
print(f"\n  rank(N mod 3) = {rank3}")
print(f"  det(N) mod 3 = {(-8) % 3} = 1  (invertible over F₃!)")
print(f"\n  Compare: rank(N mod 2) = 1  (nearly singular)")
print(f"           rank(N mod 3) = 3  (FULL RANK, invertible)")

# ── Newton polygon for p=3 ───────────────────────────────────────────
print(f"\n§2g. Newton polygon of χ(λ) for p=3:")
print("-" * 70)
# χ(λ) = λ³ − λ² − 18λ + 8
coeffs_chi = [8, -18, -1, 1]  # from λ⁰ to λ³
print(f"  χ(λ) = λ³ − λ² − 18λ + 8")
print(f"  Coefficients [λ⁰, λ¹, λ², λ³] = {coeffs_chi}")
print(f"  v₃ of coefficients: {[vp(c,3) if c != 0 else '∞' for c in coeffs_chi]}")
points = [(i, vp(coeffs_chi[i], 3)) for i in range(4)]
print(f"  Newton polygon points: {points}")
print(f"  Slopes: ", end="")
for i in range(3):
    s = (points[i+1][1] - points[i][1])
    print(f"({i}→{i+1}): {s}  ", end="")
print()
print(f"\n  All slopes ≤ 0: all eigenvalues have v₃ ≥ 0 (are 3-adic integers)")
print(f"  Slope (0→1) = {points[1][1] - points[0][1]}: one eigenvalue with v₃ = {abs(points[1][1] - points[0][1])}")

# ══════════════════════════════════════════════════════════════════════
# PART 3: COMBINED PORTRAIT — THE FULL PICTURE
# ══════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 3: COMBINED p-ADIC PORTRAIT")
print(f"{'='*70}")

print(f"\n§3a. Complete factorization: n = (-1)^ε · 2^a · 3^b · u")
print("-" * 70)
print(f"  {'key':>12s}  {'n':>4s}  {'ε':>2s}  {'2^a':>4s}  {'3^b':>4s}  {'u':>3s}  "
      f"decomposition")
for i, s in enumerate(sectors):
    for j, k in enumerate(coords):
        n = n_vals[(s,k)]
        eps = 0 if n > 0 else 1
        a = vp(n, 2)
        b = vp(n, 3)
        u = abs(n) // (2**a * 3**b)
        tag = " [F]" if (s,k) in flat else ""
        print(f"  {s:>6s}/{k}:  {n:>+4d}   {eps}  2^{a}   3^{b}   {u:>2d}  "
              f"(-1)^{eps}·2^{a}·3^{b}{tag}")

print(f"\n  ALL residual units u = 1!")
print(f"  Every n is of the form ±2^a · 3^b  (a ∈ {{0..3}}, b ∈ {{0,1}})")
print(f"  Only ONE entry has b=1: lepton/r = +2¹·3¹ = 6")
print(f"  The other 8 are PURE 2-powers with sign: n = ±2^a")

# ── The (a,b) grid ───────────────────────────────────────────────────
print(f"\n§3b. The (v₂, v₃) grid:")
print("-" * 70)
print(f"\n  Each n lives at a point (v₂, v₃) in the 'prime spectrum':")
print(f"          p         q         r")
for i, s in enumerate(sectors):
    row = ""
    for j, k in enumerate(coords):
        n = n_vals[(s,k)]
        a = vp(n, 2)
        b = vp(n, 3)
        tag = "*" if (s,k) in flat else " "
        row += f"  ({a},{b}){tag}  "
    print(f"  {s:>6s}  {row}")

print(f"\n  The ONLY point off the v₃=0 axis: lepton/r at (1,1)")
print(f"  THIS is where ln3 'leaks' into the transition matrix.")

# ── The (v₂, sign) as quantum numbers ────────────────────────────────
print(f"\n§3c. (v₂, sign) as effective quantum numbers for FREE entries:")
print("-" * 70)
print(f"  Free entries only (6 values):")
print(f"  {'key':>12s}  {'n':>4s}  {'v₂':>3s}  {'sgn':>4s}  →  bit pattern")
free_entries = []
for s in sectors:
    for k in coords:
        if (s,k) not in flat:
            n = n_vals[(s,k)]
            a = vp(n, 2)
            sgn = "+" if n > 0 else "−"
            bits = f"{'+' if n>0 else '-'}{a}"
            free_entries.append((s, k, n, a, sgn, bits))
            print(f"  {s:>6s}/{k}:  {n:>+4d}   {a:>2d}   {sgn:>3s}   →  {bits}")

print(f"\n  6 free entries encode: 1 sign bit + v₂ ∈ {{0,1,2,3}} (2 bits)")
print(f"  Total info: 6 × 3 = 18 bits")
print(f"  This is EXACTLY the 18 DOF count!")

# ── Connection to generations ─────────────────────────────────────────
print(f"\n§3d. Smith tower and the generation lattice:")
print("-" * 70)
print(f"  Smith(N) = diag(1, 2, 4) = (2⁰, 2¹, 2²)")
print(f"  Cokernel: Z³/N·Z³ ≅ Z/1 × Z/2 × Z/4 = {{0}} × Z₂ × Z₄")
print(f"  This is a group of order 8 = 2³")
print(f"  |cokernel| = |det(N)| = 8")
print()

# What are the cosets? Which generation vectors land where?
# The physical generations g=1,2,3 give position vectors x(g).
# In the 2x representation (integer), 2x(g) ∈ Z³ per sector.
# The cokernel lives in Z³/(image of N).

# For the lepton sector:
print(f"  Lepton sector (as integer vectors 2x(g)):")
for g in [1, 2, 3]:
    xvec = [int(2*x_at('lepton', k, g)) for k in coords]
    print(f"    g={g}: {xvec}")
print(f"    Δ₁₂ = N[lepton] = {list(N[0])}")
print(f"    Δ₂₃ = T₂₃[lepton] = {list(T23[0])}")

# The image of N in Z³ is generated by the 3 rows of N
# Cokernel = Z³ / (row-span of N)
# Since Smith = diag(1,2,4), we can change basis so the
# quotient map is simply (x,y,z) → (0, y mod 2, z mod 4)

print(f"\n  In Smith basis, the quotient map is:")
print(f"    φ: Z³ → Z/1 × Z/2 × Z/4")
print(f"    φ(x,y,z) = (0, y mod 2, z mod 4)")
print()

# Find transformation matrices
# We need L, R such that L·N·R = diag(1,2,4)
# Then in new basis e' = R⁻¹·e, the map is diagonal

# Let's just check: what images do the generation vectors have?
# For that we need the actual Smith transformation
# Simpler: reduce N to Smith form and track transformations

def smith_with_transform(M):
    """Return (D, L, R) with L·M·R = D."""
    n = len(M)
    M = [list(row) for row in M]
    L = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    R = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    
    for col in range(n):
        for _ in range(200):
            best = None
            for i in range(col, n):
                for j in range(col, n):
                    if M[i][j] != 0 and (best is None or abs(M[i][j]) < abs(best[0])):
                        best = (M[i][j], i, j)
            if best is None:
                break
            _, pi, pj = best
            if pi != col:
                M[pi], M[col] = M[col], M[pi]
                L[pi], L[col] = L[col], L[pi]
            if pj != col:
                for i in range(n):
                    M[i][pj], M[i][col] = M[i][col], M[i][pj]
                    R[i][pj], R[i][col] = R[i][col], R[i][pj]
            
            changed = False
            for i in range(col+1, n):
                if M[i][col] != 0:
                    q = M[i][col] // M[col][col]
                    for j in range(n):
                        M[i][j] -= q * M[col][j]
                        L[i][j] -= q * L[col][j]
                    if M[i][col] != 0:
                        changed = True
            for j in range(col+1, n):
                if M[col][j] != 0:
                    q = M[col][j] // M[col][col]
                    for i in range(n):
                        M[i][j] -= q * M[i][col]
                        R[i][j] -= q * R[i][col]
                    if M[col][j] != 0:
                        changed = True
            if not changed:
                ok = True
                for i in range(col+1, n):
                    for j in range(col+1, n):
                        if M[i][j] % M[col][col] != 0:
                            for jj in range(n):
                                M[col][jj] += M[i][jj]
                                L[col][jj] += L[i][jj]
                            ok = False
                            break
                    if not ok:
                        break
                if ok:
                    break
        if col < n and M[col][col] < 0:
            for j in range(n):
                M[col][j] = -M[col][j]
                L[col][j] = -L[col][j]
    
    return (np.array(M), np.array(L), np.array(R))

D, L, R = smith_with_transform(N.tolist())
print(f"  L·N·R = D, where D = diag{tuple(D[i,i] for i in range(3))}")
print(f"  L = {L.tolist()}")
print(f"  R = {R.tolist()}")
print(f"  det(L) = {int(round(np.linalg.det(L.astype(float))))}, "
      f"det(R) = {int(round(np.linalg.det(R.astype(float))))}")

# Verify
LNR = L @ N @ R
print(f"  Verify: L·N·R = {LNR.tolist()}")

# Now: in the new basis, the generation vectors map via L
# The cokernel class of a vector v is: (L·v) mod diag(1,2,4)
print(f"\n  Cokernel class φ(v) = (L·v) mod (1,2,4):")
print(f"  For each sector and generation:")
for s_idx, s in enumerate(sectors):
    print(f"\n  {s}:")
    for g in [1, 2, 3]:
        xvec = np.array([int(2*x_at(s, k, g)) for k in coords])
        Lx = L @ xvec
        cokernel = [int(Lx[i]) % D[i,i] for i in range(3)]
        print(f"    g={g}: 2x = {str(xvec.tolist()):>20s}  →  L·(2x) = {str(Lx.tolist()):>20s}  "
              f"→  mod(1,2,4) = {cokernel}")

# ── The BIG test: det(N + g·4A) ──────────────────────────────────────
print(f"\n§3e. 2-adic and 3-adic content of det(T) vs generation:")
print("-" * 70)
print(f"  {'g':>3s}  {'det':>8s}  {'±2^a·3^b':>12s}  {'v₂':>4s}  {'v₃':>4s}  {'residual':>8s}")
for g in range(-1, 8):
    T = N + (g-1)*A4
    det_val = int(round(np.linalg.det(T.astype(float))))
    if det_val == 0:
        print(f"  {g:>3d}  {det_val:>8d}  {'0':>12s}  {'∞':>4s}  {'∞':>4s}  {'-':>8s}")
        continue
    a2 = vp(det_val, 2)
    a3 = vp(det_val, 3)
    res = abs(det_val) // (2**a2 * 3**a3)
    sgn = "+" if det_val > 0 else "-"
    decomp = f"{sgn}2^{a2}·3^{a3}"
    if res > 1: decomp += f"·{res}"
    marker = ""
    if g in [1,2,3]: marker = " ← gen " + str(g)
    print(f"  {g:>3d}  {det_val:>8d}  {decomp:>12s}  {a2:>4d}  {a3:>4d}  {res:>8d}{marker}")

# ══════════════════════════════════════════════════════════════════════
# PART 4: THE CRUCIAL INSIGHT
# ══════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 4: WHERE 3 HIDES AND WHY 2 DOMINATES")
print(f"{'='*70}")

print(f"""
  The lattice basis is {{ln2, lnR₃, ln3}}.
  The coordinates (p, q, r) are the exponents in this basis.
  
  KEY OBSERVATION: the transition matrix N = 2·Δx₁₂ lives in Z
  (integers), NOT in the lattice Z/2.
  
  When we compute n = 2·(x(2) - x(1)):
    - The factor 2 clears ALL denominators (coords are half-integers)
    - The p-direction (ln2) contributes ONLY powers of 2 to n
    - The q-direction (lnR₃) contributes ONLY powers of 2 to n  
    - The r-direction (ln3) contributes factors of 2 AND 3
  
  The r-coordinate is the GATEWAY for 3 into the integer matrix.
  But it enters only through the curvature a_r:
    a_r(lepton) = -3/2  →  4a = -6  (factor of 3!)
    a_r(up)     = -3/2  →  4a = -6  (same!)
    a_r(down)   =  1    →  4a =  4  (no factor of 3)
  
  For N itself (g=1→2 transition):
    n_r = 6 + 2I  →  n_r(0) = 6 = 2·3   (lepton: 3 enters)
                      n_r(±1) = 8, 4       (quarks: 3 cancels!)
  
  The 3 ENTERS at g=1→2 for the lepton, then gets ABSORBED
  because subsequent transitions add 4a_r which ALSO contains 3.
  
  In the determinant: det(N) = -8 = -2³.
  The single factor of 3 in entry (lepton,r) gets cancelled by
  the matrix algebra (cofactor expansion).
  
  CONCLUSION:
  ━━━━━━━━━━
  • 2 is the STRUCTURAL prime: Smith form = (2⁰, 2¹, 2²),
    determinant = -2³, 8 of 9 entries are ±2^a.
    
  • 3 is the GEOMETRIC prime: it enters only through one lattice 
    direction (ln3), affects only one coordinate (r), and only 
    one sector at the first transition (lepton/r = 6).
    
  • The Smith tower (1, 2, 4) IS the generation hierarchy:
    the cokernel Z/2 × Z/4 has elements of order 1, 2, 4 = 2^g-1
    for g = 1, 2, 3.
    
  • The mass spectrum is fundamentally 2-ADIC in character,
    with 3 appearing only as a "geometric echo" of the lattice.
""")
