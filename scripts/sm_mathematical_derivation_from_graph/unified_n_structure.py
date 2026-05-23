#!/usr/bin/env python3
"""
unified_n_structure.py
======================
Look for a SINGLE rule governing all 9 values of n = 2Δx₁₂.
The 6 free ones are powers of 2, but lepton/r gives 6 = 2×3.
Is there a unified structure?
"""
from fractions import Fraction as F
import numpy as np

# ── data ──────────────────────────────────────────────────────────────
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
flat = {('lepton','r'), ('up','q'), ('down','p')}
sectors = ['lepton', 'up', 'down']
coords = ['p', 'q', 'r']

# ── compute key matrices ──────────────────────────────────────────────
# n = 2·Δx₁₂ = 2(x(2) - x(1)) = 2(3a + b) = 2(b + 3a)
# Also: n = 2(δb - a)  where δb = b + 4a

def x_at(s, k, g):
    a = a_data[(s,k)]
    b = b_data[(s,k)]
    c = c_data[(s,k)]
    return a*g**2 + b*g + c

print("="*70)
print("  KEY INSIGHT: n = 2·Δx₁₂ = 2(x(2) − x(1))")
print("="*70)
print()

# Build matrices as Fraction arrays
N = {}  # n matrix
Dx12 = {}  # Δx₁₂
Dx23 = {}  # Δx₂₃
A = {}  # curvature
for s in sectors:
    for k in coords:
        dx12 = x_at(s, k, 2) - x_at(s, k, 1)
        dx23 = x_at(s, k, 3) - x_at(s, k, 2)
        N[(s,k)] = 2 * dx12
        Dx12[(s,k)] = dx12
        Dx23[(s,k)] = dx23
        A[(s,k)] = a_data[(s,k)]

print("  Δx₁₂ = x(2) − x(1) matrix:")
print(f"  {'':>8s}     p       q       r")
for s in sectors:
    row = ""
    for k in coords:
        v = Dx12[(s,k)]
        tag = "*" if (s,k) in flat else " "
        row += f"  {float(v):>+5.1f}{tag} "
    print(f"  {s:>8s}  {row}")

print()
print("  Δx₂₃ = x(3) − x(2) matrix:")
print(f"  {'':>8s}     p       q       r")
for s in sectors:
    row = ""
    for k in coords:
        v = Dx23[(s,k)]
        tag = "*" if (s,k) in flat else " "
        row += f"  {float(v):>+5.1f}{tag} "
    print(f"  {s:>8s}  {row}")

print()
print("  n = 2Δx₁₂ matrix:")
print(f"  {'':>8s}     p       q       r")
for s in sectors:
    row = ""
    for k in coords:
        v = N[(s,k)]
        tag = "*" if (s,k) in flat else " "
        row += f"  {float(v):>+5.0f}{tag} "
    print(f"  {s:>8s}  {row}")

# ── matrix algebra on n ───────────────────────────────────────────────
print()
print("="*70)
print("  MATRIX ALGEBRA ON N")
print("="*70)

n_mat = np.array([[float(N[(s,k)]) for k in coords] for s in sectors])
a_mat = np.array([[float(A[(s,k)]) for k in coords] for s in sectors])

print(f"\n  N matrix:\n{n_mat}")
print(f"\n  det(N) = {np.linalg.det(n_mat):.1f}")
print(f"  tr(N)  = {np.trace(n_mat):.1f}")

eigvals = np.linalg.eigvals(n_mat)
print(f"  eigenvalues: {eigvals}")
print(f"  |eigenvalues|: {np.abs(eigvals)}")
print(f"  product of eigenvalues = det = {np.prod(eigvals):.1f}")

# Singular values
U, s_vals, Vt = np.linalg.svd(n_mat)
print(f"  singular values: {s_vals}")
print(f"  product of sing vals: {np.prod(s_vals):.1f}")

# ── factor n = something × something? ────────────────────────────────
print()
print("="*70)
print("  FACTORIZATION ATTEMPTS")
print("="*70)

# n = 2(3a + b). Can we write the 3x3 matrix as a rank-1 or structured product?
# Let's try: is N = u ⊗ v + w ⊗ z (rank 2)?
rank = np.linalg.matrix_rank(n_mat)
print(f"\n  rank(N) = {rank}")

# Is N = D_sector × M × D_coord for diagonal D?
# Try: is there a row-scaling that makes N simpler?
print("\n  Row-normalized N (divide each row by GCD):")
for i, s in enumerate(sectors):
    row = [int(N[(s,k)]) for k in coords]
    from math import gcd
    g12 = gcd(abs(row[0]), abs(row[1]))
    g = gcd(g12, abs(row[2]))
    norm_row = [x // g for x in row]
    print(f"    {s}: {row} / {g} = {norm_row}")

# ── key test: n = -4a for flat, but what is n for free? ──────────────
print()
print("="*70)
print("  WHAT DETERMINES n FOR EACH COORD?")
print("="*70)
print()
print("  For flat:  n = -4a  (forced by flatness b=-5a)")
print("  For free:  n = -4a + 2δb  (δb is the free parameter)")
print()
print("  Let's look at n + 4a = 2δb (the 'excess over flat'):")
print(f"  {'':>8s}     p       q       r")
for s in sectors:
    row = ""
    for k in coords:
        v = N[(s,k)] + 4*A[(s,k)]
        tag = "*" if (s,k) in flat else " "
        row += f"  {float(v):>+5.0f}{tag} "
    print(f"  {s:>8s}  {row}")
print("  (flat coords should be 0 = -2a·2 + 2(-a) + 4a... wait)")
print()

# Actually 2δb = n + 4a? No.
# n = 2(δb - a), so 2δb = n + 2a.
print("  CORRECTION: 2δb = n + 2a:")
print(f"  {'':>8s}     p       q       r")
db2 = {}
for s in sectors:
    row = ""
    for k in coords:
        v = N[(s,k)] + 2*A[(s,k)]
        db2[(s,k)] = v
        tag = "*" if (s,k) in flat else " "
        row += f"  {float(v):>+5.0f}{tag} "
    print(f"  {s:>8s}  {row}")

# ── THE BIG QUESTION: is there a single formula for Δx₁₂? ───────────
print()
print("="*70)
print("  DIRECT: what are the actual x values?")
print("="*70)
print()

# Print x(1), x(2), x(3) matrices
for g in [1, 2, 3]:
    names_g = {1: ['e', 'u', 'd'], 2: ['μ', 'c', 's'], 3: ['τ', 't', 'b']}
    label = names_g[g]
    print(f"  x(g={g}):   [{label[0]:>2s}, {label[1]:>2s}, {label[2]:>2s}]")
    print(f"  {'':>8s}     p       q       r")
    for i, s in enumerate(sectors):
        row = ""
        for k in coords:
            v = x_at(s, k, g)
            row += f"  {float(v):>+5.1f}  "
        print(f"  {s:>8s}  {row}")
    print()

# ── Look at 2x values (always integers) ──────────────────────────────
print("="*70)
print("  INTEGER MATRICES: 2x(g)")
print("="*70)
print()
for g in [1, 2, 3]:
    names_g = {1: ['e', 'u', 'd'], 2: ['μ', 'c', 's'], 3: ['τ', 't', 'b']}
    label = names_g[g]
    print(f"  2x(g={g}):")
    print(f"  {'':>8s}     p       q       r     Σrow")
    for i, s in enumerate(sectors):
        vals = []
        for k in coords:
            v = 2*x_at(s, k, g)
            vals.append(int(v))
        row = "".join(f"  {v:>+5d}  " for v in vals)
        print(f"  {s:>8s}  {row}  {sum(vals):>+5d}")
    # Column sums
    col_sums = []
    for k in coords:
        cs = sum(int(2*x_at(s, k, g)) for s in sectors)
        col_sums.append(cs)
    row = "".join(f"  {v:>+5d}  " for v in col_sums)
    print(f"  {'Σcol':>8s}  {row}")
    print()

# ── Factorize n through the GENERATING FUNCTION structure ────────────
print("="*70)
print("  N THROUGH THE GENERATING FUNCTION")
print("="*70)
print()
print("  The generating function gives x_k(g) as function of (2I3, Nc, g).")
print("  For the 3 sectors:")
print("    lepton: (2I3, Nc) = (0, 1)")
print("    up:     (2I3, Nc) = (1, 3)")
print("    down:   (2I3, Nc) = (-1, 3)")
print()
print("  So Δx₁₂ = F(2I3, Nc, g=2) − F(2I3, Nc, g=1)")
print()

# The x values are: x_k(g) = a_k(I3,Nc)·g² + b_k(I3,Nc)·g + c_k(I3,Nc)
# Recall from parameter_reduction:
#   a_p(2I3, Nc) = (1/4)(2I3)² + (1/4)(2I3) + (3/4)Nc − (3/4)
#   etc. (we have the exact formulas)
# So Δx₁₂ = 3a + b is also a polynomial in (2I3, Nc).

# Let I = 2I3, N = Nc
Is = [0, 1, -1]  # lepton, up, down
Ncs = [1, 3, 3]

print("  Δx₁₂ as a function of (I=2I3, N=Nc):")
print(f"  {'':>8s}  I   N     Δp     Δq     Δr")
for idx, s in enumerate(sectors):
    I = Is[idx]
    Nc = Ncs[idx]
    dps = [float(Dx12[(s,k)]) for k in coords]
    print(f"  {s:>8s}  {I:>+2d}  {Nc:>2d}  {dps[0]:>+6.2f} {dps[1]:>+6.2f} {dps[2]:>+6.2f}")

print()

# Try to fit: Δx_k = α_k + β_k·I + γ_k·N + δ_k·I² + ε_k·I·N + ...
# We have 3 data points per coordinate, so max 3 parameters
# Use (I, N) = (0,1), (1,3), (-1,3)
# Basis: {1, I, N} or {1, I, I²} since N = 1 + 2|I| but not polynomial...

# Actually N and I are NOT independent in this 3-point set.
# (I,N) = (0,1), (1,3), (-1,3). So N = 1 + 2I² (check: 1+0=1, 1+2=3, 1+2=3) ✓!
# This means N − 1 = 2I², so any function of (I,N) reduces to a function of I alone.

print("  KEY: N = 1 + 2I²  (since Nc = 1 + 2(2I3)² for our 3 sectors)")
print("  So Δx₁₂ is a function of I = 2I3 ONLY (up to artifacts of 3 points).")
print()

# Fit quadratic in I: Δx_k = α + β·I + γ·I²
from numpy.linalg import solve

Ivals = np.array([0, 1, -1])
V = np.column_stack([np.ones(3), Ivals, Ivals**2])  # Vandermonde-like

for k in coords:
    dx_vals = np.array([float(Dx12[(s,k)]) for s in sectors])
    coeffs = solve(V, dx_vals)
    # α, β, γ
    alpha, beta, gamma = coeffs
    alpha_f = F(alpha).limit_denominator(100)
    beta_f = F(beta).limit_denominator(100)
    gamma_f = F(gamma).limit_denominator(100)
    print(f"  Δx₁₂({k}) = {float(alpha_f)} + {float(beta_f)}·I + {float(gamma_f)}·I²")
    print(f"          = {alpha_f} + ({beta_f})·I + ({gamma_f})·I²")
    # Verify
    for idx, s in enumerate(sectors):
        I = Ivals[idx]
        pred = alpha + beta*I + gamma*I**2
        actual = float(Dx12[(s,k)])
        match = "✓" if abs(pred - actual) < 1e-10 else "✗"
        print(f"          {s}: I={I:>+2d} → {pred:>+6.2f} (actual {actual:>+6.2f}) {match}")
    print()

# ── Now the same for Δx₂₃ ────────────────────────────────────────────
print("  Δx₂₃ as quadratic in I:")
for k in coords:
    dx_vals = np.array([float(Dx23[(s,k)]) for s in sectors])
    coeffs = solve(V, dx_vals)
    alpha, beta, gamma = coeffs
    alpha_f = F(alpha).limit_denominator(100)
    beta_f = F(beta).limit_denominator(100)
    gamma_f = F(gamma).limit_denominator(100)
    print(f"  Δx₂₃({k}) = {alpha_f} + ({beta_f})·I + ({gamma_f})·I²")

print()

# ── PRIME FACTORIZATION of n values ──────────────────────────────────
print("="*70)
print("  PRIME FACTORIZATION OF n")
print("="*70)
print()

def factorize(n):
    if n == 0: return "0"
    if n < 0: return "-" + factorize(-n)
    if n == 1: return "1"
    factors = []
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors.append(d)
            temp //= d
        d += 1
    if temp > 1:
        factors.append(temp)
    return "×".join(str(f) for f in factors)

print(f"  {'key':>12s}  {'n':>5s}  factorization   v₂(n)  odd part")
for s in sectors:
    for k in coords:
        n = int(N[(s,k)])
        abs_n = abs(n)
        tag = " [FLAT]" if (s,k) in flat else ""
        # 2-adic valuation
        v2 = 0
        temp = abs_n
        while temp > 0 and temp % 2 == 0:
            v2 += 1
            temp //= 2
        odd = abs_n // (2**v2) if abs_n > 0 else 0
        sign = "+" if n > 0 else "-"
        print(f"  {s:>6s}/{k}:  {n:>+5d}  {sign}{factorize(abs_n):>12s}  "
              f"   {v2}      {odd}{tag}")

print()
print("  2-adic valuation v₂(n) matrix:")
print(f"  {'':>8s}     p    q    r")
for s in sectors:
    row = ""
    for k in coords:
        n = int(N[(s,k)])
        abs_n = abs(n)
        v2 = 0
        temp = abs_n
        while temp > 0 and temp % 2 == 0:
            v2 += 1
            temp //= 2
        tag = "*" if (s,k) in flat else " "
        row += f"   {v2:>2d}{tag} "
    print(f"  {s:>8s}  {row}")

print()
print("  Odd part matrix (n / 2^v₂(n)):")
print(f"  {'':>8s}     p    q    r")
for s in sectors:
    row = ""
    for k in coords:
        n = int(N[(s,k)])
        sign = 1 if n > 0 else -1
        abs_n = abs(n)
        v2 = 0
        temp = abs_n
        while temp > 0 and temp % 2 == 0:
            v2 += 1
            temp //= 2
        odd = sign * (abs_n // (2**v2)) if abs_n > 0 else 0
        tag = "*" if (s,k) in flat else " "
        row += f"  {odd:>+3d}{tag} "
    print(f"  {s:>8s}  {row}")

print()
print("  OBSERVATION: the odd part is ±1 for all FREE coords")
print("  and is +3 for the one problematic flat coord (lepton/r: 6=2×3)")
print()

# ── WHERE DOES THE 3 COME FROM? ──────────────────────────────────────
print("="*70)
print("  WHY lepton/r HAS A FACTOR OF 3")
print("="*70)
print()
print("  For flat coords: n = -4a")
print("  lepton/r: a = -3/2  →  n = -4(-3/2) = 6 = 2×3")
print("  up/q:     a = 1/2   →  n = -4(1/2)  = -2")  
print("  down/p:   a = -1    →  n = -4(-1)   = 4")
print()
print("  The 3 in n=6 comes from |a| = 3/2.")
print("  The r-coordinate curvature has 3 in numerator (from ln3 basis).")
print()

# ── Is there a SINGLE formula? ───────────────────────────────────────
print("="*70)
print("  UNIFIED VIEW: n = 2Δx₁₂ = quadratic in I = 2I₃")
print("="*70)
print()
print("  The n values are NOT subject to two different rules.")
print("  ALL 9 are determined by the SAME generating function:")
print("    n_k = 2[F_k(I, g=2) − F_k(I, g=1)]")
print()
print("  The 'power of 2' pattern in 6 free coords is a CONSEQUENCE")
print("  of the specific coefficients of F_k, not an input.")
print("  The flatness constraint doesn't change the rule — it constrains")
print("  the generating function coefficients so that Δx₂₃ = 0.")
print()

# ── Check: Δx₁₂ polynomials in I, evaluated ─────────────────────────
# We found:
# Δx₁₂(p) = α_p + β_p·I + γ_p·I²
# etc.
# The n values are 2× these. Let's see if the polynomial coefficients
# have nice structure.

print("  Polynomial coefficients (Δx₁₂ = α + β·I + γ·I²):")
print(f"  {'coord':>5s}  {'α':>8s}  {'β':>8s}  {'γ':>8s}  {'2α':>5s}  {'2β':>5s}  {'2γ':>5s}")
for k in coords:
    dx_vals = np.array([float(Dx12[(s,k)]) for s in sectors])
    coeffs = solve(V, dx_vals)
    alpha, beta, gamma = coeffs
    af = F(alpha).limit_denominator(100)
    bf = F(beta).limit_denominator(100)
    gf = F(gamma).limit_denominator(100)
    print(f"  {k:>5s}  {str(af):>8s}  {str(bf):>8s}  {str(gf):>8s}  "
          f"{int(2*alpha):>5d}  {int(2*beta):>5d}  {int(2*gamma):>5d}")

print()
print("  So the n polynomial is n_k(I) = 2α_k + 2β_k·I + 2γ_k·I²:")
print()
for k in coords:
    dx_vals = np.array([float(Dx12[(s,k)]) for s in sectors])
    coeffs = solve(V, dx_vals)
    a2 = F(2*coeffs[0]).limit_denominator(100)
    b2 = F(2*coeffs[1]).limit_denominator(100) 
    g2 = F(2*coeffs[2]).limit_denominator(100)
    print(f"  n_{k}(I) = {a2} + ({b2})·I + ({g2})·I²")
    for idx, s in enumerate(sectors):
        I = Ivals[idx]
        val = int(2*coeffs[0] + 2*coeffs[1]*I + 2*coeffs[2]*I**2)
        actual = int(N[(s,k)])
        print(f"    I={I:>+2d} ({s:>6s}): n={val:>+3d}  (actual {actual:>+3d})")
    print()

# ── The 2-adic valuation of the polynomial ───────────────────────────
print("="*70)
print("  2-ADIC STRUCTURE OF n_k(I)")
print("="*70)
print()
print("  Why are 8 of 9 values powers of 2?")
print()
print("  n_p(I) = -1 + 5I/2 + 5I²/2")
print("    I=0: -1        v₂ = 0")
print("    I=±1: -1±5/2+5/2 = 4 or -6... wait")
print()

# Let me recompute carefully with fractions
for k in coords:
    dx_vals = [Dx12[(s,k)] for s in sectors]
    # Solve: α + β·0 + γ·0 = dx[lepton]  → α = dx[lepton]
    #        α + β·1 + γ·1 = dx[up]
    #        α - β·1 + γ·1 = dx[down]
    alpha = dx_vals[0]  # I=0
    # up: α + β + γ = dx_up
    # down: α - β + γ = dx_down
    # β = (dx_up - dx_down)/2
    # γ = (dx_up + dx_down)/2 - α
    beta = (dx_vals[1] - dx_vals[2]) / 2
    gamma = (dx_vals[1] + dx_vals[2]) / 2 - alpha
    
    print(f"  Δx₁₂({k}): α={alpha}, β={beta}, γ={gamma}")
    print(f"  n_{k}(I) = {2*alpha} + ({2*beta})·I + ({2*gamma})·I²")
    print(f"    I= 0: n = {2*alpha}")
    print(f"    I=+1: n = {2*alpha + 2*beta + 2*gamma} = {2*(alpha+beta+gamma)}")
    print(f"    I=-1: n = {2*alpha - 2*beta + 2*gamma} = {2*(alpha-beta+gamma)}")
    
    # Factor each
    for I_val, s in zip([0, 1, -1], sectors):
        n_val = 2*(alpha + beta*I_val + gamma*I_val**2)
        print(f"      {s}: n = {n_val} = {factorize(abs(int(n_val)))}", end="")
        if (s,k) in flat:
            print(" [FLAT]")
        else:
            print()
    print()
