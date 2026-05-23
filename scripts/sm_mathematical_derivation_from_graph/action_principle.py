#!/usr/bin/env python3
"""
action_principle.py — Variational Principle for the Mass Tensor
================================================================

Given 27 parabolic coefficients (a,b,c) × 3 coords × 3 sectors with
18 free parameters after constraints, we search for a scalar functional
S whose critical point reproduces the observed values.

Approach:
  The generation-g trajectory x_k(g) = a_k·g² + b_k·g + c_k is treated
  as a curve in 3D (p,q,r)-space.  Natural actions:

  §1  Kinetic action:   S_K = Σ_s ∫₁³ ½|ẋ_s|² dg
  §2  Curvature action: S_A = Σ_s |a_s|²
  §3  Mixed action:     S_M = Σ_s (α|a_s|² + β a_s·b_s + γ|b_s|² + ...)
  §4  Minimum-norm:     do the 18 free params minimize a norm subject to constraints?
  §5  Full variational: construct Euler-Lagrange for L(x,ẋ,s) and check consistency
  §6  Matrix invariants: det, trace, eigenvalues of A/B/C matrices
  §7  Sector coupling:  cross-sector terms that could unify the action

Author: SS-PEG program
Date: April 2026
"""

import numpy as np
from fractions import Fraction
from collections import OrderedDict
from itertools import product
import sys, os

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "action_principle_output.md")
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
# DATA (from parameter_reduction.py)
# ════════════════════════════════════════════════════════════════

F = Fraction

PARTICLES = OrderedDict([
    ("e",   (F(0),    F(0),    F(0))),
    ("mu",  (F(-1,2), F(-4),   F(3))),
    ("tau", (F(1),    F(-7),   F(3))),
    ("u",   (F(-3,2), F(-6),   F(-1))),
    ("c",   (F(1,2),  F(-7),   F(3))),
    ("t",   (F(6),    F(-7),   F(4))),
    ("d",   (F(-1,2), F(-8),   F(-2))),
    ("s",   (F(3,2),  F(-7),   F(0))),
    ("b",   (F(3,2),  F(-6),   F(4))),
])

SECTORS = OrderedDict([
    ("lepton", ["e", "mu", "tau"]),
    ("up",     ["u", "c", "t"]),
    ("down",   ["d", "s", "b"]),
])

FLAT_COORD = {"lepton": 2, "up": 1, "down": 0}  # index: p=0, q=1, r=2

def fit_parabola(x1, x2, x3):
    x1, x2, x3 = F(x1), F(x2), F(x3)
    a = (x1 - 2*x2 + x3) / 2
    b = (-5*x1 + 8*x2 - 3*x3) / 2
    c = 3*x1 - 3*x2 + x3
    return a, b, c

# Build coefficient arrays: A[s][k], B[s][k], C[s][k] as Fractions
A, B, CC = {}, {}, {}
for sname, names in SECTORS.items():
    coords = [PARTICLES[n] for n in names]
    A[sname], B[sname], CC[sname] = [], [], []
    for k in range(3):
        a, b, c = fit_parabola(coords[0][k], coords[1][k], coords[2][k])
        A[sname].append(a)
        B[sname].append(b)
        CC[sname].append(c)

# Float versions for numerical computation
Af = {s: np.array([float(x) for x in A[s]]) for s in SECTORS}
Bf = {s: np.array([float(x) for x in B[s]]) for s in SECTORS}
Cf = {s: np.array([float(x) for x in CC[s]]) for s in SECTORS}

sector_names = ["lepton", "up", "down"]

def print_header(n, title):
    print(f"\n{'=' * 72}")
    print(f"  §{n}. {title}")
    print(f"{'=' * 72}\n")


# ════════════════════════════════════════════════════════════════
# §1. KINETIC ACTION: S_K = Σ_s ∫₁³ ½|ẋ|² dg
# ════════════════════════════════════════════════════════════════

print_header(1, "KINETIC ACTION")

print("""  For x_k(g) = a_k g² + b_k g + c_k, the velocity is ẋ_k(g) = 2a_k g + b_k.
  
  S_K = Σ_s Σ_k ∫₁³ ½(2a_k g + b_k)² dg
      = Σ_s Σ_k ½[4a_k² ∫g² + 4a_k b_k ∫g + b_k² ∫1] dg
  
  ∫₁³ g² dg = 26/3,  ∫₁³ g dg = 4,  ∫₁³ 1 dg = 2
  
  So: S_K = Σ_s Σ_k [ (52/3) a_k² + 8 a_k b_k + b_k² ]
""")

# Compute per-sector kinetic action
S_K_sectors = {}
S_K_total = F(0)
for s in sector_names:
    val = F(0)
    for k in range(3):
        a, b = A[s][k], B[s][k]
        val += F(52,3)*a**2 + 8*a*b + b**2
    S_K_sectors[s] = val
    S_K_total += val

print(f"  Per-sector kinetic action:")
for s in sector_names:
    print(f"    S_K({s}) = {S_K_sectors[s]} = {float(S_K_sectors[s]):.6f}")
print(f"    S_K(total) = {S_K_total} = {float(S_K_total):.6f}")

# ════════════════════════════════════════════════════════════════
# §2. CURVATURE / POTENTIAL / OFFSET ACTIONS
# ════════════════════════════════════════════════════════════════

print_header(2, "CURVATURE, VELOCITY-NORM, AND OFFSET ACTIONS")

# S_A = Σ |a|², S_B = Σ |b|², S_C = Σ |c|²
S_A = sum(sum(A[s][k]**2 for k in range(3)) for s in sector_names)
S_B = sum(sum(B[s][k]**2 for k in range(3)) for s in sector_names)
S_C = sum(sum(CC[s][k]**2 for k in range(3)) for s in sector_names)
S_AB = sum(sum(A[s][k]*B[s][k] for k in range(3)) for s in sector_names)
S_AC = sum(sum(A[s][k]*CC[s][k] for k in range(3)) for s in sector_names)
S_BC = sum(sum(B[s][k]*CC[s][k] for k in range(3)) for s in sector_names)

print(f"  S_A  = Σ |a|²  = {S_A} = {float(S_A):.6f}")
print(f"  S_B  = Σ |b|²  = {S_B} = {float(S_B):.6f}")
print(f"  S_C  = Σ |c|²  = {S_C} = {float(S_C):.6f}")
print(f"  S_AB = Σ a·b   = {S_AB} = {float(S_AB):.6f}")
print(f"  S_AC = Σ a·c   = {S_AC} = {float(S_AC):.6f}")
print(f"  S_BC = Σ b·c   = {S_BC} = {float(S_BC):.6f}")

# Check: is S_AB = -5 S_A?  (would mean b = -5a universally)
print(f"\n  S_AB / S_A = {float(S_AB)/float(S_A):.6f}  (would be -5 if b = -5a universally)")

# ════════════════════════════════════════════════════════════════
# §3. MATRIX INVARIANTS
# ════════════════════════════════════════════════════════════════

print_header(3, "MATRIX INVARIANTS")

print("""  Form 3×3 matrices from coefficient vectors:
    A_{sk} = a_k(sector s),  B_{sk} = b_k(s),  C_{sk} = c_k(s)
  where rows = sectors (ℓ, u, d), cols = coordinates (p, q, r).
""")

Am = np.array([[float(A[s][k]) for k in range(3)] for s in sector_names])
Bm = np.array([[float(B[s][k]) for k in range(3)] for s in sector_names])
Cm = np.array([[float(CC[s][k]) for k in range(3)] for s in sector_names])

for label, M in [("A (curvature)", Am), ("B (velocity)", Bm), ("C (offset)", Cm)]:
    det = np.linalg.det(M)
    tr = np.trace(M)
    eigs = np.linalg.eigvals(M)
    svs = np.linalg.svd(M, compute_uv=False)
    print(f"  {label}:")
    print(f"    det = {det:.6f}")
    print(f"    tr  = {tr:.6f}")
    print(f"    eigenvalues = {', '.join(f'{e.real:.4f}{'+' if e.imag>=0 else ''}{e.imag:.4f}i' for e in eigs)}")
    print(f"    singular values = {', '.join(f'{v:.4f}' for v in svs)}")
    print(f"    Frobenius norm = {np.linalg.norm(M, 'fro'):.6f}")
    print()

# Check A + B + C  (this is x(1) for each sector)
print(f"  A + B + C = x(g=1):")
X1 = Am + Bm + Cm
for i, s in enumerate(sector_names):
    print(f"    {s}: {X1[i]}")
    
# Check 4A + 2B + C  (this is x(2))
print(f"\n  4A + 2B + C = x(g=2):")
X2 = 4*Am + 2*Bm + Cm
for i, s in enumerate(sector_names):
    print(f"    {s}: {X2[i]}")

# Check 9A + 3B + C  (this is x(3))
print(f"\n  9A + 3B + C = x(g=3):")
X3 = 9*Am + 3*Bm + Cm
for i, s in enumerate(sector_names):
    print(f"    {s}: {X3[i]}")

# ════════════════════════════════════════════════════════════════
# §4. THE FLATNESS-CONSTRAINED LAGRANGIAN
# ════════════════════════════════════════════════════════════════

print_header(4, "FLATNESS-CONSTRAINED LAGRANGIAN")

print("""  The flatness condition b_{k(s)} = -5 a_{k(s)} means the velocity
  vanishes at g* = 5/2 in the flat coordinate.  This is equivalent to:
  
    ẋ_{k(s)}(5/2) = 0  ⟺  x_{k(s)} has a turning point at g = 5/2
  
  If we define the Lagrangian L = ½|ẋ|² - V(x), the Euler-Lagrange eq is:
    ẍ_k = -∂V/∂x_k
  
  For parabolic motion: ẍ_k = 2a_k = const.
  So the "force" F_k = -∂V/∂x_k = 2a_k is constant — the potential is LINEAR:
  
    V(x) = -2 Σ_k a_k x_k + const = -2 a · x
  
  The action becomes:
    S[x] = ∫₁³ [½|ẋ|² + 2 a · x] dg
  
  But this is trivially satisfied — every parabola is a solution.
  The real question is: what determines the CURVATURE VECTOR a itself?
""")

# The force (acceleration) for each sector
print(f"  Force vectors F = 2a (constant per sector):")
for s in sector_names:
    fvec = [2*A[s][k] for k in range(3)]
    print(f"    F_{s} = ({fvec[0]}, {fvec[1]}, {fvec[2]})")

# ════════════════════════════════════════════════════════════════
# §5. SECOND-LEVEL ACTION: WHAT DETERMINES a(s)?
# ════════════════════════════════════════════════════════════════

print_header(5, "SECOND-LEVEL ACTION: WHAT DETERMINES a(s)?")

print("""  The curvature matrix A_{sk} = a_k(s) is a 3×3 matrix.
  We need a SECOND-LEVEL action Ŝ[A] whose variation δŜ/δA = 0
  yields the observed curvature.  Similarly for B, C (or B, C derived from A).
  
  Candidate functionals:
    Ŝ₁ = Tr(AᵀA) = Frobenius norm² = Σ a_k²
    Ŝ₂ = det(A)
    Ŝ₃ = Tr(A) + λ det(A)  (characteristic polynomial)
    Ŝ₄ = Tr(AᵀA) + μ Tr(A²)  (with symmetry)
    Ŝ₅ = Σ_s (Σ_k a_k)²  (trace per sector)
    Ŝ₆ = Σ_s a_s · diag(w) · a_s  (weighted norm with metric)
""")

# Compute key invariants of A (exact fractions)
AF = [[A[s][k] for k in range(3)] for s in sector_names]
BF = [[B[s][k] for k in range(3)] for s in sector_names]
CF = [[CC[s][k] for k in range(3)] for s in sector_names]

# Trace
trA = sum(AF[i][i] for i in range(3))  # a_p^ℓ + a_q^u + a_r^d
print(f"  Tr(A) = a_p^ℓ + a_q^u + a_r^d = {AF[0][0]} + {AF[1][1]} + {AF[2][2]} = {trA}")

# Anti-trace (flat-coordinate diagonal)
anti_trA = AF[0][2] + AF[1][1] + AF[2][0]  # a_r^ℓ + a_q^u + a_p^d (flat coords!)
print(f"  Flat-diagonal: a_r^ℓ + a_q^u + a_p^d = {AF[0][2]} + {AF[1][1]} + {AF[2][0]} = {anti_trA}")
print(f"  ← This is the sum of curvatures in FLAT coordinates!")

# Determinant (exact)
def det3_frac(M):
    return (M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
          - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
          + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]))

detA = det3_frac(AF)
detB = det3_frac(BF)
detC = det3_frac(CF)

print(f"\n  det(A) = {detA} = {float(detA):.6f}")
print(f"  det(B) = {detB} = {float(detB):.6f}")
print(f"  det(C) = {detC} = {float(detC):.6f}")

# Tr(A²)
def tr_sq(M):
    n = len(M)
    val = F(0)
    for i in range(n):
        for j in range(n):
            val += M[i][j] * M[j][i]
    return val

trA2 = tr_sq(AF)
trB2 = tr_sq(BF)

print(f"\n  Tr(A²) = {trA2} = {float(trA2):.6f}")
print(f"  Tr(B²) = {trB2} = {float(trB2):.6f}")

# Frobenius: Σ a_{sk}²
frobA = sum(AF[i][j]**2 for i in range(3) for j in range(3))
frobB = sum(BF[i][j]**2 for i in range(3) for j in range(3))

print(f"\n  ||A||_F² = Tr(AᵀA) = {frobA} = {float(frobA):.6f}")
print(f"  ||B||_F² = Tr(BᵀB) = {frobB} = {float(frobB):.6f}")

# Characteristic polynomial of A: λ³ - Tr(A)λ² + ½(Tr²-TrA²)λ - det = 0
p2 = (trA**2 - trA2) / 2
print(f"\n  Characteristic polynomial of A:")
print(f"    λ³ - ({trA})λ² + ({p2})λ - ({detA}) = 0")

# ════════════════════════════════════════════════════════════════
# §6. MINIMUM-NORM TEST
# ════════════════════════════════════════════════════════════════

print_header(6, "MINIMUM-NORM TEST")

print("""  Question: among all coefficient sets satisfying the 7 constraints
  (electron origin ×3, flatness ×3, q-unification ×1), do the physical
  values minimize the total Frobenius norm ||A||² + ||B||² + ||C||²?
  
  Strategy: parameterize the constraint manifold and check if the
  gradient of the norm vanishes at the physical point.
""")

# Set up the 27 variables as a vector θ
# θ = [a_p^ℓ, a_q^ℓ, a_r^ℓ, a_p^u, a_q^u, a_r^u, a_p^d, a_q^d, a_r^d,
#       b_p^ℓ, ..., b_r^d, c_p^ℓ, ..., c_r^d]

theta_phys = np.zeros(27)
for i, s in enumerate(sector_names):
    for k in range(3):
        theta_phys[i*3 + k] = float(A[s][k])         # a's: idx 0-8
        theta_phys[9 + i*3 + k] = float(B[s][k])      # b's: idx 9-17
        theta_phys[18 + i*3 + k] = float(CC[s][k])     # c's: idx 18-26

# Constraint matrix: G·θ = 0 (or = h for inhomogeneous)
# C1: electron origin: a_ℓ + b_ℓ + c_ℓ = 0 (for k=p,q,r)
#     → θ[0]+θ[9]+θ[18] = 0, θ[1]+θ[10]+θ[19] = 0, θ[2]+θ[11]+θ[20] = 0
# C2: flatness: b_{k(s)} = -5·a_{k(s)}
#     ℓ→r: θ[11] + 5·θ[2] = 0   (b_r^ℓ + 5·a_r^ℓ = 0)
#     u→q: θ[13] + 5·θ[4] = 0   (b_q^u + 5·a_q^u = 0)
#     d→p: θ[15] + 5·θ[6] = 0   (b_p^d + 5·a_p^d = 0)
# C3: q-unification: 3·(b_q^u - b_q^ℓ) + (c_q^u - c_q^ℓ) = 0
#     3·(θ[13]-θ[10]) + (θ[22]-θ[19]) = 0

n_constraints = 7
G = np.zeros((n_constraints, 27))
h = np.zeros(n_constraints)

# C1: electron origin (3 eqs)
for k in range(3):
    G[k, 0*3+k] = 1    # a_k^ℓ
    G[k, 9+0*3+k] = 1  # b_k^ℓ
    G[k, 18+0*3+k] = 1 # c_k^ℓ

# C2: flatness (3 eqs)
# ℓ→r (k=2): 5·a_r^ℓ + b_r^ℓ = 0
G[3, 2] = 5      # a_r^ℓ
G[3, 9+2] = 1    # b_r^ℓ
# u→q (k=1): 5·a_q^u + b_q^u = 0
G[4, 3+1] = 5    # a_q^u
G[4, 9+3+1] = 1  # b_q^u
# d→p (k=0): 5·a_p^d + b_p^d = 0
G[5, 6+0] = 5    # a_p^d
G[5, 9+6+0] = 1  # b_p^d

# C3: q-unification: 3(b_q^u - b_q^ℓ) + (c_q^u - c_q^ℓ) = 0
G[6, 9+3+1] = 3   # +3 b_q^u
G[6, 9+0+1] = -3  # -3 b_q^ℓ
G[6, 18+3+1] = 1  # +c_q^u
G[6, 18+0+1] = -1 # -c_q^ℓ

# Verify constraints are satisfied
residual = G @ theta_phys - h
print(f"  Constraint residual at physical point: {np.max(np.abs(residual)):.2e}")

# Project gradient of ||θ||² onto constraint manifold
# Gradient of ||θ||² = 2θ
# Projected gradient = (I - Gᵀ(GGᵀ)⁻¹G) · 2θ
grad = 2 * theta_phys
GGT_inv = np.linalg.inv(G @ G.T)
proj = np.eye(27) - G.T @ GGT_inv @ G
proj_grad = proj @ grad

print(f"\n  ||θ||² at physical point: {np.dot(theta_phys, theta_phys):.6f}")
print(f"  |projected gradient|: {np.linalg.norm(proj_grad):.6f}")
print(f"  → Physical point is {'NOT ' if np.linalg.norm(proj_grad) > 1e-10 else ''}a minimum of ||θ||²")

# Test various weighted norms: ||A||² + λ||B||² + μ||C||²
print(f"\n  Searching for weights (1, λ, μ) such that ||A||² + λ||B||² + μ||C||²")
print(f"  has a critical point at the physical values...")

# Weight matrix W = diag(1,...,1, λ,...,λ, μ,...,μ) for the 9+9+9 components
# Gradient of θᵀWθ = 2Wθ
# Projected gradient = (I - Gᵀ(GGᵀ)⁻¹G) · 2Wθ  -- but this changes with W!
# Actually: projected gradient = P · 2Wθ where P = I - Gᵀ(GWGᵀ)⁻¹GW ... no.
# The correct projection for minimum of θᵀWθ s.t. Gθ=0 is:
# 2Wθ = Gᵀ·λ  ⟹  θ = ½W⁻¹Gᵀλ  and  Gθ=0  ⟹  ½GW⁻¹Gᵀλ = 0.
# For a nontrivial solution of Gθ=h, we need θ = W⁻¹Gᵀ(GW⁻¹Gᵀ)⁻¹h + null(G) part.
# So the minimum-norm solution (in the W-norm) of Gθ=0 is θ=0.
# That's trivial. The constraints Gθ=0 INCLUDE θ=0 (all zeros → all particles at origin).
# So the minimum-norm test makes more sense if we have INHOMOGENEOUS constraints.

# Let's think differently. Fix the electron origin constraint as BOUNDARY (not norm),
# and ask: among all (a, B, C) with these constraints, what minimizes an action?

# Alternative: treat a as given (9 params), derive b,c from constraints + action
# Constraints determine 7 of the 18 other params (b,c). 
# Remaining 11 b,c params need another principle.

# Let's try: minimize S_K (kinetic action) with respect to a
# while holding b as determined by constraints where possible
print(f"\n  Per-sector kinetic action breakdown:")
for s in sector_names:
    a_vec = Af[s]
    b_vec = Bf[s]
    flat_k = FLAT_COORD[s]
    
    S_flat = float(F(52,3)) * a_vec[flat_k]**2 + 8*a_vec[flat_k]*b_vec[flat_k] + b_vec[flat_k]**2
    S_nonflat = 0
    for k in range(3):
        if k == flat_k:
            continue
        S_nonflat += float(F(52,3)) * a_vec[k]**2 + 8*a_vec[k]*b_vec[k] + b_vec[k]**2
    print(f"    {s}: S_flat = {S_flat:.4f}, S_nonflat = {S_nonflat:.4f}, total = {S_flat+S_nonflat:.4f}")
    # For flat coord: b = -5a, so S = (52/3)a² + 8a(-5a) + 25a² = (52/3 - 40 + 25)a² = (52/3 - 15)a² = 7/3 a²
    S_flat_check = float(F(7,3)) * a_vec[flat_k]**2
    print(f"           S_flat check (7/3·a²): {S_flat_check:.4f} {'✓' if abs(S_flat - S_flat_check) < 1e-10 else '✗'}")

# ════════════════════════════════════════════════════════════════
# §7. GENERATIONAL EULER-LAGRANGE
# ════════════════════════════════════════════════════════════════

print_header(7, "GENERATIONAL EULER-LAGRANGE ANALYSIS")

print("""  For a general Lagrangian L(x, ẋ, g) = T(ẋ) - V(x, g), the E-L eq is:
    ẍ_k = -∂V/∂x_k
  
  Since ẍ_k = 2a_k = const, V must be linear in x: V = -2a·x.
  
  BUT: the force a = a(s) depends on the sector, not on g.
  So V = V(x, s) = -2 a(s) · x.
  
  The deeper question: is there a UNIVERSAL potential V(x) valid for all sectors?
  → Check if a(s) = -½ ∇V evaluated at the sector's trajectory.
""")

# A universal V would need ∇V to be the same function of x for all sectors.
# Since a(s) varies, V cannot be a FIXED function of x unless x(g,s) parameterizes
# different regions where ∇V is approximately constant.

# Check: positions of gen-2 particles (midpoint of trajectory)
print(f"  Gen-2 positions (midpoint of each sector trajectory):")
for s in sector_names:
    x2 = 4*Af[s] + 2*Bf[s] + Cf[s]
    print(f"    {s}: x₂ = ({x2[0]:.2f}, {x2[1]:.2f}, {x2[2]:.2f}), |a| = {np.linalg.norm(Af[s]):.4f}")

# Check: is |a| × something related to position?
print(f"\n  Force magnitude vs position norm:")
for s in sector_names:
    for g in [1, 2, 3]:
        xg = g**2 * Af[s] + g * Bf[s] + Cf[s]
        print(f"    {s} g={g}: |x| = {np.linalg.norm(xg):.4f}, |F| = {2*np.linalg.norm(Af[s]):.4f}, |F|/|x| = {2*np.linalg.norm(Af[s])/max(np.linalg.norm(xg),1e-15):.4f}")

# ════════════════════════════════════════════════════════════════
# §8. CROSS-MATRIX RELATIONS: B = f(A)?
# ════════════════════════════════════════════════════════════════

print_header(8, "CROSS-MATRIX RELATIONS")

print("""  Can we express B as a function of A?
  Test: B = α·A + β·I + γ·P (where P is the flatness permutation matrix)
""")

# The flatness permutation maps: ℓ↔r, u↔q, d↔p
# As a matrix acting on coordinates: P_{kk'} permutes columns
# For row s, the flat coordinate k(s) gets b = -5a
# Non-flat coords are free.

# Test: B = M·A for some 3×3 matrix M (acting on coordinates/columns)
# B_{sk} = Σ_j M_{kj} A_{sj}
# This is B = A·Mᵀ in matrix notation
# So Mᵀ = A⁻¹ B  (if A is invertible)

print(f"  det(A) = {np.linalg.det(Am):.6f}")
if abs(np.linalg.det(Am)) > 1e-10:
    MT = np.linalg.solve(Am, Bm)
    M = MT.T
    print(f"\n  B = A · Mᵀ  where M =")
    for row in M:
        print(f"    [{', '.join(f'{x:8.4f}' for x in row)}]")
    
    residual = Am @ MT - Bm
    print(f"  Residual ||A·Mᵀ - B||_F = {np.linalg.norm(residual):.2e}")
    
    # Express M as exact fractions
    print(f"\n  M in exact fractions:")
    AF_mat = [[A[s][k] for k in range(3)] for s in sector_names]
    BF_mat = [[B[s][k] for k in range(3)] for s in sector_names]
    
    # Solve A·Mᵀ = B using Fraction arithmetic (Cramer's rule on columns)
    def solve_3x3_frac(M_in, rhs_col):
        dM = det3_frac(M_in)
        if dM == 0:
            return None
        sol = []
        for col in range(3):
            M2 = [[M_in[r][c] if c != col else rhs_col[r] for c in range(3)] for r in range(3)]
            sol.append(det3_frac(M2) / dM)
        return sol
    
    MT_frac = []
    for j in range(3):  # column j of Mᵀ = column j of B, solved via A
        rhs = [BF_mat[i][j] for i in range(3)]
        sol = solve_3x3_frac(AF_mat, rhs)
        MT_frac.append(sol)
    
    print(f"    Mᵀ (rows correspond to coords p,q,r):")
    for j in range(3):
        coord = ["p", "q", "r"][j]
        print(f"      Mᵀ[{coord},:] = [{MT_frac[j][0]}, {MT_frac[j][1]}, {MT_frac[j][2]}]")
    
    # M[i][j] = (Mᵀ)_{j,i} = (A⁻¹B)_{j,i}: note MT_frac[j] = col j of A⁻¹B
    M_frac = [[MT_frac[i][j] for j in range(3)] for i in range(3)]
    print(f"\n    M = (Mᵀ)ᵀ  (exact):")
    for i in range(3):
        coord = ["p", "q", "r"][i]
        print(f"      M[{coord},:] = [{M_frac[i][0]}, {M_frac[i][1]}, {M_frac[i][2]}]")
    
    # Properties of M
    trM = sum(M_frac[i][i] for i in range(3))
    detM = det3_frac(M_frac)
    print(f"\n    Tr(M) = {trM} = {float(trM):.6f}")
    print(f"    det(M) = {detM} = {float(detM):.6f}")

    # Verify M_frac matches numpy M = MT.T
    M_check = np.array([[float(M_frac[i][j]) for j in range(3)] for i in range(3)])
    print(f"    Consistency ||M_frac - M_numpy||: {np.linalg.norm(M_check - M):.2e}")

    # Similarly: C = A · Nᵀ ?
    print(f"\n  C = A · Nᵀ  where:")
    CF_mat = [[CC[s][k] for k in range(3)] for s in sector_names]
    NT_frac = []
    for j in range(3):
        rhs = [CF_mat[i][j] for i in range(3)]
        sol = solve_3x3_frac(AF_mat, rhs)
        NT_frac.append(sol)
    
    N_frac = [[NT_frac[i][j] for j in range(3)] for i in range(3)]
    print(f"    N = ")
    for i in range(3):
        coord = ["p", "q", "r"][i]
        print(f"      N[{coord},:] = [{N_frac[i][0]}, {N_frac[i][1]}, {N_frac[i][2]}]")
    
    trN = sum(N_frac[i][i] for i in range(3))
    detN = det3_frac(N_frac)
    print(f"\n    Tr(N) = {trN} = {float(trN):.6f}")
    print(f"    det(N) = {detN} = {float(detN):.6f}")

    # Verify: A·Nᵀ should equal C   (Nᵀ_{j,k} = M_{k,j} → need to reconstruct Nᵀ)
    NT_mat = np.linalg.solve(Am, Cm)
    print(f"    Residual ||A·Nᵀ - C||_F = {np.linalg.norm(Am @ NT_mat - Cm):.2e}")
    N_check = np.array([[float(N_frac[i][j]) for j in range(3)] for i in range(3)])
    print(f"    Consistency ||N_frac - N_numpy||: {np.linalg.norm(N_check - NT_mat.T):.2e}")

# ════════════════════════════════════════════════════════════════
# §9. B = A · M AND THE STRUCTURE OF M
# ════════════════════════════════════════════════════════════════

print_header(9, "THE VELOCITY-CURVATURE TRANSFER MATRIX M")

print("""  We found B = A · Mᵀ, i.e., the velocity matrix is determined by
  the curvature matrix via a universal 3×3 transfer matrix M.
  
  This means: b_k(s) = Σ_j M_{kj} a_j(s)
  
  So the velocity in coordinate k depends on curvatures in ALL coords j.
  The 9 free parameters in B are encoded in the 9 entries of M.
  
  Can we understand M's structure?
""")

if abs(np.linalg.det(Am)) > 1e-10:
    # Check: is M symmetric?
    is_sym = all(M_frac[i][j] == M_frac[j][i] for i in range(3) for j in range(3))
    print(f"  Is M symmetric? {is_sym}")
    
    # Antisymmetric part
    AS = [[M_frac[i][j] - M_frac[j][i] for j in range(3)] for i in range(3)]
    print(f"  Antisymmetric part M - Mᵀ:")
    for i in range(3):
        print(f"    [{AS[i][0]}, {AS[i][1]}, {AS[i][2]}]")
    
    # Diagonal elements
    print(f"\n  Diagonal of M: [{M_frac[0][0]}, {M_frac[1][1]}, {M_frac[2][2]}]")
    print(f"  ← The flat-coord diagonal entries should relate to -5 (flatness)")
    
    # For the flat coordinate k(s): b_{k(s)} = -5 a_{k(s)}
    # But from B = A·Mᵀ: b_k(s) = Σ_j a_j(s) · M_{kj}
    # The flatness condition constrains specific COMBINATIONS of a and M,
    # not M directly. Let's check what flatness implies for M.
    
    print(f"\n  Flatness check via M:")
    for s_idx, s in enumerate(sector_names):
        fk = FLAT_COORD[s]
        coord_name = ["p", "q", "r"][fk]
        # b_{fk}(s) = Σ_j M_{fk,j} a_j(s) should equal -5 a_{fk}(s)
        lhs = sum(float(M_frac[fk][j]) * float(AF_mat[s_idx][j]) for j in range(3))
        rhs = -5 * float(AF_mat[s_idx][fk])
        print(f"    {s} ({coord_name}-flat): Σ_j M_{coord_name}j · a_j = {lhs:.6f}, -5·a_{coord_name} = {rhs:.6f} {'✓' if abs(lhs-rhs)<1e-10 else '✗'}")

    # Check: does M have eigenvalue -5?
    M_float = np.array([[float(M_frac[i][j]) for j in range(3)] for i in range(3)])
    eigvals_M = np.linalg.eigvals(M_float)
    print(f"\n  Eigenvalues of M: {', '.join(f'{e.real:.6f}{'+' if e.imag>=0 else ''}{e.imag:.6f}i' for e in eigvals_M)}")
    print(f"  Does M have eigenvalue -5? {any(abs(e + 5) < 1e-6 for e in eigvals_M)}")
    
    # Similarly for N (C = A·Nᵀ)
    N_float = np.array([[float(N_frac[i][j]) for j in range(3)] for i in range(3)])
    eigvals_N = np.linalg.eigvals(N_float)
    print(f"\n  Eigenvalues of N: {', '.join(f'{e.real:.6f}{'+' if e.imag>=0 else ''}{e.imag:.6f}i' for e in eigvals_N)}")

    # M · N relation?
    MN = M_float @ N_float
    NM = N_float @ M_float
    print(f"\n  M · N:")
    for row in MN:
        print(f"    [{', '.join(f'{x:8.4f}' for x in row)}]")
    print(f"  N · M:")
    for row in NM:
        print(f"    [{', '.join(f'{x:8.4f}' for x in row)}]")
    print(f"  Commutator [M,N] = M·N - N·M:")
    comm = MN - NM
    for row in comm:
        print(f"    [{', '.join(f'{x:8.4f}' for x in row)}]")
    print(f"  ||[M,N]|| = {np.linalg.norm(comm):.6f}")

# ════════════════════════════════════════════════════════════════
# §10. EFFECTIVE ACTION IN TERMS OF A AND M
# ════════════════════════════════════════════════════════════════

print_header(10, "EFFECTIVE ACTION S[A, M, N]")

print("""  Given B = A·Mᵀ and C = A·Nᵀ, the 27 parameters reduce to:
    9 (curvature A) + 9 (transfer M) + 9 (transfer N) = 27
  but with constraints becoming conditions on M and N.
  
  The kinetic action becomes:
    S_K = Σ_s ∫₁³ ½|2A_s g + B_s|² dg
        = Σ_s [(52/3)|a_s|² + 8 a_s·(A·Mᵀ)_s + |A·Mᵀ_s|²]
  
  If we further express A in terms of quantum numbers:
    A_{sk} = α_k · (2I₃) + β_k · Nc + γ_k
  then the entire action is a function of (α, β, γ, M, N).
""")

# Total DOF accounting:
# A is parameterized by (α, β, γ) × 3 coords: 9 params (but a_q, a_r have α=-β: 7 params)
# M and N: 9 + 9 = 18 params
# But constraints connect them:
# Electron origin: c_ℓ = -a_ℓ - b_ℓ → C = A·Nᵀ evaluated at ℓ: A_ℓ(I+Mᵀ+Nᵀ) = 0
#   Since A_ℓ ≠ 0, this means (I+Mᵀ+Nᵀ) must annihilate A_ℓ, or more precisely:
#   a_ℓ + A_ℓ·M_col + A_ℓ·N_col = 0 for each coord column.

# Let's check: I + M + N (note: C = A·Nᵀ, B = A·Mᵀ, so x(1) = A + B + C = A(I + Mᵀ + Nᵀ))
# Electron origin: x_ℓ(1) = 0 → A_ℓ(I + Mᵀ + Nᵀ) = 0 → (I + Mᵀ + Nᵀ) applied to A's row ℓ = 0.

IMN = np.eye(3) + M_float.T + N_float.T
print(f"  I + Mᵀ + Nᵀ (should annihilate A_ℓ):")
for row in IMN:
    print(f"    [{', '.join(f'{x:8.4f}' for x in row)}]")

# Check: A_ℓ · (I + Mᵀ + Nᵀ)  (row × matrix → row vector)
a_ell = Am[0]
result = a_ell @ IMN
print(f"\n  A_ℓ · (I + Mᵀ + Nᵀ) = [{', '.join(f'{x:.6f}' for x in result)}]")
print(f"  (should be zero: {np.allclose(result, 0)})")

# For other sectors:
for i, s in enumerate(["up", "down"]):
    result = Am[i+1] @ IMN
    print(f"  A_{s} · (I+Mᵀ+Nᵀ) = [{', '.join(f'{x:.6f}' for x in result)}]  (= x_{s}(g=1))")

# ════════════════════════════════════════════════════════════════
# §11. SEARCH FOR ACTION WHOSE CRITICAL POINT IS M
# ════════════════════════════════════════════════════════════════

print_header(11, "VARIATIONAL SEARCH FOR M")

print("""  Given that A is fixed by quantum numbers, we ask:
  Is there a function f(M) such that ∂f/∂M_{ij} = 0 at the physical M?
  
  Test candidates:
    f₁ = Tr(M)          → trivially not (Tr doesn't couple entries)
    f₂ = det(M)          → ∂det/∂M = cofactor matrix
    f₃ = Tr(M²)          → ∂/∂M = 2M
    f₄ = Tr(M²) + λ det(M)  → combines
    f₅ = Tr(MᵀM)         → ∂/∂M = 2M
    f₆ = ||AM||² = Tr(AᵀA MᵀM) → couples A and M
""")

# f₃ = Tr(M²): critical point requires M = 0. Not our case.
print(f"  Tr(M) = {np.trace(M_float):.6f}")
print(f"  Tr(M²) = {np.trace(M_float @ M_float):.6f}")
print(f"  det(M) = {np.linalg.det(M_float):.6f}")

# Check: is M a critical point of Tr(M²) + λ det(M)?
# ∂/∂M [Tr(M²) + λ det(M)] = 2M + λ cof(M) = 0
# → M = -λ/2 · cof(M) = -λ/2 · det(M) · M⁻ᵀ
# → M · Mᵀ = -λ/2 · det(M) · I
# → M·Mᵀ must be proportional to I
MMT = M_float @ M_float.T
print(f"\n  M · Mᵀ:")
for row in MMT:
    print(f"    [{', '.join(f'{x:8.4f}' for x in row)}]")
print(f"  M · Mᵀ proportional to I? {np.allclose(MMT, MMT[0,0]*np.eye(3), atol=0.1)}")

# How about M as critical point of S_K[M] = Σ_s ∫|2a_s g + (a_s·M)|² dg
# where b_s = Σ_j M_kj a_j(s)?
# ∂S_K/∂M_{kj} = Σ_s ∂/∂M_{kj} [(52/3)|a_s|² + 8 Σ_m a_m(s) M_{mk} a_k(s) ... ]
# This is getting complex. Let me compute it numerically.

print(f"\n  Testing: is M a saddle point of S_K[M] (kinetic action as function of M)?")
print(f"  Holding A fixed, varying M:")

def kinetic_action_of_M(M_mat, A_mat):
    """Compute S_K = Σ_s Σ_k [(52/3) a_k² + 8 a_k b_k + b_k²] where b = a·Mᵀ."""
    S = 0
    for s_idx in range(3):
        a_s = A_mat[s_idx]
        b_s = a_s @ M_mat.T
        for k in range(3):
            S += (52/3) * a_s[k]**2 + 8 * a_s[k] * b_s[k] + b_s[k]**2
    return S

# Numerical gradient
eps = 1e-7
S0 = kinetic_action_of_M(M_float, Am)
grad_M = np.zeros((3,3))
for i in range(3):
    for j in range(3):
        dM = np.zeros((3,3))
        dM[i,j] = eps
        S_plus = kinetic_action_of_M(M_float + dM, Am)
        grad_M[i,j] = (S_plus - S0) / eps

print(f"  S_K(physical M) = {S0:.6f}")
print(f"  ∂S_K/∂M:")
for row in grad_M:
    print(f"    [{', '.join(f'{x:8.4f}' for x in row)}]")
print(f"  |∇_M S_K| = {np.linalg.norm(grad_M):.6f}")
print(f"  → M is {'NOT ' if np.linalg.norm(grad_M) > 0.01 else ''}a critical point of S_K")

# What if we add a constraint penalty? S = S_K + λ·S_constraint
# The flatness constraints on M are: for each (s, flat_k),
# Σ_j M_{fk,j} a_j(s) = -5 a_{fk}(s)
# This constrains M along 3 directions in the 9D space.

# Compute constrained gradient (project out flatness-constraint directions)
# The flatness constraint is: for (s, fk): e_{fk}ᵀ M a_s + 5 a_{fk}(s) = 0
# Linearized: δ(e_{fk}ᵀ M a_s) = e_{fk}ᵀ δM a_s = Σ_j δM_{fk,j} a_j(s)
# This constrains 3 linear combos of the 9 entries of M.

# Build constraint Jacobian for M (3 × 9 matrix)
J_flat = np.zeros((3, 9))
for c_idx, (s_idx, s) in enumerate([(0, "lepton"), (1, "up"), (2, "down")]):
    fk = FLAT_COORD[s]
    for j in range(3):
        J_flat[c_idx, fk*3 + j] = Am[s_idx, j]

# Also add q-unification constraint on M (acts through B)
# q_u(3) = q_ℓ(3) → 9a_q^u + 3b_q^u + c_q^u = 9a_q^ℓ + 3b_q^ℓ + c_q^ℓ
# Since a_q^u = a_q^ℓ, this is 3(b_q^u - b_q^ℓ) + (c_q^u - c_q^ℓ) = 0
# In terms of M: b_q(s) = Σ_j M_{1j} a_j(s)
# and C = A·Nᵀ (N is a separate unknown here)
# So this constraint involves both M and N. Skip for now.

# Also add electron origin constraint: A_ℓ(I + Mᵀ + Nᵀ) = 0
# Again involves both M and N.

# Project gradient of S_K w.r.t. M onto the constraint manifold
grad_flat = grad_M.flatten()
rank_J = np.linalg.matrix_rank(J_flat)
print(f"\n  Flatness constraint rank: {rank_J}")

if rank_J > 0:
    JJT_inv = np.linalg.pinv(J_flat @ J_flat.T)
    P_perp = np.eye(9) - J_flat.T @ JJT_inv @ J_flat
    proj_grad_flat = P_perp @ grad_flat
    print(f"  |projected ∇_M S_K| (after flatness constraints) = {np.linalg.norm(proj_grad_flat):.6f}")

# ════════════════════════════════════════════════════════════════
# §12. THE FULL 27-PARAMETER ACTION LANDSCAPE
# ════════════════════════════════════════════════════════════════

print_header(12, "FULL ACTION LANDSCAPE")

print("""  The kinetic action S_K = Σ_s ∫₁³ ½|ẋ|² dg is a QUADRATIC function of (a, b).
  Since c does not appear (only ẋ = 2ag + b depends on a,b), S_K is independent of c.
  
  So S_K determines (a, b) = 18 of the 27 params, while c is determined by:
    - Electron origin: 3 eqs
    - Bridge/q-unification: further constraints
  
  Remaining: the action S_K(a,b) subject to flatness constraints.
""")

# Core insight: S_K = Σ_s Σ_k [(52/3)a² + 8ab + b²]
# Per coordinate: f(a,b) = (52/3)a² + 8ab + b²
# ∂f/∂b = 8a + 2b = 0 → b = -4a (for UNCONSTRAINED minimum!)
# But flatness says b_flat = -5a. And the unconstrained minimum says b = -4a.
# The non-flat coords have free b.

print(f"  Unconstrained minimum of f(a,b) = (52/3)a² + 8ab + b²:")
print(f"    ∂f/∂b = 8a + 2b = 0 → b* = -4a (turning point at g* = 2)")
print(f"    ∂f/∂a = (104/3)a + 8b = 0 → b* = -(13/3)a (turning point at g* = 13/6)")
print(f"    Joint optimum: a = b = 0 (trivial)")
print(f"")
print(f"  With fixed a: b* = -4a minimizes S_K → g* = -b/(2a) = 2 (exactly gen-2!)")
print(f"  Physical flat coords have b = -5a → g* = 5/2 (between gen 2 and 3)")
print(f"  Physical non-flat coords have various b/a ratios:")

for s in sector_names:
    fk = FLAT_COORD[s]
    for k in range(3):
        if k == fk:
            continue
        ratio = float(B[s][k] / A[s][k]) if A[s][k] != 0 else float('inf')
        gstar = -ratio / 2
        deviation = ratio - (-4)  # deviation from b = -4a
        print(f"    {s}/{['p','q','r'][k]}: b/a = {B[s][k]/A[s][k] if A[s][k] != 0 else '∞':>8}, g* = {gstar:.4f}, deviation from -4: {B[s][k]/A[s][k]+4 if A[s][k] != 0 else '–'}")

print(f"""
  Key insight: if b = -4a minimized S_K (turning point at g=2),
  and flatness forces b = -5a (turning point at g=5/2),
  then the flat coordinate is PUSHED from the kinetic optimum.
  
  The "cost" of flatness: ΔS = f(a,-5a) - f(a,-4a) for each flat coord
""")

for s in sector_names:
    fk = FLAT_COORD[s]
    a_val = float(A[s][fk])
    f_opt = (52/3)*a_val**2 + 8*a_val*(-4*a_val) + (-4*a_val)**2  # = (52/3-32+16)a² = (4/3)a²
    f_flat = (52/3)*a_val**2 + 8*a_val*(-5*a_val) + (-5*a_val)**2  # = (52/3-40+25)a² = (7/3)a²
    delta = f_flat - f_opt
    print(f"    {s}: a_{['p','q','r'][fk]} = {A[s][fk]}, ΔS = {F(7,3)*A[s][fk]**2 - F(4,3)*A[s][fk]**2} = a²·(7/3 - 4/3) = a² = {A[s][fk]**2}")

total_flatness_cost = sum(A[s][FLAT_COORD[s]]**2 for s in sector_names)
print(f"    Total flatness cost: Σ a²_flat = {total_flatness_cost} = {float(total_flatness_cost):.6f}")

# ════════════════════════════════════════════════════════════════
# §13. SUMMARY AND INTERPRETATION
# ════════════════════════════════════════════════════════════════

print_header(13, "DEVIATION FROM KINETIC OPTIMUM")

print("""  Key insight: with fixed a, the kinetic-optimal velocity is b* = -4a
  (turning point at g = 2).  Define the deviation:
    δb_k(s) = b_k(s) - (-4 a_k(s)) = b_k(s) + 4 a_k(s)
  
  For flat coords: δb_flat = -5a + 4a = -a (deviation = -a).
  For non-flat: δb is free.
""")

# Compute deviation vectors
print(f"  Deviation vectors δb = b + 4a:")
delta_b = {}
for s in sector_names:
    db = [B[s][k] + 4*A[s][k] for k in range(3)]
    delta_b[s] = db
    fk = FLAT_COORD[s]
    coords = [f"{'['+str(db[k])+']' if k == fk else str(db[k])}" for k in range(3)]
    print(f"    {s}: ({', '.join(coords)})  [brackets = flat coord]")

# Check: is δb_flat = -a_flat?
print(f"\n  Flatness check: δb_flat = -a_flat?")
for s in sector_names:
    fk = FLAT_COORD[s]
    db = delta_b[s][fk]
    neg_a = -A[s][fk]
    print(f"    {s}: δb_{['p','q','r'][fk]} = {db}, -a = {neg_a} {'✓' if db == neg_a else '✗'}")

# Form δB matrix (3×3)
delta_B_mat = [[delta_b[s][k] for k in range(3)] for s in sector_names]
print(f"\n  δB matrix (rows = sectors, cols = p,q,r):")
for i, s in enumerate(sector_names):
    print(f"    {s}: [{delta_B_mat[i][0]:>6}, {delta_B_mat[i][1]:>6}, {delta_B_mat[i][2]:>6}]")

# δB = B + 4A.  Is there a simple structure?
dBm = Bm + 4 * Am
print(f"\n  rank(δB) = {np.linalg.matrix_rank(dBm)}")
print(f"  det(δB) = {det3_frac(delta_B_mat)} = {float(det3_frac(delta_B_mat)):.6f}")
print(f"  Tr(δB) = {sum(delta_B_mat[i][i] for i in range(3))} = {float(sum(delta_B_mat[i][i] for i in range(3))):.6f}")

# δB columns
for k in range(3):
    col = [delta_B_mat[i][k] for i in range(3)]
    print(f"    δB[:,{['p','q','r'][k]}] = ({col[0]}, {col[1]}, {col[2]})")

# Frobenius norm of δB = "total deviation from kinetic optimum"
frob_dB = sum(delta_B_mat[i][k]**2 for i in range(3) for k in range(3))
print(f"\n  ||δB||² = {frob_dB} = {float(frob_dB):.6f}  (total deviation from kinetic optimum)")
print(f"  ||B||²  = {frobB} = {float(frobB):.6f}")
print(f"  ||δB||²/||B||² = {float(frob_dB/frobB):.6f}")

# ════════════════════════════════════════════════════════════════
# §14. CONSTRAINED S_K OPTIMIZATION: DOES PHYSICAL = MINIMUM?
# ════════════════════════════════════════════════════════════════

print_header(14, "CONSTRAINED KINETIC ACTION OPTIMIZATION")

print("""  S_K = Σ_s Σ_k f(a_k, b_k) where f(a,b) = (52/3)a² + 8ab + b²
  
  Since per-coordinate terms are independent, with fixed a:
    - Non-flat coords: min_b f(a,b) → b = -4a (kinetic optimum)
    - Flat coords: b = -5a (constrained → not at optimum)
  
  So the constrained-S_K minimum for fixed a gives:
    b_nonflat = -4 a_nonflat
    b_flat = -5 a_flat
  
  The physical non-flat b values are NOT all -4a:
""")

# Compare physical vs constrained-optimum
total_SK_physical = F(0)
total_SK_optimal = F(0)

for s in sector_names:
    fk = FLAT_COORD[s]
    for k in range(3):
        a, b = A[s][k], B[s][k]
        f_phys = F(52,3)*a**2 + 8*a*b + b**2
        total_SK_physical += f_phys
        # Optimal b
        b_opt = -5*a if k == fk else -4*a
        f_opt = F(52,3)*a**2 + 8*a*b_opt + b_opt**2
        total_SK_optimal += f_opt
        if k != fk:
            dev = b / a if a != 0 else float('inf')
            print(f"    {s}/{['p','q','r'][k]}: b/a = {dev}, physical f = {f_phys}, optimal f = {F(4,3)*a**2}, excess = {f_phys - F(4,3)*a**2}")

print(f"\n  S_K (physical) = {total_SK_physical} = {float(total_SK_physical):.6f}")
print(f"  S_K (constrained optimum) = {total_SK_optimal} = {float(total_SK_optimal):.6f}")
print(f"  Excess = S_K(phys) - S_K(opt) = {total_SK_physical - total_SK_optimal} = {float(total_SK_physical - total_SK_optimal):.6f}")
print(f"  → Physical point is {'NOT ' if total_SK_physical != total_SK_optimal else ''}the constrained S_K minimum")

# ════════════════════════════════════════════════════════════════
# §15. GENERATING-FUNCTION VARIATIONAL PRINCIPLE
# ════════════════════════════════════════════════════════════════

print_header(15, "GENERATING-FUNCTION VARIATIONAL PRINCIPLE")

print("""  The generating function (from parameter_reduction.py) gives:
    a_k(s) = α_k^a · (2I₃) + β_k^a · Nc + γ_k^a
    b_k(s) = α_k^b · (2I₃) + β_k^b · Nc + γ_k^b
    c_k(s) = α_k^c · (2I₃) + β_k^c · Nc + γ_k^c
  
  After constraints, 18 free parameters remain.
  
  Express S_K entirely in terms of the generating-function params:
    S_K = Σ_s Σ_k [(52/3)(α_k^a·X + β_k^a·N + γ_k^a)²
                    + 8(α_k^a·X + ...)(α_k^b·X + ...)
                    + (α_k^b·X + β_k^b·N + γ_k^b)²]
  where X = 2I₃(s), N = Nc(s).
  
  This is a QUADRATIC function of the 18 parameters.
  Look for critical points.
""")

# Quantum numbers for sectors
QN = {
    "lepton": (F(-1), F(1)),   # (2I₃, Nc)
    "up":     (F(1),  F(3)),
    "down":   (F(-1), F(3)),
}

# Generating function parameters from parameter_reduction.py (uniform (2I₃, Nc) basis)
GEN_PARAMS = {
    # (α, β, γ) for each coefficient, in (2I₃, Nc) basis
    "a_p": (F(11,8), F(-1), F(27,8)),
    "a_q": (F(1,4), F(-1,4), F(1)),       # α = -β
    "a_r": (F(-5,4), F(5,4), F(-4)),       # α = -β
    "b_p": (F(-33,8), F(17,4), F(-95,8)),
    "b_q": (F(-7,4), F(13,4), F(-21,2)),
    "b_r": (F(19,4), F(-17,4), F(33,2)),
    "c_p": (F(9,4), F(-7,2), F(33,4)),
    "c_q": (F(5,2), F(-7), F(29,2)),
    "c_r": (F(-3), F(2), F(-11)),
}

# Verify generating function reproduces physical coefficients
print(f"  Verification of generating function:")
labels = ["a_p", "a_q", "a_r", "b_p", "b_q", "b_r", "c_p", "c_q", "c_r"]
coeff_maps = {
    "a_p": lambda s: A[s][0], "a_q": lambda s: A[s][1], "a_r": lambda s: A[s][2],
    "b_p": lambda s: B[s][0], "b_q": lambda s: B[s][1], "b_r": lambda s: B[s][2],
    "c_p": lambda s: CC[s][0], "c_q": lambda s: CC[s][1], "c_r": lambda s: CC[s][2],
}
all_ok = True
for lbl in labels:
    alpha, beta, gamma = GEN_PARAMS[lbl]
    for s in sector_names:
        X, N_ = QN[s]
        pred = alpha * X + beta * N_ + gamma
        actual = coeff_maps[lbl](s)
        if pred != actual:
            print(f"    MISMATCH: {lbl}({s}) = {pred}, expected {actual}")
            all_ok = False
print(f"    All 27 coefficients verified: {'✓' if all_ok else '✗'}")

# Now express S_K as function of the 18 generating-function parameters.
# S_K = Σ_s Σ_k f(a_k(s), b_k(s))
# where f(a,b) = (52/3)a² + 8ab + b²
# and a_k(s) = α_k^a · X_s + β_k^a · N_s + γ_k^a, similarly for b.

# Order the 27 params as a vector: [α_a, β_a, γ_a for p,q,r; α_b, ...; α_c, ...]
# But only a and b appear in S_K (27-9=18 relevant params).
# After constraints: 18 params total. Let's use all 18 a+b params and check gradient.

# The 18 a+b params:  α,β,γ for a_p, a_q, a_r, b_p, b_q, b_r
param_names = []
param_values = []
for coord in ["p", "q", "r"]:
    for ab in ["a", "b"]:
        key = f"{ab}_{coord}"
        alpha, beta, gamma = GEN_PARAMS[key]
        param_names.extend([f"α_{key}", f"β_{key}", f"γ_{key}"])
        param_values.extend([float(alpha), float(beta), float(gamma)])

param_values = np.array(param_values)

def SK_from_gen_params(params):
    """Compute S_K from the 18 generating-function parameters (a,b only)."""
    # Extract params
    gen = {}
    idx = 0
    for coord in ["p", "q", "r"]:
        for ab in ["a", "b"]:
            key = f"{ab}_{coord}"
            gen[key] = (params[idx], params[idx+1], params[idx+2])
            idx += 3
    
    total = 0.0
    for s in sector_names:
        X_s = float(QN[s][0])
        N_s = float(QN[s][1])
        for k, coord in enumerate(["p", "q", "r"]):
            a_alpha, a_beta, a_gamma = gen[f"a_{coord}"]
            b_alpha, b_beta, b_gamma = gen[f"b_{coord}"]
            
            a_val = a_alpha * X_s + a_beta * N_s + a_gamma
            b_val = b_alpha * X_s + b_beta * N_s + b_gamma
            
            total += (52/3) * a_val**2 + 8 * a_val * b_val + b_val**2
    return total

# Verify S_K at physical values
SK_gen = SK_from_gen_params(param_values)
print(f"\n  S_K from gen params: {SK_gen:.6f}")
print(f"  S_K direct:         {float(S_K_total):.6f}")
print(f"  Match: {'✓' if abs(SK_gen - float(S_K_total)) < 1e-10 else '✗'}")

# Numerical gradient of S_K w.r.t. the 18 gen params
eps = 1e-8
grad_gen = np.zeros(18)
for i in range(18):
    dp = np.zeros(18)
    dp[i] = eps
    grad_gen[i] = (SK_from_gen_params(param_values + dp) - SK_gen) / eps

print(f"\n  Gradient of S_K w.r.t. generating-function parameters:")
for i in range(18):
    if abs(grad_gen[i]) > 0.01:
        print(f"    ∂S_K/∂{param_names[i]:>8} = {grad_gen[i]:+10.4f}")

print(f"  |∇S_K| = {np.linalg.norm(grad_gen):.6f}")
print(f"  → Physical point is {'NOT ' if np.linalg.norm(grad_gen) > 0.01 else ''}a critical point of S_K")

# Now test: if we COULD freely optimize the b params (9 of 18),
# what would the optimal b params be?
# For each coord k: ∂f/∂b = 8a + 2b = 0 → b = -4a
# In gen-param space: b_k(s) = -4·a_k(s) for all s
# → α_b = -4α_a, β_b = -4β_a, γ_b = -4γ_a

b_opt_params = param_values.copy()
for k, coord in enumerate(["p", "q", "r"]):
    # a params are at index k*6, k*6+1, k*6+2
    # b params are at index k*6+3, k*6+4, k*6+5
    for j in range(3):
        b_opt_params[k*6 + 3 + j] = -4 * param_values[k*6 + j]

SK_opt = SK_from_gen_params(b_opt_params)
print(f"\n  S_K with b = -4a (fully optimal): {SK_opt:.6f}")
print(f"  S_K physical:                      {float(S_K_total):.6f}")
print(f"  Ratio S_K(phys)/S_K(opt) = {float(S_K_total)/SK_opt:.6f}")

# What constraints prevent b = -4a?
# Flatness: for the flat coord k(s), b_flat = -5a (not -4a)
# This means α_b^flat ≠ -4α_a^flat, etc. in general.
# But in generating function, flatness is per-sector, not per-parameter.
# The gen-function gives b_k(s) as a LINEAR function of QN.
# We need b_{flat(s)}(s) = -5 a_{flat(s)}(s) for each specific (s, flat(s)).
# This gives 3 equations in the 9 b-generating params.

print(f"\n  Constraints on b-generating params from flatness:")
print(f"    (These are linear constraints mixing α,β,γ of different b_k)")
for s in sector_names:
    fk = FLAT_COORD[s]
    coord = ["p", "q", "r"][fk]
    X, N_ = QN[s]
    # b_{coord}(s) = -5 a_{coord}(s)
    # α_b X + β_b N + γ_b = -5(α_a X + β_a N + γ_a)
    a_alpha, a_beta, a_gamma = GEN_PARAMS[f"a_{coord}"]
    b_alpha, b_beta, b_gamma = GEN_PARAMS[f"b_{coord}"]
    lhs = b_alpha * X + b_beta * N_ + b_gamma
    rhs = -5 * (a_alpha * X + a_beta * N_ + a_gamma)
    print(f"    {s}/{coord}: b_{coord}({s}) = {lhs}, -5·a_{coord}({s}) = {rhs} {'✓' if lhs == rhs else '✗'}")
    # Deviation from -4a:
    opt = -4 * (a_alpha * X + a_beta * N_ + a_gamma)
    print(f"           kinetic-opt: {opt}, deviation from opt: {lhs - opt}")

# ════════════════════════════════════════════════════════════════
# §16. LAGRANGIAN-CONSTRAINED S_K MINIMUM (FULL COMPUTATION)
# ════════════════════════════════════════════════════════════════

print_header(16, "LAGRANGIAN-CONSTRAINED S_K MINIMUM")

print("""  Minimize S_K(A,B params) subject to:
    C1: flatness: b_{flat(s)}(s) = -5·a_{flat(s)}(s) for 3 (s, flat(s)) pairs
    C2: electron origin: (a + b + c)_ℓ = 0 → c_ℓ = -(a+b)_ℓ (doesn't affect S_K)
    C3: q-unification: 3(b_q^u-b_q^ℓ) + (c_q^u-c_q^ℓ) = 0 (involves c → doesn't affect S_K)
  
  Only C1 constrains (a,b).  Since a and b enter S_K quadratically,
  and the flatness constraints are LINEAR in (a,b), this is a standard
  quadratic program.
  
  For fixed a: S_K is quadratic in b with per-coordinate separation.
  min_{b} S_K s.t. b_{flat} = -5a_{flat}:
    → b_nonflat = -4a  (unconstrained)
    → b_flat = -5a     (constrained)
""")

# Compute the ACTUAL constrained optimum
print(f"  Comparison: physical vs constrained-S_K-optimal b:")
print(f"  {'Sector/coord':>14} {'a':>8} {'b(phys)':>10} {'b(opt)':>10} {'b(phys)/a':>10} {'b(opt)/a':>10} {'δb':>8}")
total_excess = F(0)
for s in sector_names:
    fk = FLAT_COORD[s]
    for k in range(3):
        a_val = A[s][k]
        b_phys = B[s][k]
        b_opt = (-5 if k == fk else -4) * a_val
        delta = b_phys - b_opt
        f_phys = F(52,3)*a_val**2 + 8*a_val*b_phys + b_phys**2
        f_opt = F(52,3)*a_val**2 + 8*a_val*b_opt + b_opt**2
        excess = f_phys - f_opt
        total_excess += excess
        label = f"{s}/{['p','q','r'][k]}{'*' if k == fk else ' '}"
        ratio_p = f"{b_phys/a_val}" if a_val != 0 else "∞"
        ratio_o = "-5" if k == fk else "-4"
        print(f"  {label:>14} {str(a_val):>8} {str(b_phys):>10} {str(b_opt):>10} {ratio_p:>10} {ratio_o:>10} {str(delta):>8}")

print(f"\n  Total S_K excess above constrained minimum: {total_excess}")
print(f"  = {float(total_excess):.6f}")

# Break down excess by non-flat coords
print(f"\n  Excess per non-flat coordinate (= δb² = (b - b_opt)² contribution to S_K):")
print(f"  Note: for fixed a, f(a,b) - f(a, -4a) = (b + 4a)² since")
print(f"        f(a,b) = (52/3)a² + 8ab + b² and f(a,-4a) = (4/3)a²")
print(f"        so f - f_opt = 8a(b+4a) + (b² - 16a²) = (b+4a)(8a+b+4a)")
print(f"             hmm, let's compute directly:")
for s in sector_names:
    fk = FLAT_COORD[s]
    for k in range(3):
        if k == fk:
            continue
        a_val = A[s][k]
        b_phys = B[s][k]
        db = b_phys + 4*a_val
        excess = F(52,3)*a_val**2 + 8*a_val*b_phys + b_phys**2 - F(4,3)*a_val**2
        print(f"    {s}/{['p','q','r'][k]}: a={a_val}, δb = b+4a = {db}, excess = {excess} (= δb(δb+8a) hmm)")
        # Actually f(a,b) - f(a,-4a) = 8a(b-(-4a)) + (b²-16a²) = 8a·δb + (b-4a)(b+4a)
        # = 8a·δb + (δb-8a)·δb = δb(8a + δb - 8a) = δb² ... no
        # f(a,b) = (52/3)a² + 8ab + b²
        # f(a, -4a) = (52/3)a² - 32a² + 16a² = (52/3 - 16)a² = (52-48)/3 · a² = 4/3 a²
        # f - f_opt = 8a(b+4a) + (b²-16a²) = 8aδb + (b-4a)(b+4a) = 8aδb + δb(b-4a)
        # hmm, b-4a = (b+4a) - 8a = δb - 8a.
        # = 8aδb + δb(δb - 8a) = 8aδb + δb² - 8aδb = δb²
        check = db**2
        print(f"           excess = δb² = ({db})² = {check}  {'✓' if excess == check else '✗'}")

# ════════════════════════════════════════════════════════════════
# §17. SUMMARY AND INTERPRETATION
# ════════════════════════════════════════════════════════════════

print_header(17, "SUMMARY AND INTERPRETATION")

print(f"""  1. PARABOLIC TRAJECTORIES: each sector traces x(g) = a·g² + b·g + c
     in (p,q,r)-space with g = 1,2,3 (generations).

  2. KINETIC ACTION: S_K = Σ_s ∫₁³ ½|ẋ|² dg is quadratic in (a, b),
     independent of c.  Total physical S_K = {S_K_total} = {float(S_K_total):.4f}.

  3. KINETIC OPTIMUM: with fixed a, minimizing S_K gives b* = -4a,
     i.e., turning point at generation g* = 2.  This is the "least action"
     for the velocity of generational transitions.

  4. FLATNESS PENALTY: the physicaly imposed flatness b_flat = -5a
     pushes g* from 2 → 5/2, costing ΔS = a² per flat coordinate.
     Total flatness cost: {total_flatness_cost} = {float(total_flatness_cost):.4f}.

  5. EXCESS ACTION: the physical point is NOT the constrained-S_K minimum.
     Total excess = {total_excess} = {float(total_excess):.4f}.
     The excess equals Σ (δb_nonflat)² where δb = b + 4a.

  6. MATRIX FACTORIZATION: B = A · Mᵀ and C = A · Nᵀ decompose all 27
     coefficients into a curvature matrix A multiplied by transfer matrices.
     det(M) = {detM}, Tr(M) = {trM}.
     det(N) = {detN}, Tr(N) = {trN}.

  7. PARAMETER HIERARCHY:
     a) Curvature A is fixed by quantum numbers (2I₃, Nc): 7 params
     b) Flat-coord velocities are fixed by flatness: b_flat = -5a
     c) Optimal non-flat velocities would be: b = -4a (least kinetic action)
     d) Physical non-flat b DEVIATE from the optimum by δb ≠ 0
     e) THE δb's ARE THE 6 UNEXPLAINED PARAMETERS
     f) c is then fixed by electron origin (3 params) + q-unification (1 param)
     g) Remaining 5 c-params are also free

  8. NEXT QUESTIONS:
     a) What principle determines the 6 non-flat δb values?
     b) Do the 6 δb values minimize a potential V(x) or bridge condition?
     c) Is the CKM/PMNS mixing encoded in the δb structure?
     d) Can we unify δb with a single additional constraint?
""")


# Close tee
sys.stdout = _orig_stdout
_outfile.close()
print(f"Output saved to: {OUTPUT_PATH}")
