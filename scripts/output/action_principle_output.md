
========================================================================
  §1. KINETIC ACTION
========================================================================

  For x_k(g) = a_k g² + b_k g + c_k, the velocity is ẋ_k(g) = 2a_k g + b_k.
  
  S_K = Σ_s Σ_k ∫₁³ ½(2a_k g + b_k)² dg
      = Σ_s Σ_k ½[4a_k² ∫g² + 4a_k b_k ∫g + b_k² ∫1] dg
  
  ∫₁³ g² dg = 26/3,  ∫₁³ g dg = 4,  ∫₁³ 1 dg = 2
  
  So: S_K = Σ_s Σ_k [ (52/3) a_k² + 8 a_k b_k + b_k² ]

  Per-sector kinetic action:
    S_K(lepton) = 233/12 = 19.416667
    S_K(up) = 1343/48 = 27.979167
    S_K(down) = 41/3 = 13.666667
    S_K(total) = 977/16 = 61.062500

========================================================================
  §2. CURVATURE, VELOCITY-NORM, AND OFFSET ACTIONS
========================================================================

  S_A  = Σ |a|²  = 177/16 = 11.062500
  S_B  = Σ |b|²  = 3437/16 = 214.812500
  S_C  = Σ |c|²  = 505/2 = 252.500000
  S_AB = Σ a·b   = -691/16 = -43.187500
  S_AC = Σ a·c   = 53/2 = 26.500000
  S_BC = Σ b·c   = -675/4 = -168.750000

  S_AB / S_A = -3.903955  (would be -5 if b = -5a universally)

========================================================================
  §3. MATRIX INVARIANTS
========================================================================

  Form 3×3 matrices from coefficient vectors:
    A_{sk} = a_k(sector s),  B_{sk} = b_k(s),  C_{sk} = c_k(s)
  where rows = sectors (ℓ, u, d), cols = coordinates (p, q, r).

  A (curvature):
    det = -0.375000
    tr  = 2.500000
    eigenvalues = 2.5889+0.0000i, -0.4276+0.0000i, 0.3387+0.0000i
    singular values = 3.2882, 0.4202, 0.2714
    Frobenius norm = 3.326034

  B (velocity):
    det = -125.500000
    tr  = -7.000000
    eigenvalues = -11.9767+0.0000i, 2.4883+2.0705i, 2.4883-2.0705i
    singular values = 13.8918, 4.1279, 2.1885
    Frobenius norm = 14.656483

  C (offset):
    det = 128.000000
    tr  = -3.500000
    eigenvalues = 9.3306+0.0000i, -1.1772+0.0000i, -11.6535+0.0000i
    singular values = 12.2906, 10.0180, 1.0396
    Frobenius norm = 15.890249

  A + B + C = x(g=1):
    lepton: [0. 0. 0.]
    up: [-1.5 -6.  -1. ]
    down: [-0.5 -8.  -2. ]

  4A + 2B + C = x(g=2):
    lepton: [-0.5 -4.   3. ]
    up: [ 0.5 -7.   3. ]
    down: [ 1.5 -7.   0. ]

  9A + 3B + C = x(g=3):
    lepton: [ 1. -7.  3.]
    up: [ 6. -7.  4.]
    down: [ 1.5 -6.   4. ]

========================================================================
  §4. FLATNESS-CONSTRAINED LAGRANGIAN
========================================================================

  The flatness condition b_{k(s)} = -5 a_{k(s)} means the velocity
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

  Force vectors F = 2a (constant per sector):
    F_lepton = (2, 1, -3)
    F_up = (7/2, 1, -3)
    F_down = (-2, 0, 2)

========================================================================
  §5. SECOND-LEVEL ACTION: WHAT DETERMINES a(s)?
========================================================================

  The curvature matrix A_{sk} = a_k(s) is a 3×3 matrix.
  We need a SECOND-LEVEL action Ŝ[A] whose variation δŜ/δA = 0
  yields the observed curvature.  Similarly for B, C (or B, C derived from A).
  
  Candidate functionals:
    Ŝ₁ = Tr(AᵀA) = Frobenius norm² = Σ a_k²
    Ŝ₂ = det(A)
    Ŝ₃ = Tr(A) + λ det(A)  (characteristic polynomial)
    Ŝ₄ = Tr(AᵀA) + μ Tr(A²)  (with symmetry)
    Ŝ₅ = Σ_s (Σ_k a_k)²  (trace per sector)
    Ŝ₆ = Σ_s a_s · diag(w) · a_s  (weighted norm with metric)

  Tr(A) = a_p^ℓ + a_q^u + a_r^d = 1 + 1/2 + 1 = 5/2
  Flat-diagonal: a_r^ℓ + a_q^u + a_p^d = -3/2 + 1/2 + -1 = -2
  ← This is the sum of curvatures in FLAT coordinates!

  det(A) = -3/8 = -0.375000
  det(B) = -251/2 = -125.500000
  det(C) = 128 = 128.000000

  Tr(A²) = 7 = 7.000000
  Tr(B²) = 589/4 = 147.250000

  ||A||_F² = Tr(AᵀA) = 177/16 = 11.062500
  ||B||_F² = Tr(BᵀB) = 3437/16 = 214.812500

  Characteristic polynomial of A:
    λ³ - (5/2)λ² + (-3/8)λ - (-3/8) = 0

========================================================================
  §6. MINIMUM-NORM TEST
========================================================================

  Question: among all coefficient sets satisfying the 7 constraints
  (electron origin ×3, flatness ×3, q-unification ×1), do the physical
  values minimize the total Frobenius norm ||A||² + ||B||² + ||C||²?
  
  Strategy: parameterize the constraint manifold and check if the
  gradient of the norm vanishes at the physical point.

  Constraint residual at physical point: 0.00e+00

  ||θ||² at physical point: 478.375000
  |projected gradient|: 43.743571
  → Physical point is NOT a minimum of ||θ||²

  Searching for weights (1, λ, μ) such that ||A||² + λ||B||² + μ||C||²
  has a critical point at the physical values...

  Per-sector kinetic action breakdown:
    lepton: S_flat = 5.2500, S_nonflat = 14.1667, total = 19.4167
           S_flat check (7/3·a²): 5.2500 ✓
    up: S_flat = 0.5833, S_nonflat = 27.3958, total = 27.9792
           S_flat check (7/3·a²): 0.5833 ✓
    down: S_flat = 2.3333, S_nonflat = 11.3333, total = 13.6667
           S_flat check (7/3·a²): 2.3333 ✓

========================================================================
  §7. GENERATIONAL EULER-LAGRANGE ANALYSIS
========================================================================

  For a general Lagrangian L(x, ẋ, g) = T(ẋ) - V(x, g), the E-L eq is:
    ẍ_k = -∂V/∂x_k
  
  Since ẍ_k = 2a_k = const, V must be linear in x: V = -2a·x.
  
  BUT: the force a = a(s) depends on the sector, not on g.
  So V = V(x, s) = -2 a(s) · x.
  
  The deeper question: is there a UNIVERSAL potential V(x) valid for all sectors?
  → Check if a(s) = -½ ∇V evaluated at the sector's trajectory.

  Gen-2 positions (midpoint of each sector trajectory):
    lepton: x₂ = (-0.50, -4.00, 3.00), |a| = 1.8708
    up: x₂ = (0.50, -7.00, 3.00), |a| = 2.3585
    down: x₂ = (1.50, -7.00, 0.00), |a| = 1.4142

  Force magnitude vs position norm:
    lepton g=1: |x| = 0.0000, |F| = 3.7417, |F|/|x| = 3741657386773941.0000
    lepton g=2: |x| = 5.0249, |F| = 3.7417, |F|/|x| = 0.7446
    lepton g=3: |x| = 7.6811, |F| = 3.7417, |F|/|x| = 0.4871
    up g=1: |x| = 6.2650, |F| = 4.7170, |F|/|x| = 0.7529
    up g=2: |x| = 7.6322, |F| = 4.7170, |F|/|x| = 0.6180
    up g=3: |x| = 10.0499, |F| = 4.7170, |F|/|x| = 0.4694
    down g=1: |x| = 8.2614, |F| = 2.8284, |F|/|x| = 0.3424
    down g=2: |x| = 7.1589, |F| = 2.8284, |F|/|x| = 0.3951
    down g=3: |x| = 7.3655, |F| = 2.8284, |F|/|x| = 0.3840

========================================================================
  §8. CROSS-MATRIX RELATIONS
========================================================================

  Can we express B as a function of A?
  Test: B = α·A + β·I + γ·P (where P is the flatness permutation matrix)

  det(A) = -0.375000

  B = A · Mᵀ  where M =
    [  0.3333,   8.3333,   5.3333]
    [  4.0000,  -4.0000,   5.0000]
    [  1.3333,  13.3333,   0.3333]
  Residual ||A·Mᵀ - B||_F = 1.72e-15

  M in exact fractions:
    Mᵀ (rows correspond to coords p,q,r):
      Mᵀ[p,:] = [1/3, 25/3, 16/3]
      Mᵀ[q,:] = [4, -4, 5]
      Mᵀ[r,:] = [4/3, 40/3, 1/3]

    M = (Mᵀ)ᵀ  (exact):
      M[p,:] = [1/3, 25/3, 16/3]
      M[q,:] = [4, -4, 5]
      M[r,:] = [4/3, 40/3, 1/3]

    Tr(M) = -10/3 = -3.333333
    det(M) = 1004/3 = 334.666667
    Consistency ||M_frac - M_numpy||: 3.13e-15

  C = A · Nᵀ  where:
    N = 
      N[p,:] = [-10/3, -71/6, -47/6]
      N[q,:] = [-12, -29, -21]
      N[r,:] = [-8/3, -62/3, -14/3]

    Tr(N) = -37 = -37.000000
    det(N) = -1024/3 = -341.333333
    Residual ||A·Nᵀ - C||_F = 9.47e-15
    Consistency ||N_frac - N_numpy||: 8.84e-15

========================================================================
  §9. THE VELOCITY-CURVATURE TRANSFER MATRIX M
========================================================================

  We found B = A · Mᵀ, i.e., the velocity matrix is determined by
  the curvature matrix via a universal 3×3 transfer matrix M.
  
  This means: b_k(s) = Σ_j M_{kj} a_j(s)
  
  So the velocity in coordinate k depends on curvatures in ALL coords j.
  The 9 free parameters in B are encoded in the 9 entries of M.
  
  Can we understand M's structure?

  Is M symmetric? False
  Antisymmetric part M - Mᵀ:
    [0, 13/3, 4]
    [-13/3, 0, -25/3]
    [-4, 25/3, 0]

  Diagonal of M: [1/3, -4, 1/3]
  ← The flat-coord diagonal entries should relate to -5 (flatness)

  Flatness check via M:
    lepton (r-flat): Σ_j M_rj · a_j = 7.500000, -5·a_r = 7.500000 ✓
    up (q-flat): Σ_j M_qj · a_j = -2.500000, -5·a_q = -2.500000 ✓
    down (p-flat): Σ_j M_pj · a_j = 5.000000, -5·a_p = 5.000000 ✓

  Eigenvalues of M: 10.363955+0.000000i, -3.026012+0.000000i, -10.671276+0.000000i
  Does M have eigenvalue -5? False

  Eigenvalues of N: -44.941980+0.000000i, 1.112007+0.000000i, 6.829973+0.000000i

  M · N:
    [-115.3333, -355.8333, -202.5000]
    [ 21.3333, -34.6667,  29.3333]
    [-165.3333, -409.3333, -292.0000]
  N · M:
    [-58.8889, -84.8889, -79.5556]
    [-148.0000, -264.0000, -216.0000]
    [-89.7778,  -1.7778, -119.1111]
  Commutator [M,N] = M·N - N·M:
    [-56.4444, -270.9444, -122.9444]
    [169.3333, 229.3333, 245.3333]
    [-75.5556, -407.5556, -172.8889]
  ||[M,N]|| = 659.445937

========================================================================
  §10. EFFECTIVE ACTION S[A, M, N]
========================================================================

  Given B = A·Mᵀ and C = A·Nᵀ, the 27 parameters reduce to:
    9 (curvature A) + 9 (transfer M) + 9 (transfer N) = 27
  but with constraints becoming conditions on M and N.
  
  The kinetic action becomes:
    S_K = Σ_s ∫₁³ ½|2A_s g + B_s|² dg
        = Σ_s [(52/3)|a_s|² + 8 a_s·(A·Mᵀ)_s + |A·Mᵀ_s|²]
  
  If we further express A in terms of quantum numbers:
    A_{sk} = α_k · (2I₃) + β_k · Nc + γ_k
  then the entire action is a function of (α, β, γ, M, N).

  I + Mᵀ + Nᵀ (should annihilate A_ℓ):
    [ -2.0000,  -8.0000,  -1.3333]
    [ -3.5000, -32.0000,  -7.3333]
    [ -2.5000, -16.0000,  -3.3333]

  A_ℓ · (I + Mᵀ + Nᵀ) = [0.000000, 0.000000, 0.000000]
  (should be zero: True)
  A_up · (I+Mᵀ+Nᵀ) = [-1.500000, -6.000000, -1.000000]  (= x_up(g=1))
  A_down · (I+Mᵀ+Nᵀ) = [-0.500000, -8.000000, -2.000000]  (= x_down(g=1))

========================================================================
  §11. VARIATIONAL SEARCH FOR M
========================================================================

  Given that A is fixed by quantum numbers, we ask:
  Is there a function f(M) such that ∂f/∂M_{ij} = 0 at the physical M?
  
  Test candidates:
    f₁ = Tr(M)          → trivially not (Tr doesn't couple entries)
    f₂ = det(M)          → ∂det/∂M = cofactor matrix
    f₃ = Tr(M²)          → ∂/∂M = 2M
    f₄ = Tr(M²) + λ det(M)  → combines
    f₅ = Tr(MᵀM)         → ∂/∂M = 2M
    f₆ = ||AM||² = Tr(AᵀA MᵀM) → couples A and M

  Tr(M) = -3.333333
  Tr(M²) = 230.444444
  det(M) = 334.666667

  M · Mᵀ:
    [ 98.0000,  -5.3333, 113.3333]
    [ -5.3333,  57.0000, -46.3333]
    [113.3333, -46.3333, 179.6667]
  M · Mᵀ proportional to I? False

  Testing: is M a saddle point of S_K[M] (kinetic action as function of M)?
  Holding A fixed, varying M:
  S_K(physical M) = 61.062500
  ∂S_K/∂M:
    [ 12.1250,   4.2500, -10.7500]
    [-10.7500,  -4.0000,  14.0000]
    [  5.7500,   4.0000,  -6.0000]
  |∇_M S_K| = 26.329938
  → M is NOT a critical point of S_K

  Flatness constraint rank: 3
  |projected ∇_M S_K| (after flatness constraints) = 6.097020

========================================================================
  §12. FULL ACTION LANDSCAPE
========================================================================

  The kinetic action S_K = Σ_s ∫₁³ ½|ẋ|² dg is a QUADRATIC function of (a, b).
  Since c does not appear (only ẋ = 2ag + b depends on a,b), S_K is independent of c.
  
  So S_K determines (a, b) = 18 of the 27 params, while c is determined by:
    - Electron origin: 3 eqs
    - Bridge/q-unification: further constraints
  
  Remaining: the action S_K(a,b) subject to flatness constraints.

  Unconstrained minimum of f(a,b) = (52/3)a² + 8ab + b²:
    ∂f/∂b = 8a + 2b = 0 → b* = -4a (turning point at g* = 2)
    ∂f/∂a = (104/3)a + 8b = 0 → b* = -(13/3)a (turning point at g* = 13/6)
    Joint optimum: a = b = 0 (trivial)

  With fixed a: b* = -4a minimizes S_K → g* = -b/(2a) = 2 (exactly gen-2!)
  Physical flat coords have b = -5a → g* = 5/2 (between gen 2 and 3)
  Physical non-flat coords have various b/a ratios:
    lepton/p: b/a =     -7/2, g* = 1.7500, deviation from -4: 1/2
    lepton/q: b/a =      -11, g* = 5.5000, deviation from -4: -7
    up/p: b/a =    -13/7, g* = 0.9286, deviation from -4: 15/7
    up/r: b/a =    -17/3, g* = 2.8333, deviation from -4: -5/3
    down/q: b/a =        ∞, g* = -inf, deviation from -4: –
    down/r: b/a =       -1, g* = 0.5000, deviation from -4: 3

  Key insight: if b = -4a minimized S_K (turning point at g=2),
  and flatness forces b = -5a (turning point at g=5/2),
  then the flat coordinate is PUSHED from the kinetic optimum.
  
  The "cost" of flatness: ΔS = f(a,-5a) - f(a,-4a) for each flat coord

    lepton: a_r = -3/2, ΔS = 9/4 = a²·(7/3 - 4/3) = a² = 9/4
    up: a_q = 1/2, ΔS = 1/4 = a²·(7/3 - 4/3) = a² = 1/4
    down: a_p = -1, ΔS = 1 = a²·(7/3 - 4/3) = a² = 1
    Total flatness cost: Σ a²_flat = 7/2 = 3.500000

========================================================================
  §13. DEVIATION FROM KINETIC OPTIMUM
========================================================================

  Key insight: with fixed a, the kinetic-optimal velocity is b* = -4a
  (turning point at g = 2).  Define the deviation:
    δb_k(s) = b_k(s) - (-4 a_k(s)) = b_k(s) + 4 a_k(s)
  
  For flat coords: δb_flat = -5a + 4a = -a (deviation = -a).
  For non-flat: δb is free.

  Deviation vectors δb = b + 4a:
    lepton: (1/2, -7/2, [3/2])  [brackets = flat coord]
    up: (15/4, [-1/2], 5/2)  [brackets = flat coord]
    down: ([1], 1, 3)  [brackets = flat coord]

  Flatness check: δb_flat = -a_flat?
    lepton: δb_r = 3/2, -a = 3/2 ✓
    up: δb_q = -1/2, -a = -1/2 ✓
    down: δb_p = 1, -a = 1 ✓

  δB matrix (rows = sectors, cols = p,q,r):
    lepton: [   1/2,   -7/2,    3/2]
    up: [  15/4,   -1/2,    5/2]
    down: [     1,      1,      3]

  rank(δB) = 3
  det(δB) = 35 = 35.000000
  Tr(δB) = 3 = 3.000000
    δB[:,p] = (1/2, 15/4, 1)
    δB[:,q] = (-7/2, -1/2, 1)
    δB[:,r] = (3/2, 5/2, 3)

  ||δB||² = 741/16 = 46.312500  (total deviation from kinetic optimum)
  ||B||²  = 3437/16 = 214.812500
  ||δB||²/||B||² = 0.215595

========================================================================
  §14. CONSTRAINED KINETIC ACTION OPTIMIZATION
========================================================================

  S_K = Σ_s Σ_k f(a_k, b_k) where f(a,b) = (52/3)a² + 8ab + b²
  
  Since per-coordinate terms are independent, with fixed a:
    - Non-flat coords: min_b f(a,b) → b = -4a (kinetic optimum)
    - Flat coords: b = -5a (constrained → not at optimum)
  
  So the constrained-S_K minimum for fixed a gives:
    b_nonflat = -4 a_nonflat
    b_flat = -5 a_flat
  
  The physical non-flat b values are NOT all -4a:

    lepton/p: b/a = -7/2, physical f = 19/12, optimal f = 4/3, excess = 1/4
    lepton/q: b/a = -11, physical f = 151/12, optimal f = 1/3, excess = 49/4
    up/p: b/a = -13/7, physical f = 871/48, optimal f = 49/12, excess = 225/16
    up/r: b/a = -17/3, physical f = 37/4, optimal f = 3, excess = 25/4
    down/q: b/a = inf, physical f = 1, optimal f = 0, excess = 1
    down/r: b/a = -1, physical f = 31/3, optimal f = 4/3, excess = 9

  S_K (physical) = 977/16 = 61.062500
  S_K (constrained optimum) = 73/4 = 18.250000
  Excess = S_K(phys) - S_K(opt) = 685/16 = 42.812500
  → Physical point is NOT the constrained S_K minimum

========================================================================
  §15. GENERATING-FUNCTION VARIATIONAL PRINCIPLE
========================================================================

  The generating function (from parameter_reduction.py) gives:
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

  Verification of generating function:
    All 27 coefficients verified: ✓

  S_K from gen params: 61.062500
  S_K direct:         61.062500
  Match: ✓

  Gradient of S_K w.r.t. generating-function parameters:
    ∂S_K/∂   α_a_p =   +22.6667
    ∂S_K/∂   β_a_p =  +126.6667
    ∂S_K/∂   γ_a_p =   +46.6667
    ∂S_K/∂   α_b_p =    +4.5000
    ∂S_K/∂   β_b_p =   +29.5000
    ∂S_K/∂   γ_b_p =   +10.5000
    ∂S_K/∂   α_a_q =   +16.0000
    ∂S_K/∂   β_a_q =   -10.6667
    ∂S_K/∂   γ_a_q =   -21.3333
    ∂S_K/∂   α_b_q =    +4.0000
    ∂S_K/∂   β_b_q =    -4.0000
    ∂S_K/∂   γ_b_q =    -6.0000
    ∂S_K/∂   α_a_r =   -18.6667
    ∂S_K/∂   β_a_r =  +136.0000
    ∂S_K/∂   γ_a_r =   +50.6667
    ∂S_K/∂   α_b_r =    -4.0000
    ∂S_K/∂   β_b_r =   +36.0000
    ∂S_K/∂   γ_b_r =   +14.0000
  |∇S_K| = 208.684441
  → Physical point is NOT a critical point of S_K

  S_K with b = -4a (fully optimal): 14.750000
  S_K physical:                      61.062500
  Ratio S_K(phys)/S_K(opt) = 4.139831

  Constraints on b-generating params from flatness:
    (These are linear constraints mixing α,β,γ of different b_k)
    lepton/r: b_r(lepton) = 15/2, -5·a_r(lepton) = 15/2 ✓
           kinetic-opt: 6, deviation from opt: 3/2
    up/q: b_q(up) = -5/2, -5·a_q(up) = -5/2 ✓
           kinetic-opt: -2, deviation from opt: -1/2
    down/p: b_p(down) = 5, -5·a_p(down) = 5 ✓
           kinetic-opt: 4, deviation from opt: 1

========================================================================
  §16. LAGRANGIAN-CONSTRAINED S_K MINIMUM
========================================================================

  Minimize S_K(A,B params) subject to:
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

  Comparison: physical vs constrained-S_K-optimal b:
    Sector/coord        a    b(phys)     b(opt)  b(phys)/a   b(opt)/a       δb
       lepton/p         1       -7/2         -4       -7/2         -4      1/2
       lepton/q       1/2      -11/2         -2        -11         -4     -7/2
       lepton/r*     -3/2       15/2       15/2         -5         -5        0
           up/p       7/4      -13/4         -7      -13/7         -4     15/4
           up/q*      1/2       -5/2       -5/2         -5         -5        0
           up/r      -3/2       17/2          6      -17/3         -4      5/2
         down/p*       -1          5          5         -5         -5        0
         down/q         0          1          0          ∞         -4        1
         down/r         1         -1         -4         -1         -4        3

  Total S_K excess above constrained minimum: 685/16
  = 42.812500

  Excess per non-flat coordinate (= δb² = (b - b_opt)² contribution to S_K):
  Note: for fixed a, f(a,b) - f(a, -4a) = (b + 4a)² since
        f(a,b) = (52/3)a² + 8ab + b² and f(a,-4a) = (4/3)a²
        so f - f_opt = 8a(b+4a) + (b² - 16a²) = (b+4a)(8a+b+4a)
             hmm, let's compute directly:
    lepton/p: a=1, δb = b+4a = 1/2, excess = 1/4 (= δb(δb+8a) hmm)
           excess = δb² = (1/2)² = 1/4  ✓
    lepton/q: a=1/2, δb = b+4a = -7/2, excess = 49/4 (= δb(δb+8a) hmm)
           excess = δb² = (-7/2)² = 49/4  ✓
    up/p: a=7/4, δb = b+4a = 15/4, excess = 225/16 (= δb(δb+8a) hmm)
           excess = δb² = (15/4)² = 225/16  ✓
    up/r: a=-3/2, δb = b+4a = 5/2, excess = 25/4 (= δb(δb+8a) hmm)
           excess = δb² = (5/2)² = 25/4  ✓
    down/q: a=0, δb = b+4a = 1, excess = 1 (= δb(δb+8a) hmm)
           excess = δb² = (1)² = 1  ✓
    down/r: a=1, δb = b+4a = 3, excess = 9 (= δb(δb+8a) hmm)
           excess = δb² = (3)² = 9  ✓

========================================================================
  §17. SUMMARY AND INTERPRETATION
========================================================================

  1. PARABOLIC TRAJECTORIES: each sector traces x(g) = a·g² + b·g + c
     in (p,q,r)-space with g = 1,2,3 (generations).

  2. KINETIC ACTION: S_K = Σ_s ∫₁³ ½|ẋ|² dg is quadratic in (a, b),
     independent of c.  Total physical S_K = 977/16 = 61.0625.

  3. KINETIC OPTIMUM: with fixed a, minimizing S_K gives b* = -4a,
     i.e., turning point at generation g* = 2.  This is the "least action"
     for the velocity of generational transitions.

  4. FLATNESS PENALTY: the physicaly imposed flatness b_flat = -5a
     pushes g* from 2 → 5/2, costing ΔS = a² per flat coordinate.
     Total flatness cost: 7/2 = 3.5000.

  5. EXCESS ACTION: the physical point is NOT the constrained-S_K minimum.
     Total excess = 685/16 = 42.8125.
     The excess equals Σ (δb_nonflat)² where δb = b + 4a.

  6. MATRIX FACTORIZATION: B = A · Mᵀ and C = A · Nᵀ decompose all 27
     coefficients into a curvature matrix A multiplied by transfer matrices.
     det(M) = 1004/3, Tr(M) = -10/3.
     det(N) = -1024/3, Tr(N) = -37.

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

