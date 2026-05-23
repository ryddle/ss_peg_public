========================================================================
  §1. CURVATURE TENSOR DERIVATION
========================================================================

  Measured parabolic coefficients (a, b, c) per sector per coordinate:
  Sector   Coord         a        b        c
  ──────── ────── ──────── ──────── ────────
  lepton   p       +1.0000  -3.5000  +2.5000
  lepton   q       +0.5000  -5.5000  +5.0000
  lepton   r       -1.5000  +7.5000  -6.0000

  up       p       +1.7500  -3.2500  +0.0000
  up       q       +0.5000  -2.5000  -4.0000
  up       r       -1.5000  +8.5000  -8.0000

  down     p       -1.0000  +5.0000  -4.5000
  down     q       +0.0000  +1.0000  -9.0000
  down     r       +1.0000  -1.0000  -2.0000

  Curvature vectors:
    a_lepton: [ 1.   0.5 -1.5]
    a_up: [ 1.75  0.5  -1.5 ]
    a_down: [-1.00000000e+00  2.22044605e-16  1.00000000e+00]

  Solving a_p = α·(2I₃) + β·Nc + γ:
    System matrix:
      lepton: [  -1.0,   +1,   +1] · (α, β, γ) = +1.0000
      up: [  +1.0,   +3,   +1] · (α, β, γ) = +1.7500
      down: [  -1.0,   +3,   +1] · (α, β, γ) = -1.0000

    Solution: α = 1.375000, β = -1.000000, γ = 3.375000
    → a_p = 1.3750·(2I₃) + -1.0000·Nc + 3.3750
    Check lepton: pred=+1.0000, actual=+1.0000, diff=0.00e+00
    Check up: pred=+1.7500, actual=+1.7500, diff=2.22e-16
    Check down: pred=-1.0000, actual=-1.0000, diff=0.00e+00

  Solving a_q = α·(2I₃) + β·Nc + γ:
    Solution: α = 0.250000, β = -0.250000, γ = 1.000000
    → a_q = 0.2500·(2I₃) + -0.2500·Nc + 1.0000

  Solving a_r = α·(2I₃) + β·Nc + γ:
    Solution: α = -1.250000, β = 1.250000, γ = -4.000000
    → a_r = -1.2500·(2I₃) + 1.2500·Nc + -4.0000

========================================================================
  §2. VELOCITY FROM FLATNESS CONSTRAINT
========================================================================

  Flatness principle: v_{k(s)}(2.5) = 5·a_{k(s)} + b_{k(s)} = 0
  → b_{k(s)} = -5·a_{k(s)} for the flat coordinate k of sector s

  lepton (r-flat): b_r = +7.5000, predicted = +7.5000  ✓
  up (q-flat): b_q = -2.5000, predicted = -2.5000  ✓
  down (p-flat): b_p = +5.0000, predicted = +5.0000  ✓

  Free velocity components (NOT constrained by flatness):
    lepton: b_p, b_q are free
      b_p = -3.5000
      b_q = -5.5000
    up: b_p, b_r are free
      b_p = -3.2500
      b_r = +8.5000
    down: b_q, b_r are free
      b_q = +1.0000
      b_r = -1.0000

  Solving b_p = α·(2I₃) + β·Nc + γ:
    Solution: α = -4.125000, β = 4.250000, γ = -11.875000

  Solving b_q = α·(2I₃) + β·Nc + γ:
    Solution: α = -1.750000, β = 3.250000, γ = -10.500000

  Solving b_r = α·(2I₃) + β·Nc + γ:
    Solution: α = 4.750000, β = -4.250000, γ = 16.500000

========================================================================
  §3. OFFSET DERIVATION (ELECTRON ORIGIN + BRIDGE)
========================================================================
  Measured offset vectors:
    c_lepton: [ 2.5  5.  -6. ]
    c_up: [ 0. -4. -8.]
    c_down: [-4.5 -9.  -2. ]

  Electron origin constraint:
    a_ℓ + b_ℓ + c_ℓ = (+0.0000, +0.0000, -0.0000)
    → c_ℓ = -a_ℓ - b_ℓ = (+2.5000, +5.0000, -6.0000)
    Measured c_ℓ: [ 2.5  5.  -6. ]
    Consistent: True

  q-unification constraint: q_ℓ(3) = q_u(3)
    q_ℓ(3) = -7.0
    q_u(3) = -7.0
    Match: True

  BVP q-unification constraint: 3·Δb_q + Δc_q = 0
    Δb_q = b_q^(u) - b_q^(ℓ) = +3.0000
    Δc_q = c_q^(u) - c_q^(ℓ) = -9.0000
    3·Δb_q + Δc_q = +0.0000

  b/τ bridge constraint:
    τ coords: [ 1. -7.  3.]
    b coords: [ 1.5 -6.   4. ]
    Bridge vector (b - τ): (+0.5, +1, +1)

  Solving c_p = α·(2I₃) + β·Nc + γ:
    Solution: α = 2.250000, β = -3.500000, γ = 8.250000

========================================================================
  §4. FULL VERIFICATION — 36 FERMION COORDINATES
========================================================================

  Reproducing all 36 coordinates (9 fermions × 3 coords + 3 bosons × 3 coords)
  from the parabolic model x_s(g) = a·g² + b·g + c

  Particle   g   k    Model   Actual          Δ
  ──────── ─── ─── ──────── ──────── ──────────
  e          1   p    +0.00    +0.00   4.44e-16 
  e          1   q    +0.00    +0.00   8.88e-16 
  e          1   r    -0.00    +0.00   8.88e-16 
  mu         2   p    -0.50    -0.50   4.44e-16 
  mu         2   q    -4.00    -4.00   0.00e+00 
  mu         2   r    +3.00    +3.00   0.00e+00 
  tau        3   p    +1.00    +1.00   4.44e-16 
  tau        3   q    -7.00    -7.00   1.78e-15 
  tau        3   r    +3.00    +3.00   1.78e-15 
  u          1   p    -1.50    -1.50   2.22e-16 
  u          1   q    -6.00    -6.00   0.00e+00 
  u          1   r    -1.00    -1.00   1.78e-15 
  c          2   p    +0.50    +0.50   0.00e+00 
  c          2   q    -7.00    -7.00   0.00e+00 
  c          2   r    +3.00    +3.00   0.00e+00 
  t          3   p    +6.00    +6.00   0.00e+00 
  t          3   q    -7.00    -7.00   0.00e+00 
  t          3   r    +4.00    +4.00   1.78e-15 
  d          1   p    -0.50    -0.50   8.88e-16 
  d          1   q    -8.00    -8.00   0.00e+00 
  d          1   r    -2.00    -2.00   0.00e+00 
  s          2   p    +1.50    +1.50   8.88e-16 
  s          2   q    -7.00    -7.00   0.00e+00 
  s          2   r    +0.00    +0.00   4.44e-16 
  b          3   p    +1.50    +1.50   8.88e-16 
  b          3   q    -6.00    -6.00   0.00e+00 
  b          3   r    +4.00    +4.00   1.78e-15 

  Total absolute error: 1.53e-14
  Mean error per coordinate: 5.67e-16
  → Parabolic model is EXACT ✓

========================================================================
  §5. CONSTRAINT ANALYSIS & DEGREES OF FREEDOM
========================================================================

  Total parameters: 3 sectors × 3 vectors (a,b,c) × 3 coords = 27

  Constraints:

  Constraint                                         Count
  ────────────────────────────────────────────────── ─────
  Electron origin: c_ℓ = -(a_ℓ + b_ℓ)                    3
  Flatness: b_{k(s)} = -5·a_{k(s)}, 1 per sector         3
  q-unification: q_ℓ(3) = q_u(3)                         1
  Half-integer quantization (discrete)                   0
  ────────────────────────────────────────────────── ─────
  Total continuous constraints                           7
  Total parameters                                      27
  Free parameters                                       20

========================================================================
  §6. BOSON EXTENSION
========================================================================

  The bosons W, Z, H are connected via spanning tree: t→W→H→Z.
  Can the parabolic structure extend to bosons?

  t→W: Δ(p,q,r) = (-0.5, -3, -2)
           Δ(ln m) = -0.7598
           m_W/m_t = 0.467756

  W→H: Δ(p,q,r) = (+1.5, +1, +0)
           Δ(ln m) = +0.4451
           m_H/m_W = 1.560580

  H→Z: Δ(p,q,r) = (+1.0, -2, -2)
           Δ(ln m) = -0.3148
           m_Z/m_H = 0.729970

  Can bosons be assigned a 'generation' g on any fermion parabola?

  Testing lepton parabola:
    W:   p→g=['4.212', '-0.712']  q→g=['6.000', '5.000']  r→g=['1.543', '3.457']
    Z:   p→g=['4.676', '-1.176']  q→no real g  r→g=['1.000', '4.000']
    H:   p→g=['4.500', '-1.000']  q→g=['7.000', '4.000']  r→g=['1.543', '3.457']

  Testing up parabola:
    W:   p→g=['2.930', '-1.073']  q→no real g  r→g=['1.667', '4.000']
    Z:   p→g=['3.260', '-1.402']  q→no real g  r→g=['1.192', '4.475']
    H:   p→g=['3.134', '-1.276']  q→no real g  r→g=['1.667', '4.000']

  Testing down parabola:
    W:   p→no real g  q→g=['-1.000']  r→g=['2.562', '-1.562']
    Z:   p→no real g  q→g=['-2.000']  r→g=['2.000', '-1.000']
    H:   p→no real g  q→g=['0.000']  r→g=['2.562', '-1.562']

========================================================================
  §7. UNIVERSALITY ANALYSIS
========================================================================

  Question: Which components of (a, b, c) are "universal" (same across sectors)
  and which are sector-specific?


  a (curvature):
    p: ℓ=+1.0000, u=+1.7500, d=-1.0000  
    q: ℓ=+0.5000, u=+0.5000, d=+0.0000  ℓ=u
    r: ℓ=-1.5000, u=-1.5000, d=+1.0000  ℓ=u

  b (velocity):
    p: ℓ=-3.5000, u=-3.2500, d=+5.0000  
    q: ℓ=-5.5000, u=-2.5000, d=+1.0000  
    r: ℓ=+7.5000, u=+8.5000, d=-1.0000  

  c (offset):
    p: ℓ=+2.5000, u=+0.0000, d=-4.5000  
    q: ℓ=+5.0000, u=-4.0000, d=-9.0000  
    r: ℓ=-6.0000, u=-8.0000, d=-2.0000  

  Koide-like mass ratios from scalar curvature:
    lepton: A = -1.252103
      m3·m1/m2² = e^(2A) = 0.081741
    up: A = -0.732242
      m3·m1/m2² = e^(2A) = 0.231197
    down: A = +0.405465
      m3·m1/m2² = e^(2A) = 2.250000


  Predicted ln(m/m_e) from parabolic model:
  Particle  ln(m/m_e)_model     ln(m/m_e)_3D          Δ
  ──────── ──────────────── ──────────────── ──────────
  e               -0.000000        +0.000000   1.20e-15
  mu              +5.327915        +5.327915   0.00e+00
  tau             +8.151625        +8.151625   3.55e-15
  u               +1.429645        +1.429645   2.00e-15
  c               +7.805052        +7.805052   0.00e+00
  t              +12.715973       +12.715973   1.78e-15
  d               +2.213506        +2.213506   8.88e-16
  s               +5.202362        +5.202362   0.00e+00
  b               +9.002148        +9.002148   0.00e+00

========================================================================
  SUMMARY
========================================================================

  §1: Curvature tensor a_k(I₃, Nc)
      a_p: unique solution α=1.3750, 
           β=-1.0000, 
           γ=3.3750

  §2: Flatness → velocity
      3 of 9 velocity components fixed: b_k(s) = -5·a_k(s)
      All verified exact ✓

  §3: Electron origin → offset
      c_ℓ fully determined from a_ℓ, b_ℓ
      q-unification verified: 3·Δb_q + Δc_q = 0.00

  §4: Full verification
      All 27 fermion coordinates reproduced exactly

  §5: DOF count
      27 total - 7 constraints = 20 free parameters

  §6: Boson extension
      Spanning tree chain t→W→H→Z analyzed
      Bosons do not lie on any fermion parabola (as expected)

  §7: Universality
      a_q and a_r partially universal (ℓ = u pattern)
      a_p fully sector-specific

