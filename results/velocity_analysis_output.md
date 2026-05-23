========================================================================
  VELOCITY ANALYSIS — Integration & Boundary Conditions
========================================================================

§1. FORMAL INTEGRATION — From Acceleration to Mass
------------------------------------------------------------------------

  The mass parabola: ln(m(g)/mₑ) = A·g² + B·g + C

  Reading this as kinematics in "generation time" g:
    Position:     x(g) = A·g² + B·g + C     [= ln(m/mₑ)]
    Velocity:     v(g) = 2A·g + B            [= d/dg ln(m/mₑ)]
    Acceleration: a(g) = 2A = constant       [= d²/dg² ln(m/mₑ)]

  The acceleration 2A is CONSTANT (uniform) — each sector has a
  fixed "generational acceleration" in mass-space.

  INTEGRATION PROCESS:
  ────────────────────
  Given: acceleration 2A (from graph topology)
  
  First integration:
    v(g) = ∫ 2A dg = 2A·g + B
    → B is the FIRST integration constant = initial velocity v(0)
    
  Second integration:
    x(g) = ∫ v(g) dg = A·g² + B·g + C
    → C is the SECOND integration constant = initial position x(0)

  To fix B and C we need TWO boundary conditions per sector.


§2. BOUNDARY CONDITIONS — What Fixes the Integration Constants?
------------------------------------------------------------------------

  Each sector needs 2 boundary conditions. Natural choices:

  OPTION A: "Fix gen-1 and gen-2" (know the 2 lightest masses)
    x(1) = ln(m₁/mₑ)  →  A + B + C = known
    x(2) = ln(m₂/mₑ)  →  4A + 2B + C = known

    Solving: B = ln(m₂/mₑ) - 3A - ln(m₁/mₑ) = ln(m₂/m₁) - 3A
             C = ln(m₁/mₑ) - A - B

    → gen-3 mass is PREDICTED:
      ln(m₃/mₑ) = 9A + 3B + C = 2A + 2·ln(m₂/mₑ) - ln(m₁/mₑ)

    Equivalently:  ln(m₃/m₁) = 2A + 2·ln(m₂/m₁)
    Or:            m₃/m₁ = e^(2A) · (m₂/m₁)²

  ★ PREDICTION FORMULA: m₃/m₁ = e^(2A) · (m₂/m₁)²

  Sector              A     e^(2A)     (m₂/m₁)²    Predicted   Actual m₃/m₁    Error
  ──────────────────────────────────────────────────────────────────────
  lepton      -1.252103   0.081741     42752.64      3494.62        3477.22  +0.5005%
  up          -0.732242   0.231197    345700.45     79924.94       79981.48  -0.0707%
  down        +0.405465   2.250000       400.00       900.00         895.07  +0.5502%

  The formula m₃/m₁ = e^(2A) · (m₂/m₁)² is EXACT by construction
  (3 points determine a unique quadratic). The errors are from
  the lattice discretization (~0.3%).

  BUT: this formula tells us something profound.
  Rearranging: m₃ · m₁ / m₂² = e^(2A)

  ★ The "Koide-like" ratio m₃·m₁/m₂² is FIXED by the curvature A!
    It equals e^(2A) — a pure function of graph topology.

  Koide-like geometric ratio m₃·m₁/m₂²:
  lepton    : m₃·m₁/m₂² = 0.081333,  e^(2A) = 0.081741
  up        : m₃·m₁/m₂² = 0.231361,  e^(2A) = 0.231197
  down      : m₃·m₁/m₂² = 2.237687,  e^(2A) = 2.250000


§3. CURVATURE A FROM GRAPH — Decomposition
------------------------------------------------------------------------

  A = aₚ·ln2 + a_q·lnR₃ + a_r·ln3

  Constraints determine a_q, a_r:
    ℓ/u: a_q = +0.5, a_r = -1.5  (isocurvature + conservation)
    d:   a_q = 0,    a_r = +1.0  (linearity + conservation)

  Only aₚ differs per sector:
    aₚ(ℓ) = 1.0,  aₚ(u) = 1.75,  aₚ(d) = -1.0

  So: A = aₚ·ln2 + FIXED_PART
  where FIXED_PART depends only on the sector type (normal/inverted).

  Fixed part (ℓ/u): 0.5·lnR₃ − 1.5·ln3 = -1.945250
  Fixed part (d):   0.0·lnR₃ + 1.0·ln3 = 1.098612

  A_l = +1.00·ln2 + -1.945250 = 0.693147 + -1.945250 = -1.252103
  A_u = +1.75·ln2 + -1.945250 = 1.213008 + -1.945250 = -0.732242
  A_d = -1.00·ln2 + 1.098612 = -0.693147 + 1.098612 = 0.405465

  ★ The ONLY free parameter in the curvature A is aₚ (per sector).
    If we can derive aₚ from quantum numbers, A is FULLY determined.

  Can aₚ = f(I₃, Nₒ)?
    ℓ: I₃ = −1/2, Nₒ = 1  → aₚ = 1.00
    u: I₃ = +1/2, Nₒ = 3  → aₚ = 1.75
    d: I₃ = −1/2, Nₒ = 3  → aₚ = −1.00

  Linear fit: aₚ = 1.3750·(2I₃) + -1.0000·Nₒ + 3.3750

  As fractions: aₚ = (11/8)·(2I₃) + (-1)·Nₒ + (27/8)

  ℓ: aₚ = 11/8·(-1) + (-1)·1 + (27/8) = 1.0000  [actual: 1.0000]  ✓
  u: aₚ = 11/8·(+1) + (-1)·3 + (27/8) = 1.7500  [actual: 1.7500]  ✓
  d: aₚ = 11/8·(-1) + (-1)·3 + (27/8) = -1.0000  [actual: -1.0000]  ✓

  This works but has 3 parameters for 3 data points (trivially exact).
  
  HOWEVER: the coefficients 11/8, -1, 27/8 are rational.
  Let's check if they relate to our building blocks:
    11/8 = ?
    -1 = ?  
    27/8 = ?

  Alternative: aₚ vs electric charge Q:
    ℓ: Q = -1.00, aₚ = +1.0000
    u: Q = +0.67, aₚ = +1.7500
    d: Q = -0.33, aₚ = -1.0000
  Linear fit: aₚ = 0.6316·Q + 0.7237
  Residual norm: 1.8655463430037693

  Full decomposition with (2I₃) and Nₒ:

  aₚ = 11/8·(2I₃) + (-1)·Nₒ + 27/8

  The (2I₃) term = 11/8 could be 11/8
  The Nₒ term = -1 (negative: more color → more negative → lighter gen-2!)
  The constant = 27/8


§4. VELOCITY B — What Fixes the First Integration Constant?
------------------------------------------------------------------------

  B = bₚ·ln2 + b_q·lnR₃ + b_r·ln3

  Velocity in mass-space (= initial log-mass growth rate):

  Sector      B (mass v.)       v(g=1)       v(g=2)       v(g=3)
  ───────────────────────────────────────────────────────
  lepton        +9.084224      +6.5800      +4.0758      +1.5716
  up            +8.572134      +7.1076      +5.6432      +4.1787
  down          +1.772461      +2.5834      +3.3943      +4.2053

  v(g) = 2A·g + B = generational mass velocity at generation g.
  
  v(g) > 0 means mass INCREASES with generation (always true for g=1,2,3).
  v decreases for ℓ and u (negative A → decelerating!).
  v increases for down (positive A → accelerating!).

  Lattice velocity vectors b = (bₚ, b_q, b_r):
  Sector           bₚ      b_q      b_r        |b|
  ────────────────────────────────────────
  lepton        -3.50    -5.50    +7.50     9.9373
  up            -3.25    -2.50    +8.50     9.4373
  down          +5.00    +1.00    -1.00     5.1962

  Angles between velocity vectors:
    lepton↔up: 18.6°
    lepton↔down: 126.2°
    up↔down: 123.8°


§5. MINIMAL INPUT — How Many Masses Must We Know?
------------------------------------------------------------------------

  SCENARIO 1: A from graph, know gen-1 and gen-2 masses → predict gen-3
  ─────────────────────────────────────────────────────────────────────
  Input:  3 curvatures (A) + 3 gen-1 masses + 3 gen-2 masses = 9
  Output: 3 gen-3 masses (predicted)
  Status: TRIVIALLY EXACT (3 points → unique quadratic)

  SCENARIO 2: A from graph, know gen-1 masses → predict gen-2 AND gen-3
  ─────────────────────────────────────────────────────────────────────
  Input:  3 curvatures + 3 gen-1 masses = 6 parameters
  Output: 6 masses (3 gen-2 + 3 gen-3)
  Status: NEED 3 velocity values B → 3 additional constraints
  
  SCENARIO 3: A from graph, B from regularity → only gen-1 masses needed
  ─────────────────────────────────────────────────────────────────────
  Input:  3 curvatures + 3 velocities + 3 gen-1 masses = 9 (but A,B from graph!)
  Output: 6 masses from gen-1 anchors alone
  Status: IDEAL — but what regularity principle fixes B?

  SCENARIO 4: Everything from graph (no experimental input)
  ─────────────────────────────────────────────────────────
  Input:  A from graph + B from regularity + C from quantum numbers
  Output: ALL 9 masses (relative to mₑ)
  Status: ULTIMATE GOAL

  SCENARIO 2 in practice — B as function of gen-1 mass:
  Sector              A    ln(m₁/mₑ)          B     B + 3A
  ───────────────────────────────────────────────────────
  lepton      -1.252103    +0.000000  +9.084224  +5.327915
  up          -0.732242    +1.441494  +8.572134  +6.375407
  down        +0.405465    +2.212545  +1.772461  +2.988856

  B + 3A = v(1) + A = velocity at g=1 plus half the acceleration.
  This equals ln(m₂/m₁) = logarithmic mass jump from gen-1 to gen-2.

  If we could derive B from A (or from graph invariants alone),
  knowing the gen-1 mass ratio would give ALL masses.


§6. VELOCITY REGULARITY — Testing Hypotheses for B
------------------------------------------------------------------------

  Hypothesis 1: Is B proportional to A?
  (Constant ratio → same "damping" across sectors)

  lepton    : B/A = 9.084224/-1.252103 = -7.2552
  up        : B/A = 8.572134/-0.732242 = -11.7067
  down      : B/A = 1.772461/0.405465 = 4.3714


  Hypothesis 2: Is the arc-length per generation constant?
  Arc length = ∫|γ'(g)|dg in lattice space

  lepton    : arc₁₂ = 5.0249, arc₂₃ = 3.3541, ratio = 0.6675
  up        : arc₁₂ = 4.5826, arc₂₃ = 5.5902, ratio = 1.2199
  down      : arc₁₂ = 3.0000, arc₂₃ = 4.1231, ratio = 1.3744


  Hypothesis 3: Minimal curvature energy (∫κ²ds)?
  For a quadratic Bézier, curvature is related to |a|.

  lepton    : |a|² = 3.5000, |a| = 1.8708
  up        : |a|² = 5.5625, |a| = 2.3585
  down      : |a|² = 2.0000, |a| = 1.4142


  Hypothesis 4: Does B live on the lattice?
  Check if B-components are expressible in the basis {ln2, lnR₃, ln3}.

  b-coefficients (should be rational if lattice-native):
  Sector             bₚ        b_q        b_r
  ──────────────────────────────────────────
  lepton          -3.50      -5.50      +7.50
  up              -3.25      -2.50      +8.50
  down            +5.00      +1.00      -1.00

  Smallest |b| = 1.0
  All b as multiples of 0.25:
  lepton    :    -14/4,    -22/4,    +30/4
  up        :    -13/4,    -10/4,    +34/4
  down      :    +20/4,     +4/4,     -4/4


§7. VELOCITY FROM CURVATURE — Are b and a Related?
------------------------------------------------------------------------

  The curvature (a) and velocity (b) vectors per sector:
  If the lattice imposes a relationship b = f(a), we can derive B from A.

  Sector         aₚ    a_q    a_r │       bₚ      b_q      b_r │  b/a ratio
  ───────────────────────────────────────────────────────────────────────────
  lepton      +1.00  +0.50  -1.50 │    -3.50    -5.50    +7.50 │ -3.50 -11.00 -5.00
  up          +1.75  +0.50  -1.50 │    -3.25    -2.50    +8.50 │ -1.86 -5.00 -5.67
  down        -1.00  +0.00  +1.00 │    +5.00    +1.00    -1.00 │ -5.00     ∞ -1.00

  Testing b + k·a for small integers k:
  k=-6: ℓ=( -9.50, -8.50,+16.50)  u=(-13.75, -5.50,+17.50)  d=(+11.00, +1.00, -7.00)
  k=-5: ℓ=( -8.50, -8.00,+15.00)  u=(-12.00, -5.00,+16.00)  d=(+10.00, +1.00, -6.00)
  k=-4: ℓ=( -7.50, -7.50,+13.50)  u=(-10.25, -4.50,+14.50)  d=( +9.00, +1.00, -5.00)
  k=-3: ℓ=( -6.50, -7.00,+12.00)  u=( -8.50, -4.00,+13.00)  d=( +8.00, +1.00, -4.00)
  k=-2: ℓ=( -5.50, -6.50,+10.50)  u=( -6.75, -3.50,+11.50)  d=( +7.00, +1.00, -3.00)
  k=-1: ℓ=( -4.50, -6.00, +9.00)  u=( -5.00, -3.00,+10.00)  d=( +6.00, +1.00, -2.00)
  k=+0: ℓ=( -3.50, -5.50, +7.50)  u=( -3.25, -2.50, +8.50)  d=( +5.00, +1.00, -1.00)
  k=+1: ℓ=( -2.50, -5.00, +6.00)  u=( -1.50, -2.00, +7.00)  d=( +4.00, +1.00, +0.00)
  k=+2: ℓ=( -1.50, -4.50, +4.50)  u=( +0.25, -1.50, +5.50)  d=( +3.00, +1.00, +1.00)
  k=+3: ℓ=( -0.50, -4.00, +3.00)  u=( +2.00, -1.00, +4.00)  d=( +2.00, +1.00, +2.00)
  k=+4: ℓ=( +0.50, -3.50, +1.50)  u=( +3.75, -0.50, +2.50)  d=( +1.00, +1.00, +3.00)
  k=+5: ℓ=( +1.50, -3.00, +0.00)  u=( +5.50, +0.00, +1.00)  d=( +0.00, +1.00, +4.00)
  k=+6: ℓ=( +2.50, -2.50, -1.50)  u=( +7.25, +0.50, -0.50)  d=( -1.00, +1.00, +5.00)

  Tangent at g=3 (unification point): v(3) = 6a + b
  lepton    : v(3) = ( +2.50,  -2.50,  -1.50)
  up        : v(3) = ( +7.25,  +0.50,  -0.50)
  down      : v(3) = ( -1.00,  +1.00,  +5.00)

  Tangent at g=0 (extrapolated pre-generation): v(0) = b
  lepton    : v(0) = ( -3.50,  -5.50,  +7.50)
  up        : v(0) = ( -3.25,  -2.50,  +8.50)
  down      : v(0) = ( +5.00,  +1.00,  -1.00)

  Position at g=0 (extrapolated): γ(0) = c
  lepton    : γ(0) = ( +2.50,  +5.00,  -6.00)
  up        : γ(0) = ( +0.00,  -4.00,  -8.00)
  down      : γ(0) = ( -4.50,  -9.00,  -2.00)


§8. VELOCITY CONSTRAINT — Unification at g=3 Restricts B
------------------------------------------------------------------------

  We know q_ℓ(3) = q_u(3) = −7 (exact unification at gen 3).
  
  This constrains the q-velocity:
    q(g) = a_q·g² + b_q·g + c_q
    q(3) = 9·a_q + 3·b_q + c_q = −7
    
  For ℓ: 9(0.5) + 3(-5.5) + 5.0 = 4.5 - 16.5 + 5.0 = -7.0 ✓
  For u: 9(0.5) + 3(-2.5) +(-4.0) = 4.5 - 7.5 - 4.0 = -7.0 ✓
  
  These share a_q = 0.5 but have DIFFERENT b_q: -5.5 vs -2.5.
  The difference b_q(u) − b_q(ℓ) = 3.0 is compensated by 
  the difference c_q(u) − c_q(ℓ) = -9.0.
  
  Constraint: 3·Δb_q + Δc_q = 0 → Δc_q = −3·Δb_q
  (The offset absorbs 3× the velocity difference to meet at g=3.)

  Δb_q = +3.0, Δc_q = -9.0, −3·Δb_q = -9.0  ✓

  ★ The q-unification at g=3 couples velocity and offset:
    "faster start" (larger b_q) requires "deeper offset" (more negative c_q)
    to arrive at the same point at generation 3.
    
  This is a BOUNDARY VALUE problem, not an initial value problem!
  The sectors are like beams anchored at both ends (g=1 and g=3),
  not free projectiles.


§9. RETHINKING — From IVP to BVP (Boundary Value Problem)
------------------------------------------------------------------------

  We've been thinking "initial value": know acceleration + initial conditions.
  But the LATTICE suggests "boundary value": know endpoints + curvature.

  For each sector, the curve passes through 3 lattice points.
  These are NOT arbitrary — they are constrained by:
  
  1. INTEGER/HALF-INTEGER coordinates (lattice quantization)
  2. Curvature conservation |a_q + a_r| = 1
  3. q-unification at g=3 (at least ℓ↔u)
  4. Electron at origin (gen-1 lepton)
  
  So the CORRECT formulation is:
  
  ┌─────────────────────────────────────────────────────────┐
  │  Given: lattice {ln2, lnR₃, ln3}                       │
  │         curvature law |a_q + a_r| = 1                   │
  │         quantization: coordinates ∈ ½ℤ                  │
  │         electron = (0,0,0)                              │
  │         unification q_ℓ(3) = q_u(3)                     │
  │                                                         │
  │  Find: all lattice positions that satisfy these AND      │
  │        produce physical mass ordering m₁ < m₂ < m₃      │
  └─────────────────────────────────────────────────────────┘
  
  This is a DISCRETE optimization: find integer points on
  quadratic Bézier curves with fixed curvature constraints!

  LATTICE QUANTIZATION CHECK:
  Are all coordinates in ½ℤ?

  b    : ( +1.5,  -6.0,  +4.0)  ✓
  c    : ( +0.5,  -7.0,  +3.0)  ✓
  d    : ( -0.5,  -8.0,  -2.0)  ✓
  e    : ( +0.0,  +0.0,  +0.0)  ✓
  mu   : ( -0.5,  -4.0,  +3.0)  ✓
  s    : ( +1.5,  -7.0,  +0.0)  ✓
  t    : ( +6.0,  -7.0,  +4.0)  ✓
  tau  : ( +1.0,  -7.0,  +3.0)  ✓
  u    : ( -1.5,  -6.0,  -1.0)  ✓

  All ∈ ½ℤ: YES

  Finer check — coordinate types:
  Types found: {'integer', 'half-integer'}


§10. DISCRETE SEARCH — Counting Compatible Lattice Configurations
------------------------------------------------------------------------

  Setup: Fix the curvatures (a_p, a_q, a_r) per sector.
  Then for each sector, the gen-2 particle is determined by gen-1 and gen-3:
    γ(1) + γ(3) − 2·γ(2) = 2a  ⇒  coord_2 = (coord_1 + coord_3)/2 − a
    
  So the FREE choices are the gen-1 and gen-3 positions (gen-2 follows).
  
  Constraints:
    - All coordinates in ½ℤ (after applying curvature)
    - Mass ordering: m₁ < m₂ < m₃ within each sector
    - Lepton gen-1 = (0,0,0)
    - q_ℓ(3) = q_u(3) (unification)
    
  Search space: gen-1 and gen-3 positions for each sector.
  Let's enumerate all valid configurations in a bounded region
  and see how constrained the system really is.

  LEPTON sector: gen-1 = (0,0,0), gen-3 q = −7
  Curvature: a = (+1.0, +0.5, -1.5)
  Search range: p₃ ∈ [-3.0, 10.0], r₃ ∈ [-3.0, 10.0]
  Valid configurations (m₁ < m₂ < m₃, half-integer): 182

  ★ ACTUAL lepton configuration:
    gen-1: (np.float64(0.0), np.float64(0.0), np.float64(0.0)) → m = 0.511 MeV
    gen-2: (np.float64(-0.5), np.float64(-4.0), np.float64(3.0)) → m = 105.270 MeV
    gen-3: (np.float64(1.0), np.float64(-7.0), np.float64(3.0)) → m = 1772.665 MeV

  Solutions with m₃ within 10% of m_τ (1776.9 MeV): 2
  Solutions with BOTH m₂ ≈ m_μ AND m₃ ≈ m_τ (10%): 2

  These are:
    gen-3 = (+1.0, −7, +3.0), m₂ = 105.3 MeV, m₃ = 1772.7 MeV
    gen-3 = (+9.0, −7, -2.0), m₂ = 108.0 MeV, m₃ = 1867.5 MeV

  ★ CONSTRAINT POWER:
    Before constraints: 27 × 27 = 729 lattice points
    After mass ordering: 182 valid
    After mass matching (10%): 2 valid
    
    The lattice + curvature law reduces the solution space dramatically.


§11. VELOCITY DECOMPOSITION — Can We Express b?
------------------------------------------------------------------------

  Like aₚ = f(I₃, Nc), can we express bₚ, b_q, b_r in terms
  of quantum numbers and/or graph invariants?
  
  The b-coefficients:
    Lepton: (−3.5, −5.5, +7.5)
    Up:     (−3.25, −2.5, +8.5)
    Down:   (+5.0, +1.0, −1.0)
    
  Differences (relative to lepton):
    Δb(u−ℓ): (+0.25, +3.0, +1.0)
    Δb(d−ℓ): (+8.5, +6.5, −8.5)

  Δa(u−ℓ) = (+0.75, +0.00, +0.00)
  Δb(u−ℓ) = (+0.25, +3.00, +1.00)

  r_u(3) − r_ℓ(3) = +1.0
  (ℓ and u are close at gen 3 in r too — gap of only 1)

  EXACT unification test at g=3:
  p: ℓ=+1.0, u=+6.0, d=+1.5, Δ(u−ℓ)=+5.0, Δ(d−ℓ)=+0.5
  q: ℓ=-7.0, u=-7.0, d=-6.0, Δ(u−ℓ)=+0.0, Δ(d−ℓ)=+1.0
  r: ℓ=+3.0, u=+4.0, d=+4.0, Δ(u−ℓ)=+1.0, Δ(d−ℓ)=+1.0

  At g = 3: q-unification is EXACT (Δ=0), r-gap = 1, p-gap very large.
  
  The q-unification is a HARD constraint: 3·Δb_q + Δc_q = 0 (ℓ↔u).
  The r near-unification is SOFT: gap = 1.
  
  If we ALSO imposed r-unification, Δb_r would be further constrained.
  Currently Δb_r(u−ℓ) = 1.0, Δc_r(u−ℓ) = −2.0, giving gap = 3(1)+(-2) = 1.
  For exact: need Δc_r = −3·Δb_r → Δb_r = Δb_r, Δc_r = −3·Δb_r
  With current Δb_r = 1: need Δc_r = −3 (currently −2).
  → Missing 1 unit of offset.

  ★ The r-gap of exactly 1 at gen 3 may be the TOPOLOGICAL UNIT
    that prevents full C¹ continuity — the "quantum" of spectral gap.


========================================================================
  §12. SUMMARY — The Status of Mass Derivation
========================================================================

  WHAT WE HAVE (from graph topology alone):
  ──────────────────────────────────────────
  ✓ 3D lattice basis: {ln2, lnR₃, ln3}
  ✓ Curvature law: |a_q + a_r| = 1
  ✓ Isocurvature: a_q(ℓ) = a_q(u), a_r(ℓ) = a_r(u)
  ✓ Down linearity: a_q(d) = 0
  ✓ q-unification: q_ℓ(3) = q_u(3) = −7
  ✓ Electron = origin
  ✓ Quantization: all coords ∈ ½ℤ

  WHAT THIS DETERMINES:
  ─────────────────────
  ✓ Curvature vector a per sector (given aₚ = f(I₃,Nc))
  ✓ Relationship: gen-2 = (gen-1 + gen-3)/2 − a  (Bézier midpoint)
  ✓ BVP coupling: Δc_q = −3·Δb_q between ℓ and u

  WHAT REMAINS FREE:
  ──────────────────
  ? aₚ per sector: expressible as 11/8·(2I₃)+(-1)·Nc+27/8
    but not yet derived from first principles
  ? Velocity b per sector per component (9 params, reduced by constraints)
  ? Gen-1 positions for up/down sectors (4 coords, after e=origin)
  ? Gen-3 positions (constrained by unification for q)

  THE CORRECT PICTURE:
  ────────────────────
  This is NOT an initial value problem (know start → predict end).
  This is a BOUNDARY VALUE problem on a discrete lattice:
  
  Find half-integer coordinates that satisfy:
    1. Curvature constraints (from gauge structure)
    2. Mass ordering m₁ < m₂ < m₃ (from stability)
    3. q-unification at gen 3 (from spectral convergence)
    4. Bézier midpoint relation (quadratic structure)
    
  The solution space is HIGHLY constrained:
  for leptons, only 2 configuration(s) match
  the observed masses within 10%.
  
  NEXT STEP: Run the discrete search for ALL THREE sectors
  simultaneously, with inter-sector constraints from CKM mixing
  (the C⁰ crossing points). This should further reduce solutions.

