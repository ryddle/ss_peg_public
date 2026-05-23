========================================================================
  §1. MASSLESS BOSONS: WHY γ AND GLUONS ARE OUTSIDE THE LATTICE
========================================================================

  The mass lattice encodes POLE MASSES of observable particles via:
    ln(m_f / m_e) = p·ln2 + q·lnR₃ + r·ln3

  For a massless particle, m = 0 ⟹ ln(m/m_e) = −∞.
  This requires (p, q, r) → −∞ in some direction — no finite lattice point.

  ┌──────────┬────────┬──────────────┬──────────────────────────────────┐
  │ Particle │ Mass   │ Lattice?     │ Reason                           │
  ├──────────┼────────┼──────────────┼──────────────────────────────────┤
  │ Photon γ │ 0      │ NO (∞ point) │ Unbroken U(1)_EM gauge boson     │
  │ Gluons g │ 0*     │ NO (∞ point) │ Unbroken SU(3)_c, confinement    │
  │ W±       │ 80 GeV │ YES          │ Broken SU(2)_L via Higgs         │
  │ Z⁰       │ 91 GeV │ YES          │ Broken SU(2)_L×U(1)_Y via Higgs │
  │ H⁰       │ 125 GeV│ YES          │ Higgs scalar                     │
  └──────────┴────────┴──────────────┴──────────────────────────────────┘

  * Gluon mass = 0 in the Lagrangian. They are never observed as free
    particles (confinement). The "mass" concept is ill-defined for gluons.
    ΛQCD ≈ 200 MeV is a scale, not a pole mass.

  CONCLUSION: The lattice naturally separates:
  - Massive particles (broken gauge symmetry) → finite lattice points
  - Massless particles (unbroken gauge symmetry) → at infinity
  
  This is consistent with the Higgs mechanism: mass ↔ lattice membership.

========================================================================
  §2. BOSON QUANTUM NUMBERS AND PARABOLA ASSIGNMENT
========================================================================

  Massive bosons and their SM quantum numbers:

  ┌──────────┬────────┬────┬──────┬─────┬─────────────────────────────┐
  │ Particle │ Mass   │ I  │ I₃   │ Nc  │ Notes                       │
  ├──────────┼────────┼────┼──────┼─────┼─────────────────────────────┤
  │ W⁺       │ 80.4   │ 1  │ +1   │ 1   │ SU(2) adjoint, color singlet│
  │ W⁻       │ 80.4   │ 1  │ −1   │ 1   │ SU(2) adjoint, color singlet│
  │ Z⁰       │ 91.2   │ 0* │ 0    │ 1   │ Mixed W³/B, color singlet   │
  │ H⁰       │ 125.1  │ ½  │ −½   │ 1   │ SU(2) doublet, color singlet│
  └──────────┴────────┴────┴──────┴─────┴─────────────────────────────┘
  * Z is a mixture, effective I₃ = 0.

  Key: H⁰ has (2I₃=−1, Nc=1) — SAME as lepton sector!
  But leptons have g ∈ {1,2,3} and H has much larger mass.
  
  Boson lattice coordinates:

    W: (p, q, r) = (11/2, -10, 2)
    Z: (p, q, r) = (8, -11, 0)
    H: (p, q, r) = (7, -9, 2)

========================================================================
  §3. EFFECTIVE GENERATION ON FERMION PARABOLAS
========================================================================

  For each boson B with coordinates (p_B, q_B, r_B), and each fermion
  sector s, solve for g_eff from:
    x_k(g) = a_k(s)·g² + b_k(s)·g + c_k(s) = x_B,k
  for each k ∈ {p, q, r}. If all 3 coords give the same g, the boson
  lies on that fermion parabola.

  ── W: target = (11/2, -10, 2) ──
    lepton parabola (2I₃=-1, Nc=1):
      p: g = 4.212214, -0.712214
      q: g = 6.000000, 5.000000
      r: g = 1.542573, 3.457427
      → No consistent g. Closest: (4.2122, 5.0000, 3.4574), spread=1.5426
    up parabola (2I₃=1, Nc=3):
      p: g = 2.929847, -1.072704
      q: g = no real roots
      r: g = 1.666667, 4.000000
      → No consistent g (some coords have no real roots)
    down parabola (2I₃=-1, Nc=3):
      p: g = no real roots
      q: g = -1
      r: g = 2.561553, -1.561553
      → No consistent g (some coords have no real roots)

  ── Z: target = (8, -11, 0) ──
    lepton parabola (2I₃=-1, Nc=1):
      p: g = 4.676175, -1.176175
      q: g = no real roots
      r: g = 1.000000, 4.000000
      → No consistent g (some coords have no real roots)
    up parabola (2I₃=1, Nc=3):
      p: g = 3.259596, -1.402453
      q: g = no real roots
      r: g = 1.191857, 4.474810
      → No consistent g (some coords have no real roots)
    down parabola (2I₃=-1, Nc=3):
      p: g = no real roots
      q: g = -2
      r: g = 2, -1
      → No consistent g (some coords have no real roots)

  ── H: target = (7, -9, 2) ──
    lepton parabola (2I₃=-1, Nc=1):
      p: g = 4.500000, -1.000000
      q: g = 7.000000, 4.000000
      r: g = 1.542573, 3.457427
      → No consistent g. Closest: (4.5000, 4.0000, 3.4574), spread=1.0426
    up parabola (2I₃=1, Nc=3):
      p: g = 3.133621, -1.276478
      q: g = no real roots
      r: g = 1.666667, 4.000000
      → No consistent g (some coords have no real roots)
    down parabola (2I₃=-1, Nc=3):
      p: g = no real roots
      q: g = 0
      r: g = 2.561553, -1.561553
      → No consistent g (some coords have no real roots)

========================================================================
  §4. GENERATING FUNCTION WITH BOSONIC QUANTUM NUMBERS
========================================================================

  Instead of using fermion parabolas, evaluate the generating function
  directly with bosonic quantum numbers.
  
  For each boson B with SM quantum numbers (2I₃, Nc), compute
    x_k(g) = a_k(2I₃,Nc)·g² + b_k(2I₃,Nc)·g + c_k(2I₃,Nc)
  and solve for g from each coordinate.


  ── W⁺: (2I₃=2, Nc=1), target=(11/2, -10, 2) ──
    p(g) = 41/8·g² + -127/8·g + 37/4
    q(g) = 5/4·g² + -43/4·g + 25/2
    r(g) = -21/4·g² + 87/4·g + -15
    p = 11/2 → g = 2.839909, 0.257652
    q = -10 → g = 5.000000, 3.600000
    r = 2 → g = 1.045407, 3.097451
    → No consistent g. Best: (2.8399, 3.6000, 3.0975), spread=0.7601

  ── Z⁰: (2I₃=0, Nc=1), target=(8, -11, 0) ──
    p(g) = 19/8·g² + -61/8·g + 19/4
    q(g) = 3/4·g² + -29/4·g + 15/2
    r(g) = -11/4·g² + 49/4·g + -9
    p = 8 → g = 3.591539, -0.381012
    q = -11 → g = no real roots
    r = 0 → g = 0.928036, 3.526510

  ── H⁰: (2I₃=-1, Nc=1), target=(7, -9, 2) ──
    p(g) = 1·g² + -7/2·g + 5/2
    q(g) = 1/2·g² + -11/2·g + 5
    r(g) = -3/2·g² + 15/2·g + -6
    p = 7 → g = 4.500000, -1.000000
    q = -9 → g = 7.000000, 4.000000
    r = 2 → g = 1.542573, 3.457427
    → No consistent g. Best: (4.5000, 4.0000, 3.4574), spread=1.0426


  ── Systematic scan: try (2I₃, Nc) with Nc ∈ {0,1,2,...,8} ──
  For each boson, find (2I₃, Nc, g) that reproduces all 3 coordinates.

  W: target = (11/2, -10, 2)
    No exact (2I₃ ∈ [-4,4], Nc ∈ [0,8], rational g) solution found.

  Z: target = (8, -11, 0)
    No exact (2I₃ ∈ [-4,4], Nc ∈ [0,8], rational g) solution found.

  H: target = (7, -9, 2)
    No exact (2I₃ ∈ [-4,4], Nc ∈ [0,8], rational g) solution found.

========================================================================
  §5. EXTENDED SCAN — HALF-INTEGER & FRACTIONAL QUANTUM NUMBERS
========================================================================

  Extend the search to half-integer 2I₃ and fractional Nc (including
  Nc = 1/3, 8/3 etc. that might arise from representation theory).


  W: target = (11/2, -10, 2)
    No solution with half-integer 2I₃ and Nc/3.

  Z: target = (8, -11, 0)
    No solution with half-integer 2I₃ and Nc/3.

  H: target = (7, -9, 2)
    No solution with half-integer 2I₃ and Nc/3.

========================================================================
  §6. BOSON CHAIN ANALYSIS
========================================================================

  The bosons in the spanning tree form a chain connected to fermions:
    t → W → H → Z
  
  Analyze the difference vectors along this chain.

  Chain: t → W → H → Z
  Step           Δp       Δq       Δr        |Δ|
  ──────── ──────── ──────── ──────── ──────────
  t→W        -1/2       -3       -2     3.6401
  W→H         3/2        1        0     1.8028
  H→Z           1       -2       -2     3.0000

  Total t→Z: Δ = (2, -4, -4)
  |Δ_total| = 6.0000

  Bridge connections (boson ↔ fermion):
    tau→b: (1/2, 1, 1)
    b→t: (9/2, -1, 0)

========================================================================
  §7. BOSON AS FERMION GENERATION EXTRAPOLATION
========================================================================

  Alternative approach: forget quantum numbers.
  For each boson, find which fermion parabola it is closest to,
  and at what (possibly non-integer) generation.

  Minimize: Σ_k [x_k^(boson) - (a_k·g² + b_k·g + c_k)]² over g


  ── W: (11/2, -10, 2) ──
    lepton: g_opt = 3.9344, residual = 2.4096
      pred: (4.209, -8.899, 0.289)
      actual: (5.500, -10.000, 2.000)
    up: g_opt = 2.9148, residual = 3.5926
      pred: (5.395, -7.039, 4.032)
      actual: (5.500, -10.000, 2.000)
    down: g_opt = 2.3900, residual = 5.1093
      pred: (1.738, -6.610, 1.322)
      actual: (5.500, -10.000, 2.000)

  ── Z: (8, -11, 0) ──
    lepton: g_opt = 4.3718, residual = 2.9448
      pred: (6.311, -9.489, -1.880)
      actual: (8.000, -11.000, 0.000)
    up: g_opt = 3.2859, residual = 5.6122
      pred: (8.216, -6.816, 3.734)
      actual: (8.000, -11.000, 0.000)
    down: g_opt = 2.1014, residual = 7.6154
      pred: (1.591, -6.899, 0.315)
      actual: (8.000, -11.000, 0.000)

  ── H: (7, -9, 2) ──
    lepton: g_opt = 4.0496, residual = 3.1839
      pred: (4.726, -9.073, -0.227)
      actual: (7.000, -9.000, 2.000)
    up: g_opt = 3.1406, residual = 2.8178
      pred: (7.054, -6.920, 3.900)
      actual: (7.000, -9.000, 2.000)
    down: g_opt = 2.4430, residual = 5.8129
      pred: (1.747, -6.557, 1.525)
      actual: (7.000, -9.000, 2.000)

========================================================================
  §8. INTER-BOSON STRUCTURE
========================================================================

  Even if bosons don't sit on fermion parabolas, they may have their
  own internal structure. Test: do bosons form their own parabola?
  
  Assign W=gen1, H=gen2, Z=gen3 (ordered by distance from top):
    t→W: nearest, t→H: intermediate, t→Z: farthest


  Ordering: W→H→Z (by mass)
    g=1: W = (11/2, -10, 2)
    g=2: H = (7, -9, 2)
    g=3: Z = (8, -11, 0)
    p(g) = -1/4·g² + 9/4·g + 7/2
    q(g) = -3/2·g² + 11/2·g + -14
    r(g) = -1·g² + 3·g + 0
    Flatness check:
      p: a=-1/4, b=9/4, b/a=-9 
      q: a=-3/2, b=11/2, b/a=-11/3 
      r: a=-1, b=3, b/a=-3 

  Ordering: W→Z→H (by q coord)
    g=1: W = (11/2, -10, 2)
    g=2: Z = (8, -11, 0)
    g=3: H = (7, -9, 2)
    p(g) = -7/4·g² + 31/4·g + -1/2
    q(g) = 3/2·g² + -11/2·g + -6
    r(g) = 2·g² + -8·g + 8
    Flatness check:
      p: a=-7/4, b=31/4, b/a=-31/7 
      q: a=3/2, b=-11/2, b/a=-11/3 
      r: a=2, b=-8, b/a=-4 

  Ordering: H→W→Z (by p coord)
    g=1: H = (7, -9, 2)
    g=2: W = (11/2, -10, 2)
    g=3: Z = (8, -11, 0)
    p(g) = 2·g² + -15/2·g + 25/2
    q(g) = 0·g² + -1·g + -8
    r(g) = -1·g² + 3·g + 0
    Flatness check:
      p: a=2, b=-15/2, b/a=-15/4 
      q: a=0, b=-1 (linear) 
      r: a=-1, b=3, b/a=-3 

  Ordering: Z→H→W (reverse mass)
    g=1: Z = (8, -11, 0)
    g=2: H = (7, -9, 2)
    g=3: W = (11/2, -10, 2)
    p(g) = -1/4·g² + -1/4·g + 17/2
    q(g) = -3/2·g² + 13/2·g + -16
    r(g) = -1·g² + 5·g + -4
    Flatness check:
      p: a=-1/4, b=-1/4, b/a=1 
      q: a=-3/2, b=13/2, b/a=-13/3 
      r: a=-1, b=5, b/a=-5 → FLAT ✓

========================================================================
  SUMMARY
========================================================================

  1. MASSLESS BOSONS (γ, gluons): Live at infinity in the mass lattice.
     They are naturally excluded — the lattice describes massive particles
     from broken gauge symmetry only.

  2. FERMION PARABOLA FIT: Massive bosons (W, Z, H) do NOT sit on any
     fermion parabola at integer or simple rational generation.

  3. QUANTUM NUMBER SCAN: Systematic search over (2I₃, Nc) with extended
     ranges finds no combination that places bosons exactly on the
     generating function at any consistent generation.

  4. NEAREST APPROACH: The up-quark parabola provides the closest
     extrapolation for W and H (expected: t is the heaviest fermion).

  5. INTER-BOSON STRUCTURE: The 3 massive bosons may form their own
     "generation" sequence with separate parabolic parameters.

  6. BOSON CHAIN: The spanning tree chain t→W→H→Z has rational
     difference vectors, suggesting algebraic (not accidental) placement.

  CONCLUSION: Bosons require a SEPARATE generating function from fermions.
  The fermion function with SM quantum numbers cannot extend to bosons.
  This is physically natural: fermions and bosons obey different statistics
  and couple to the Higgs differently (Yukawa vs gauge coupling).


========================================================================
  §9. GENERATING FIGURES
========================================================================
  Figure 1 saved: D:\Desarrollo\projects\pyshics\ss_peg\output\boson_fermion_parabolas_components.png
  Figure 2 saved: D:\Desarrollo\projects\pyshics\ss_peg\output\boson_fermion_parabolas_3d.png
  Figure 3 saved: D:\Desarrollo\projects\pyshics\ss_peg\output\boson_flatness_highlight.png
  Figure 4 saved: D:\Desarrollo\projects\pyshics\ss_peg\output\boson_chain_3d.png
  Figure 5 saved: D:\Desarrollo\projects\pyshics\ss_peg\output\boson_fermion_flatness_comparison.png

  All figures generated ✓
