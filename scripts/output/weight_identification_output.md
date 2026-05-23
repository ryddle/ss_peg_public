========================================================================
  §1. ROOT SYSTEM RECONSTRUCTION
========================================================================

  Strategy: In a Lie algebra representation, the difference between
  consecutive weights along a root chain is a root α. We extract
  all gen→gen+1 differences as candidate roots.

  lepton:
    α₁₂ = e→mu = (-0.5, -4, +3)
    α₂₃ = mu→tau = (+1.5, -3, +0)
  up:
    α₁₂ = u→c = (+2.0, -1, +4)
    α₂₃ = c→t = (+5.5, +0, +1)
  down:
    α₁₂ = d→s = (+2.0, +1, +2)
    α₂₃ = s→b = (+0.0, +1, +4)

  Root relationships:

  Cross-sector roots (same generation):
    e→u (gen 1): (-1.5, -6, -1)
    e→d (gen 1): (-0.5, -8, -2)
    u→d (gen 1): (+1.0, -2, -1)
    mu→c (gen 2): (+1.0, -3, +0)
    mu→s (gen 2): (+2.0, -3, -3)
    c→s (gen 2): (+1.0, +0, -3)
    tau→t (gen 3): (+5.0, +0, +1)
    tau→b (gen 3): (+0.5, +1, +1)
    t→b (gen 3): (-4.5, +1, +0)

========================================================================
  §2. DYNKIN EMBEDDING TEST
========================================================================

  Testing: can we find a linear map L: R³→R³ such that L(pqr)
  are Dynkin labels (integers) of some rank-3 Lie algebra?

  For A₃: weights are (a₁, a₂, a₃) with aᵢ ∈ Z
  For B₃: similar but with different norms for short/long roots


  Sector-wise Cartan matrix construction:
  Using α₁₂ and α₂₃ as candidate simple roots per sector

  lepton:
    α₁ = [-0.5 -4.   3. ], |α₁| = 5.0249
    α₂ = [ 1.5 -3.   0. ], |α₂| = 3.3541
    |α₁|/|α₂| = 1.4981
    ⟨α₁,α₂⟩/⟨α₂,α₂⟩ × 2 = 2.0000
    ⟨α₁,α₂⟩/⟨α₁,α₁⟩ × 2 = 0.8911
    Cartan-like: [[2, 2.00], [0.89, 2]]
    Integer Cartan? False

  up:
    α₁ = [ 2. -1.  4.], |α₁| = 4.5826
    α₂ = [5.5 0.  1. ], |α₂| = 5.5902
    |α₁|/|α₂| = 0.8198
    ⟨α₁,α₂⟩/⟨α₂,α₂⟩ × 2 = 0.9600
    ⟨α₁,α₂⟩/⟨α₁,α₁⟩ × 2 = 1.4286
    Cartan-like: [[2, 0.96], [1.43, 2]]
    Integer Cartan? False

  down:
    α₁ = [2. 1. 2.], |α₁| = 3.0000
    α₂ = [0. 1. 4.], |α₂| = 4.1231
    |α₁|/|α₂| = 0.7276
    ⟨α₁,α₂⟩/⟨α₂,α₂⟩ × 2 = 1.0588
    ⟨α₁,α₂⟩/⟨α₁,α₁⟩ × 2 = 2.0000
    Cartan-like: [[2, 1.06], [2.00, 2]]
    Integer Cartan? True

========================================================================
  §3. SHIFTED LATTICE & WEYL SYMMETRY
========================================================================

  Weyl character formula uses shifted weights: λ + ρ where ρ = half
  sum of positive roots. Test symmetry of shifted weights.

  Fermion centroid: (0.8889, -5.7778, 1.5556)

  Shifted coordinates (w - centroid):
  Particle       Δp       Δq       Δr      |Δ|
  ──────── ──────── ──────── ──────── ────────
  e         -0.8889  +5.7778  -1.5556   6.0492
  mu        -1.3889  +1.7778  +1.4444   2.6788
  tau       +0.1111  -1.2222  +1.4444   1.8954
  u         -2.3889  -0.2222  -2.5556   3.5053
  c         -0.3889  -1.2222  +1.4444   1.9317
  t         +5.1111  -1.2222  +2.4444   5.7959
  d         -1.3889  -2.2222  -3.5556   4.4169
  s         +0.6111  -1.2222  -1.5556   2.0705
  b         +0.6111  -0.2222  +2.4444   2.5295

  Reflection symmetry test:
  For each shifted weight w, check if -w is also present (up to tolerance)
    e ↔ -e closest to tau (dist=4.6228) ✗
    mu ↔ -mu closest to s (dist=0.9623) ✗
    tau ↔ -tau closest to s (dist=2.5513) ✗
    u ↔ -u closest to b (dist=1.8359) ✗
    c ↔ -c closest to s (dist=2.4570) ✗
    t ↔ -t closest to u (dist=3.0837) ✗
    d ↔ -d closest to b (dist=2.7955) ✗
    s ↔ -s closest to mu (dist=0.9623) ✗
    b ↔ -b closest to u (dist=1.8359) ✗

  Per-sector shifted coordinates (relative to sector centroid):

  lepton (centroid: [ 0.16666667 -3.66666667  2.        ]):
    e: (-0.1667, +3.6667, -2.0000)  |w|=4.1800
    mu: (-0.6667, -0.3333, +1.0000)  |w|=1.2472
    tau: (+0.8333, -3.3333, +1.0000)  |w|=3.5785
    gen1 + gen3 - 2·gen2 deviation: 1.247219

  up (centroid: [ 1.66666667 -6.66666667  2.        ]):
    u: (-3.1667, +0.6667, -3.0000)  |w|=4.4127
    c: (-1.1667, -0.3333, +1.0000)  |w|=1.5723
    t: (+4.3333, -0.3333, +2.0000)  |w|=4.7842
    gen1 + gen3 - 2·gen2 deviation: 1.572330

  down (centroid: [ 0.83333333 -7.          0.66666667]):
    d: (-1.3333, -1.0000, -2.6667)  |w|=3.1447
    s: (+0.6667, +0.0000, -0.6667)  |w|=0.9428
    b: (+0.6667, +1.0000, +3.3333)  |w|=3.5434
    gen1 + gen3 - 2·gen2 deviation: 0.942809

========================================================================
  §4. SECTOR-WISE REPRESENTATION DECOMPOSITION
========================================================================

  Each sector has 3 particles (generations 1, 2, 3). In SU(2)
  language, 3 states = spin-1 (triplet). Test if the generational
  structure is consistent with SU(2)_generation.


  lepton: gen1→gen3 direction = [ 1. -7.  3.], |d| = 7.6811
    Projections: ['0.0000', '4.7519', '7.6811']
    Expected (equal): ['0.0000', '3.8406', '7.6811']
    Gen2 midpoint deviation: 0.911322
    → NOT an SU(2) triplet ✗
    Gen2 perpendicular component: [-1.11864407  0.33050847  1.1440678 ] (|⊥| = 1.633858)

  up: gen1→gen3 direction = [ 7.5 -1.   5. ], |d| = 9.0692
    Projections: ['0.0000', '3.9695', '9.0692']
    Expected (equal): ['0.0000', '4.5346', '9.0692']
    Gen2 midpoint deviation: 0.565101
    → NOT an SU(2) triplet ✗
    Gen2 perpendicular component: [-1.28267477 -0.56231003  1.81155015] (|⊥| = 2.289795)

  down: gen1→gen3 direction = [2. 2. 6.], |d| = 6.6332
    Projections: ['0.0000', '2.7136', '6.6332']
    Expected (equal): ['0.0000', '3.3166', '6.6332']
    Gen2 midpoint deviation: 0.603023
    → NOT an SU(2) triplet ✗
    Gen2 perpendicular component: [ 1.18181818  0.18181818 -0.45454545] (|⊥| = 1.279204)

========================================================================
  §5. BOSONIC EXTENSION
========================================================================

  The 3 bosons (W, Z, H) are connected to the fermion lattice via
  the spanning tree. Test their position relative to the fermion
  weight structure.

  Boson positions relative to fermion centroid [ 0.88888889 -5.77777778  1.55555556]:
    W: pqr=[  5.5 -10.    2. ], shifted=(+4.6111, -4.2222, +0.4444)
    Z: pqr=[  8. -11.   0.], shifted=(+7.1111, -5.2222, -1.5556)
    H: pqr=[ 7. -9.  2.], shifted=(+6.1111, -3.2222, +0.4444)

  Convex hull test (are bosons inside the fermion polytope?):
    Fermion convex hull volume: 58.7500
    W: outside ✗ (ΔV = 20.333333)
    Z: outside ✗ (ΔV = 56.666667)
    H: outside ✗ (ΔV = 28.333333)

  Nearest fermion to each boson:
    W ↔ t: dist=3.6401, Δ=(+0.5, +3, +2)
    Z ↔ t: dist=6.0000, Δ=(-2.0, +4, +4)
    H ↔ t: dist=3.0000, Δ=(-1.0, +2, +2)

  Boson chains (from spanning tree):
    t→W: (-0.5, -3, -2)  |Δ|=3.6401
    W→H: (+1.5, +1, +0)  |Δ|=1.8028
    H→Z: (+1.0, -2, -2)  |Δ|=3.0000

  Full 12-point coordinate matrix analysis:
    Rank: 3
    Singular values: [28.18920656  8.34940688  6.73840029]
    Principal directions:
      PC1: (-0.4062, +0.8948, -0.1851)  explains 87.3%
      PC2: (-0.8771, -0.4386, -0.1955)  explains 7.7%
      PC3: (+0.2561, -0.0829, -0.9631)  explains 5.0%

========================================================================
  SUMMARY
========================================================================

  §1: Root system reconstruction
      6 generational roots extracted (2 per sector)
      Cross-sector roots also tabulated

  §2: Dynkin embedding
      Cartan-like matrix computed per sector
      Integer Cartan condition tested

  §3: Weyl symmetry
      Fermion centroid: [ 0.88888889 -5.77777778  1.55555556]
      Reflection matches evaluated

  §4: SU(2)_gen decomposition
      Tested equally-spaced triplet hypothesis per sector

  §5: Bosonic extension
      Spanning tree chain t→W→H→Z analyzed
      Boson positions relative to fermion convex hull determined
      12-point coordinate matrix rank: 3

