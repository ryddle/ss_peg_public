========================================================================
  BOSONIC SECTOR SIGNIFICANCE TEST
========================================================================

§1. 12-PARTICLE LATTICE RECONSTRUCTION
--------------------------------------------------

  Building lattice from spanning tree + b/τ bridge
  (same procedure as consistent_lattice.py)

  Particle      p      q      r    m_pred(MeV)   m_exp(MeV)   Error%  Half-int?
  ────────────────────────────────────────────────────────────────────────
  e          +0.0     +0     +0           0.51         0.51   0.000%       ✓
  mu         -0.5     -4     +3         105.27       105.66   0.368%       ✓
  tau        +1.0     -7     +3        1772.66      1776.86   0.236%       ✓
  u          -1.5     -6     -1           2.13         2.16   1.178%       ✓
  c          +0.5     -7     +3        1253.46      1270.00   1.302%       ✓
  t          +6.0     -7     +4      170175.50    172500.00   1.348%       ✓
  d          -0.5     -8     -2           4.67         4.70   0.543%       ✓
  s          +1.5     -7     +0          92.85        93.40   0.590%       ✓
  b          +1.5     -6     +4        4149.57      4180.00   0.728%       ✓
  W          +5.5    -10     +2       79600.56     80379.00   0.968%       ✓  BOSON
  Z          +8.0    -11     +0       90679.16     91187.60   0.558%       ✓  BOSON
  H          +7.0     -9     +2      124223.07    125100.00   0.701%       ✓  BOSON

  Fermion mean |err|: 0.699%   max: 1.348%
  Boson mean |err|:   0.742%   max: 0.968%
  Full lattice mean:  0.710%


§2. BOSON COORDINATES — Half-integer quality
--------------------------------------------------

  W: (p, q, r) = (+5.50, -10.00, +2.00)
      p = +5.5000  →  nearest ½-int: +5.5  deviation: 0.0000  ✓
      q = -10.0000  →  nearest ½-int: -10.0  deviation: 0.0000  ✓
      r = +2.0000  →  nearest ½-int: +2.0  deviation: 0.0000  ✓

  Z: (p, q, r) = (+8.00, -11.00, +0.00)
      p = +8.0000  →  nearest ½-int: +8.0  deviation: 0.0000  ✓
      q = -11.0000  →  nearest ½-int: -11.0  deviation: 0.0000  ✓
      r = +0.0000  →  nearest ½-int: +0.0  deviation: 0.0000  ✓

  H: (p, q, r) = (+7.00, -9.00, +2.00)
      p = +7.0000  →  nearest ½-int: +7.0  deviation: 0.0000  ✓
      q = -9.0000  →  nearest ½-int: -9.0  deviation: 0.0000  ✓
      r = +2.0000  →  nearest ½-int: +2.0  deviation: 0.0000  ✓

  NOTE: Boson coordinates arise from the sum of integer/half-integer
  ratio vectors, so they are ALGEBRAICALLY half-integer by construction
  of the spanning tree. The physical question is: are the ratio vectors
  themselves accurate? → See §3.


§3. BOSON MASS RATIOS — Lattice vs experiment
--------------------------------------------------

  Ratio                 Exp    Lattice   Error%                    5D vector
  ──────────────────────────────────────────────────────────────────────────
  m_t/m_W            2.1461     2.1379   0.383%              (+0,+3,-1,+2,+0)
  m_H/m_W            1.5564     1.5606   0.270%              (+0,+1,-1,+0,+1)
  m_H/m_Z            1.3719     1.3699   0.144%              (+1,+2,+0,+2,+0)
  m_H/m_t            0.7252     0.7300   0.655%              (+0,-2,+0,-2,+1)
  m_W/m_Z            0.8815     0.8778   0.413%              (+1,+1,+1,+2,-1)
  m_Z/m_t            0.5286     0.5329   0.801%              (-1,-4,+0,-4,+1)
  m_W/m_b           19.2294    19.1828   0.242%              (+0,-4,+0,-2,+4)
  m_Z/m_b           21.8152    21.8527   0.172%              (-1,-5,-1,-4,+5)
  m_W/m_tau         45.2365    44.9045   0.734%              (+0,-3,-1,-1,+4)

  Mean |err| (boson ratios): 0.424%
  Max |err|:                 0.801%


§4. BOSON SECTOR GEOMETRY
--------------------------------------------------

  W: (+5.5, -10.0, +2.0)
  Z: (+8.0, -11.0, +0.0)
  H: (+7.0, -9.0, +2.0)

  Z − W: (+2.5, -1.0, -2.0)
  H − W: (+1.5, +1.0, +0.0)
  H − Z: (-1.0, +2.0, +2.0)

  Boson positions relative to gen-3 fermions:
    W − t: (-0.5, -3.0, -2.0)  |Δ| = 3.64
    W − b: (+4.0, -4.0, -2.0)  |Δ| = 6.00
    W − τ: (+4.5, -3.0, -1.0)  |Δ| = 5.50

    Z − t: (+2.0, -4.0, -4.0)  |Δ| = 6.00
    Z − b: (+6.5, -5.0, -4.0)  |Δ| = 9.12
    Z − τ: (+7.0, -4.0, -3.0)  |Δ| = 8.60

    H − t: (+1.0, -2.0, -2.0)  |Δ| = 3.00
    H − b: (+5.5, -3.0, -2.0)  |Δ| = 6.58
    H − τ: (+6.0, -2.0, -1.0)  |Δ| = 6.40

  Boson triangle area: 5.39
  Normal to boson plane: (+0.371, -0.557, +0.743)


§5. LATTICE DENSITY & PROPER SIGNIFICANCE TEST
--------------------------------------------------

  The 5D formula space collapses to 3D: ln(f) = p·ln2 + q·lnR₃ + r·ln3
  With |exponents| ≤ 4, the effective (p,q,r) triples give ~3,300
  distinct formula values. The significance test must account for
  this lattice density.

  5D formulas: 9^5 = 59,049
  Unique 3D values: 3,323

  Random ratio targets: N = 500,000
  Range: [0.50, 200.0]
  Mean best-match error:   0.1596%
  Median best-match error: 0.1392%

  Per-ratio p-values (boson sector):
  Name             Target   Error%          p_i
  ──────────────────────────────────────────────
  m_t/m_W          2.1461   0.383%  0.967462
  m_H/m_W          1.5564   0.270%  0.805834
  m_H/m_Z          1.3719   0.144%  0.516462
  m_H/m_t          0.7252   0.655%  1.000000

  Joint boson probability (4 ratios):
    log₁₀(P) = -0.40
    P ≈ 10^-0.4

  Per-ratio p-values (ALL 14 ratios):
  Name             Target   Error%          p_i   Sector
  ────────────────────────────────────────────────────────
  m_mu/m_e       206.7683   0.368%  0.951880   fermion
  m_tau/m_mu      16.8170   0.132%  0.476240   fermion
  m_tau/m_e     3477.2283   0.542%  0.994974   fermion
  m_c/m_u        587.9630   0.126%  0.455918   fermion
  m_t/m_c        135.8268   0.046%  0.192996   fermion
  m_s/m_d         19.8723   0.047%  0.199074   fermion
  m_b/m_s         44.7537   0.139%  0.498762   fermion
  m_t/m_b         41.2679   0.151%  0.538382   fermion
  m_d/m_u          2.1759   0.643%  1.000000   fermion
  m_b/m_tau        2.3525   0.493%  0.988146   fermion
  m_t/m_W          2.1461   0.383%  0.967462     boson
  m_H/m_W          1.5564   0.270%  0.805834     boson
  m_H/m_Z          1.3719   0.144%  0.516462     boson
  m_H/m_t          0.7252   0.655%  1.000000     boson

  JOINT PROBABILITIES:
    Fermion only (10 ratios): log₁₀(P) = -2.68  →  P ≈ 10^-2.7
    All 14 ratios:            log₁₀(P) = -3.07  →  P ≈ 10^-3.1
    Bosons add:               0.4 decades of evidence
    Significance: 3.1σ



§6. BUILDING-BLOCK RANDOMIZATION — Are the 5 invariants special?
--------------------------------------------------

  Ultimate test: replace R₂, R₃, R_EE, h∨₃, h∨₂ with 5 random
  numbers in comparable ranges. How often do ALL 14 mass ratios
  (10 fermion + 4 boson) get matched within 1%?

  Real building blocks: mean err = 0.296%, max = 0.655%

  Random building blocks test (10,000,000 trials):
    P(max err ≤ 0.655%) = 0.000000  (0/10,000,000)
    P(mean err ≤ 0.296%) = 0.000000  (0/10,000,000)
    P < 1.0e-07  →  significance > 5.2σ
    P < 1.0e-07  →  significance > 5.2σ

  Comparison — building-block randomization:
    Fermion only (10 ratios, max ≤ 0.643%): P = 0.000000  (0/10,000,000)
    All 14 ratios (max ≤ 0.655%):            P = 0.000000  (0/10,000,000)


========================================================================
  §7. SUMMARY — BOSONIC SECTOR ON THE LATTICE
========================================================================

  The W, Z, and H bosons are placed on the same 3D lattice
  {ln2, lnR₃, ln3} as all 9 fermions, via a spanning tree of
  integer-exponent ratio formulas.

  KEY RESULTS:
  ┌─────────────────────────────────────────────────────────┐
  │  All 3 boson coordinates are (half-)integer                 │
  │  Boson mass accuracy: mean 0.742%, max 0.968%          │
  │  4 boson ratio formulas: all within 1%                      │
  │  12-particle lattice: mean 0.710%, max 1.348%          │
  │  Per-ratio p-value test: P ≈ 10^-3 (all 14 ratios)    │
  │  Building-block randomization: see §6 for P-values          │
  └─────────────────────────────────────────────────────────────┘

  SIGNIFICANCE:
    The per-ratio density test (§5) correctly accounts for lattice
    density: ~3,323 distinct formula values for ~59,049 nominal 5D vectors.
    The building-block randomization (§6) is the definitive test:
    it asks whether ANY 5 random numbers could reproduce 14 mass
    ratios with the SAME exponent vectors at comparable accuracy.

  CONCLUSION:
    The bosonic sector is FULLY CONSISTENT with the fermion lattice.
    No new building blocks needed. Same 5 graph invariants, same
    3D lattice structure, same sub-1% accuracy.

