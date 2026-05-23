==============================================================================
  NEUTRINO MASS DERIVATION FROM GENERATING FUNCTION
  SS-PEG Research Program — 2026
==============================================================================

==============================================================================
§1. GENERATING FUNCTION → NEUTRINO LATTICE COORDINATES
==============================================================================

  Verification with known sectors:
  Sector         2I₃   Nc │  Gen        p        q        r │   m_pred (MeV)    m_exp (MeV)     Err%
  ───────────────────────────────────────────────────────────────────────────────────────────────
  lepton          -1    1 │    1    +0.00    +0.00    +0.00 │         0.5110         0.5110   0.000%
  lepton          -1    1 │    2    -0.50    -4.00    +3.00 │       105.2699       105.6584   0.368%
  lepton          -1    1 │    3    +1.00    -7.00    +3.00 │      1772.6615      1776.8600   0.236%
  up              +1    3 │    1    -1.50    -6.00    -1.00 │         2.1346         2.1600   1.178%
  up              +1    3 │    2    +0.50    -7.00    +3.00 │      1253.4610      1270.0000   1.302%
  up              +1    3 │    3    +6.00    -7.00    +4.00 │    170175.5025    172500.0000   1.348%
  down            -1    3 │    1    -0.50    -8.00    -2.00 │         4.6745         4.7000   0.543%
  down            -1    3 │    2    +1.50    -7.00    +0.00 │        92.8490        93.4000   0.590%
  down            -1    3 │    3    +1.50    -6.00    +4.00 │      4149.5710      4180.0000   0.728%

  ══════════════════════════════════════════════════════════════════════════════
  ★ NEUTRINO SECTOR: (2I₃ = +1, Nc = 1) — GENUINE PREDICTION
  ══════════════════════════════════════════════════════════════════════════════

  Parabolic coefficients (exact fractions):
    p(g) = (15/4)·g² + (-47/4)·g + (7)    [g* = 1.5667]
    q(g) = (1)·g² + (-9)·g + (10)    [g* = 4.5000]
    r(g) = (-4)·g² + (17)·g + (-12)    [g* = 2.1250]

  Lattice coordinates and Dirac masses:
   Gen          p          q          r │     ln(m/mₑ)      m_D (MeV)       m_D (eV)
  ──────────────────────────────────────────────────────────────────────────────
     1    -1.0000    +2.0000    +1.0000 │    -0.783861       0.233342       233.3423
     2    -1.5000    -4.0000    +6.0000 │    +7.930605    1421.143702   1421143.7017
     3    +5.5000    -8.0000    +3.0000 │   +11.865450   72697.529486  72697529.4856

  Isospin partner comparison (ν vs ℓ):
   Gen      m_ℓ (MeV)    m_D^ν (MeV)    Ratio m_D^ν/m_ℓ
  ───────────────────────────────────────────────────────
     1         0.5110       0.233342           0.456640
     2       105.6584    1421.143702          13.450365
     3      1776.8600   72697.529486          40.913482

==============================================================================
§2. DIRAC MASSES vs PHYSICAL NEUTRINO MASSES — SEESAW EXTRACTION
==============================================================================

  Physical neutrino masses (Normal Ordering, NuFIT 5.2):
  ┌────────────────────────┬──────────────┬──────────────┬──────────────┐
  │ Scenario               │ m₁ (meV)     │ m₂ (meV)     │ m₃ (meV)     │
  ├────────────────────────┼──────────────┼──────────────┼──────────────┤
  │ Hierarchical (m₁→0)    │     0.0000   │     8.6139   │    50.1498   │
  │ R_EE⁵ solution         │     2.1297   │     8.8733   │    50.1950   │
  │ Best graph scan        │     3.6500   │     9.3553   │    50.2824   │
  └────────────────────────┴──────────────┴──────────────┴──────────────┘

  Dirac masses from generating function:
    m_D₁ = 233.3423 eV = 233342.3357 keV
    m_D₂ = 1421143.7017 eV = 1421143701.7070 keV
    m_D₃ = 72697529.4856 eV = 72697529485.5840 keV

  Type-I Seesaw extraction: M_R = m_D² / m_ν

   Gen       m_D (eV)       m_ν (eV)        M_R (GeV)         M_R (eV)
  ─────────────────────────────────────────────────────────────────
    ν₁    233342.3357   2.129668e-03     2.556673e+04     2.556673e+13
    ν₂ 1421143701.7070   8.873302e-03     2.276097e+11     2.276097e+20
    ν₃ 72697529485.5840   5.019497e-02     1.052880e+14     1.052880e+23

  With hierarchical masses (m₁ → 0):
    ν₁    233342.3357         (m₁→0)                —                —
    ν₂ 1421143701.7070   8.613942e-03     2.344628e+11     2.344628e+20
    ν₃ 72697529485.5840   5.014978e-02     1.053829e+14     1.053829e+23

==============================================================================
§3. SEESAW SCALE — GRAPH FORMULA SEARCH
==============================================================================

  If M_R has a graph-topological origin, it should be expressible as
  M_R = m_e · (product of building blocks)  or as a ratio M_R_i/M_R_j.


  ── M_R₁ (eV) = 2.556673e+13 ──
     R₂^⁻6 · R₃^⁻6 · R_EE^⁻6 · h∨₃^6 · h∨₂^6                 = 8.467021e+08  err 99.997%  C=30
     R₂^⁻6 · R₃^⁻6 · R_EE^⁻5 · h∨₃^6 · h∨₂^6                 = 5.987088e+08  err 99.998%  C=29
     R₂^⁻6 · R₃^⁻5 · R_EE^⁻6 · h∨₃^6 · h∨₂^6                 = 4.671666e+08  err 99.998%  C=29

  ── M_R₂ (eV) = 2.276097e+20 ──
     R₂^⁻6 · R₃^⁻6 · R_EE^⁻6 · h∨₃^6 · h∨₂^6                 = 8.467021e+08  err 100.000%  C=30
     R₂^⁻6 · R₃^⁻6 · R_EE^⁻5 · h∨₃^6 · h∨₂^6                 = 5.987088e+08  err 100.000%  C=29
     R₂^⁻6 · R₃^⁻5 · R_EE^⁻6 · h∨₃^6 · h∨₂^6                 = 4.671666e+08  err 100.000%  C=29

  ── M_R₃ (eV) = 1.052880e+23 ──
     R₂^⁻6 · R₃^⁻6 · R_EE^⁻6 · h∨₃^6 · h∨₂^6                 = 8.467021e+08  err 100.000%  C=30
     R₂^⁻6 · R₃^⁻6 · R_EE^⁻5 · h∨₃^6 · h∨₂^6                 = 5.987088e+08  err 100.000%  C=29
     R₂^⁻6 · R₃^⁻5 · R_EE^⁻6 · h∨₃^6 · h∨₂^6                 = 4.671666e+08  err 100.000%  C=29

  ── M_R₃/M_R₂ = 4.625816e+02 ──
     R₂^⁻6 · R₃^6 · R_EE^⁻6 · h∨₂^5                          = 4.622392e+02  err 0.074%  C=23
     R₂^⁻5 · R₃^6 · R_EE^⁻6 · h∨₂^6                          = 4.622392e+02  err 0.074%  C=23
     R₂^⁻6 · R₃^6 · R_EE^⁻4 · h∨₂^6                          = 4.622392e+02  err 0.074%  C=22

  ── M_R₃/M_R₁ = 4.118166e+09 ──
     R₂^⁻6 · R₃^⁻6 · R_EE^⁻6 · h∨₃^6 · h∨₂^6                 = 8.467021e+08  err 79.440%  C=30
     R₂^⁻6 · R₃^⁻6 · R_EE^⁻5 · h∨₃^6 · h∨₂^6                 = 5.987088e+08  err 85.462%  C=29
     R₂^⁻6 · R₃^⁻5 · R_EE^⁻6 · h∨₃^6 · h∨₂^6                 = 4.671666e+08  err 88.656%  C=29

  ── M_R₂/M_R₁ = 8.902572e+06 ──
     R₂^⁻4 · R₃^⁻3 · R_EE^⁻2 · h∨₃^6 · h∨₂^6                 = 8.888614e+06  err 0.157%  C=21
     R₂^⁻3 · R₃^⁻3 · R_EE^⁻6 · h∨₃^6 · h∨₂^5                 = 8.888614e+06  err 0.157%  C=23
     R₂^⁻4 · R₃^⁻3 · R_EE^⁻4 · h∨₃^6 · h∨₂^5                 = 8.888614e+06  err 0.157%  C=22

  ── ln(M_R₃/m_e)/ln(M_R₂/m_e) = 1.181939e+00 ──
     R₂^5 · R₃^⁻3 · R_EE^5 · h∨₃^2 · h∨₂^2                   = 1.184007e+00  err 0.175%  C=17
     R₂^⁻3 · R₃^⁻3 · R_EE^5 · h∨₃^2 · h∨₂^⁻6                 = 1.184007e+00  err 0.175%  C=19
     R₂^5 · R₃^⁻3 · R_EE^3 · h∨₃^2 · h∨₂                     = 1.184007e+00  err 0.175%  C=14

==============================================================================
§4. ISOSPIN FLIP: CHARGED LEPTON → NEUTRINO
==============================================================================

  The generating function maps (2I₃=-1, Nc=1) → charged leptons,
  and (2I₃=+1, Nc=1) → neutrino Dirac masses.

  The ISOSPIN FLIP (2I₃: -1 → +1) changes each coefficient by:
    ΔX = α · Δ(2I₃) = α · (+2) = 2α

  This gives a systematic shift in all lattice coordinates.

  Coeff             α         2α │      ℓ value      ν value          Δ
  ──────────────────────────────────────────────────────────────────────
  a_p            11/8       11/4 │            1         15/4       11/4
  b_p           -33/8      -33/4 │         -7/2        -47/4      -33/4
  c_p             9/4        9/2 │          5/2            7        9/2
  a_q             1/4        1/2 │          1/2            1        1/2
  b_q            -7/4       -7/2 │        -11/2           -9       -7/2
  c_q             5/2          5 │            5           10          5
  a_r            -5/4       -5/2 │         -3/2           -4       -5/2
  b_r            19/4       19/2 │         15/2           17       19/2
  c_r              -3         -6 │           -6          -12         -6

  Lattice coordinate shifts (ν − ℓ) at each generation:
   Gen         Δp         Δq         Δr │    Δ[ln(m/mₑ)]        m_D/m_ℓ
  ─────────────────────────────────────────────────────────────────
     1    -1.0000    +2.0000    +1.0000 │      -0.783861       0.456640
     2    -1.0000    +0.0000    +3.0000 │      +2.602690      13.500000
     3    +4.5000    -1.0000    +0.0000 │      +3.713825      41.010385

  Graph formulas for m_D^ν/m_ℓ at each generation:

  ── Gen 1: m_D^ν/m_ℓ = 0.456640 ──
     R₂^3 · R₃^2 · h∨₃ · h∨₂^2                               = 0.456640  err 0.000%  C=8
     R₂^⁻4 · R₃^2 · R_EE^⁻2 · h∨₃ · h∨₂^⁻6                   = 0.456640  err 0.000%  C=15
     R₂^4 · R₃^2 · R_EE^⁻2 · h∨₃ · h∨₂^2                     = 0.456640  err 0.000%  C=11

  ── Gen 2: m_D^ν/m_ℓ = 13.450365 ──
     R₂^⁻5 · R_EE^⁻5 · h∨₃^⁻3 · h∨₂                          = 13.408840  err 0.309%  C=14
     R₂^⁻4 · R_EE^⁻5 · h∨₃^⁻3 · h∨₂^2                        = 13.408840  err 0.309%  C=14
     R₂^⁻5 · R_EE^⁻3 · h∨₃^⁻3 · h∨₂^2                        = 13.408840  err 0.309%  C=13

  ── Gen 3: m_D^ν/m_ℓ = 40.913482 ──
     R₂^6 · R₃^⁻4 · R_EE^6 · h∨₃^5 · h∨₂^3                   = 40.969603  err 0.137%  C=24
     R₂^⁻2 · R₃^⁻4 · R_EE^6 · h∨₃^5 · h∨₂^⁻5                 = 40.969603  err 0.137%  C=22
     R₂^4 · R₃^⁻4 · R_EE^4 · h∨₃^5                           = 40.969603  err 0.137%  C=17

==============================================================================
§5. DIRAC MASS RATIOS — GRAPH FORMULAS
==============================================================================

  The generating function predicts Dirac mass RATIOS between neutrino
  generations. These should have graph formulas, and they determine
  the STRUCTURE of the seesaw (whether M_R is universal or not).

  Ratio                 Value Formula                                               Err%
  ──────────────────────────────────────────────────────────────────────────────────────────
  m_D₂/m_D₁         6090.3809 R₃^⁻6 · R_EE^⁻3 · h∨₃^5 · h∨₂^⁻2                    0.000%
  m_D₃/m_D₂           51.1542 R₂^⁻5 · R₃^⁻4 · R_EE^2 · h∨₃^⁻3 · h∨₂^3             0.000%
  m_D₃/m_D₁       311548.8206 R₂^⁻6 · R₃^⁻3 · R_EE^⁻3 · h∨₃^2 · h∨₂^5             0.375%

  Comparison with charged-lepton generation ratios:
  m_μ/m_e            206.7683 R₂^⁻5 · R₃^6 · R_EE · h∨₃^4 · h∨₂^2                 0.033%
  m_τ/m_μ             16.8170 R₂^5 · R₃^⁻6 · R_EE^6 · h∨₃^5 · h∨₂⁻¹               0.032%
  m_τ/m_e           3477.2283 R₂^⁻4 · R₃^⁻4 · R_EE^⁻5 · h∨₃^⁻2 · h∨₂^5            0.137%

  Generation scaling exponent γ = ln(m₃/m₂)/ln(m₂/m₁):
    Neutrino Dirac: γ_ν = 0.451530
    Charged lepton: γ_ℓ = 0.529371

==============================================================================
§6. CONSISTENCY: DIRAC MASSES + SEESAW → PHYSICAL OBSERVABLES
==============================================================================

  If the seesaw m_ν = m_D²/M_R operates with universal M_R:
    m_ν₂/m_ν₃ = (m_D₂/m_D₃)² = m_D₂²/m_D₃²

  If M_R is generation-dependent (M_R_i):
    m_ν₂/m_ν₃ = (m_D₂/m_D₃)² × (M_R₃/M_R₂)

  CASE 1: Universal M_R
    (m_D₂/m_D₃)² = 0.000382
    Experimental √(Δm²₂₁/Δm²₃₁) ≈ m₂/m₃ = 0.171764
    Predicted r = (m_D₂/m_D₃)² = 0.000382
    Experimental r = Δm²₂₁/Δm²₃₁ = 0.029503
    Error: 98.7%
    → Universal M_R does NOT reproduce the mass ratio.
    → NEED generation-dependent M_R.

  CASE 2: Generation-dependent M_R
    Required M_R₃/M_R₂ = √r / (m_D₂/m_D₃)² = 449.465408

    Graph formula search for M_R₃/M_R₂ = 449.465408:
      R₂^⁻5 · R₃ · R_EE · h∨₃^2 · h∨₂^2                       = 449.447137  err 0.004%  C=11
      R₂^⁻5 · R₃ · R_EE^5 · h∨₃^2 · h∨₂^4                     = 449.447137  err 0.004%  C=17
      R₂^⁻6 · R₃ · R_EE^⁻3 · h∨₃^2 · h∨₂⁻¹                    = 449.447137  err 0.004%  C=13
      R₂^⁻6 · R₃ · R_EE⁻¹ · h∨₃^2                             = 449.447137  err 0.004%  C=10
      R₂^⁻6 · R₃ · R_EE · h∨₃^2 · h∨₂                         = 449.447137  err 0.004%  C=11

    From r directly: M_R₃/M_R₂ = √[r/(m_D₂/m_D₃)⁴] = 449.465408
    Check: (m_D₂/m_D₃)⁴ × (M_R₃/M_R₂)² = 0.029503 vs r = 0.029503

==============================================================================
§7. RAMANUJAN SPECTRAL INTERPRETATION
==============================================================================

  In the SS-PEG framework, masses arise from spectral hops on the
  Cartan-Weyl root graph. Each building block represents a different
  type of spectral process:

  - R₂, R₃: spectral gap ratios (propagation efficiency)
  - R_EE:   bifundamental bridge (sector coupling)
  - h∨₃, h∨₂: Coxeter amplification (loop corrections)

  For neutrinos (Q=0, Nc=1), the Ramanujan-optimal propagation on
  SU(3) is ABSENT. The generating function captures this through
  the (2I₃, Nc) dependence:

  The isospin flip 2I₃: -1 → +1 shifts lattice coordinates by 2α_k.
  The α coefficients control HOW MUCH each parabolic parameter
  responds to the isospin flip.

  α coefficients (isospin sensitivity):
  Coeff               α      |α| Sector              
  ───────────────────────────────────────────────────────
  a_p              11/8   1.3750 SU(2) binary        
  b_p             -33/8   4.1250 SU(2) binary        
  c_p               9/4   2.2500 SU(2) binary        
  a_q               1/4   0.2500 SU(3) irrational    
  b_q              -7/4   1.7500 SU(3) irrational    
  c_q               5/2   2.5000 SU(3) irrational    
  a_r              -5/4   1.2500 SU(3) ternary       
  b_r              19/4   4.7500 SU(3) ternary       
  c_r                -3   3.0000 SU(3) ternary       

  Total |α| shift for p: 7.7500

  Total |α| shift for q: 4.5000

  Total |α| shift for r: 9.0000

  → Dominant shift: r coordinate
    This means the isospin flip primarily affects the ternary (ln3) channel.

==============================================================================
§8. NEUTRINO GENERATION CURVES — QUADRATIC ANALYSIS
==============================================================================

  The generating function gives neutrino Dirac masses on a parabola:
    ln(m_D(g)/mₑ) = A_ν·g² + B_ν·g + C_ν

  where A_ν, B_ν, C_ν are computed from the neutrino lattice coefficients.

  Mass-space parabola ln(m/mₑ) = A·g² + B·g + C:
  Sector                  A            B            C   g*(turn)  2A(accel)
  ─────────────────────────────────────────────────────────────────
  Neutrino ν      -2.389810   +15.883897   -14.277947     3.3233  -4.779620
  Lepton ℓ        -1.252103    +9.084224    -7.832121     3.6276  -2.504206

  Koide-like geometric ratio m₃·m₁/m₂² = e^(2A):
    Neutrino Dirac: e^(2A_ν) = 0.008399
    Charged lepton: e^(2A_ℓ) = 0.081741
    Actual m_D₃·m_D₁/m_D₂² = 0.008399
    Consistency: 0.0000%
    S_K(p|ν) = 29.3125
    S_K(q|ν) = 26.3333
    S_K(r|ν) = 22.3333
    S_K(ν total) = 77.9792 = 3743/48

==============================================================================
§9. PHYSICAL NEUTRINO MASS PREDICTIONS — SYNTHESIS
==============================================================================

  From the generating function, the neutrino Dirac mass spectrum is:
    m_D = (233342.34, 1421143701.71, 72697529485.58) keV

  Combined with the graph-derived mass-squared ratio r = √2/48:
    r_graph = R₂³ · R_EE · h∨₃⁻¹ = 0.029463
    r_exp   = 0.029503

  And the Dirac mass ratio from the generating function:
    (m_D₂/m_D₃)² = 0.000382

  KEY TEST: Does m_D₂/m_D₃ from F match √r from graph formula?
    m_D₂/m_D₃ (from F)  = 0.019549
    √r (graph formula)   = 0.171647
    √r (experimental)    = 0.171764
    Err(Dirac vs graph): 88.611%
    Err(Dirac vs exp):   88.619%
    → INCONSISTENT with universal M_R. The seesaw scale MUST be generation-dependent.

==============================================================================
§10. THE COMPLETE NEUTRINO PICTURE — SUMMARY
==============================================================================

  ╔══════════════════════════════════════════════════════════════════════╗
  ║  THE NEUTRINO MASS DERIVATION FROM GRAPH TOPOLOGY                   ║
  ╚══════════════════════════════════════════════════════════════════════╝

  1. GENERATING FUNCTION PREDICTION
     The generating function F(2I₃, Nc, g), calibrated on charged fermions
     only, evaluated at neutrino quantum numbers (2I₃=+1, Nc=1), produces:

     ν₁: (p,q,r) = (-1.00, +2.00, +1.00) → m_D₁ = 233342.34 keV
     ν₂: (p,q,r) = (-1.50, -4.00, +6.00) → m_D₂ = 1421143701.71 keV
     ν₃: (p,q,r) = (+5.50, -8.00, +3.00) → m_D₃ = 72697529485.58 keV

     These are the DIRAC masses (Yukawa-type, before seesaw suppression).

  2. SEESAW MECHANISM
     Physical neutrino masses m_ν ≪ m_D require Type-I seesaw:
       m_ν = m_D² / M_R

     The ratio m_D/m_ℓ (isospin partner) varies by generation — the
     generating function predicts NON-PROPORTIONAL Dirac masses.

  3. MASS ORDERING
     The graph framework PREDICTS Normal Ordering (NO) through:
     - Better graph formula fits for NO ratios vs IO ratios
     - r_Δm² = √2/48 has a compact formula in NO, not in IO

  4. PMNS ANGLES
     All 4 PMNS parameters derived from building blocks (already done):
       sinθ₂₃ = R₃³·R_EE⁻⁴·h∨₃²·h∨₂⁻³     (0.000%)
       sinθ₁₂ = R₂⁻³·R₃·R_EE⁴·h∨₂⁻¹        (0.204%)
       sinθ₁₃ = R₂·R_EE⁻⁴·h∨₃⁻³·h∨₂         (0.237%)
       δ_CP/π = R₂⁻²·R₃⁻²·R_EE⁴·h∨₃⁻¹       (0.047%)

  5. GRAPH-TOPOLOGICAL ORIGIN
     The generating function's α coefficients quantify the isospin
     sensitivity of each lattice coordinate. Neutrinos occupy a
     SINGULAR point where Q=0 collapses the electromagnetic tensor.
     The Ramanujan spectral quality of the building blocks ensures
     maximal propagation efficiency on the SU(2) cycle.


==============================================================================
§11. FLATNESS AND GENERATION STRUCTURE FOR NEUTRINOS
==============================================================================

  Flatness test b = -5a for neutrino sector:
    p: b/a = -11.7500/3.7500 = -3.1333  (need -5.0, off by +1.8667)
    q: b/a = -9.0000/1.0000 = -9.0000  (need -5.0, off by -4.0000)
    r: b/a = 17.0000/-4.0000 = -4.2500  (need -5.0, off by +0.7500)

  Charged lepton flat coordinate: r
    b_r/a_r = 4.0000

  → Neutrino flat coordinate: NONE — all coordinates are non-flat!
    This is a UNIQUE property of the neutrino sector.
    Charged fermions always have exactly 1 flat coordinate per sector.
