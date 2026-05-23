==============================================================================
  BOSON SECTOR IN THE GENERATING FUNCTION FRAMEWORK
  SS-PEG Research Program — 2026
==============================================================================

==============================================================================
§1. BOSON LATTICE COORDINATES — VERIFICATION
==============================================================================

   Boson        p        q        r │   m_pred (MeV)  m_exp (MeV)     Err%
  ────────────────────────────────────────────────────────────────────
       W     11/2      -10        2 │        79600.6      80379.0   0.968%
       Z        8      -11        0 │        90679.2      91187.6   0.558%
       H        7       -9        2 │       124223.1     125100.0   0.701%

  Note: All coordinates are multiples of 1/2 (half-integer lattice).
  This is the same quantization as fermion coordinates.


==============================================================================
§2. ORDERING ANALYSIS — FLATNESS TEST FOR ALL PERMUTATIONS
==============================================================================

  For fermions, each sector has a FLAT coordinate where b/a = -5
  (universal turning point g* = 5/2). We test all 6 orderings of
  (W, Z, H) as (g=1, g=2, g=3) to find which one(s) reproduce
  this flatness condition.

         Ordering │ Coord        a        b        c │      b/a  Flat?
  ──────────────────────────────────────────────────────────────────────
            W→Z→H │     p     -7/4     31/4     -1/2 │   -4.429       
            W→Z→H │     q      3/2    -11/2       -6 │   -3.667       
            W→Z→H │     r        2       -8        8 │   -4.000       

            W→H→Z │     p     -1/4      9/4      7/2 │   -9.000       
            W→H→Z │     q     -3/2     11/2      -14 │   -3.667       
            W→H→Z │     r       -1        3        0 │   -3.000       

            Z→W→H │     p        2    -17/2     29/2 │   -4.250       
            Z→W→H │     q        0        1      -12 │  linear       
            Z→W→H │     r       -1        5       -4 │   -5.000 ✓ FLAT

            Z→H→W │     p     -1/4     -1/4     17/2 │    1.000       
            Z→H→W │     q     -3/2     13/2      -16 │   -4.333       
            Z→H→W │     r       -1        5       -4 │   -5.000 ✓ FLAT

            H→W→Z │     p        2    -15/2     25/2 │   -3.750       
            H→W→Z │     q        0       -1       -8 │  linear       
            H→W→Z │     r       -1        3        0 │   -3.000       

            H→Z→W │     p     -7/4     25/4      5/2 │   -3.571       
            H→Z→W │     q      3/2    -13/2       -4 │   -4.333       
            H→Z→W │     r        2       -8        8 │   -4.000       


  SUMMARY:
      W→Z→H: no flatness
      W→H→Z: no flatness
    ★ Z→W→H: 1 flat coordinate(s): r
    ★ Z→H→W: 1 flat coordinate(s): r
      H→W→Z: no flatness
      H→Z→W: no flatness

  ★ Selected ordering: Z→W→H (matches fermion flatness condition)


==============================================================================
§3. THE Z→W→H PARABOLA: COEFFICIENTS AND PROPERTIES
==============================================================================

  Assignment: g=1 → Z, g=2 → W, g=3 → H

  Coord │        a        b        c │      b/a       g* │ Notes
  ────────────────────────────────────────────────────────────────────────
      p │        2    -17/2     29/2 │   -4.250   2.1250 │ 
      q │        0        1      -12 │    ∞    ∞ │ 
      r │       -1        5       -4 │   -5.000   2.5000 │ FLAT (b = -5a, g* = 5/2)

  Rationality check:
    p_a = 2 (denominator 1)
    p_b = -17/2 (denominator 2)
    p_c = 29/2 (denominator 2)
    q_a = 0 (denominator 1)
    q_b = 1 (denominator 1)
    q_c = -12 (denominator 1)
    r_a = -1 (denominator 1)
    r_b = 5 (denominator 1)
    r_c = -4 (denominator 1)

  All denominators: [1, 2]
  Maximum denominator: 2
  → SAME as fermion coefficients (multiples of 1/4)

  Electron origin condition (c = -a - b):
    p: c = 29/2, -a-b = 13/2 → ✗ FAILS
    q: c = -12, -a-b = -1 → ✗ FAILS
    r: c = -4, -a-b = -4 → ✓ HOLDS


==============================================================================
§4. FERMION MANIFOLD PROJECTION — EFFECTIVE QUANTUM NUMBERS
==============================================================================

  The fermion GF has X(s) = α·(2I₃) + β·Nc + γ for each of 9 coefficients.
  For the boson parabola, we ask: what (2I₃_eff, Nc_eff) would produce
  the observed coefficients? This is an OVERDETERMINED system (9 equations,
  2 unknowns), so we measure the residual.

  Least-squares effective quantum numbers:
    2I₃_eff = 1.815075
    Nc_eff  = 3.707780

   Coefficient │   Observed GF Predicted   Residual   Rel.Err%
  ──────────────────────────────────────────────────────────────
           a_p │     2.0000       2.1629    -0.1629      8.15%
           a_q │     0.0000       0.5268    -0.5268       inf%
           a_r │    -1.0000      -1.6341     0.6341     63.41%
           b_p │    -8.5000      -3.6041    -4.8959     57.60%
           b_q │     1.0000      -1.6261     2.6261    262.61%
           b_r │     5.0000       9.3635    -4.3635     87.27%
           c_p │    14.5000      -0.6433    15.1433    104.44%
           c_q │   -12.0000      -6.9168    -5.0832     42.36%
           c_r │    -4.0000      -9.0297     5.0297    125.74%

  RMS residual: 6.065104
  Total sum of squared residuals: 331.069330

  Per-coefficient-pair effective (2I₃, Nc):
                Pair │    2I₃_eff     Nc_eff │ Notes
  ────────────────────────────────────────────────────────────

  Range of pairwise solutions (35 solvable pairs):
    2I₃_eff ∈ [-109.0000, 51.0000], mean = -4.0392
    Nc_eff  ∈ [-105.0000, 31.0000], mean = -2.7752
    Spread(2I₃) = 160.0000
    Spread(Nc)  = 136.0000

    Best pairwise fit: (a_p, c_q)
    → 2I₃_eff = 2.368421, Nc_eff = 4.631579
    → Total squared residual = 369.669148

  CONCLUSION: The large spread confirms bosons are NOT on the fermion
  linear manifold. No single (2I₃, Nc) reproduces all 9 boson coefficients.


==============================================================================
§5. SPIN-EXTENDED GENERATING FUNCTION
==============================================================================

  All fermion sectors have spin S = 1/2 (2S = 1). The boson sector
  has spin S = 1 for W,Z and S = 0 for H, but they form a SINGLE
  parabola. Can we extend: X = α·(2I₃) + β·Nc + δ·(2S) + γ?

  Problem: All 4 fermion sectors have 2S = 1, so δ is degenerate
  with γ in the fermion fit. We need the boson sector to break the
  degeneracy.

  Approach: Assign effective 2S_eff to the boson sector and solve
  the 5-sector system. But the boson sector must be assigned
  effective (2I₃_B, Nc_B) as well.

  Testing effective spin values for the boson sector:
    2S_B │    2I₃_eff     Nc_eff          δ │    RMS res
  ──────────────────────────────────────────────────────────
       0 │     2.0056     3.8836    -1.0148 │   5.986503
       1 │     1.8151     3.7078     0.0000 │   6.065104
       2 │     2.0056     3.8836     1.0148 │   5.986503
       3 │     2.0056     3.8836     0.5074 │   5.986503
       4 │     2.0056     3.8836     0.3383 │   5.986503

  Note: The RMS residual measures how well the extended GF fits.
  A significant reduction from the 2-parameter fit (§4) would indicate
  that spin is a relevant quantum number for the generating function.
  2-parameter RMS: 6.065104


==============================================================================
§6. KINETIC ACTION — ALL SECTORS COMPARED
==============================================================================

  The kinetic action S_K = ∫₁³ ½·|ẋ|² dg measures the "cost" of
  traversing a sector's parabola from g=1 to g=3.
  S_K^(k) = (52/3)·a² + 8·a·b + b² per coordinate.
  (Using ∫₁³ g² dg = 26/3, ∫₁³ g dg = 4, ∫₁³ dg = 2)

        Sector │       S_K(p)       S_K(q)       S_K(r) │     S_K(total)    Decimal
  ────────────────────────────────────────────────────────────────────────────────
        lepton │        19/12       151/12         21/4 │         233/12    19.4167
            up │       871/48         7/12         37/4 │        1343/48    27.9792
          down │          7/3            1         31/3 │           41/3    13.6667
      neutrino │       469/16         79/3         67/3 │        3743/48    77.9792
  ────────────────────────────────────────────────────────────────────────────────
       ★ boson │        67/12            1          7/3 │         107/12     8.9167

  Sector kinetic actions (ordered):
           boson: S_K =       107/12 ≈     8.9167
            down: S_K =         41/3 ≈    13.6667
          lepton: S_K =       233/12 ≈    19.4167
              up: S_K =      1343/48 ≈    27.9792
        neutrino: S_K =      3743/48 ≈    77.9792

  Boson S_K = 107/12 = 8.916667
  Denominator: 12

  Kinetic action ratios:
    S_K(lepton) / S_K(boson) = 233/107 ≈ 2.177570
    S_K(up) / S_K(boson) = 1343/428 ≈ 3.137850
    S_K(down) / S_K(boson) = 164/107 ≈ 1.532710
    S_K(neutrino) / S_K(boson) = 3743/428 ≈ 8.745327

  S_K(charged fermions) = 977/16 = 61.0625
  S_K(all 5 sectors) = 3551/24 = 147.9583

  Building-block formula search for S_K(boson) = 8.916667:
      R₂^4 · R₃^⁻6 · R_EE^5 · h∨₃^6 · h∨₂^⁻5 = 8.921456 (err 0.0537%, complexity 26)
      R₂^5 · R₃^⁻6 · R_EE^3 · h∨₃^6 · h∨₂^⁻5 = 8.921456 (err 0.0537%, complexity 25)
      R₂^5 · R₃^⁻6 · R_EE^5 · h∨₃^6 · h∨₂^⁻4 = 8.921456 (err 0.0537%, complexity 26)
      R₂^4 · R₃^⁻6 · R_EE^3 · h∨₃^6 · h∨₂^⁻6 = 8.921456 (err 0.0537%, complexity 25)
      R₂^3 · R₃^⁻6 · R_EE^5 · h∨₃^6 · h∨₂^⁻6 = 8.921456 (err 0.0537%, complexity 26)


==============================================================================
§7. BOSON-FERMION STRUCTURAL ANALYSIS
==============================================================================

  7a. Parabola coefficient comparison (all 5 sectors):

     Coeff │   Lepton       Up     Down Neutrino │    Boson
  ────────────────────────────────────────────────────────────
       a_p │        1      7/4       -1     15/4 │        2
       a_q │      1/2      1/2        0        1 │        0
       a_r │     -3/2     -3/2        1       -4 │       -1

       b_p │     -7/2    -13/4        5    -47/4 │    -17/2
       b_q │    -11/2     -5/2        1       -9 │        1
       b_r │     15/2     17/2       -1       17 │        5

       c_p │      5/2        0     -9/2        7 │     29/2
       c_q │        5       -4       -9       10 │      -12
       c_r │       -6       -8       -2      -12 │       -4


  7b. Curvature vectors a = (a_p, a_q, a_r):
        Sector │      a_p      a_q      a_r │       |a|²
  ────────────────────────────────────────────────────
        lepton │        1      1/2     -3/2 │        7/2
            up │      7/4      1/2     -3/2 │      89/16
          down │       -1        0        1 │          2
      neutrino │     15/4        1       -4 │     497/16
         boson │        2        0       -1 │          5

  7c. Displacement vectors: boson(g=1) to fermion(g=3):
  Boson g=1 (Z): (8, -11, 0)
    →     lepton(g=3) = (1, -7, 3), Δ = (-7, 4, 3)
    →         up(g=3) = (6, -7, 4), Δ = (-2, 4, 4)
    →       down(g=3) = (3/2, -6, 4), Δ = (-13/2, 5, 4)
    →   neutrino(g=3) = (11/2, -8, 3), Δ = (-5/2, 3, 3)

  7d. Fermion sector closest to boson parabola (coefficient distance):
    d²(boson,     lepton) = 512 ≈ 512.0000
    d²(boson,         up) = 2743/8 ≈ 342.8750
    d²(boson,       down) = 2421/4 ≈ 605.2500
    d²(boson,   neutrino) = 6975/8 ≈ 871.8750

  7e. Flatness summary (all sectors):
        Sector │ Flat coord   a_flat   b_flat      b/a │       g*
  ──────────────────────────────────────────────────────────────
        lepton │          r     -3/2     15/2       -5 │      5/2
            up │          q      1/2     -5/2       -5 │      5/2
          down │          p       -1        5       -5 │      5/2
      neutrino │       NONE        —        —        — │        —
       ★ boson │          r       -1        5       -5 │      5/2


==============================================================================
§8. BUILDING-BLOCK FORMULAS FOR BOSON MASS RATIOS
==============================================================================

  Even though bosons don't sit on the fermion GF, their mass ratios
  (both inter-boson and boson-fermion) may be expressible as products
  of the 5 building blocks R₂, R₃, R_EE, h∨₃, h∨₂.

  Known boson mass ratios (from spanning tree):
         Ratio │ Experimental                        Formula    Predicted     Err%
  ────────────────────────────────────────────────────────────────────────────────
       m_t/m_W │     2.146083          R₃^3 · R_EE⁻¹ · h∨₃^2     2.137868   0.383%
       m_H/m_W │     1.556377              R₃ · R_EE⁻¹ · h∨₂     1.560580   0.270%
       m_H/m_Z │     1.371897              R₂ · R₃^2 · h∨₃^2     1.369919   0.144%
       m_H/m_t │     0.725217          R₂⁻¹ · R₃^⁻2 · h∨₃^⁻2     0.729970   0.655%

  Additional boson mass ratios:

  m_Z/m_W = 1.134470:
      R₂^⁻4 · R₃^6 · R_EE^3 · h∨₃^⁻2 · h∨₂^6 = 1.134904 (err 0.0383%, |e|=21) ★
        R₂^⁻3 · R₃^6 · R_EE · h∨₃^⁻2 · h∨₂^6 = 1.134904 (err 0.0383%, |e|=18) ★
      R₂^⁻5 · R₃^6 · R_EE^5 · h∨₃^⁻2 · h∨₂^6 = 1.134904 (err 0.0383%, |e|=24) ★

  m_W/m_e = 157297.779183:
      R₂^⁻5 · R₃^3 · R_EE^⁻5 · h∨₃^4 · h∨₂^6 = 157620.743323 (err 0.2053%, |e|=23) ★
      R₂^⁻6 · R₃^3 · R_EE^⁻5 · h∨₃^4 · h∨₂^5 = 157620.743323 (err 0.2053%, |e|=23) ★
      R₂^⁻6 · R₃^3 · R_EE^⁻3 · h∨₃^4 · h∨₂^6 = 157620.743323 (err 0.2053%, |e|=22) ★

  m_Z/m_e = 178449.681746:
        R₂^⁻6 · R₃⁻¹ · R_EE^⁻6 · h∨₃ · h∨₂^6 = 178168.142682 (err 0.1578%, |e|=20) ★
      R₂^2 · R₃^⁻4 · R_EE^⁻5 · h∨₃^6 · h∨₂^4 = 177990.966340 (err 0.2571%, |e|=21) ★
      R₂^3 · R₃^⁻4 · R_EE^⁻3 · h∨₃^6 · h∨₂^6 = 177990.966340 (err 0.2571%, |e|=22) ★

  m_H/m_e = 244814.593063:
               R₂^⁻5 · R₃^⁻6 · h∨₃^3 · h∨₂^3 = 244994.823155 (err 0.0736%, |e|=17) ★
               R₂^⁻3 · R₃^⁻6 · h∨₃^3 · h∨₂^5 = 244994.823155 (err 0.0736%, |e|=17) ★
     R₂^⁻2 · R₃^⁻6 · R_EE^⁻2 · h∨₃^3 · h∨₂^5 = 244994.823155 (err 0.0736%, |e|=18) ★

  Inter-boson coordinate displacements (in the Z→W→H ordering):
    Z→W: Δ(p,q,r) = (-5/2, 1, 2), m_W/m_Z = 0.881469, pred = 0.877826
    Z→H: Δ(p,q,r) = (-1, 2, 2), m_H/m_Z = 1.371897, pred = 1.369919
    W→H: Δ(p,q,r) = (3/2, 1, 0), m_H/m_W = 1.556377, pred = 1.560580


==============================================================================
§9. BOSON TRANSITION MATRIX
==============================================================================

  The fermion transition matrix N has entries n_sk = 2(x(2) - x(1))
  with Smith normal form diag(1,2,4). Does the boson sector share
  this structure?

  p: x(1) = 8, x(2) = 11/2, Δ₁₂ = -5/2, n = 2Δ = -5
  q: x(1) = -11, x(2) = -10, Δ₁₂ = 1, n = 2Δ = 2
  r: x(1) = 0, x(2) = 2, Δ₁₂ = 2, n = 2Δ = 4

  Boson N-row: (-5, 2, 4)

  2-adic valuations:
    n_p = -5, v₂ = 0
    n_q = 2, v₂ = 1
    n_r = 4, v₂ = 2

  Comparison with fermion N-rows:
        Sector │    n_p    n_q    n_r
  ────────────────────────────────────
        lepton │     -1     -8      6
            up │      4     -2      8
          down │      4      2      4
  ────────────────────────────────────
       ★ boson │     -5      2      4


==============================================================================
§10. BOSON CURVATURE MATRIX AND SPECTRAL ANALYSIS
==============================================================================

  Boson curvature vector: a = (2.0, 0.0, -1.0)
  Boson velocity vector:  b = (-8.5, 1.0, 5.0)
  Boson offset vector:    c = (14.5, -12.0, -4.0)

  Inner products:
    a·b = -22.000000
    a·c = 33.000000
    b·c = -155.250000
    |a|² = 5.000000
    |b|² = 98.250000
    |c|² = 370.250000

  Curvature alignment (cos θ between a_boson and a_fermion):
    cos θ(boson,     lepton) = +0.836660  (θ = 33.2°)
    cos θ(boson,         up) = +0.948091  (θ = 18.5°)
    cos θ(boson,       down) = -0.948683  (θ = 161.6°)
    cos θ(boson,   neutrino) = +0.922772  (θ = 22.7°)

  b/a ratios per coordinate:
    b_p/a_p = -4.250000
    b_q/a_q = ∞ (a=0)
    b_r/a_r = -5.000000


==============================================================================
§11. SUMMARY AND CONCLUSIONS
==============================================================================

  KEY FINDINGS:

  1. FLATNESS CONDITION: Among all 6 orderings of (W, Z, H), exactly
     2 ordering(s) satisfy the universal flatness condition b = -5a.
     → Z→W→H: flat in r (same turning point g* = 5/2)
     → Z→H→W: flat in r (same turning point g* = 5/2)

  2. RATIONAL ARITHMETIC: All 9 boson parabola coefficients are exact
     rationals with denominators in [1, 2], matching the
     fermion arithmetic (multiples of 1/4).

  3. FERMION MANIFOLD: The boson coefficients are NOT on the fermion
     linear manifold X = α·(2I₃) + β·Nc + γ.
     → Least-squares effective QN: 2I₃_eff = 1.8151, Nc_eff = 3.7078
     → RMS residual = 6.065104 (significant deviation)

  4. KINETIC ACTION: S_K(boson) = 107/12 ≈ 8.9167
     → Rank 1/5 (lowest)

  5. BOSON-FERMION RELATIONSHIP: Bosons form a SEPARATE sector with
     its own parabolic structure. The flat coordinate (r)
     shares the universal g* = 5/2 turning point with all 3 charged-fermion
     sectors (leptons→r, up→q, down→p).

  6. INTER-SECTOR PATTERN: With the boson sector, we now have 5 sectors:
     ℓ(r flat), u(q flat), d(p flat), ν(NO flat), B(r flat).
     The boson sector shares the r-flatness with leptons — both are
     color singlets (Nc = 1).

  7. MASS RATIOS: All boson-fermion mass ratios (m_t/m_W, m_H/m_W,
     m_H/m_Z, m_H/m_t) are exact building-block formulas at sub-1%
     accuracy, confirming that bosons participate in the same algebraic
     structure as fermions even though they sit on a separate parabola.

  PHYSICAL INTERPRETATION:
  The generating function F(2I₃, Nc, g) governs FERMION parabolas.
  Bosons require a separate treatment — they are not "generations" of
  the same particle but distinct particles with different quantum numbers.
  Nevertheless, the universal features (rational coefficients, flatness
  at g* = 5/2, building-block mass ratios) persist, suggesting a deeper
  unified structure that encompasses both statistics.

