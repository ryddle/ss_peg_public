========================================================================
  §1. W-H-Z TRIANGLE: RIGHT ANGLE VERIFICATION
========================================================================

  Boson coordinates (p, q, r):
    W = (11/2, -10, 2)
    H = (7, -9, 2)
    Z = (8, -11, 0)

  Side vectors:
    WH = H − W = (+1.5, +1.0, +0.0)
    WZ = Z − W = (+2.5, -1.0, -2.0)
    HZ = Z − H = (+1.0, -2.0, -2.0)

  Side lengths:
    |WH| = 1.802776
    |WZ| = 3.354102
    |HZ| = 3.000000

  Exact squared lengths:
    |WH|² = 13/4 = 3.2500
    |WZ|² = 45/4 = 11.2500
    |HZ|² = 9 = 9.0000

  Pythagorean test (a² + b² = c²?):
    |WH|² + |HZ|² = 49/4 = 12.2500  vs  |WZ|² = 11.2500
    |WH|² + |WZ|² = 29/2 = 14.5000  vs  |HZ|² = 9.0000
    |HZ|² + |WZ|² = 81/4 = 20.2500  vs  |WH|² = 3.2500

  ∠W: dot product = 11/4, cos(θ) = 0.454794, θ = 62.95°

  ∠H: dot product = 1/2, cos(θ) = 0.092450, θ = 84.70°

  ∠Z: dot product = 17/2, cos(θ) = 0.844737, θ = 32.36°

  Sum of angles: 180.00° (should be 180°)

  ► The angle at H is 84.70° — nearly a RIGHT ANGLE!
    Deviation from 90°: 5.30°

  In mass space (coords × basis):
    ∠H (mass space) = 90.25°

========================================================================
  §2. THE 18 FREE COORDINATE PARAMETERS
========================================================================

  Instead of (a, b, c) parametrization, work directly with coordinates.
  Each (p, q, r) ∈ (1/2)Z. Total: 9 fermions × 3 coords = 27.
  
  The 9 constraints:
    C1: Electron origin → (p_e, q_e, r_e) = (0, 0, 0)     [3 eqs]
    C2: Flatness → x_flat(2) = x_flat(3)                    [3 eqs]
        r_μ = r_τ,  q_c = q_t,  p_s = p_b
    C3: q-unification → q_τ = q_t                           [1 eq]
    C4: Curvature universality → a_q^(ℓ) = a_q^(u),
        a_r^(ℓ) = a_r^(u)                                   [2 eqs]
    
  Total: 3 + 3 + 1 + 2 = 9 constraints.
  Free coordinates: 27 − 3 (electron) − 6 (C2+C3+C4) = 18

  Derived (constrained) quantities:
    r_τ = r_μ                          (from C2: lepton r-flat)
    q_c = q_τ                          (from C2+C3: up q-flat + q-unif)
    q_t = q_τ                          (from C3: q-unification)
    p_b = p_s                          (from C2: down p-flat)
    q_u = 2q_τ − 2q_μ                 (from C4: a_q universality)
    r_t = 2r_c − r_u − r_μ            (from C4: a_r universality)

  Verification with physical values:
    r_τ = r_μ: 3 = 3  ✓
    q_c = q_τ: -7 = -7  ✓
    q_t = q_τ: -7 = -7  ✓
    p_b = p_s: 3/2 = 3/2  ✓
    q_u = 2q_τ − 2q_μ: -6 = -6  ✓
    r_t = 2r_c − r_u − r_μ: 4 = 4  ✓
    All verified: YES ✓

  The 18 free parameters (all ∈ (1/2)Z):
    # Name   Description                           Physical
  ─── ────── ─────────────────────────────────── ──────────
    1 p_μ    lepton gen-2 p                            -1/2
    2 q_μ    lepton gen-2 q                              -4
    3 r_μ    lepton gen-2 r (= r_τ)                       3
    4 p_τ    lepton gen-3 p                               1
    5 q_τ    lepton gen-3 q (= q_c = q_t)                -7
    6 p_u    up gen-1 p                                -3/2
    7 r_u    up gen-1 r                                  -1
    8 p_c    up gen-2 p                                 1/2
    9 r_c    up gen-2 r                                   3
   10 p_t    up gen-3 p                                   6
   11 p_d    down gen-1 p                              -1/2
   12 q_d    down gen-1 q                                -8
   13 r_d    down gen-1 r                                -2
   14 p_s    down gen-2 p (= p_b)                       3/2
   15 q_s    down gen-2 q                                -7
   16 r_s    down gen-2 r                                 0
   17 q_b    down gen-3 q                                -6
   18 r_b    down gen-3 r                                 4

  Parity constraints on derived quantities:
    q_u = 2(q_τ − q_μ) → always integer (since 2 × half-integer-diff)
    Physical: q_u = -6 (integer ✓)
    r_t = 2r_c − r_u − r_μ → ∈ (1/2)Z always, but:
    If r_u + r_μ ∈ Z: r_t ∈ Z
    If r_u + r_μ ∈ Z + 1/2: r_t ∈ Z + 1/2
    Physical: r_u + r_μ = -1 + 3 = 2 (integer → r_t ∈ Z)
    Physical r_t = 4 (integer ✓)

========================================================================
  §3. DECOMPOSITION INTO 7 INDEPENDENT GROUPS
========================================================================

  F_total = Σ_{s,k} [a_k(s)² + b_k(s)² + c_k(s)²]
  
  Each (s, k) pair depends on at most 3 free parameters.
  The dependencies cluster into 7 INDEPENDENT groups:

  G1: Lepton p-coordinates
       Free params (2): p_μ, p_τ
       Affects: lepton/p
       F_G1(physical) = 39/2 = 19.5000

  G2: Lepton+Up q-coordinates (coupled via universality)
       Free params (2): q_μ, q_τ
       Affects: lepton/q, up/q
       F_G2(physical) = 78 = 78.0000

  G3: Lepton+Up r-coordinates (coupled via universality)
       Free params (3): r_μ, r_u, r_c
       Affects: lepton/r, up/r
       F_G3(physical) = 233 = 233.0000

  G4: Up p-coordinates
       Free params (3): p_u, p_c, p_t
       Affects: up/p
       F_G4(physical) = 109/8 = 13.6250

  G5: Down p-coordinates (p is flat)
       Free params (2): p_d, p_s
       Affects: down/p
       F_G5(physical) = 185/4 = 46.2500

  G6: Down q-coordinates
       Free params (3): q_d, q_s, q_b
       Affects: down/q
       F_G6(physical) = 82 = 82.0000

  G7: Down r-coordinates
       Free params (3): r_d, r_s, r_b
       Affects: down/r
       F_G7(physical) = 6 = 6.0000

  F_total = F_G1 + ... + F_G7 = 3827/8 = 478.3750
  Verification: 478.375 = 478.3750 ✓

========================================================================
  §4. QUADRATIC FORMS AND ELLIPSOIDAL BOUNDS
========================================================================

  For each group, express F_i as a quadratic form in the free params:
    F_i(x) = x^T · A_i · x
  
  The eigenvalues of A_i determine the bounding ellipsoid.
  Within each ellipsoid, count the lattice points (step = 1/2).


  G1 (Lepton p-coordinates):
    Quadratic form matrix A_G1:
    [  +26.0000     -9.5000]
    [   -9.5000     +3.5000]
    Eigenvalues: ['0.0254', '29.4746']
    F_G1(physical) = 19.5000
    Ellipsoid semi-axes (√(F/λ)):
      λ_1 = 0.0254 → semi-axis = 27.68 → ~111 half-integer steps
        direction: -0.3435·p_μ, -0.9392·p_τ
      λ_2 = 29.4746 → semi-axis = 0.81 → ~4 half-integer steps
        direction: -0.9392·p_μ, +0.3435·p_τ
    Verification: x_phys^T A x_phys = 19.5000 (should be 19.5000)

  G2 (Lepton+Up q-coordinates (coupled via universality)):
    Quadratic form matrix A_G2:
    [  +88.0000    -46.5000]
    [  -46.5000    +26.0000]
    Eigenvalues: ['1.1140', '112.8860']
    F_G2(physical) = 78.0000
    Ellipsoid semi-axes (√(F/λ)):
      λ_1 = 1.1140 → semi-axis = 8.37 → ~34 half-integer steps
        direction: -0.4719·q_μ, -0.8817·q_τ
      λ_2 = 112.8860 → semi-axis = 0.83 → ~4 half-integer steps
        direction: -0.8817·q_μ, +0.4719·q_τ
    Verification: x_phys^T A x_phys = 78.0000 (should be 78.0000)

  G3 (Lepton+Up r-coordinates (coupled via universality)):
    Quadratic form matrix A_G3:
    [  +14.0000     -3.5000     +2.5000]
    [   -3.5000     +5.0000     -3.0000]
    [   +2.5000     -3.0000     +2.0000]
    Eigenvalues: ['0.1385', '4.8459', '16.0155']
    F_G3(physical) = 233.0000
    Ellipsoid semi-axes (√(F/λ)):
      λ_1 = 0.1385 → semi-axis = 41.01 → ~165 half-integer steps
        direction: -0.0258·r_μ, +0.5115·r_u, +0.8589·r_c
      λ_2 = 4.8459 → semi-axis = 6.93 → ~28 half-integer steps
        direction: +0.4237·r_μ, +0.7838·r_u, -0.4540·r_c
      λ_3 = 16.0155 → semi-axis = 3.81 → ~16 half-integer steps
        direction: -0.9054·r_μ, +0.3522·r_u, -0.2369·r_c
    Verification: x_phys^T A x_phys = 233.0000 (should be 233.0000)

  G4 (Up p-coordinates):
    Quadratic form matrix A_G4:
    [  +15.5000    -19.5000     +7.0000]
    [  -19.5000    +26.0000     -9.5000]
    [   +7.0000     -9.5000     +3.5000]
    Eigenvalues: ['0.0088', '0.6393', '44.3519']
    F_G4(physical) = 13.6250
    Ellipsoid semi-axes (√(F/λ)):
      λ_1 = 0.0088 → semi-axis = 39.31 → ~158 half-integer steps
        direction: +0.1324·p_u, +0.4264·p_c, +0.8948·p_t
      λ_2 = 0.6393 → semi-axis = 4.62 → ~19 half-integer steps
        direction: +0.8014·p_u, +0.4852·p_c, -0.3498·p_t
      λ_3 = 44.3519 → semi-axis = 0.55 → ~3 half-integer steps
        direction: -0.5833·p_u, +0.7634·p_c, -0.2775·p_t
    Verification: x_phys^T A x_phys = 13.6250 (should be 13.6250)

  G5 (Down p-coordinates (p is flat)):
    Quadratic form matrix A_G5:
    [  +15.5000    -12.5000]
    [  -12.5000    +10.5000]
    Eigenvalues: ['0.2525', '25.7475']
    F_G5(physical) = 46.2500
    Ellipsoid semi-axes (√(F/λ)):
      λ_1 = 0.2525 → semi-axis = 13.54 → ~55 half-integer steps
        direction: -0.6340·p_d, -0.7733·p_s
      λ_2 = 25.7475 → semi-axis = 1.34 → ~6 half-integer steps
        direction: -0.7733·p_d, +0.6340·p_s
    Verification: x_phys^T A x_phys = 46.2500 (should be 46.2500)

  G6 (Down q-coordinates):
    Quadratic form matrix A_G6:
    [  +15.5000    -19.5000     +7.0000]
    [  -19.5000    +26.0000     -9.5000]
    [   +7.0000     -9.5000     +3.5000]
    Eigenvalues: ['0.0088', '0.6393', '44.3519']
    F_G6(physical) = 82.0000
    Ellipsoid semi-axes (√(F/λ)):
      λ_1 = 0.0088 → semi-axis = 96.44 → ~386 half-integer steps
        direction: +0.1324·q_d, +0.4264·q_s, +0.8948·q_b
      λ_2 = 0.6393 → semi-axis = 11.33 → ~46 half-integer steps
        direction: +0.8014·q_d, +0.4852·q_s, -0.3498·q_b
      λ_3 = 44.3519 → semi-axis = 1.36 → ~6 half-integer steps
        direction: -0.5833·q_d, +0.7634·q_s, -0.2775·q_b
    Verification: x_phys^T A x_phys = 82.0000 (should be 82.0000)

  G7 (Down r-coordinates):
    Quadratic form matrix A_G7:
    [  +15.5000    -19.5000     +7.0000]
    [  -19.5000    +26.0000     -9.5000]
    [   +7.0000     -9.5000     +3.5000]
    Eigenvalues: ['0.0088', '0.6393', '44.3519']
    F_G7(physical) = 6.0000
    Ellipsoid semi-axes (√(F/λ)):
      λ_1 = 0.0088 → semi-axis = 26.09 → ~105 half-integer steps
        direction: +0.1324·r_d, +0.4264·r_s, +0.8948·r_b
      λ_2 = 0.6393 → semi-axis = 3.06 → ~13 half-integer steps
        direction: +0.8014·r_d, +0.4852·r_s, -0.3498·r_b
      λ_3 = 44.3519 → semi-axis = 0.37 → ~2 half-integer steps
        direction: -0.5833·r_d, +0.7634·r_s, -0.2775·r_b
    Verification: x_phys^T A x_phys = 6.0000 (should be 6.0000)

========================================================================
  §5. LATTICE POINT ENUMERATION PER GROUP
========================================================================

  For each group, enumerate ALL half-integer lattice points satisfying
  F_i ≤ F_i(physical). This gives the complete set of "candidates"
  that could potentially be part of a solution with F_total ≤ F_phys.


  G1 (F ≤ 19.5000, 2D):
    Total lattice points: 291
    Top 15 (lowest F):
       #                          Coords           F   Phys?
    ────  ──────────────────────────────  ──────────  ──────
       1  (                        0, 0)      0.0000  
       2  (                  -1/2, -3/2)      0.1250  
       3  (                    1/2, 3/2)      0.1250  
       4  (                    -1, -5/2)      0.3750  
       5  (                      1, 5/2)      0.3750  
       6  (                    -3/2, -4)      0.5000  
       7  (                      -1, -3)      0.5000  
       8  (                    -1/2, -1)      0.5000  
       9  (                      1/2, 1)      0.5000  
      10  (                        1, 3)      0.5000  
      11  (                      3/2, 4)      0.5000  
      12  (                   -2, -11/2)      0.8750  
      13  (                     0, -1/2)      0.8750  
      14  (                      0, 1/2)      0.8750  
      15  (                     2, 11/2)      0.8750  

    Physical solution rank: 285 / 291

  G2 (F ≤ 78.0000, 2D):
    Total lattice points: 87
    Top 15 (lowest F):
       #                          Coords           F   Phys?
    ────  ──────────────────────────────  ──────────  ──────
       1  (                        0, 0)      0.0000  
       2  (                    -1/2, -1)      1.5000  
       3  (                      1/2, 1)      1.5000  
       4  (                  -1/2, -1/2)      5.2500  
       5  (                    1/2, 1/2)      5.2500  
       6  (                      -1, -2)      6.0000  
       7  (                        1, 2)      6.0000  
       8  (                     0, -1/2)      6.5000  
       9  (                      0, 1/2)      6.5000  
      10  (                    -1, -3/2)      7.0000  
      11  (                      1, 3/2)      7.0000  
      12  (                  -1/2, -3/2)     10.7500  
      13  (                    1/2, 3/2)     10.7500  
      14  (                  -3/2, -5/2)     11.7500  
      15  (                    3/2, 5/2)     11.7500  

    Physical solution rank: 86 / 87

  G3 (F ≤ 233.0000, 3D):
    Total lattice points: 36401
    Top 15 (lowest F):
       #                          Coords           F   Phys?
    ────  ──────────────────────────────  ──────────  ──────
       1  (                     0, 0, 0)      0.0000  
       2  (                 0, -1/2, -1)      0.2500  
       3  (               0, -1/2, -1/2)      0.2500  
       4  (                 0, 1/2, 1/2)      0.2500  
       5  (                   0, 1/2, 1)      0.2500  
       6  (                 0, -1, -3/2)      0.5000  
       7  (                  0, 0, -1/2)      0.5000  
       8  (                   0, 0, 1/2)      0.5000  
       9  (                   0, 1, 3/2)      0.5000  
      10  (                   0, -1, -2)      1.0000  
      11  (                   0, -1, -1)      1.0000  
      12  (                     0, 1, 1)      1.0000  
      13  (                     0, 1, 2)      1.0000  
      14  (               0, -3/2, -5/2)      1.2500  
      15  (                 0, -3/2, -2)      1.2500  

    Physical solution rank: 36379 / 36401

  G4 (F ≤ 13.6250, 3D):
    Total lattice points: 3435
    Top 15 (lowest F):
       #                          Coords           F   Phys?
    ────  ──────────────────────────────  ──────────  ──────
       1  (                     0, 0, 0)      0.0000  
       2  (              -1/2, -3/2, -3)      0.1250  
       3  (               0, -1/2, -3/2)      0.1250  
       4  (                 0, 1/2, 3/2)      0.1250  
       5  (                 1/2, 3/2, 3)      0.1250  
       6  (              -1/2, -2, -9/2)      0.2500  
       7  (              -1/2, -1, -3/2)      0.2500  
       8  (            -1/2, -1/2, -1/2)      0.2500  
       9  (               1/2, 1/2, 1/2)      0.2500  
      10  (                 1/2, 1, 3/2)      0.2500  
      11  (                 1/2, 2, 9/2)      0.2500  
      12  (                -1, -2, -7/2)      0.3750  
      13  (                -1/2, -1, -2)      0.3750  
      14  (                  -1/2, 0, 1)      0.3750  
      15  (                 0, -1, -5/2)      0.3750  

    Physical solution rank: 3410 / 3435

  G5 (F ≤ 46.2500, 2D):
    Total lattice points: 233
    Top 15 (lowest F):
       #                          Coords           F   Phys?
    ────  ──────────────────────────────  ──────────  ──────
       1  (                        0, 0)      0.0000  
       2  (                  -1/2, -1/2)      0.2500  
       3  (                    1/2, 1/2)      0.2500  
       4  (                      -1, -1)      1.0000  
       5  (                        1, 1)      1.0000  
       6  (                    -1, -3/2)      1.6250  
       7  (                      1, 3/2)      1.6250  
       8  (                    -3/2, -2)      1.8750  
       9  (                    -1/2, -1)      1.8750  
      10  (                      1/2, 1)      1.8750  
      11  (                      3/2, 2)      1.8750  
      12  (                  -3/2, -3/2)      2.2500  
      13  (                    3/2, 3/2)      2.2500  
      14  (                    -2, -5/2)      2.6250  
      15  (                     0, -1/2)      2.6250  

    Physical solution rank: 231 / 233

  G6 (F ≤ 82.0000, 3D):
    Total lattice points: 49729
    Top 15 (lowest F):
       #                          Coords           F   Phys?
    ────  ──────────────────────────────  ──────────  ──────
       1  (                     0, 0, 0)      0.0000  
       2  (              -1/2, -3/2, -3)      0.1250  
       3  (               0, -1/2, -3/2)      0.1250  
       4  (                 0, 1/2, 3/2)      0.1250  
       5  (                 1/2, 3/2, 3)      0.1250  
       6  (              -1/2, -2, -9/2)      0.2500  
       7  (              -1/2, -1, -3/2)      0.2500  
       8  (            -1/2, -1/2, -1/2)      0.2500  
       9  (               1/2, 1/2, 1/2)      0.2500  
      10  (                 1/2, 1, 3/2)      0.2500  
      11  (                 1/2, 2, 9/2)      0.2500  
      12  (                -1, -2, -7/2)      0.3750  
      13  (                -1/2, -1, -2)      0.3750  
      14  (                  -1/2, 0, 1)      0.3750  
      15  (                 0, -1, -5/2)      0.3750  

    Physical solution rank: 49699 / 49729

  G7 (F ≤ 6.0000, 3D):
    Total lattice points: 985
    Top 15 (lowest F):
       #                          Coords           F   Phys?
    ────  ──────────────────────────────  ──────────  ──────
       1  (                     0, 0, 0)      0.0000  
       2  (              -1/2, -3/2, -3)      0.1250  
       3  (               0, -1/2, -3/2)      0.1250  
       4  (                 0, 1/2, 3/2)      0.1250  
       5  (                 1/2, 3/2, 3)      0.1250  
       6  (              -1/2, -2, -9/2)      0.2500  
       7  (              -1/2, -1, -3/2)      0.2500  
       8  (            -1/2, -1/2, -1/2)      0.2500  
       9  (               1/2, 1/2, 1/2)      0.2500  
      10  (                 1/2, 1, 3/2)      0.2500  
      11  (                 1/2, 2, 9/2)      0.2500  
      12  (                -1, -2, -7/2)      0.3750  
      13  (                -1/2, -1, -2)      0.3750  
      14  (                  -1/2, 0, 1)      0.3750  
      15  (                 0, -1, -5/2)      0.3750  

    Physical solution rank: 970 / 985

========================================================================
  §6. COMBINED SEARCH SPACE ANALYSIS
========================================================================

  Group sizes:
    G1: 291 lattice points
    G2: 87 lattice points
    G3: 36401 lattice points
    G4: 3435 lattice points
    G5: 233 lattice points
    G6: 49729 lattice points
    G7: 985 lattice points

  Naive product: 36,128,853,103,183,744,904,775
  (all combinations of group solutions with each F_i ≤ F_i(phys))

  Key insight: since we bound each F_i ≤ F_i^(phys),
  ALL 36,128,853,103,183,744,904,775 combinations satisfy F_total ≤ F_total^(phys).
  This is an UPPER BOUND on the search space.

========================================================================
  §7. NON-TRIVIALITY CONSTRAINTS
========================================================================

  The trivial solution (all coordinates = 0) gives F_total = 0.
  To select the physical solution, we need non-triviality constraints.
  
  Options:
    NT1: All 9 fermions have DISTINCT coordinate vectors
    NT2: Within each sector, masses are strictly ordered (gen hierarchy)
    NT3: Within each sector, at least 2 generations differ
    NT4: The curvature a ≠ 0 for at least one coordinate per sector
  
  Let's see which constraints are needed and their effect.

  Trivial solutions per group (all params = 0, F = 0):
    G1: zero solution exists = True
    G2: zero solution exists = True
    G3: zero solution exists = True
    G4: zero solution exists = True
    G5: zero solution exists = True
    G6: zero solution exists = True
    G7: zero solution exists = True

  Non-trivial solutions per group (at least one coord ≠ 0):
    G1: 290 non-trivial (of 291 total)
    G2: 86 non-trivial (of 87 total)
    G3: 36400 non-trivial (of 36401 total)
    G4: 3434 non-trivial (of 3435 total)
    G5: 232 non-trivial (of 233 total)
    G6: 49728 non-trivial (of 49729 total)
    G7: 984 non-trivial (of 985 total)

  'Generation-distinct' solutions per group:
  (all generations give different coordinate values)
    G1: 280 generation-distinct solutions
    G2: 64 generation-distinct solutions
    G3: 32398 generation-distinct solutions
    G4: 3302 generation-distinct solutions
    G5: 206 generation-distinct solutions
    G6: 48854 generation-distinct solutions
    G7: 928 generation-distinct solutions

========================================================================
  §8. LATTICE-MINIMAL SOLUTIONS PER GROUP
========================================================================

  Since F_total = F_G1 + ... + F_G7 and groups are independent,
  the global F_total minimum is the SUM of individual group minima.
  
  For each group, find the minimum non-trivial F value and the
  corresponding coordinates. Then compare with physical values.


  G1 (Lepton p-coordinates):
    Minimum non-trivial F = 0.1250
    Physical F = 19.5000
    Number of solutions at minimum: 2
    Physical IS the minimum: NO
    Minimum solutions:
      (-1/2, -3/2)
      (1/2, 3/2)
    Physical: (-1/2, 1)

  G2 (Lepton+Up q-coordinates (coupled via universality)):
    Minimum non-trivial F = 1.5000
    Physical F = 78.0000
    Number of solutions at minimum: 2
    Physical IS the minimum: NO
    Minimum solutions:
      (-1/2, -1)
      (1/2, 1)
    Physical: (-4, -7)

  G3 (Lepton+Up r-coordinates (coupled via universality)):
    Minimum non-trivial F = 0.2500
    Physical F = 233.0000
    Number of solutions at minimum: 4
    Physical IS the minimum: NO
    Minimum solutions:
      (0, -1/2, -1)
      (0, -1/2, -1/2)
      (0, 1/2, 1/2)
      (0, 1/2, 1)
    Physical: (3, -1, 3)

  G4 (Up p-coordinates):
    Minimum non-trivial F = 0.1250
    Physical F = 13.6250
    Number of solutions at minimum: 4
    Physical IS the minimum: NO
    Minimum solutions:
      (-1/2, -3/2, -3)
      (0, -1/2, -3/2)
      (0, 1/2, 3/2)
      (1/2, 3/2, 3)
    Physical: (-3/2, 1/2, 6)

  G5 (Down p-coordinates (p is flat)):
    Minimum non-trivial F = 0.2500
    Physical F = 46.2500
    Number of solutions at minimum: 2
    Physical IS the minimum: NO
    Minimum solutions:
      (-1/2, -1/2)
      (1/2, 1/2)
    Physical: (-1/2, 3/2)

  G6 (Down q-coordinates):
    Minimum non-trivial F = 0.1250
    Physical F = 82.0000
    Number of solutions at minimum: 4
    Physical IS the minimum: NO
    Minimum solutions:
      (-1/2, -3/2, -3)
      (0, -1/2, -3/2)
      (0, 1/2, 3/2)
      (1/2, 3/2, 3)
    Physical: (-8, -7, -6)

  G7 (Down r-coordinates):
    Minimum non-trivial F = 0.1250
    Physical F = 6.0000
    Number of solutions at minimum: 4
    Physical IS the minimum: NO
    Minimum solutions:
      (-1/2, -3/2, -3)
      (0, -1/2, -3/2)
      (0, 1/2, 3/2)
      (1/2, 3/2, 3)
    Physical: (-2, 0, 4)

  =============================================================
  COMBINED RESULTS:
  ============================================================
  Sum of group non-trivial minima: 2.5000
  Physical F_total:                478.3750
  Ratio (phys / min):              191.3500
  Physical solution is group-wise minimal: NO

  Groups where physical ≠ minimum:
    G1: phys F = 19.5000, min F = 0.1250, gap = 19.3750
    G2: phys F = 78.0000, min F = 1.5000, gap = 76.5000
    G3: phys F = 233.0000, min F = 0.2500, gap = 232.7500
    G4: phys F = 13.6250, min F = 0.1250, gap = 13.5000
    G5: phys F = 46.2500, min F = 0.2500, gap = 46.0000
    G6: phys F = 82.0000, min F = 0.1250, gap = 81.8750
    G7: phys F = 6.0000, min F = 0.1250, gap = 5.8750

========================================================================
  SUMMARY
========================================================================

  §1: The W-H-Z triangle has an angle of 84.70° at H — 
      nearly a right angle (deviation 5.30°).
      
  §2: 18 free parameters identified in coordinate language,
      all constrained to (1/2)Z.
      
  §3: F_total decomposes into 7 INDEPENDENT groups.
      Each group can be bounded and enumerated separately.
      
  §4-5: Ellipsoidal bounds from F_i ≤ F_i^(phys) yield finite
         lattice point counts per group.

  Group lattice point counts:
    G1 (2D, F≤19.5): 291 points
    G2 (2D, F≤78.0): 87 points
    G3 (3D, F≤233.0): 36401 points
    G4 (3D, F≤13.6): 3435 points
    G5 (2D, F≤46.2): 233 points
    G6 (3D, F≤82.0): 49729 points
    G7 (3D, F≤6.0): 985 points

  Total naive search space: 36,128,853,103,183,744,904,775

  Physical F_total = 478.3750
  Sum of non-trivial group minima = 2.5000
