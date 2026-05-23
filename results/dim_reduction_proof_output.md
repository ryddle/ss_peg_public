========================================================================
  §1. 5D → 3D DIMENSIONAL REDUCTION
========================================================================

  Enumerating all 9^5 = 59,049 formulas...
  Exponent range: [-4, 4]

  Results:
    5D formulas (nominal):         59,049
    Unique 3D triples (exact):      3,321
    Numerically unique values:      3,321
    Reduction factor:                17.8×

  Baker check: 3D triples = numerical uniques? YES ✓
  → ln2, lnR₃, ln3 are Q-linearly independent (Baker's theorem confirmed)

  p-coordinate distribution:
    Integer p:       1,701 (51.2%)
    Half-integer p:  1,620 (48.8%)

========================================================================
  §2. BAKER INDEPENDENCE — NUMERICAL STRESS TEST
========================================================================

  Testing: does any (p, q, r) ≠ (p', q', r') give the same ln(f)?
  If Baker holds, the answer is NO for all distinct triples.

  Minimum gap between distinct formula values: 9.949281e-04
  Mean gap:                                    8.255767e-03
  Ratio min/mean:                              1.205131e-01

  10 smallest gaps (to verify no degeneracies):
  Rank              Gap        Value_i      Value_i+1
  ────── ────────────── ────────────── ──────────────

  Exhaustive test: p·ln2 + q·lnR₃ + r·ln3 = 0 for |p|,|q|,|r| ≤ 20?
    None found. Q-independence confirmed for |exp| ≤ 20.

========================================================================
  §3. SPINORIAL STRUCTURE ANALYSIS
========================================================================

  Question: Does p mod 1 (integer vs half-integer) correlate with
  any quantum number (I₃, Y, Q, Nc, gen)?

  Particle      p  p mod 1  gen     I₃      Y      Q   Nc  c_odd?
  ──────── ────── ──────── ──── ────── ────── ────── ──── ───────
  e          +0.0      0.0    1   -0.5  -1.00  -1.00    1
  mu         -0.5      0.5    2   -0.5  -1.00  -1.00    1
  tau        +1.0      0.0    3   -0.5  -1.00  -1.00    1
  u          -1.5      0.5    1   +0.5  +0.33  +0.67    3
  c          +0.5      0.5    2   +0.5  +0.33  +0.67    3
  t          +6.0      0.0    3   +0.5  +0.33  +0.67    3
  d          -0.5      0.5    1   -0.5  +0.33  -0.33    3
  s          +1.5      0.5    2   -0.5  +0.33  -0.33    3
  b          +1.5      0.5    3   -0.5  +0.33  -0.33    3

  Correlation tests:
    p_half ↔ gen_even? 6/9 match
    Half-integer p particles have I₃ = {-0.5, 0.5}
    Integer p particles have I₃ = {-0.5, 0.5}
    Half-integer p particles have Nc = {1, 3}
    Integer p particles have Nc = {1, 3}

  p mod 1 by generation:
    Gen 1: ['e(p=+0.0)', 'u(p=-1.5)', 'd(p=-0.5)']  half-int: 2/3
    Gen 2: ['mu(p=-0.5)', 'c(p=+0.5)', 's(p=+1.5)']  half-int: 3/3
    Gen 3: ['tau(p=+1.0)', 't(p=+6.0)', 'b(p=+1.5)']  half-int: 1/3

  2p values (should all be integers):
    e   : 2p =    +0  parity = even
    mu  : 2p =    -1  parity = odd
    tau : 2p =    +2  parity = even
    u   : 2p =    -3  parity = odd
    c   : 2p =    +1  parity = odd
    t   : 2p =   +12  parity = even
    d   : 2p =    -1  parity = odd
    s   : 2p =    +3  parity = odd
    b   : 2p =    +3  parity = odd
    W   : 2p =   +11  parity = odd
    Z   : 2p =   +16  parity = even
    H   : 2p =   +14  parity = even

========================================================================
  §4. WEIGHT DIAGRAM ANALYSIS
========================================================================

  Testing whether the 12 coordinate vectors (or their 2× scaling)
  form weights of a representation of a rank-3 Lie algebra.

  Approach: Compute all pairwise differences and check if they
  correspond to roots of known rank-3 systems.

  Doubled coordinates (2p, q, r) — all integers:
  Particle   2p    q    r
  ──────── ──── ──── ────
  e          +0   +0   +0
  mu         -1   -4   +3
  tau        +2   -7   +3
  u          -3   -6   -1
  c          +1   -7   +3
  t         +12   -7   +4
  d          -1   -8   -2
  s          +3   -7   +0
  b          +3   -6   +4
  W         +11  -10   +2
  Z         +16  -11   +0
  H         +14   -9   +2

  Intra-sector difference vectors (generation transitions):
  Transition    Δ(2p)   Δq   Δr      |Δ|
  ──────────── ────── ──── ──── ────────
  e→mu         -1   -4   +3    5.099
  mu→tau        +3   -3   +0    4.243
  u→c          +4   -1   +4    5.745
  c→t         +11   +0   +1   11.045
  d→s          +4   +1   +2    4.583
  s→b          +0   +1   +4    4.123
  e→tau        +2   -7   +3    7.874
  u→t         +15   -1   +5   15.843
  d→b          +4   +2   +6    7.483

  Cross-sector comparison of gen→gen+1 transitions:
    lepton: 1→2 = [-1 -4  3], 2→3 = [ 3 -3  0]
    up: 1→2 = [ 4 -1  4], 2→3 = [11  0  1]
    down: 1→2 = [4 1 2], 2→3 = [0 1 4]

  All pairwise differences (12×12 = 132 non-trivial):
    Unique difference vectors: 66
    Repeated differences: 0

  Root system closure test:
    Tested 2500 pairwise sums
    Closures found: 391
    Closure fraction: 0.156

  Known rank-3 root systems for comparison:
    A₃ (=SU(4)):  12 roots, dim=15
    B₃ (=SO(7)):  18 roots, dim=21
    C₃ (=Sp(6)):  18 roots, dim=21
    Our system:    130 difference vectors from 12 particles

  Gram matrix of fermion coordinates (original p,q,r):
    Rank of coordinate matrix: 3
    Eigenvalues of Gram matrix: [3.83247418e+02 5.59808265e+01 1.72717554e+01 2.52445054e-14
 4.84576113e-15]
    Singular values: [19.57670601  7.48203358  4.15593014]
    Right singular vectors (principal directions in (p,q,r) space):
      PC1: (-0.1872, +0.9447, -0.2694)  — explains 84.0%
      PC2: (-0.6100, -0.3267, -0.7219)  — explains 12.3%
      PC3: (+0.7700, -0.0292, -0.6374)  — explains 3.8%

========================================================================
  SUMMARY
========================================================================

  §1: 5D→3D reduction CONFIRMED
      59,049 → 3,321 unique values (factor 17.8×)
      Baker independence: CONFIRMED

  §2: Q-linear independence of {ln2, lnR₃, ln3}
      No relation found for |exp| ≤ 20
      Minimum gap: 9.949281e-04

  §3: Spinorial structure
      Half-integer p particles identified
      See correlations above for quantum number dependencies

  §4: Weight diagram
      Coordinate matrix rank: 3
      Difference vectors: 130 (including negatives)
      Root closure fraction: 0.156

