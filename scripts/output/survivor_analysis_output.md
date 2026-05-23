========================================================================
  SURVIVOR ANALYSIS — Anatomy of the 72 Configs at 2%
========================================================================

§1. REBUILDING JOINT CONFIGURATION POOL
------------------------------------------------------------------------
  Joint configs with gen-3 gap structure: 2500
  At 10%: 2500
  At  5%: 625
  At  3%: 224
  At  2%: 72


§2. WHERE IS THE DEGENERACY? — Fixed vs Variable Coordinates
------------------------------------------------------------------------

  Analyzing 72 configs at 2% tolerance.

  === LEPTON SECTOR ===
    gen-1: 1 unique values
      ( +0.0,  +0.0,  +0.0)  ×72 (100%) ★ ACTUAL
    gen-2: 1 unique values
      ( -0.5,  -4.0,  +3.0)  ×72 (100%) ★ ACTUAL
    gen-3: 1 unique values
      ( +1.0,  -7.0,  +3.0)  ×72 (100%) ★ ACTUAL

  === UP SECTOR ===
    gen-1: 8 unique values
      ( -4.5,  -4.0,  +2.0)  ×9 (12%)
      ( -2.5,  -9.0,  -2.0)  ×9 (12%)
      ( -1.5,  -6.0,  -1.0)  ×9 (12%) ★ ACTUAL
      ( -0.5,  -3.0,  +0.0)  ×9 (12%)
      ( +0.5,  +0.0,  +1.0)  ×9 (12%)
      ( +4.5,  +1.0,  -1.0)  ×9 (12%)
      ( +8.5, -11.0, -10.0)  ×9 (12%)
      ( +9.5,  -8.0,  -9.0)  ×9 (12%)
    gen-2: 8 unique values
      ( -1.0,  -6.0,  +4.5)  ×9 (12%)
      ( +0.0,  -8.5,  +2.5)  ×9 (12%)
      ( +0.5,  -7.0,  +3.0)  ×9 (12%) ★ ACTUAL
      ( +1.0,  -5.5,  +3.5)  ×9 (12%)
      ( +1.5,  -4.0,  +4.0)  ×9 (12%)
      ( +3.5,  -3.5,  +3.0)  ×9 (12%)
      ( +5.5,  -9.5,  -1.5)  ×9 (12%)
      ( +6.0,  -8.0,  -1.0)  ×9 (12%)
    gen-3: 1 unique values
      ( +6.0,  -7.0,  +4.0)  ×72 (100%) ★ ACTUAL

  === DOWN SECTOR ===
    gen-1: 9 unique values
      ( -2.5, -14.0,  -4.0)  ×8 (11%)
      ( -1.5, -11.0,  -3.0)  ×8 (11%)
      ( -0.5,  -8.0,  -2.0)  ×8 (11%) ★ ACTUAL
      ( +0.5,  -5.0,  -1.0)  ×8 (11%)
      ( +1.5,  -2.0,  +0.0)  ×8 (11%)
      ( +3.5,  -7.0,  -4.0)  ×8 (11%)
      ( +4.5,  -4.0,  -3.0)  ×8 (11%)
      ( +5.5,  -1.0,  -2.0)  ×8 (11%)
      ( +6.5,  +2.0,  -1.0)  ×8 (11%)
    gen-2: 9 unique values
      ( +0.5, -10.0,  -1.0)  ×8 (11%)
      ( +1.0,  -8.5,  -0.5)  ×8 (11%)
      ( +1.5,  -7.0,  +0.0)  ×8 (11%) ★ ACTUAL
      ( +2.0,  -5.5,  +0.5)  ×8 (11%)
      ( +2.5,  -4.0,  +1.0)  ×8 (11%)
      ( +3.5,  -6.5,  -1.0)  ×8 (11%)
      ( +4.0,  -5.0,  -0.5)  ×8 (11%)
      ( +4.5,  -3.5,  +0.0)  ×8 (11%)
      ( +5.0,  -2.0,  +0.5)  ×8 (11%)
    gen-3: 1 unique values
      ( +1.5,  -6.0,  +4.0)  ×72 (100%) ★ ACTUAL


§3. COORDINATE-BY-COORDINATE — Which components vary?
------------------------------------------------------------------------

  === LEPTON ===
    gen-1:
      p: FIXED at +0.0
      q: FIXED at +0.0
      r: FIXED at +0.0
    gen-2:
      p: FIXED at -0.5
      q: FIXED at -4.0
      r: FIXED at +3.0
    gen-3:
      p: FIXED at +1.0
      q: FIXED at -7.0
      r: FIXED at +3.0

  === UP ===
    gen-1:
      p: range [-4.5, +9.5], span=14.0, 8 values
      q: range [-11.0, +1.0], span=12.0, 8 values
      r: range [-10.0, +2.0], span=12.0, 7 values
    gen-2:
      p: range [-1.0, +6.0], span=7.0, 8 values
      q: range [-9.5, -3.5], span=6.0, 8 values
      r: range [-1.5, +4.5], span=6.0, 7 values
    gen-3:
      p: FIXED at +6.0
      q: FIXED at -7.0
      r: FIXED at +4.0

  === DOWN ===
    gen-1:
      p: range [-2.5, +6.5], span=9.0, 9 values
      q: range [-14.0, +2.0], span=16.0, 9 values
      r: range [-4.0, +0.0], span=4.0, 5 values
    gen-2:
      p: range [+0.5, +5.0], span=4.5, 9 values
      q: range [-10.0, -2.0], span=8.0, 9 values
      r: range [-1.0, +1.0], span=2.0, 5 values
    gen-3:
      p: FIXED at +1.5
      q: FIXED at -6.0
      r: FIXED at +4.0


§4. GEN-1 ANATOMY — Decomposing the Degeneracy
------------------------------------------------------------------------

  Up gen-1 p values:     [np.float64(-4.5), np.float64(-2.5), np.float64(-1.5), np.float64(-0.5), np.float64(0.5), np.float64(4.5), np.float64(8.5), np.float64(9.5)]
  Up gen-1 q values:     [np.float64(-11.0), np.float64(-9.0), np.float64(-8.0), np.float64(-6.0), np.float64(-4.0), np.float64(-3.0), np.float64(0.0), np.float64(1.0)]
  Up gen-1 r values:     [np.float64(-10.0), np.float64(-9.0), np.float64(-2.0), np.float64(-1.0), np.float64(0.0), np.float64(1.0), np.float64(2.0)]
  Down gen-1 p values:   [np.float64(-2.5), np.float64(-1.5), np.float64(-0.5), np.float64(0.5), np.float64(1.5), np.float64(3.5), np.float64(4.5), np.float64(5.5), np.float64(6.5)]
  Down gen-1 q values:   [np.float64(-14.0), np.float64(-11.0), np.float64(-8.0), np.float64(-7.0), np.float64(-5.0), np.float64(-4.0), np.float64(-2.0), np.float64(-1.0), np.float64(2.0)]
  Down gen-1 r values:   [np.float64(-4.0), np.float64(-3.0), np.float64(-2.0), np.float64(-1.0), np.float64(0.0)]

  Unique (up_gen1, down_gen1) pairs: 72
  Unique up gen-1: 8
  Unique down gen-1: 9
  Product (independent): 72
  Actual pairs: 72
  → Up and Down gen-1 vary INDEPENDENTLY (no cross-constraint)


§5. VELOCITY B DISTRIBUTION — Across Survivors
------------------------------------------------------------------------

  === LEPTON VELOCITY B ===
    B_p: FIXED at -1.50
    B_q: FIXED at -4.50
    B_r: FIXED at +4.50
    Unique B vectors: 1
    Actual B: (-1.50, -4.50, +4.50)

  === UP VELOCITY B ===
    B_p: 8 values, range [-5.25, +1.75]
        +1.75  ×9
        +0.75  ×9
        +0.25  ×9 ★
        -0.25  ×9
        -0.75  ×9
    B_q: 8 values, range [-5.00, +1.00]
        -2.50  ×9
        +0.00  ×9
        -1.50  ×9 ★
        -3.00  ×9
        -4.50  ×9
    B_r: 7 values, range [+4.00, +10.00]
        +5.50  ×18 ★
        +4.00  ×9
        +6.00  ×9
        +5.00  ×9
        +4.50  ×9
    Unique B vectors: 8
    Actual B: (+0.25, -1.50, +5.50)

  === DOWN VELOCITY B ===
    B_p: 9 values, range [-0.50, +4.00]
        +4.00  ×8
        +3.50  ×8
        +3.00  ×8 ★
        +2.50  ×8
        +2.00  ×8
    B_q: 9 values, range [-4.00, +4.00]
        +4.00  ×8
        +2.50  ×8
        +1.00  ×8 ★
        -0.50  ×8
        -2.00  ×8
    B_r: 5 values, range [+0.00, +2.00]
        +2.00  ×16
        +1.50  ×16
        +1.00  ×16 ★
        +0.50  ×16
        +0.00  ×8
    Unique B vectors: 9
    Actual B: (+3.00, +1.00, +1.00)


§6. CURVATURE A — Consistency Check
------------------------------------------------------------------------
   lepton: A = -1.252103  (fixed by construction)
       up: A = -0.732242  (fixed by construction)
     down: A = +0.405465  (fixed by construction)

  ✓ Curvature A is IDENTICAL for all survivors (it's an input, not a result).


§7. DOWN gen-1 PATTERNS — Hidden Arithmetic?
------------------------------------------------------------------------

  9 unique down gen-1 coordinates:

    ( -2.5, -14.0,  -4.0)  m_d=4.602 MeV  p+q+r=-20.5  p-r= +1.5  2p+r= -9.0  p·ln2+q·lnR3+r·ln3=+2.1980
    ( -1.5, -11.0,  -3.0)  m_d=4.638 MeV  p+q+r=-15.5  p-r= +1.5  2p+r= -6.0  p·ln2+q·lnR3+r·ln3=+2.2057
    ( -0.5,  -8.0,  -2.0)  m_d=4.674 MeV  p+q+r=-10.5  p-r= +1.5  2p+r= -3.0  p·ln2+q·lnR3+r·ln3=+2.2135 ★ ACTUAL
    ( +0.5,  -5.0,  -1.0)  m_d=4.711 MeV  p+q+r= -5.5  p-r= +1.5  2p+r= +0.0  p·ln2+q·lnR3+r·ln3=+2.2213
    ( +1.5,  -2.0,  +0.0)  m_d=4.748 MeV  p+q+r= -0.5  p-r= +1.5  2p+r= +3.0  p·ln2+q·lnR3+r·ln3=+2.2290
    ( +3.5,  -7.0,  -4.0)  m_d=4.585 MeV  p+q+r= -7.5  p-r= +7.5  2p+r= +3.0  p·ln2+q·lnR3+r·ln3=+2.1942
    ( +4.5,  -4.0,  -3.0)  m_d=4.621 MeV  p+q+r= -2.5  p-r= +7.5  2p+r= +6.0  p·ln2+q·lnR3+r·ln3=+2.2020
    ( +5.5,  -1.0,  -2.0)  m_d=4.657 MeV  p+q+r= +2.5  p-r= +7.5  2p+r= +9.0  p·ln2+q·lnR3+r·ln3=+2.2097
    ( +6.5,  +2.0,  -1.0)  m_d=4.693 MeV  p+q+r= +7.5  p-r= +7.5  2p+r=+12.0  p·ln2+q·lnR3+r·ln3=+2.2175

  Looking for linear constraints on down gen-1...

  8 unique up gen-1 coordinates:

    ( -4.5,  -4.0,  +2.0)  m_u=2.193 MeV  p+q+r= -6.5  p-r= -6.5
    ( -2.5,  -9.0,  -2.0)  m_u=2.118 MeV  p+q+r=-13.5  p-r= -0.5
    ( -1.5,  -6.0,  -1.0)  m_u=2.135 MeV  p+q+r= -8.5  p-r= -0.5 ★ ACTUAL
    ( -0.5,  -3.0,  +0.0)  m_u=2.151 MeV  p+q+r= -3.5  p-r= -0.5
    ( +0.5,  +0.0,  +1.0)  m_u=2.168 MeV  p+q+r= +1.5  p-r= -0.5
    ( +4.5,  +1.0,  -1.0)  m_u=2.127 MeV  p+q+r= +4.5  p-r= +5.5
    ( +8.5, -11.0, -10.0)  m_u=2.172 MeV  p+q+r=-12.5  p-r=+18.5
    ( +9.5,  -8.0,  -9.0)  m_u=2.189 MeV  p+q+r= -7.5  p-r=+18.5

  Looking for linear constraints on up gen-1...


§8. INTER-SECTOR RELATION — u↔d Gen-1 Distance
------------------------------------------------------------------------

  For each survivor, compute Δ₁ = d(gen-1) - u(gen-1):
  (This is the gen-1 displacement between up and down sectors)

  Unique Δ₁ = d₁ - u₁ vectors: 36

    Δ₁=(-12.0, -6.0, +5.0)  |Δ|_lattice=14.32  Δ(lnm)=+0.7433  ×1
    Δ₁=(-11.0, -3.0, +6.0)  |Δ|_lattice=12.88  Δ(lnm)=+0.7510  ×2
    Δ₁=(-10.0, +0.0, +7.0)  |Δ|_lattice=12.21  Δ(lnm)=+0.7588  ×2
    Δ₁=( -9.0, +3.0, +8.0)  |Δ|_lattice=12.41  Δ(lnm)=+0.7666  ×2
    Δ₁=( -8.0, +6.0, +9.0)  |Δ|_lattice=13.45  Δ(lnm)=+0.7744  ×2
    Δ₁=( -7.0,-15.0, -3.0)  |Δ|_lattice=16.82  Δ(lnm)=+0.7721  ×1
    Δ₁=( -7.0, +9.0,+10.0)  |Δ|_lattice=15.17  Δ(lnm)=+0.7821  ×1
    Δ₁=( -6.0,-12.0, -2.0)  |Δ|_lattice=13.56  Δ(lnm)=+0.7798  ×1
    Δ₁=( -6.0, +1.0, +5.0)  |Δ|_lattice=7.87  Δ(lnm)=+0.7395  ×1
    Δ₁=( -5.0, -9.0, -1.0)  |Δ|_lattice=10.34  Δ(lnm)=+0.7876  ×1
    Δ₁=( -5.0, +4.0, +6.0)  |Δ|_lattice=8.77  Δ(lnm)=+0.7473  ×2
    Δ₁=( -4.0, -6.0, +0.0)  |Δ|_lattice=7.21  Δ(lnm)=+0.7954  ×1
    Δ₁=( -4.0, +7.0, +7.0)  |Δ|_lattice=10.68  Δ(lnm)=+0.7551  ×2
    Δ₁=( -3.0,-14.0, -5.0)  |Δ|_lattice=15.17  Δ(lnm)=+0.7528  ×1
    Δ₁=( -3.0, -3.0, +1.0)  |Δ|_lattice=4.36  Δ(lnm)=+0.8032  ×1

  ln(m_d/m_u) at gen-1:
    from coords: [np.float64(0.7375), np.float64(0.7395), np.float64(0.7413), np.float64(0.7433), np.float64(0.7453), np.float64(0.7473), np.float64(0.749), np.float64(0.751), np.float64(0.7528), np.float64(0.753), np.float64(0.7551), np.float64(0.7568), np.float64(0.7588), np.float64(0.7605), np.float64(0.7608), np.float64(0.7628), np.float64(0.7646), np.float64(0.7666), np.float64(0.7683), np.float64(0.7706), np.float64(0.7721), np.float64(0.7723), np.float64(0.7744), np.float64(0.7761), np.float64(0.7798), np.float64(0.7801), np.float64(0.7821), np.float64(0.7839), np.float64(0.7876), np.float64(0.7879), np.float64(0.7916), np.float64(0.7954), np.float64(0.7956), np.float64(0.7994), np.float64(0.8032), np.float64(0.8072)]
    experiment:   +0.7711


§9. RESIDUAL PATTERN — Error Sign Structure (+ → −)
------------------------------------------------------------------------

  For each survivor, compute the fractional error per generation.

   Config       e       μ       τ       u       c       t       d       s       b
  ───────────────────────────────────────────────────────────────────────────
      1   +0.0%  -0.4%  -0.2%  +1.5%  +0.0%  -1.5%  -1.4%  -1.4%  -0.7%
      2   +0.0%  -0.4%  -0.2%  +1.5%  +0.0%  -1.5%  -0.7%  -1.0%  -0.7%
      3   +0.0%  -0.4%  -0.2%  +1.5%  +0.0%  -1.5%  +0.1%  -0.6%  -0.7%
      4   +0.0%  -0.4%  -0.2%  +1.5%  +0.0%  -1.5%  +0.9%  -0.2%  -0.7%
      5   +0.0%  -0.4%  -0.2%  +1.5%  +0.0%  -1.5%  +1.7%  +0.2%  -0.7%
      6   +0.0%  -0.4%  -0.2%  +1.5%  +0.0%  -1.5%  -1.8%  -1.5%  -0.7%
      7   +0.0%  -0.4%  -0.2%  +1.5%  +0.0%  -1.5%  -1.1%  -1.2%  -0.7%
      8   +0.0%  -0.4%  -0.2%  +1.5%  +0.0%  -1.5%  -0.3%  -0.8%  -0.7%
      9   +0.0%  -0.4%  -0.2%  +1.5%  +0.0%  -1.5%  +0.5%  -0.4%  -0.7%
     10   +0.0%  -0.4%  -0.2%  -1.9%  -1.7%  -1.5%  -1.4%  -1.4%  -0.7%
     11   +0.0%  -0.4%  -0.2%  -1.9%  -1.7%  -1.5%  -0.7%  -1.0%  -0.7%
     12   +0.0%  -0.4%  -0.2%  -1.9%  -1.7%  -1.5%  +0.1%  -0.6%  -0.7%
     13   +0.0%  -0.4%  -0.2%  -1.9%  -1.7%  -1.5%  +0.9%  -0.2%  -0.7%
     14   +0.0%  -0.4%  -0.2%  -1.9%  -1.7%  -1.5%  +1.7%  +0.2%  -0.7%
     15   +0.0%  -0.4%  -0.2%  -1.9%  -1.7%  -1.5%  -1.8%  -1.5%  -0.7%
     21 ★   +0.0%  -0.4%  -0.2%  -1.2%  -1.3%  -1.5%  +0.1%  -0.6%  -0.7%
  (actual config is #21)


  SIGN PATTERN across generations (within each sector):
  ───────────────────────────────────────────────────

       Pattern (ℓ|u|d)  Count  %
          +--|---|---|     20  27.8%
          +--|+--|---|     15  20.8%
          +--|---|+--|     12  16.7%
          +--|+--|+--|      9  12.5%
          +--|++-|---|      5  6.9%
          +--|---|++-|      4  5.6%
          +--|++-|+--|      3  4.2%
          +--|+--|++-|      3  4.2%
          +--|++-|++-|      1  1.4%


§10. RESIDUAL DECOMPOSITION — Linear (β) vs Quadratic (γ)
------------------------------------------------------------------------

  For each sector, decompose the residual δ(g) into:
    δ(g) = α + β·(g-2) + γ·(g-2)²
  where:
    γ = (δ₁ - 2δ₂ + δ₃)/2  → curvature error
    β = (δ₃ - δ₁)/2         → velocity error
    α = δ₂                   → offset error

  === LEPTON ===
    γ (curvature err):  mean=+0.249%  std=0.000%  range=[+0.249%, +0.249%]
    β (velocity err):   mean=-0.118%  std=0.000%  range=[-0.118%, -0.118%]
    α (offset err):     mean=-0.367%  std=0.000%  range=[-0.367%, -0.367%]

  === UP ===
    γ (curvature err):  mean=-0.031%  std=0.004%  range=[-0.035%, -0.024%]
    β (velocity err):   mean=-0.667%  std=0.614%  range=[-1.515%, +0.224%]
    α (offset err):     mean=-0.798%  std=0.610%  range=[-1.685%, +0.043%]

  === DOWN ===
    γ (curvature err):  mean=+0.274%  std=0.003%  range=[+0.272%, +0.282%]
    β (velocity err):   mean=-0.245%  std=0.536%  range=[-1.196%, +0.545%]
    α (offset err):     mean=-0.757%  std=0.533%  range=[-1.544%, +0.186%]


§11. MINIMAL VELOCITY ERROR — Ranking by |β|
------------------------------------------------------------------------

  If the sign-flip pattern (+ → −) is a velocity error β,
  the "best" config should minimize |β| across all 3 sectors.

  Rank  Config    |β_ℓ|    |β_u|    |β_d|      Σ|β|  Actual?
  ─────────────────────────────────────────────────────────────────
     1      47   0.118%   0.026%   0.025%    0.169%  
     2      20   0.118%   0.159%   0.025%    0.302%  
     3      52   0.118%   0.026%   0.162%    0.306%  
     4      11   0.118%   0.224%   0.025%    0.366%  
     5      53   0.118%   0.026%   0.224%    0.369%  
     6      25   0.118%   0.159%   0.162%    0.439%  
     7      26   0.118%   0.159%   0.224%    0.501%  
     8      16   0.118%   0.224%   0.162%    0.503%  
     9      46   0.118%   0.026%   0.360%    0.504%  
    10      48   0.118%   0.026%   0.412%    0.556%  
    11      17   0.118%   0.224%   0.224%    0.566%  
    12      19   0.118%   0.159%   0.360%    0.637%  
    13      29   0.118%   0.544%   0.025%    0.687%  
    14      21   0.118%   0.159%   0.412%    0.689%  ★
    15      51   0.118%   0.026%   0.545%    0.689%  
    16      10   0.118%   0.224%   0.360%    0.701%  
    17      12   0.118%   0.224%   0.412%    0.754%  
    18      54   0.118%   0.026%   0.613%    0.758%  
    19      24   0.118%   0.159%   0.545%    0.822%  
    20      34   0.118%   0.544%   0.162%    0.824%  

  → ACTUAL config rank by |β|: #14 of 72


§12. THE CONSTRAINT THAT KILLS — Actual vs Others
------------------------------------------------------------------------

  Actual gen-1:  u=(-1.5,-6.0,-1.0)  d=(-0.5,-8.0,-2.0)
  Actual Δ₁:    (+1.0,-2.0,-1.0)

  Testing distinguishing properties:
    q₂ = q₃ in up: 9/72 configs  (actual: YES)
    q₂ = q₃ in down: 0/72 configs  (actual: NO)

  Coordinate type at gen-1 (integer vs half-integer):
    up gen-1 type pattern:
      ('H', 'I', 'I')  ×72 (100%) ★ ACTUAL
    down gen-1 type pattern:
      ('H', 'I', 'I')  ×72 (100%) ★ ACTUAL

  Sum of gen-1 coordinates p+q+r:
    up: values = [np.float64(-13.5), np.float64(-12.5), np.float64(-8.5), np.float64(-7.5), np.float64(-6.5), np.float64(-3.5), np.float64(1.5), np.float64(4.5)]
    actual sum = -8.5
    down: values = [np.float64(-20.5), np.float64(-15.5), np.float64(-10.5), np.float64(-7.5), np.float64(-5.5), np.float64(-2.5), np.float64(-0.5), np.float64(2.5), np.float64(7.5)]
    actual sum = -10.5

  Is velocity B a half-integer vector?
    up: B=(+0.25,-1.50,+5.50)  half-int: False
    down: B=(+3.00,+1.00,+1.00)  half-int: True

  How many survivors have half-integer B?
    up: 0/72 have half-int B
    down: 72/72 have half-int B


========================================================================
  §13. SUMMARY — Anatomy of Degeneracy
========================================================================

  FROM 2500 joint configs (gen-3 gap structure, 10% tolerance),
  72 survive at 2%.

  FIXED (identical across all survivors):
    - Lepton sector: fully fixed (gen-1 = origin, unique gen-3)
    - Gen-3 for all sectors: fully fixed by gap structure
    - Curvature A: fixed by construction
    - Gen-2: determined by midpoint formula

  VARIABLE:
    - Up gen-1 coordinates
    - Down gen-1 coordinates
    → The entire degeneracy is in WHERE gen-1 quarks sit on the lattice.

  The velocity B inherits this degeneracy (B depends on gen-1).
  The mass errors show a systematic + → − sign flip with generation,
  confirming the residual is dominated by a VELOCITY correction (β),
  not a curvature correction (γ ≈ 0).

  To break the degeneracy, we need a constraint on gen-1 positions.
  Candidates: CKM angles, velocity regularity, integer-type restrictions.

