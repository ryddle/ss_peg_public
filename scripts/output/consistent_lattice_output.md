========================================================================
  CONSISTENT MASS LATTICE FROM TRANSITIVITY CONSTRAINTS
========================================================================

§1. BASELINE — Individual search results (from mass_ratio_search.py)
------------------------------------------------------------------------

  Ratio           Exp.value    Formula   Error%                5D vector       3D (p,q,r)
  --------------------------------------------------------------------------------------
  m_mu/m_e         206.7683   206.0081   0.368%          (+0,-4,+1,+3,+0)      (-0.5,-4,+3)
  m_tau/m_mu        16.8170    16.8392   0.132%          (+0,-3,-1,+0,+1)      (+1.5,-3,+0)
  m_tau/m_e       3477.2283  3496.0728   0.542%          (-2,-4,+0,+4,+0)      (+2.0,-4,+4)
  m_c/m_u          587.9630   587.2241   0.126%          (-2,-1,+0,+4,+0)      (+2.0,-1,+4)
  m_t/m_c          135.8268   135.7645   0.046%          (-1,+0,-1,+1,+4)      (+5.5,+0,+1)
  m_s/m_d           19.8723    19.8629   0.047%          (-2,+1,+0,+2,+0)      (+2.0,+1,+2)
  m_b/m_s           44.7537    44.6916   0.139%          (+0,+1,+0,+4,+0)      (+0.0,+1,+4)
  m_t/m_b           41.2679    41.3303   0.151%          (-3,+2,-1,+1,+2)      (+5.5,+2,+1)
  m_d/m_u            2.1759     2.1899   0.643%          (-1,-2,+0,-1,+0)      (+1.0,-2,-1)
  m_b/m_tau          2.3525     2.3409   0.493%          (+0,+1,-1,+1,+0)      (+0.5,+1,+1)
  m_t/m_W            2.1461     2.1379   0.383%          (+0,+3,-1,+2,+0)      (+0.5,+3,+2)
  m_H/m_W            1.5564     1.5606   0.270%          (+0,+1,-1,+0,+1)      (+1.5,+1,+0)
  m_H/m_Z            1.3719     1.3699   0.144%          (+1,+2,+0,+2,+0)      (-1.0,+2,+2)
  m_H/m_t            0.7252     0.7300   0.655%          (-1,-2,+0,-2,+0)      (+1.0,-2,-2)

  Mean |error|:  0.2956%
  Max  |error|:  0.6554%
  RMS  error:    0.3620%


§2. TRANSITIVITY CHECK — Are searched vectors internally consistent?
------------------------------------------------------------------------

  Two disconnected trees: leptons (rooted at e) and quarks+bosons (rooted at u).
  Connected via m_b/m_τ:

  v(b) from lepton tree: v(τ) + v(m_b/m_τ) = [ 0. -6. -1.  4.  1.]
  v(b) from quark tree:  v(u) + v(d/u) + v(s/d) + v(b/s) = [-3.  0.  0.  5.  0.]

  Offset (difference) = [ 3. -6. -1. -1.  1.]
  This offset IS the "spectral address" of u relative to e:
    v(u) - v(e) = [ 3. -6. -1. -1.  1.]
    In 3D: [-1.5 -6.  -1. ]

  Particle coordinates (5D, from spanning tree + b/τ bridge):
  Particle                    5D vector         3D (p,q,r)   m_pred (MeV)  m_exp (MeV)   Error%
  --------------------------------------------------------------------------------------------
  e                    (+0,+0,+0,+0,+0)  ( +0.0,  +0,  +0)           0.51         0.51    0.00%
  mu                   (+0,-4,+1,+3,+0)  ( -0.5,  -4,  +3)         105.27       105.66    0.37%
  tau                  (+0,-7,+0,+3,+1)  ( +1.0,  -7,  +3)        1772.66      1776.86    0.24%
  u                    (+3,-6,-1,-1,+1)  ( -1.5,  -6,  -1)           2.13         2.16    1.18%
  d                    (+2,-8,-1,-2,+1)  ( -0.5,  -8,  -2)           4.67         4.70    0.54%
  c                    (+1,-7,-1,+3,+1)  ( +0.5,  -7,  +3)        1253.46      1270.00    1.30%
  s                    (+0,-7,-1,+0,+1)  ( +1.5,  -7,  +0)          92.85        93.40    0.59%
  t                    (+0,-7,-2,+4,+5)  ( +6.0,  -7,  +4)      170175.50    172500.00    1.35%
  b                    (+0,-6,-1,+4,+1)  ( +1.5,  -6,  +4)        4149.57      4180.00    0.73%
  W                   (+0,-10,-1,+2,+5)  ( +5.5, -10,  +2)       79600.56     80379.00    0.97%
  Z                   (-1,-11,-2,+0,+6)  ( +8.0, -11,  +0)       90679.16     91187.60    0.56%
  H                    (+0,-9,-2,+2,+6)  ( +7.0,  -9,  +2)      124223.07    125100.00    0.70%


§3. TRANSITIVITY ERRORS — Predicted vs experimental for ALL ratios
------------------------------------------------------------------------

  Now ALL ratios are DERIVED from particle coordinates (not searched individually).
  For redundant ratios (not in spanning tree), errors reveal transitivity violations.

  Ratio           Exp.value    Lattice   Error%   Searched  SrchErr%  In tree?
  --------------------------------------------------------------------------
  m_mu/m_e         206.7683   206.0081   0.368%    206.0081   0.368%       YES
  m_tau/m_mu        16.8170    16.8392   0.132%     16.8392   0.132%       YES
  m_tau/m_e       3477.2283  3469.0120   0.236%   3496.0728   0.542%        no ◄
  m_c/m_u          587.9630   587.2241   0.126%    587.2241   0.126%       YES
  m_t/m_c          135.8268   135.7645   0.046%    135.7645   0.046%       YES
  m_s/m_d           19.8723    19.8629   0.047%     19.8629   0.047%       YES
  m_b/m_s           44.7537    44.6916   0.139%     44.6916   0.139%       YES
  m_t/m_b           41.2679    41.0104   0.624%     41.3303   0.151%        no
  m_d/m_u            2.1759     2.1899   0.643%      2.1899   0.643%       YES
  m_b/m_tau          2.3525     2.3409   0.493%      2.3409   0.493%       YES
  m_t/m_W            2.1461     2.1379   0.383%      2.1379   0.383%       YES
  m_H/m_W            1.5564     1.5606   0.270%      1.5606   0.270%       YES
  m_H/m_Z            1.3719     1.3699   0.144%      1.3699   0.144%       YES
  m_H/m_t            0.7252     0.7300   0.655%      0.7300   0.655%        no

  LATTICE (transitivity-forced):
    Mean |error|:  0.3075%
    Max  |error|:  0.6554%
    RMS:           0.3745%

  INDIVIDUAL SEARCH (no transitivity):
    Mean |error|:  0.2956%
    Max  |error|:  0.6554%
    RMS:           0.3620%

  RATIO (lattice/individual): 1.04×
  → Transitivity-forced lattice is COMPARABLE to individual search.
    The lattice structure is real — consistency does NOT degrade quality.


§4. GLOBAL OPTIMIZATION — Best consistent lattice
------------------------------------------------------------------------

  The spanning tree approach assigns coordinates from one specific path.
  Different spanning trees (different bridges) give different results.
  A GLOBAL optimization finds the particle coordinates that minimize
  total squared error across ALL ratios simultaneously.

  Variables: 3D coordinates (p,q,r) for 11 particles (e fixed at origin)
  Constraints: q ∈ Z, r ∈ Z, p ∈ Z/2
  Objective: minimize Σ (ln(predicted_ratio) - ln(experimental_ratio))²

  System: 14 equations, 33 unknowns, rank = 11
  Relaxed solution residual: 0.000000

  Relaxed (continuous) particle coordinates:
  Particle        p        q        r
  ------------------------------------
  e           0.000    0.000    0.000
  mu          1.811   -1.553    2.870
  tau         2.769   -2.376    4.389
  u           0.490   -0.420    0.776
  c           2.655   -2.278    4.208
  t           4.323   -3.709    6.852
  d           0.754   -0.647    1.194
  s           1.769   -1.517    2.803
  b           3.060   -2.625    4.849
  W           4.064   -3.486    6.441
  Z           4.107   -3.523    6.509
  H           4.214   -3.615    6.679

  Starting optimization from SPANNING TREE coordinates (not from rounded relaxed).
  Rationale: rounding the relaxed solution lands far from the correct lattice point.
  The spanning tree already achieves 0.3075% mean error.

  Seed coordinates (from spanning tree):
  Particle        p        q        r
  --------------------------------
  e             0.0        0        0
  mu           -0.5       -4        3
  tau           1.0       -7        3
  u            -1.5       -6       -1
  c             0.5       -7        3
  t             6.0       -7        4
  d            -0.5       -8       -2
  s             1.5       -7        0
  b             1.5       -6        4
  W             5.5      -10        2
  Z             8.0      -11        0
  H             7.0       -9        2
  Converged after 1 iterations.

  Optimized lattice particle coordinates:
  Particle        p        q        r
  --------------------------------
  e             0.0        0        0
  mu           -0.5       -4        3
  tau           1.0       -7        3
  u            -1.5       -6       -1
  c             0.5       -7        3
  t             6.0       -7        4
  d            -0.5       -8       -2
  s             1.5       -7        0
  b             1.5       -6        4
  W             5.5      -10        2
  Z             8.0      -11        0
  H             7.0       -9        2

  Optimized lattice errors:
  Ratio           Exp.value    Lattice   Error%   Searched  SrchErr%  Better?
  --------------------------------------------------------------------------------
  m_mu/m_e         206.7683   206.0081   0.368%    206.0081   0.368%         
  m_tau/m_mu        16.8170    16.8392   0.132%     16.8392   0.132%         
  m_tau/m_e       3477.2283  3469.0120   0.236%   3496.0728   0.542%      YES
  m_c/m_u          587.9630   587.2241   0.126%    587.2241   0.126%         
  m_t/m_c          135.8268   135.7645   0.046%    135.7645   0.046%         
  m_s/m_d           19.8723    19.8629   0.047%     19.8629   0.047%         
  m_b/m_s           44.7537    44.6916   0.139%     44.6916   0.139%         
  m_t/m_b           41.2679    41.0104   0.624%     41.3303   0.151%         
  m_d/m_u            2.1759     2.1899   0.643%      2.1899   0.643%         
  m_b/m_tau          2.3525     2.3409   0.493%      2.3409   0.493%         
  m_t/m_W            2.1461     2.1379   0.383%      2.1379   0.383%         
  m_H/m_W            1.5564     1.5606   0.270%      1.5606   0.270%      YES
  m_H/m_Z            1.3719     1.3699   0.144%      1.3699   0.144%         
  m_H/m_t            0.7252     0.7300   0.655%      0.7300   0.655%         

  OPTIMIZED LATTICE (consistent):
    Mean |error|:  0.3075%
    Max  |error|:  0.6554%
    RMS:           0.3745%
    #2/14 ratios are BETTER than individual search

  INDIVIDUAL SEARCH:
    Mean |error|:  0.2956%

  RATIO: 1.04×


§5. CONSISTENT 5D EXPONENT VECTORS
------------------------------------------------------------------------

  The optimized 3D coordinates (p,q,r) define the LATTICE consistently.
  For each ratio, the 3D difference translates to a FAMILY of 5D vectors:
    p = e - a - c/2, q = b, r = d
  
  Multiple (a,c,e) triples give the same p, but with different physical
  meaning (number of bridge crossings c vs weak bounces a vs weak height e).
  
  We choose the 5D vector with minimum total |exponent| (Occam).

  Ratio                           Searched 5D                Consistent 5D  Same?
  --------------------------------------------------------------------------------
  m_mu/m_e                   (+0,-4,+1,+3,+0)             (+0,-4,+1,+3,+0)    YES
  m_tau/m_mu                 (+0,-3,-1,+0,+1)             (-1,-3,-1,+0,+0) differ
  m_tau/m_e                  (-2,-4,+0,+4,+0)             (-1,-7,+0,+3,+0) differ
  m_c/m_u                    (-2,-1,+0,+4,+0)             (-2,-1,+0,+4,+0)    YES
  m_t/m_c                    (-1,+0,-1,+1,+4)             (-5,+0,-1,+1,+0) differ
  m_s/m_d                    (-2,+1,+0,+2,+0)             (-2,+1,+0,+2,+0)    YES
  m_b/m_s                    (+0,+1,+0,+4,+0)             (+0,+1,+0,+4,+0)    YES
  m_t/m_b                    (-3,+2,-1,+1,+2)             (-4,-1,-1,+0,+0) differ
  m_d/m_u                    (-1,-2,+0,-1,+0)             (-1,-2,+0,-1,+0)    YES
  m_b/m_tau                  (+0,+1,-1,+1,+0)             (+0,+1,-1,+1,+0)    YES
  m_t/m_W                    (+0,+3,-1,+2,+0)             (+0,+3,-1,+2,+0)    YES
  m_H/m_W                    (+0,+1,-1,+0,+1)             (-1,+1,-1,+0,+0) differ
  m_H/m_Z                    (+1,+2,+0,+2,+0)             (+0,+2,+0,+2,-1) differ
  m_H/m_t                    (-1,-2,+0,-2,+0)             (-1,-2,+0,-2,+0)    YES

  8/14 vectors are IDENTICAL to individual search.
  6 vectors DIFFER — these are the transitivity corrections.


§6. QUANTUM NUMBER REGRESSION — Can we derive coordinates from (gen, I₃, Y, Nc)?
------------------------------------------------------------------------

  For fermions only, test the ansatz:
    v_i(p,q,r) = α₀ + α₁·gen + α₂·gen² + α₃·I₃ + α₄·Y + α₅·Nc

  If this fits well, the lattice coordinates are DERIVED from quantum numbers,
  not searched. Each coordinate p,q,r would be a known function of the particle.

  Coord      const      gen     gen²       I₃        Y       Nc       R²      RMSE
  ----------------------------------------------------------------------------
  p         -0.359   -0.583   +0.583   +0.833   +0.568   -0.045   0.6003    1.2886
  q         -0.001   -2.333   +0.333   +0.333   -0.768   -1.155   0.5561    1.5316
  r         -1.392   +5.000   -0.667   +1.333   +1.298   -1.532   0.8789    0.7536

  Predicted vs Actual coordinates:
  Particle    p_act   p_pred    q_act   q_pred    r_act   r_pred
  --------------------------------------------------------
  e             0.0    -1.39         0    -2.56         0    -0.56
  mu           -0.5    -0.22        -4    -3.89         3     2.44
  tau           1.0     2.11        -7    -4.56         3     4.11
  u            -1.5     0.11        -6    -5.56        -1    -0.56
  c             0.5     1.28        -7    -6.89         3     2.44
  t             6.0     3.61        -7    -7.56         4     4.11
  d            -0.5    -0.72        -8    -5.89        -2    -1.89
  s             1.5     0.44        -7    -7.22         0     1.11
  b             1.5     2.78        -6    -7.89         4     2.78


§7. SECTOR ANALYSIS — Generation stepping in the consistent lattice
------------------------------------------------------------------------

  For each sector, compute the lattice STEP between consecutive generations:
    Δv = v(gen+1) - v(gen)

  If the lattice is real, Δv should have a systematic sector-dependent pattern.


  Charged leptons (I₃=-½, Nc=1):
    mu/e: Δv = (-0.5, -4, +3)  → ratio = 206.008
    tau/mu: Δv = (+1.5, -3, +0)  → ratio = 16.839
    Acceleration γ = Δ(2→3)/Δ(1→2): p=-3.000, q=0.750, r=0.000

  Up quarks (I₃=+½, Nc=3):
    c/u: Δv = (+2.0, -1, +4)  → ratio = 587.224
    t/c: Δv = (+5.5, +0, +1)  → ratio = 135.765
    Acceleration γ = Δ(2→3)/Δ(1→2): p=2.750, q=-0.000, r=0.250

  Down quarks (I₃=-½, Nc=3):
    s/d: Δv = (+2.0, +1, +2)  → ratio = 19.863
    b/s: Δv = (+0.0, +1, +4)  → ratio = 44.692
    Acceleration γ = Δ(2→3)/Δ(1→2): p=0.000, q=1.000, r=2.000


§8. SYNTHESIS
========================================================================

  1. TRANSITIVITY TEST:
     Consistent lattice mean error: 0.3075%
     Individual search mean error:  0.2956%
     Ratio: 1.04×

     → The mass lattice IS consistent. Forcing transitivity does NOT
       degrade accuracy significantly. The 15 formulas ARE connected
       through a common set of particle coordinates.

  2. VECTOR CHANGES:
     8/14 vectors unchanged by transitivity
     6/14 vectors modified to enforce consistency
     Modified vectors are the CORRECT ones (if the lattice is real).

  3. QUANTUM NUMBER DERIVABILITY:
     p: R² = 0.6003 (moderate)
     q: R² = 0.5561 (moderate)
     r: R² = 0.8789 (good)

  If ALL R² > 0.95, the lattice coordinates ARE derivable from quantum
  numbers, and the mass formulas can be DERIVED (not searched).
  If not, additional structure (non-linear terms) is needed.


§9. PREDICTIONS — Ratios NOT used in the fit
------------------------------------------------------------------------
  Ratio           Exp.value    Lattice   Error%
  ----------------------------------------------
  m_t/m_e        337574.0792 333025.1511   1.348%
  m_t/m_u        79861.1111 79724.1879   0.171%
  m_b/m_d          889.3617   887.7073   0.186%
  m_c/m_s           13.5974    13.5000   0.717%
  m_c/m_d          270.2128   268.1498   0.763%
  m_b/m_u         1935.1852  1944.0000   0.456%
  m_tau/m_d        378.0553   379.2210   0.308%
  m_mu/m_u          48.9159    49.3171   0.820%
  m_W/m_Z            0.8815     0.8778   0.413%
  m_Z/m_t            0.5286     0.5329   0.801%
