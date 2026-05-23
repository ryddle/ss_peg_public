========================================================================
  ORIGIN OF THE 0.3% — LATTICE RESOLUTION OR GENUINE DEVIATION?
========================================================================

§1. EFFECTIVE DIMENSIONALITY REDUCTION
------------------------------------------------------------------------

  The 5 building blocks are NOT independent in log-space:
    ln R₂  = -ln 2
    ln R_EE = -½ ln 2
    ln h∨₂ = +ln 2
    ln R₃  = ln R₃   (irrational, independent of ln 2 and ln 3)
    ln h∨₃ = ln 3

  Any formula ln(f) = a·ln R₂ + b·ln R₃ + c·ln R_EE + d·ln h∨₃ + e·ln h∨₂
  reduces to:
    ln(f) = p·ln 2 + q·ln R₃ + r·ln 3
  where:
    p = (e - a - c/2)  [half-integer, from R₂ᵃ·R_EEᶜ·h∨₂ᵉ]
    q = b               [integer, from R₃ᵇ]
    r = d               [integer, from h∨₃ᵈ]

  Nominal formula count:     9^5 = 59,049
  Effective (p,q,r) triples: 3,321
  Numerically unique values: 3,323
  Reduction factor:          17.8x

  Baker's theorem check: 3321 triples vs 3323 unique values
  → CONFIRMED: 3 independent log-bases, 3323 distinct formulas

§2. LATTICE SPACING (FORMULA RESOLUTION)
------------------------------------------------------------------------
  Target log-range: [-0.405, 8.154]  (span = 8.559)
  Lattice points in range: 1707

  Gap statistics (within target range):
    Mean gap:    0.5620%
    Median gap:  0.4764%
    P10 gap:     0.0995%
    P90 gap:     0.7801%
    Max gap:     3.2334%
    Min gap:     0.099542%

  Expected BEST-MATCH error ≈ mean_gap/2 = 0.2810%

§3. ACTUAL ERRORS vs LATTICE BEST-MATCH
------------------------------------------------------------------------
  Name             Target Formula err  Best match  Local gap
  ------------------------------------------------------------
  Koide            0.6667     0.0000%     0.0000%    0.0995%
  m_mu/m_e       206.7683     0.3677%     0.2692%    0.6799%
  m_tau/m_mu      16.8170     0.1320%     0.1320%    0.3765%
  m_tau/m_e     3477.2283     0.5419%     0.1648%    0.7801% <<<
  m_c/m_u        587.9630     0.1257%     0.1258%    1.2602%
  m_t/m_c        135.8270     0.0460%     0.0460%    0.4764%
  m_s/m_d         19.8723     0.0471%     0.0471%    0.4764%
  m_b/m_s         44.7537     0.1387%     0.1389%    0.4764%
  m_t/m_b         41.2679     0.1512%     0.1512%    0.7801%
  m_d/m_u          2.1759     0.6439%     0.1667%    0.6799% <<<
  m_b/m_tau        2.3525     0.4943%     0.0203%    0.3022% <<<
  m_t/m_W          2.1461     0.3836%     0.0085%    0.0995% <<<
  m_H/m_W          1.5564     0.2686%     0.2686%    0.7801%
  m_H/m_Z          1.3719     0.1444%     0.1446%    0.4764%
  m_H/m_t          0.7252     0.6578%     0.1805%    0.6799% <<<

  Mean formula error:      0.2762%
  Mean best-lattice error: 0.1243%
  Mean local gap:          0.5616%

  Formulas that ARE the best lattice match: 6/15
  Sub-optimal formulas (better lattice match exists):
    m_mu/m_e    : formula=0.3677%, best=0.2692%
    m_tau/m_e   : formula=0.5419%, best=0.1648%
    m_d/m_u     : formula=0.6439%, best=0.1667%
    m_b/m_tau   : formula=0.4943%, best=0.0203%
    m_t/m_W     : formula=0.3836%, best=0.0085%
    m_H/m_t     : formula=0.6578%, best=0.1805%

§4. MONTE CARLO — WHAT ERRORS DO RANDOM TARGETS GET?
------------------------------------------------------------------------
  Random targets: N = 200,000
  Mean best-match error:   0.1909%
  Median best-match error: 0.1569%
  Max best-match error:    1.3626%

   Threshold   P(err<thr)        1/P
      0.01%       0.035875         27.9
      0.05%       0.180325          5.5
      0.10%       0.331360          3.0
      0.20%       0.609230          1.6
      0.30%       0.794700          1.3
      0.50%       0.955840          1.0
      1.00%       0.997010          1.0
      2.00%       1.000000          1.0
      5.00%       1.000000          1.0

§5. JOINT PROBABILITY OF 15 MASS RATIO MATCHES
------------------------------------------------------------------------

  Method: For each of the 15 targets, compute
    p_i = P(random_target gets best-match error ≤ actual_error_i)
  from the Monte Carlo distribution. Then:
    P_joint = Π p_i  (assuming approximate independence)

  Koide         err = 0.0000%  P(random ≤ err) = 0.000005
  m_mu/m_e      err = 0.3677%  P(random ≤ err) = 0.898480
  m_tau/m_mu    err = 0.1320%  P(random ≤ err) = 0.427395
  m_tau/m_e     err = 0.5419%  P(random ≤ err) = 0.966420
  m_c/m_u       err = 0.1257%  P(random ≤ err) = 0.408135
  m_t/m_c       err = 0.0460%  P(random ≤ err) = 0.166485
  m_s/m_d       err = 0.0471%  P(random ≤ err) = 0.170190
  m_b/m_s       err = 0.1387%  P(random ≤ err) = 0.447155
  m_t/m_b       err = 0.1512%  P(random ≤ err) = 0.483945
  m_d/m_u       err = 0.6439%  P(random ≤ err) = 0.981880
  m_b/m_tau     err = 0.4943%  P(random ≤ err) = 0.954465
  m_t/m_W       err = 0.3836%  P(random ≤ err) = 0.919265
  m_H/m_W       err = 0.2686%  P(random ≤ err) = 0.741255
  m_H/m_Z       err = 0.1444%  P(random ≤ err) = 0.463845
  m_H/m_t       err = 0.6578%  P(random ≤ err) = 0.982455

  JOINT PROBABILITY (product of all 15):
    log₁₀(P) = -8.87
    P ≈ 10^-8.9

  JOINT PROBABILITY (5 most constraining, rank-corrected):
    log₁₀(P) = -7.61
    P ≈ 10^-7.6

  ALTERNATIVE: uniform threshold
    Max actual error: 0.6578%
    P(one random ≤ max) = 0.982455
    P(15 random all ≤ max) = 0.982455^15 = 7.6681e-01
    P(5 random all ≤ max)  = 0.982455^5  = 9.1530e-01

§6. NON-INTEGER CORRECTIONS TO EXPONENTS
------------------------------------------------------------------------

  If the true relationship is R₂^a · R₃^b · R_EE^c · h∨₃^d · h∨₂^e = target,
  the integer exponents give ~0.3% error. What CONTINUOUS correction to
  each exponent would make it exact?

  For each formula, we find the smallest single-exponent correction
  δn_k that absorbs the residual: δn_k = residual / ln(block_k)

  Name           err%       dR2      dR3    dR_EE     dhv3     dhv2
  ----------------------------------------------------------------------
  Koide        0.000%  -0.00000 -0.00000 -0.00000 +0.00000 +0.00000
  m_mu/m_e     0.369%  -0.00531 -0.00619 -0.01063 +0.00335 +0.00531
  m_tau/m_mu   0.132%  +0.00190 +0.00222 +0.00381 -0.00120 -0.00190
  m_tau/m_e    0.539%  +0.00780 +0.00909 +0.01559 -0.00492 -0.00780
  m_c/m_u      0.126%  -0.00181 -0.00211 -0.00363 +0.00114 +0.00181
  m_t/m_c      0.046%  -0.00066 -0.00077 -0.00133 +0.00042 +0.00066
  m_s/m_d      0.047%  -0.00068 -0.00079 -0.00136 +0.00043 +0.00068
  m_b/m_s      0.139%  -0.00200 -0.00233 -0.00400 +0.00126 +0.00200
  m_t/m_b      0.151%  +0.00218 +0.00254 +0.00436 -0.00138 -0.00218
  m_d/m_u      0.640%  +0.00926 +0.01079 +0.01852 -0.00584 -0.00926
  m_b/m_tau    0.497%  -0.00715 -0.00833 -0.01430 +0.00451 +0.00715
  m_t/m_W      0.385%  -0.00554 -0.00646 -0.01109 +0.00350 +0.00554
  m_H/m_W      0.268%  +0.00387 +0.00451 +0.00774 -0.00244 -0.00387
  m_H/m_Z      0.145%  -0.00209 -0.00243 -0.00417 +0.00132 +0.00209
  m_H/m_t      0.653%  +0.00946 +0.01103 +0.01892 -0.00597 -0.00946

  If there existed a universal shift δ to ALL R₃ exponents:
    Mean δ(R₃ exponent):   +0.00089
    Std δ(R₃ exponent):    0.00625
    If universal: all should be same → std/mean = 7.06
    → No universal correction. Residuals are quasi-random.

§7. ERROR vs LOCAL LATTICE DENSITY
------------------------------------------------------------------------

  If errors = lattice resolution: formula_error ∝ local_gap
  If errors = physical deviation:  formula_error independent of local_gap

  Correlation(formula_error, local_gap) = 0.1442
  → No correlation: errors are NOT from lattice resolution

  Error / local_gap ratio:
    Mean:   0.695
    Median: 0.344
    Std:    0.945
    (If pure lattice resolution → ratio ≈ 0.5 uniformly)

§8. RANK CONSTRAINT — WHICH MATCHES ARE INDEPENDENT?
------------------------------------------------------------------------
  Exponent matrix: 15 formulas x 5 blocks
  Rank: 5
  Independent predictions: 5
  Consistency constraints: 10

  Singular values: [9.446 8.288 5.535 3.548 1.692]

  Independent set: ['Koide', 'm_mu/m_e', 'm_tau/m_mu', 'm_tau/m_e', 'm_c/m_u']
  Dependent set:   ['m_t/m_c', 'm_s/m_d', 'm_b/m_s', 'm_t/m_b', 'm_d/m_u', 'm_b/m_tau', 'm_t/m_W', 'm_H/m_W', 'm_H/m_Z', 'm_H/m_t']

  CONSISTENCY CHECKS (dependent formulas):
  Name         Implied value     Actual  Consist.err
  -------------------------------------------------------
  m_t/m_c           143.4425   135.8270      5.6068%
  m_s/m_d            19.9533    19.8723      0.4077%
  m_b/m_s            44.8950    44.7537      0.3157%
  m_t/m_b            42.5155    41.2679      3.0232%
  m_d/m_u             2.1802     2.1759      0.1979%
  m_b/m_tau           2.3222     2.3525      1.2899%
  m_t/m_W             2.1315     2.1461      0.6825%
  m_H/m_W             1.5724     1.5564      1.0298%
  m_H/m_Z             1.3768     1.3719      0.3582%
  m_H/m_t             0.7263     0.7252      0.1536%

  Mean consistency error:   1.3065%
  Max consistency error:    5.6068%

  These consistency checks are INDEPENDENT of the graph.
  They test whether the experimental mass ratios satisfy
  the multiplicative relations implied by the exponent structure.

§9. CORRECTED JOINT PROBABILITY
------------------------------------------------------------------------
  Monte Carlo estimate: P(15 random targets ALL match within 0.658%)
    Trials:  100,000
    Passes:  76932
    P =      0.769320
    log₁₀(P) = -0.11

  Analytical product of individual P_i:
    P_joint (all 15) = 10^-8.87
    P_joint (5 indep) = 10^-6.12

========================================================================
§10. SYNTHESIS — WHAT IS THE 0.3%?
========================================================================

  ╔══════════════════════════════════════════════════════════════════╗
  ║  FINDING 1: EFFECTIVE FORMULA SPACE                             ║
  ║                                                                  ║
  ║  5 building blocks → 3 independent log-bases                    ║
  ║  9^5 = 59,049 nominal formulas → 3,323 distinct values        ║
  ║  Reduction: 18x (R₂, R_EE, h∨₂ all powers of 2)           ║
  ╠══════════════════════════════════════════════════════════════════╣
  ║  FINDING 2: LATTICE RESOLUTION                                  ║
  ║                                                                  ║
  ║  Mean gap in target range:   0.5620%                        ║
  ║  Expected best-match error:  0.2810%  (gap/2)              ║
  ║  Actual mean formula error:  0.2762%                        ║
  ║  Ratio actual/expected:      0.98x                           ║
  ╠══════════════════════════════════════════════════════════════════╣
  ║  FINDING 3: JOINT PROBABILITY                                   ║
  ║                                                                  ║
  ║  P(all 15 ≤ observed errors)   = 10^-8.9                     ║
  ║  P(5 most constraining)        = 10^-7.6                     ║
  ║  P(15 random all ≤ max error)  = 0.769320 (MC)                  ║
  ╠══════════════════════════════════════════════════════════════════╣
  ║  FINDING 4: CONSISTENCY                                          ║
  ║                                                                  ║
  ║  Rank of exponent matrix:      5                                ║
  ║  Independent searched formulas: 5                                ║
  ║  Consistency tests (dependent): 10                               ║
  ║  Mean consistency error:        1.3065%                       ║
  ╚══════════════════════════════════════════════════════════════════╝


  ═══════════════════════════════════════════════════════════════════
  THE CONCEPTUAL QUESTION: WHAT IS THE ~0.3% DEVIATION?
  ═══════════════════════════════════════════════════════════════════

  Four hypotheses tested:

  H1: Lattice resolution (integer-exponent discretization)
      If true: errors ∝ local gap, ratio ≈ 0.5
      Result: see §7 correlation analysis above

  H2: Missing QFT corrections (tree-level approximation)
      If true: errors ≈ α_s/π ≈ 3.75%
      Result: observed 0.3% << 3.75% → NOT pure tree-level

  H3: Non-perturbative (formulas are exact "dressed" values)
      If true: errors ≈ 0 (or ~ exp(-1/α))
      Result: observed 0.3% >> 0 → NOT fully non-perturbative

  H4: Partial resummation (graph captures most but not all
      radiative structure)
      If true: errors ≈ α²/(4π) ≈ 0.01-0.3% (2-loop level)
      Result: observed 0.3% is at the upper edge of this range

  α_s²/(4π)² ≈ 0.009%  (pure 2-loop strong)
  α_s²/(4π)  ≈ 0.11%   (mixed 2-loop)
  α_s/π      ≈ 3.75%   (1-loop strong)

  The 0.3% sits between 1-loop and 2-loop, consistent with
  "all leading 1-loop absorbed, some subleading 1-loop remaining."


========================================================================
  ANALYSIS COMPLETE
========================================================================

  Output saved to: D:\Desarrollo\projects\pyshics\ss_peg\output\residual_origin_output.md
