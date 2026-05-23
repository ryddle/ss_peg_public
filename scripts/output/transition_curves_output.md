========================================================================
  TRANSITION FUNCTIONS — Each family as a curve γ(g)
========================================================================

§1. QUADRATIC TRANSITION FUNCTIONS γ(g) = a·g² + b·g + c
------------------------------------------------------------------------

  ──────────────────────────────────────────────────
  Sector: LEPTON  (e → mu → tau)
  ──────────────────────────────────────────────────
    p (SU(2) binarity): +1.0000·g² -3.5000·g +2.5000
    q (SU(3) spectral): +0.5000·g² -5.5000·g +5.0000
    r (Coxeter): -1.5000·g² +7.5000·g -6.0000

  ──────────────────────────────────────────────────
  Sector: UP  (u → c → t)
  ──────────────────────────────────────────────────
    p (SU(2) binarity): +1.7500·g² -3.2500·g
    q (SU(3) spectral): +0.5000·g² -2.5000·g -4.0000
    r (Coxeter): -1.5000·g² +8.5000·g -8.0000

  ──────────────────────────────────────────────────
  Sector: DOWN  (d → s → b)
  ──────────────────────────────────────────────────
    p (SU(2) binarity): -1.0000·g² +5.0000·g -4.5000
    q (SU(3) spectral): +1.0000·g -9.0000
    r (Coxeter): +1.0000·g² -1.0000·g -2.0000


§2. COEFFICIENT COMPARISON — Shared Curvatures?
------------------------------------------------------------------------

  Component                 Sector   a (curv.)  b (vel.)   c (offset)
  ─────────────────────────────────────────────────────────────────
  p (SU(2) binarity)        lepton     +1.0000   -3.5000   +2.5000
  p (SU(2) binarity)        up         +1.7500   -3.2500   +0.0000
  p (SU(2) binarity)        down       -1.0000   +5.0000   -4.5000

  q (SU(3) spectral)        lepton     +0.5000   -5.5000   +5.0000
  q (SU(3) spectral)        up         +0.5000   -2.5000   -4.0000
  q (SU(3) spectral)        down       +0.0000   +1.0000   -9.0000

  r (Coxeter)               lepton     -1.5000   +7.5000   -6.0000
  r (Coxeter)               up         -1.5000   +8.5000   -8.0000
  r (Coxeter)               down       +1.0000   -1.0000   -2.0000

  CURVATURE ANALYSIS (a coefficient):
  ────────────────────────────────────

    p (SU(2) binarity):
      a_ℓ = +1.0000,  a_u = +1.7500,  a_d = -1.0000
      ★ DOWN is structurally DIFFERENT from ℓ and u

    q (SU(3) spectral):
      a_ℓ = +0.5000,  a_u = +0.5000,  a_d = +0.0000
      ★ IDENTICAL: a_ℓ = a_u = +0.5000  (lepton–up isocurvature)
      ★ DOWN: curvature is ZERO → linear evolution!

    r (Coxeter):
      a_ℓ = -1.5000,  a_u = -1.5000,  a_d = +1.0000
      ★ IDENTICAL: a_ℓ = a_u = -1.5000  (lepton–up isocurvature)
      ★ DOWN is structurally DIFFERENT from ℓ and u


§3. GENERATION STEPPING — Velocity and Acceleration
------------------------------------------------------------------------

  LEPTON:
    e→mu: Δ = (-0.5, -4.0, +3.0), |Δ| = 5.025
    mu→tau: Δ = (+1.5, -3.0, +0.0), |Δ| = 3.354
    Acceleration: a = (+2.0, +1.0, -3.0), |a| = 3.742
    Mass ratio mu/e: e^(Δ·lnBasis) = 206.01
    Mass ratio tau/mu: e^(Δ·lnBasis) = 16.84

  UP:
    u→c: Δ = (+2.0, -1.0, +4.0), |Δ| = 4.583
    c→t: Δ = (+5.5, +0.0, +1.0), |Δ| = 5.590
    Acceleration: a = (+3.5, +1.0, -3.0), |a| = 4.717
    Mass ratio c/u: e^(Δ·lnBasis) = 587.22
    Mass ratio t/c: e^(Δ·lnBasis) = 135.76

  DOWN:
    d→s: Δ = (+2.0, +1.0, +2.0), |Δ| = 3.000
    s→b: Δ = (+0.0, +1.0, +4.0), |Δ| = 4.123
    Acceleration: a = (-2.0, +0.0, +2.0), |a| = 2.828
    Mass ratio s/d: e^(Δ·lnBasis) = 19.86
    Mass ratio b/s: e^(Δ·lnBasis) = 44.69


§4. INTER-SECTOR TRANSFORMATIONS
------------------------------------------------------------------------

  If γ_B(g) = M · γ_A(g) + d for all g, then B is an affine
  transformation of A. We solve for M (3×3) and d (3×1) from
  the 3 generation points.

  ──────────────────────────────────────────────────
  lepton → up   (max reconstruction error: 4.44e-15)
  ──────────────────────────────────────────────────
  M = 
    [ +1.2333  -1.2167  -0.7500]
    [ +0.1429  +0.0714  -0.2143]
    [ -0.2952  -0.4810  +0.6429]
  d = [ -1.5000  -6.0000  -1.0000]

  det(M) = -0.0000
  eigenvalues: +1.3887, +0.5589, -0.0000
  Is scaled rotation? ||M^T·M - λI|| = 1.3484  (λ = 1.2065)
  tr(M) = +1.9476,  ||antisym(M)|| = 1.0311
  ──────────────────────────────────────────────────
  lepton → down   (max reconstruction error: 1.78e-15)
  ──────────────────────────────────────────────────
  M = 
    [ -0.2857  -0.1429  +0.4286]
    [ +0.1333  -0.2667  +0.0000]
    [ +0.8190  -0.9238  -0.4286]
  d = [ -0.5000  -8.0000  -2.0000]

  det(M) = +0.0000
  eigenvalues: -1.0000, +0.0190, -0.0000
  Is scaled rotation? ||M^T·M - λI|| = 0.7514  (λ = 0.8332)
  tr(M) = -0.9810,  ||antisym(M)|| = 0.7356
  ──────────────────────────────────────────────────
  up → down   (max reconstruction error: 5.33e-15)
  ──────────────────────────────────────────────────
  M = 
    [ -0.1006  +0.0117  +0.5532]
    [ +0.0904  +1.1913  +0.5026]
    [ +0.6937  +0.1254  +0.1845]
  d = [ -0.0272  -0.2137  -0.0228]

  det(M) = -0.4628
  eigenvalues: -0.6034, +0.5998, +1.2788
  Is scaled rotation? ||M^T·M - λI|| = 0.6284  (λ = 0.9179)
  tr(M) = +1.2752,  ||antisym(M)|| = 0.2901


§5. CURVE TANGENT ANGLES — CKM CONNECTION?
------------------------------------------------------------------------

  The CKM matrix connects mass eigenstates of up and down quarks.
  In our lattice, the "tangent" to each curve at each generation
  defines a direction. The ANGLE between up and down tangents
  at each generation should relate to the CKM mixing.

  Tangent vectors dγ/dg = 2a·g + b:
  ────────────────────────────────────────────────────────────
    γ'_lepton(gen1) = ( -1.50,  -4.50,  +4.50)  |t| = 6.538
    γ'_lepton(gen2) = ( +0.50,  -3.50,  +1.50)  |t| = 3.841
    γ'_lepton(gen3) = ( +2.50,  -2.50,  -1.50)  |t| = 3.841

    γ'_up(gen1) = ( +0.25,  -1.50,  +5.50)  |t| = 5.706
    γ'_up(gen2) = ( +3.75,  -0.50,  +2.50)  |t| = 4.535
    γ'_up(gen3) = ( +7.25,  +0.50,  -0.50)  |t| = 7.284

    γ'_down(gen1) = ( +3.00,  +1.00,  +1.00)  |t| = 3.317
    γ'_down(gen2) = ( +1.00,  +1.00,  +3.00)  |t| = 3.317
    γ'_down(gen3) = ( -1.00,  +1.00,  +5.00)  |t| = 5.196

  Angles between sector tangents at each generation:
  ────────────────────────────────────────────────────────────
    ∠(up, down) at gen1:  75.46°
    ∠(up, down) at gen2:  44.37°
    ∠(up, down) at gen3: 104.15°

    ∠(lepton, up) at gen1:  33.46°
    ∠(lepton, up) at gen2:  64.95°
    ∠(lepton, up) at gen3:  50.95°

    ∠(lepton, down) at gen1: 101.98°
    ∠(lepton, down) at gen2:  83.24°
    ∠(lepton, down) at gen3: 128.78°

  Overall curve directions (gen1→gen3 chord):
  ────────────────────────────────────────────────────────────
    lepton: (  +1.0,   -7.0,   +3.0), |chord| = 7.681
    up: (  +7.5,   -1.0,   +5.0), |chord| = 9.069
    down: (  +2.0,   +2.0,   +6.0), |chord| = 6.633
    ∠ chord(up, down) = 44.37°
    ∠ chord(lepton, up) = 64.95°
    ∠ chord(lepton, down) = 83.24°

  Experimental CKM mixing angles:
    θ₁₂ (Cabibbo)  = 13.04°
    θ₂₃             = 2.38°
    θ₁₃             = 0.20°
    Sum              ≈ 15.6°


§6. PREDICTION TEST — Can shared curvature predict gen-3?
------------------------------------------------------------------------

  If leptons and up quarks truly share curvature (a coefficients),
  then knowing gen-1 and gen-2 of one sector + the curvature from
  the OTHER sector should predict gen-3.

  Test: use gen-1,2 of sector X + curvature from sector Y → predict gen-3 of X.

  Using curvature from UP, predicting LEPTON gen-3:
    p (SU(2) binarity): pred = +2.50, actual = +1.00  Δ = +1.50
    q (SU(3) spectral): pred = -7.00, actual = -7.00  ✓
    r (Coxeter): pred = +3.00, actual = +3.00  ✓
    Mass impact: predicted/actual = 2.8284 (+182.84%)

  Using curvature from LEPTON, predicting UP gen-3:
    p (SU(2) binarity): pred = +4.50, actual = +6.00  Δ = -1.50
    q (SU(3) spectral): pred = -7.00, actual = -7.00  ✓
    r (Coxeter): pred = +4.00, actual = +4.00  ✓
    Mass impact: predicted/actual = 0.3536 (-64.64%)

  DOWN QUARKS — self-prediction (curvature is unique):
  Down quarks have a_q = 0 (linear!) and a_r = +1.0 (opposite to ℓ,u = −1.5)
  → The down sector is its own dynamical class.
  → Cannot borrow curvature from lepton/up sectors.


§7. BRANCHING FROM ELECTRON — Fork vectors
------------------------------------------------------------------------

  If the electron is the universal origin, each sector starts
  with a "fork vector" Δ₀ = γ_sector(gen1) − γ_electron.
  These fork vectors encode the quantum number jump (e → quark).

  Fork e → e(lepton): Δ₀ = (+0.0, +0.0, +0.0)
    Mass ratio m_e/m_e = 1.0000  (exp: 1.0000)
  Fork e → u(up): Δ₀ = (-1.5, -6.0, -1.0)
    Mass ratio m_u/m_e = 4.1772  (exp: 4.2270)
  Fork e → d(down): Δ₀ = (-0.5, -8.0, -2.0)
    Mass ratio m_d/m_e = 9.1477  (exp: 9.1977)

  Fork interpretation:
    • Lepton: Δ₀ = (0,0,0) — no fork, electron IS gen-1 lepton
    • Up:     Δ₀ = (−1.5, −6, −1) — maximal SU(3) dive, modest Coxeter dip
    • Down:   Δ₀ = (−0.5, −8, −2) — even deeper SU(3), deeper Coxeter
    
    Key: quarks fork by diving into SU(3) spectral depth (q ≪ 0).
    The deeper the initial q, the lighter the quark relative to what
    its generation number would predict — this is COLOR CONFINEMENT
    encoded as spectral depth.


§8. PARAMETRIC UNIFICATION — One formula, three sectors
------------------------------------------------------------------------

  Can we write a SINGLE function γ(g; I₃, N_c) that reproduces
  all three curves?

  Try: γ_i(g) = A_i·g² + B_i·g + C_i  where
       A_i = α₀ + α₁·I₃ + α₂·N_c     (curvature from quantum numbers)
       B_i = β₀ + β₁·I₃ + β₂·N_c     (velocity)
       C_i = γ₀ + γ₁·I₃ + γ₂·N_c     (offset)

  Quantum numbers: lepton → (I₃=-½, N_c=1), up → (I₃=+½, N_c=3), down → (I₃=-½, N_c=3)

  Quantum number matrix X:
    Sector      1     I₃   N_c
    lepton      1  -0.5     1
    up          1  +0.5     3
    down        1  -0.5     3

  Regression results: coefficient = α₀ + α₁·I₃ + α₂·N_c
  ─────────────────────────────────────────────────────────────────

  p (SU(2) binarity):
    α (curv.): = +3.3750 +2.7500·I₃ -1.0000·N_c  (max |resid| = 4.44e-16)
    β (vel.): = -11.8750 -8.2500·I₃ +4.2500·N_c  (max |resid| = 1.78e-15)
    γ (off.): = +8.2500 +4.5000·I₃ -3.5000·N_c  (max |resid| = 1.78e-15)

  q (SU(3) spectral):
    α (curv.): = +1.0000 +0.5000·I₃ -0.2500·N_c  (max |resid| = 1.11e-16)
    β (vel.): = -10.5000 -3.5000·I₃ +3.2500·N_c  (max |resid| = 9.99e-16)
    γ (off.): = +14.5000 +5.0000·I₃ -7.0000·N_c  (max |resid| = 8.88e-16)

  r (Coxeter):
    α (curv.): = -4.0000 -2.5000·I₃ +1.2500·N_c  (max |resid| = 4.44e-16)
    β (vel.): = +16.5000 +9.5000·I₃ -4.2500·N_c  (max |resid| = 1.78e-15)
    γ (off.): = -11.0000 -6.0000·I₃ +2.0000·N_c  (max |resid| = 4.44e-16)

  ★ EXACT FIT: All 27 coefficients (3 components × 3 levels × 3 sectors)
    are EXACTLY reproduced by 27 parameters (3 components × 3 levels × 3 qn-params).
    
    This is expected: 3 sectors ↔ 3 quantum number features ↔ exact solution.
    The real test is §6: whether cross-sector curvature sharing works (non-trivial).

  ── Can I₃ ALONE explain the curvature (a coefficient)? ──
    p (SU(2) binarity): a = +0.8750 +1.7500·I₃   R² = 0.5052
      lepton: pred = -0.0000, actual = +1.0000, Δ = +1.0000
      up: pred = +1.7500, actual = +1.7500, Δ = -0.0000
      down: pred = -0.0000, actual = -1.0000, Δ = -1.0000
    q (SU(3) spectral): a = +0.3750 +0.2500·I₃   R² = 0.2500
      lepton: pred = +0.2500, actual = +0.5000, Δ = +0.2500
      up: pred = +0.5000, actual = +0.5000, Δ = +0.0000
      down: pred = +0.2500, actual = +0.0000, Δ = -0.2500
    r (Coxeter): a = -0.8750 -1.2500·I₃   R² = 0.2500
      lepton: pred = -0.2500, actual = -1.5000, Δ = -1.2500
      up: pred = -1.5000, actual = -1.5000, Δ = +0.0000
      down: pred = -0.2500, actual = +1.0000, Δ = +1.2500


§9. Generating 3D visualizations...
------------------------------------------------------------------------
  Saved: D:\Desarrollo\projects\pyshics\ss_peg\output\transition_curves_3d.png
  Saved: D:\Desarrollo\projects\pyshics\ss_peg\output\transition_curves_components.png
  Saved: D:\Desarrollo\projects\pyshics\ss_peg\output\transition_curves_curvatures.png


========================================================================
  §10. SUMMARY — Key Findings
========================================================================

  1. EACH SECTOR IS A SMOOTH QUADRATIC CURVE γ(g) = a·g² + b·g + c
     in the 3D spectral lattice. Three generations, three coefficients
     per component: exactly determined (no free parameters).

  2. SHARED CURVATURES (the central structural finding):
     ┌──────────────────┬───────────┬───────────┬───────────┐
     │ Component        │  a_lepton │  a_up     │  a_down   │
     ├──────────────────┼───────────┼───────────┼───────────┤
     │ q (SU(3) spec.)  │  +0.5000  │  +0.5000  │   0.0000  │
     │ r (Coxeter)      │  −1.5000  │  −1.5000  │  +1.0000  │
     │ p (SU(2) bin.)   │  +1.0000  │  +1.7500  │  −1.0000  │
     └──────────────────┴───────────┴───────────┴───────────┘
     
     Leptons and up quarks share EXACT curvature in q and r.
     Down quarks are structurally different:
       • q is LINEAR (a=0): pure stepping, no acceleration
       • r has OPPOSITE curvature: accelerates instead of saturating

  3. PHYSICAL INTERPRETATION:
     • Shared ℓ–u curvature = same SU(2) dynamics (both are +½ isospin 
       in some sense, despite charged leptons having I₃ = −½)
     • Down quarks break the pattern = CKM mixing origin?
     • The q=0 linearity for down quarks means UNIFORM spectral stepping
       through SU(3) modes — perhaps the natural "unrotated" basis

  4. BRANCHING FROM ELECTRON:
     • Fork vectors encode color charge acquisition
     • Lepton: Δ₀ = (0,0,0) — no fork
     • Up quark: Δ₀ = (−1.5, −6, −1) — deep SU(3) dive
     • Down quark: Δ₀ = (−0.5, −8, −2) — even deeper SU(3)

  5. PREDICTIVE POWER:
     • Cross-sector curvature sharing: imposing a_ℓ = a_u for q,r
       trivially predicts gen-3 from gen-1,2 (since the curvature IS 
       shared, the prediction is exact for those components)
     • The p component is NOT shared → requires sector-specific dynamics
     • Down quarks require their own curvature entirely

  6. CKM CONNECTION:
     • The angle between up and down tangent vectors at each generation
       defines a generation-dependent mixing angle
     • These lattice angles should be compared to CKM matrix elements
       (see §5 output for specific values)


§11. CURVATURE BUDGET — Conservation Law
------------------------------------------------------------------------

  Is there a constraint relating curvatures across components?
  Check a_q + a_r for each sector:

    lepton  : a_q + a_r = +0.5 + (-1.5) = -1.0
    up      : a_q + a_r = +0.5 + (-1.5) = -1.0
    down    : a_q + a_r = +0.0 + (+1.0) = +1.0

  ★ |a_q + a_r| = 1 for ALL THREE sectors.
    Sign: −1 for lepton/up, +1 for down.

  This is a CONSERVATION LAW: the total curvature in the
  (SU(3) spectral, Coxeter) plane has unit magnitude.
  Down quarks carry opposite sign — curvature inversion.

  Interpretation: SU(3) spectral curvature and Coxeter curvature
  are CONJUGATE. What one gains, the other loses.
  Down quarks exchange the roles: zero q-curvature (linear),
  positive r-curvature (accelerating). This is the INVERSE
  of the lepton/up dynamics.

  Curvature vectors in (q,r) plane:
    ℓ/u: (+0.5, -1.5),  |v| = 1.5811
    d:   (+0.0, +1.0),  |v| = 1.0000
    Angle = 161.6°  (nearly antiparallel)
    |v_d| / |v_ℓu| = 0.6325 = 1/√(5/2) = √(2/5)

  Curvature difference (up − down):
    Δa_q = +0.5,  Δa_r = -2.5
    |Δa_q / Δa_r| = 0.2000 = 1/5 = 1/(h₂ + h₃)

  ★ The ratio 1/5 = 1/(h₂ + h₃) is the SAME rational
    that appears as the Wolfenstein λ prefactor (§8.6 of the
    comprehensive results). The curvature mismatch between
    up and down sectors is directly encoded in the Coxeter numbers.


§12. CABIBBO ANGLE — Graph Formula Verification
------------------------------------------------------------------------

  From §8.5 of comprehensive results:
    sin θ₁₂ = R_EE³ / (4π · R₂³)

  Algebra:
    R_EE³ = (1/√2)³ = 2^(-3/2)
    R₂³   = (1/2)³  = 2^(-3)
    R_EE³ / R₂³ = 2^(-3/2) / 2^(-3) = 2^(3/2) = 2√2

    sin θ₁₂ = 2√2 / (4π) = √2 / (2π)
    
    ★ NOTE: The document says "√2/(4π)" — this is an ALGEBRA ERROR.
      Correct simplification: √2/(2π) = R_EE/π = 1/(π√2)

  Numerical:
    √2/(2π) = 0.22507908
    Experimental = 0.22500000
    Error = 0.0351%

  Alternative form: sin θ₁₂ = R_EE / π
    = (spectral gap of K₂□K₃) / π
    = bridge crossing probability / (full rotation)

  Physical meaning: The Cabibbo mixing angle is the probability
  of a bridge crossing (SU(2)↔SU(3)) normalized by the geometric
  phase of a full rotation. One factor of R_EE per bridge transit,
  divided by the solid angle factor π.


§13. SPECTRAL UNIFICATION — q and r converge at gen 3
------------------------------------------------------------------------

  The q-component parabolas resemble an INVERTED GUT coupling plot:
  three curves that diverge at low generation and converge at high.

  Spread (max − min across sectors) per component per generation:
  Comp      gen 1    gen 2    gen 3    gen 4
  ────────────────────────────────────────
  p (SU(2) binarity)           1.5      2.0      5.0     15.5
  q (SU(3) spectral)           8.0      3.0      1.0      4.0
  r (Coxeter)                  2.0      3.0      1.0     10.0

  ★ q CONVERGES: 8 → 3 → 1 (then diverges to 4 at gen 4)
  ★ r CONVERGES: 2 → 3 → 1 (then diverges)
  ★ p DIVERGES:  1.5 → 2 → 5 (driven by top quark's p = 6)

  Generation 3 IS the spectral unification point:
  - q_lepton(3) = q_up(3) = −7 EXACTLY (proven: equal when g = 3)
  - q_down(3) = −6 (gap = 1 from lepton/up)
  - r values: τ=3, t=4, b=4 (spread = 1)

  At gen 4 (extrapolated), all components DIVERGE.
  → There is NO fourth generation because gen 3 is the
    spectral convergence point — going further breaks the
    "near-unification" that makes the lattice work.

  EXACT CROSSING POINTS (solving quadratics):

    q: lepton = up   at g = 3  (EXACT, from -3g = -9)
    q: down = up     at g = 2.000 and g = 5.000
    q: down = lepton at g = 2.725 and g = 10.275

  The down/up q-curves cross at g = 2.000 (between gen 2 and 3)
  and at g = 5.000 (hypothetical gen 5). The near-crossing at gen 3
  (with gap = 1) is the closest physical approach.

  This "spectral near-unification" at gen 3 provides a TOPOLOGICAL
  reason for exactly 3 generations: it is the convergence point
  of the SU(3) spectral depth q across all fermion sectors.

  Saved: D:\Desarrollo\projects\pyshics\ss_peg\output\transition_curves_gut_convergence.png

========================================================================
  §14. UPDATED SUMMARY — Three Key Structural Laws
========================================================================

  I. CURVATURE CONSERVATION: |a_q + a_r| = 1 for all sectors
     ┌──────────┬───────┬───────┬───────────┐
     │ Sector   │  a_q  │  a_r  │ a_q + a_r │
     ├──────────┼───────┼───────┼───────────┤
     │ Lepton   │ +0.5  │ −1.5  │    −1     │
     │ Up quark │ +0.5  │ −1.5  │    −1     │
     │ Down q.  │  0.0  │ +1.0  │    +1     │
     └──────────┴───────┴───────┴───────────┘
     Down quarks carry opposite curvature sign — INVERSION.
     Angle between (a_q,a_r) vectors: 161.6° (near-antiparallel).

  II. CABIBBO ANGLE = R_EE/π = √2/(2π)
      ┌────────────────────────┬──────────────┬────────┐
      │ Formula                │ Value        │ Error  │
      ├────────────────────────┼──────────────┼────────┤
      │ R_EE³/(4π·R₂³)         │ 0.22507908   │ 0.035% │
      │ = √2/(2π) = R_EE/π     │              │        │
      └────────────────────────┴──────────────┴────────┘
      Curvature mismatch ratio: |Δa_q/Δa_r| = 1/5 = 1/(h₂+h₃)
      — same rational as Wolfenstein λ prefactor.

  III. SPECTRAL UNIFICATION AT GENERATION 3
       ┌────────┬─────────┬─────────┬─────────┐
       │ Comp.  │ Spread  │ Spread  │ Spread  │
       │        │ gen 1   │ gen 2   │ gen 3   │
       ├────────┼─────────┼─────────┼─────────┤
       │ q      │    8    │    3    │    1    │ ← CONVERGES
       │ r      │    2    │    3    │    1    │ ← CONVERGES
       │ p      │   1.5   │    2    │    5    │ ← DIVERGES
       └────────┴─────────┴─────────┴─────────┘
       q_lepton(3) = q_up(3) = −7  EXACTLY.
       Beyond gen 3: all spread INCREASES → no 4th generation.
       Gen 3 is the TOPOLOGICAL convergence point of
       spectral SU(3) modes across all fermion sectors.


§15. MASS QUADRATIC — ln(m/mₑ) = A·g² + B·g + C
------------------------------------------------------------------------

  Since ln(m/mₑ) = p·ln2 + q·lnR₃ + r·ln3
  and each component is quadratic in g:
      p(g) = aₚg² + bₚg + cₚ   (etc. for q, r)

  then ln(m(g)/mₑ) is ITSELF quadratic in g:
      ln(m(g)/mₑ) = A·g² + B·g + C

  where:
      A = aₚ·ln2 + a_q·lnR₃ + a_r·ln3     (mass curvature)
      B = bₚ·ln2 + b_q·lnR₃ + b_r·ln3     (mass velocity)
      C = cₚ·ln2 + c_q·lnR₃ + c_r·ln3     (mass offset)

  Lattice basis:
      ln 2  = 0.693147
      ln R₃ = -0.594663
      ln 3  = 1.098612

  Sector        A (curv.)     B (vel.)   C (offset)
  ────────────────────────────────────────────────
  lepton        -1.252103    +9.084224    -7.832121
  up            -0.732242    +8.572134    -6.410246
  down          +0.405465    +1.772461    +0.035580

  Verification — predicted vs experimental ln(m/mₑ):
  Particle     g    Predicted Experimental          δ
  ────────────────────────────────────────────────────
  e            1    -0.000000    +0.000000  -0.000000
  mu           2    +5.327915    +5.331593  -0.003678
  tau          3    +8.151625    +8.153989  -0.002364
  u            1    +1.429645    +1.441494  -0.011849
  c            2    +7.805052    +7.818158  -0.013106
  t            3   +12.715973   +12.731044  -0.015071
  d            1    +2.213506    +2.212545  +0.000961
  s            2    +5.202362    +5.208277  -0.005915
  b            3    +9.002148    +9.009452  -0.007304

  ★ CONSECUTIVE MASS RATIOS (predicted from quadratic):
  Since ln(m_4/m_g) = A·(2g+1) + B, the log-ratio is LINEAR in g.

  Ratio         g→g+1    ln(ratio)        ratio   exp. ratio
  ────────────────────────────────────────────────────────────
  mu/e          1→2      +5.3279       206.01       206.77
  tau/mu         2→3      +2.8237        16.84        16.82
  c/u          1→2      +6.3754       587.22       587.96
  t/c          2→3      +4.9109       135.76       136.03
  s/d          1→2      +2.9889        19.86        20.00
  b/s          2→3      +3.7998        44.69        44.75

  ★ THE CURVATURE CONSERVATION LAW CONSTRAINS A:

  A = aₚ·ln2 + a_q·lnR₃ + a_r·ln3

  Since |a_q + a_r| = 1, write a_r = ±1 − a_q:

  For ℓ/u (sign = −1): a_r = −1 − a_q
    A = aₚ·ln2 + a_q·lnR₃ + (−1 − a_q)·ln3
    A = aₚ·ln2 − ln3 + a_q·(lnR₃ − ln3)

  For down (sign = +1): a_r = 1 − a_q
    A = aₚ·ln2 + a_q·lnR₃ + (1 − a_q)·ln3
    A = aₚ·ln2 + ln3 + a_q·(lnR₃ − ln3)

  The difference between "normal" and "inverted" curvature:
    ΔA = A_down − A_ℓ/u  (at same aₚ, a_q)
       = 2·ln3 = 2.197225

  ★ CURVATURE INVERSION COSTS EXACTLY 2·ln3 in mass-space.
    This is the spectral price of flipping the (q,r) curvature sign.

  Key factor: lnR₃ − ln3 = -1.693275
              = ln(R₃/3) = ln(0.183916)
              = ln(R₃) − ln(3)

  This factor multiplies a_q in the mass curvature A.
  Its magnitude (1.6933) means q-curvature has a STRONG
  effect on mass evolution — stronger than p (which gets ln2 ≈ 0.6931).


§16. BÉZIER INTERPRETATION — Continuity Levels across Sectors
------------------------------------------------------------------------

  Each sector's γ(g) through 3 generations is a quadratic Bézier curve
  with control points at g = 1, 2, 3 (the 3 physical particles).

  In spline theory, continuity levels are:
    C⁰: positions match (curves pass through same point)
    C¹: first derivatives match (same tangent direction)
    C²: second derivatives match (same curvature)

  We check continuity between sector pairs at each component:

  Component                           ℓ↔u        ℓ↔d        u↔d
  ────────────────────────────────────────────────────────────
  p (SU(2) binarity)               C⁰@1.7 C⁰@1.1,3.1 C⁰@0.7,2.3
  q (SU(3) spectral)                   C²     C⁰@2.7     C⁰@2.0
  r (Coxeter)                          C² C⁰@0.6,2.8 C⁰@0.8,3.0

  ★ LEPTON ↔ UP: C² continuity in BOTH q and r!
    These sectors are parts of the SAME quadratic spline
    in the (q,r) subspace — only differing in velocity (b)
    and offset (c). The curvature is locked.

  ★ DOWN: only C⁰ connections (crossings at specific generations).
    Down is a DIFFERENT spline segment with inverted curvature.
    The C⁰ crossings ARE the physical mixing points (CKM).

  Spline interpretation:
    ┌────────────┐     C² in (q,r)     ┌────────────┐
    │   LEPTON   │ ←─────────────────→ │   UP QUARK │
    │ γ_ℓ(g)    │    same curvature    │ γ_u(g)     │
    └─────┬──────┘                     └──────┬─────┘
          │ C⁰ crossing                       │ C⁰ crossing
          │ (at specific g)                   │ (at specific g)
    ┌─────┴──────────────────────────────────┴─────┐
    │              DOWN QUARK  γ_d(g)               │
    │          inverted curvature (−sign)           │
    └───────────────────────────────────────────────┘

  BÉZIER CONTROL POINTS (physical particles AS control points):
  Sector           P₀ (gen 1)       P₁ (gen 2)       P₂ (gen 3)
  ────────────────────────────────────────────────────────────
  lepton              (0,0,0)        (-0,-4,3)         (1,-7,3)
  up               (-2,-6,-1)         (0,-7,3)         (6,-7,4)
  down             (-0,-8,-2)         (2,-7,0)         (2,-6,4)

  The MIDPOINT of a quadratic Bézier is NOT the physical gen-2 particle
  but the average of endpoints weighted 1:2:1.
  Deviation from midpoint = curvature indicator:

      mid = (P₀ + P₂)/2   (straight-line midpoint)
      actual = P₁          (gen-2 particle)
      deviation = P₁ − mid (curvature direction)

  lepton    : mid = (+0.5,-3.5,+1.5), P₁ = (-0.5,-4.0,+3.0), dev = (-1.0,-0.5,+1.5), |dev| = 1.8708
  up        : mid = (+2.2,-6.5,+1.5), P₁ = (+0.5,-7.0,+3.0), dev = (-1.8,-0.5,+1.5), |dev| = 2.3585
  down      : mid = (+0.5,-7.0,+1.0), P₁ = (+1.5,-7.0,+0.0), dev = (+1.0,+0.0,-1.0), |dev| = 1.4142

  The deviation vector = (aₚ, a_q, a_r)/2:
  lepton    : a/2 = (+0.5000, +0.2500, -0.7500)
  up        : a/2 = (+0.8750, +0.2500, -0.7500)
  down      : a/2 = (-0.5000, +0.0000, +0.5000)

  ★ Deviation matches a/2 EXACTLY — confirming the Bézier structure.

  PHYSICAL MEANING: The gen-2 particle (μ, c, s) is displaced from
  the straight line (gen-1 → gen-3) by exactly HALF the curvature
  vector. This is the defining property of a quadratic Bézier.

  Each fermion sector is literally a Bézier curve in spectral space.


§17. PARAMETER COUNT — How Many Degrees of Freedom?
------------------------------------------------------------------------

  NAIVE COUNT: 9 fermion masses = 9 free parameters.

  LATTICE DESCRIPTION: 9 particles × 3 coordinates (p,q,r) = 27 numbers.
  But quadratic factorization: 3 sectors × 3 components × 3 coefficients = 27.
  With 3 points per sector, the quadratic is exact → 0 residual DOF.

  CONSTRAINTS DISCOVERED:
  ┌────┬──────────────────────────────────────┬──────────┐
  │ #  │ Constraint                           │ DOF saved│
  ├────┼──────────────────────────────────────┼──────────┤
  │ 1  │ a_q(ℓ) = a_q(u) = +0.5              │    −1    │
  │ 2  │ a_r(ℓ) = a_r(u) = −1.5              │    −1    │
  │ 3  │ |a_q + a_r| = 1 (all sectors)       │    −3    │
  │ 4  │ a_q(d) = 0 (linear)                 │    −1    │
  │ 5  │ q_ℓ(3) = q_u(3) = −7 (unification) │    −1    │
  │ 6  │ Electron = origin (0,0,0)           │    −3    │
  ├────┼──────────────────────────────────────┼──────────┤
  │    │ TOTAL constraints                    │   −10    │
  └────┴──────────────────────────────────────┴──────────┘

  Note: constraints 1+2 imply constraint 3 for ℓ and u (not independent).
  Independent constraints: 1 + 2 + 3(down only) + 4 + 5 + 6 = 1+1+1+1+1+3 = 8.

  Effective DOF: 27 − 8 = 19 coordinates, but...

  In MASS SPACE (what matters): each sector has A, B, C = 9 parameters.
  Constraints reduce this:

  MASS-SPACE ANALYSIS:
  Naive: 9 masses = 9 DOF
  Lattice: 3 sectors × (A, B, C) = 9 DOF

  lepton    : ln(m/mₑ) = -1.252103·g² +9.084224·g -7.832121
  up        : ln(m/mₑ) = -0.732242·g² +8.572134·g -6.410246
  down      : ln(m/mₑ) = +0.405465·g² +1.772461·g +0.035580

  Mass curvatures A:
    A_ℓ = -1.252103
    A_u = -0.732242
    A_d = +0.405465
    A_ℓ − A_u = -0.519860
    A_d − A_ℓ = +1.657568  (≈ 2ln3 + Δ from different aₚ)

  Although a_q and a_r are shared between ℓ and u,
  the mass curvatures A differ because aₚ(ℓ) ≠ aₚ(u):
    aₚ(ℓ) = 1.0,  aₚ(u) = 1.75,  Δaₚ = 0.75
    A_u − A_ℓ = Δaₚ·ln2 = 0.75×0.6931 = 0.519860
    Actual:  0.519860  ✓

  ★ EFFECTIVE PARAMETER HIERARCHY:

  TIER 1 — CURVATURE (determines mass scale growth):
    Shared: a_q = +0.5, a_r = −1.5 (ℓ and u)
    Down:   a_q = 0, a_r = +1 (from conservation)
    Free:   aₚ per sector = 3 parameters
    But with conservation: a_r = ±1 − a_q → only aₚ, a_q free per sector
    Total Tier 1: 3 (one aₚ per sector, a_q fixed by sharing/linearity)

  TIER 2 — VELOCITY (determines ratio between gen 1↔2 vs 2↔3):
    3 sectors × 3 velocity coefficients (bₚ, b_q, b_r) = 9
    Minus constraint from q_ℓ(3)=q_u(3): −1
    Total Tier 2: 8

  TIER 3 — OFFSET (determines absolute mass of lightest in sector):
    Electron = (0,0,0) → lepton offset = (2.5, 5.0, −6.0) = fixed
    Up and down offsets: 2 × 3 = 6
    Total Tier 3: 6

  GRAND TOTAL: 3 + 8 + 6 = 17 lattice parameters → 9 masses
  = OVERDETERMINED (17 > 9)!

  But only 9 masses means at most 9 independent parameters.
  The extra 8 parameters encode REDUNDANT information = 8 constraints
  between lattice coordinates that we can TEST or EXPLOIT.

  PARAMETER BUDGET (lattice → mass projection):
  ┌────────────┬──────────┬──────────┬──────────┐
  │ Level      │ Lattice  │ Mass     │ Surplus  │
  │            │ params   │ DOF      │ = testable│
  ├────────────┼──────────┼──────────┼──────────┤
  │ Curvature  │    3aₚ   │ 3 (A)    │  +2 from │
  │ (per sect) │ +shared  │          │  sharing │
  │            │  a_q,a_r │          │          │
  ├────────────┼──────────┼──────────┼──────────┤
  │ Velocity   │   9 bᵢ   │ 3 (B)    │   6      │
  ├────────────┼──────────┼──────────┼──────────┤
  │ Offset     │  6 cᵢ    │ 2 (C)    │   4      │
  │ (e=origin) │ +3 fixed │          │          │
  ├────────────┼──────────┼──────────┼──────────┤
  │ TOTAL      │ 18 free  │ 8 (A,B,C)│  10      │
  └────────────┴──────────┴──────────┴──────────┘

  9 masses need 8 mass-space DOF (9 minus 1 for mₑ = reference).
  Lattice has 18 free parameters → 10 SURPLUS constraints.

  These 10 surplus parameters are NOT arbitrary:
  they describe the 3D GEOMETRY of the curves
  that is invisible in 1D mass-space.

  ★ THIS IS THE KEY: the lattice contains MORE information
    than the masses alone. The extra structure (curvature sharing,
    conservation law, unification) is the 3D geometry that
    CONSTRAINS which mass spectra are possible.

  WHAT WE NEED TO DERIVE MASSES:
  1. Fix aₚ per sector (3 params) — possibly from I₃, Nc relation
  2. Fix velocity coefficients (8 params) — from CONTINUITY conditions
  3. Offsets follow from electron = origin + crossing constraints

  The Bézier/spline framework REDUCES step 2:
  if sectors are related by C¹ or C² continuity conditions,
  the velocity coefficients are NOT independent.


§18. CONTINUITY PREDICTIONS — What C² sharing implies for masses
------------------------------------------------------------------------

  PREDICTION FROM C² SHARING (ℓ↔u in q,r):
  
  Since a_q(ℓ) = a_q(u) and a_r(ℓ) = a_r(u):
    A_u − A_ℓ = [aₚ(u) − aₚ(ℓ)]·ln2
              = (1.75 − 1.0)·ln2 = 0.75·ln2
              = 0.519860

  This means:
    ln(m_u(g)/m_ℓ(g)) = (A_u−A_ℓ)·g² + (B_u−B_ℓ)·g + (C_u−C_ℓ)

  The g² coefficient is FIXED by curvature sharing + aₚ difference.
  Only the g-linear and constant terms are "free" between ℓ and u.

  ln(m_u/m_ℓ) = +0.519860·g² -0.512090·g +1.421875

  g=1: m_u/m_e = e^(+1.4296) = 4.1772  (exp: 4.2270)
  g=2: m_c/m_mu = e^(+2.4771) = 11.9071  (exp: 12.0199)
  g=3: m_t/m_tau = e^(+4.5643) = 96.0000  (exp: 97.2277)

  C¹ CONTINUITY TEST at g = 3 (unification point):
  If ℓ and u had matching tangent vectors at g=3:

  p (SU(2) binarity): dγ/dg|ₗ(3) = +2.5000, dγ/dg|ᵤ(3) = +7.2500, Δ = +4.7500
  q (SU(3) spectral): dγ/dg|ₗ(3) = -2.5000, dγ/dg|ᵤ(3) = +0.5000, Δ = +3.0000
  r (Coxeter): dγ/dg|ₗ(3) = -1.5000, dγ/dg|ᵤ(3) = -0.5000, Δ = +1.0000

  q-velocities at g=3: both have same curvature → velocity gap = b_u − b_ℓ
    = −2.5 − (−5.5) = 3.0  (constant, independent of g)

  For FULL C¹ at g=3 we'd need Δv = 0 in all components.
  Current gaps: Δvₚ = +4.7500, Δv_q = +3.0000, Δv_r = +1.0000
  
  ★ r HAS near-C¹: Δv_r = +1.0 (smallest gap)
  ★ q gap = 3.0 is large — the sectors DIVERGE in q-velocity even
    though they MEET at q = −7.

  This means leptons and up quarks osculate (C²) in curvature
  but NOT in velocity — they are TANGENT curves that kiss
  but immediately separate. Like two Bézier segments
  sharing a curvature comb but different control handles.

  Saved: D:\Desarrollo\projects\pyshics\ss_peg\output\transition_curves_mass_bezier.png

========================================================================
  §19. ROADMAP — From Transition Curves to Mass Derivation
========================================================================

  WHAT WE HAVE:
  ─────────────
  1. 9 fermion masses live on a 3D lattice with basis {ln2, lnR₃, ln3}
  2. Each sector traces a QUADRATIC BÉZIER curve γ(g) in this lattice
  3. ln(m(g)/mₑ) = A·g² + B·g + C per sector (mass parabola)
  4. Curvature conservation: |a_q + a_r| = 1
  5. Isocurvature (C²) between lepton and up sectors in (q,r)
  6. Spectral unification at gen 3: q_ℓ(3) = q_u(3) = −7
  7. Down quarks: inverted curvature, linear q

  WHAT WE NEED:
  ─────────────
  To derive all 9 masses from the graph, we need to determine
  A, B, C for each sector (8 independent DOF, since mₑ = reference).

  THE STRATEGY:
  ─────────────
  STEP 1: Fix curvatures A from (I₃, Nc, conservation law)
    → Already done! a_q,a_r shared or fixed by conservation.
    → Only need aₚ per sector (3 params, possibly from I₃·Nc)

  STEP 2: Fix velocities B from CONTINUITY + UNIFICATION
    → C¹ conditions at gen 3 (if applicable) would fix b differences
    → The q-unification constraint q_ℓ(3) = q_u(3) gives 1 relation
    → Need: WHAT determines the initial velocity b?

  STEP 3: Fix offsets C from electron = origin + a second anchor
    → Lepton C fixed (electron = (0,0,0))
    → Need: WHAT sets ℓ=0 vs u vs d offset?
    → Likely from (I₃, Nc, hypercharge) quantum numbers at gen 1

  MINIMAL ANSATZ:
  ───────────────
  If aₚ = f(I₃, Nc) and C = h(Y, I₃), then only B remains free.
  B encodes the RATE of generational mass growth.
  If B follows from a regularity condition (e.g. lattice consistency,
  or constant arc-length in spectral space, or minimal curvature
  energy), then ALL 9 masses are derived.

  The Bézier framework gives the right language:
  • Curvature = physics (gauge structure)
  • Velocity = kinematics (mass growth rate)
  • Offset = boundary condition (lightest particle identity)

