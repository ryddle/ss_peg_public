  ===========================================================================
  PART 1: THE DATA — N AND ITS 2-ADIC STRUCTURE
  ===========================================================================

  The transition matrix N:
     -1   -8   +6   ← lepton (ℓ)
     +4   -2   +8   ← up (u)
     +4   +2   +4   ← down (d)
       p     q     r

  Entry-by-entry decomposition:
      Sector       Coord     n   |n|   2^v   v₂   odd  flat?  pow2?
    ───────────────────────────────────────────────────────────────────────────
    lepton (ℓ)    p (ln 2)    -1     1  2^0     0    -1   free      ✓
    lepton (ℓ)   q (ln R₃)    -8     8  2^3     3    -1   free      ✓
    lepton (ℓ)    r (ln 3)    +6     6  2^1     1    +3   FLAT  ✗ (3^1)
      up (u)    p (ln 2)    +4     4  2^2     2    +1   free      ✓
      up (u)   q (ln R₃)    -2     2  2^1     1    -1   FLAT      ✓
      up (u)    r (ln 3)    +8     8  2^3     3    +1   free      ✓
    down (d)    p (ln 2)    +4     4  2^2     2    +1   FLAT      ✓
    down (d)   q (ln R₃)    +2     2  2^1     1    +1   free      ✓
    down (d)    r (ln 3)    +4     4  2^2     2    +1   free      ✓

  Free entries: 6, all |n| = 2^v: YES ✓
  Flat entries: 3, all |n| = 2^v: NO ✗
    Exception: n_l,r = 6 = (1, 1, 1, 1)

  ===========================================================================
  PART 2: {2,3}-SMOOTHNESS OF ALL ENTRIES
  ===========================================================================

  All 9 entries of N:
     n_l,p =  -1 =         −2^0   ✓ {2,3}-smooth
     n_l,q =  -8 =         −2^3   ✓ {2,3}-smooth
     n_l,r =  +6 =     +2^1·3^1   ✓ {2,3}-smooth
     n_u,p =  +4 =         +2^2   ✓ {2,3}-smooth
     n_u,q =  -2 =         −2^1   ✓ {2,3}-smooth
     n_u,r =  +8 =         +2^3   ✓ {2,3}-smooth
     n_d,p =  +4 =         +2^2   ✓ {2,3}-smooth
     n_d,q =  +2 =         +2^1   ✓ {2,3}-smooth
     n_d,r =  +4 =         +2^2   ✓ {2,3}-smooth

  ALL entries {2,3}-smooth: YES ✓

  This means: n = ±2^a · 3^b for all entries, with a ≥ 0, b ≥ 0.
  Since the lattice basis is (ln 2, ln R₃, ln 3), only the primes
  2 and 3 appear — the number system is exactly {2,3}-smooth.

  Therefore: |n| = 2^v  ⟺  3 ∤ n  (for {2,3}-smooth integers).
  
  To prove S4, we only need to show: 3 DOES NOT DIVIDE any free entry.

  ===========================================================================
  PART 3: THE MOD-3 PROOF — S4 IS A THEOREM OF S1
  ===========================================================================

  The generating polynomials n_k(I) = α_k + β_k·I + γ_k·I²:
    n_p(I) = -1 +0·I +5·I²
    n_q(I) = -8 -2·I +8·I²
    n_r(I) = +6 +2·I +0·I²

  Evaluation table (all 9 entries):
    Column    I    Sector     n   n mod 3   3|n?   Type
    ───────────────────────────────────────────────────────
         p   +0    lepton    -1         2  no  ✓   free
         p   +1        up    +4         1  no  ✓   free
         p   -1      down    +4         1  no  ✓   FLAT
         q   +0    lepton    -8         1  no  ✓   free
         q   +1        up    -2         1  no  ✓   FLAT
         q   -1      down    +2         2  no  ✓   free
         r   +0    lepton    +6         0  YES ✗   FLAT
         r   +1        up    +8         2  no  ✓   free
         r   -1      down    +4         1  no  ✓   free

  ╔══════════════════════════════════════════════════════════════════════╗
  ║  THEOREM (S4 — mod 3 proof):                                       ║
  ║                                                                     ║
  ║  Given the generating polynomials n_k(I) from S1 (§12.3):          ║
  ║    n_p(I) = −1 + 5I²                                               ║
  ║    n_q(I) = −8 − 2I + 8I²                                          ║
  ║    n_r(I) = 6 + 2I                                                  ║
  ║                                                                     ║
  ║  Reduce mod 3:                                                      ║
  ║    n_p(I) ≡ 2 + 2I²  (mod 3)                                       ║
  ║    n_q(I) ≡ 1 + I + 2I²  (mod 3)                                   ║
  ║    n_r(I) ≡ 2I  (mod 3)                                             ║
  ║                                                                     ║
  ║  Evaluating at I ∈ {0, ±1}:                                         ║
  ║    n_p mod 3: {2, 1, 1}     — never 0                               ║
  ║    n_q mod 3: {1, 1, 2}     — never 0                               ║
  ║    n_r mod 3: {0, 2, 1}     — 0 only at I=0 (lepton = FLAT entry)  ║
  ║                                                                     ║
  ║  ⟹ No FREE entry is divisible by 3.                                ║
  ║  ⟹ Since all entries are {2,3}-smooth: |n_free| = 2^v.  ∎          ║
  ╚══════════════════════════════════════════════════════════════════════╝

  ===========================================================================
  PART 4: WHY THE FACTOR OF 3 FALLS ON THE FLAT ENTRY
  ===========================================================================

  The r-column polynomial is n_r(I) = 6 + 2I.

  At I = 0 (leptons): n_ℓr = 6 = 2·3.
  This is the FLAT entry — forced by the flatness constraint b = −5a.

  The r-coordinate corresponds to ln 3 (the SU(3) Coxeter direction).
  The curvature in this direction for leptons is a_ℓr = −3/2, so:
    n_flat = −4a = −4·(−3/2) = 6 = 2·3.

  The factor of 3 enters because the lattice basis vector in the
  r-direction IS ln 3 — the Coxeter number h₃ = 3 is baked into
  the lattice geometry.

  But this entry is FLAT (anti-diagonal), so it doesn't count as
  a free parameter. The flatness permutation principle (S2) places
  this 3-contaminated entry exactly where it doesn't affect S4.

  Anti-diagonal (flat) assignments:
    (ℓ, r): n = +6 = 2^1·3^1  ← contains factor 3
    (u, q): n = -2 = −2^1  ← pure power of 2
    (d, p): n = +4 = 2^2  ← pure power of 2

  The flatness principle (S2) assigns the r-column flat entry to the
  LEPTON sector. Since leptons have I₃ = 0, the polynomial n_r(I=0) = 6
  picks up the factor of 3 from the lattice. But S2 ensures this
  lands on the anti-diagonal, making it flat — hence harmless for S4.

  ===========================================================================
  PART 5: THREE MECHANISMS — WHY EACH COLUMN GIVES POWERS OF 2
  ===========================================================================

  ─── Column p: n_p(I) = −1 + 5I² ───

  Mechanism: CONSTANT + QUADRATIC CANCELLATION
  • α_p = −1 = −2⁰ (the lepton value)
  • γ_p = 5 (quark-lepton splitting)
  • At I² = 1: α + γ = −1 + 5 = 4 = 2²

  Why γ_p = 5?
    γ_p = [n_p(+1) + n_p(−1)]/2 − n_p(0) = (4+4)/2 − (−1) = 5
    = h₂ + h₃ = 2 + 3 (Coxeter numbers of SU(2) and SU(3))

  So the key relation is: α_p + γ_p = −1 + (h₂ + h₃) = −1 + 5 = 4 = 2²
  This works because h₂ + h₃ = 5 ≡ 2 mod 3, and −1 ≡ 2 mod 3,
  so −1 + 5 = 4 ≡ 1 mod 3 (not divisible by 3).

  ─── Column q: n_q(I) = −8 − 2I + 8I² ───

  Mechanism: QUADRATIC CANCELS CONSTANT (γ = −α)
  • α_q = −8 = −2³
  • γ_q = +8 = 2³ = −α_q
  • At I = ±1: n = α + γ ± β = 0 ± (−2) = ∓2

  The quadratic exactly annihilates the constant:
    n_q(±1) = ±β_q = ∓2 = ∓2¹

  This happens because the UP and DOWN quark q-displacements average to zero:
    n_q(+1) + n_q(−1) = −2 + 2 = 0
  ⟹ The q-direction has ISOSPIN antisymmetry for quarks.

  Why γ_q = −α_q?
    n_q values: ['-8', '-2', '2']
    Average of quark n_q: (-2 + 2)/2 = 0
    This equals α_q + γ_q = -8 + 8 = 0

  ─── Column r: n_r(I) = 6 + 2I ───

  Mechanism: LINEAR SHIFT FROM {2,3}-SMOOTH BASE
  • α_r = 6 = 2·3 (flat entry, contains factor 3)
  • β_r = 2 (the linear term)
  • At I = ±1: n = α ± β = 6 ± 2 = {8, 4}

  The pair (6, 2) has the special property:
    6 + 2 = 8 = 2³
    6 − 2 = 4 = 2²
  Both are powers of 2 despite 6 itself containing factor 3.

  This works because 6 ≡ 0 mod 3 and 2 ≡ 2 mod 3, so:
    6 ± 2 ≡ ±2 ≡ {2, 1} mod 3 — neither is 0.

  Geometrically: the linear term β = 2 is the isospin-driven
  displacement that REMOVES the factor of 3 from the quark entries
  by shifting away from the 3-divisible base value.


  ===========================================================================
  PART 6: CONNECTION TO COXETER NUMBERS h₂ = 2, h₃ = 3
  ===========================================================================

  SU(2): h₂ = 2,  SU(3): h₃ = 3
  h₂ + h₃ = 5
  h₂ · h₃ = 6

  Polynomial coefficient origins:
    γ_p = 5 = h₂ + h₃ = 2 + 3
    γ_q = 8 = 2³  (pure power of 2)
    α_r = 6 = 2·h₃ = 2·3
    β_r = 2 = h₂ = 2
    β_q = -2 = −h₂ = −2

  The three column mechanisms in terms of Coxeter numbers:
    p: α + γ = −1 + (h₂+h₃) = −1 + 5 = 4 = 2²
    q: γ = −α = 2³, so quark values = ±β = ±h₂ = ±2
    r: α ± β = 2h₃ ± h₂ = 2·3 ± 2 = {8, 4} = {8, 4}

  Key number-theoretic relations:
    −1 + (h₂ + h₃) = 4 = 2²  ✓ (power of 2)
    2h₃ + h₂ = 8 = 2³  ✓ (power of 2)
    2h₃ − h₂ = 4 = 2²  ✓ (power of 2)
    h₂ = 2 = 2¹  ✓ (power of 2)
    h₃ = 3 = 3¹  (NOT a power of 2 — but appears only in flat entry)

  These relations hold because h₂ = 2 is itself a power of 2,
  and h₃ = 3 is the 'other prime'. The combinations 2h₃ ± h₂ = 6 ± 2
  land on powers of 2 because 3·2 ± 2 = 2(3±1) = 2·{4,2} = {8,4}.

  ╔══════════════════════════════════════════════════════════════════════╗
  ║  The power-of-2 pattern emerges because:                            ║
  ║  (1) h₂ = 2 is a power of 2, providing 'clean' linear shifts       ║
  ║  (2) h₃ = 3 appears only as 2·h₃ = 6 in the r-direction          ║
  ║  (3) 2h₃ ± h₂ = 2(h₃ ± 1) = 2·{4,2} — the ±1 shift kills the 3  ║
  ║  (4) The flatness principle places the lone 3-divisible entry       ║
  ║      on the anti-diagonal, where it's constrained (not free)        ║
  ╚══════════════════════════════════════════════════════════════════════╝

  ===========================================================================
  PART 7: THE ODD-PART MATRIX — S4 IN ONE LINE
  ===========================================================================

  The odd part u(n) = n / 2^{v₂(n)} strips all factors of 2:

  u(N) =
      -1   -1   +3*   ← lepton (ℓ)
      +1   -1*  +1    ← up (u)
      +1*  +1   +1    ← down (d)
       p    q    r
  (* = flat entry)

  S4 states: ALL free entries have u(n) = ±1.
  Equivalently: the odd-part matrix restricted to free entries is ±1.

  Free odd parts: [np.int64(-1), np.int64(-1), np.int64(1), np.int64(1), np.int64(1), np.int64(1)]
  All ±1: YES ✓
  Flat odd parts: [np.int64(3), np.int64(-1), np.int64(1)]

  ===========================================================================
  PART 8: THE v₂ MATRIX AND EXPONENT DISTRIBUTION
  ===========================================================================

  v₂(N) =
      0   3   1*   ← lepton (ℓ)
      2   1*  3    ← up (u)
      2*  1   2    ← down (d)
     p  q  r

  Row sums:   [np.int64(4), np.int64(6), np.int64(5)] → sectors ℓ, u, d
  Col sums:   [np.int64(4), np.int64(5), np.int64(6)] → coords p, q, r
  Trace:      3 (= v₂(n_ℓp) + v₂(n_uq) + v₂(n_dr))
  Anti-diag:  4 (flat entries)
  Total:      15

  Free v₂ values (sorted): [np.int64(0), np.int64(1), np.int64(2), np.int64(2), np.int64(3), np.int64(3)]
  Multiset: {0, 1, 2, 2, 3, 3}, sum = 11
  Mean: 1.83, Median: 2

  ===========================================================================
  PART 9: S4 + S1 → CONSTRAINTS ON THE SMITH FORM
  ===========================================================================

  If S4 holds (all free |n| = 2^v), what can we say about Smith(N)?

  The 9 entries of N are:
    Free: {±1, ±2, ±4, ±8} (powers of 2)
    Flat: {-4a} where a comes from the curvature

  Smith form d₁ = gcd of all entries.
  Since n_ℓp = −1 (the unique v₂ = 0 entry), gcd = 1.
  ⟹ d₁ = 1  ✓ (first Smith invariant is forced)

  All 2×2 minors: [34, -32, -52, 30, -28, -44, 16, -16, -24]
  gcd of all 2×2 minors: 2
  ⟹ d₁·d₂ = 2, so d₂ = 2
  det(N) = -8
  d₁·d₂·d₃ = |det| = 8
  ⟹ d₃ = 4
  Smith(N) = (1, 2, 4) = (1, 2, 4)  ✓

  Can we derive d₂ = 2 from S4 (power-of-2 structure)?

  A 2×2 minor has the form: n₁n₄ − n₂n₃.
  If all nᵢ = ±2^{vᵢ}, then the minor = ±2^{v₁+v₄} ± 2^{v₂+v₃}.

  For the gcd to be exactly 2, we need at least one minor
  with v₂(minor) = 1 and no minor with v₂ = 0.

  v₂ of each 2×2 minor:
    rows(0,1) cols(0,1): minor =   +34, v₂ = 1
    rows(0,1) cols(0,2): minor =   -32, v₂ = 5
    rows(0,1) cols(1,2): minor =   -52, v₂ = 2
    rows(0,2) cols(0,1): minor =   +30, v₂ = 1
    rows(0,2) cols(0,2): minor =   -28, v₂ = 2
    rows(0,2) cols(1,2): minor =   -44, v₂ = 2
    rows(1,2) cols(0,1): minor =   +16, v₂ = 4
    rows(1,2) cols(0,2): minor =   -16, v₂ = 4
    rows(1,2) cols(1,2): minor =   -24, v₂ = 3

  Minimum v₂ among all minors: 1
  This gives d₂ = 2^1 = 2

  The n_ℓp = −1 entry (v₂ = 0) participates in 4 minors:
  Any minor involving row 0, col 0 has terms like:
    (−1)·n_{ij} − n_{0j}·n_{i0}
  where at least one factor has v₂ = 0.

  For such a minor: v₂(minor) = v₂(−n_{ij} − n_{0j}·n_{i0})
  Since n_{ij} and n_{0j}·n_{i0} are both integers with v₂ ≥ 0,
  the minor's v₂ depends on whether the two terms have the same
  or different v₂. This requires case analysis of the specific entries.

  ===========================================================================
  PART 10: HOW SPECIAL IS S4? RANDOM POLYNOMIAL ANALYSIS
  ===========================================================================

  Q: For a quadratic polynomial n(I) = α + βI + γI² with integer
     coefficients, how often do ALL three values at I = {0,+1,−1}
     come out as ±2^v (i.e., ALL have trivial odd part)?

  Range: α, β, γ ∈ [−10, +10]
  Total polynomials (no zero values): 8220
  All 3 values are ±2^v: 400
  Fraction: 400/8220 = 4.8662%

  At least 2 of 3 values are ±2^v: 2208
  Fraction: 2208/8220 = 26.8613%

  For our specific generating polynomials:
    Column p: values = [-1, 4, 4], free = [-1, 4], all free pow2: True
    Column q: values = [-8, -2, 2], free = [-8, 2], all free pow2: True
    Column r: values = [6, 8, 4], free = [8, 4], all free pow2: True

  The probability that ALL THREE COLUMNS simultaneously satisfy S4
  is roughly (4.9%)³ ≈ 1.2e-04 if the columns were independent —
  showing S4 is highly non-trivial from a random perspective.
  
  But from the generating function perspective, the polynomial
  coefficients are DETERMINED by the graph topology, not random.
  S4 is a necessary consequence of the specific Coxeter structure
  of SU(2)×SU(3).

  ===========================================================================
  PART 11: WHY (2, 3)-SMOOTH? THE LATTICE BASIS
  ===========================================================================

  The mass formula is: ln(m_f / m_e) = p·ln 2 + q·ln R₃ + r·ln 3.

  The lattice Λ has basis {ln 2, ln R₃, ln 3} and all coordinates
  x(g) ∈ (Z/2)³ — half-integer multiples of the basis vectors.

  The displacement n = 2(x(2)−x(1)) has INTEGER coordinates in this basis.
  These integers involve only the LATTICE STRUCTURE, not the irrational
  values R₃, ln 2, ln 3 themselves.

  From the parabolic formula: n = 6a + 2b, where a, b are exact rationals.
  The parabolic parameters involve ratios of graph invariants (§10.5).
  Since the graph invariants of SU(2)×SU(3) involve only:
    • Integers 1, 2, 3, 4, 5, 6, 7, 8, 9, ...
    • Fractions with denominators dividing 4 or 8

  the resulting n values are always {2,3}-smooth integers.

  Explicit computation: n = 6a + 2b for each entry:
    Sector  Coord         a         b        6a        2b     n  2,3-sm
    ─────────────────────────────────────────────────────────────────
         ℓ      p         1      -7/2         6        -7    -1       ✓  
         ℓ      q       1/2     -11/2         3       -11    -8       ✓  
         ℓ      r      -3/2      15/2        -9        15    +6       ✓ *
         u      p       7/4     -13/4      21/2     -13/2    +4       ✓  
         u      q       1/2      -5/2         3        -5    -2       ✓ *
         u      r      -3/2      17/2        -9        17    +8       ✓  
         d      p        -1         5        -6        10    +4       ✓ *
         d      q         0         1         0         2    +2       ✓  
         d      r         1        -1         6        -2    +4       ✓  

  The denominators in a, b are at most 4 (from §10.1).
  So 6a has denominator ≤ 4/gcd(6,4) = 2, and 2b has denominator ≤ 2.
  Their sum n = 6a + 2b is always an integer (verified above).

  The numerators involve only small prime factors because the
  graph invariants of SU(2)×SU(3) are all {2,3,5,7}-smooth.
  But the specific combinations 6a + 2b turn out to be {2,3}-smooth
  for ALL 9 entries — eliminating factors of 5 and 7.

  ===========================================================================
  PART 12: COUNTERFACTUAL — WHAT IF h₃ ≠ 3?
  ===========================================================================

  Would S4 survive if the gauge group had different Coxeter numbers?
  Test: replace h₃ by hypothetical values and check the power-of-2 pattern.

  The key relations involving h₃:
    • α_r = 2·h₃ (flat lepton-r entry)
    • β_r = h₂ = 2
    • n_r(±1) = 2h₃ ± 2
    • γ_p ≈ h₂ + h₃ → α_p + γ_p = −1 + h₂ + h₃

    h₃    2h₃+2    2h₃−2    −1+2+h₃   All pow2?
  ──────────────────────────────────────────
     1        4  0 (ZERO)          2        ✗ no
     2        6        2          3        ✗ no
     3        8        4          4       ✓ YES ★ (SU(3))
     4       10        6          5        ✗ no
     5       12        8          6        ✗ no
     6       14       10          7        ✗ no
     7       16       12          8        ✗ no
     8       18       14          9        ✗ no
     9       20       16         10        ✗ no
    10       22       18         11        ✗ no
    11       24       20         12        ✗ no
    12       26       22         13        ✗ no

  Note: the p-column condition −1 + h₂ + h₃ = 2^v alone admits Mersenne-like
  values h₃ ∈ {1, 3, 7, 15, ...} (2^k − 1). But the r-column requires BOTH
  2h₃ + 2 and 2h₃ − 2 to be powers of 2 (or zero). This is far more restrictive:
  • h₃ = 1: 2h₃ − 2 = 0  (degenerate)
  • h₃ = 7: 2h₃ − 2 = 12  (NOT a power of 2)
  • h₃ = 15: 2h₃ − 2 = 28  (NOT a power of 2)

  All three conditions simultaneously hold ONLY for h₃ = 3 (SU(3)).
  This is a sharp selection criterion on the gauge group itself!

  ===========================================================================
  PART 13: DOES S1 DETERMINE S5 (SMITH FORM)?
  ===========================================================================

  S1 (the generating function) determines the polynomial coefficients,
  which determine all 9 entries of N, which determine Smith(N).

  This means S5 is ALSO a consequence of S1:
  S1 → N → Smith(N) = (1, 2, 4).

  The only question is whether S5 can be understood WITHOUT
  computing the full matrix — i.e., from structural properties alone.

  From the v₂ matrix:
    v₂(N) = [[0, 3, 1], [2, 1, 3], [2, 1, 2]]
    det(N) = -8 = −2³
    d₁ · d₂ · d₃ = |det| = 8 = 2³
    d₁ | d₂ | d₃ (divisibility chain)

  With d₁ = 1 (forced by n_ℓp = −1 being odd):
    d₂ · d₃ = 8, with d₂ | d₃.
    Options: (d₂, d₃) ∈ {(1,8), (2,4)}

    From gcd of 2×2 minors = 2: d₂ = 2, d₃ = 4.

  Can gcd(minors) = 2 be derived from S4 alone?

  2×2 minors involve products of pairs. With the v₂ structure:
    min(v₂ of row_i/col_j product − row_i/col_k product)

  The minor with minimum v₂:
    rows (0,1), cols (0,1):
    = n_00·n_11 − n_01·n_10
    = (-1)·(-2) − (-8)·(4)
    = 2 − -32
    = 34
    v₂ = 1

  ===========================================================================
  PART 14: S4 VIA THE DISPLACEMENT IDENTITY
  ===========================================================================

  From §12.7: n = 2(x(2) − x(1)), so |n| = 2|Δ₁₂|.
  S4 says |n_free| = 2^v, equivalently |Δ₁₂| = 2^{v−1} ∈ Z/2.

  Displacement matrix Δ₁₂ = N/2:
        -0.5     -4.0     +3.0*
        +2.0     -1.0*    +4.0 
        +2.0*    +1.0     +2.0 

  The displacements are:
    Δ₁₂(ℓ,p) =   -0.5 → |Δ| =   0.5  ✓ 2^-1
    Δ₁₂(ℓ,q) =   -4.0 → |Δ| =   4.0  ✓ 2^2
    Δ₁₂(ℓ,r) =   +3.0 → |Δ| =   3.0  ✗ (FLAT)
    Δ₁₂(u,p) =   +2.0 → |Δ| =   2.0  ✓ 2^1
    Δ₁₂(u,q) =   -1.0 → |Δ| =   1.0  ✓ 2^0 (FLAT)
    Δ₁₂(u,r) =   +4.0 → |Δ| =   4.0  ✓ 2^2
    Δ₁₂(d,p) =   +2.0 → |Δ| =   2.0  ✓ 2^1 (FLAT)
    Δ₁₂(d,q) =   +1.0 → |Δ| =   1.0  ✓ 2^0
    Δ₁₂(d,r) =   +2.0 → |Δ| =   2.0  ✓ 2^1

  ALL displacements are of the form ±2^j (j ∈ {−1, 0, 1, 2}) EXCEPT
  Δ₁₂(ℓ, r) = 3 — which is the FLAT entry.

  The displacement in each non-flat direction is a BINARY step:
  moving by ½, 1, 2, or 4 half-integer lattice units.

  ===========================================================================
  PART 15: WHY n_ℓp = −1 IS THE UNIQUE ODD ENTRY
  ===========================================================================

  n_ℓp = −1 = −2⁰ is the ONLY entry with v₂ = 0. Why?

  1. ELECTRON-AT-ORIGIN: x_ℓ(1) = (0, 0, 0), so n_ℓk = 2·x_μ,k.
     For n_ℓk to have v₂ = 0, we need x_μ,k ∈ Z/2 \ Z 
     (a half-integer that's NOT an integer).

  2. MUON COORDINATES: x_μ = (−1/2, −4, 3).
     Only the p-coordinate is a non-integer half-integer (−1/2).
     The q and r coordinates are integers (−4 and 3).

  3. Therefore n_ℓp = 2·(−1/2) = −1 is the unique odd entry.

  4. Why is x_μ,p = −1/2?
     The p-coordinate is ln 2, associated with SU(2).
     From n = 6a + 2b: n_ℓp = 6(1) + 2(−7/2) = 6 − 7 = −1.
     Since n_ℓp = 2·x_μ,p, we get x_μ,p = −1/2.
     The parabolic parameter a_ℓp = 1 is the integer curvature
     for this cell; combined with b_ℓp = −7/2, they yield the
     unique odd entry in the entire matrix.

  5. MINIMALITY: |n_ℓp| = 1 = 2⁰ is the smallest possible |n| > 0
     for an integer, making this the minimum-displacement entry.

  Curvatures |a| for all 9 entries:
    a_ℓ,p =      1 =  +1.0000  |a| = 1.0000
    a_ℓ,q =    1/2 =  +0.5000  |a| = 0.5000
    a_ℓ,r =   -3/2 =  -1.5000  |a| = 1.5000
    a_u,p =    7/4 =  +1.7500  |a| = 1.7500
    a_u,q =    1/2 =  +0.5000  |a| = 0.5000
    a_u,r =   -3/2 =  -1.5000  |a| = 1.5000
    a_d,p =     -1 =  -1.0000  |a| = 1.0000
    a_d,q =      0 =  +0.0000  |a| = 0.0000
    a_d,r =      1 =  +1.0000  |a| = 1.0000

  Minimum nonzero |a|: 0.5 (at ℓ,q, u,q)

  ===========================================================================
  PART 16: THE COMPLETE PROOF — S4 AS A THEOREM OF S1
  ===========================================================================


  ╔══════════════════════════════════════════════════════════════════════╗
  ║                    THEOREM: S4 IS A THEOREM OF S1                   ║
  ╠══════════════════════════════════════════════════════════════════════╣
  ║                                                                     ║
  ║  CLAIM: Every free entry of N satisfies |n| = 2^v (v ∈ {0,1,2,3}). ║
  ║                                                                     ║
  ║  PROOF:                                                             ║
  ║                                                                     ║
  ║  Step 1. (Generating function — S1, §12.3)                         ║
  ║  The columns of N are quadratic polynomials in I = 2I₃:            ║
  ║    n_p(I) = −1 + 5I²                                               ║
  ║    n_q(I) = −8 − 2I + 8I²                                          ║
  ║    n_r(I) = 6 + 2I                                                  ║
  ║  evaluated at I ∈ {0, +1, −1} (lepton, up, down).                  ║
  ║                                                                     ║
  ║  Step 2. ({2,3}-smoothness)                                        ║
  ║  Direct evaluation: all 9 values are                                ║
  ║    {−1, −8, 6, 4, −2, 8, 4, 2, 4}                                  ║
  ║  which are all {2,3}-smooth (only prime factors 2 and 3).           ║
  ║                                                                     ║
  ║  Step 3. (Mod-3 analysis)                                          ║
  ║  Reduce the polynomials mod 3:                                      ║
  ║    n_p ≡ 2 + 2I²  (mod 3)                                          ║
  ║    n_q ≡ 1 + I + 2I²  (mod 3)                                      ║
  ║    n_r ≡ 2I  (mod 3)                                                ║
  ║                                                                     ║
  ║  Evaluate at I ∈ {0, ±1}:                                           ║
  ║    p: {2, 1, 1} — no zeros                                          ║
  ║    q: {1, 1, 2} — no zeros                                          ║
  ║    r: {0, 2, 1} — zero ONLY at I=0 (lepton/r = FLAT entry)         ║
  ║                                                                     ║
  ║  Step 4. (Conclusion)                                               ║
  ║  No free entry is divisible by 3.                                   ║
  ║  Combined with {2,3}-smoothness: |n_free| = 2^v.                   ║
  ║  The exponents are v ∈ {0, 1, 2, 2, 3, 3}.  ∎                     ║
  ║                                                                     ║
  ║  COROLLARY: S5 (Smith form) is also determined by S1.              ║
  ║  Since S1 fixes all 9 entries, Smith(N) = (1,2,4) follows.         ║
  ║                                                                     ║
  ║  STATUS OF THE DERIVATION CHAIN:                                    ║
  ║  S1 (generating function) alone determines N completely.            ║
  ║  S2 (flatness), S3 (lattice), S4 (powers of 2), S5 (Smith form),  ║
  ║  and S6a (sign rule) are ALL consequences of S1.                    ║
  ║                                                                     ║
  ║  Only S6b (min ||N−I||²) remains as a potentially independent      ║
  ║  selection principle — choosing among the 6 sign-compatible         ║
  ║  candidates generated by the generating function.                   ║
  ╚══════════════════════════════════════════════════════════════════════╝

  ===========================================================================
  PART 17: REVISED DERIVATION ARCHITECTURE
  ===========================================================================

  BEFORE (§13 of 00_DERIVATION.md):
  
    262,144  →[S4+S5: Smith(1,2,4)]→  3,056  →[S6a: sign]→  6  →[S6b: max tr]→  1

  AFTER (incorporating this investigation + yesterday's S6a proof):

    The generating function S1 determines the polynomial coefficients,
    which fix all 9 entries of N. Rules S2–S6a are checkpoints, not inputs:

    ┌─────── S1: Generating Function F(2I₃, N_c, g) ───────┐
    │                                                        │
    │  Polynomial coefficients → evaluate at I ∈ {0,±1}     │
    │                                                        │
    │  ⟹ S2 (flatness):   anti-diagonal entries = −4a       │  THEOREM
    │  ⟹ S3 (lattice):    all x(g) ∈ Z/2                   │  THEOREM
    │  ⟹ S4 (powers of 2): free |n| = 2^v (mod 3 proof)    │  THEOREM
    │  ⟹ S5 (Smith form):  Smith(N) = (1,2,4) (computed)    │  THEOREM
    │  ⟹ S6a (sign rule): sign = (−1)^{1+Nc} (§13.2)      │  THEOREM
    │                                                        │
    │  These determine 6 CANDIDATE matrices.                │
    │  S6b selects the physical one:  min ||N−I||²_F        │  OBSERVED
    └────────────────────────────────────────────────────────┘

  The derivation has collapsed from 7 independent rules to:
    1 theorem (S1) + 1 axiom (S6b)  +  5 consequences.

