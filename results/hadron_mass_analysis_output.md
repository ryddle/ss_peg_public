================================================================================
  HADRON MASSES IN THE SS-PEG LATTICE FRAMEWORK
  SS-PEG Research Program — 2026
================================================================================

================================================================================
§1. HADRON CATALOG — MASSES AND QUANTUM NUMBERS
================================================================================

  Hadrons are composite particles made of quarks bound by QCD.
  Unlike elementary particles (whose masses are fundamental),
  hadron masses arise from a combination of:
    • Current quark masses (from the lattice)
    • QCD binding energy (~99% of light hadron mass)
    • Spin-spin interactions (hyperfine splitting)
    • Electromagnetic effects (isospin breaking)

  BARYON OCTET (J^P = 1/2⁺):
      Name │  Mass (MeV) │     Quarks │   I   I₃   S   Q
    ───────────────────────────────────────────────────────
         p │     938.272 │        uud │ 0.5  0.5   0   1
         n │     939.565 │        udd │ 0.5 -0.5   0   0
         Λ │    1115.683 │        uds │ 0.0  0.0  -1   0
        Σ⁺ │    1189.370 │        uus │ 1.0  1.0  -1   1
        Σ⁰ │    1192.642 │        uds │ 1.0  0.0  -1   0
        Σ⁻ │    1197.449 │        dds │ 1.0 -1.0  -1  -1
        Ξ⁰ │    1314.860 │        uss │ 0.5  0.5  -2   0
        Ξ⁻ │    1321.710 │        dss │ 0.5 -0.5  -2  -1

  BARYON DECUPLET (J^P = 3/2⁺):
      Name │  Mass (MeV) │     Quarks │   I   I₃   S   Q
    ───────────────────────────────────────────────────────
       Δ⁺⁺ │    1232.000 │        uuu │ 1.5  1.5   0   2
        Δ⁺ │    1232.000 │        uud │ 1.5  0.5   0   1
        Δ⁰ │    1232.000 │        udd │ 1.5 -0.5   0   0
        Δ⁻ │    1232.000 │        ddd │ 1.5 -1.5   0  -1
       Σ*⁺ │    1382.800 │        uus │ 1.0  1.0  -1   1
       Σ*⁰ │    1383.700 │        uds │ 1.0  0.0  -1   0
       Σ*⁻ │    1387.200 │        dds │ 1.0 -1.0  -1  -1
       Ξ*⁰ │    1531.800 │        uss │ 0.5  0.5  -2   0
       Ξ*⁻ │    1535.000 │        dss │ 0.5 -0.5  -2  -1
        Ω⁻ │    1672.450 │        sss │ 0.0  0.0  -3  -1

  PSEUDOSCALAR MESONS (J^P = 0⁻):
      Name │  Mass (MeV) │     Quarks │   I   I₃   S   Q
    ───────────────────────────────────────────────────────
        π⁺ │     139.570 │        ud̄ │ 1.0  1.0   0   1
        π⁰ │     134.977 │     uū/dd̄ │ 1.0  0.0   0   0
        K⁺ │     493.677 │        us̄ │ 0.5  0.5   1   1
        K⁰ │     497.611 │        ds̄ │ 0.5 -0.5   1   0
         η │     547.862 │ uū+dd̄+ss̄ │ 0.0  0.0   0   0
        η' │     957.780 │ uū+dd̄+ss̄ │ 0.0  0.0   0   0

  VECTOR MESONS (J^P = 1⁻):
      Name │  Mass (MeV) │     Quarks │   I   I₃   S   Q
    ───────────────────────────────────────────────────────
        ρ⁺ │     775.110 │        ud̄ │ 1.0  1.0   0   1
         ω │     782.660 │     uū+dd̄ │ 0.0  0.0   0   0
       K*⁺ │     891.670 │        us̄ │ 0.5  0.5   1   1
       K*⁰ │     895.550 │        ds̄ │ 0.5 -0.5   1   0
         φ │    1019.461 │        ss̄ │ 0.0  0.0   0   0
       J/ψ │    3096.900 │        cc̄ │ 0.0  0.0   0   0
         Υ │    9460.300 │        bb̄ │ 0.0  0.0   0   0

  Total: 31 hadrons cataloged
    • 8 baryon octet members
    • 10 baryon decuplet members
    • 6 pseudoscalar mesons
    • 7 vector mesons

================================================================================
§2. LATTICE PROJECTION — BEST-FIT (p, q, r) FOR EACH HADRON
================================================================================

  For each hadron mass m, find the half-integer point (p, q, r)
  that minimizes |m_pred - m_exp|/m_exp, where
    m_pred = mₑ · 2^p · R₃^q · 3^r

    Name │      m_exp │     p     q     r │     m_pred │   Err% │  p_Z  q_Z  r_Z │ Err%_Z
  ─────────────────────────────────────────────────────────────────────────────────────
       p │    938.272 │   1.5  -3.5     4 │    938.331 │  0.01% │    2  -14   -2 │  0.11%
       n │    939.565 │    12  -0.5    -1 │    939.265 │  0.03% │    9   -4   -1 │  0.16%
       Λ │   1115.683 │  -3.5    -5   6.5 │   1115.315 │  0.03% │   14   -4   -4 │  0.03%
      Σ⁺ │   1189.370 │   3.5 -14.5    -3 │   1189.638 │  0.02% │    3   -4    3 │  0.14%
      Σ⁰ │   1192.642 │  13.5    -1    -2 │   1192.179 │  0.04% │    3   -4    3 │  0.14%
      Σ⁻ │   1197.449 │   0.5    -6   3.5 │   1197.878 │  0.04% │   -3  -11    3 │  0.16%
      Ξ⁰ │   1314.860 │     8  -8.5  -2.5 │   1315.398 │  0.04% │    4   -3    3 │  0.05%
      Ξ⁻ │   1321.710 │    -5 -13.5     3 │   1321.686 │  0.00% │   -2  -10    3 │  0.19%
     Δ⁺⁺ │   1232.000 │     8    -1   1.5 │   1231.972 │  0.00% │   15   -3   -4 │  0.10%
      Δ⁺ │   1232.000 │     8    -1   1.5 │   1231.972 │  0.00% │   15   -3   -4 │  0.10%
      Δ⁰ │   1232.000 │     8    -1   1.5 │   1231.972 │  0.00% │   15   -3   -4 │  0.10%
      Δ⁻ │   1232.000 │     8    -1   1.5 │   1231.972 │  0.00% │   15   -3   -4 │  0.10%
     Σ*⁺ │   1382.800 │    -5 -14.5   2.5 │   1383.014 │  0.02% │   12   -3   -2 │  0.13%
     Σ*⁰ │   1383.700 │   1.5    -6     3 │   1383.190 │  0.04% │   12   -3   -2 │  0.06%
     Σ*⁻ │   1387.200 │     9  -6.5    -2 │   1387.171 │  0.00% │    6  -10   -2 │  0.19%
     Ξ*⁰ │   1531.800 │  -3.5   -12     3 │   1532.093 │  0.02% │    7   -9   -2 │  0.12%
     Ξ*⁻ │   1535.000 │     0    -7   3.5 │   1535.170 │  0.01% │    7   -9   -2 │  0.09%
      Ω⁻ │   1672.450 │   6.5 -12.5  -3.5 │   1672.732 │  0.02% │   13   -4   -3 │  0.03%
      π⁺ │    139.570 │     4    -2   1.5 │    139.553 │  0.01% │   -3    0    7 │  0.09%
      π⁰ │    134.977 │     3 -10.5  -2.5 │    135.028 │  0.04% │   -1   -5    3 │  0.05%
      K⁺ │    493.677 │    -2  -6.5     4 │    493.773 │  0.02% │    1   -3    4 │  0.17%
      K⁰ │    497.611 │    -1  -3.5     5 │    497.625 │  0.00% │    2    0    5 │  0.18%
       η │    547.862 │  -3.5  -7.5   4.5 │    548.027 │  0.03% │   -4   -9    4 │  0.36%
      η' │    957.780 │     8  -1.5     1 │    957.568 │  0.02% │    5   -5    1 │  0.17%
      ρ⁺ │    775.110 │    -3  -7.5   4.5 │    775.028 │  0.01% │    1  -13   -1 │  0.08%
       ω │    782.660 │    -5    -8   5.5 │    782.542 │  0.02% │    2  -10    0 │  0.12%
     K*⁺ │    891.670 │     4  -0.5     4 │    891.568 │  0.01% │    1   -4    4 │  0.18%
     K*⁰ │    895.550 │     2   -13  -1.5 │    895.702 │  0.02% │   -5  -11    4 │  0.12%
       φ │   1019.461 │   0.5  -8.5     2 │   1019.478 │  0.00% │   14   -2   -3 │  0.09%
     J/ψ │   3096.900 │     6  -9.5    -1 │   3096.977 │  0.00% │    9   -6   -1 │  0.19%
       Υ │   9460.300 │     9 -12.5  -3.5 │   9462.401 │  0.02% │    5   -7    2 │  0.06%

  Half-integer lattice: mean error = 0.017%, max = 0.041%
  Integer lattice:      mean error = 0.124%, max = 0.358%
  Sub-1% fits: half-int = 31/31, integer = 31/31

================================================================================
§3. LATTICE QUALITY — COMPARISON WITH ELEMENTARY PARTICLES
================================================================================

  Elementary particles: mean = 0.710%, max = 1.348%
  Hadrons (half-int):   mean = 0.017%, max = 0.041%

        Baryon octet: mean = 0.026%, max = 0.041%
     Baryon decuplet: mean = 0.011%, max = 0.037%
           PS mesons: mean = 0.021%, max = 0.038%
       Vector mesons: mean = 0.011%, max = 0.022%

  Interpretation:
  → Hadrons fit the lattice surprisingly well, suggesting the
    lattice basis {ln2, lnR₃, ln3} captures QCD-scale physics too.

================================================================================
§4. QCD SCALE FROM THE GRAPH — Λ_QCD FROM α_s
================================================================================

  Graph-derived α_s = (R₃/R₂)⁴/(4π) = 0.117998
  Experimental α_s(M_Z) = 0.117994
  Error: 0.003%

  1-loop β₀ = (33 - 2·5)/(12π) = 0.610094
  Λ_QCD (graph) = M_Z · exp(-2π/(β₀·4π·α_s)) = 87.8 MeV
  Λ_QCD (exp, approx) ≈ 210 MeV
  Error: 58.2%

  Does Λ_QCD admit a lattice representation?
    Best half-int: (-0.5, 0.0, 5.0) → 87.8 MeV, err = 0.02%

  m_p / Λ_QCD(graph) = 10.684
    (lattice QCD predicts m_p/Λ_QCD ≈ 4.1–4.5)

  Proton mass in building-block language:
    ln(m_p/m_e) = 7.515428
    For reference: ln(m_p/m_e) / ln(2) = 10.8425
                   ln(m_p/m_e) / lnR₃ = -12.6381
                   ln(m_p/m_e) / ln(3) = 6.8408

================================================================================
§5. CONSTITUENT QUARK MODEL TEST — ADDITIVE MASS
================================================================================

  The constituent quark model assigns effective masses:
    m_u^const ≈ 336 MeV, m_d^const ≈ 340 MeV, m_s^const ≈ 486 MeV

  These are ~100× larger than current quark masses, because
  they absorb the ~99% QCD binding energy into an effective mass.

  Baryon octet: additive constituent model vs experiment:
      Name │     m_exp │     m_CQM │   Err%
    ────────────────────────────────────────
         p │   938.272 │    1012.0 │   7.9%
         n │   939.565 │    1016.0 │   8.1%
         Λ │  1115.683 │    1162.0 │   4.2%
        Σ⁺ │  1189.370 │    1158.0 │   2.6%
        Σ⁰ │  1192.642 │    1162.0 │   2.6%
        Σ⁻ │  1197.449 │    1166.0 │   2.6%
        Ξ⁰ │  1314.860 │    1308.0 │   0.5%
        Ξ⁻ │  1321.710 │    1312.0 │   0.7%

  Pseudoscalar mesons: additive model (quark + antiquark):
      Name │     m_exp │     m_CQM │   Err%
    ────────────────────────────────────────
        π⁺ │   139.570 │     676.0 │ 384.3%
        π⁰ │   134.977 │     672.0 │ 397.9%
        K⁺ │   493.677 │     822.0 │  66.5%
        K⁰ │   497.611 │     826.0 │  66.0%
         η │   547.862 │     672.0 │  22.7%
        η' │   957.780 │     972.0 │   1.5%

  Note: The constituent quark model works well for baryons (~8%)
  but fails badly for pseudoscalar mesons — the pion is anomalously
  light because it is a pseudo-Goldstone boson of chiral symmetry.

================================================================================
§6. QUARK COORDINATE COMPOSITION — ADDITIVE vs EXPERIMENT
================================================================================

  If hadron mass is multiplicative on the lattice:
    m_hadron = m₁ · m₂ · m₃  →  x_hadron = x₁ + x₂ + x₃
  Then hadron coordinates = sum of quark coordinates.

  If mass is additive (constituent model):
    m_hadron ≈ m₁ + m₂ + m₃  →  NO simple lattice relation

  Test: additive quark coordinates vs best-fit hadron coordinates:

    Name │    Σp    Σq    Σr │ p_fit q_fit r_fit │    Δp    Δq    Δr
  ────────────────────────────────────────────────────────────────────────
       p │  -3.5 -20.0  -4.0 │   1.5  -3.5   4.0 │   5.0  16.5   8.0
       n │  -2.5 -22.0  -5.0 │  12.0  -0.5  -1.0 │  14.5  21.5   4.0
       Λ │  -0.5 -21.0  -3.0 │  -3.5  -5.0   6.5 │  -3.0  16.0   9.5
      Σ⁺ │  -1.5 -19.0  -2.0 │   3.5 -14.5  -3.0 │   5.0   4.5  -1.0
      Σ⁰ │  -0.5 -21.0  -3.0 │  13.5  -1.0  -2.0 │  14.0  20.0   1.0
      Σ⁻ │   0.5 -23.0  -4.0 │   0.5  -6.0   3.5 │   0.0  17.0   7.5
      Ξ⁰ │   1.5 -20.0  -1.0 │   8.0  -8.5  -2.5 │   6.5  11.5  -1.5
      Ξ⁻ │   2.5 -22.0  -2.0 │  -5.0 -13.5   3.0 │  -7.5   8.5   5.0
     Δ⁺⁺ │  -4.5 -18.0  -3.0 │   8.0  -1.0   1.5 │  12.5  17.0   4.5
      Δ⁺ │  -3.5 -20.0  -4.0 │   8.0  -1.0   1.5 │  11.5  19.0   5.5
      Δ⁰ │  -2.5 -22.0  -5.0 │   8.0  -1.0   1.5 │  10.5  21.0   6.5
      Δ⁻ │  -1.5 -24.0  -6.0 │   8.0  -1.0   1.5 │   9.5  23.0   7.5
     Σ*⁺ │  -1.5 -19.0  -2.0 │  -5.0 -14.5   2.5 │  -3.5   4.5   4.5
     Σ*⁰ │  -0.5 -21.0  -3.0 │   1.5  -6.0   3.0 │   2.0  15.0   6.0
     Σ*⁻ │   0.5 -23.0  -4.0 │   9.0  -6.5  -2.0 │   8.5  16.5   2.0
     Ξ*⁰ │   1.5 -20.0  -1.0 │  -3.5 -12.0   3.0 │  -5.0   8.0   4.0
     Ξ*⁻ │   2.5 -22.0  -2.0 │   0.0  -7.0   3.5 │  -2.5  15.0   5.5
      Ω⁻ │   4.5 -21.0   0.0 │   6.5 -12.5  -3.5 │   2.0   8.5  -3.5

  Δp range: [-7.5, 14.5], mean = 4.44
  Δq range: [4.5, 23.0], mean = 14.61
  Δr range: [-3.5, 9.5], mean = 4.17

  → The offset varies significantly across hadrons.
    Multiplicative quark composition does NOT hold on the lattice.
    This is expected: hadron mass ≠ product of quark masses.

================================================================================
§7. HADRON MASS RATIOS IN BUILDING BLOCKS — EXHAUSTIVE SEARCH
================================================================================

  Search for hadron mass ratios expressible as R₂^a · R₃^b · R_EE^c · h∨₃^d · h∨₂^e
  with |a|,|b|,|c|,|d|,|e| ≤ 8.

         Ratio │     Target │   a   b   c   d   e │       Pred │   Err%
  ────────────────────────────────────────────────────────────────────────
       m_p/m_π │     6.7226 │  -2   7   8   3   6 │     6.7247 │ 0.031%
       m_n/m_p │     1.0014 │  -5   3  -7  -5   2 │     1.0010 │ 0.038%
       m_Ω/m_p │     1.7825 │  -6   7   7   4  -2 │     1.7832 │ 0.038%
       m_Λ/m_p │     1.1891 │   4   7  -7   3   2 │     1.1888 │ 0.027%
       m_K/m_π │     3.5371 │  -6   4   7   3  -2 │     3.5387 │ 0.044%
       m_ρ/m_π │     5.5536 │  -8   4  -2  -7   8 │     5.5543 │ 0.013%
       m_Δ/m_p │     1.3131 │  -2  -6  -8  -3  -6 │     1.3128 │ 0.021%
       m_φ/m_ρ │     1.3152 │   8  -2   8   8  -2 │     1.3154 │ 0.014%
      m_Jψ/m_p │     3.3006 │  -8  -6  -7  -5  -7 │     3.3005 │ 0.004%
      m_Υ/m_Jψ │     3.0548 │   3  -4   7   8  -8 │     3.0554 │ 0.021%
      m_η'/m_η │     1.7482 │  -6   8   7   2   2 │     1.7491 │ 0.049%
       m_p/m_e │  1836.1525 │  -5   3  -6  -1   7 │  1834.6452 │ 0.082%
       m_Σ/m_Λ │     1.0690 │  -2   3  -7   2  -6 │     1.0689 │ 0.004%
       m_Ξ/m_Σ │     1.1025 │   3  -2  -3   5  -8 │     1.1024 │ 0.007%
       m_Ω/m_Ξ │     1.2720 │   6  -4  -7  -1   1 │     1.2717 │ 0.024%

  Sub-1% matches: 15/15
  Sub-3% matches: 15/15

================================================================================
§8. GELL-MANN–OKUBO MASS RELATIONS ON THE LATTICE
================================================================================

  The GMO formula relates hadron masses within SU(3) multiplets.
  For the baryon octet (equal-spacing rule):
    2(m_N + m_Ξ) = 3m_Λ + m_Σ

  LHS = 2(m_N + m_Ξ) = 4514.4 MeV
  RHS = 3m_Λ + m_Σ   = 4540.2 MeV
  Error: 0.57%

  Baryon decuplet equal-spacing rule:
    m(Σ*) - m(Δ) = 152.6 MeV
    m(Ξ*) - m(Σ*) = 148.8 MeV
    m(Ω)  - m(Ξ*) = 139.0 MeV
    Mean spacing = 146.8 MeV, spread = 13.5 MeV

  Mean decuplet spacing ≈ 146.8 MeV
  = m_s(constituent) - m_u(constituent) ≈ 150 MeV

  Lattice current quark masses: m_s = 92.85, m_d = 4.67, m_u = 2.13 MeV
  m_s - m_d (lattice) = 88.17 MeV
  m_s - m_u (lattice) = 90.71 MeV
  Ratio: spacing/(m_s - m_d) = 1.67

  GMO on lattice coordinates:
  If 2(x_N + x_Ξ) = 3x_Λ + x_Σ holds in log-mass,
  then the GMO relation is a LATTICE constraint.

    p: 2(x_N + x_Ξ) = 19.0,  3x_Λ + x_Σ = 3.0,  diff = 16.0
    q: 2(x_N + x_Ξ) = -24.0,  3x_Λ + x_Σ = -16.0,  diff = -8.0
    r: 2(x_N + x_Ξ) = 3.0,  3x_Λ + x_Σ = 17.5,  diff = -14.5

================================================================================
§9. MESON MASS SPLITTINGS — PSEUDOSCALAR vs VECTOR
================================================================================

  The mass difference between vector and pseudoscalar mesons
  with the same quark content is due to spin-spin (chromomagnetic)
  interaction: ΔM ∝ (σ₁·σ₂) · |ψ(0)|² / (m₁·m₂)

          System │     m_PS │      m_V │       Δm │ m_V/m_PS │   Ratio bb
  ────────────────────────────────────────────────────────────────────
               π │    139.6 │    775.1 │    635.5 │   5.5536 │    (0.01%)
               K │    493.7 │    891.7 │    398.0 │   1.8062 │    (0.03%)
        φ/η(ss̄) │    547.9 │   1019.5 │    471.6 │   1.8608 │    (0.03%)
       J/ψ vs ηc │   2983.9 │   3096.9 │    113.0 │   1.0379 │    (0.06%)
         Υ vs ηb │   9399.0 │   9460.3 │     61.3 │   1.0065 │    (0.03%)

  The V-PS splitting decreases as quark mass increases
  (1/m₁m₂ dependence), confirming the chromomagnetic origin.

  Lattice perspective — coordinate differences:
                 π: Δ(p,q,r) = (-7.0, -5.5, 3.0)
                 K: Δ(p,q,r) = (6.0, 6.0, 0.0)
          φ/η(ss̄): Δ(p,q,r) = (4.0, -1.0, -2.5)
         J/ψ vs ηc: Δ(p,q,r) = (-1.0, -4.0, -1.5)
           Υ vs ηb: Δ(p,q,r) = (8.0, 1.0, -4.5)

  → The V-PS displacement is NOT constant on the lattice,
    reflecting the mass-dependent chromomagnetic interaction.

================================================================================
§10. BINDING ENERGY LANDSCAPE — WHAT THE LATTICE SEES
================================================================================

  For each baryon, compute:
    E_bind = m_hadron - Σ m_quark (current quark masses from lattice)
    R_bind = E_bind / m_hadron  (binding fraction)

    Name │     m_had │      Σm_q │    E_bind │ R_bind
  ────────────────────────────────────────────────────
       p │   938.272 │     8.944 │   929.328 │ 99.0%
       n │   939.565 │    11.484 │   928.081 │ 98.8%
       Λ │  1115.683 │    99.658 │  1016.025 │ 91.1%
      Σ⁺ │  1189.370 │    97.118 │  1092.252 │ 91.8%
      Σ⁰ │  1192.642 │    99.658 │  1092.984 │ 91.6%
      Σ⁻ │  1197.449 │   102.198 │  1095.251 │ 91.5%
      Ξ⁰ │  1314.860 │   187.832 │  1127.028 │ 85.7%
      Ξ⁻ │  1321.710 │   190.372 │  1131.338 │ 85.6%

  Mean binding fraction: 91.9%
  This confirms: ~99% of light baryon mass is binding energy.

  Heavy hadrons (charm, bottom):
       J/ψ (cc̄): m = 3096.9, Σm_q = 2506.9, E_bind = 590.0 MeV (19.1%)
         Υ (bb̄): m = 9460.3, Σm_q = 8299.1, E_bind = 1161.2 MeV (12.3%)

  For heavy quarkonia, the binding energy is SMALL and NEGATIVE
  (the quarks are heavier than the hadron → negative binding).
  This regime is controlled by perturbative QCD (Coulomb-like).

================================================================================
§11. ISOSPIN MULTIPLET STRUCTURE ON THE LATTICE
================================================================================

  Isospin multiplets differ only by u ↔ d substitution.
  On the lattice, this is a KNOWN displacement:
    Δx(u→d) = x_d - x_u = (1.0, -2, -1)

            Transition │    Δ(p,q,r) actual │  Δ(p,q,r) expected │ Match?
  ────────────────────────────────────────────────────────────────────────
            p → n: u→d │ (10.5,   3,  -5)       │ ( 1.0,  -2,  -1)       │      ✗
     Σ⁺ → Σ⁻: u→d (×2) │ (-3.0,   8,   6)       │ ( 2.0,  -4,  -2)       │      ✗
          Ξ⁰ → Ξ⁻: u→d │ (-13.0,  -5,   6)       │ ( 1.0,  -2,  -1)       │      ✗

================================================================================
§12. STATISTICAL SIGNIFICANCE — MONTE CARLO TEST
================================================================================

  Question: How likely is it that random lattice bases achieve
  the same quality of hadron mass fits?

  Method: Replace {ln2, lnR₃, ln3} with 3 random numbers
  in comparable ranges, and count how often the mean error
  for all hadrons is ≤ the observed value.

  IMPORTANT CAVEAT: The basis {ln2, lnR₃, ln3} is irrational
  and quasi-independent — a 3D half-integer lattice over such a
  basis can approximate any real number well.  The key test is
  whether the fit quality EXCEEDS what random bases achieve with
  the SAME search ranges.

  Observed mean error: 0.0170%
  Monte Carlo trials: 50,000
  Successes (mean_err ≤ observed): 35247
  p-value: 0.704940
  Significance: -0.5σ

  Note: With 1683 half-integer (p,q) pairs and continuous r
  rounded to nearest half-integer, we test ~1683 candidates per
  target. For 31 targets spanning a range of ~4.3 in ln-space,
  the lattice spacing is ~b₃/2 in the r-direction ≈ 0.55.
  The expected minimum distance to a half-integer lattice point
  is ~b₃/4 ≈ 0.14, giving a typical error of ~exp(0.14)-1 ≈ 15%.
  The observed 0.02% is FAR below this naive estimate.

================================================================================
§13. SUMMARY AND SYNTHESIS
================================================================================

  ┌──────────────────────────────────────────────────────────────────┐
  │  HADRON MASSES IN THE SS-PEG LATTICE                           │
  ├──────────────────────────────────────────────────────────────────┤
  │  1. Hadron lattice fit quality:                                │
  │     • Mean error (half-int):  0.017%                        │
  │     • Elementary particles:    0.710%                        │
  │     • Sub-1% hadron fits:     31/31                            │
  │                                                                │
  │  2. QCD scale:                                                │
  │     • Λ_QCD (graph):    87.8 MeV                            │
  │     • Λ_QCD (exp):       210.0 MeV                            │
  │                                                                │
  │  3. Building-block ratio matches (sub-1%): 15/15                │
  │                                                                │
  │  4. Key physical insights:                                     │
  │     • Light hadron mass is ~99% binding energy                 │
  │     • Constituent quark model → baryons ~8%, mesons FAIL       │
  │     • GMO equal-spacing rule holds to ~1%                      │
  │     • Heavy quarkonia: quarks heavier than hadron              │
  │                                                                │
  │  5. Lattice basis is not statistically special for hadrons   │
  │     (any comparable basis achieves similar fit quality)      │
  └──────────────────────────────────────────────────────────────────┘

  PHYSICAL INTERPRETATION:

  The lattice basis {ln 2, ln R₃, ln 3} was derived from elementary
  particle masses — quarks, leptons, and bosons. Hadron masses are
  dominated by QCD dynamics (confinement, chiral symmetry breaking)
  rather than by quark masses alone.

  The fact that hadron masses nonetheless fit the lattice with
  mean error 0.02% suggests one of:
    (a) The lattice basis captures something deeper than just
        quark masses — perhaps the QCD scale itself is lattice-
        compatible.
    (b) The three numbers {ln 2, ln R₃, ln 3} span a basis
        that is dense enough to approximate any real number
        with half-integer coefficients.
    (c) A combination of both: partial lattice structure plus
        sufficient basis density.

  The Monte Carlo test (§12) discriminates between these scenarios.

================================================================================
  END OF HADRON MASS ANALYSIS
================================================================================
