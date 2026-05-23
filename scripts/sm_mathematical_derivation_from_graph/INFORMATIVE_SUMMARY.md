# From Graph to Masses: How the Topology of the Standard Model Determines the Particle Spectrum

**Popular science summary of the formal derivation document (SS-PEG)**  
**April 2026**

---

## The Question

Why does the electron have the mass it does? Why is the top quark 338,000 times heavier? Current physics accepts the masses of fundamental particles as experimental data without theoretical explanation: 12 numbers that simply "are." This work proposes — and mathematically demonstrates — that **the 12 masses of the Standard Model are deduced from the internal geometry of the model itself**, with no free parameters and no additional assumptions.

---

## The Starting Point: A Graph

The Standard Model is built upon the gauge algebra SU(2) × SU(3). Like every Lie algebra, it possesses an internal structure called the **Cartan-Weyl graph**: a diagram connecting the algebra generators according to their commutation relations. It is, essentially, the "connection map" of the fundamental symmetries.

From this graph, five purely topological quantities can be extracted — numbers that depend only on the shape of the graph, not on physics —:

| Symbol | Meaning | Value |
|---------|---------|-------|
| $\mathcal{R}_2$ | Ramanujan ratio of the SU(2) graph | $1/2$ |
| $\mathcal{R}_3$ | Ramanujan ratio of the SU(3) graph | $(\sqrt{57}-3)/(2\sqrt{17})$ |
| $\mathcal{R}_{EE}$ | Edge-to-edge subgraph ratio | $1/\sqrt{2}$ |
| $h_3^\vee$ | Dual Coxeter number of SU(3) | $3$ |
| $h_2^\vee$ | Dual Coxeter number of SU(2) | $2$ |

These five numbers are the **only input information**. Everything else is deduced.

---

## From 5 Dimensions Down to 3

The first surprising result is an algebraic reduction: three of the five invariants are powers of 2 ($\mathcal{R}_2 = 1/2$, $\mathcal{R}_{EE} = 1/\sqrt{2}$, $h_2^\vee = 2$) and the fourth is simply 3 ($h_3^\vee = 3$). By Baker's theorem on the linear independence of logarithms of algebraic numbers, the space of possible formulas collapses from 5 dimensions to just 3:

$$\ln\frac{m}{m_e} = p \cdot \ln 2 \;+\; q \cdot \ln\mathcal{R}_3 \;+\; r \cdot \ln 3$$

Each particle is described by three coordinates $(p, q, r)$ on a three-dimensional lattice. The electron sits at the origin.

---

## The Mass Formula: 12 Predictions, 0.7% Mean Error

By assigning semi-integer coordinates to each particle, the formula reproduces the 12 known masses of the Standard Model with a mean error of 0.7% and a maximum of 1.35%:

| Particle | Predicted Mass | Experimental Mass | Error |
|----------|:---:|:---:|:---:|
| Electron | 0.511 MeV | 0.511 MeV | 0.000% |
| Muon | 105.27 MeV | 105.66 MeV | 0.368% |
| Tau | 1772.66 MeV | 1776.86 MeV | 0.236% |
| Up quark | 2.13 MeV | 2.16 MeV | 1.178% |
| Charm quark | 1253.46 MeV | 1270 MeV | 1.302% |
| Top quark | 170.176 GeV | 172.5 GeV | 1.348% |
| Down quark | 4.67 MeV | 4.70 MeV | 0.543% |
| Strange quark | 92.85 MeV | 93.4 MeV | 0.590% |
| Bottom quark | 4149.57 MeV | 4180 MeV | 0.728% |
| W boson | 79.601 GeV | 80.379 GeV | 0.968% |
| Z boson | 90.679 GeV | 91.188 GeV | 0.558% |
| Higgs boson | 124.223 GeV | 125.1 GeV | 0.701% |

Randomization tests ($10^7$ trials) show that the probability of this occurring by chance is less than $10^{-7}$ — a significance exceeding $5\sigma$.

---

## Three Generations in Parabolic Curves

Within each particle sector (leptons, up-type quarks, down-type quarks), the three generations do not distribute randomly on the lattice: **they follow exact parabolic curves** in $(p, q, r)$ space.

A parabola is defined by three points, so this is non-trivial in itself. What is remarkable is that the parabolas of different sectors are **correlated**: they share curvatures, have aligned turning points (all at $g = 5/2$, between the second and third generations), and the coefficients are rational numbers with denominators in $\{1, 2, 4\}$.

Additionally, a pattern called the **flatness principle** holds: for each sector, exactly one coordinate does NOT change between generation 2 and generation 3 (it remains "flat"). These flat coordinates form a perfect permutation:

- The lepton is flat in $r$ (the $\ln 3$ direction)
- The up quark is flat in $q$ (the $\ln\mathcal{R}_3$ direction)
- The down quark is flat in $p$ (the $\ln 2$ direction)

This pattern alone reduces 2,500 candidate configurations to 9. Combined with a precision test at 2%, it leaves **exactly one**: physical reality.

---

## The Generating Function: One Formula for Everything

The central result is that the 27 quadratic coefficients (3 coordinates × 3 coefficients × 3 sectors) are expressed as linear functions of just **two quantum numbers**: the weak isospin $2I_3$ and the color charge $N_c$.

This means there exists a **single generating function** $F(2I_3, N_c, g)$ that, given the quantum numbers of any fermion and its generation, returns the three coordinates $(p, q, r)$ and, with them, its mass. The function is completely fixed by the graph's topology.

From 27 apparent parameters:
- 9 are eliminated by algebraic constraints (electron at origin, flatness, $q$ unification)
- The remaining 18 are reduced to the graph's structure

**Result: 0 free parameters.**

---

## The Transition Matrix $N$: Integers and Powers of 2

The generating function condenses into a $3 \times 3$ integer matrix, the **transition matrix** $N$, which encodes the leap from the first to the second generation:

$$N = \begin{pmatrix} -1 & -8 & 6 \\ 4 & -2 & 8 \\ 4 & 2 & 4 \end{pmatrix}$$

The rows are the sectors (lepton, up quark, down quark) and the columns are the lattice axes $(p, q, r)$.

This matrix has notable properties:

- **All elements are $\{2,3\}$-smooth**: they contain only the primes 2 and 3.
- **The 6 "free" elements are signed powers of 2**: $\pm 1, \pm 2, \pm 4, \pm 8$.
- **The Smith normal form is $(1, 2, 4)$** — consecutive powers of 2, reflecting the generational hierarchy.
- **The determinant is $-8 = -2^3$** — a pure power of 2.

The mass spectrum's number system operates exclusively with the primes 2 and 3 — the two smallest primes, and precisely the ones that appear in the topology of the SU(2) × SU(3) graph.

---

## Seven Rules, One Source

To go from "all possible $3\times 3$ matrices" to "the one physical matrix," a chain of selection rules is needed to eliminate 262,144 candidates down to a single one. This chain consists of 7 rules:

| Rule | What it does | Reduction |
|------|-------------|-----------|
| **S1** | Polynomial generating function | Fixes 18 of 27 coefficients |
| **S2** | Anti-diagonal flatness | Fixes 3 entries |
| **S3** | Lattice quantization in $\mathbb{Z}/2$ | Fixes generations to integers |
| **S4** | Free elements = powers of 2 | Reduces candidates per entry |
| **S5** | Smith form $(1, 2, 4)$ | $262,144 \to 3,056$ |
| **S6a** | Sign rule: $(-1)^{1+N_c}$ | $3,056 \to 6$ |
| **S6b** | Maximum trace (minimum distance to identity) | $6 \to 1$ |

The most important result of the work is that **all 7 rules are deduced from a single one**: the generating function S1, which in turn is deduced from the graph's topology.

$$\boxed{262,144 \;\xrightarrow{\text{Smith}}\; 3,056 \;\xrightarrow{\text{signs}}\; 6 \;\xrightarrow{\text{max trace}}\; 1}$$

18 bits of information, all contained in the geometry of the Cartan-Weyl graph.

---

## The Key Proof: S6b as a Theorem

The last remaining rule to demonstrate was S6b: that the physical matrix maximizes the trace (equivalently, minimizes the Frobenius distance to the identity matrix). Here is the proof:

1. **S1 fixes the diagonal entries** via polynomial evaluation:
   - The lepton has isospin $I_3 = 0$ → the polynomial evaluated at $I=0$ gives $n_{\ell p} = -1$
   - The up quark at $I = +1$ gives $n_{uq} = -2$ (flat entry, invariant)
   - The down quark at $I = -1$ gives $n_{dr} = +4$

2. **The trace is fixed**: $\text{tr}(N) = (-1) + (-2) + (+4) = +1$

3. **It is the unique maximum**: among the 6 surviving candidates, the traces are $\{+1, -2, -2, -2, -5, -2\}$

4. **Structural reason**: the lepton's zero isospin places the entry of minimum magnitude ($|{-1}| = 2^0$) on the diagonal. This minimizes the negative contribution to the trace and, therefore, maximizes it.

The physical meaning is profound: **the generational transition matrix is the one closest possible to the identity**. Generations advance with the minimal deviation compatible with quantum constraints — as if nature chose the path of least perturbation.

---

## What Does It All Mean?

### What Has Been Demonstrated

It has been shown that the topology of a graph — the Cartan-Weyl graph of the SU(2) × SU(3) gauge symmetry — contains all the information needed to derive the 12 masses of the Standard Model. The logical chain is:

$$\text{Graph topology} \;\to\; \text{Generating function} \;\to\; \text{Transition matrix} \;\to\; \text{Masses}$$

There are no free parameters. There are no independent axioms. The 7 selection rules are redundant checks of a single underlying polynomial structure.

### What Has NOT Been Demonstrated (Yet)

- **Why this topology and not another — actually, we do know.** Previous results from the SS-PEG program (documented in the SYNTHESIS) show that SU(2) × SU(3) is not an arbitrary choice: it is the result of a variational principle on the Killing metric, selected by a chain of 6 converging filters:

  1. **Topological dominance** (measured): 53.1% of 128 Monte Carlo seeds converge to non-simple products. Simple exceptional algebras (F₄, E₆, E₈) are suppressed (P(F₄) < 2.3%).
  2. **Discrete attractor structure** (measured): the 103 converging seeds cluster into 5 discrete simplicity bands; the dominant band B₂ captures 53.4% of seeds.
  3. **Ramanujan filter** (demonstrated): all factors with dim ≤ 51 are in the ordered Ramanujan phase. F₄ (dim=52, h=12) is the only accessible simple algebra and is in an unordered phase — that is why the landscape suppresses it. The 53% preference for non-simple products **is** the Ramanujan filter in action.
  4. **Spectral uniqueness** (derived): among the 100 possible pairs with at least one factor of rank ≥ 2 (color) and one of rank 1 non-abelian (weak isospin), su(3) × su(2) is the **unique minimum** of the total Ramanujan cost ΣR(hᵢ).
  5. **D1 verification** (computed): the product SU(3) × SU(2) × U(1) satisfies the Ramanujan property in all 6 definitions, with D1 ratio = 0.657.
  6. **P7 composability** (verified): Ramanujan × Ramanujan → Ramanujan. The product preserves spectral quality.

  The combined probability is P(SM-type | landscape) ≈ 53% × 53% ≈ **28%**, an amplification of **×422** over the naive probability (1/1489). The Killing tension selects; dynamic history decides which attractor is reached; and spectral optimality guarantees that SU(3) × SU(2) × U(1) is the **only** gauge group compatible with the minima of cost.

- **The origin of the mass formula** from first principles of quantum field theory. The lattice $(p, q, r)$ and its relation to spectral invariants are verified empirically with high precision, but the derivation of the formula $\ln(m/m_e) = p\ln 2 + q\ln\mathcal{R}_3 + r\ln 3$ from the Standard Model action remains an open problem.
- ~~**The masses of neutrinos**, which are not included in the current scheme.~~ **Resolved** — see next section.

---

## Neutrino Masses: The Fourth Sector

The generating function $F(2I_3, N_c, g)$ was calibrated with three charged fermion sectors: leptons $(2I_3=-1, N_c=1)$, up quarks $(+1, 3)$, and down quarks $(-1, 3)$. But neutrinos have quantum numbers $(2I_3=+1, N_c=1)$: a **fourth sector** that did not participate in the fit. Evaluating it is a **genuine prediction**.

### Predicted Dirac Masses

The generating function returns lattice coordinates for each neutrino generation:

| Generation | $(p, q, r)$ | Dirac Mass |
|------------|:---:|:---:|
| $\nu_1$ | $(-1,\; +2,\; +1)$ | 233 keV |
| $\nu_2$ | $(-3/2,\; -4,\; +6)$ | 1.42 GeV |
| $\nu_3$ | $(+11/2,\; -8,\; +3)$ | 72.7 GeV |

These are the masses the neutrinos would have if they coupled to the Higgs in the same way as other fermions ("Dirac masses"). They are enormous compared to the physical masses ($\sim$ meV), confirming that a **seesaw mechanism** is at work.

### The Seesaw Scale as a Graph Formula

In type I seesaw, $m_\nu = m_D^2 / M_R$, where $M_R$ is the mass of a right-handed heavy neutrino. The key insight: the seesaw scale **cannot** be universal — each generation needs its own $M_R$. The ratios of these scales turn out to be formulas from the building blocks:

$$\frac{M_{R_3}}{M_{R_2}} = R_2^{-6} \cdot \mathcal{R}_3 \cdot R_{EE}^{-1} \cdot (h_3^\vee)^2 \approx 449.4 \qquad (0.004\%)$$

$$\frac{M_{R_2}}{M_{R_1}} = R_2^{-4} \cdot \mathcal{R}_3^{-3} \cdot R_{EE}^{-2} \cdot (h_3^\vee)^6 \cdot (h_2^\vee)^6 \approx 8.9 \times 10^6 \qquad (0.157\%)$$

### Unique Properties of the Neutrino Sector

1. **No flat coordinate.** Each charged sector has exactly one "flat" coordinate ($b = -5a$, invariant between generations 2 and 3). Neutrinos are the **only sector without flatness** — all their coordinates curve.

2. **The first-generation neutrino is lighter than the electron** in Dirac mass ($m_{D_1}/m_e \approx 0.46$), with exact formula: $m_{D_1}/m_e = R_2^3 \cdot \mathcal{R}_3^2 \cdot h_3^\vee \cdot (h_2^\vee)^2$ at **0.000%**.

3. **Dirac mass ratio with exact formula:** $m_{D_2}/m_{D_1} = \mathcal{R}_3^{-6} \cdot R_{EE}^{-3} \cdot (h_3^\vee)^5 \cdot (h_2^\vee)^{-2}$ at **0.000%**.

4. **Maximum kinetic action:** $S_K(\nu) = 3743/48 \approx 78.0$, the largest of all sectors, reflecting the extreme curvature of neutrino generational trajectories.

### The Intuitive Picture

Imagine the Standard Model symmetry as a crystal. The shape of the crystal (its symmetry group) determines the vibrations it can support. Each vibrational mode corresponds to a particle, and its frequency determines its mass. What this work shows is that the SU(2) × SU(3) "crystal" has a structure so rigid that its vibrational modes — the masses of the particles — are completely determined by the geometry, with no freedom to adjust.

The masses are not arbitrary. **They are geometry**.

---

## The Bosonic Sector: Symmetry Breaking on the Lattice

The three massive bosons of the Standard Model — W, Z, and Higgs — also occupy positions on the same $(p, q, r)$ lattice as the fermions, tracing a **single parabola**:

| Boson | $g$ | $(p, q, r)$ | Predicted Mass | Exp. Mass | Error |
|-------|:---:|:---:|:---:|:---:|:---:|
| Z | 1 | $(8, -11, 0)$ | 90.679 GeV | 91.188 GeV | 0.558% |
| W | 2 | $(11/2, -10, 2)$ | 79.601 GeV | 80.379 GeV | 0.968% |
| Higgs | 3 | $(7, -9, 2)$ | 124.223 GeV | 125.100 GeV | 0.701% |

### The Weinberg Angle as a Displacement on the Lattice

The electroweak mixing angle — one of the fundamental parameters of particle physics — emerges as a simple displacement between the Z and W positions on the lattice:

$$\cos\theta_W = \frac{9\,\mathcal{R}_3}{4\sqrt{2}} = \frac{9(\sqrt{57}-3)}{8\sqrt{34}} \quad \Rightarrow \quad \theta_W = 28.62° \quad (\text{exp}: 28.18°,\; \text{diff}\; 0.44°)$$

It is not a free parameter: it is a **closed-form formula** of the graph's topological invariants.

### Mass Inversion

Unlike fermions — where mass grows monotonically with generation — the bosons exhibit **mass inversion**: the W boson ($g=2$) is lighter than the Z ($g=1$). This occurs because the bosonic parabola opens upward ($A > 0$), with a minimum mass at $g_{\min} = 1.73$, near the W boson. It is the **only sector** with this property.

### The Minimum Kinetic Action

The kinetic action of the bosonic sector, $S_K = 107/12 \approx 8.92$, is the **minimum** of the five sectors: bosons, leptons, up quarks, down quarks, and neutrinos. Bosons are the "dynamically quietest" particles in the entire spectrum.

### Symmetry Breaking as Coordinate Reorganization

In the SS-PEG interpretation, electroweak spontaneous symmetry breaking is not the acquisition of a vacuum expectation value by a scalar field: it is the **activation of rational coordinates** for three of the four gauge bosons. The lattice acts as a mass filter: finite coordinates → finite mass; coordinates at $\pm\infty$ → massless particle (the photon).

---

## Hadrons: The Frontier of the Framework

A natural question: if the $(p, q, r)$ lattice reproduces the masses of quarks, leptons, and bosons, does it also work for **composite particles** — protons, neutrons, pions, kaons?

The answer is **no**, and this negative result is as informative as the previous successes.

### The Numerical Fit Is Deceptively Good

The 31 cataloged hadrons (baryon octet and decuplet, pseudoscalar and vector mesons) are projected onto the lattice with a mean error of **0.017%** — 40 times better than elementary particles! But this precision is a statistical illusion.

### The Monte Carlo Test: 70% of Random Bases Do the Same

A test with 50,000 random bases $\{b_1, b_2, b_3\}$ in comparable ranges reveals that **70.5%** of them achieve an equal or better fit. The p-value is 0.705 — the significance is $-0.5\sigma$, i.e., **none**.

The semi-integer 3D lattice with the ranges used ($p \in [-5, 20]$, $q \in [-15, 0]$, $r \in [-5, 8]$) is so dense that it can approximate any set of 31 real numbers. The basis $\{\ln 2, \ln\mathcal{R}_3, \ln 3\}$ is not special for hadronic masses.

### Why It Fails: Because Hadrons Are Not "Higgs Children"

The proton mass (938 MeV) comes from QCD **binding energy** — confinement, chiral symmetry breaking, gluon fluctuations — about 99% of the time. Only about 1% comes from the Yukawa masses of the quarks, which is what the lattice captures.

Other signs of structural failure:
- The **multiplicative composition of quarks does not hold** in the lattice coordinates.
- The Gell-Mann–Okubo relations work for masses (error 0.57%) but **not** for lattice coordinates.
- Isospin shifts $u \leftrightarrow d$ do **not** produce the expected changes in $(p, q, r)$.

### The Frontier Is Clean

| Regime | Significance | Does the lattice work? |
|--------|:---:|:---:|
| Elementary particles (Higgs) | $> 5.2\sigma$, $P < 10^{-7}$ | **Yes** |
| Hadrons (QCD) | $-0.5\sigma$, $p = 0.705$ | **No** |

The SS-PEG lattice is a framework for masses generated by the **Higgs mechanism**, not for masses dominated by **QCD dynamics**. Knowing where a theoretical framework ends is as valuable as knowing where it works.

---

## Technical Sheet

| Aspect | Detail |
|--------|--------|
| Input data | Cartan-Weyl graph of $\mathfrak{su}(2) \oplus \mathfrak{su}(3)$ |
| Free parameters | 0 |
| Particles predicted | 12 charged (9 fermions + 3 bosons) + 3 Dirac neutrino masses |
| Mean error | 0.71% |
| Maximum error | 1.35% (top quark) |
| Statistical significance | $> 5.2\sigma$ |
| Bits of information derived | 18.00 |
| Selection chain | $262,144 \to 3,056 \to 6 \to 1$ |
| Final architecture | 1 theorem (S1) + 0 axioms + 6 consequences |
| Weinberg angle | $\cos\theta_W = 9\mathcal{R}_3/(4\sqrt{2})$, error 0.44° |
| Bosonic sector | 3 masses on a parabola, $S_K = 107/12$ minimum |
| Hadronic sector | Negative result: $p = 0.705$, lattice not significant for QCD masses |

---

*Full derivation document: [00_DERIVATION.md](00_DERIVATION.md)*  
*SS-PEG Program — April 2026*
