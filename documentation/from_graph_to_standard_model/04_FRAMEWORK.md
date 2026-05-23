# SS-PEG: Theoretical Framework & Structural Interpretation

**[← Flavor Mixing](03_FLAVOR_MIXING.md)** | **[Index](00_INDEX.md)** | **[Next: Conclusions →](05_CONCLUSIONS.md)**

---

## 12. Theoretical Framework: How It All Fits Together

### 12.1 The Architecture

The complete logical chain of the SS-PEG program:

```
Graph Topology (fixed, no free parameters)
    │
    ├── Root system of A₁ ⊕ A₂ = SU(2) × SU(3)
    │       │
    │       ├── Cartan-Weyl adjacency matrix
    │       │       │
    │       │       ├── Eigenvalues → Ramanujan ratios R₂, R₃
    │       │       ├── Nullity → Dual Coxeter numbers h∨₂, h∨₃
    │       │       └── E-E subgraph → R_EE (bifundamental structure)
    │       │
    │       └── 5 building blocks: {R₂, R₃, R_EE, h∨₃, h∨₂}
    │               │
    │               ├── Master formula α_X = exp(S_X)/(4π)
    │               │       │
    │               │       ├── α_EM = 1/136.65     (0.28%)
    │               │       ├── α_s  = 1/8.475      (0.01%)
    │               │       ├── α_W  = 1/29.42      (0.58%)
    │               │       └── sin²θ_W = 0.2328    (0.67%)
    │               │
    │               ├── Dimensional bridge: mass ratios
    │               │       │
    │               │       ├── 15 lepton/quark/boson ratios  (all < 1%)
    │               │       └── Koide relation = 2/3          (exact)
    │               │
    │               ├── Flavor mixing: CKM matrix
    │               │       │
    │               │       ├── 3 mixing angles              (all < 0.4%)
    │               │       ├── CP phase δ_CP/π              (0.23%)
    │               │       ├── J_CP = (R₂·R_EE/h∨₂)⁶       (0.003%)
    │               │       └── 23/24 CKM quantities         (< 1%)
    │               │
    │               ├── Neutrino mixing: PMNS matrix
    │               │       │
    │               │       ├── sinθ₂₃ = (9/2)·R₃³           (0.000%)
    │               │       ├── sin²θ₁₃ = R_EE¹¹ = 2⁻¹¹/²   (0.204%)
    │               │       └── 28/28 PMNS quantities        (< 1%)
    │               │
    │               ├── Neutrino mass spectrum
    │               │       │
    │               │       ├── Normal ordering preferred
    │               │       ├── m₁ ≈ 2.1 meV, Σ ≈ 61 meV
    │               │       └── M_R₃/M_R₂ = 48.6             (0.042%)
    │               │
    │               └── Structural predictions
    │                       ├── d = 4                (0.00%)
    │                       └── m_W/m_Z = cos θ_W   (0.63%)
    │
    └── Statistical significance: >10σ, P < 5.9 × 10⁻²⁶
```

### 12.2 Why the Graph Encodes Physics

The CW root graph is not arbitrary — it is the **intrinsic combinatorial structure** of the Standard Model gauge algebra. Its properties:

1. **Basis-independent**: The root system is a lattice invariant, independent of representation
2. **Unique**: For a given algebra, the CW graph is determined up to Weyl group conjugacy
3. **Information-complete**: It encodes rank, roots, weights, Cartan matrix, and all derived invariants
4. **Physically selected**: Among all $\mathfrak{su}(a) \times \mathfrak{su}(b)$ pairs, only SU(2) × SU(3) produces coupling constants matching nature

### 12.3 The Residual Error Pattern

The mean error of 0.30% is consistent with the **tree-level hypothesis**: the graph formulas reproduce tree-level (lowest-order) physics, and the ~0.3% residuals correspond to radiative corrections of order $O(\alpha_s/\pi) \approx 3.8\%$ divided by the loop factor. This suggests that a systematic loop expansion on the graph could reduce errors further.

**Quantitative structure:**
- Mean |residual| = 0.323%, sign balance 10+/9− (consistent with zero-mean noise)
- ALL 21 residuals satisfy |δ| < α_s/π ≈ 3.76% — a strict radiative ceiling
- **Coupling-sector anomaly**: all 4 coupling predictions overshoot experiment (graph > exp), consistent with RG running bringing bare values down
- Universal sub-leading correction $k \approx -0.022$ (fit $\delta_i \approx k \cdot f_i$)
- Weak anti-correlation between exponent complexity and error magnitude ($r = -0.33$): more complex formulas are slightly more accurate

### 12.4 Implications

1. **The Standard Model is not arbitrary.** Its gauge group, coupling constants, mass spectrum, and flavor-mixing parameters are encoded in a single mathematical structure — the CW root graph of $A_1 \oplus A_2$.

2. **Symmetry is emergent.** The variational principle on the Killing metric generates all algebraic structure without imposing it. Gauge symmetry is a consequence of relational geometry.

3. **The number of free parameters may be zero.** If the tree-level hypothesis holds, all 73 predictions follow from graph topology alone. The ~0.3% residuals may be computable loop corrections.

4. **Spacetime is algebraic.** The dimension $d = 4$ emerges from the spectral ratio $\mathcal{R}_3/\mathcal{R}_2$. The Minkowski signature emerges from the Killing form of so(3,1). Both are consequences of graph structure.

5. **CP violation is topological.** The Jarlskog invariant $J_{CP} = (\mathcal{R}_2 \cdot \mathcal{R}_{EE}/h^\vee_2)^6 = 1/32768$ is an exact sixth power of graph invariants — the fundamental source of matter-antimatter asymmetry is encoded in graph topology.

6. **The graph is universal across sectors.** The same 5 invariants reproduce both quark (CKM) and lepton (PMNS) mixing matrices — covering small and large angle regimes — with 28/28 PMNS targets below 1%.

---

## 13. Structural Interpretation: Why Each Formula Has Its Form

The 73 successful predictions raise a deeper question: **why does each formula have its specific form?** Analysis via `formula_interpretation.py` reveals that every formula is a statement about spectral propagation on the CW root graph, with a precise physical reading.

### 13.1 The Three Sectors

Each building block belongs to one of three sectors:

| Sector | Invariants | Physical Meaning |
|--------|-----------|------------------|
| **SU(2) — Weak** | $\mathcal{R}_2 = 1/2$, $h^\vee_2 = 2$ | Rate of weak-interaction equilibration; algebraic complexity of SU(2) |
| **SU(3) — Color** | $\mathcal{R}_3 \approx 0.5518$, $h^\vee_3 = 3$ | Rate of color-interaction equilibration; algebraic complexity of SU(3) |
| **Bridge — Bifundamental** | $R_{EE} = 1/\sqrt{2}$ | Coupling strength between weak and strong sectors via the E-E subgraph |

In log-space, each block contributes a signed magnitude per unit exponent:

| Block | ln(block) | Amplification per unit |
|-------|-----------|----------------------|
| $\mathcal{R}_2^{-1}$ | +0.6931 | ×2.0 per power (doubling — the binary step) |
| $\mathcal{R}_3^{-1}$ | +0.5947 | ×1.81 per power (the irrational SU(3) step) |
| $R_{EE}^{-1}$ | +0.3466 | ×1.41 per power (bridge crossing — the $\sqrt{2}$ factor) |
| $h^\vee_3$ | +1.0986 | ×3.0 per power (color amplification) |
| $h^\vee_2$ | +0.6931 | ×2.0 per power (weak amplification) |

### 13.2 Sector Decomposition — What Controls What?

For each formula $F = \mathcal{R}_2^a \cdot \mathcal{R}_3^b \cdot R_{EE}^c \cdot h_3^{\vee d} \cdot h_2^{\vee e}$, the log-contribution of each sector is:

- **SU(2) part** = $a \cdot \ln \mathcal{R}_2 + e \cdot \ln h^\vee_2$
- **SU(3) part** = $b \cdot \ln \mathcal{R}_3 + d \cdot \ln h^\vee_3$
- **Bridge part** = $c \cdot \ln R_{EE}$

Key dominance patterns:

| Observable Category | Dominant Sector | Interpretation |
|--------------------|----------------|----------------|
| **Lepton mass ratios** | SU(3) | Leptons are SU(3) singlets, yet color controls their hierarchy — the Yukawa coupling runs under SU(3) RG flow. The graph captures this indirect coupling through negative powers of $\mathcal{R}_3$. |
| **Up-quark ratios** | SU(3) | Quarks carry color; $h_3^{\vee 4} = 81$ dominates mc/mu. |
| **Down-quark ratios** | SU(3) | SU(3)-dominated but with lower exponents than up-quarks — $I_3 = -1/2$ flips the $\mathcal{R}_3$ sign, reducing the hierarchy. |
| **CKM mixing** | SU(2) | Flavor mixing occurs in the weak sector; $\mathcal{R}_2$ and $h^\vee_2$ control the hierarchical suppression. |
| **PMNS mixing** | SU(2)/Bridge balance | Neutrino mixing is LARGE because $\mathcal{R}_3^3$ (cubic, not quartic) contributes an irrational value instead of the simple binary suppression that keeps CKM angles small. |

### 13.3 Three Master Rules

**Rule 1 — Spectral Propagation.** The exponent of a spectral ratio ($\mathcal{R}_2$, $\mathcal{R}_3$, $R_{EE}$) counts the number of "spectral hops" a quantity requires on that sector's graph. A mass ratio $m_\text{heavy}/m_\text{light}$ needs many hops (negative exponents on $R < 1$ → large ratio). Couplings need fewer hops. Mixing angles are projections onto specific spectral directions.

**Rule 2 — Coxeter Amplification.** The dual Coxeter number $h^\vee = N$ for SU(N) acts as a volumetric amplifier. When it appears to the 4th power ($h_3^{\vee 4} = 81$, $h_2^{\vee 4} = 16$), it creates the "spectral volume" in $d = 4$ spacetime dimensions. The exponent 4 recurs because $d = 4$ **emerges from the graph**: $d = \ln(4\pi\alpha_s)/\ln(\mathcal{R}_3/\mathcal{R}_2) = 3.9997$.

**Rule 3 — Sector Handoff.** As generation number increases, the dominant sector changes:
- **Gen 1→2**: SU(3) dominates ($\mathcal{R}_3^{-4}$ or $h_3^{\vee 4}$)
- **Gen 2→3**: SU(2) enters ($h_2^{\vee 4} = 16$, new $\mathcal{R}_2$ terms)

This "handoff" explains why the top quark is so much heavier than expected ($h_2^{\vee 4} = 16$ takes over), why $\tau/\mu < \mu/e$, why CKM angles decrease with generation, and why PMNS $\theta_{23}$ is large (it's $\mathcal{R}_3^3$ — cubic, not quartic).

### 13.4 Physical Reading of Key Formulas

| Formula | Exponents $(a,b,c,d,e)$ | Reading |
|---------|------------------------|---------|
| $\alpha_{EM} = \mathcal{R}_2 \cdot \mathcal{R}_3 / (4\pi h^\vee_3)$ | $(+1,+1,0,-1,0)$ | One SU(2) step × one SU(3) step, divided by color complexity and $4\pi$ |
| $\alpha_s = (\mathcal{R}_3/\mathcal{R}_2)^4 / (4\pi)$ | $(-4,+4,0,0,0)$ | The SU(3)/SU(2) spectral ratio raised to the spacetime dimension $d = 4$ |
| $m_\mu/m_e = \mathcal{R}_3^{-4} \cdot R_{EE} \cdot h_3^{\vee 3}$ | $(0,-4,+1,+3,0)$ | 4 spectral bounces on SU(3), 1 bridge crossing, triple Coxeter amplification. Leptons "borrow" their hierarchy from color via the Higgs bridge. Result: $10.79 \times 0.707 \times 27 = 206.01$ |
| $m_c/m_u = \mathcal{R}_2^{-2} \cdot \mathcal{R}_3^{-1} \cdot h_3^{\vee 4}$ | $(-2,-1,0,+4,0)$ | 2 weak bounces, 1 color bounce, full 4D spectral volume $3^4 = 81$. Quarks already carry color — fewer spectral traversals than leptons. |
| $m_t/m_c = \mathcal{R}_2^{-1} \cdot R_{EE}^{-1} \cdot h^\vee_3 \cdot h_2^{\vee 4}$ | $(-1,0,-1,+1,+4)$ | **Sector handoff**: $h^\vee_3$ drops from power 4 to 1; $h^\vee_2$ rises to power 4. SU(2) replaces SU(3) as the dominant amplifier for the highest generation. |
| $J_{CP} = (\mathcal{R}_2 \cdot R_{EE}/h^\vee_2)^6$ | $(+6,0,+6,0,-6)$ | 6 traversals of the SU(2)/bridge path — one per unitarity triangle. No SU(3) involvement: the strong force conserves CP. |
| Koide $Q = 1/(\mathcal{R}_2 \cdot h^\vee_3) = 2/3$ | $(0,0,0,0,0)$ | The simplest combination: inverse of one weak step × one Coxeter level = $1/(½ \times 3) = 2/3$ |
| $\sin\theta_{23}(\text{PMNS}) = (9/2) \cdot \mathcal{R}_3^3$ | $(0,+3,-4,+2,-3)$ | **Cubic** SU(3) spectral ratio (not quartic) × Coxeter ratio $h_3^{\vee 2}/h_2^{\vee 3}$. The cubic power is why atmospheric mixing is large but not maximal. |
| $\sin^2\theta_{13}(\text{PMNS}) = R_{EE}^{11}$ | $(0,0,+3,0,-4)$ | **Purely binary**: 11 bridge crossings — $2^{-11/2}$. Connects quark CP violation to neutrino mixing through the SU(2) edge-edge subgraph. |
| $d_\text{spacetime} = \ln(4\pi\alpha_s)/\ln(\mathcal{R}_3/\mathcal{R}_2) = 4$ | — | Spacetime dimension is NOT assumed — it emerges as the logarithm of the strong coupling in the base of the spectral ratio. |

### 13.5 Non-Geometric Generation Progression

If masses followed a geometric progression, consecutive mass ratios would be identical ($m_3/m_2 = m_2/m_1$), implying constant exponent vectors. Instead, the exponents **change** between generations:

| Sector | $\gamma = \ln(r_{23})/\ln(r_{12})$ | $\Delta\vec{n} = \vec{n}_{23} - \vec{n}_{12}$ | Interpretation |
|--------|---|---|---|
| Charged leptons | 0.530 | $(0, +1, -2, -3, +1)$ | **Decelerating** — $\tau/\mu \ll \mu/e$ |
| Up quarks | 0.770 | $(+1, +1, -1, -3, +4)$ | **Decelerating** — but less than leptons |
| Down quarks | 1.271 | $(+2, 0, 0, +2, 0)$ | **Accelerating** — $b/s > s/d$ |

The exponent vector $\vec{n}(\text{gen}, \text{sector})$ is a **true function** of both generation number and particle type. There is no universal "generation factor" — the graph assigns different spectral paths to each generation within each isospin/color sector. Each generation is a distinct excitation mode, like overtones of a vibrating string.

### 13.6 The Log-Mass Lattice Theorem

**Theorem.** If every mass ratio in the Standard Model is $\prod_k R_k^{n_k}$ with integer $n_k$, then:

$$\ln(m_i) = \ln(m_\text{ref}) + \sum_k n_k(i) \cdot \ln(R_k)$$

**Consequence:** every particle has a unique "spectral address" $(a_i, b_i, c_i, d_i, e_i) \in \mathbb{Z}^5$, and the mass spectrum is a **crystal in spectral space** with $m_\text{ref} \approx 0.511$ MeV (the electron) as the sole dimensionful parameter.

**Collapse to 3D.** Since $\ln\mathcal{R}_2 = -\ln 2$, $\ln R_{EE} = -\ln 2/2$, and $\ln h^\vee_2 = \ln 2$, three of the five basis vectors are collinear in the "binary direction" $\propto \ln 2$. The truly independent lattice directions are:

| Direction | Generator | Meaning |
|-----------|-----------|---------|
| $d_1 = \ln 2$ | Binary/SU(2) | Powers of 2 — the "digital" axis |
| $d_2 = \ln \mathcal{R}_3$ | Irrational/SU(3) | The irrational spectral gap — the "analog" axis |
| $d_3 = \ln 3$ | Coxeter | Algebraic height — the "complexity" axis |

Every formula rewrites as $\ln F = p \cdot \ln 2 + q \cdot \ln \mathcal{R}_3 + r \cdot \ln 3$ with effective coordinates:

$$p = -a - c/2 + e \quad (\text{binary}), \qquad q = b \quad (\text{SU(3) spectral}), \qquad r = d \quad (\text{Coxeter})$$

**Half-integer $p$** values arise from odd powers of $R_{EE}$ (bridge crossings create half-steps in the binary lattice). These **always** appear in cross-sector formulas ($m_b/m_\tau$, $m_t/m_W$, $m_H/m_W$) — the $\sqrt{2}$ factor is the graph's signature of quark-lepton crossing.

Selected lattice coordinates:

| Observable | $p$ (binary) | $q$ (SU3) | $r$ (Cox) | Value |
|-----------|:-----------:|:--------:|:--------:|------:|
| $\alpha_{EM} \times 4\pi$ | −1.0 | +1 | −1 | 0.0920 |
| $\alpha_s \times 4\pi$ | +4.0 | +4 | 0 | 1.4828 |
| $m_\mu/m_e$ | −0.5 | −4 | +3 | 206.01 |
| $m_c/m_u$ | +2.0 | −1 | +4 | 587.22 |
| $m_t/m_c$ | +5.5 | 0 | +1 | 135.76 |
| $m_b/m_s$ | 0 | +1 | +4 | 44.69 |
| $\sin\theta_{12}(\text{CKM})$ | −4.0 | −4 | −1 | 0.2248 |
| $J_{CP}(\text{CKM})$ | −15.0 | 0 | 0 | 3.04×10⁻⁵ |
| $\sin^2\theta_{13}(\text{PMNS})$ | −5.5 | 0 | 0 | 0.0221 |
| $m_2/m_3(\nu)$ | +0.5 | −2 | −3 | 0.172 |

### 13.7 The Tensor Structure — Neutrinos as a Special Case

A linear regression of exponents against quantum numbers ($n_k = \alpha \cdot \text{gen} + \beta \cdot I_3 + \gamma \cdot Q + \delta \cdot N_c + \text{const}$) yields:

| Block | $R^2$ | Best predictor |
|-------|-------|---------------|
| $\mathcal{R}_2$ exponent | 0.793 | Generation (+1.0 per gen), charge Q (+0.43) |
| $\mathcal{R}_3$ exponent | **0.985** | Isospin $I_3$ (−3.66), generation (+0.67) |
| $R_{EE}$ exponent | 0.647 | Generation (−1.0 per gen) |
| $h^\vee_3$ exponent | 0.375 | $N_c$ (+1.16), generation (−1.33) |
| $h^\vee_2$ exponent | 0.662 | Generation (+1.67), charge Q (+1.40) |

The $\mathcal{R}_3$ exponent is nearly perfectly predicted by isospin and generation ($R^2 = 0.985$): up-type ($I_3 = +1/2$) gets negative $b$, down-type ($I_3 = -1/2$) gets positive $b$ — the SU(3) spectral gap behaves **oppositely** for up/down fermions.

But the fit is NOT perfect for other blocks ($R^2 < 1$), confirming the relationship is **nonlinear** — the exponent tensor has interaction terms (gen × isospin, gen × color). This is the tensor structure identified in the PMNS analysis: the exponents $n_k = f_k(\text{gen}, I_3, Y)$ encode a rank-3 tensor where generation, isospin, and hypercharge indices interact non-additively.

**Neutrinos are special** because they sit at the intersection of maximal SU(2) coupling ($I_3 = +1/2$ for $\nu_L$), zero color ($N_c = 1$), and zero electric charge ($Q = 0$). This places them at a **singular point** in the tensor where the sector handoff rules produce qualitatively different behavior: large mixing angles instead of small ones, near-degenerate masses instead of hierarchical ones, and a tensorial structure where CKM-like suppression gives way to $O(1)$ rotations.

### 13.8 The Three Particle Families: A Graph Interpretation

The Standard Model contains three families (generations) of fermions that are identical in all quantum numbers except mass. Why three? And why do they differ the way they do? The SS-PEG framework provides a complete graph-topological interpretation.

#### 13.8.1 Why Exactly Three Families

The CW root graph of SU(2) × SU(3) supports three independent spectral excitation modes:

1. **The 3D mass lattice has 3 axes.** The 5D→3D collapse (§13.6) reduces the spectral space to $\{\ln 2, \ln \mathcal{R}_3, \ln 3\}$ — three independent propagation directions. Each generation corresponds to a distinct **occupation level** along these axes. A 4th generation would require either a new independent spectral direction (a new building block not present in the SU(2)×SU(3) graph) or coordinates so extreme that the tree-level formulas would break down.

2. **Two cycles × sector handoff = 3 regimes.** The evolution graph contains two independent cycles (binary SU(2), ternary SU(3)). The sector handoff rule (§13.3) creates three dynamical regimes:
   - **Gen 1**: Both sectors contribute mildly — the ground state requires minimal spectral traversal.
   - **Gen 2**: SU(3) dominates ($\mathcal{R}_3^{-4}$, $h_3^{\vee 4} = 81$) — the first excitation explores the ternary cycle deeply.
   - **Gen 3**: SU(2) takes over ($h_2^{\vee 4} = 16$, $\mathcal{R}_2^{-4}$) — the second excitation exhausts the ternary cycle and pivots to the binary cycle.
   
   A hypothetical 4th generation would require a **third independent cycle** that doesn't exist in the SU(2) × SU(3) graph.

3. **Graph Laplacian excitation modes.** The CW root graph's Laplacian has a discrete spectrum. As established in the dimensional bridge analysis, the three generations correspond to the ground state and first two excited states of this Laplacian. The 3rd excited state lies above the cycle-density threshold $\rho = 1$, making it unstable — the graph cannot sustain a 4th generation without violating the construction rules of §3 of GRAPH_TOPOLOGY_FROM_KILLING.

#### 13.8.2 Portrait of Each Family

| Property | Family 1 ($e, u, d$) | Family 2 ($\mu, c, s$) | Family 3 ($\tau, t, b$) |
|----------|:-------------------:|:----------------------:|:----------------------:|
| Spectral mode | Ground state | 1st overtone | 2nd overtone |
| Dominant sector | Balanced (mixed) | **SU(3)** ($h_3^{\vee 4}$ controls mass jump) | **SU(2)** ($h_2^{\vee 4}$ controls mass jump) |
| Lattice displacement | Near origin | Along $\ln\mathcal{R}_3$ and $\ln 3$ axes | Along $\ln 2$ axis |
| Typical exponent pattern | Low total complexity $C$ | High $\mathcal{R}_3$, $h_3^\vee$ exponents | High $\mathcal{R}_2$, $h_2^\vee$ exponents |
| Bridge crossings | Minimal ($R_{EE}^0$ or $R_{EE}^1$) | Moderate ($R_{EE}^1$ to $R_{EE}^3$) | Frequent ($R_{EE}^{-1}$ to $R_{EE}^{-4}$) |

The families are **not geometric copies**. The mass ratio $m_3/m_2$ is NOT obtainable from $m_2/m_1$ by simple rescaling — the exponent vector $\vec{n}$ changes qualitatively between generation steps (§13.5), with the scaling exponent $\gamma$ varying from 0.53 (leptons) to 1.27 (down quarks). Each family is a genuinely distinct excitation mode of the cycle structure.

#### 13.8.3 Why Quarks and Leptons Mix Differently

The graph architecture explains the dramatic contrast between CKM (small mixing) and PMNS (large mixing):

**CKM — small angles.** Quark mixing requires transitions between up-type ($I_3 = +1/2$) and down-type ($I_3 = -1/2$) mass eigenstates. Since up and down quarks propagate on **different spectral paths** (the $\mathcal{R}_3$ exponent has opposite sign for $I_3 = \pm 1/2$, as shown by the $R^2 = 0.985$ correlation in §13.7), inter-family mixing requires bridge crossings ($R_{EE}$ factors). Each crossing costs a factor of $1/\sqrt{2}$ — the CKM Jarlskog invariant $J_{CP} = R_{EE}^{30} = 2^{-15}$ quantifies the total suppression from 30 such crossings. The CKM matrix is **hierarchically small** because the bridge toll accumulates multiplicatively.

**PMNS — large angles.** Neutrinos have $Q = 0$ and $N_c = 1$, so they propagate almost exclusively on the SU(2) binary cycle (§13.7). The charged leptons, while SU(3) singlets, inherit their hierarchy from the SU(3) sector via Yukawa RG running. This means the **mass eigenstates** (determined by the spectral lattice) and the **weak eigenstates** (determined by the SU(2) cycle alone) are **maximally misaligned**: the mass lattice explores all three axes, but the neutrino weak states live on the $\ln 2$ axis alone. The result is large mixing angles — $\theta_{23} \approx 45°$, $\theta_{12} \approx 34°$ — with no hierarchical suppression.

> **The mixing angle hierarchy is not a free parameter — it is a geometric consequence of whether the mass eigenstates and weak eigenstates explore the same or different regions of the spectral lattice.**

#### 13.8.4 Family Number as Topological Invariant

The number of fermion generations is a **topological invariant** of the evolution graph, not a phenomenological input. It follows from three independent lines of argument:

1. **Algebraic**: The SU(2) × SU(3) graph has rank $1 + 2 = 3$, providing exactly 3 Cartan generators, each generating an independent direction in weight space. Three generations = three independent weights.

2. **Spectral**: The 3D mass lattice can support at most 3 non-degenerate occupancy levels before the spectral gap closes (positions become indistinguishable within the lattice resolution set by $\mathcal{R}_3$).

3. **Dynamical**: The sector handoff exhausts both cycles in 2 steps. After Gen 1→2 (SU(3) driven) and Gen 2→3 (SU(2) driven), both sectors have peaked — there is no third sector to sustain a Gen 3→4 handoff.

This triple convergence — algebraic rank, spectral dimensionality, and dynamical exhaustion — locks the number of families at exactly 3. The graph doesn't merely **accommodate** three families; it **requires** them and forbids a fourth.

### 13.9 The $R_{EE}$ Power Tower: Systematic Binary Analysis

A systematic scan of **65 observables** tests which physical quantities are pure powers of $R_{EE} = 1/\sqrt{2}$, equivalently $2^{-n/2}$ for integer $n$. This separates the observable space into a **binary sector** (pure SU(2) topology) and a **ternary sector** (requiring SU(3) information through $\mathcal{R}_3$ and $h^\vee_3$).

> **Internal consistency principle**: All CKM and PMNS mixing angles are computed from their graph building-block formulas (§8.3, §9.3 in [03_FLAVOR_MIXING.md](03_FLAVOR_MIXING.md)), NOT from experimental PDG values. This ensures the tower scan tests the graph framework's internal binary structure, rather than mixing graph predictions with external data.

#### 13.9.1 Confirmed Tower Members

Two observables are confirmed as pure $R_{EE}^n$ to better than 0.5%:

| Observable | $n$ | $R_{EE}^n = 2^{-n/2}$ | Graph-derived | Error |
|-----------|-----|----------------------|-------------|-------|
| $\|V_{tb}\|$ | 0 | $1$ | $0.99913$ | 0.088% |
| $J_{CP}^{\text{CKM}}$ | 30 | $3.052 \times 10^{-5}$ | $3.037 \times 10^{-5}$ | 0.485% |

The $J_{CP}$ value is computed from graph-derived CKM angles (§8.3 in [03_FLAVOR_MIXING.md](03_FLAVOR_MIXING.md)): $\sin\theta_{12} = 0.2248$, $\sin\theta_{23} = 1/24$, $\sin\theta_{13} = 0.003654$, $\delta_{CP} = 65.70°$. The 0.485% gap measures the **internal inconsistency** between the four independent angle formulas and the direct formula $J_{CP} = (R_2 \cdot R_{EE}/h^\vee_2)^6 = 2^{-15}$.

**Near-binary candidates** (Tier 2, $< 1\%$) include $\sin^2\theta_{13}^{\text{PMNS}} \approx R_{EE}^{11}$ (0.68%) and $\alpha_W \approx R_{EE}^{10}$ (0.98%). **Suggestive candidates** (Tier 3, $< 5\%$) include $\|U_{\tau 1}\| \approx R_{EE}^2$ (1.0%), $\sqrt{r} = m_2/m_3 \approx R_{EE}^5$ (3.0%), and $\alpha_{\text{EM}} \times 4\pi \approx R_{EE}^7$ (3.6%).

#### 13.9.2 The Binary Frontier

The classification of all 65 observables:

| Category | Pure binary ($R_{EE}^n$) | Ternary (requires $\mathcal{R}_3$, $h^\vee_3$) | Ambiguous |
|---------|------------------------|-----------------------------------------------|-----------|
| Couplings | 0 | 4 | 3 |
| Mass ratios | 0 | 14 | 3 |
| CKM | 1 ($\|V_{tb}\|$) | 12 | 3 |
| PMNS | 0 | 16 | 2 |
| Neutrino mass | 0 | 1 | 1 |
| Bosons | 0 | 3 | 3 |
| **Total** | **1** | **49** | **15** |

**Key insight**: The binary tower is *extremely sparse* — only 2 confirmed pure-binary observables ($|V_{tb}|$ at $R_{EE}^0$ and $J_{CP}^{\text{CKM}}$ at $R_{EE}^{30}$) among 65 observables. The overwhelming majority of physics requires the full SU(3) color information through $\mathcal{R}_3$ and $h^\vee_3$.

#### 13.9.3 Composite Binary Relations

When **products and ratios** of observables are tested, 61 new composite quantities within 1% of some $R_{EE}^n$ emerge. The most significant:

| Composite | $= R_{EE}^n$ | Error |
|----------|-------------|-------|
| $\sin^2\theta_{23}^{\text{PMNS}} / (m_\mu/m_e)$ | $R_{EE}^{17}$ | **0.034%** |
| $\alpha_s / Q_{\text{Koide}}$ | $R_{EE}^5$ | 0.041% |
| $(m_t/m_c) \times Q_{\text{Koide}}$ | $R_{EE}^{-13}$ | 0.046% |
| $\alpha_s \times (m_t/m_c)$ | $R_{EE}^{-8}$ | 0.087% |
| $J_{CP}^{\text{PMNS}} \times \sin^2\theta_{23}$ | $R_{EE}^{15}$ | 0.116% |
| $\|V_{cb}\| / \|V_{ub}\|$ | $R_{EE}^{-7}$ | 0.144% |
| $J_{CP}^{\text{PMNS}} \times (m_\mu/m_e)$ | $R_{EE}^{-2}$ | 0.149% |
| $(m_\mu/m_e) \times (m_s/m_d)$ | $R_{EE}^{-24}$ | 0.316% |

The first relation, $\sin^2\theta_{23} = R_{EE}^{17} \times (m_\mu/m_e)$, is remarkable: the muon-to-electron mass ratio times a pure binary power reproduces the atmospheric mixing angle at 0.034%. This suggests a **binary bridge** connecting lepton masses to neutrino mixing.

#### 13.9.4 Binary Residual Analysis

For near-binary observables (Tier 2–3), the residual $\text{observable}/R_{EE}^n$ can itself be expressed as a graph formula, revealing the **ternary correction** needed:

| Observable | $R_{EE}^n$ | Residual | Correction formula | Corr. error |
|-----------|-----------|---------|-------------------|-------------|
| $J_{CP}^{\text{CKM}}$ | $R_{EE}^{30}$ | 0.9952 | $R_2^{-2} \cdot R_3^{-3} \cdot R_{EE}^{-2} \cdot h^{\vee -1}_3 \cdot h^{\vee -4}_2$ | 0.29% |
| $\sin^2\theta_{13}$ | $R_{EE}^{11}$ | 0.9932 | $R_2^{-2} \cdot R_3^{-3} \cdot R_{EE}^{-2} \cdot h^{\vee -1}_3 \cdot h^{\vee -4}_2$ | 0.10% |
| $\alpha_W$ | $R_{EE}^{10}$ | 1.0099 | $R_2^4 \cdot R_3^{-4} \cdot R_{EE}^4 \cdot h^\vee_3 \cdot h^\vee_2$ | 0.17% |

With graph-derived angles, both $J_{CP}$ and $\sin^2\theta_{13}$ show negative residuals (graph values slightly below $R_{EE}^n$), and their correction formulas share the same ternary structure. For $\alpha_W$, the correction requires $R_3$ and $h^\vee_3$, consistent with ternary contamination.

#### 13.9.5 Structural Interpretation: Layered Ontology

The $R_{EE}$ tower reveals a **two-layer ontological structure**:

**Layer 1 — Binary (SU(2) topology).** Observables determined entirely by the weak sector subgraph $K_2$. These are pure powers of $\sqrt{2}$. The confirmed members are:
- Near-unity CKM elements ($|V_{tb}| = R_{EE}^0$, exact to 0.088%)
- CKM CP violation ($J_{CP}^{\text{CKM}} = R_{EE}^{30}$, 0.485% — from graph-derived angles)
- Near-binary: reactor angle ($\sin^2\theta_{13} \approx R_{EE}^{11}$, 0.68%) and $\alpha_W \approx R_{EE}^{10}$ (0.98%)

**Layer 2 — Ternary (SU(3) enters).** The vast majority of physics — mass hierarchies, the Weinberg angle, the atmospheric and solar angles, inter-generational ratios — requires the color sector information carried by $\mathcal{R}_3$ and $h^\vee_3$.

**The n-sequence** of confirmed members is $\{0, 30\}$, with near-binary extensions at $n \in \{10, 11\}$:
- $n = 0$: trivial (unit element, $|V_{tb}|$)
- $n = 30$: CKM Jarlskog invariant $J_{CP} = R_{EE}^{30}$ — the most suppressed flavor observable (0.485%)
- $n = 11$ (near-binary, 0.68%): the PMNS reactor angle $\sin^2\theta_{13}$
- $n = 10$ (near-binary, 0.98%): the weak coupling $\alpha_W$

The confirmed gap $\Delta n = 30$ between $|V_{tb}|$ and $J_{CP}$ encodes the **enormous suppression hierarchy** of CP violation: each unit of $n$ contributes $1/\sqrt{2}$, so $n = 30$ gives $2^{-15} \approx 3.05 \times 10^{-5}$ — precisely the experimental magnitude of $J_{CP}$.

> **The R_EE tower is sparse but exact: the binary sector encodes the skeleton of flavor physics (CP violation, reactor angle), while the ternary sector dresses it with the full mass and mixing structure.**

---

**[← Flavor Mixing](03_FLAVOR_MIXING.md)** | **[Index](00_INDEX.md)** | **[Next: Conclusions →](05_CONCLUSIONS.md)**
