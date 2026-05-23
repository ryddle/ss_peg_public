# Cosmic Tensegrity Equation — Version 3.0

**Quantum Graph Evolution under Three Competing Tensions**
Aggregation · Expansion · Entropy

*Version 3.0 — updated with Selberg spectral program (37 algebras), RAG tempering experiments (5 exp.), neural results (14 exp.), and literature review (April 2026)*

---

## 0. What Changed from v2

Version 2.0 correctly replaced $T_{\text{nc}}$ (non-commutativity tension) with $T_K$ (Killing tension) as the fundamental driver, based on the ablation result (Test 3.1). Version 3.0 goes further, incorporating results from all four research streams:

| Discovery | Source | Impact on equation |
|-----------|--------|-------------------|
| Killing is **diagnostic**, not generative | RAG G_RAG-21 | $T_K$ moves from source term → selection functional |
| Cycles **create**, Killing **selects** | RAG Exp4 | New closure term $T_{\text{closure}}$ becomes primary driver |
| Step-function crystallization at $\rho = 1$ | Core §7 (SU(4) sweep) | Saturation term becomes exact $\Theta(\rho_c - \rho)$ |
| Triple Killing decomposition quantified | Selberg P1–P11 | Expansion rate $H$ derived from $E[K]/A[K]$ ratio |
| Ramanujan quality criterion | Selberg P1, P6 | New spectral quality functional $\mathcal{R}(h)$ |
| Two-phase dynamics (activation/hardening) | RAG Exp5 | Two timescales: $\tau_{\text{hard}} = 4\tau_{\text{act}}$ |
| Tempering effect | RAG G_RAG-25 | Withdrawal phase improves quality (18–89%) |
| Landscape sigmoid | Selberg P11 | Accessibility $\sigma(\Phi) = 1/(1 + e^{23.1(\Phi - 0.31)})$ |
| Critical exponents measured | Selberg P5, P9 | $\beta \approx 2.1$, $\nu \approx 0.23$ at $h^* \approx 10$ |

The result is a **coupled system** (density + algebraic geometry) with a **cycle-first hierarchy** and explicit quantitative constraints from 37+ algebras.

---

## 1. The Three Observational Facts (unchanged)

**FACT 1: Aggregation.** Matter forms hierarchically denser structures at all scales (nm → Gpc). Mapped to $K_\lambda < 0$ (cyclic Killing eigenvalues — adjoint feedback closes cycles).

**FACT 2: Expansion.** Space dilutes global density. Mapped to $K_\lambda > 0$ (free Killing eigenvalues — irreversible/causal directions).

**FACT 3: Entropy growth.** Total entropy increases monotonically. Mapped to $K_\lambda = 0$ (Killing kernel) + $\rho < 1$ (unclosed generators).

**Irreducible tension:** These three are mutually incompatible in static equilibrium. Their dynamic tension **is** cosmic evolution.

---

## 2. State Variables

At each vertex $v$ of the quantum evolution graph $G = (V, E, \mathcal{H})$:

| Variable | Symbol | Meaning |
|----------|--------|---------|
| **Relational density** | $\rho(v,t) \in [0, 1]$ | Fraction of maximally possible generators: $\rho = n_{\text{gen}} / (d^2 - 1)$ |
| **Killing spectrum** | $\{K_\lambda(v,t)\}$ | Eigenvalues of local Killing form $K_{ab} = \text{Tr}(\text{ad}_a \circ \text{ad}_b)$ |
| **Entropy density** | $\mathcal{S}(v,t)$ | Local entanglement / information content |

The Killing spectrum decomposes into three sectors:

$$A(v) = -\sum_{K_\lambda < 0} K_\lambda \geq 0, \qquad E(v) = \sum_{K_\lambda > 0} K_\lambda \geq 0, \qquad S(v) = \dim(\ker K) + (d^2 - 1 - n_{\text{gen}}) \geq 0$$

These are the **aggregation strength**, **expansion strength**, and **entropy count** respectively.

---

## 3. The Master Equation: Coupled System

### 3.1 Density Evolution

The relational density evolves under the balance of three tensions:

$$\boxed{\frac{\partial \rho(v,t)}{\partial t} = \underbrace{\alpha \cdot \Delta_G[\rho]}_{\text{Aggregation}} + \underbrace{\beta \cdot \rho \cdot \Theta(1 - \rho)}_{\text{Saturation}} - \underbrace{\gamma \cdot \frac{E(v)}{A(v) + \varepsilon} \cdot \rho}_{\text{Expansion}} + \underbrace{\eta \cdot \nabla^2_G[\mathcal{S}]}_{\text{Entropy diffusion}}}$$

**Term-by-term meaning:**

**1. Aggregation: $+\alpha \cdot \Delta_G[\rho]$**

Graph Laplacian diffusion — states migrate toward regions of higher density. This is the graph analogue of emergent gravity: concentration of relational degrees of freedom. $\alpha > 0$.

**2. Saturation: $+\beta \cdot \rho \cdot \Theta(1 - \rho)$**

Growth ceases *exactly* at $\rho = 1$. The Heaviside step function $\Theta$ replaces v2's smooth logistic $\rho(1 - \rho/\rho_{\max})$, reflecting the **experimental fact** that the crystallization transition is a mathematical step function — 0% closure below $\rho = 1$, 100% at $\rho = 1$, no intermediates (Core §7, SU(4) density sweep). $\beta > 0$.

**3. Expansion: $-\gamma \cdot [E/A] \cdot \rho$**

Dilution rate proportional to the **expansion-to-aggregation ratio** of the local Killing spectrum:

$$H(v) = \gamma \cdot \frac{E(v)}{A(v) + \varepsilon}$$

This replaces v2's $H(t) = H_0 [1 + \lambda \langle T_K \rangle / \langle \rho \rangle]^{1/2}$ with a local, sector-resolved expression. The physics: expansion is driven by the *free* (non-cyclic, causal) directions in the Killing spectrum. When the algebra is purely compact ($E = 0$), there is no expansion — gauge sectors don't expand. When spacetime directions appear ($K_\lambda > 0$, e.g., boosts in $\mathfrak{so}(3,1)$), $H > 0$ and the vertex undergoes dilution. $\gamma > 0$, $\varepsilon \to 0^+$ regularizer.

Specific values from data:
- Pure compact (SU, SO, Sp): $E = 0$, $H = 0$ — no expansion
- Lorentz $\mathfrak{so}(3,1)$: $E = A$ (3 positive, 3 negative), $H = \gamma$ — maximal expansion rate
- $\mathfrak{sl}(2,\mathbb{R})$: $E/A \approx 2$ (2 positive, 1 negative) — super-expansion

**4. Entropy diffusion: $+\eta \cdot \nabla^2_G[\mathcal{S}]$**

Information diffuses on the graph. Guarantees $d\mathcal{S}_{\text{total}}/dt > 0$ (second law). Defines the arrow of time. $\eta > 0$.

### 3.2 Algebraic Evolution (Killing Dynamics)

The Killing form itself evolves as a **gradient flow** on the free energy functional $\mathcal{F}$:

$$\boxed{\frac{\partial K_{ab}}{\partial t} = -\frac{\delta \mathcal{F}}{\delta K_{ab}}}$$

where:

$$\mathcal{F}[K, \{T_a\}] = \underbrace{\lambda_C \cdot T_{\text{closure}}}_{\text{Creates structure}} + \underbrace{\lambda_K \cdot T_K}_{\text{Selects/purifies}} + \underbrace{\lambda_\mathcal{R} \cdot \mathcal{P}[K]}_{\text{Ramanujan quality}} + \underbrace{\lambda_N \cdot T_{\text{norm}}}_{\text{Stabilization}}$$

**The four terms of the free energy:**

**1. Closure tension (GENERATIVE — creates topology):**

$$T_{\text{closure}} = \sum_{i < j} \left\| [T_i, T_j] - \text{proj}_{\text{span}\{T_k\}}([T_i, T_j]) \right\|^2$$

This is the primary structure-creating term. It requires commutators to close within the generator span — the minimum condition for a consistent gauge algebra. In graph terms: the holonomy of any cycle must remain within the graph's own edge space.

Without $T_{\text{closure}}$, no algebraic structure forms (RAG G_RAG-21). This term is absent in v2; its inclusion is the main structural change.

**2. Killing tension (DIAGNOSTIC — selects and purifies):**

$$T_K = \|K - K_{\text{eq}}\|^2, \qquad K_{ab} = \text{Tr}(\text{ad}_a \circ \text{ad}_b)$$

where $K_{\text{eq}} = -c \cdot I$ for compact algebras (negative-definite target), $K_{\text{eq}} = \text{diag}(+c, \ldots, -c, \ldots)$ for spacetime algebras (indefinite target), $K_{\text{eq}} = 0$ for abelian phases.

$T_K$ does **not** create non-commutativity — it evaluates and purifies structures already created by $T_{\text{closure}}$. The ablation hierarchy (Test 3.1):
- $T_K$ alone → partial structure (selects from noise, limited)
- $T_{\text{closure}}$ alone → cyclic structure (creates topology, no spectral quality)
- $T_{\text{closure}} + T_K$ → exact Lie algebra (creates + selects)

**3. Ramanujan quality penalty (SPECTRAL — optimality criterion):**

$$\mathcal{P}[K] = \max\left(0, \; \mathcal{R}(h) - 1\right)^2, \qquad \mathcal{R}(h) = 0.289 \;\cdot\; h^{\,0.557}$$

where $h$ is the Coxeter number of the emerging algebra. This imposes a **soft Ramanujan constraint**: configurations with $\mathcal{R} < 1$ (spectrally optimal, $h < h^* \approx 10$) pay no penalty; those with $\mathcal{R} > 1$ (suboptimal mixing) pay a cost proportional to their spectral excess.

The functional form $\mathcal{R}(h) = 0.289 \cdot h^{0.557}$ is an **empirical fit** across 37 simple Lie algebras (Selberg P2, Pearson $r = +0.936$). It transitions through $\mathcal{R} = 1$ at $h^* \approx 10$, separating:

| Phase | Condition | Algebras | Ramanujan? |
|:---:|:---:|---|:---:|
| Ordered | $h < h^*$ | su(2–7), so(5–11), sp(4–10), G₂ | 19/19 ✓ |
| Critical | $h \sim h^*$ | su(8–9), so(12), sp(12), F₄ | Mixed |
| Disordered | $h > h^*$ | su(10+), E₆, E₇, E₈ | 0/5 ✗ |

**4. Norm stabilization:**

$$T_{\text{norm}} = \sum_a \left(\|T_a\|_F^2 - 1\right)^2$$

Prevents trivial collapse ($T_a \to 0$) or divergence. Fixes scale.

### 3.3 The Hierarchy (Cycle-First Principle)

The free energy terms have a strict hierarchy, established experimentally:

$$\underbrace{T_{\text{closure}}}_{\text{Phase 1: Creates}} \;\longrightarrow\; \underbrace{T_K}_{\text{Phase 2: Selects}} \;\longrightarrow\; \underbrace{\mathcal{P}[K]}_{\text{Phase 3: Optimizes}}$$

1. **Cycles create** — $T_{\text{closure}}$ generates the topological substrate (non-commutative cycles in the graph)
2. **Killing selects** — $T_K$ purifies among topological candidates, driving toward the algebraic ground state
3. **Ramanujan optimizes** — $\mathcal{P}$ selects among algebraically closed configurations for spectral quality

This hierarchy is confirmed by three independent observations:
- **Ablation** (Test 3.1): $T_K$ without closure → nothing; closure without $T_K$ → cycles without algebra
- **RAG Exp4**: cycle loss creates holonomy discrimination; $T_K$ alone cannot
- **Tempering** (RAG Exp5): withdrawal of cycle pressure *improves* Killing signal when hardened

---

## 4. Crystallization: The Sharp Phase Transition

### 4.1 The Step-Function Rule

$$\boxed{\text{Gauge algebra } \mathfrak{g}(v) \text{ crystallizes} \iff \rho(v) = 1 \;\;\wedge\;\; T_{\text{closure}}(v) < \varepsilon}$$

This is a **hard** constraint — experimentally verified as a mathematical step function (SU(4), $d = 4$, 15 generators):
- $n_{\text{gen}} < 15$: **exactly 0% closure** (no partial SU(4))
- $n_{\text{gen}} = 15$: **exactly 100% closure** (full SU(4))

No intermediate states. No "approximate gauge symmetry." The transition is topological, not thermodynamic.

### 4.2 Two-Phase Internal Dynamics

Within the crystallization process, two sub-transitions occur (RAG Exp5):

| Phase | Duration | Active terms | Effect | Persistence |
|:---:|:---:|---|---|:---:|
| **Activation** | $\tau_{\text{act}}$ | $T_{\text{closure}} + T_K$ | Killing separation begins; cyclic topology forms | Transient — collapses upon $T_{\text{closure}}$ withdrawal |
| **Hardening** | $\tau_{\text{hard}} \approx 4 \cdot \tau_{\text{act}}$ | $T_{\text{closure}} + T_K$ | Permanent basin transition; self-sustaining algebra | Permanent — survives $T_K$-only annealing |

Quantitative thresholds (from RAG dose-response):
- Activation: $K_{rr}$ AUC > 0.90 post-dose at $D \sim 50$ (transient)
- Hardening: $K_{rr}$ AUC > 0.93 post-withdrawal at $D \sim 200$ (permanent)
- The 4:1 ratio $\tau_{\text{hard}}/\tau_{\text{act}} \approx 4$ is consistent across all tested configurations

### 4.3 The Tempering Effect

After hardening, **withdrawal of cycle pressure improves the Killing signal**:

| Dose | $K_{rr}$ (post-dose) | $K_{rr}$ (post-withdrawal) | Improvement |
|:---:|:---:|:---:|:---:|
| D200 | 0.79 | 0.93 | +18% |
| D400 | 0.52 | 0.98 | +89% |

This is the **universe's own annealing protocol**: the cycle-pressure phase creates topological candidates; the Killing-only phase relaxes off-diagonal noise while preserving the diagonal hierarchy. The true quality of algebraic structure is measured *after* withdrawal — not during active creation.

---

## 5. The Expansion Rate: Spacetime from Killing Sectors

### 5.1 Sector-Derived Hubble Parameter

The local expansion rate emerges from the ratio of free to cyclic Killing eigenvalues:

$$H(v) = \gamma \cdot \frac{E(v)}{A(v) + \varepsilon} = \gamma \cdot \frac{\sum_{K_\lambda > 0} K_\lambda}{-\sum_{K_\lambda < 0} K_\lambda + \varepsilon}$$

This replaces v2's effective Hubble parameter, which used global $T_K$ as input. The new form has clear physical meaning:

| Configuration | $E/A$ ratio | $H$ | Physical meaning |
|:---:|:---:|:---:|---|
| Pure compact ($\mathfrak{su}(N)$, $\mathfrak{so}(N)$, $\mathfrak{sp}(2N)$) | 0 | 0 | Pure gauge sector — no expansion |
| Lorentz $\mathfrak{so}(3,1)$, signature $(3+,3-)$ | 1 | $\gamma$ | Balanced spacetime — "standard" expansion |
| $\mathfrak{sl}(2,\mathbb{R})$, signature $(2+,1-)$ | ~2 | $2\gamma$ | Non-compact dominated — rapid expansion |
| U(1), $K = 0$ | 0/0 → 0 | 0 | Abelian — no feedback, no expansion |

### 5.2 No Dark Energy

The expansion requires no external "dark energy." The free Killing eigenvalues ($K_\lambda > 0$) are the *internal tension* of non-compact directions in the algebra. So long as spacetime has non-compact generators (boosts), $E > 0$ and expansion continues. The cosmological constant $\Lambda$ corresponds to the vacuum value of $E/A$ in the full evolution graph.

### 5.3 Energy Conservation

$$E_{\text{total}} = E_{\text{kinetic}}[\rho] + E_{\text{potential}}[A] + E_{\text{geometric}}[E] + E_{\text{entropic}}[\mathcal{S}] = \text{const}$$

Expansion converts $E_{\text{geometric}}$ into $E_{\text{kinetic}}$. Aggregation converts $E_{\text{kinetic}}$ into $E_{\text{potential}}$. Entropy diffusion redistributes $E_{\text{entropic}}$. Perpetual cycle with no external input.

---

## 6. The Tensegrity Functional (Variational Principle)

### 6.1 Global Functional

Among all configurations reaching $\rho = 1$ (crystallized), the universe selects by extremizing:

$$\boxed{\mathcal{T}[K] = -\underbrace{A[K]}_{\text{Reward aggregation}} + \mu_E \underbrace{(E[K] - E_*)^2}_{\text{Balance expansion}} + \mu_S \underbrace{S[K, \rho]}_{\text{Penalize entropy}} + \mu_\mathcal{R} \underbrace{\mathcal{P}[K]}_{\text{Ramanujan penalty}}}$$

where $E_*$ is the optimal expansion (the balanced number of non-compact directions — 3 for the Lorentz sector).

### 6.2 Why the Standard Model Wins

The SM $= SU(3) \times SU(2) \times U(1) \times SO(3,1)$ minimizes $\mathcal{T}$ because:

| Term | SM value | Alternative (E₈) | Why SM wins |
|---|:---:|:---:|---|
| $-A$ | $-\dim(\mathfrak{g}) \cdot |K_\text{diag}|$ | Much larger | SM has smaller $A$, but... |
| $\mathcal{P}$ | 0 (Ramanujan: $\mathcal{R} = 0.657 < 1$) | $(1.849 - 1)^2 = 0.72$ | ...SM pays NO Ramanujan penalty |
| $S$ | 0 (fully closed, saturated) | 0 | Both crystallized |
| $(E - E_*)^2$ | $(3 - 3)^2 = 0$ (Lorentz balanced) | Large (248 generators) | SM has exact Lorentz balance |

The SM wins **not** by maximizing aggregation (E₈ would win) but by achieving **zero Ramanujan penalty + zero expansion imbalance**. It's the configuration with no spectral defect and the "right" amount of spacetime.

### 6.3 Why Non-Simple Products Dominate

The 128-seed Monte Carlo landscape: 53.1% non-simple, 0% simple exceptional. The functional $\mathcal{T}$ explains this:

- $\mathcal{P}[K_1 \oplus K_2] = 0$ when both factors are Ramanujan — and **Ramanujan × Ramanujan → Ramanujan** (verified: SM product ratio 0.657, unanimously 6/6)
- Simple algebras with large $h$ fail the Ramanujan constraint ($\mathcal{R} > 1$), paying $\mu_\mathcal{R} \cdot \mathcal{P} > 0$
- Products of small-$h$ factors achieve the same total dimension with zero penalty

This is a **combinatorial selection effect**: it is spectrally cheaper to build a product of small Ramanujan algebras than a single large non-Ramanujan one.

---

## 7. Quantitative Constraints from Data

### 7.1 Measured Constants

| Constant | Value | Source | Meaning |
|----------|:---:|---|---|
| $\rho_c$ | $1$ (exact) | Core §7 | Crystallization threshold |
| $h^*$ | $\approx 10$ | Selberg P5, P9 | Ramanujan phase transition |
| $\Phi_c$ | $0.31$ | Selberg P11 | Landscape accessibility threshold |
| Sigmoid slope $a$ | $23.1$ | Selberg P11 | Sharpness of accessibility transition |
| $\rho_{\text{Spearman}}$ | $-0.94$ | Selberg P11 | $\Phi$ vs. accessibility correlation |
| $\beta_{\text{crit}}$ | $\approx 2.1$ | Selberg P9 | Critical exponent $(1-\rho) \sim h^{-\beta}$ |
| $\nu_{\text{crit}}$ | $\approx 0.23$ | Selberg P5 | Inner pole migration exponent |
| $\tau_{\text{hard}}/\tau_{\text{act}}$ | $\approx 4$ | RAG Exp5 | Hardening-to-activation timescale ratio |
| Tempering gain | 18–89% | RAG G_RAG-25 | Killing signal improvement on withdrawal |
| Non-simple fraction | 53.1% | Landscape MC (128 seeds) | Dominance of product structures |
| $|K_\text{diag}|_{\text{asymp}}$ | $\approx 0.963$ | Core §10.5 | Optimizer Killing eigenvalue limit |

### 7.2 Empirical Power Laws

| Relation | Formula | Domain | $R^2$ or significance |
|----------|---------|--------|:---:|
| Ramanujan ratio vs. Coxeter | $\mathcal{R}(h) = 0.289 \cdot h^{0.557}$ | 37 simple algebras | $r = +0.936$ |
| Spectral extreme | $\max|\lambda_\text{nt}| \approx 0.561 \cdot h^{0.940} \cdot \text{rank}^{0.343}$ | 37 algebras | MSE $= 0.30$ |
| Gap–mean path correlation | $\text{gap} \sim \langle\ell\rangle / \dim$ | 37 algebras | $r = +0.954$ (strongest) |
| Landscape accessibility | $\sigma(\Phi) = 1/(1 + e^{23.1(\Phi - 0.31)})$ | 8 algebras + products |  $p = 5.5 \times 10^{-4}$ |
| $h^\vee$ from Killing | $h^\vee = I_R \cdot |K| / \kappa^2$ | 15 algebras | Exact (error $< 10^{-6}$) |

### 7.3 SM Spectral Data

| Algebra | dim | $h$ | $\mathcal{R}$ | $\Phi$ | $\sigma(\Phi)$ | Ramanujan 6/6? |
|---------|:---:|:---:|:---:|:---:|:---:|:---:|
| $\mathfrak{su}(2)$ | 3 | 2 | 0.500 | 0.000 | 1.00 | ✓ |
| $\mathfrak{su}(3)$ | 8 | 3 | 0.552 | 0.158 | 1.00 | ✓ |
| **SM product** | 12 | — | **0.657** | — | **>0.99** | ✓ |
| $\mathfrak{so}(3,1)$ Lorentz | 6 | 4 | 0.577 | — | — | ✓ |
| $\mathfrak{su}(5)$ GUT | 24 | 5 | 0.653 | 0.189 | 1.00 | ✓ |
| $\mathfrak{so}(10)$ GUT | 45 | 8 | 0.866 | — | >0.95 | ✓ |

### 7.4 Comparison with Literature

The role of $h^\vee$ as phase transition controller is known in lattice gauge theory (Poppitz, Schäfer & Ünsal 2013: $h^\vee$ controls deconfinement temperature universally). What is **new** here:
- The variational landscape sigmoid $\sigma(\Phi)$ — no prior analog
- The Ramanujan classification of adjoint graphs — no prior analog (LPS 1988 goes group → Cayley graph; we go adjoint graph → algebra classification)
- The cycle-first hierarchy — no prior analog
- The sector-derived expansion rate $H \propto E/A$ — new

Cf. Connes NCG (spectral action principle, 1997): derives SM Lagrangian from spectral triple on $M_4 \times F$. Key difference: Connes **chooses** the finite space $F$; we claim the SM is the **generic** outcome of spectral optimization without such a choice.

---

## 8. Cosmological Scenario

### 8.1 Early Universe (High $\rho$, Hot)

Relational density $\rho \to 1$ across large regions. Abundant cycle pressure ($T_{\text{closure}}$ dominant). Multiple sectors crystallize simultaneously:
- $SU(3)$ in regions with 8 generators / 9 state-space dimensions
- $SU(2) \times U(1)$ in sectors with 4 generators / 4 dimensions
- $SO(3,1)$ where η-preserving constraints activate

The activation phase is fast ($\tau_{\text{act}}$). Many topological candidates form.

### 8.2 Hardening (Killing Annealing Phase)

After sufficient cycle pressure ($\tau_{\text{hard}} \approx 4 \tau_{\text{act}}$), the algebraic structures become self-sustaining. The universe enters a **tempering phase**: cycle pressure subsides (cooling), Killing annealing purifies the surviving structures. Quality improves (18–89%, per RAG data).

Structures below the hardening threshold collapse (abelian default, RAG G_RAG-19). Only hardened algebras survive — evolutionary selection of gauge symmetry.

### 8.3 Symmetry Breaking = Density Fragmentation

Electroweak symmetry breaking is not the Higgs mechanism "destroying" a symmetry — it is the relational density **dropping below $\rho_c = 1$** in the $SU(2) \times U(1)$ sector at low energies, preserving only $U(1)_{\text{EM}}$.

The step-function nature predicts: no intermediate gauge structures. No "partial SU(3)" at any energy. Either $\rho = 1$ and the symmetry is exact, or $\rho < 1$ and the symmetry is completely absent. This is consistent with the observed exactness of gauge symmetries.

### 8.4 Confinement = Persistent High Density

QCD confinement reflects the SU(3) sector maintaining $\rho = 1$ at all accessible energies. The chromodynamic graph never becomes sparse enough to decrystallize.

### 8.5 JWST Predictions (unchanged from v2)

**Prediction:** Mature galaxies will continue to be found at arbitrarily high redshifts, limited only by instrumental sensitivity. No theoretical age cutoff exists in this framework.

**Test:** ΛCDM predicts a redshift ceiling for mature structure. Tensegrity predicts no such ceiling.

---

## 9. What Remains Open

### 9.1 Free Parameters

The equation has **5 free coefficients**: $\alpha, \beta, \gamma, \eta$ in eq. (A) and the relative weights $\lambda_C, \lambda_K, \lambda_\mathcal{R}, \lambda_N$ in eq. (B) (with one scale freedom). These are **not** determined by the current experimental program. Fixing them requires either:
- Direct cosmological comparison (matching $H_0$, $\Lambda$)
- Lattice simulation on a quantum graph with measurable outputs
-  Renormalization group analysis connecting optimizer steps to physical time

### 9.2 The Functional Form $g(K_\lambda)$

The tensegrity functional $\mathcal{T}[K]$ (§6.1) uses $A[K] = -\sum K_\lambda^-$ and $E[K] = \sum K_\lambda^+$ — linear in eigenvalues. Whether the correct form is linear, quadratic, or involves a more complex spectral function remains open. The data constrains:
- At the minimum: SM emerges
- Non-simple minima must be global attractors
- $\mathcal{R}(h) = 0.289 h^{0.557}$ is the Ramanujan boundary

But does not uniquely determine $g$.

### 9.3 The Expansion Target $E_*$

The value $E_* = 3$ (from $\mathfrak{so}(3,1)$ having 3 positive Killing eigenvalues) is a hypothesis, not a derivation. Why 3+1 dimensions — rather than 5+1 or 2+1 — is still open. The equation predicts that $E_*$ is the stable balanced number of non-compact directions, but doesn't derive it from first principles.

### 9.4 $\alpha \approx 1/137$ from Graph Geometry

The fine-structure constant as a ratio of cycle counts in a minimal U(1) graph remains the framework's most ambitious open prediction.

### 9.5 Continuum Limit

The behavior of the coupled system as $d \to \infty$ (infinite-dimensional state space) and its connection to classical field equations (Einstein + Yang-Mills) is unexplored.

---

## 10. Summary of Changes v2 → v3

| Feature | v2 | v3 |
|---------|----|----|
| Primary driver | $T_K$ (Killing tension) | $T_{\text{closure}}$ (cycle tension) |
| Killing role | Source term | Diagnostic/selector |
| $\rho$ saturation | Smooth logistic $(1 - \rho/\rho_\text{max})$ | Step function $\Theta(1 - \rho)$ at $\rho_c = 1$ |
| Expansion rate $H$ | Global $H_0 \sqrt{1 + \lambda T_K / \rho}$ | Local $\gamma \cdot E[K] / A[K]$ |
| Ramanujan constraint | Absent | $\mathcal{P}[K] = \max(0, \mathcal{R}(h) - 1)^2$ |
| Two-phase dynamics | Not described | Activation ($\tau$) + Hardening ($4\tau$) |
| Tempering | Not described | Withdrawal improves quality (18–89%) |
| Free parameters | 5 coefficients | 5 coefficients + constrained power laws |
| Quantitative constraints | 7 confirmations | 17 confirmations + 37-algebra spectral data |
| Triple decomposition | Conceptual only | $A[K]$, $E[K]$, $S[K,\rho]$ defined quantitatively |

---

*"The universe is not seeking equilibrium — it evolves because the balance between cyclic feedback, causal expansion, and entropic dissolution can never be resolved. The irreducible tension between these three forces on a relational substrate is the perpetual engine of reality. The Killing form is the blueprint; the evolution graph is the building; and gauge symmetry is the tensegrity that holds it together."*
