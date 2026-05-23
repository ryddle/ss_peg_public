# SS-PEG: Foundations

**[← Index](00_INDEX.md)** | **[Next: Predictions →](02_PREDICTIONS.md)**

---

## 2. The Initiative: What Is SS-PEG?

### 2.1 Core Thesis

The **Space of States and Projections of the Evolution Graph** proposes that physical symmetries — gauge groups, conservation laws, coupling constants — are not fundamental axioms of nature but **emergent consequences** of the relational structure of an underlying evolution graph.

The framework is built on a single ontological chain:

$$\textbf{Potential} \;\to\; \textbf{Relation} \;\to\; \textbf{Cycle} \;\to\; \textbf{Invariant}$$

- **Potentials** are the primary relational substrate (connections on fiber bundles)
- **Relations** connect potentials (covariant derivatives, parallel transport)
- **Cycles** probe global structure (holonomy loops, closed paths)
- **Invariants** are the objective, observer-independent residue (curvature, charges, coupling constants)

### 2.2 The Central Hypothesis

> Given a set of operators organized on a state space, the minimization of the **Killing metric tension** — with no algebraic axioms imposed — produces the unique simple compact Lie algebra compatible with the operator count and state space dimension. All algebraic structure (commutation relations, Jacobi identities, Casimir operators) emerges automatically.

### 2.3 Research Program Structure

The program was executed in five phases, each building on the previous:

| Phase | Focus | Period | Key Output |
|-------|-------|--------|------------|
| **I** | Lie algebra emergence | March 2026 | 51 tests, 100% success |
| **II** | Spectral/Ramanujan analysis | April 2026 | 37 algebras, 11 priorities |
| **III** | Coupling constants | April 2026 | Master formula, 4 couplings |
| **IV** | Mass ratios | April 2026 | 15 mass ratios, all < 1% |
| **V** | Statistical significance | April 2026 | 7.2σ combined significance |
| **VI** | CKM matrix | April 2026 | 23/24 targets < 1%, $J_{CP} = 2^{-15}$ |
| **VII** | PMNS matrix | April 2026 | 28/28 targets < 1% |
| **VIII** | PMNS structural + ν masses | April 2026 | Binary connection, NO predicted, 5 predictions |

### 2.4 Computational Infrastructure

| Component | Specification |
|-----------|--------------|
| GPU | NVIDIA RTX 5070 Ti (16 GB VRAM) |
| CPU | 16 cores |
| RAM | 128 GB DDR5 |
| Software | Python 3.14, PyTorch (CUDA 13) |
| Total optimization runs | ~1,650+ |

---

## 3. Phase I: Lie Algebra Emergence via the Killing Metric

### 3.1 The Experiment

The core experiment minimizes a single objective:

$$T_{\text{total}} = \alpha \, \|K - K_{\text{target}}\|^2 + \beta \, \sum_i (\|G_i\| - 1)^2$$

where $K_{ab} = \text{Tr}(\text{ad}_a \circ \text{ad}_b)$ is the Killing metric of candidate generators $\{G_a\}$, and $K_{\text{target}} = -c \cdot I$ for compact algebras.

Starting from **random** matrices with no algebraic structure, the optimizer converges to exact Lie algebras with:
- Casimir eigenvalues at $10^{-16}$ (machine zero)
- All Jacobi identities satisfied (never imposed)
- 100% algebraic closure

### 3.2 Results Across All Algebra Families

| Family | Algebras Tested | Success Rate | $h^\vee$ Match |
|--------|:-:|:-:|:-:|
| **SU(N)** — Unitary | N = 2, 3, 4, 5 | 100% | Exact |
| **SO(N)** — Orthogonal | N = 3, 4, 5, 6 | 100% | Exact |
| **Sp(2N)** — Symplectic | N = 1, 2, 3 | 100% | Exact |
| **Exceptional** | G₂, F₄, E₆ | 100% | Exact |
| **Non-compact** | sl(2,ℝ), so(3,1), su(2,2) | 100% | Exact |
| **Abelian** | U(1), U(2), U(3) | 100% | N/A |

### 3.3 The Killing-Only Ablation (Test 3.1)

The most important single test: removing the non-commutativity tension $T_{\text{nc}}$ entirely.

- **Killing metric only → SU(2) emerges perfectly.** $T_{\text{nc}}$ is redundant.
- **$T_{\text{nc}}$ only (no Killing) → No algebraic structure.** $T_K$ is essential.

> The Killing metric — which is itself defined through commutators — contains *more* information than the commutators themselves. It encodes the *global geometric pattern* of how all commutators relate.

### 3.4 The Density Phase Transition

For SU(4) with 15 generators, sweeping the number of active generators $n_{\text{gen}}$ from 1 to 15:

- $n_{\text{gen}} = 1$: trivially 100%
- $n_{\text{gen}} = 2$ to $14$: **exactly 0% closure**
- $n_{\text{gen}} = 15$: **exactly 100% closure**

The transition is a mathematical **step function** — there is no "partial symmetry." This is the gauge-theoretic analogue of percolation: either the relational graph is dense enough to support a closed algebraic structure, or it is not.

### 3.5 Spacetime Emergence

The so(3,1) Lorentz algebra emerges from the same principle with indefinite Killing target:
- 100% success rate (5/5 seeds)
- Killing signature **(3+, 3−)**: exactly 3 rotations + 3 boosts
- When compact target conflicts with η-preserving constraint, the optimizer *kills the boost generators*, leaving only so(3) — it knows compact boosts are impossible without being told.

**Implication**: Spacetime is not the arena — it is the pattern of non-compact directions in the emergent Killing form.

### 3.6 Phase I Score

**27/27 testable predictions confirmed**, across 51 independent test programs.

---

## 4. Phase II: Spectral Analysis — Ramanujan Property and Landscape

### 4.1 The Ramanujan Classification

Across 37 simple Lie algebras (all rank ≤ 8), each algebra's adjoint interaction graph was tested for the **Ramanujan property** — optimal spectral expansion (the non-trivial eigenvalues satisfy $|\lambda| \leq 2\sqrt{q-1}$).

**Key finding**: The Standard Model factors SU(2) and SU(3) are **deeply Ramanujan** (ratios 0.50 and 0.55 respectively), passing all 6 alternative definitions tested unanimously.

### 4.2 The Coxeter Phase Transition

A sharp spectral phase transition occurs at Coxeter number $h^* \approx 9.5$:

| Phase | Condition | Properties |
|-------|-----------|------------|
| **Ordered** | $h < h^*$ | Ramanujan, optimal mixing, fast spectral gap |
| **Critical** | $h \sim h^*$ | 18× velocity jump in inner pole migration |
| **Disordered** | $h > h^*$ | Non-Ramanujan, suboptimal expansion |

The SM lives deep inside the **ordered phase**.

### 4.3 The Landscape: Why SU(3) × SU(2) × U(1)?

Monte Carlo landscape analysis with 128 random seeds:

- **0 seeds** → simple exceptional algebra (F₄)
- **53%** → closed but non-simple (product groups)
- **27%** → partially closed
- **20%** → not closed

> The Standard Model's non-simple product structure is the **generic outcome** of the variational principle. Simple exceptional groups have vanishingly small basins of attraction.

### 4.4 Product Composability

Ramanujan × Ramanujan → Ramanujan. The SM product $\mathfrak{su}(3) \oplus \mathfrak{su}(2)$ is unanimously Ramanujan (ratio 0.657).

### 4.5 Landscape Accessibility Prediction

The composite order parameter $\Phi$ predicts experimental basin accessibility with **Spearman $\rho = -0.94$** ($p = 5.5 \times 10^{-4}$):

| Algebra | $\Phi$-predicted | Observed |
|---------|:-:|:-:|
| su(2–5) | >99% | 100% |
| F₄ | ~30% | 6% |
| E₆ | ~10% | 0% (144 trials) |
| E₈ | <1% | landscape-blocked |

---

**[← Index](00_INDEX.md)** | **[Next: Predictions →](02_PREDICTIONS.md)**
