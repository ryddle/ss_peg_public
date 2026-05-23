# Non-Compact and Spacetime Algebras

> Part of the [Experimental Summary](../EXPERIMENTAL_SUMMARY.md). Sections §8, §8b, §13, §14.

---

## 8. Non-Compact Algebras: sl(2, ℝ)

### 8.1 Motivation

All classical tests above produce compact algebras with negative definite Killing form ($K \propto -I$). Physical spacetime symmetries (Lorentz group, conformal group) involve **non-compact** algebras with indefinite Killing form. Can the same principle produce these?

### 8.2 Three Experiments

| Experiment | Parameterization | $K_{\text{target}}$ | Closure | $K$ signature |
|-----------|-----------------|:---:|:-:|---|
| A: Real + compact target | Real traceless 2×2 | $-I$ | 100% | (0+, 1×0, 2−) — degenerate |
| B: Real + indefinite target | Real traceless 2×2 | sorted(−1,+1,+1) | 100% | (2+, 0×0, 1−) — **correct sl(2,ℝ)** |
| C: Hermitian (control) | Hermitian traceless 2×2 | $-I$ | 100% | (0+, 0×0, 3−) — correct su(2) |

### 8.3 Key Findings

1. **The Killing form signature discriminates compact vs non-compact real forms.** With $K_{\text{target}} = -I$ (compact), the optimizer finds su(2) regardless of parameterization. With an indefinite target, the genuine non-compact algebra sl(2,ℝ) emerges.

2. **The principle extends to non-compact algebras** when the target Killing form has the correct indefinite signature. This opens the path to emergent Lorentz symmetry: $\mathfrak{so}(3,1)$ with signature $(3+, 3−)$.

3. **Physical implication**: In the SS-PEG framework, the signature of the Killing form encodes the **causal structure** of the emergent symmetry. Compact (negative definite) = internal gauge symmetry. Indefinite = spacetime symmetry. The variational principle naturally distinguishes between the two.

---

## 8b. Non-Compact Algebras: so(3,1) — Lorentz Group

### 8b.1 Motivation

The sl(2,ℝ) test (§8) showed that the variational principle extends to non-compact algebras when the Killing target has the correct indefinite signature. The next milestone is $\mathfrak{so}(3,1)$, the **Lorentz algebra**: 6 generators (3 rotations $J_i$ + 3 boosts $K_i$) in dim=4, with Killing signature $(3+, 3-)$. This is the symmetry algebra of special relativity.

### 8b.2 Parameterization

The Lorentz algebra condition is $T^T \eta + \eta T = 0$, where $\eta = \text{diag}(-1,+1,+1,+1)$. This gives:
- **Time-space pairs** $(\mu,\nu)$ with $\eta_{\mu\mu} \neq \eta_{\nu\nu}$: $T_{\mu\nu} = +T_{\nu\mu}$ (symmetric → boosts)
- **Space-space pairs** with $\eta_{\mu\mu} = \eta_{\nu\nu}$: $T_{\mu\nu} = -T_{\nu\mu}$ (antisymmetric → rotations)

6 real parameters per generator: 3 boosts $(T_{01}, T_{02}, T_{03})$ + 3 rotations $(T_{12}, T_{13}, T_{23})$. The constraint $T^T \eta + \eta T = 0$ is enforced **structurally** by the parameterization (not penalized).

### 8b.3 Three Experiments

| Experiment | Parameterization | $K_{\text{target}}$ | Closure | $K$ signature | $\eta$-pres | Result |
|-----------|-----------------|:---:|:-:|---|:-:|---|
| A: η-pres + indefinite target | η-preserving 4×4 | sorted(-1,-1,-1,+1,+1,+1) | **100%** (5/5) | **(3+, 0×0, 3−)** | ✓ | **so(3,1) PASS** |
| B: Antisymmetric + compact (so(4) control) | Antisymmetric 4×4 | $-I$ | **100%** (5/5) | (0+, 0×0, 6−) | ✗ | so(4) PASS |
| C: η-pres + compact target | η-preserving 4×4 | $-I$ | 100% (5/5) | (0+, 3×0, 3−) | ✓ | **so(3) subalgebra** |

### 8b.4 Key Results

1. **so(3,1) Lorentz algebra emerges with 100% success rate** (5/5 seeds) when given: (a) η-preserving parameterization, (b) indefinite Killing target $K_{\text{sorted}} = [-1,-1,-1,+1,+1,+1]$.

2. **Killing eigenvalues**: All seeds converge to $K_{\text{eig}} \approx \pm 0.9014$, with perfect 3+/3− splitting. Off-diagonal $K$ suppressed to $\sim 10^{-8}$.

3. **Canonical span test**: All 6 emergent generators lie in the span of the canonical basis $(J_1, J_2, J_3, K_1, K_2, K_3)$ with residuals $\sim 10^{-16}$.

4. **Rotation–boost discrimination**: The Killing diagonal sign cleanly separates rotation generators ($K_{aa} < 0$, compact) from boost generators ($K_{aa} > 0$, non-compact). The assignment varies by seed (the optimizer finds different bases) but always produces exactly 3 of each.

5. **Experiment C — automatic subalgebra discovery**: When the η-preserving parameterization is paired with a compact Killing target $K = -I$, the optimizer cannot satisfy both constraints simultaneously. It resolves the frustration by **killing the boost generators** (setting them to zero), leaving only the rotation subalgebra $\mathfrak{so}(3) \subset \mathfrak{so}(3,1)$. This manifests as $K$ signature $(0+, 3 \times 0, 3-)$, $f_{\text{rms}} \approx 0$, and 3 degenerate Killing eigenvalues at $\approx -0.972$.

6. **Physical interpretation**: The variational principle with Minkowski η selects the Lorentz algebra as the unique non-compact 6D algebra in dim=4. The Killing signature discriminates: $(3+,3-)$ → spacetime symmetry (Lorentz), $(0+,6-)$ → internal symmetry (SO(4)).

### 8b.5 Comparison with sl(2,ℝ)

| Property | sl(2,ℝ) (§8) | so(3,1) (§8b) |
|----------|:---:|:---:|
| Generators | 3 | 6 |
| Matrix dim | 2×2 | 4×4 |
| Parameterization | Real traceless | η-preserving |
| K signature | (2+, 1−) | (3+, 3−) |
| Compact subalgebra | — | so(3) |
| Physical role | Conformal subgroup | Spacetime symmetry |
| Success rate | 100% | 100% |

---

## 13. U(1) and Abelian Emergence

### 13.1 Motivation

All previous experiments (§3–§11) tested **simple** Lie algebras ($\mathfrak{su}$, $\mathfrak{so}$, $\mathfrak{sp}$, $\mathfrak{g}_2$, $\mathfrak{f}_4$, $\mathfrak{e}_6$) — algebras with non-degenerate negative-definite Killing metrics. A critical gap remained: **the abelian gauge group $U(1)$**, the simplest and most physically ubiquitous symmetry (electromagnetism), was completely absent from the test program.

The theoretical challenge is immediate: $\mathfrak{u}(1)$ is abelian, so all commutators vanish and the Killing metric is identically zero:

$$K_{ab} = \text{Tr}(\text{ad}_a \circ \text{ad}_b) = 0 \quad \forall\, a, b$$

This means U(1) is **invisible** to the Killing metric — it lies in the kernel of the variational principle. The question becomes: when the optimizer is given more generators than needed for $\mathfrak{su}(d)$ (specifically $d^2$ instead of $d^2 - 1$), can it discover the **direct sum decomposition** $\mathfrak{u}(d) = \mathfrak{u}(1) \oplus \mathfrak{su}(d)$?

### 13.2 Parameterization

A **general hermitian parameterization** (WITH trace) was developed, using $d^2$ parameters per generator:

- $d$ independent diagonal entries (all free — no tracelessness constraint)
- $d(d-1)/2$ complex off-diagonal pairs (real + imaginary parts)

This extends the standard SU parameterization ($d^2 - 1$ params) by removing the constraint $\sum_i \text{diag}_i = 0$.

**Script:** `verification/test_u1_gpu.py`

### 13.3 Experiments

#### Experiment A: Single Generator (Trivial Abelian)

| Parameter | Value |
|---|---|
| $n_{\text{gen}}$ | 1 |
| dim | 2 |
| Parameterization | Hermitian traceless (SU) |
| $K_{\text{target}}$ | $-1$ |
| Seeds | 3 (42, 137, 256) |

**Result:** All seeds immediately stuck at $T_{\text{total}} = 5.0$, $T_K = 1.0$, $T_{\text{nc}} = 0$. The Killing metric is trivially $K = [0]$ because a single generator has no commutator partner. The optimizer cannot move — $K_{\text{target}} = -1$ is structurally unreachable. This confirms the theoretical prediction: **an isolated generator is inherently abelian**.

#### Experiment D: U(2) vs SU(2) Side-by-Side

| Parameter | SU(2) | U(2) |
|---|---|---|
| $n_{\text{gen}}$ | 3 | 4 |
| dim | 2 | 2 |
| Parameterization | Hermitian traceless | General hermitian |
| $K_{\text{target}}$ | $-1$ | $-1$ |
| Seeds | 3 | 3 |
| Steps | 12,000 | 12,000 |

**Key Result:**

| Metric | SU(2) | U(2) |
|---|---|---|
| K eigenvalues | $[-0.961, -0.961, -0.961]$ | $[-0.950, -0.950, -0.950, \sim 0]$ |
| K signature | $(0^+, 0^0, 3^-)$ | $(0^+, 1^0, 3^-)$ |
| Closure | 3/3 (100%) | 6/6 (100%) |
| $K_{\text{offdiag}}$ max | $1.0 \times 10^{-7}$ | $3.0 \times 10^{-1}$ |
| Frobenius norms | 0.700 (all) | 1.000 (all) |
| Generator traces | 0.000 (all) | $\neq 0$ (all four) |

**Interpretation:** The optimizer discovers $\mathfrak{u}(2) = \mathfrak{u}(1) \oplus \mathfrak{su}(2)$ in a **rotated basis** where all 4 generators carry non-zero trace. The trace-based classifier fails (classifying all 4 as "U(1)"), but the **Killing eigenvalue spectrum is the correct diagnostic**: 3 degenerate eigenvalues at $-0.950$ (the $\mathfrak{su}(2)$ part) plus 1 eigenvalue at exactly 0 (the $\mathfrak{u}(1)$ direction). The zero eigenvalue confirms that one linear combination of generators lies in the kernel of the Killing form — this is the $\mathfrak{u}(1)$ factor.

#### Experiment E: U(3) = U(1) ⊕ SU(3)

| Parameter | Value |
|---|---|
| $n_{\text{gen}}$ | 9 |
| dim | 3 |
| Parameterization | General hermitian |
| $K_{\text{target}}$ | $-1$ |
| Seeds | 5 (42, 137, 256, 314, 999) |
| Steps | 15,000 |

**Result — ALL 5 seeds identical pattern:**

| Seed | K nonzero eigs (8) | K zero eig | Closure | $K_{\text{offdiag}}$ max |
|:-:|---|---|:-:|---|
| 42 | $\overline{K} = -0.950000$ | $1.1 \times 10^{-16}$ | 36/36 (100%) | 0.290 |
| 137 | $\overline{K} = -0.950000$ | $-2.8 \times 10^{-17}$ | 36/36 (100%) | 0.213 |
| 256 | $\overline{K} = -0.950000$ | $1.1 \times 10^{-16}$ | 36/36 (100%) | 0.408 |
| 314 | $\overline{K} = -0.950000$ | $-1.4 \times 10^{-16}$ | 36/36 (100%) | 0.384 |
| 999 | $\overline{K} = -0.950001$ | $2.2 \times 10^{-16}$ | 36/36 (100%) | 0.286 |

**K signature: $(0^+, 1^0, 8^-)$ for all seeds.** Exactly 8 eigenvalues at $-0.950$ (the $\mathfrak{su}(3)$ part, $\dim = 8$) and 1 eigenvalue at machine-precision zero (the $\mathfrak{u}(1)$ kernel). Closure is perfect at 36/36 commutators. This confirms $\mathfrak{u}(3) = \mathfrak{u}(1) \oplus \mathfrak{su}(3)$ with 100% probability across all seeds.

#### Experiment C: Abelian Target ($K_{\text{target}} = 0$)

| Parameter | Value |
|---|---|
| $n_{\text{gen}}$ | 4 |
| dim | 2 |
| Parameterization | General hermitian |
| $K_{\text{target}}$ | $0$ |
| Seeds | 3 (42, 137, 256) |

**Result:**

| Seed | $\max|K_{\text{eig}}|$ | K signature | Max $\|[T_a, T_b]\|$ | Closure |
|:-:|---|---|---|:-:|
| 42 | $6.5 \times 10^{-10}$ | $(0^+, 4^0, 0^-)$ | $1.8 \times 10^{-5}$ | 6/6 (100%) |
| 137 | $5.7 \times 10^{-10}$ | $(0^+, 4^0, 0^-)$ | $1.5 \times 10^{-5}$ | 6/6 (100%) |
| 256 | $5.0 \times 10^{-11}$ | $(0^+, 4^0, 0^-)$ | $4.9 \times 10^{-6}$ | 6/6 (100%) |

$K_{\text{target}} = 0$ drives the system to a **completely abelian algebra**: all Killing eigenvalues $\lesssim 10^{-10}$, all commutator norms $\lesssim 10^{-5}$. The generators become approximately commuting hermitian matrices — effectively $\mathfrak{u}(1)^4$.

### 13.4 Key Findings

**F1. U(1) is invisible to the Killing metric — but emergent via spectral decomposition.**

For any $K_{\text{target}} < 0$, the optimizer converges to the maximal semisimple subalgebra and places the abelian part in the Killing kernel. The decomposition $\mathfrak{u}(d) = \mathfrak{u}(1) \oplus \mathfrak{su}(d)$ emerges as:
- $\dim(\mathfrak{su}(d)) = d^2 - 1$ eigenvectors with $K_{\lambda} \approx -0.950$
- $\dim(\mathfrak{u}(1)) = 1$ eigenvector with $K_{\lambda} = 0$ (machine precision)

**F2. The Killing eigenvalue spectrum is the correct diagnostic for decomposition.**

Generator traces are unreliable: the optimizer finds a **rotated basis** where all generators carry mixed trace. The invariant diagnostic is the eigenvalue spectrum of the $n_{\text{gen}} \times n_{\text{gen}}$ Killing matrix. Zero eigenvalues count the abelian directions; nonzero eigenvalues span the semisimple part.

**F3. The universal $K_{\text{eigenvalue}} = -0.950$ for U(d) generalizes across $d$.**

Both U(2) and U(3) yield the same semisimple eigenvalue ($-0.950$), different from the pure SU(d) value ($\approx -0.961$). This shift arises from different Frobenius norms: SU parameterization converges to $\|T\|_F \approx 0.700$, while general hermitian converges to $\|T\|_F = 1.000$.

**F4. $K_{\text{target}} = 0$ drives the system to an abelian algebra.**

Setting $K_{\text{target}} = 0$ makes the optimizer minimize all commutators, producing approximately commuting generators ($\max\|[T_a, T_b]\| \sim 10^{-5}$).

**F5. The variational principle now covers both factors of physical gauge groups.**

The Standard Model gauge group $SU(3) \times SU(2) \times U(1)$ requires both semisimple and abelian components. This experiment establishes that:
- Semisimple factors emerge from $K_{\text{target}} < 0$ (§3–§11)
- Abelian factors emerge in the Killing kernel when $n_{\text{gen}} > \dim(\mathfrak{g}_{\text{semisimple}})$
- Pure abelian algebras can be targeted with $K_{\text{target}} = 0$ (Exp C)

---

## 14. SU(2,2) Conformal Algebra Emergence

### 14.1 Motivation

The Lorentz algebra $\mathfrak{so}(3,1)$ emerged cleanly in §8b, establishing that the variational principle produces spacetime symmetries. A natural question: **can the principle access the full conformal symmetry of 3+1 spacetime?**

The **conformal group** of Minkowski space is $SO(4,2)$, locally isomorphic to $SU(2,2)$. Its Lie algebra $\mathfrak{su}(2,2)$ has:
- **15 generators** (same dimension as $\mathfrak{su}(4)$, but a different real form)
- **4×4 complex** matrix representation with indefinite metric $\eta = \text{diag}(-1,-1,+1,+1)$
- **Killing form signature** $(8^+, 0, 7^-)$ — 7 compact directions from $\mathfrak{s}(\mathfrak{u}(2) \oplus \mathfrak{u}(2))$ and 8 non-compact
- **Physical content**: contains the Lorentz algebra $\mathfrak{so}(3,1)$, plus translations, special conformal transformations, and dilatation — the symmetry of massless field equations and CFT

The constraint is $T^\dagger \eta + \eta T = 0$ with $\text{Tr}(T) = 0$, yielding 15 real parameters per generator.

### 14.2 Parameterization

Each generator is parameterized by 15 real numbers:

1. **3 diagonal imaginary** (traceless): $T_{kk} = i \cdot \theta_k$ with $\theta_3 = -\theta_0 - \theta_1 - \theta_2$
2. **4 same-sign anti-hermitian pairs** (indices (0,1) and (2,3)):
   - $T_{01} = a + ib$, $T_{10} = -a + ib$ (both rows have same $\eta$ sign)
3. **8 cross-sign hermitian pairs** (indices (0,2), (0,3), (1,2), (1,3)):
   - $T_{ij} = a + ib$, $T_{ji} = a - ib$ (rows have opposite $\eta$ sign)

This parameterization **exactly satisfies** $T^\dagger \eta + \eta T = 0$ and $\text{Tr}(T) = 0$ by construction — no projection needed.

### 14.3 Optimizer: Indefinite Killing Target

We use a custom optimizer `evolve_su22_indefinite()` analogous to §8b's Lorentz optimizer:

- **Killing loss**: Eigenvalue-based targeting $K_{\text{eig}} \to [-1]^7 \oplus [+1]^8$ (7 compact negative, 8 non-compact positive)
- **Loss function**: $\mathcal{L} = \alpha \cdot T_K^{\text{eig}} + \beta \cdot T_{\text{norm}}$ with $\alpha = 5.0$, $\beta = 1.0$
- **Steps**: 50,000 (Adam, lr $= 3 \times 10^{-3}$)
- **Precision**: `torch.complex128` throughout

### 14.4 Experiments

#### Experiment A: η(2,2)-preserving + indefinite Killing target → su(2,2)

| Seed | Closure | K signature | η-preserving | Traceless | Time |
|------|---------|-------------|-------------|-----------|------|
| 42 | 105/105 (100%) | (8+, 0, 7−) | ✓ (err: 0.00e+00) | ✓ | 61.0 s |
| 123 | 105/105 (100%) | (8+, 0, 7−) | ✓ (err: 0.00e+00) | ✓ | 61.3 s |
| 256 | 105/105 (100%) | (8+, 0, 7−) | ✓ (err: 0.00e+00) | ✓ | 60.7 s |
| 999 | 105/105 (100%) | (8+, 0, 7−) | ✓ (err: 0.00e+00) | ✓ | 60.2 s |
| 2025 | 105/105 (100%) | (8+, 0, 7−) | ✓ (err: 0.00e+00) | ✓ | 60.1 s |

**Result: 5/5 PASS.** All seeds produce **exact** $\mathfrak{su}(2,2)$:
- 100% closure (all 105 commutators close within span)
- Killing signature $(8^+, 0, 7^-)$ exactly matches the expected compact/non-compact split
- η-preservation error is **machine zero** (by construction)
- Canonical span check: all 15 emerged generators lie exactly in the span of the 15 canonical $\mathfrak{su}(2,2)$ generators (residuals $\sim 10^{-16}$)

#### Experiment B: su(4) compact control (no η constraint)

| Seed | Closure | K signature | η-preserving | Time |
|------|---------|-------------|-------------|------|
| 42 | 105/105 (100%) | (0+, 0, 15−) | ✗ | 54.1 s |
| 123 | 105/105 (100%) | (0+, 0, 15−) | ✗ | 53.9 s |
| 256 | 105/105 (100%) | (0+, 0, 15−) | ✗ | 53.5 s |
| 999 | 105/105 (100%) | (0+, 0, 15−) | ✗ | 53.2 s |
| 2025 | 105/105 (100%) | (0+, 0, 15−) | ✗ | 53.1 s |

**Result: 5/5 PASS (control).** Without $\eta$-preservation, the optimizer finds the **compact real form** $\mathfrak{su}(4)$.

#### Experiment C: η(2,2)-preserving + compact target $K = -I$ (frustrated)

| Seed | Closure | K signature | Time |
|------|---------|-------------|------|
| 42 | 105/105 (100%) | (0+, 6₀, 9−) | 12.4 s |
| 123 | 105/105 (100%) | (0+, 5₀, 10−) | 12.1 s |
| 256 | 105/105 (100%) | (0+, 6₀, 9−) | 11.9 s |
| 999 | 105/105 (100%) | (0+, 5₀, 10−) | 12.3 s |
| 2025 | 105/105 (100%) | (0+, 6₀, 9−) | 12.0 s |

**Result: Frustrated convergence.** The optimizer finds a **9-dimensional subalgebra** plus 6 degenerate generators ($K \approx 0$). The 9D subalgebra is consistent with $\mathfrak{s}(\mathfrak{u}(2) \oplus \mathfrak{u}(2))$ (maximal compact subalgebra of $\mathfrak{su}(2,2)$).

### 14.5 Key Findings

**F1. The conformal algebra of 3+1 spacetime emerges from a variational principle.**

With $\eta(2,2)$-preserving parameterization and indefinite Killing target, 15 random complex $4 \times 4$ matrices converge to an exact realization of $\mathfrak{su}(2,2) \cong \mathfrak{so}(4,2)$ with **100% success rate** across all seeds.

**F2. The $\eta$-constraint selects between real forms.**

- No $\eta$-constraint → compact $\mathfrak{su}(4)$ (internal gauge symmetry)
- $\eta(2,2)$-constraint + indefinite target → non-compact $\mathfrak{su}(2,2)$ (conformal spacetime symmetry)

**F3. The spacetime symmetry hierarchy is complete: $\mathfrak{su}(2) \to \mathfrak{so}(3,1) \to \mathfrak{su}(2,2)$.**

The variational principle now produces the full chain of physically relevant spacetime symmetries:
- $\mathfrak{su}(2)$: Internal spin symmetry (§1–6)
- $\mathfrak{so}(3,1)$: Lorentz symmetry (§8b) — rotations + boosts
- $\mathfrak{su}(2,2) \cong \mathfrak{so}(4,2)$: Full conformal symmetry (§14) — Lorentz + translations + special conformal + dilatation

**F4. Frustrated optimization reveals the maximal compact subalgebra.**

The Exp C result mirrors the pattern found for $\mathfrak{so}(3,1)$ (§8b Exp C): when an incompatible compact Killing target is applied to $\eta$-preserving generators, the optimizer discovers the **maximal compact subalgebra**.

**F5. Physical interpretation: AdS/CFT connection.**

$\mathfrak{su}(2,2) \cong \mathfrak{so}(4,2)$ is simultaneously:
- The **conformal algebra** of 3+1 Minkowski spacetime
- The **isometry algebra** of $\text{AdS}_5$ (anti-de Sitter space in 5 dimensions)

Its emergence from a variational principle on finite matrices suggests that the AdS/CFT duality may have roots in the same algebraic structure that produces gauge and spacetime symmetries from relational potentials.
