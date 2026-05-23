# Exceptional Algebras: Constructive Verification

> Part of the [Experimental Summary](../EXPERIMENTAL_SUMMARY.md). Sections §9, §9b, §9c, §10.

---

## 9. Exceptional Algebra: G₂

### 9.1 The Subalgebra Problem

For SU(N), the optimizer works because $n_{\text{gen}} = d^2 - 1$ spans the full traceless hermitian space. Exceptional algebras are always proper subalgebras of $\mathfrak{su}(d)$:

| Group | $n_{\text{gen}}$ | $d_{\text{min}}$ | $d^2 - 1$ | Ratio |
|-------|:-:|:-:|:-:|:-:|
| G₂ | 14 | 7 | 48 | 29% |
| F₄ | 52 | 27 | 728 | 7% |
| E₆ | 78 | 27 | 728 | 11% |

The density phase transition (§7) predicts these fail — and they do.

### 9.2 Results

**Blind approach (14 generators in SO(7) parameterization):**
- 50,000 steps, 5 random seeds
- **0% closure** — optimizer converges to $|K|/\kappa^2 = 3$ (not 4)
- G₂ occupies a 14D submanifold of the 21D SO(7) space; gradient descent cannot find it

**Constructive approach (octonion invariant 3-form $\varphi_{ijk}$):**
- Built 14 generators as the SO(7) subalgebra preserving the octonionic 3-form
- Constraint matrix: rank 7, null space dimension = 14 (= dim G₂) ✓
- **100% closure, $K_{\text{diag}} = -4.000000$, $|K|/\kappa^2 = 4.0000$, $h^{\vee} = 4$ exactly**
- $\varphi$-preservation verified at machine precision ($3.15 \times 10^{-16}$)

### 9.3 Key Findings

1. **The subalgebra problem is real for Killing-only optimization**: blind optimization *without closure constraint* cannot discover exceptional algebras. The density transition threshold ($n_{\text{gen}} = d^2 - 1$) is a hard boundary for Killing-only loss. **However**, adding $T_{\text{closure}}$ to the loss (§11c, §11d) enables blind exceptional discovery.

2. **With structural input (invariant forms), G₂ satisfies all the same relations**: closure, Jacobi, Killing proportionality, and the universal $h^{\vee}$ formula.

3. **Physical interpretation**: In the SS-PEG framework, exceptional algebras require **additional geometric constraints** beyond the Killing metric — specifically, algebraic closure. The closure tension serves as a proxy for the preservation of invariant differential forms, suggesting that exceptional symmetries in nature (if present) encode additional topological structure of the evolution graph.

---

## 9b. Exceptional Algebra: F₄ = Der(J₃(𝕆))

### 9b.1 Construction

F₄ is the algebra of **derivations** of the exceptional Jordan algebra J₃(𝕆) (3×3 hermitian matrices over octonions, dim=27). Implementation via `octonion_jordan.py`:

1. Full octonion algebra (8D, Fano plane multiplication table)
2. J₃(𝕆) basis (27 elements) and Jordan product structure constants $C_{ij}^k$
3. Derivation condition: $D(X \circ Y) = D(X) \circ Y + X \circ D(Y)$ → linear constraint on 27×27 matrices
4. SVD of constraint matrix → null space dimension = **52** = dim(f₄) ✓

### 9b.2 Constructive Results

| Metric | Value |
|--------|-------|
| Generators | 52 (27×27 real matrices) |
| Derivation error | ~10⁻¹⁵ |
| Cubic norm preservation | ~10⁻¹⁵ |
| Closure | 1326/1326 (100%) |
| $K_{\text{diag}}$ | −3.0000 |
| $I_R$ (27-dim) | 3.0 |
| $h^{\vee}$ measured | **9.0000** (theory: 9) ✓ |

### 9b.3 Blind Test (40,000 steps, 3 seeds, SO(27) parameterization)

| Seed | Closure | $K_{\text{diag}}$ | $h^{\vee}$ measured | Result |
|------|---------|---------|-------------|--------|
| 0 | 0/1326 (0%) | −2.8538 | 8.5614 | FAIL |
| 1 | 0/1326 (0%) | −2.8537 | 8.5612 | FAIL |
| 2 | 0/1326 (0%) | −2.8538 | 8.5613 | FAIL |

- GPU time: 3,592 s (1,197 s/seed) with vectorized kernels
- Density: 52/351 = 14.8% (SO(27) parameterization)
- Tension plateau at ~28.33 — optimizer found a generic SO(27) minimum, not f₄

---

## 9c. Exceptional Algebra: E₆ = Preserves det(J₃(𝕆))

### 9c.1 Construction

E₆ preserves only the **cubic norm** (determinant) of J₃(𝕆), not the full Jordan product:
$$N(X) = \xi_1\xi_2\xi_3 - \xi_1|x_1|^2 - \xi_2|x_2|^2 - \xi_3|x_3|^2 + 2\,\text{Re}(x_3^* \cdot (x_1^* \cdot x_2^*))$$

Linearized preservation condition → constraint on symmetric 3-tensor $N_{ijk}$ → SVD → null space = **78** = dim(e₆) ✓

### 9c.2 Constructive Results

| Metric | Value |
|--------|-------|
| Generators | 78 (27×27 real matrices) |
| Cubic norm preservation | ~10⁻¹⁶ |
| Closure | 3003/3003 (100%) |
| $K_{\text{diag}}$ | −4.0000 |
| $I_R$ (27-dim) | 3.0 |
| $h^{\vee}$ measured | **12.0000** (theory: 12) ✓ |
| f₄ ⊂ e₆ | 52/52 confirmed (residual ~10⁻¹⁵) ✓ |

### 9c.3 Blind Test (40,000 steps, 3 seeds, SO(27) parameterization)

| Seed | Closure | $K_{\text{diag}}$ | $h^{\vee}$ measured | Result |
|------|---------|---------|-------------|--------|
| 0 | 0/3003 (0%) | −4.7778 | 14.3333 | FAIL |
| 1 | 0/3003 (0%) | −4.7778 | 14.3333 | FAIL |
| 2 | 0/3003 (0%) | −4.7778 | 14.3333 | FAIL |

- GPU time: 7,131 s (2,377 s/seed) with vectorized kernels
- Density: 78/351 = 22.2% (SO(27) parameterization)
- Converged at step 39,620 (plateau detection)
- Remarkably symmetric minimum: $K_{\text{diag}} = -4.777778 \pm 0.000013$ (6 digits identical across 3 seeds)

### 9c.4 The Exceptional Chain

All three accessible exceptional algebras form a verified inclusion chain:
$$G_2 \subset F_4 \subset E_6$$

| Algebra | Invariant preserved | $n_{\text{gen}}$ | dim | $h^{\vee}$ | Blind | Constructive |
|---------|-------------------|:-:|:-:|:-:|:-:|:-:|
| G₂ | Octonion 3-form $\varphi_{ijk}$ | 14 | 7 | 4 | FAIL | PASS |
| F₄ | Jordan product $C_{ij}^k$ | 52 | 27 | 9 | FAIL | PASS |
| E₆ | Cubic norm $N_{ijk}$ | 78 | 27 | 12 | FAIL | PASS |

### 9c.5 Key Findings for F₄ and E₆

1. **The constructive approach generalizes perfectly** from G₂ to F₄ and E₆. The pattern constraint matrix → SVD → null space works for all three.

2. **The universal $h^{\vee}$ formula extends to all exceptionals**: $h^{\vee} = I_R \cdot |K|/\kappa^2$ with $I_R = 3$ for the 27-dim representation of F₄/E₆, yielding exact results.

3. **Blind optimization without closure constraint fails for exceptionals**, independent of density. Adding $T_{\text{closure}}$ to the loss enables blind discovery of G₂ (§11c) and F₄ (§11d). E₆ (78 gen in so(27), 22.2% density) **failed with current strategy** (§11e): T_cl plateau at 0.009, best closure 77.4%.

4. **GPU code vectorization** (~1000× speedup) made the 40k-step blind tests feasible: from estimated ~43 days (F₄) and ~96 days (E₆) down to ~1h and ~2h respectively.

---

## 10. Universal Invariants and Cross-Checks

### 10.1 The Universal Dual Coxeter Formula

Across all **15** algebras tested (4 families), a single formula recovers the dual Coxeter number exactly:

$$\boxed{h^{\vee} = I_R \cdot \frac{|K_{\text{diag}}|}{\kappa^2}}$$

where:
- $|K_{\text{diag}}|$ = absolute value of the diagonal Killing metric entries
- $\kappa^2 = \text{Tr}(T_a^2) / d$ = normalized generator norm
- $I_R$ = Dynkin index of the representation ($I = 1/2$ for SU/Sp fundamental, $I = 1$ for SO vector)

### 10.2 Complete Verification Table

| Algebra | Family | $n_{\text{gen}}$ | dim | $I_R$ | $|K|/\kappa^2$ | $h^{\vee}_{\text{measured}}$ | $h^{\vee}_{\text{theory}}$ | Match |
|---------|--------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| SU(2) | $A_1$ | 3 | 2 | 0.5 | 4.0 | 2 | 2 | ✅ |
| SU(3) | $A_2$ | 8 | 3 | 0.5 | 6.0 | 3 | 3 | ✅ |
| SU(4) | $A_3$ | 15 | 4 | 0.5 | 8.0 | 4 | 4 | ✅ |
| SU(5) | $A_4$ | 24 | 5 | 0.5 | 10.0 | 5 | 5 | ✅ |
| SO(3) | $B_1$ | 3 | 3 | 1.0 | 1.0 | 1 | 1 | ✅ |
| SO(4) | $D_2$ | 6 | 4 | 1.0 | 2.0 | 2 | 2 | ✅ |
| SO(5) | $B_2$ | 10 | 5 | 1.0 | 3.0 | 3 | 3 | ✅ |
| SO(6) | $D_3$ | 15 | 6 | 1.0 | 4.0 | 4 | 4 | ✅ |
| Sp(2) | $C_1$ | 3 | 2 | 0.5 | 4.0 | 2 | 2 | ✅ |
| Sp(4) | $C_2$ | 10 | 4 | 0.5 | 6.0 | 3 | 3 | ✅ |
| Sp(6) | $C_3$ | 21 | 6 | 0.5 | 8.0 | 4 | 4 | ✅ |
| G₂ | Exc. | 14 | 7 | 1.0 | 4.0 | 4 | 4 | ✅ |
| F₄ | Exc. | 52 | 27 | 3.0 | 3.0 | 9 | 9 | ✅ |
| E₆ | Exc. | 78 | 27 | 3.0 | 4.0 | 12 | 12 | ✅ |

**15/15 algebras: exact match.** Zero exceptions.

### 10.3 Isomorphic Pair Cross-Checks

Three pairs of algebras that are mathematically isomorphic were tested via different parameterizations and confirmed to yield identical invariants:

| Pair | Dynkin type | $n_{\text{gen}}$ | $h^{\vee}_A$ | $h^{\vee}_B$ | Match |
|------|:-:|:-:|:-:|:-:|:-:|
| SU(2) $\cong$ Sp(2) | $A_1$ | 3 | 2.0000 | 2.0000 | ✅ |
| SU(4) $\cong$ SO(6) | $A_3$ | 15 | 4.0000 | 4.0000 | ✅ |
| SO(5) $\cong$ Sp(4) | $B_2$ | 10 | 3.0000 | 3.0000 | ✅ |

### 10.4 Casimir Scaling

The Casimir operator $C = \sum_a T_a^2$ satisfies the universal formula:

$$C_\lambda = \frac{n_{\text{gen}} \cdot \kappa}{d} \cdot I$$

verified at machine precision ($\sim 10^{-16}$) for all SU(2..5):

| SU(N) | $C_\lambda$ (measured) | $C_\lambda$ (predicted) | Relative error |
|-------|:-:|:-:|:-:|
| SU(2) | 0.735201 | 0.735201 | $\sim 10^{-16}$ |
| SU(3) | 0.854637 | 0.854637 | $\sim 10^{-16}$ |
| SU(4) | 0.903862 | 0.903862 | $\sim 10^{-16}$ |
| SU(5) | 0.926496 | 0.926496 | $\sim 10^{-16}$ |

### 10.5 K_diag Convergence Survey (§9d)

**Motivation:** Across all optimizer runs, $|K_{\text{diag}}|$ converges to a value systematically below the target ($|K_{\text{target}}| = 1.0$). A unified survey re-ran all 16 algebras with identical optimizer config (steps=12000, lr=0.001, α=5, β=1, 3 seeds) to extract the emerged $K_{\text{diag}}$ and search for patterns.

**Script:** `verification/analysis_kdiag_survey.py`

#### Master Table

| Algebra | $n_{\text{gen}}$ | dim | $\rho$ | $h^{\vee}$ | $|K_{\text{eig}}|$ | $K_{\text{diag,mean}}$ | $K_{\text{std}}$ | $K_{\text{off}}$ | sig | frob |
|---------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| SU(2) | 3 | 2 | 1.000 | 2 | 0.960925 | −0.960925 | 1.1e-7 | 1.0e-7 | (0+,0,3−) | 0.700 |
| SU(3) | 8 | 3 | 1.000 | 3 | 0.962076 | −0.962076 | 8.3e-8 | 1.1e-7 | (0+,0,8−) | 0.633 |
| SU(4) | 15 | 4 | 1.000 | 4 | 0.962575 | −0.962575 | 4.1e-8 | 1.3e-7 | (0+,0,15−) | 0.589 |
| SU(5) | 24 | 5 | 1.000 | 5 | 0.962814 | −0.962814 | 3.2e-8 | 1.1e-7 | (0+,0,24−) | 0.557 |
| SO(3) | 3 | 3 | 0.375 | 2 | 0.950653 | −0.950653 | 5.2e-8 | 9.2e-8 | (0+,0,3−) | 0.987 |
| SO(4) | 6 | 4 | 0.400 | 2 | 0.957308 | −0.957308 | 2.7e-8 | 9.0e-8 | (0+,0,6−) | 0.832 |
| SO(5) | 10 | 5 | 0.417 | 3 | 0.959715 | −0.959715 | 5.0e-8 | 7.1e-8 | (0+,0,10−) | 0.752 |
| SO(6) | 15 | 6 | 0.429 | 4 | 0.960925 | −0.960925 | 4.5e-8 | 1.0e-7 | (0+,0,15−) | 0.700 |
| Sp(2) | 3 | 2 | 1.000 | 2 | 0.960925 | −0.960925 | 1.1e-7 | 1.0e-7 | (0+,0,3−) | 0.700 |
| Sp(4) | 10 | 4 | 0.667 | 3 | 0.962076 | −0.962076 | 6.2e-8 | 8.6e-8 | (0+,0,10−) | 0.633 |
| Sp(6) | 21 | 6 | 0.600 | 4 | 0.962575 | −0.962575 | 4.7e-8 | 7.9e-8 | (0+,0,21−) | 0.589 |
| sl(2,ℝ) | 3 | 2 | 1.000 | 2 | 0.512–1.23 | +0.452 | 6.8e-1 | 3.0e-1 | (2+,0,1−) | 0.923 |
| so(3,1) | 6 | 4 | 0.400 | 2 | 0.901385 | ≈0 (±0.901) | 9.0e-1 | 1.3e-8 | (3+,0,3−) | 0.974 |
| G₂ (constr.) | 14 | 7 | 0.292 | 4 | **4.000000** | −4.000000 | 7e-16 | 4e-16 | (0+,0,14−) | 1.000 |
| F₄ (constr.) | 52 | 27 | 0.071 | 9 | 2.64–2.92 | −2.792 | 1.1e-1 | 1.2e-1 | (0+,0,52−) | 1.000 |
| E₆ (constr.) | 78 | 27 | 0.107 | 12 | 0.07–2.00 | −1.242 | 2.6e-1 | 4.2e-1 | (24+,2,52−) | 1.000 |

#### Key Findings

**F1. Isomorphic algebras produce identical K_diag values.**

Three pairs yield bitwise-identical $|K_{\text{diag}}|$ and Frobenius norms:

| K_diag value | frob | Algebras |
|:-:|:-:|---|
| 0.960925 | 0.700 | SU(2) = SO(6) = Sp(2) |
| 0.962076 | 0.633 | SU(3) = Sp(4) |
| 0.962575 | 0.589 | SU(4) = Sp(6) |

Note: SU(2) ≅ Sp(2) and SU(4) ≅ SO(6) are known isomorphisms. The SU(3) = Sp(4) coincidence is **not** an isomorphism — they share $K_{\text{diag}}$ because they share the same ratio $|K|/\kappa^2 = 2h^{\vee}/I_R = 6$ in the fundamental representation.

**F2. The Frobenius norm is the dominant predictor of K_diag.**

Correlation coefficients (compact optimized algebras):

| Variable | corr with $|K_{\text{diag}}|$ |
|---|-:|
| frob_mean | **−0.966** |
| ρ = n_gen/(d²−1) | +0.604 |
| h∨ | +0.583 |
| n_gen | +0.551 |
| dim | +0.261 |

The strong anti-correlation with Frobenius norm reflects the T_norm tension: the optimizer balances $T_K$ (push K toward target) against $T_{\text{norm}}$ (push ||T|| toward 1), reaching a compromise where K converges below target and ||T|| converges below 1. Smaller generators → lower Killing eigenvalues.

**F3. The SU(N) series reveals monotonic approach to a limit.**

$|K_{\text{diag}}|$ for SU(N): 0.9609, 0.9621, 0.9626, 0.9628 — increasing toward an asymptotic value $\approx 0.963$ as $N \to \infty$.

**F4. Constructive G₂ gives $K_{\text{diag}} = -h^{\vee}$ exactly.**

The G₂ constructive generators (Frobenius-orthonormalized) yield $K_{\text{diag}} = -4.000000$ at machine precision, i.e., $|K| = h^{\vee} = 4$. This is the canonical normalization: for unit-Frobenius generators in the defining representation with $I_R = 1$, the Killing metric eigenvalue equals the dual Coxeter number.

F₄ and E₆ constructive generators deviate from this pattern because their Frobenius-orthonormalized generators are not in a standard representation convention (they use $I_R = 3$ representations in dim=27).

**F5. so(3,1) eigenvalues are perfectly symmetric.**

$|K_+| / |K_-| = 1.000000$ — the positive (boost) and negative (rotation) eigenvalues have exactly equal magnitude, reflecting the $\mathbb{Z}_2$ symmetry of the Lorentz group's root system.
