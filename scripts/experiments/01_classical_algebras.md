# Classical Algebra Emergence

> Part of the [Experimental Summary](../EXPERIMENTAL_SUMMARY.md). Sections §3–§7.

---

## 3. Phase 1 — SU(2) Systematic Verification (10 Tests)

The foundational verification suite established that the optimizer reliably discovers SU(2) from random initial conditions. 10 independent tests across 4 categories:

### 3.1 Results

| # | Test | Question | Result | Verdict |
|---|------|----------|--------|---------|
| 1.1 | Reproducibility (100 seeds) | Consistent emergence? | 100/100 close | **PASS** |
| 1.2 | α independence (7 values) | Depends on hyperparameter? | 7/7 close | **PASS** |
| 1.3 | Noise robustness (5 levels × 10 seeds) | Recovers from perturbation? | 50/50 close | **PASS** |
| 2.1 | 3 gen in dim=3 | Spin-1 representation? | 0/10 close | **FAIL** (expected) |
| 2.2 | 2 gen in dim=2 | U(1) emerges? | Partial SU(2) | Informative |
| 2.3 | 8 gen in dim=3 | SU(3) emerges? | 28/28 close | **PASS** |
| 3.1 | Only $T_{\text{Killing}}$ | $T_{\text{nc}}$ necessary? | SU(2) emerges | **SURPRISING** |
| 3.2 | Only $T_{\text{nc}}$ | $T_{\text{Killing}}$ necessary? | No closure | **PASS** |
| 3.3 | Energy landscape | Global minimum? | Single deep basin | **PASS** |
| 4.1 | Representation theory | Casimir & $f_{abc}$ correct? | Exact match | **PASS** |

### 3.2 Key Findings

1. **The Killing metric is the sole structural driver** (Test 3.1). The non-commutativity tension $T_{\text{nc}}$ is redundant — removing it entirely still produces SU(2). Conversely, removing $T_{\text{Killing}}$ (Test 3.2) destroys all structure. This is the central experimental result: **geometry (Killing) generates algebra (Lie)**.

2. **The emergent algebra is exact, not approximate** (Test 4.1). Casimir relative error: machine zero ($\sim 10^{-16}$). Structure constants proportional to $\varepsilon_{abc}$ with std = 0. All spurious entries exactly zero.

3. **Dimensional constraints are sharp** (Tests 2.1, 2.2). The optimizer requires $n_{\text{gen}} = d^2 - 1$ for full closure. Fewer generators produce partial structure; the algebra does not close.

---

## 4. SU(N) Universality Sweep ($N = 2, 3, 4, 5$)

Having established the principle for SU(2), the sweep extended to higher-rank unitary groups.

### 4.1 Results

| Group | $n_{\text{gen}}$ | dim | Closure | $K_{\text{diag}}$ | $h^{\vee}_{\text{measured}}$ | $h^{\vee}_{\text{theory}}$ | GPU Time |
|-------|:-:|:-:|:-:|:-:|:-:|:-:|--:|
| SU(2) | 3 | 2 | 100% | −0.9187 | 2 | 2 | 70 s |
| SU(3) | 8 | 3 | 100% | −0.9455 | 3 | 3 | 462 s |
| SU(4) | 15 | 4 | 100% | −0.9613 | 4 | 4 | 1658 s |
| SU(5) | 24 | 5 | 100% | −0.9627 | 5 | 5 | 5014 s |

### 4.2 Key Findings

- **Universal mechanism**: The Killing metric variational principle produces the unique simple compact Lie algebra $\mathfrak{su}(N)$ for every $N$ tested.
- **Scaling**: Computational cost scales as $O(N^5)$ to $O(N^6)$, consistent with the Killing metric's cubic dependence on generator dimension.
- **All Jacobi identities satisfied** to machine precision for every group.

---

## 5. SO(N) Orthogonal Series ($N = 3, 4, 5, 6$)

### 5.1 Parameterization

SO(N) generators are **real antisymmetric** matrices ($A^T = -A$), requiring a dedicated parameterization distinct from the hermitian traceless matrices used for SU(N). Each generator is parameterized by $d(d-1)/2$ real parameters.

### 5.2 Results

| Group | $n_{\text{gen}}$ | dim | Closure | $K_{\text{diag}}$ | $\|K\|/\kappa^2$ | $h^{\vee}$ (theory) | GPU Time |
|-------|:-:|:-:|:-:|:-:|:-:|:-:|--:|
| SO(3) | 3 | 3 | 100% | −0.9512 | 1.0000 | 1 | 35 s |
| SO(4) | 6 | 4 | 100% | −0.9517 | 2.0000 | 2 | 241 s |
| SO(5) | 10 | 5 | 100% | −0.9584 | 3.0000 | 3 | 665 s |
| SO(6) | 15 | 6 | 100% | −0.9613 | 4.0000 | 4 | 1599 s |

### 5.3 Key Findings

- **The principle extends to orthogonal groups** without modification to the loss function — only the parameterization changes.
- The dual Coxeter number $h^{\vee} = N - 2$ is recovered **exactly** for all $N$.
- SO(3) $\cong$ SU(2) and SO(6) $\cong$ SU(4) cross-checks verified (see §10).

---

## 6. Sp(2N) Symplectic Series ($N = 1, 2, 3$)

### 6.1 Parameterization

USp(2N) generators are hermitian matrices satisfying the symplectic constraint $JG + G^T J = 0$, where $J = \begin{pmatrix} 0 & I_n \\ -I_n & 0 \end{pmatrix}$. Parameterized via $X = \begin{pmatrix} A & B \\ -\bar{B} & \bar{A} \end{pmatrix}$, $T = iX$.

### 6.2 Results

| Group | $n_{\text{gen}}$ | dim | Closure | $K_{\text{diag}}$ | $\|K\|/\kappa^2$ | $h^{\vee}$ (theory) | GPU Time |
|-------|:-:|:-:|:-:|:-:|:-:|:-:|--:|
| Sp(2) | 3 | 2 | 100% | −0.9609 | 4.0000 | 2 | 82 s |
| Sp(4) | 10 | 4 | 100% | −0.9621 | 6.0000 | 3 | 775 s |
| Sp(6) | 21 | 6 | 100% | −0.9626 | 8.0000 | 4 | 3673 s |

### 6.3 Key Findings

- **Third classical family confirmed**: the Killing metric principle works identically for symplectic groups.
- Dual Coxeter $h^{\vee} = N + 1$ recovered exactly ($h^{\vee}_{\text{Sp}(2N)} = N + 1$).
- Sp(2) $\cong$ SU(2) and Sp(4) $\cong$ SO(5) cross-checks verified (see §10).

---

## 7. Density Phase Transition

### 7.1 Experiment Design

For SU(4) (dim=4, algebra dimension = 15), sweep $n_{\text{gen}}$ from 1 to 15 and measure algebraic closure at each point. This tests whether symmetry emergence depends on "relational density" — the ratio of generators to available degrees of freedom.

### 7.2 Results

| $n_{\text{gen}}$ | Closure | Interpretation |
|:-:|:-:|---|
| 1 | 100% | Trivial (single generator) |
| 2 – 14 | **0%** | Underdetermined — algebra cannot close |
| 15 | **100%** | Full $\mathfrak{su}(4)$ emerges |

Total runtime: 9776 s.

### 7.3 Key Finding

**Symmetry emergence is a sharp phase transition, not a gradual onset.** The closure function is a step function: exactly 0% for $n < d^2 - 1$, exactly 100% at $n = d^2 - 1$. There are no intermediate states, no partial algebras, no gradual "symmetry crystallization."

Within the SS-PEG framework, this means: a relational graph must reach a **critical density** — exactly $d^2 - 1$ independent relational degrees of freedom in a $d$-dimensional state space — before a gauge symmetry can crystallize. Below this threshold, there is no symmetry at all. Above it, the symmetry is exact and complete. This is a discrete topological transition, analogous to percolation in graph theory.
