# Literature Review & Novelty Assessment

## SS-PEG Claims vs. Existing Published Work

**Date:** April 2026  
**Purpose:** Systematic comparison of SS-PEG experimental results against published literature to assess novelty, identify overlaps, and locate the framework within the existing body of knowledge.

---

## Overview of Search Methodology

Systematic Google Scholar searches were conducted across 8 thematic axes:
1. Ramanujan graphs + Lie algebra spectral properties
2. Cayley graphs of adjoint representations + expander properties
3. ML / optimization for Lie algebra and symmetry discovery
4. Emergent gauge symmetry from graph/network dynamics
5. Standard Model gauge group uniqueness / selection principles
6. Dual Coxeter number as phase transition controller
7. Connes' Noncommutative Geometry and the spectral action principle
8. Quantum Graphity and related graph-based quantum gravity

The existing workspace bibliography (34 refs in .tex, plus marco_matematico_formal.md) was also inventoried. **Notable gap identified:** the current bibliography contains NO references to Ramanujan graphs, spectral graph theory on Lie algebras, ML for symmetry discovery, or Connes' NCG program.

---

## Claim-by-Claim Assessment

### CLAIM 1: Killing metric alone generates Lie algebras (variational principle)

**SS-PEG result:** Minimizing $T_K = \|K - K_{\text{target}}\|^2$ on random matrices produces the unique simple compact Lie algebra of the correct dimension, with 100% success. Removing $T_K$ and keeping only non-commutativity loss produces nothing.

**Closest prior work:**

| Reference | Cites | What they do | Difference from SS-PEG |
|-----------|:-----:|-------------|----------------------|
| **Dehmamy, Walters, Liu (NeurIPS 2021)** "Automatic symmetry discovery with Lie algebra convolutional network" | 153 | Discovers symmetry transformations of data by optimizing Noether current conservation | Discovers symmetries *of data/equations*, not algebra *from the Killing form*. The optimization target is Noether loss, not Killing tension. |
| **Yang, Walters, Dehmamy (ICML 2023)** "Generative adversarial symmetry discovery" | 71 | Uses GAN loss to discover Lie algebra bases including gauge symmetry | Discovers symmetries of a learned manifold via adversarial training. Algebra emerges from data, not from pure $K$-optimization. |
| **Forestano et al. (Phys. Lett. B, 2023)** "Discovering sparse representations of Lie groups with ML" | 12 | Deep learning for symmetry transformations and algebra generators | Similar spirit but fundamentally different mechanism: learns from data, not from Killing form as sole input. |
| **Karjol et al. (2025)** "Automatic Discovery of One-Parameter Subgroups" | — | Discovers SO(n) and SL(n) subgroups from optimization | Discovers subgroups of known groups, not algebras from scratch. |
| **Shumaylov et al. (2024)** "Lie algebra canonicalization" | 18 | Equivariant neural operators under arbitrary Lie groups | Engineering tool, not discovery tool. |

**Novelty assessment: ★★★★☆ HIGHLY NOVEL**

The specific claim — that the Killing metric *alone*, without data, without Noether currents, without adversarial training, generates Lie algebras from random matrices — has no direct precedent. The ML symmetry discovery literature discovers symmetries *of physical systems/data*; SS-PEG claims the Killing geometry is sufficient *without any external data*. The Dehmamy et al. line of work is the closest conceptual neighbor but operates in a fundamentally different regime: physics data in → symmetry out, vs. Killing target in → algebra out.

**Risk of triviality:** LOW. The ablation result (Test 3.1: $T_K$ essential, $T_{\text{nc}}$ redundant) is a genuinely novel finding with no published counterpart.

---

### CLAIM 2: SM gauge group is "spectrally optimal" (Ramanujan property)

**SS-PEG result:** The adjoint interaction graphs of $\mathfrak{su}(2)$ and $\mathfrak{su}(3)$ satisfy the Ramanujan bound. SM products remain Ramanujan under composition. Algebras are classified by spectral gap ratio relative to the Ramanujan bound.

**Closest prior work:**

| Reference | Cites | What they do | Difference from SS-PEG |
|-----------|:-----:|-------------|----------------------|
| **Lubotzky, Phillips, Sarnak (1988)** "Ramanujan graphs" | 1400+ | Construct $(p+1)$-regular Ramanujan graphs from $PGL(2, \mathbb{Z}/q\mathbb{Z})$ | The LPS construction goes *from* a Lie group quotient *to* a Cayley graph. SS-PEG goes *from* the adjoint representation structure *to* a Ramanujan classification. The graph objects are completely different. |
| **Lubotzky (1994)** *Discrete groups, expanding graphs and invariant measures* | 1068 | Connects property (T) of Lie groups to expansion of Cayley graphs | Studies Cayley graphs of the groups themselves, not adjoint interaction graphs of their algebras. |
| **Tao (2015)** *Expansion in finite simple groups of Lie type* | 78 | Expansion properties of Cayley graphs of $SL_n(\mathbb{F}_p)$ etc. | Cayley graphs of finite groups of Lie type, not spectral classification of adjoint interaction graphs. |
| **Fontana et al. (2023)** "The adjoint is all you need" | 68 | Uses Killing form norm and spectral gap in quantum circuit ansätze to characterize barren plateaus | Uses adjoint + Killing form + spectral gap, but for quantum algorithm design, not for classifying Lie algebras. Closest structural analog. |
| **Dalal et al. (2026)** "Ramanujan Complexes from Unitary Groups" | — | Golden gates for PU(5) from Ramanujan complexes | Constructs quantum gates, not spectral classification of algebras. |

**Novelty assessment: ★★★★★ NOVEL**

No one has applied the Ramanujan graph criterion to *adjoint interaction graphs* of Lie algebras as a classification/selection principle. The entire literature on Ramanujan graphs and Lie groups operates in the *opposite direction* (Lie group → Cayley graph → prove Ramanujan). The SS-PEG approach of asking "is the adjoint-representation interaction graph Ramanujan?" and using this to classify algebras appears to be entirely new.

Fontana et al. (2023) is the closest structural analog (adjoint + Killing + spectral gap), but their purpose is quantum computing, not algebra classification.

**Risk of triviality:** VERY LOW. This is a genuinely novel application of spectral graph theory to Lie algebra classification.

---

### CLAIM 3: Dual Coxeter number $h^\vee$ as the master parameter

**SS-PEG result:** $h^\vee = I_R \cdot |K_{\text{diag}}| / \kappa^2$ universally. $h^\vee$ controls the spectral phase transition at $h^* \approx 9.5$, landscape accessibility via $\Phi(h) \approx 0.289 \cdot h^{0.557}$, and the Ramanujan classification.

**Closest prior work:**

| Reference | Cites | What they do | Difference from SS-PEG |
|-----------|:-----:|-------------|----------------------|
| **Poppitz, Schäfer, Ünsal (2013)** "Universal mechanism of (semi-classical) deconfinement" | 147 | Show $h^\vee$ controls the deconfinement temperature universally for all simple gauge groups in YM theory | $h^\vee$ as universal controller of phase transitions is KNOWN in lattice gauge theory. However, they study thermal deconfinement, not spectral phase transitions in an optimization landscape. |
| **Buisseret (2011)** "Structure of Yang-Mills spectrum for arbitrary simple gauge algebras" | 8 | Phase transitions in YM spectrum governed by $h^\vee$ | Same: $h^\vee$ controls transitions, but in a physical (lattice QCD) context, not an optimization landscape. |
| **Gepner (1987)** CFT associated with Lie algebras | 302 | $h^\vee$ appears in CFT partition functions and central charges | Standard result in CFT: $c = k \cdot \dim(\mathfrak{g}) / (k + h^\vee)$. |
| **Korff (2000)** Lie algebraic structures in integrable models | — | Coxeter / dual Coxeter control phase transitions in integrable systems | Mathematical physics, different context. |

**Novelty assessment: ★★★☆☆ PARTIALLY NOVEL**

⚠️ **The role of $h^\vee$ as a universal organizer of phase transitions is WELL ESTABLISHED in lattice gauge theory (Poppitz et al.) and CFT (Gepner).** This is the area of greatest overlap with existing literature.

However, the *specific* SS-PEG contributions are novel:
- The formula $h^\vee = I_R \cdot |K_{\text{diag}}| / \kappa^2$ from variational optimization (not from theory)
- The spectral phase transition in an *optimization landscape* (not thermal/deconfinement)
- The sigmoid accessibility function $\Phi(h)$ and the existence of $\Phi_c = 0.31$ as a landscape threshold
- The connection of $h^\vee$ to Ramanujan classification

**Risk of triviality:** MODERATE for the general role of $h^\vee$. LOW for the specific landscape results ($\Phi(h)$, $\Phi_c$, sigmoid, spectral phase transition in optimization).

**Recommendation:** The paper MUST cite Poppitz et al. (2013) and Gepner (1987). It should frame the $h^\vee$ results as: "The known role of $h^\vee$ in physical phase transitions (Poppitz et al.) manifests in our variational landscape as..." This prevents the appearance of claiming novelty for the general importance of $h^\vee$.

---

### CLAIM 4: SM gauge group from spectral/geometric principles (uniqueness)

**SS-PEG result:** The SM group $SU(3) \times SU(2) \times U(1)$ emerges as the spectrally optimal configuration — Ramanujan, landscape-dominant, composable.

**Closest prior work:**

| Reference | Cites | What they do | Difference from SS-PEG |
|-----------|:-----:|-------------|----------------------|
| **Chamseddine & Connes (1996, 1997)** "The spectral action principle" | ~984 | Derive the full SM Lagrangian from a spectral action on product geometry $M_4 \times F$ within noncommutative geometry | **THE MAIN COMPETITOR.** Also derives SM from spectral principles. But uses a fundamentally different mechanism: inner automorphisms of an almost-commutative algebra on a spectral triple, not variational optimization of the Killing metric on random matrices. |
| **Connes (2006)** "Noncommutative geometry and the standard model with neutrino mixing" | 400+ | Derives SM with massive neutrinos from NCG | Extended Connes program. Different formalism, same spirit. |
| **Besnard (2019)** NCG + B-L extension | — | Extends Connes program to include B-L symmetry | Extension of Connes, same family. |
| **Axelkrans (2025)** "SM Gauge Symmetries from Projection Invariance" | — | Derives SM gauge group from projection invariance | Different geometric principle, same goal. Very recent. |
| **Livolsi** "Emergent Gauge Fields from Universal $\Psi$ Lagrangian" | — | Derives $SU(3) \times SU(2) \times U(1)$ without fundamental gauge symmetries | Different mechanism (Lagrangian-based, not spectral). |
| **Slansky (1981)** "Group theory for unified model building" | 2071 | Comprehensive review of group selection for unification | Surveys options but does NOT derive SM uniqueness. |

**Novelty assessment: ★★★★☆ NOVEL (but must acknowledge Connes)**

The Connes NCG program (spectral action principle) is the most important comparison point. Both approaches share the philosophy that the SM gauge group is *derived* from spectral/geometric principles rather than postulated. Key differences:

| Aspect | Connes NCG | SS-PEG |
|--------|-----------|--------|
| Starting point | Almost-commutative algebra on $M_4 \times F$ | Random matrices + Killing tension |
| Mechanism | Inner automorphisms of spectral triple | Variational optimization of $T_K$ |
| What's derived | Full SM Lagrangian (gauge + Higgs + fermion masses) | Gauge algebra classification + spectral optimality |
| Spacetime | Input ($M_4$) | Emergent (from indefinite $K$ signature) |
| Finite space $F$ | Input (chosen to give SM) | Not needed — algebra emerges from dimensionality |

**The critical difference:** Connes must *choose* the finite space $F$ (essentially encoding the SM content), and then the spectral action produces the SM Lagrangian. SS-PEG claims the SM is the *generic* outcome of variational optimization (landscape dominance + Ramanujan optimality), requiring no such choice.

**Risk of triviality:** LOW. The mechanism is genuinely different from Connes. But incomplete: SS-PEG derives algebra classification, not the full Lagrangian. Connes derives more.

**Recommendation:** The paper MUST discuss Connes extensively. It should frame SS-PEG as complementary: "While Connes' NCG derives the SM Lagrangian from a chosen spectral triple, our approach derives the *selection* of the gauge algebra as the spectrally optimal outcome of a variational principle, without choosing a finite geometry."

---

### CLAIM 5: Emergent gauge symmetry from graph dynamics / relational density

**SS-PEG result:** Gauge groups crystallize at a sharp density threshold $\rho_c = 1$. No intermediate gauge structures exist (step function). The graph relational density determines which symmetry emerges.

**Closest prior work:**

| Reference | Cites | What they do | Difference from SS-PEG |
|-----------|:-----:|-------------|----------------------|
| **Konopka, Markopoulou, Severini (2008)** "Quantum Graphity" | 263 | Emergent U(1) gauge from graph-state phase transitions | Already in workspace bibliography. Emerges U(1) only, from specific Hamiltonian on graph states. SS-PEG claims universal Lie algebra emergence from Killing optimization. |
| **Van Damme, Halimeh, Hauke (2020)** "Gauge-symmetry violation quantum phase transition" | 26 | Phase transitions in lattice gauge theories via gauge violation | Studies existing gauge theories, not emergence of new ones. |
| **Halimeh et al. (2024)** "Spin exchange-enabled quantum simulator for large-scale non-Abelian gauge theories" | 35 | Simulates non-Abelian gauge theories on spin systems | Simulation, not derivation of gauge emergence. |
| **F. Michael (2026)** "Emergent Geometry, Gauge Structure, and Unique Four-Dimensional Vacuum" | — | Derives gauge fields from discrete state updates on lattice/graph | Very recent. Similar spirit. Need to compare mechanism in detail. |

**Novelty assessment: ★★★★☆ HIGHLY NOVEL**

Quantum Graphity is the closest predecessor, but it achieves only U(1) emergence and uses a very specific Hamiltonian. SS-PEG's claim of *universal* gauge algebra emergence from a single variational principle ($T_K$) with sharp percolation-like phase transition is substantially more general.

The step-function transition with no intermediate states is a strong, testable, and apparently novel prediction.

**Risk of triviality:** LOW. The specific mechanism and universality claim are new.

---

### CLAIM 6: Triple balance / tensegrity (aggregation, expansion, entropy)

**SS-PEG result:** The evolution graph architecture is determined by the equilibrium of three sectors of Killing eigenvalues: aggregation ($K_\lambda < 0$), expansion ($K_\lambda > 0$), and entropy ($K_\lambda = 0$).

**Closest prior work:** No direct analog found in the literature.

**Novelty assessment: ★★★★★ NOVEL**

The decomposition of Killing eigenvalue sectors into competing "tensions" with physical interpretations (aggregation vs. expansion vs. entropy) appears to have no precedent. This is an internal SS-PEG construction with no established counterpart.

**Risk of triviality:** LOW, but also UNTESTED — the tensegrity equation form is listed as "Open" in the predictions table.

---

### CLAIM 7: Killing form is diagnostic not generative (RAG result)

**SS-PEG result:** Cycle loss creates algebraic structure; Killing tension selects/purifies. They play complementary roles: creation vs. selection.

**Closest prior work:** No direct analog in the literature. This is an empirical finding about the optimization process itself.

**Novelty assessment: ★★★★★ NOVEL**

This is a meta-result about the optimization dynamics, not a physics claim per se. No one else has studied this decomposition because no one else has the SS-PEG variational framework.

**Risk of triviality:** LOW.

---

### CLAIM 8: Spacetime signature from Killing form signature

**SS-PEG result:** The indefinite Killing form of $\mathfrak{so}(3,1)$ emerges from η-preserving parameterization + indefinite target. The (3+, 3−) signature maps to Minkowski spacetime structure.

**Closest prior work:**

| Reference | Cites | What they do |
|-----------|:-----:|-------------|
| Classical result | — | The Killing form of $\mathfrak{so}(3,1)$ having signature (3,3) is a standard mathematical fact |
| Standard gauge theory texts | — | The distinction between compact (internal) and non-compact (spacetime) algebras via Killing signature is textbook material |

**Novelty assessment: ★★☆☆☆ PARTIALLY KNOWN**

⚠️ **The connection between Killing form signature and compact/non-compact classification is standard Lie theory.** What is potentially new is the *emergence* of the correct signature from variational optimization (not assumed from the start), and the philosophical interpretation that "spacetime is the pattern of non-compact directions."

**Risk of triviality:** MODERATE for the mathematical fact. LOW for the variational emergence mechanism.

**Recommendation:** Frame carefully — the paper should NOT claim that "Killing signature classifies compact vs. non-compact" as a discovery. It should say: "The *known* relationship between Killing signature and spacetime structure is realized *dynamically* in our variational framework..."

---

## Summary Table

| # | Claim | Novelty | Triviality Risk | Key Comparison |
|:-:|-------|:-------:|:---------------:|---------------|
| 1 | Killing metric alone generates Lie algebras | ★★★★☆ | LOW | Dehmamy et al. (NeurIPS 2021) — different mechanism |
| 2 | SM as Ramanujan (spectrally optimal) | ★★★★★ | VERY LOW | LPS (1988) — opposite direction |
| 3 | $h^\vee$ as master parameter | ★★★☆☆ | MODERATE | **Poppitz et al. (2013)** — $h^\vee$ in phase transitions is KNOWN |
| 4 | SM from spectral principles | ★★★★☆ | LOW | **Connes NCG (1996–2006)** — different mechanism, similar philosophy |
| 5 | Emergent gauge from graph density | ★★★★☆ | LOW | Quantum Graphity (2008) — U(1) only |
| 6 | Triple balance / tensegrity | ★★★★★ | LOW | No precedent found |
| 7 | Killing diagnostic not generative | ★★★★★ | LOW | No precedent found |
| 8 | Spacetime from Killing signature | ★★☆☆☆ | MODERATE | Standard Lie theory classification |

---

## Critical References to Add to Bibliography

### Must-cite (first priority)

1. **Chamseddine, A. H., & Connes, A. (1997).** *The spectral action principle.* Comm. Math. Phys. 186, 731–750. [~984 cites] — THE main competitor approach deriving SM from spectral principles.

2. **Poppitz, E., Schäfer, T., & Ünsal, M. (2013).** *Universal mechanism of (semi-)classical deconfinement and θ dependence for all simple groups.* JHEP 03, 087. [147 cites] — $h^\vee$ as universal phase transition controller in gauge theory.

3. **Dehmamy, N., Walters, R., & Liu, Y.-Y. (2021).** *Automatic symmetry discovery with Lie algebra convolutional network.* NeurIPS. [153 cites] — ML for symmetry discovery.

4. **Lubotzky, A., Phillips, R., & Sarnak, P. (1988).** *Ramanujan graphs.* Combinatorica 8(3), 261–277. [1400+ cites] — Original Ramanujan graph construction from Lie groups.

5. **Lubotzky, A. (1994).** *Discrete groups, expanding graphs and invariant measures.* Birkhäuser. [1068 cites] — Property (T), expansion, and Lie groups.

### Should-cite (second priority)

6. **Gepner, D. (1987).** *On the spectrum of 2D conformal field theories associated with Lie algebras.* Nuclear Physics B. [302 cites] — $h^\vee$ in CFT.

7. **Fontana, E., et al. (2023).** *The adjoint is all you need.* [68 cites] — Killing form + spectral gap in quantum computing.

8. **Yang, J., Walters, R., & Dehmamy, N. (2023).** *Generative adversarial symmetry discovery.* ICML. [71 cites] — GAN-based Lie algebra discovery.

9. **Connes, A. (2006).** *Noncommutative geometry and the standard model with neutrino mixing.* JHEP. [400+ cites] — Extended NCG program.

10. **Tao, T. (2015).** *Expansion in finite simple groups of Lie type.* Graduate Studies in Mathematics, AMS. [78 cites] — Expansion of Cayley graphs.

### Nice-to-cite (context)

11. **Buisseret, F. (2011).** *Structure of the Yang-Mills spectrum for arbitrary simple gauge algebras.* [8 cites]
12. **Forestano, K. T., et al. (2023).** *Discovering sparse representations of Lie groups with ML.* Phys. Lett. B. [12 cites]
13. **Axelkrans, J. (2025).** *Standard Model Gauge Symmetries from Projection Invariance.*
14. **Sawicki, A., & Karnas, K. (2017).** *Universality of single-qudit gates.* [44 cites] — Killing form + spectral gap
15. **Méliot, P.-L. (2019).** *Asymptotic representation theory and spectrum of random geometric graph on compact Lie group.*

---

## Conclusions

### What is NOT novel (known results in different clothing)
- The general importance of $h^\vee$ for phase transitions → Poppitz et al. (2013), Gepner (1987)
- Killing form signature classifying compact vs. non-compact → standard Lie theory
- ML can discover symmetries of data → Dehmamy et al. (2021), Yang et al. (2023)

### What IS genuinely novel
- **Killing form as sole generative principle** (no data, no Noether currents, no algebraic axioms)
- **Ramanujan classification of adjoint interaction graphs** and SM as spectrally optimal
- **Sharp density phase transition** with step-function character (not thermal)
- **Triple balance / tensegrity** decomposition of Killing eigenvalue sectors
- **Killing diagnostic not generative** (creation vs. selection duality in optimization)
- **Landscape sigmoid** $\Phi(h)$ and $\Phi_c = 0.31$ as spectral accessibility threshold
- **Combined framework** deriving gauge algebra selection from variational optimization (vs. Connes' NCG which requires choosing the finite space $F$)

### Overall verdict

**SS-PEG is NOT duplicating known results.** The core mechanism (Killing-only variational principle) and its application to algebra classification via spectral graph theory are genuinely new. The framework sits in a novel intersection:

```
Spectral Graph Theory (Ramanujan/LPS)
    ↕
SS-PEG Framework  ←→  Connes NCG (spectral action)
    ↕
ML Symmetry Discovery (Dehmamy et al.)
```

The main risk of appearing unoriginal comes from **Claim 3** ($h^\vee$ as master parameter) and **Claim 8** (Killing signature), which overlap with well-known results. These must be carefully framed as "the known property X manifests in our variational framework as Y."

The framework is NOT complete (does not derive the full SM Lagrangian, as Connes does), but it addresses a *different* question: not "what is the Lagrangian?" but "why THIS gauge group rather than another?".
