# SS-PEG: Space of States and Projections of the Evolution Graph

**Experimental Program on Emergent Gauge Symmetries from Killing Metric Geometry**

---

## Overview

SS-PEG is an experimental program investigating whether gauge symmetries (SU(3), SU(2), U(1), SO(N), Sp(2N), exceptional algebras) emerge from a purely geometric variational principle acting on the Killing metric, without imposing any algebraic axioms.

### Central Result

The entire algebraic structure of gauge symmetry — commutation relations, Jacobi identities, Casimir operators, dual Coxeter numbers, representation invariants — emerges as the unique output of a geometric variational principle on the Killing metric. No algebraic axioms are required.

---

## Publications

| Paper | Title | arXiv |
|-------|-------|-------|
| **Paper I** | Lie Algebras from the Killing Metric: Gauge Symmetries as Emergent Consequences of Relational Geometry | [pending] |
| **Paper II** | Graph Topology & Spectral Analysis: Ramanujan Properties and the Selberg Zeta Function | [pending] |
| **Paper III** | Standard Model from Graph Topology: Coupling Constants, Mass Spectrum, and Flavor Mixing | [pending] |

---

## Repository Structure

```
sspeg-public/
├── data/                    # Experimental results (npz, logs, json)
│   ├── verification/        # Landscape, bifurcation, E6/E7/E8 verification
│   ├── training_ai/         # Neural algebraic learning experiments
│   ├── algebraic_graph_rag/ # Algebraic graph RAG experiments
│   └── selberg_zeta/        # Selberg/Ihara spectral analysis
│
├── scripts/                 # All Python scripts (reproducible experiments)
│   ├── verification/        # 80+ test/analysis scripts
│   ├── algebraic_graph_rag/ # 18 scripts + utils
│   ├── selberg_zeta/        # 20 scripts
│   ├── from_graph_to_standard_model/  # 40+ scripts
│   ├── sm_mathematical_derivation_from_graph/  # 20 scripts
│   └── training_ai_test/    # 14 scripts + utils
│
├── figures/                 # All generated figures (PNG, SVG)
│   ├── fig1-edge_classes.png
│   ├── fig2-archetype_I_pure_cycle.png
│   ├── ...
│   └── fig20-emergence_minima.png
│
├── experiments/             # Experiment documentation
├── results/                 # Script output (markdown, text)
├── documentation/           # Theoretical documents and synthesis
│   ├── SYNTHESIS.md
│   ├── marco_formal_completo.md
│   ├── tensegrity_ecuation_v3.md
│   ├── from_graph_to_standard_model/
│   ├── selberg_zeta/
│   ├── algebraic_graph_rag/
│   └── training_ai_test/
│
├── README.md
├── LICENSE
└── .gitignore
```

---

## Reproducing Experiments

### Prerequisites

- Python 3.10+
- PyTorch (CUDA 12+)
- NumPy, SciPy, Matplotlib
- NetworkX

### Quick Start

```bash
# Clone the repository
git clone https://github.com/ryddle/sspeg-public.git
cd sspeg-public

# Install dependencies
pip install -r requirements.txt  # if available

# Run all verification tests
python scripts/verification/run_all_tests.py

# Reproduce Paper I results (Killing metric emergence)
python scripts/verification/test_3_1_only_killing.py

# Reproduce landscape analysis
python scripts/verification/_landscape_montecarlo.py
```

### Key Experiments

| Experiment | Script | Description |
|------------|--------|-------------|
| Killing-only ablation | `test_3_1_only_killing.py` | Proves $T_K$ generates algebra, $T_{nc}$ is redundant |
| Density phase transition | `test_density_transition_gpu.py` | Step function at $ho = 1$ |
| Landscape analysis | `_landscape_montecarlo.py` | 128 seeds in so(27) |
| F4 bifurcation | `_test_bifurcation_v2.py` | History-dependent symmetry selection |
| E6 blind discovery | `test_e6_gpu.py` | E6 from random initialization |
| Neural algebraic learning | `training_ai_test/exp1_killing_mlp.py` | 14 neural experiments |
| Selberg spectral analysis | `selberg_zeta/ramanujan_check.py` | Ramanujan property verification |
| Algebraic graph RAG | `algebraic_graph_rag/exp1_holonomy_consistency.py` | 7 RAG experiments |

---

## Key Results Summary

- **51 core tests**, 17 Lie algebras, 4 classical families, 3 exceptional algebras
- **100% success rate** across all classical and exceptional families
- **Density phase transition**: Step function at $ho = 1$ (0% closure for $ho < 1$, 100% at $ho = 1$)
- **Killing-only ablation**: $T_K$ alone generates algebra; $T_{nc}$ is redundant
- **F4 bifurcation**: History-dependent symmetry selection via Adam optimizer state
- **Landscape analysis**: Non-simple products dominate (53.1% of seeds)
- **Dark matter prediction**: Connection without particles → filamentary topology, cored profiles
- **Gravity from entanglement**: Jacobson (1995) derivation in evolution graph framework

---

## License

- **Documentation**: CC-BY-4.0
- **Code**: MIT License

---

## Author

José M. Moreno  
Independent Researcher  
SS-PEG Project

---

## Zenodo DOI

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

---

## Acknowledgments

This work was conducted as an independent research project. Special thanks to the open-source community for PyTorch, NumPy, SciPy, and all computational tools used in this program.
