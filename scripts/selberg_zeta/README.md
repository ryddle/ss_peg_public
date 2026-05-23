# SS-PEG ↔ Selberg / Ihara Zeta Connection

Experimental scripts exploring the connection between the SS-PEG framework
(graph-based emergent gauge algebras) and the Selberg/Ihara zeta function.

## Core idea

The SS-PEG evolution graph — where Killing eigenvalues classify edges as cyclic,
free, or decoupled — naturally defines a graph zeta function (Ihara zeta).
The spectral properties of this zeta encode the same information as the
trace-over-cycles formula, connecting:

```
spectrum ↔ closed cycles ↔ holonomies ↔ gauge structure
```

## Scripts

| Script | Description |
|--------|-------------|
| `toy_model_ihara.py` | 3-node cyclic graph: Ihara zeta, poles, trace formula |
| `holonomy_zeta_su2.py` | Extension with SU(2) operators on edges (non-abelian holonomy) |
| `ramanujan_check.py` | Test Ramanujan property for graphs from SS-PEG algebras |

## Key questions

1. Do SS-PEG emergent algebras produce Ramanujan (or near-Ramanujan) graphs?
2. Is the density transition at ρ=1 visible as a phase transition in the zeta?
3. Can the F₄ bifurcation be reinterpreted as a pole-crossing in Z(u)?

## Dependencies

- numpy, matplotlib, torch (same venv as `verification/`)
