"""
Topology to Algebra Converter
==============================
Maps a weighted, typed graph (extracted from real data) to a set of
algebraic relation operators Q_r ∈ O(n_gen) in the SU(d) adjoint representation.

This is the core bridge between arbitrary data modalities and the Killing
classifier pipeline (Exp4/Exp5).

The pipeline:
  1. Graph topology (edge types, edge weights, adjacency structure)
  2. Spectral embedding of each edge-type subgraph → direction in su(d)
  3. Angle scaling by edge-type weight → generator A_r
  4. Matrix exponential → unitary U_r = exp(i · A_r) ∈ SU(d)
  5. Adjoint representation → Q_r = Ad(U_r) ∈ O(n_gen)

Design principle: each distinct *relation type* in the data becomes one operator.
The more relation types, the richer the Killing classifier signal (G_RAG-19/21).
"""

import sys
import os
import numpy as np
from typing import List, Dict, Optional, Tuple, Any
from dataclasses import dataclass, field
from scipy.linalg import expm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'training_ai_test'))
from killing_utils import _su_basis


# ============================================================
#  DATA STRUCTURES
# ============================================================

@dataclass
class RelationalGraph:
    """
    A typed, weighted graph extracted from a data source.

    Nodes represent entities (concepts, tokens, pixels, etc.).
    Edges are typed (each type = one algebraic operator).
    """
    nodes: List[str]                                    # node labels
    edges: List[Tuple[int, int, str, float]]           # (src, dst, type, weight)

    @property
    def edge_types(self) -> List[str]:
        types = []
        seen = set()
        for _, _, t, _ in self.edges:
            if t not in seen:
                types.append(t)
                seen.add(t)
        return types

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def n_edges(self) -> int:
        return len(self.edges)

    def adjacency_for_type(self, edge_type: str) -> np.ndarray:
        """Weighted adjacency matrix for a single edge type."""
        n = self.n_nodes
        A = np.zeros((n, n), dtype=np.float64)
        for src, dst, t, w in self.edges:
            if t == edge_type:
                A[src, dst] += w
                A[dst, src] += w     # symmetric: undirected
        return A


@dataclass
class AlgebraicRelation:
    """
    One relation operator Q_r ∈ O(n_gen) together with its metadata.
    """
    name: str
    operator: np.ndarray             # (n_gen, n_gen) real orthogonal matrix
    angle_scale: float               # effective rotation magnitude
    spectral_gap: float              # Fiedler value of the subgraph
    n_edges_in_subgraph: int
    source_modality: str = "unknown"
    extra: Dict[str, Any] = field(default_factory=dict)


# ============================================================
#  ADJOINT REPRESENTATION  (same as kg_utils._adjoint_matrix)
# ============================================================

def _adjoint_matrix(U: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """
    Adjoint representation Ad(U) ∈ O(n_gen) of U ∈ SU(d).

    Ad(U)_{ab} = (1/2) Re Tr(T_a U T_b U†)
    with GML basis satisfying Tr(T_a T_b) = 2 δ_{ab}.
    """
    n_gen = basis.shape[0]
    U_dag = U.conj().T
    Ad = np.zeros((n_gen, n_gen), dtype=np.float64)
    for a in range(n_gen):
        for b in range(n_gen):
            val = np.trace(basis[a] @ U @ basis[b] @ U_dag)
            Ad[a, b] = 0.5 * val.real
    return Ad


# ============================================================
#  SPECTRAL GRAPH → ALGEBRA DIRECTION
# ============================================================

def _spectral_direction(A: np.ndarray, n_gen: int) -> Tuple[np.ndarray, float]:
    """
    Extract a unit direction in R^{n_gen} from a weighted adjacency matrix.

    Method:
    1. Compute normalised Laplacian L_sym = D^{-1/2}(D-A)D^{-1/2}
    2. Eigendecompose L_sym: eigenvectors carry graph topology
    3. Use the k = min(n_gen, n_nodes) Fiedler-and-above eigenvectors
    4. Map the spectral content to a direction vector in R^{n_gen}:
       direction[k] = Σ_{i} λ_i · v_i[k % n_nodes]  (spectral projection)
    5. Normalise to unit vector

    Returns:
      direction : (n_gen,) unit vector, real
      fiedler   : second-smallest eigenvalue (Fiedler value = spectral gap)
    """
    n = A.shape[0]
    if n == 0:
        return np.zeros(n_gen), 0.0

    # Degree vector
    deg = A.sum(axis=1)
    deg = np.where(deg > 0, deg, 1.0)   # avoid division by zero

    # Normalised Laplacian
    D_inv_sqrt = np.diag(1.0 / np.sqrt(deg))
    L = np.eye(n) - D_inv_sqrt @ A @ D_inv_sqrt

    # Eigendecompose (symmetric → real eigenvalues)
    eigvals, eigvecs = np.linalg.eigh(L)

    fiedler = eigvals[1] if n > 1 else 0.0

    # Build direction in R^{n_gen} via spectral projection
    k_use = min(n_gen, n)
    direction = np.zeros(n_gen)
    for k in range(k_use):
        # Map eigenvector index k to generator index (cyclic if n < n_gen)
        gen_idx = k % n_gen
        # Spectral weight: mix eigenvalue magnitude and eigenvector energy
        lam = eigvals[k]
        v_energy = np.sum(eigvecs[:, k] ** 2)   # always 1 for normalised, but keep for clarity
        direction[gen_idx] += np.abs(lam) * v_energy

    # If n_gen > k_use, fill remaining entries with higher-frequency components
    if n_gen > k_use:
        for i in range(k_use, n_gen):
            j = i % k_use
            direction[i] = direction[j] * (0.5 ** (i // k_use))

    norm = np.linalg.norm(direction)
    if norm > 1e-12:
        direction = direction / norm
    else:
        # Fallback: uniform direction
        direction = np.ones(n_gen) / np.sqrt(n_gen)

    return direction, fiedler


def direction_to_operator(
    direction: np.ndarray,
    angle_scale: float,
    d: int,
    basis: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Convert a unit direction vector in R^{n_gen} to an adjoint operator Q ∈ O(n_gen).

    Steps:
    1. Build algebra element A = Σ_k direction[k] * T_k   (Hermitian, traceless)
    2. U = exp(i · angle_scale · A) ∈ SU(d)
    3. Q = Ad(U) ∈ O(n_gen)

    Parameters
    ----------
    direction   : (n_gen,) unit vector
    angle_scale : rotation magnitude in radians
    d           : SU(d) dimension (n_gen = d²-1)
    basis       : (n_gen, d, d) optional pre-computed GML basis
    """
    if basis is None:
        basis = _su_basis(d, 'cpu').numpy()

    n_gen = d * d - 1
    assert len(direction) == n_gen, f"direction length {len(direction)} ≠ n_gen {n_gen}"

    # Algebra element: Hermitian traceless matrix
    A = sum(float(direction[k]) * basis[k] for k in range(n_gen))

    # Unitary: U = exp(i · angle_scale · A)
    U = expm(1j * angle_scale * A)

    # Adjoint representation
    Q = _adjoint_matrix(U, basis)

    return Q


# ============================================================
#  MAIN CONVERTER
# ============================================================

def graph_to_operators(
    graph: RelationalGraph,
    d: int = 3,
    default_angle_scale: float = 0.25,
    modality: str = "unknown",
    min_edges_per_type: int = 1,
    seed: int = 42,
) -> List[AlgebraicRelation]:
    """
    Convert a RelationalGraph to a list of AlgebraicRelation operators.

    One operator per edge type.
    Edge types with fewer than min_edges_per_type edges are skipped.

    Parameters
    ----------
    graph              : input relational graph
    d                  : SU(d) algebra dimension (n_gen = d²-1 = 8 for d=3)
    default_angle_scale: fallback angle scale if edge weights are uniform
    modality           : label for the source data type
    min_edges_per_type : discard relation types with too few edges
    """
    basis = _su_basis(d, 'cpu').numpy()
    n_gen = d * d - 1

    relations: List[AlgebraicRelation] = []

    for edge_type in graph.edge_types:
        # Count edges of this type
        type_edges = [(s, t, w) for s, t, et, w in graph.edges if et == edge_type]
        if len(type_edges) < min_edges_per_type:
            continue

        # Build adjacency matrix for this type
        A = graph.adjacency_for_type(edge_type)

        # Mean edge weight as angle scale modifier
        mean_weight = np.mean([w for _, _, w in type_edges])
        angle_scale = default_angle_scale * np.clip(mean_weight, 0.05, 5.0)

        # Spectral direction
        direction, fiedler = _spectral_direction(A, n_gen)

        # Operator
        Q = direction_to_operator(direction, angle_scale, d, basis)

        relations.append(AlgebraicRelation(
            name=edge_type,
            operator=Q,
            angle_scale=angle_scale,
            spectral_gap=fiedler,
            n_edges_in_subgraph=len(type_edges),
            source_modality=modality,
        ))

    return relations


def operators_to_array(relations: List[AlgebraicRelation]) -> np.ndarray:
    """Stack operators into (K, n_gen, n_gen) array for Killing analysis."""
    if not relations:
        raise ValueError("Empty relations list")
    return np.stack([r.operator for r in relations], axis=0)
