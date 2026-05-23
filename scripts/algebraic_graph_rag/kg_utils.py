"""
Algebraic Knowledge Graph Utilities
====================================
Synthetic KG generation, holonomy computation, and algebraic layers
for the Algebraic Graph RAG Phase 1 experiments.

Reuses the su(d) basis from training_ai_test/killing_utils.py.
"""

import sys, os
import torch
import numpy as np
from typing import Tuple, List, Optional
from scipy.linalg import expm

# Allow importing from training_ai_test
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'training_ai_test'))
from killing_utils import (
    _su_basis,
    killing_metric,
    killing_tension,
    DEVICE,
)

# ============================================================
#  GROUND-TRUTH RELATION OPERATORS (Adjoint representation)
# ============================================================

def _adjoint_matrix(U: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """
    Compute the adjoint representation matrix of a unitary U ∈ SU(d).

    Ad(U)_{ab} = 2 * Re(Tr(T_a · U · T_b · U†))

    where {T_a} is the GML basis with Tr(T_a T_b) = 2 δ_{ab}.

    U     : (d, d) complex unitary
    basis : (n_gen, d, d) complex hermitian traceless
    Returns: (n_gen, n_gen) real orthogonal matrix
    """
    n_gen = basis.shape[0]
    U_dag = U.conj().T
    Ad = np.zeros((n_gen, n_gen), dtype=np.float64)
    for a in range(n_gen):
        for b in range(n_gen):
            # Tr(T_a U T_b U†) — factor of 2 from normalization Tr(T_a T_b) = 2δ
            val = np.trace(basis[a] @ U @ basis[b] @ U_dag)
            Ad[a, b] = 0.5 * val.real  # The 2 and 1/2 cancel with norm convention
    return Ad


def generate_relation_operators(
    n_relations: int,
    d: int = 3,
    angle_scale: float = 0.8,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate ground-truth relation operators as adjoint-representation matrices
    of random SU(d) elements.

    Each relation r is generated as:
      1. Pick random algebra element: A_r = Σ_k α_rk T_k  (random coefficients)
      2. Compute unitary: U_r = exp(i · angle_scale · A_r) ∈ SU(d)
      3. Compute adjoint: Q_r = Ad(U_r) ∈ O(n_gen)

    Returns:
      Q_ops  : (n_relations, n_gen, n_gen) real orthogonal matrices
      alphas : (n_relations, n_gen) algebra coefficients (ground truth)
    """
    rng = np.random.default_rng(seed)
    n_gen = d * d - 1
    basis_np = _su_basis(d, 'cpu').numpy()  # (n_gen, d, d) complex

    Q_ops = np.zeros((n_relations, n_gen, n_gen), dtype=np.float64)
    alphas = np.zeros((n_relations, n_gen), dtype=np.float64)

    for r in range(n_relations):
        alpha = rng.standard_normal(n_gen)
        alpha = alpha / np.linalg.norm(alpha)  # unit direction
        alphas[r] = alpha

        # Algebra element A = Σ_k α_k T_k
        A = sum(alpha[k] * basis_np[k] for k in range(n_gen))
        # Unitary U = exp(i · ε · A)  — A is Hermitian, so iA is anti-Hermitian → U is unitary
        U = expm(1j * angle_scale * A)
        # Adjoint representation
        Q_ops[r] = _adjoint_matrix(U, basis_np)

    return Q_ops, alphas


# ============================================================
#  SYNTHETIC KG GENERATION
# ============================================================

def generate_synthetic_kg(
    n_entities: int = 200,
    n_relations: int = 6,
    n_train_triples: int = 5000,
    n_val_triples: int = 1000,
    d: int = 3,
    angle_scale: float = 0.8,
    noise_std: float = 0.05,
    seed: int = 42,
) -> dict:
    """
    Generate a synthetic algebraic knowledge graph.

    Entities are random unit vectors in R^{n_gen}.
    Relations are adjoint-representation operators from SU(d).
    Valid triples: (h, r, t) where e_t ≈ Q_r · e_h.
    Invalid triples: random (h, r, t) not satisfying the relation.

    Returns dict with:
      entity_embeddings : (n_entities, n_gen) ground truth
      relation_ops      : (n_relations, n_gen, n_gen) ground truth
      train_triples     : (n_train, 3) int [head_idx, rel_idx, tail_idx]
      train_labels      : (n_train,) int [0=invalid, 1=valid]
      val_triples       : (n_val, 3)
      val_labels        : (n_val,)
    """
    rng = np.random.default_rng(seed)
    n_gen = d * d - 1

    # Entity embeddings: random unit vectors in R^{n_gen}
    entities = rng.standard_normal((n_entities, n_gen))
    entities = entities / np.linalg.norm(entities, axis=1, keepdims=True)

    # Relation operators
    Q_ops, alphas = generate_relation_operators(n_relations, d, angle_scale, seed)

    # Precompute all transformed entities: Q_r · e_h for all (r, h)
    # transformed[r, h] = Q_r @ entities[h]
    transformed = np.einsum('rij,hj->rhi', Q_ops, entities)  # (n_rel, n_ent, n_gen)

    # Find nearest entity for each (r, h) → candidate valid tail
    # nearest_tail[r, h] = argmin_t ||transformed[r,h] - entities[t]||
    # Exclude self-loops (t != h)
    nearest_tail = np.zeros((n_relations, n_entities), dtype=np.int64)
    nearest_dist = np.zeros((n_relations, n_entities), dtype=np.float64)

    for r in range(n_relations):
        for h in range(n_entities):
            diffs = entities - transformed[r, h]  # (n_ent, n_gen)
            dists = np.linalg.norm(diffs, axis=1)
            dists[h] = 1e10  # exclude self
            nearest_tail[r, h] = np.argmin(dists)
            nearest_dist[r, h] = dists[nearest_tail[r, h]]

    def _generate_split(n_triples, rng_local):
        triples = []
        labels = []
        n_valid = n_triples // 2
        n_invalid = n_triples - n_valid

        # Valid triples: pick (h, r) with small nearest_dist
        for _ in range(n_valid):
            r = rng_local.integers(n_relations)
            h = rng_local.integers(n_entities)
            t = nearest_tail[r, h]
            triples.append([h, r, t])
            labels.append(1)

        # Invalid triples: random (h, r, t) where t != nearest_tail[r,h]
        for _ in range(n_invalid):
            r = rng_local.integers(n_relations)
            h = rng_local.integers(n_entities)
            t = rng_local.integers(n_entities)
            # Ensure it's not accidentally valid
            while t == nearest_tail[r, h] or t == h:
                t = rng_local.integers(n_entities)
            triples.append([h, r, t])
            labels.append(0)

        perm = rng_local.permutation(len(triples))
        triples = np.array(triples)[perm]
        labels = np.array(labels)[perm]
        return triples, labels

    train_triples, train_labels = _generate_split(n_train_triples, rng)
    val_rng = np.random.default_rng(seed + 1000)
    val_triples, val_labels = _generate_split(n_val_triples, val_rng)

    return {
        'entity_embeddings': entities,
        'relation_ops': Q_ops,
        'relation_alphas': alphas,
        'nearest_tail': nearest_tail,
        'nearest_dist': nearest_dist,
        'train_triples': train_triples,
        'train_labels': train_labels,
        'val_triples': val_triples,
        'val_labels': val_labels,
        'd': d,
        'n_gen': n_gen,
    }


# ============================================================
#  CYCLE GENERATION (Consistent & Inconsistent Triangles)
# ============================================================

def generate_cycles(
    relation_ops: np.ndarray,
    n_consistent: int = 200,
    n_inconsistent: int = 200,
    seed: int = 42,
) -> dict:
    """
    Generate consistent and inconsistent triangular cycles.

    A cycle is a triple of relation indices (r1, r2, r3).
    Consistent: Q_{r3} · Q_{r2} · Q_{r1} ≈ I  (holonomy ≈ 0)
    Inconsistent: random triple → holonomy typically >> 0

    To create consistent cycles:
      Pick r1, r2 randomly, then set Q_{r3} = (Q_{r2} · Q_{r1})^{-1} = Q_{r1}^T · Q_{r2}^T.
      Since we only have K discrete relation types, we find the r3 whose Q_{r3}
      is closest to the required inverse. The residual holonomy is the unavoidable
      discretization error.

    Returns dict with:
      consistent_cycles     : (n_con, 3) int relation indices
      consistent_holonomy   : (n_con,) float ground-truth ||H - I||_F
      inconsistent_cycles   : (n_incon, 3) int
      inconsistent_holonomy : (n_incon,) float
    """
    rng = np.random.default_rng(seed)
    n_rel = relation_ops.shape[0]
    n_gen = relation_ops.shape[1]
    I_n = np.eye(n_gen)

    def _holonomy_norm(r1, r2, r3):
        H = relation_ops[r3] @ relation_ops[r2] @ relation_ops[r1]
        return np.linalg.norm(H - I_n, 'fro')

    # --- Consistent cycles ---
    # Strategy: for each (r1, r2), find the r3 that minimizes ||Q_r3 Q_r2 Q_r1 - I||
    # This is equivalent to finding r3 closest to (Q_r2 Q_r1)^{-1} = Q_r1^T Q_r2^T
    best_closures = []
    for r1 in range(n_rel):
        for r2 in range(n_rel):
            target = (relation_ops[r2] @ relation_ops[r1]).T  # inverse
            best_r3, best_dist = -1, 1e10
            for r3 in range(n_rel):
                dist = np.linalg.norm(relation_ops[r3] - target, 'fro')
                if dist < best_dist:
                    best_dist = dist
                    best_r3 = r3
            hol = _holonomy_norm(r1, r2, best_r3)
            best_closures.append((r1, r2, best_r3, hol))

    # Sort by holonomy and take the n_consistent best
    best_closures.sort(key=lambda x: x[3])
    consistent_cycles = []
    consistent_holonomy = []
    seen = set()
    for r1, r2, r3, hol in best_closures:
        key = (r1, r2, r3)
        if key not in seen:
            seen.add(key)
            consistent_cycles.append([r1, r2, r3])
            consistent_holonomy.append(hol)
            if len(consistent_cycles) >= n_consistent:
                break

    # Pad if not enough (repeat with noise)
    while len(consistent_cycles) < n_consistent:
        idx = rng.integers(len(consistent_cycles))
        consistent_cycles.append(consistent_cycles[idx])
        consistent_holonomy.append(consistent_holonomy[idx])

    # --- Inconsistent cycles ---
    # Enumerate ALL possible cycles and rank by holonomy (finite pool: n_rel^3)
    all_cycles_hol = []
    for r1 in range(n_rel):
        for r2 in range(n_rel):
            for r3 in range(n_rel):
                hol = _holonomy_norm(r1, r2, r3)
                all_cycles_hol.append((r1, r2, r3, hol))

    # Sort descending by holonomy → top are most inconsistent
    all_cycles_hol.sort(key=lambda x: -x[3])

    # Consistent holonomy median as separation reference
    con_med = np.median(consistent_holonomy) if consistent_holonomy else 0.5

    inconsistent_cycles = []
    inconsistent_holonomy = []
    for r1, r2, r3, hol in all_cycles_hol:
        key = (r1, r2, r3)
        if key in seen:
            continue  # skip cycles already in the consistent set
        if hol > con_med * 1.2:  # minimal separation from consistent set
            inconsistent_cycles.append([r1, r2, r3])
            inconsistent_holonomy.append(hol)
            seen.add(key)
            if len(inconsistent_cycles) >= n_inconsistent:
                break

    # If still not enough, take any remaining cycles not in consistent set
    if len(inconsistent_cycles) < n_inconsistent:
        for r1, r2, r3, hol in all_cycles_hol:
            key = (r1, r2, r3)
            if key not in seen:
                inconsistent_cycles.append([r1, r2, r3])
                inconsistent_holonomy.append(hol)
                seen.add(key)
                if len(inconsistent_cycles) >= n_inconsistent:
                    break

    return {
        'consistent_cycles': np.array(consistent_cycles[:n_consistent]),
        'consistent_holonomy': np.array(consistent_holonomy[:n_consistent]),
        'inconsistent_cycles': np.array(inconsistent_cycles[:n_inconsistent]),
        'inconsistent_holonomy': np.array(inconsistent_holonomy[:n_inconsistent]),
    }


# ============================================================
#  HOLONOMY COMPUTATION (on learned operators)
# ============================================================

def compute_holonomy(
    relation_matrices: torch.Tensor,
    cycles: np.ndarray,
) -> torch.Tensor:
    """
    Compute holonomy Δ_info = ||Π θ_r - I||_F for each cycle.

    relation_matrices : (n_relations, n_gen, n_gen) learned operators
    cycles            : (n_cycles, k) int, relation indices per cycle
    Returns           : (n_cycles,) holonomy norms
    """
    n_gen = relation_matrices.shape[1]
    I_n = torch.eye(n_gen, dtype=relation_matrices.dtype, device=relation_matrices.device)
    holonomies = []

    for cycle in cycles:
        H = I_n.clone()
        for r_idx in cycle:
            H = relation_matrices[r_idx] @ H
        delta = torch.norm(H - I_n, p='fro')
        holonomies.append(delta)

    return torch.stack(holonomies)


# ============================================================
#  KILLING TENSION FOR RELATION OPERATORS
# ============================================================

def generate_relation_operators_with_closures(
    n_base: int = 12,
    n_closures: int = 12,
    d: int = 3,
    angle_scale: float = 0.25,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, List[Tuple[int, int, int]]]:
    """
    Generate relation operators with explicit closure relations.

    Creates n_base random SU(d) adjoint operators, then for selected pairs
    (r_i, r_j) adds closure relations r_c where Q_c = (Q_j @ Q_i)^{-1},
    guaranteeing that cycle (r_i, r_j, r_c) has EXACTLY zero holonomy:
      Q_c · Q_j · Q_i = I.

    Returns:
      Q_ops          : (n_base + n_closures, n_gen, n_gen) orthogonal matrices
      alphas         : (n_base, n_gen) algebra coefficients (base only)
      closure_cycles : list of (r_i, r_j, r_c) tuples with GT holonomy = 0
    """
    # Generate base relations
    Q_base, alphas = generate_relation_operators(n_base, d, angle_scale, seed)
    n_gen = Q_base.shape[1]

    # Pick pairs for closures
    rng = np.random.default_rng(seed + 7777)
    pairs = []
    seen_pairs = set()
    while len(pairs) < n_closures:
        ri = rng.integers(n_base)
        rj = rng.integers(n_base)
        if (ri, rj) not in seen_pairs:
            seen_pairs.add((ri, rj))
            pairs.append((ri, rj))

    # Generate closure operators
    Q_closures = np.zeros((n_closures, n_gen, n_gen), dtype=np.float64)
    closure_cycles = []
    for c, (ri, rj) in enumerate(pairs):
        # Q_c = (Q_j @ Q_i)^{-1} = (Q_j @ Q_i)^T  (since they're orthogonal)
        Q_closures[c] = (Q_base[rj] @ Q_base[ri]).T
        closure_cycles.append((ri, rj, n_base + c))

    Q_ops = np.concatenate([Q_base, Q_closures], axis=0)
    return Q_ops, alphas, closure_cycles


def generate_cycles_v2(
    relation_ops: np.ndarray,
    closure_cycles: List[Tuple[int, int, int]],
    n_consistent: int = 200,
    n_inconsistent: int = 200,
    seed: int = 42,
) -> dict:
    """
    Generate consistent/inconsistent cycles using closure information.

    Consistent cycles: drawn from the guaranteed-zero-holonomy closure_cycles
      (padded with repetitions if fewer than n_consistent).
    Inconsistent cycles: random triples NOT in the closure set, ranked by holonomy.

    Returns dict with:
      consistent_cycles     : (n_con, 3) int
      consistent_holonomy   : (n_con,) float (should be ~0.0)
      inconsistent_cycles   : (n_incon, 3) int
      inconsistent_holonomy : (n_incon,) float
    """
    rng = np.random.default_rng(seed)
    n_rel = relation_ops.shape[0]
    n_gen = relation_ops.shape[1]
    I_n = np.eye(n_gen)

    def _holonomy_norm(r1, r2, r3):
        H = relation_ops[r3] @ relation_ops[r2] @ relation_ops[r1]
        return np.linalg.norm(H - I_n, 'fro')

    # --- Consistent cycles from closures ---
    closure_set = set()
    consistent_cycles = []
    consistent_holonomy = []
    for ri, rj, rc in closure_cycles:
        hol = _holonomy_norm(ri, rj, rc)
        consistent_cycles.append([ri, rj, rc])
        consistent_holonomy.append(hol)
        closure_set.add((ri, rj, rc))

    # Pad by repeating if needed
    base_n = len(consistent_cycles)
    while len(consistent_cycles) < n_consistent:
        idx = rng.integers(base_n)
        consistent_cycles.append(consistent_cycles[idx])
        consistent_holonomy.append(consistent_holonomy[idx])

    # --- Inconsistent cycles: sample random triples not in closure set ---
    candidate_pool = []
    # Sample a large pool of random triples
    n_candidates = min(n_rel ** 3, 50000)
    seen = set(closure_set)
    attempts = 0
    while len(candidate_pool) < n_candidates and attempts < n_candidates * 5:
        r1, r2, r3 = rng.integers(n_rel, size=3)
        key = (int(r1), int(r2), int(r3))
        if key not in seen:
            hol = _holonomy_norm(r1, r2, r3)
            candidate_pool.append((r1, r2, r3, hol))
            seen.add(key)
        attempts += 1

    # Sort descending by holonomy → take worst closers
    candidate_pool.sort(key=lambda x: -x[3])

    inconsistent_cycles = []
    inconsistent_holonomy = []
    for r1, r2, r3, hol in candidate_pool:
        inconsistent_cycles.append([int(r1), int(r2), int(r3)])
        inconsistent_holonomy.append(hol)
        if len(inconsistent_cycles) >= n_inconsistent:
            break

    return {
        'consistent_cycles': np.array(consistent_cycles[:n_consistent]),
        'consistent_holonomy': np.array(consistent_holonomy[:n_consistent]),
        'inconsistent_cycles': np.array(inconsistent_cycles[:n_inconsistent]),
        'inconsistent_holonomy': np.array(inconsistent_holonomy[:n_inconsistent]),
    }


# ============================================================
#  ENTITY-GROUNDED CYCLES (for path holonomy evaluation)
# ============================================================

def generate_entity_grounded_cycles(
    nearest_tail: np.ndarray,
    closure_cycles,
    n_entities: int,
    n_rel: int,
    n_consistent: int = 200,
    n_inconsistent: int = 200,
    seed: int = 42,
):
    """
    Generate entity-grounded triangular cycles for path holonomy evaluation.

    Each cycle is (h1, r1, h2, r2, h3, r3) specifying both entities and
    relations along the path: h1 --r1--> h2 --r2--> h3 --r3--> h1.

    Consistent: closure cycle (ri, rj, rc) from random starting entity.
    Inconsistent: same path but closure relation replaced with random.

    nearest_tail   : (n_rel, n_entities) int
    closure_cycles : list of (ri, rj, rc) tuples with zero GT holonomy

    Returns:
      consistent_paths   : list of (h1, r1, h2, r2, h3, r3) tuples
      inconsistent_paths : list of (h1, r1, h2, r2, h3, r3) tuples
    """
    rng = np.random.default_rng(seed)
    n_closures = len(closure_cycles)

    consistent_paths = []
    for _ in range(n_consistent):
        ci = rng.integers(n_closures)
        ri, rj, rc = closure_cycles[ci]
        h1 = int(rng.integers(n_entities))
        h2 = int(nearest_tail[ri, h1])
        h3 = int(nearest_tail[rj, h2])
        consistent_paths.append((h1, ri, h2, rj, h3, rc))

    inconsistent_paths = []
    for _ in range(n_inconsistent):
        ci = rng.integers(n_closures)
        ri, rj, rc = closure_cycles[ci]
        h1 = int(rng.integers(n_entities))
        h2 = int(nearest_tail[ri, h1])
        h3 = int(nearest_tail[rj, h2])
        r_rand = int(rng.integers(n_rel))
        while r_rand == rc:
            r_rand = int(rng.integers(n_rel))
        inconsistent_paths.append((h1, ri, h2, rj, h3, r_rand))

    return consistent_paths, inconsistent_paths


# ============================================================
#  KILLING TENSION FOR RELATION OPERATORS (original position)
# ============================================================

def relation_killing_tension(relation_matrices: torch.Tensor) -> torch.Tensor:
    """
    Compute T_K for a set of relation operator matrices.

    Treats each relation matrix θ_r ∈ R^{n×n} as a "generator" and
    computes the Killing form of the set {θ_r}.

    For this to make physical sense, we interpret the relation operators
    as elements near the identity of the group, so their logs are
    approximately proportional to algebra elements.

    relation_matrices : (n_relations, n_gen, n_gen)
    Returns           : scalar T_K
    """
    # Use the matrices directly as generators for the Killing metric
    # K_{ab} = Σ_c Tr([θ_a, θ_c][θ_b, θ_c])
    return killing_tension(relation_matrices)


def relation_ortho_tension(relation_matrices: torch.Tensor) -> torch.Tensor:
    """
    T_ortho = Σ_r ||θ_r θ_r^T - I||_F^2

    Forces each relation operator to be orthogonal.

    relation_matrices : (n_relations, n_gen, n_gen)
    Returns           : scalar
    """
    n_gen = relation_matrices.shape[1]
    I_n = torch.eye(n_gen, dtype=relation_matrices.dtype,
                    device=relation_matrices.device)
    loss = torch.tensor(0.0, dtype=relation_matrices.dtype,
                        device=relation_matrices.device)
    for r in range(relation_matrices.shape[0]):
        theta = relation_matrices[r]
        loss = loss + ((theta @ theta.T - I_n) ** 2).sum()
    return loss
