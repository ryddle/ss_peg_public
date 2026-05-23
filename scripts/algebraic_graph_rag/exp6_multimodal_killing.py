"""
Experiment 6 — Killing Classifier on Multimodal Operators
==========================================================
Extension of Exp4 using REAL operators derived from multimodal data
(text, SVG, audio) instead of synthetic SU(3) adjoint operators.

Research question:
  Do operators extracted from real data via spectral-graph embedding
  exhibit Killing structure separable from random SU(d)-adjoint operators —
  and can cycle loss enhance that separation just as in the synthetic case?

Setup:
  - Type M (multimodal): operators from text + SVG + audio encoders
    → directions derived from spectral structure of real data graphs
  - Type B (baseline): random SU(d)-adjoint operators with random directions
    → same algebraic subspace as Type M, but no real-data structure

Key design choice (fix over original Exp6 attempt):
  Both Type M and Type B operators now live in the adjoint image of SU(d)
  inside SO(n_gen).  This ensures a fair comparison: the only difference is
  where the generator direction comes from (real data vs uniform random).
  Using random SO(n_gen) operators (Exp4 Type B) was unfair because those
  span a larger subspace (so(8), dim 28) than adjoint SU(3) (dim 8),
  causing K_rr polarity inversion unrelated to data structure.

Protocol (identical to Exp4):
  Phase 0: Extract multimodal operators; mine approximate closure cycles
  Phase 1: Train with L_triple + T_K + T_ortho (NO cycle loss)
  Phase 2: Killing spectrum analysis → classify M vs B
  Phase 3: Selective cycle loss on Type M subset

New vs Exp4:
  - Operators are NOT algebraically exact closures (only approximate)
  - Closure cycles are mined from the operators themselves using a
    nearest-product search (|Q_k - Q_j @ Q_i|_F minimisation)
  - Same AlgKGModel, train_phase, analyze_killing_spectrum, evaluate_holonomy

Motivated by G_RAG-19/21 and the multimodal encoder pipeline:
  "If real structure maps to O(n_gen) operators, the Killing metric
   should distinguish them from generic noise just as it does for
   algebraic operators."
"""

import sys, os, time, copy
import torch
import torch.nn as nn
import numpy as np
from scipy.linalg import expm
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(__file__))

from kg_utils import (
    generate_cycles_v2,
    generate_entity_grounded_cycles,
    compute_holonomy,
    relation_killing_tension,
    relation_ortho_tension,
    DEVICE,
)
# Reuse all machinery from Exp4
from exp4_killing_classifier import (
    AlgKGModel,
    train_phase,
    analyze_killing_spectrum,
    evaluate_holonomy,
)
# Multimodal encoders
from relational_encoder import MultimodalEncoder
from topology_to_algebra import operators_to_array, direction_to_operator, _su_basis

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'training_ai_test'))
from killing_utils import killing_metric

sys.stdout.reconfigure(encoding='utf-8')


# ============================================================
#  SAMPLE MULTIMODAL DATA (self-contained, no external files)
# ============================================================

SAMPLE_TEXTS = [
    # Particle physics — dense lexical overlap for high-connectivity cooccurrence graph
    """Quarks carry colour charge and are confined inside hadrons by gluons.
    Gluons are the mediators of the strong nuclear force between colour-charged quarks.
    Protons and neutrons are hadrons composed of three quarks each bound by gluons.
    The strong force confines quarks inside protons, neutrons, and other hadrons.
    Mesons consist of a quark and an antiquark bound together by gluons.
    Pions are the lightest mesons and mediate the residual nuclear binding force.
    Kaons contain strange quarks and exhibit CP violation in their weak decay.
    Baryons are hadrons with three quarks, the most familiar being protons and neutrons.
    Colour charge of quarks comes in red, green, and blue varieties under SU(3).
    Gluons themselves carry colour charge and interact with each other via gluon vertices.
    Asymptotic freedom means quarks behave as nearly free particles at very short distances.
    Confinement ensures that isolated free quarks are never observed in nature.
    Quantum chromodynamics is the gauge field theory of quarks and gluons.
    The strong coupling constant governs the interaction strength between quarks and gluons.
    Lattice QCD computes hadron masses from first principles on a discrete spacetime lattice.
    Flavour symmetry relates up, down, and strange quarks in the approximate quark model.
    Isospin symmetry connects protons and neutrons as an SU(2) doublet of the strong force.
    The quark model classifies all hadrons by their quark content, spin, and parity.
    Deep inelastic scattering experiments revealed the quark substructure inside protons.
    Colour confinement is a direct consequence of the non-Abelian gauge structure of gluons.""",

    # Symmetry breaking and gauge theory — overlapping vocabulary for dense cooccurrence
    """Symmetry breaking gives mass to gauge bosons via the Higgs mechanism.
    The Higgs field permeates all of space with a non-zero vacuum expectation value.
    Gauge bosons acquire mass by absorbing Goldstone bosons during spontaneous symmetry breaking.
    The W and Z bosons are the massive gauge bosons of the electroweak interaction.
    Photons remain massless because electromagnetic gauge symmetry remains unbroken.
    The Standard Model gauge group is SU(3) cross SU(2) cross U(1).
    Spontaneous symmetry breaking occurs when the ground state of a field breaks a symmetry.
    Goldstone bosons arise whenever a continuous global symmetry is spontaneously broken.
    Electroweak unification merges electromagnetism and the weak nuclear force above the Z mass.
    The Higgs boson is the physical particle arising from the Higgs field after symmetry breaking.
    Yukawa couplings between fermions and the Higgs field generate the fermion mass spectrum.
    The fermion masses in the Standard Model all arise from Higgs Yukawa interactions.
    SU(2) isospin symmetry relates up and down quarks in electroweak doublet representations.
    The weak force couples exclusively to left-handed fermion doublets through W and Z bosons.
    Parity violation distinguishes left-handed from right-handed fermion fields in weak decay.
    Anomaly cancellation requires quarks and leptons to appear in complete Standard Model generations.
    Running coupling constants vary with the energy scale due to radiative quantum corrections.
    Renormalization group equations describe how coupling constants flow with the energy scale.
    Grand unified theories extend the Standard Model gauge group into a single larger gauge group.
    Supersymmetry introduces a symmetry that relates every boson to a fermionic superpartner.""",

    # Spectral graph theory, algebraic topology, and holonomy — rich shared vocabulary
    """A graph consists of vertices connected by weighted edges encoding relational structure.
    The graph Laplacian matrix encodes the full connectivity structure of a weighted graph.
    Spectral graph theory studies graphs through the eigenvalues and eigenvectors of the Laplacian.
    The Fiedler value is the second smallest eigenvalue of the normalised graph Laplacian.
    The Fiedler vector partitions graph vertices into two clusters by spectral bisection.
    Connected graphs have a strictly positive Fiedler value and a single connected component.
    Algebraic connectivity measures how well the vertices of a graph are linked by edges.
    Random walks on graphs converge to stationary distributions determined by graph eigenvalues.
    Spectral clustering uses leading eigenvectors of the Laplacian to partition graph vertices.
    Graph Laplacians generalise the continuous Laplace-Beltrami operator from smooth manifolds.
    Holonomy measures the failure of parallel transport of a vector around a closed graph cycle.
    A connection on a fibre bundle determines how vectors are transported along graph edges.
    Curvature is the infinitesimal holonomy of a connection around an infinitesimally small loop.
    The Killing metric on a Lie algebra is defined via the trace of the adjoint representation.
    Lie algebras are vector spaces equipped with an antisymmetric bilinear Lie bracket operation.
    The adjoint representation maps each Lie algebra element to a linear operator on the algebra.
    Cycle spaces of graphs correspond to the first homology groups of the underlying graph topology.
    Persistent homology tracks topological features of a graph across a filtration of subcomplexes.
    Graph neural networks propagate vertex features along edges using spectral graph convolutions.
    The graph eigenvalue spectrum encodes geometric properties such as diameter and vertex expansion.""",
]

SAMPLE_SVG = """<svg xmlns="http://www.w3.org/2000/svg" width="500" height="420">
  <defs>
    <style>
      .node { stroke-width: 2; }
      .edge { stroke: #888; stroke-width: 1.5; fill: none; }
      .label { font-size: 10px; text-anchor: middle; fill: #fff; font-weight: bold; }
      .group-a { fill: #e74c3c; stroke: #c0392b; }
      .group-b { fill: #2ecc71; stroke: #27ae60; }
      .group-c { fill: #3498db; stroke: #2980b9; }
      .hub { stroke-width: 4; }
    </style>
  </defs>
  <g id="layer_edges">
    <line class="edge" x1="80"  y1="80"  x2="160" y2="140"/>
    <line class="edge" x1="80"  y1="80"  x2="50"  y2="210"/>
    <line class="edge" x1="160" y1="140" x2="50"  y2="210"/>
    <line class="edge" x1="50"  y1="210" x2="110" y2="340"/>
    <line class="edge" x1="160" y1="140" x2="110" y2="340"/>
    <line class="edge" x1="250" y1="50"  x2="350" y2="110"/>
    <line class="edge" x1="250" y1="50"  x2="270" y2="200"/>
    <line class="edge" x1="350" y1="110" x2="270" y2="200"/>
    <line class="edge" x1="270" y1="200" x2="300" y2="340"/>
    <line class="edge" x1="350" y1="110" x2="300" y2="340"/>
    <line class="edge" x1="430" y1="80"  x2="460" y2="190"/>
    <line class="edge" x1="460" y1="190" x2="420" y2="330"/>
    <line class="edge" x1="430" y1="80"  x2="420" y2="330"/>
    <line class="edge" x1="160" y1="140" x2="250" y2="50"/>
    <line class="edge" x1="270" y1="200" x2="460" y2="190"/>
    <line class="edge" x1="110" y1="340" x2="300" y2="340"/>
    <line class="edge" x1="300" y1="340" x2="420" y2="330"/>
    <line class="edge" x1="80"  y1="80"  x2="250" y2="50"/>
    <line class="edge" x1="430" y1="80"  x2="350" y2="110"/>
  </g>
  <g id="layer_group_a">
    <circle id="n1" class="node group-a hub" cx="80"  cy="80"  r="22"/>
    <circle id="n2" class="node group-a"     cx="160" cy="140" r="17"/>
    <circle id="n3" class="node group-a"     cx="50"  cy="210" r="17"/>
    <circle id="n4" class="node group-a"     cx="110" cy="340" r="17"/>
  </g>
  <g id="layer_group_b">
    <circle id="n5" class="node group-b hub" cx="250" cy="50"  r="22"/>
    <circle id="n6" class="node group-b"     cx="350" cy="110" r="17"/>
    <circle id="n7" class="node group-b"     cx="270" cy="200" r="17"/>
    <circle id="n8" class="node group-b"     cx="300" cy="340" r="17"/>
  </g>
  <g id="layer_group_c">
    <circle id="n9"  class="node group-c hub" cx="430" cy="80"  r="22"/>
    <circle id="n10" class="node group-c"     cx="460" cy="190" r="17"/>
    <circle id="n11" class="node group-c"     cx="420" cy="330" r="17"/>
    <rect   id="n12" class="node group-c"     x="395"  y="175"  width="28" height="18" rx="4"/>
  </g>
  <g id="layer_labels">
    <text class="label" x="80"  y="84">A1</text>
    <text class="label" x="160" y="144">A2</text>
    <text class="label" x="50"  y="214">A3</text>
    <text class="label" x="110" y="344">A4</text>
    <text class="label" x="250" y="54">B1</text>
    <text class="label" x="350" y="114">B2</text>
    <text class="label" x="270" y="204">B3</text>
    <text class="label" x="300" y="344">B4</text>
    <text class="label" x="430" y="84">C1</text>
    <text class="label" x="460" y="194">C2</text>
    <text class="label" x="420" y="334">C3</text>
    <text class="label" x="409" y="187">C4</text>
  </g>
</svg>"""


def _make_synthetic_audio(duration_sec: float = 3.0, sr: int = 22050, seed: int = 42):
    """Generate a multi-frequency synthetic audio signal."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, duration_sec, int(duration_sec * sr))
    # Chord-like: fundamental + harmonics with phase offsets
    signal = (
        1.0  * np.sin(2 * np.pi * 220 * t + 0.0) +
        0.7  * np.sin(2 * np.pi * 330 * t + 0.3) +
        0.5  * np.sin(2 * np.pi * 440 * t + 0.7) +
        0.3  * np.sin(2 * np.pi * 660 * t + 1.1) +
        0.15 * np.sin(2 * np.pi * 880 * t + 1.8) +
        0.05 * rng.standard_normal(len(t))
    )
    return signal / np.abs(signal).max(), sr


# ============================================================
#  TYPE B: RANDOM SU(d)-ADJOINT OPERATORS
# ============================================================

def generate_random_adjoint_operators(
    n_ops: int,
    n_gen: int,
    d: int = 3,
    angle_scale: float = 0.25,
    seed: int = 9999,
) -> np.ndarray:
    """
    Generate random operators that live in the SAME algebraic subspace as
    multimodal operators: the adjoint image of SU(d) inside SO(n_gen).

    Each operator is built from a uniformly random unit direction in R^{n_gen}
    via the same `direction_to_operator` pipeline used by the encoders.
    This ensures Type M vs Type B differ only in the *origin* of the
    direction vector (real data structure vs random noise), not in the
    ambient algebraic space.
    """
    basis = _su_basis(d, 'cpu').numpy()
    rng = np.random.default_rng(seed)
    ops = np.zeros((n_ops, n_gen, n_gen), dtype=np.float64)
    for r in range(n_ops):
        direction = rng.standard_normal(n_gen)
        direction = direction / np.linalg.norm(direction)
        ops[r] = direction_to_operator(direction, angle_scale, d, basis)
    return ops


# ============================================================
#  APPROXIMATE CLOSURE MINING
# ============================================================

def mine_closure_cycles(
    ops: np.ndarray,
    n_cycles: int = 50,
    max_closure_hol: float | None = None,
    seed: int = 42,
) -> list:
    """
    Mine approximate algebraic closure cycles from a set of operators.

    A closure cycle (i, j, k) satisfies Q_k ≈ Q_j @ Q_i, i.e.
    the holonomy ||Q_k @ Q_j @ Q_i - I||_F is small.

    When max_closure_hol is None (default), the threshold is set automatically
    to the 25th percentile of all triple holonomies so that at most 25% of
    triples qualify as "consistent" — guaranteeing that inconsistent triples
    exist regardless of the number of operators.

    Returns list of (i, j, k) tuples sorted by ascending holonomy.
    """
    n = ops.shape[0]
    I = np.eye(ops.shape[1])
    cycles = []

    for i in range(n):
        for j in range(n):
            if j == i:
                continue
            product = ops[j] @ ops[i]
            for k in range(n):
                if k == i or k == j:
                    continue
                # Forward closure: Q_k^{-1} @ Q_j @ Q_i ≈ I
                # Since Q_k is orthogonal, Q_k^{-1} = Q_k^T
                residual = ops[k].T @ product
                hol = np.linalg.norm(residual - I, 'fro')
                cycles.append((i, j, k, hol))

    # Sort by holonomy (best closures first)
    cycles.sort(key=lambda x: x[3])

    # Auto threshold: 25th percentile so at most 25% qualify, leaving room
    # for inconsistent triples in build_cycle_data_from_mined.
    if max_closure_hol is None:
        all_hols = [h for _, _, _, h in cycles]
        max_closure_hol = float(np.percentile(all_hols, 25)) if all_hols else 2.5

    # Keep only cycles below threshold and limit count
    good = [(i, j, k) for i, j, k, h in cycles if h < max_closure_hol]
    if not good:
        # Fallback: take top n_cycles // 4 regardless of threshold
        good = [(i, j, k) for i, j, k, _ in cycles[:max(n_cycles // 4, 1)]]

    # Deduplicate: each triple (i,j,k) appears at most once
    seen = set()
    unique = []
    for cyc in good:
        key = tuple(cyc)
        if key not in seen:
            seen.add(key)
            unique.append(cyc)

    return unique[:n_cycles]


def build_cycle_data_from_mined(
    ops: np.ndarray,
    closure_cycles: list,
    n_consistent: int = 150,
    n_inconsistent: int = 150,
    seed: int = 42,
) -> dict:
    """
    Build consistent/inconsistent cycle tensors from mined closures.

    Consistent:   cycles from mine_closure_cycles (low holonomy)
    Inconsistent: random triple permutations (high expected holonomy)
    """
    rng = np.random.default_rng(seed)
    n = ops.shape[0]
    I = np.eye(ops.shape[1])

    # Consistent: sample from mined cycles
    if len(closure_cycles) >= n_consistent:
        idxs = rng.choice(len(closure_cycles), n_consistent, replace=False)
        con = [closure_cycles[i] for i in idxs]
    else:
        con = closure_cycles * (n_consistent // max(len(closure_cycles), 1) + 1)
        con = con[:n_consistent]

    # Inconsistent: random triples that are NOT in the closure list
    closure_set = set(map(tuple, closure_cycles))
    incon = []
    attempts = 0
    while len(incon) < n_inconsistent and attempts < 50_000:
        tri = tuple(rng.choice(n, 3, replace=False))
        if tri not in closure_set:
            incon.append(list(tri))
        attempts += 1

    # Compute holonomies
    con_arr = np.array(con, dtype=np.int64)
    incon_arr = np.array(incon, dtype=np.int64)

    def _hol(tri_arr):
        hols = np.zeros(len(tri_arr))
        for idx, (i, j, k) in enumerate(tri_arr):
            H = ops[k] @ ops[j] @ ops[i]
            hols[idx] = np.linalg.norm(H - I, 'fro')
        return hols

    return {
        'consistent_cycles': con_arr,
        'inconsistent_cycles': incon_arr,
        'consistent_holonomy': _hol(con_arr),
        'inconsistent_holonomy': _hol(incon_arr),
    }


# ============================================================
#  CONFIG
# ============================================================

D = 3
ANGLE_SCALE = 0.5   # 0.25 caused K_rr polarity flip (eff. angle ~0.093 → K_rr ≈1e-4 for all)
N_ENTITIES = 200
N_TRAIN, N_VAL = 8000, 2000
PHASE1_EPOCHS = 800
PHASE3_EPOCHS = 400
LR = 1e-3
BATCH_SIZE = 256
HEAD_HIDDEN = 32
LAMBDA_K, LAMBDA_ORTHO = 0.5, 1.0
LAMBDA_CYCLE = 1.0
MARGIN = 1.0
GRAD_CLIP = 1.0
SEED = 42

N_RANDOM_B = 8     # random Type B operators (adjoint SU(d), same subspace as M)
MIN_GAP = 0.0      # keep all operators (gap=0 → uniform direction, still valid for volume comparison)
PHASE3_LOG_EVERY = 50  # evaluate K_rr AUC every N epochs during Phase 3 (feature B)
SR = 22050


# ============================================================
#  MAIN EXPERIMENT
# ============================================================

def run_experiment(seed: int = SEED, quick: bool = False,
                   gap_filter: float | None = None,
                   angle_scale: float | None = None):
    global PHASE1_EPOCHS, PHASE3_EPOCHS, N_TRAIN, N_VAL, ANGLE_SCALE
    if quick:
        PHASE1_EPOCHS = 200
        PHASE3_EPOCHS = 100
        N_TRAIN = 2000
        N_VAL = 500
    if angle_scale is not None:
        ANGLE_SCALE = angle_scale
    _gap = MIN_GAP if gap_filter is None else gap_filter

    print("=" * 70)
    print("  Exp6 — Killing Classifier on Multimodal Operators")
    print("=" * 70)
    print(f"  Device : {DEVICE}")
    print(f"  Seed   : {seed}")
    print(f"  d      : {D}  (n_gen = {D*D-1})")
    print(f"  Type M : operators from text + SVG + audio encoders")
    print(f"  Type B : {N_RANDOM_B} random SO(n_gen) operators")

    t0 = time.time()
    n_gen = D * D - 1

    # ========================================================
    #  PHASE 0A: EXTRACT MULTIMODAL OPERATORS
    # ========================================================

    print("\n--- Extracting multimodal operators ---")
    enc = MultimodalEncoder(d=D, angle_scale=ANGLE_SCALE, seed=seed)

    for i, txt in enumerate(SAMPLE_TEXTS):
        enc.encode_text(txt, use_spacy=False, window=4, min_edges_per_type=2)
        print(f"  text[{i}]: +{enc.n_relations} relations so far")

    n_after_text = enc.n_relations

    enc.encode_svg(SAMPLE_SVG, use_geometric=True, min_edges_per_type=1)
    print(f"  svg:      +{enc.n_relations - n_after_text} new relations "
          f"({enc.n_relations} total)")

    n_after_svg = enc.n_relations

    signal, sr = _make_synthetic_audio(duration_sec=3.0, sr=SR, seed=seed)
    enc.encode_audio(signal, sr, n_fft=512, hop_length=256,
                     n_freq_bands=16, max_frames=64,
                     use_chroma=True, min_edges_per_type=3)
    print(f"  audio:    +{enc.n_relations - n_after_svg} new relations "
          f"({enc.n_relations} total)")

    # Optional: filter out operators whose spectral direction is effectively random.
    # With _gap=0.0 all operators are kept (gap=0 → uniform direction, still valid
    # because angle_scale is calibrated to the edge weight distribution).
    if _gap > 0:
        before = enc.n_relations
        enc._relations = [r for r in enc._relations if r.spectral_gap >= _gap]
        removed = before - enc.n_relations
        if removed:
            print(f"  [gap filter] Removed {removed} operators with gap < {_gap} "
                  f"({enc.n_relations} remaining)")

    if enc.n_relations < 2:
        raise RuntimeError("Not enough multimodal operators extracted — "
                           "check encoder outputs.")

    Q_M_arr, relations_M = enc.to_kg_operators()
    N_M = Q_M_arr.shape[0]
    print(f"\n  Total multimodal operators : {N_M}")
    print(f"  Source modalities : "
          f"{set(r.source_modality for r in relations_M)}")
    print(f"  Relation names    :")
    for rel in relations_M:
        print(f"    [{rel.source_modality:8s}]  {rel.name:30s}  "
              f"gap={rel.spectral_gap:.4f}  edges={rel.n_edges_in_subgraph}")

    # ========================================================
    #  PHASE 0B: TYPE B OPERATORS
    # ========================================================

    # Match type B angle_scale to mean angle_scale of actual M operators so that
    # ||Q_B - I|| ≈ ||Q_M - I|| on average — fair algebraic comparison.
    mean_angle_scale_M = float(np.mean([r.angle_scale for r in relations_M]))
    print(f"\n--- Generating {N_RANDOM_B} Type B (random SU(d)-adjoint) operators ---")
    print(f"  Using mean M angle_scale = {mean_angle_scale_M:.4f} for type B")
    Q_B = generate_random_adjoint_operators(
        N_RANDOM_B, n_gen, d=D, angle_scale=mean_angle_scale_M, seed=seed + 5000
    )
    I_np = np.eye(n_gen)
    dist_M = np.mean([np.linalg.norm(Q_M_arr[r] - I_np, 'fro') for r in range(N_M)])
    dist_B = np.mean([np.linalg.norm(Q_B[r] - I_np, 'fro') for r in range(N_RANDOM_B)])
    print(f"  ||Q - I|| mean: Type M={dist_M:.4f}, Type B={dist_B:.4f}")
    print(f"  (Both operators live in adjoint image of SU({D}) ⊆ SO({n_gen}), matched magnitude)")

    # Concatenate
    Q_all = np.concatenate([Q_M_arr, Q_B], axis=0)
    N_REL = Q_all.shape[0]
    type_labels = np.array([0] * N_M + [1] * N_RANDOM_B)
    print(f"  Total relations: {N_REL}  (M={N_M}, B={N_RANDOM_B})")

    # ========================================================
    #  PHASE 0C: MINE CLOSURE CYCLES FROM TYPE M
    # ========================================================

    print("\n--- Mining approximate closure cycles from Type M operators ---")
    # Use auto percentile threshold: only the best 25% of triples are "consistent",
    # ensuring there is always room for inconsistent triples regardless of N_M.
    raw_closures = mine_closure_cycles(
        Q_M_arr,
        n_cycles=300,
        max_closure_hol=None,  # auto: use 25th percentile of all holonomies
        seed=seed,
    )
    print(f"  Raw mined closure candidates: {len(raw_closures)}")

    # Show top 10 best closures
    hols_mined = []
    for (i, j, k) in raw_closures[:10]:
        H = Q_M_arr[k] @ Q_M_arr[j] @ Q_M_arr[i]
        hol = np.linalg.norm(H - I_np, 'fro')
        hols_mined.append(hol)
    if hols_mined:
        print(f"  Best 10 closure holonomies: min={min(hols_mined):.4f}, "
              f"max={max(hols_mined):.4f}, mean={np.mean(hols_mined):.4f}")
    else:
        print("  WARNING: No closure cycles mined — using random triples")

    # Compare: random triples holonomy
    rng0 = np.random.default_rng(seed)
    rand_hols = []
    for _ in range(50):
        i, j, k = rng0.choice(N_M, 3, replace=False)
        H = Q_M_arr[k] @ Q_M_arr[j] @ Q_M_arr[i]
        rand_hols.append(np.linalg.norm(H - I_np, 'fro'))
    print(f"  Random triple holonomies (baseline): mean={np.mean(rand_hols):.4f}")

    # Shift indices: in Q_all, Type M is 0..N_M-1 (no change needed)
    cycle_data = build_cycle_data_from_mined(
        Q_M_arr,
        raw_closures,
        n_consistent=150,
        n_inconsistent=150,
        seed=seed,
    )
    print(f"  Consistent cycles:   {len(cycle_data['consistent_cycles'])}, "
          f"mean hol={cycle_data['consistent_holonomy'].mean():.4f}")
    print(f"  Inconsistent cycles: {len(cycle_data['inconsistent_cycles'])}, "
          f"mean hol={cycle_data['inconsistent_holonomy'].mean():.4f}")
    hol_sep_gt = (cycle_data['inconsistent_holonomy'].mean() /
                  max(cycle_data['consistent_holonomy'].mean(), 1e-10))
    print(f"  GT holonomy separation: {hol_sep_gt:.2f}x")

    con_cycles_t = torch.tensor(
        cycle_data['consistent_cycles'], dtype=torch.long, device=DEVICE)
    incon_cycles_t = torch.tensor(
        cycle_data['inconsistent_cycles'], dtype=torch.long, device=DEVICE)

    # ========================================================
    #  PHASE 0D: GROUND-TRUTH KILLING ANALYSIS
    # ========================================================

    print("\n--- Ground truth Killing analysis (on Q_all) ---")
    Q_all_t = torch.tensor(Q_all, dtype=torch.float64, device=DEVICE)
    with torch.no_grad():
        K_gt = killing_metric(Q_all_t).cpu().numpy()
    K_gt_diag = np.diag(K_gt)
    print(f"  GT K_rr  Type M: mean={K_gt_diag[:N_M].mean():.4f}, "
          f"std={K_gt_diag[:N_M].std():.4f}")
    print(f"  GT K_rr  Type B: mean={K_gt_diag[N_M:].mean():.4f}, "
          f"std={K_gt_diag[N_M:].std():.4f}")
    try:
        gt_k_auc = roc_auc_score(type_labels, K_gt_diag)
    except ValueError:
        gt_k_auc = float('nan')
    print(f"  GT K_rr AUC (M vs B): {gt_k_auc:.4f}")

    # ========================================================
    #  PHASE 0E: SYNTHETIC KG GENERATION
    # ========================================================

    print("\n--- Generating synthetic KG (mixed Type M + B) ---")
    rng = np.random.default_rng(seed)
    entities = rng.standard_normal((N_ENTITIES, n_gen))
    entities = entities / np.linalg.norm(entities, axis=1, keepdims=True)

    transformed = np.einsum('rij,hj->rhi', Q_all, entities)
    nearest_tail = np.zeros((N_REL, N_ENTITIES), dtype=np.int64)
    nearest_dist = np.zeros((N_REL, N_ENTITIES), dtype=np.float64)
    for r in range(N_REL):
        for h in range(N_ENTITIES):
            diffs = entities - transformed[r, h]
            dists = np.linalg.norm(diffs, axis=1)
            dists[h] = 1e10
            nearest_tail[r, h] = np.argmin(dists)
            nearest_dist[r, h] = dists[nearest_tail[r, h]]

    print(f"  NN dist: M mean={nearest_dist[:N_M].mean():.4f}, "
          f"B mean={nearest_dist[N_M:].mean():.4f}")

    def _gen_split(n_triples, rng_local):
        triples, labels = [], []
        n_valid = n_triples // 2
        for _ in range(n_valid):
            r = rng_local.integers(N_REL)
            h = rng_local.integers(N_ENTITIES)
            triples.append([h, r, nearest_tail[r, h]])
            labels.append(1)
        for _ in range(n_triples - n_valid):
            r = rng_local.integers(N_REL)
            h = rng_local.integers(N_ENTITIES)
            t = rng_local.integers(N_ENTITIES)
            while t == nearest_tail[r, h] or t == h:
                t = rng_local.integers(N_ENTITIES)
            triples.append([h, r, t])
            labels.append(0)
        perm = rng_local.permutation(len(triples))
        return np.array(triples)[perm], np.array(labels)[perm]

    train_triples, train_labels = _gen_split(N_TRAIN, np.random.default_rng(seed + 100))
    val_triples, val_labels = _gen_split(N_VAL, np.random.default_rng(seed + 200))
    print(f"  Train: {len(train_triples)}  |  Val: {len(val_triples)}")

    # Entity-grounded cycle paths
    print("\n--- Generating entity-grounded cycle paths ---")
    con_paths, incon_paths = generate_entity_grounded_cycles(
        nearest_tail, raw_closures,
        N_ENTITIES, N_REL,
        n_consistent=200, n_inconsistent=200, seed=seed + 500,
    )
    gt_con_hol = np.zeros(len(con_paths))
    gt_incon_hol = np.zeros(len(incon_paths))
    for i, (h1, r1, h2, r2, h3, r3) in enumerate(con_paths):
        v = Q_all[r3] @ Q_all[r2] @ Q_all[r1] @ entities[h1]
        gt_con_hol[i] = np.linalg.norm(v - entities[h1])
    for i, (h1, r1, h2, r2, h3, r3) in enumerate(incon_paths):
        v = Q_all[r3] @ Q_all[r2] @ Q_all[r1] @ entities[h1]
        gt_incon_hol[i] = np.linalg.norm(v - entities[h1])

    gt_path_sep = gt_incon_hol.mean() / max(gt_con_hol.mean(), 1e-10)
    print(f"  Consistent paths  : {len(con_paths)}, GT hol={gt_con_hol.mean():.4f}")
    print(f"  Inconsistent paths: {len(incon_paths)}, GT hol={gt_incon_hol.mean():.4f}")
    print(f"  GT path separation: {gt_path_sep:.2f}x")

    # ========================================================
    #  PHASE 1: TRIPLE-ONLY TRAINING
    # ========================================================

    print(f"\n{'='*70}")
    print(f"  PHASE 1: Triple-only Training ({PHASE1_EPOCHS} ep)")
    print(f"  L_triple + T_K + T_ortho  (NO cycle loss)")
    print(f"{'='*70}")

    torch.manual_seed(seed)
    model = AlgKGModel(N_ENTITIES, N_REL, n_gen, HEAD_HIDDEN).to(DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)

    hist_p1 = train_phase(
        model, optimizer,
        train_triples, train_labels, val_triples, val_labels,
        epochs=PHASE1_EPOCHS, batch_size=BATCH_SIZE,
        lambda_k=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
        lambda_cycle=0.0,
        con_cycles_t=con_cycles_t, incon_cycles_t=incon_cycles_t,
        grad_clip=GRAD_CLIP, label="Phase1", log_every=100,
    )

    print(f"\n  Phase 1 holonomy evaluation (without cycle loss):")
    hol_p1 = evaluate_holonomy(
        model, cycle_data, con_paths, incon_paths,
        gt_con_hol, gt_incon_hol, label="Phase1"
    )

    # ========================================================
    #  PHASE 2: KILLING SPECTRUM ANALYSIS
    # ========================================================

    print(f"\n{'='*70}")
    print(f"  PHASE 2: Killing Spectrum Analysis")
    print(f"{'='*70}")

    analysis = analyze_killing_spectrum(model, type_labels, N_M, N_RANDOM_B)

    print(f"\n  Eigenvalue spectrum of learned Killing matrix:")
    eigs = analysis['eigenvalues']
    n_neg = (eigs < -1e-6).sum()
    n_zero = (np.abs(eigs) < 1e-6).sum()
    n_pos = (eigs > 1e-6).sum()
    print(f"    neg={n_neg}  zero={n_zero}  pos={n_pos}")
    print(f"    Top 5 (abs): {eigs[:5]}")

    print(f"\n  K_rr (Killing diagonal) per type:")
    print(f"    Type M (multimodal): mean={analysis['K_diag_A'].mean():.4f}  "
          f"std={analysis['K_diag_A'].std():.4f}  "
          f"range=[{analysis['K_diag_A'].min():.4f}, {analysis['K_diag_A'].max():.4f}]")
    print(f"    Type B (random):     mean={analysis['K_diag_B'].mean():.4f}  "
          f"std={analysis['K_diag_B'].std():.4f}  "
          f"range=[{analysis['K_diag_B'].min():.4f}, {analysis['K_diag_B'].max():.4f}]")

    print(f"\n  Classification AUCs (Phase 1 — key result):")
    print(f"    K_rr ROC-AUC:         {analysis['k_diag_auc']:.4f}")
    print(f"    |K_rr| ROC-AUC:       {analysis['k_diag_abs_auc']:.4f}")
    print(f"    CommTotal ROC-AUC:     {analysis['comm_total_auc']:.4f}")

    print(f"\n  Commutator block structure (mean ||[θ_a, θ_b]||_F):")
    print(f"    M↔M={analysis['comm_AA']:.6f}   "
          f"M↔B={analysis['comm_AB']:.6f}   "
          f"B↔B={analysis['comm_BB']:.6f}")

    print(f"\n  Per-relation detail:")
    print(f"    {'Rel':>4s}  {'Modality':>10s}  {'Name':>28s}  "
          f"{'K_rr':>10s}  {'CommTotal':>10s}")
    for r in range(N_REL):
        if r < N_M:
            mod = relations_M[r].source_modality
            name = relations_M[r].name[:28]
        else:
            mod = "random"
            name = f"B_{r - N_M}"
        print(f"    {r:4d}  {mod:>10s}  {name:>28s}  "
              f"{analysis['K_diag'][r]:10.4f}  "
              f"{analysis['comm_total'][r]:10.4f}")

    state_after_p1 = copy.deepcopy(model.state_dict())
    adam_after_p1 = copy.deepcopy(optimizer.state_dict())

    # ========================================================
    #  PHASE 3: SELECTIVE CYCLE LOSS
    # ========================================================

    print(f"\n{'='*70}")
    print(f"  PHASE 3: Selective Cycle Loss ({PHASE3_EPOCHS} ep, Type M closures,")
    print(f"           K_rr tracked every {PHASE3_LOG_EVERY} ep — feature B)")
    print(f"{'='*70}")

    model.load_state_dict(copy.deepcopy(state_after_p1))
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    optimizer.load_state_dict(adam_after_p1)

    # Chunked Phase 3: run in blocks of PHASE3_LOG_EVERY epochs and snapshot
    # K_rr AUC after each block to map the K_rr decay phenomenon (feature B).
    k_auc_p3_trace = []
    hist_p3 = {'val_acc': [], 'tk': []}
    n_chunks = max(1, PHASE3_EPOCHS // PHASE3_LOG_EVERY)
    remainder = PHASE3_EPOCHS % PHASE3_LOG_EVERY
    analysis_p3 = None

    for chunk_idx in range(n_chunks + (1 if remainder > 0 else 0)):
        chunk_ep = PHASE3_LOG_EVERY if chunk_idx < n_chunks else remainder
        ep_end   = chunk_idx * PHASE3_LOG_EVERY + chunk_ep
        hist_chunk = train_phase(
            model, optimizer,
            train_triples, train_labels, val_triples, val_labels,
            epochs=chunk_ep, batch_size=BATCH_SIZE,
            lambda_k=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
            lambda_cycle=LAMBDA_CYCLE,
            con_cycles_t=con_cycles_t, incon_cycles_t=incon_cycles_t,
            margin=MARGIN,
            grad_clip=GRAD_CLIP, label=f"P3@{ep_end:4d}", log_every=chunk_ep,
        )
        for key in hist_p3:
            hist_p3[key].extend(hist_chunk.get(key, []))
        analysis_p3 = analyze_killing_spectrum(model, type_labels, N_M, N_RANDOM_B)
        k_auc_p3_trace.append({
            'epoch': ep_end,
            'k_auc': analysis_p3['k_diag_auc'],
            'comm_auc': analysis_p3['comm_total_auc'],
        })

    print(f"\n  Phase 3 holonomy evaluation:")
    hol_p3 = evaluate_holonomy(
        model, cycle_data, con_paths, incon_paths,
        gt_con_hol, gt_incon_hol, label="Phase3"
    )

    print(f"\n  Phase 3 Killing spectrum:")
    print(f"    K_rr  M: mean={analysis_p3['K_diag_A'].mean():.4f}  "
          f"std={analysis_p3['K_diag_A'].std():.4f}")
    print(f"    K_rr  B: mean={analysis_p3['K_diag_B'].mean():.4f}  "
          f"std={analysis_p3['K_diag_B'].std():.4f}")
    print(f"    K_rr AUC (P3): {analysis_p3['k_diag_auc']:.4f}")
    print(f"    Comm  M↔M={analysis_p3['comm_AA']:.6f}  "
          f"M↔B={analysis_p3['comm_AB']:.6f}  "
          f"B↔B={analysis_p3['comm_BB']:.6f}")

    # ========================================================
    #  SUMMARY
    # ========================================================

    t_total = time.time() - t0

    print(f"\n{'='*70}")
    print(f"  SUMMARY — Exp6: Killing Classifier on Multimodal Operators")
    print(f"{'='*70}")

    print(f"\n  Multimodal operators: {N_M}")
    print(f"    Sources: "
          f"text={sum(1 for r in relations_M if r.source_modality=='text')}, "
          f"svg={sum(1 for r in relations_M if r.source_modality=='svg')}, "
          f"audio={sum(1 for r in relations_M if r.source_modality=='audio')}")
    print(f"    Mean spectral gap:  "
          f"{np.mean([r.spectral_gap for r in relations_M]):.4f}")

    print(f"\n  Ground truth Killing (Q_M vs Q_B):")
    print(f"    GT K_rr  M: {K_gt_diag[:N_M].mean():.4f}  "
          f"B: {K_gt_diag[N_M:].mean():.4f}  "
          f"AUC: {gt_k_auc:.4f}")

    print(f"\n  Closure mining:")
    print(f"    Mined cycles: {len(raw_closures)}")
    print(f"    Consistent mean hol:   {cycle_data['consistent_holonomy'].mean():.4f}")
    print(f"    Inconsistent mean hol: {cycle_data['inconsistent_holonomy'].mean():.4f}")
    print(f"    GT holonomy sep: {hol_sep_gt:.2f}x")

    print(f"\n  Phase 1 (triple-only, {PHASE1_EPOCHS}ep):")
    print(f"    Val acc:   {hist_p1['val_acc'][-1]:.4f}")
    print(f"    T_K final: {hist_p1['tk'][-1]:.4f}")
    print(f"    path_AUC:  {hol_p1['path_auc']:.4f}")

    print(f"\n  Killing Classification (Phase 1 — KEY RESULT):")
    print(f"    K_rr ROC-AUC:     {analysis['k_diag_auc']:.4f}")
    print(f"    |K_rr| ROC-AUC:   {analysis['k_diag_abs_auc']:.4f}")
    print(f"    CommTotal ROC-AUC: {analysis['comm_total_auc']:.4f}")
    print(f"    NOTE: type B is now adjoint SU({D}) — same subspace as type M")

    print(f"\n  Phase 3 (+ cycle loss, {PHASE3_EPOCHS}ep):")
    print(f"    Val acc:   {hist_p3['val_acc'][-1]:.4f}")
    print(f"    T_K final: {hist_p3['tk'][-1]:.4f}")
    print(f"    path_AUC:  {hol_p3['path_auc']:.4f}")
    print(f"    K_rr AUC:  {analysis_p3['k_diag_auc']:.4f}")

    print(f"\n  K_rr AUC decay trace (Phase 3, every {PHASE3_LOG_EVERY} ep):")
    print(f"    {'Epoch':>6s}  {'K_rr AUC':>10s}  {'CommTotal AUC':>14s}")
    for entry in k_auc_p3_trace:
        print(f"    {entry['epoch']:6d}  {entry['k_auc']:10.4f}  {entry['comm_auc']:14.4f}")

    print(f"\n  Total time: {t_total:.1f}s")

    # ========================================================
    #  COMPARISON TABLE (vs Exp4 reference values)
    # ========================================================

    print(f"\n{'='*70}")
    print(f"  COMPARISON vs Exp4 (synthetic SU(3) operators)")
    print(f"{'='*70}")
    print(f"  Metric                     Exp6-P1    Exp6-P3")
    print(f"  ------------------------   -------    -------")
    print(f"  GT K_rr AUC (M vs B_adj)   {gt_k_auc:.4f}     —")
    print(f"  K_rr AUC (Phase 1)         {analysis['k_diag_auc']:.4f}     —")
    print(f"  K_rr AUC (Phase 3)         —          {analysis_p3['k_diag_auc']:.4f}")
    print(f"  path_AUC (Phase 1)         {hol_p1['path_auc']:.4f}     —")
    print(f"  path_AUC (Phase 3)         —          {hol_p3['path_auc']:.4f}")
    print(f"  CommTotal AUC (P1)         {analysis['comm_total_auc']:.4f}     —")
    print(f"  CommTotal AUC (P3)         —          {analysis_p3['comm_total_auc']:.4f}")
    print(f"  Type B space               adjoint SU({D}) — same subspace as type M")

    return {
        'n_multimodal': N_M,
        'n_random_b': N_RANDOM_B,
        'n_mined_closures': len(raw_closures),
        'gt_k_auc': gt_k_auc,
        'hol_sep_gt': hol_sep_gt,
        'p1_val_acc': hist_p1['val_acc'][-1],
        'p1_tk': hist_p1['tk'][-1],
        'p1_path_auc': hol_p1['path_auc'],
        'p1_k_auc': analysis['k_diag_auc'],
        'p1_comm_auc': analysis['comm_total_auc'],
        'p3_val_acc': hist_p3['val_acc'][-1],
        'p3_path_auc': hol_p3['path_auc'],
        'p3_k_auc': analysis_p3['k_diag_auc'],
        'p3_comm_auc': analysis_p3['comm_total_auc'],
        'k_auc_p3_trace': k_auc_p3_trace,
        'time_s': t_total,
    }


# ============================================================
#  ENTRY POINT
# ============================================================

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=SEED,
                        help="Single seed (ignored if --seeds is provided)")
    parser.add_argument("--seeds", type=int, nargs="+", default=None,
                        help="One or more seeds for multi-seed runs (e.g. --seeds 42 123 777)")
    parser.add_argument("--quick", action="store_true",
                        help="Short run for smoke-testing (200+100 epochs)")
    parser.add_argument("--gap-filter", type=float, default=None, dest="gap_filter",
                        help="Discard operators with spectral_gap < threshold (ablation D)")
    parser.add_argument("--angle-scale", type=float, default=None, dest="angle_scale",
                        help="Override ANGLE_SCALE (default 0.25). Try 0.5, 1.0, 1.5 to boost K_rr signal.")
    args = parser.parse_args()

    seeds = args.seeds if args.seeds is not None else [args.seed]
    all_results = []
    for s in seeds:
        res = run_experiment(seed=s, quick=args.quick,
                             gap_filter=args.gap_filter,
                             angle_scale=args.angle_scale)
        all_results.append(res)

    if len(all_results) > 1:
        print(f"\n{'='*70}")
        print(f"  MULTI-SEED AGGREGATE  ({len(seeds)} seeds: {seeds})")
        print(f"{'='*70}")
        keys = [
            'gt_k_auc', 'p1_k_auc', 'p1_comm_auc',
            'p1_path_auc', 'p3_k_auc', 'p3_comm_auc', 'p3_path_auc',
        ]
        print(f"  {'Metric':<32s}  {'Mean':>8s}  {'Std':>8s}  {'Min':>8s}  {'Max':>8s}")
        print(f"  {'-'*32}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")
        for key in keys:
            vals = [r[key] for r in all_results]
            print(f"  {key:<32s}  {np.mean(vals):8.4f}  {np.std(vals):8.4f}  "
                  f"{min(vals):8.4f}  {max(vals):8.4f}")
        p1_k_vals = [f'{r["p1_k_auc"]:.4f}' for r in all_results]
        p3_c_vals = [f'{r["p3_comm_auc"]:.4f}' for r in all_results]
        print(f"\n  Per-seed K_rr AUC P1:  {p1_k_vals}")
        print(f"  Per-seed CommTotal P3: {p3_c_vals}")
