"""
Experiment 7 — Wilson Loop Curvature of the KG Connection
==========================================================
Extension of Exp6 measuring the LOCAL CURVATURE of the gauge connection
defined by the learned KG operators.

Research question:
  Under cycle loss (Phase 3), do Type M (multimodal, physics domain)
  operators form a FLATTER gauge connection — lower Wilson loop curvature
  F_{r1,r2} — compared to Type B operators from a DIFFERENT real dataset?

  This is a REAL vs OTHER REAL comparison (not real vs random).
  Type M: physics texts (QCD, SM gauge, spectral graph) + structured SVG + audio
  Type B: non-physics texts (genomics, economics, music theory) + tree SVG + audio
  Both encode real data through the SAME encoder pipeline and SAME angle_scale.

Key insight (from Exp6):
  T_K (Killing metric) is Ad-invariant → K_rr ≡ C·||θ||² for all
  elements of the same simple algebra.  Comparing random vs real (Exp6)
  still confounded by different spectral structures.  Here BOTH sets
  come from real data, so K_rr has essentially NO discriminating power.

  Wilson loop curvature F_{r1,r2} = ||Q_{r1}^T Q_{r2}^T Q_{r1} Q_{r2} - I||_F
  measures the group commutator (parallel transport around an
  infinitesimal plaquette).  Physics-domain operators (cycle-structured)
  should abelianize under cycle loss → lower curvature.  Non-physics
  operators (no cycle structure) should retain higher curvature.

Resolution hierarchy metaphor:
  1. K          — satellite: detects algebra class
  2. T_closure  — aerial:    detects specific algebra
  3. Holonomy   — drone:     which loops close
  4. path_AUC   — microscope: local flatness
  5. Curvature  — below microscope: F_{r1,r2}  ← THIS EXPERIMENT

Protocol:
  Same as Exp6 Phases 0–3, PLUS:
  - After each Phase 3 chunk, compute per-relation mean curvature C_r
  - Track C_r^M mean, C_r^B mean, curve_auc (ROC-AUC of M vs B by C_r)
  - Track T_closure(M), T_closure(B) co-evolution during Phase 3
  - Report alongside K_rr AUC and CommTotal AUC for comparison

New metrics:
  curve_auc       : ROC-AUC of (M vs B) classification by C_r
                    Hypothesis: C_r^M < C_r^B → AUC > 0.5
  closure_M_final : T_closure(θ_M block) at end of Phase 3
  closure_B_final : T_closure(θ_B block) at end of Phase 3
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
from exp4_killing_classifier import (
    AlgKGModel,
    train_phase,
    analyze_killing_spectrum,
    evaluate_holonomy,
)
from relational_encoder import MultimodalEncoder
from topology_to_algebra import operators_to_array, direction_to_operator, _su_basis

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'training_ai_test'))
from killing_utils import killing_metric, closure_tension

sys.stdout.reconfigure(encoding='utf-8')


# ============================================================
#  SAMPLE MULTIMODAL DATA (identical to Exp6 — self-contained)
# ============================================================

SAMPLE_TEXTS = [
    # Text 0: QCD / hadrons
    """
    Quarks are confined within hadrons by the strong nuclear force mediated
    by gluons. The color charge of quarks comes in three types: red, green,
    and blue. Gluons carry color charge and interact with each other as well
    as with quarks. Confinement prevents free quarks from being observed in
    isolation. The strong force is described by quantum chromodynamics, or
    QCD. Asymptotic freedom means the coupling constant decreases at high
    energies. Hadrons include both baryons and mesons as composite particles.
    Baryons consist of three quarks while mesons contain a quark-antiquark
    pair. The proton is the lightest stable baryon with two up and one down
    quark. Neutrons contain two down quarks and one up quark bound together.
    Pions are the lightest mesons and mediate nuclear forces at low energies.
    The mass gap problem asks why the strong force has a short range despite
    gluons being massless in the Lagrangian formulation. Color singlet states
    are the only observable combinations due to confinement constraints.
    Deep inelastic scattering revealed the parton structure inside nucleons.
    The running coupling constant of QCD decreases logarithmically with scale.
    Chiral symmetry breaking generates most of the hadron masses dynamically.
    The QCD vacuum has nontrivial topological structure with instanton effects.
    Lattice QCD provides nonperturbative calculations of hadron properties.
    The confinement mechanism remains one of the deepest unsolved problems.
    """,
    # Text 1: symmetry breaking / SM gauge
    """
    The Standard Model unifies electromagnetic and weak interactions through
    electroweak symmetry breaking. The Higgs boson was discovered at the
    Large Hadron Collider in 2012 confirming the Higgs mechanism. Spontaneous
    symmetry breaking occurs when the vacuum state does not respect the full
    symmetry of the Lagrangian. The W and Z bosons acquire mass through the
    Higgs mechanism while the photon remains massless. Gauge invariance is
    the fundamental principle underlying all interactions in the Standard
    Model. The electroweak symmetry group is SU(2) times U(1) broken to U(1).
    Fermion masses arise from Yukawa couplings to the Higgs field in the
    ground state. The Cabibbo-Kobayashi-Maskawa matrix describes quark flavor
    mixing in weak interactions. CP violation in kaon and B meson systems
    arises from complex phases in the CKM matrix. Neutrino oscillations
    demonstrate that neutrinos have small but nonzero masses. The seesaw
    mechanism provides a natural explanation for the smallness of neutrino
    masses. Grand unified theories attempt to unify all three gauge couplings
    at a high energy scale. Supersymmetry proposes that every fermion has
    a bosonic superpartner and vice versa. The hierarchy problem asks why
    the Higgs mass is so much smaller than the Planck scale. Renormalization
    group equations describe how coupling constants run with energy scale.
    Anomaly cancellation constrains the fermion content of gauge theories.
    The strong CP problem concerns why QCD preserves CP symmetry so well.
    Axions were proposed to solve the strong CP problem dynamically.
    String theory provides a framework that may unify gravity with gauge forces.
    """,
    # Text 2: spectral graph / holonomy
    """
    Graph Laplacians encode the connectivity structure of networks through
    their eigenvalue spectra. The spectral gap determines mixing times and
    expander properties of random walks on graphs. Holonomy measures the
    failure of parallel transport to return to its starting value around
    closed loops. Connections on principal bundles generalize the notion
    of covariant derivatives to curved spaces. Curvature tensors measure
    how much parallel transport around infinitesimal loops deviates from
    identity. Flat connections have zero curvature and trivial holonomy
    for contractible loops. The Atiyah-Singer index theorem relates analytic
    and topological invariants of differential operators on manifolds.
    Spectral clustering partitions graphs by using eigenvectors of the
    Laplacian as coordinates. The Cheeger inequality bounds the spectral
    gap by the isoperimetric ratio of graph cuts. Ramanujan graphs achieve
    optimal spectral gaps making them excellent expanders. Heat kernels on
    graphs connect spectral theory to random walk probabilities over time.
    Graph neural networks learn hierarchical representations through message
    passing along edges. The Weisfeiler-Leman algorithm tests graph isomorphism
    through iterative neighborhood aggregation. Persistent homology tracks
    topological features across multiple scales in a filtration. Discrete
    exterior calculus provides a framework for differential geometry on meshes.
    Algebraic topology studies spaces through chain complexes and homology groups.
    The fundamental group classifies loops up to continuous deformation in a space.
    Covering spaces correspond to subgroups of the fundamental group via lifting.
    Morse theory connects critical points of smooth functions to topology.
    """,
]

SAMPLE_SVG = """<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="550" height="420"
     viewBox="0 0 550 420">
  <defs>
    <style>
      .node-A { fill: #4a90d9; stroke: #2c6fad; stroke-width: 2; }
      .node-B { fill: #e85d4a; stroke: #b33a28; stroke-width: 2; }
      .node-C { fill: #5cb85c; stroke: #3d8f3d; stroke-width: 2; }
      .edge   { stroke: #888; stroke-width: 1.5; fill: none; }
      .bridge { stroke: #bbb; stroke-width: 1; stroke-dasharray: 4 3; fill: none; }
      .label  { font-size: 11px; font-family: sans-serif; fill: #333; }
    </style>
  </defs>

  <!-- Edges (drawn first) -->
  <g id="layer_edges">
    <!-- Cluster A internal -->
    <line class="edge" x1="90"  y1="80"  x2="165" y2="140"/>
    <line class="edge" x1="90"  y1="80"  x2="55"  y2="210"/>
    <line class="edge" x1="165" y1="140" x2="55"  y2="210"/>
    <line class="edge" x1="55"  y1="210" x2="115" y2="340"/>
    <line class="edge" x1="165" y1="140" x2="115" y2="340"/>

    <!-- Cluster B internal -->
    <line class="edge" x1="255" y1="50"  x2="355" y2="110"/>
    <line class="edge" x1="255" y1="50"  x2="275" y2="200"/>
    <line class="edge" x1="355" y1="110" x2="275" y2="200"/>
    <line class="edge" x1="275" y1="200" x2="305" y2="340"/>
    <line class="edge" x1="355" y1="110" x2="305" y2="340"/>

    <!-- Cluster C internal -->
    <line class="edge" x1="435" y1="80"  x2="465" y2="190"/>
    <line class="edge" x1="435" y1="80"  x2="425" y2="330"/>
    <line class="edge" x1="465" y1="190" x2="425" y2="330"/>
    <line class="edge" x1="465" y1="190" x2="414" y2="183"/>
    <line class="edge" x1="425" y1="330" x2="414" y2="183"/>

    <!-- Inter-cluster bridges -->
    <line class="bridge" x1="165" y1="140" x2="255" y2="50"/>
    <line class="bridge" x1="355" y1="110" x2="435" y2="80"/>
    <line class="bridge" x1="115" y1="340" x2="305" y2="340"/>
    <line class="bridge" x1="305" y1="340" x2="425" y2="330"/>
  </g>

  <!-- Nodes -->
  <g id="layer_nodes">
    <!-- Cluster A (blue) -->
    <circle class="node-A" cx="90"  cy="80"  r="14"/>
    <circle class="node-A" cx="165" cy="140" r="14"/>
    <circle class="node-A" cx="55"  cy="210" r="14"/>
    <circle class="node-A" cx="115" cy="340" r="14"/>
    <!-- Cluster B (red) -->
    <circle class="node-B" cx="255" cy="50"  r="14"/>
    <circle class="node-B" cx="355" cy="110" r="14"/>
    <circle class="node-B" cx="275" cy="200" r="14"/>
    <circle class="node-B" cx="305" cy="340" r="14"/>
    <!-- Cluster C (green) -->
    <circle class="node-C" cx="435" cy="80"  r="14"/>
    <circle class="node-C" cx="465" cy="190" r="14"/>
    <circle class="node-C" cx="425" cy="330" r="14"/>
    <circle class="node-C" cx="414" cy="183" r="14"/>
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


# ============================================================
#  SAMPLE DATA FOR TYPE B: DIFFERENT REAL DOMAINS
# ============================================================

SAMPLE_TEXTS_B = [
    # Text B0: Molecular biology / genomics
    """
    DNA stores genetic information as sequences of nucleotide base pairs in
    the double helix structure. Transcription converts DNA into messenger RNA
    through RNA polymerase binding at promoter regions. Translation occurs at
    ribosomes where transfer RNA molecules match codons to amino acids.
    Proteins fold into three-dimensional structures determined by their amino
    acid sequences and chemical interactions. Enzymes catalyze biochemical
    reactions by lowering activation energy barriers in metabolic pathways.
    Gene regulation controls which proteins are expressed in different cell
    types and developmental stages. The epigenome includes chemical
    modifications to DNA and histones that influence gene expression.
    CRISPR-Cas9 enables precise editing of genomic sequences through guided
    RNA mechanisms. Horizontal gene transfer moves genetic material between
    organisms outside normal reproduction processes. Evolution proceeds through
    natural selection acting on heritable variation in populations over time.
    Mutations in DNA arise from copying errors or chemical damage to nucleotide
    bases during replication. Genetic drift causes random changes in allele
    frequencies especially in small isolated populations. Developmental biology
    studies how a single fertilized cell generates complex multicellular
    organisms. Cell signaling pathways transmit information from surface
    receptors to nuclear gene regulatory machinery. The proteome represents
    all proteins expressed by a cell under given environmental conditions.
    Comparative genomics reveals evolutionary relationships and functional
    elements across species genomes. Non-coding RNA molecules perform
    regulatory functions beyond their roles as messenger intermediaries.
    Chromatin structure determines the accessibility of genomic regions to
    transcription factors and regulatory proteins. Metabolic networks integrate
    thousands of enzymatic reactions to sustain cellular life processes.
    Systems biology models biological networks to understand emergent behaviors.
    """,
    # Text B1: Macroeconomics / behavioral economics
    """
    Markets allocate resources through the interaction of supply and demand
    at equilibrium prices determined by buyer and seller decisions. Monetary
    policy controls the money supply and interest rates to stabilize output.
    Fiscal policy uses government spending and taxation to influence aggregate
    demand in economies during expansions and recessions. Behavioral economics
    studies how cognitive biases lead individuals to deviate from rational
    choice predictions. Inflation erodes purchasing power when the money supply
    grows faster than real productive output capacity. Unemployment arises from
    structural mismatches between worker skills and available job opportunities.
    The business cycle describes recurring expansions and contractions in
    aggregate economic activity over medium-term horizons. International trade
    exploits comparative advantages that allow nations to specialize production
    efficiently. Financial markets aggregate dispersed information about future
    economic conditions into observable asset prices. Game theory models
    strategic interactions where outcomes depend on the simultaneous decisions
    of multiple agents. Central banks use open market operations to influence
    short-term interest rates and overall liquidity conditions. The wealth
    effect describes how changes in asset values influence consumer spending
    behavior. Keynesian economics emphasizes the role of aggregate demand in
    determining short-run economic output levels. Credit markets channel
    savings toward productive investments through financial intermediation.
    Risk aversion leads investors to demand higher returns from assets with
    greater uncertainty and volatility. Portfolio diversification reduces
    exposure to idiosyncratic risk by combining uncorrelated asset classes.
    Labor productivity growth drives improvements in living standards over
    long time horizons through technological progress. Institutional economics
    examines how rules and enforcement mechanisms shape economic behavior.
    Market failures occur when prices fail to reflect full social costs of
    transactions, justifying regulatory intervention. Economic growth models
    identify capital accumulation and technological progress as primary drivers.
    """,
    # Text B2: Music theory / harmonic analysis
    """
    Harmony arises from simultaneous pitch combinations that create consonant
    or dissonant sound qualities perceived by listeners. Counterpoint governs
    the combination of independent melodic lines according to voice leading
    rules developed over centuries of practice. Functional tonality organizes
    pitches around a central tonic using hierarchical chord progression
    relationships. Scales are ordered collections of pitches that define the
    tonal vocabulary of a musical passage or composition. Rhythm organizes
    sound in time through patterns of duration stress and metrical grouping.
    The circle of fifths maps relationships between major and minor keys
    through chromatic distance and common tones. Modulation shifts the tonal
    center of a piece to a different key creating contrast and forward motion.
    Chromaticism uses pitches outside the prevailing diatonic scale to create
    expressive harmonic tension and color. Spectral analysis decomposes musical
    timbres into component frequencies and their amplitudes over time.
    The Fourier series represents periodic waveforms as sums of sinusoidal
    harmonics at integer multiples of the fundamental. The overtone series
    determines natural harmonic relationships between simultaneously sounding
    pitches in acoustic instruments. Voice leading minimizes melodic motion
    between chord tones to produce smooth harmonic progressions. Polyphony
    combines multiple independent voices each carrying a complete melodic line
    simultaneously through contrapuntal techniques. Cadences mark phrase
    endings through characteristic harmonic progressions that establish closure.
    Serialism organizes pitch duration dynamics and articulation through
    pre-composed ordering sequences. Modal scales including Dorian Phrygian
    and Lydian provide alternatives to major and minor tonality. Microtonality
    divides the octave into intervals smaller than the Western semitone.
    Form in music emerges from patterns of repetition variation and contrast
    across sections of a composition. Orchestration assigns musical material
    to instruments based on timbral characteristics and register. Temporal
    hierarchies in music span from individual note durations to large sections.
    """,
]

SAMPLE_SVG_B = """<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="550" height="420"
     viewBox="0 0 550 420">
  <defs>
    <style>
      .node-root { fill: #9b59b6; stroke: #7d3c98; stroke-width: 2; }
      .node-mid  { fill: #e67e22; stroke: #ca6f1e; stroke-width: 2; }
      .node-leaf { fill: #27ae60; stroke: #1e8449; stroke-width: 2; }
      .edge      { stroke: #888; stroke-width: 1.5; fill: none; }
      .cross     { stroke: #bbb; stroke-width: 1; stroke-dasharray: 4 3; fill: none; }
      .label     { font-size: 11px; font-family: sans-serif; fill: #333; }
    </style>
  </defs>

  <!-- Hierarchical tree topology (3-level, 12 nodes) -->
  <!-- Different from 3-clique topology of SAMPLE_SVG -->
  <g id="layer_edges">
    <!-- Root to level-1 -->
    <line class="edge" x1="275" y1="50"  x2="110" y2="160"/>
    <line class="edge" x1="275" y1="50"  x2="275" y2="160"/>
    <line class="edge" x1="275" y1="50"  x2="440" y2="160"/>
    <!-- N1 subtree -->
    <line class="edge" x1="110" y1="160" x2="50"  y2="290"/>
    <line class="edge" x1="110" y1="160" x2="170" y2="290"/>
    <!-- N2 subtree -->
    <line class="edge" x1="275" y1="160" x2="210" y2="290"/>
    <line class="edge" x1="275" y1="160" x2="275" y2="290"/>
    <line class="edge" x1="275" y1="160" x2="340" y2="290"/>
    <!-- N3 subtree -->
    <line class="edge" x1="440" y1="160" x2="380" y2="290"/>
    <line class="edge" x1="440" y1="160" x2="500" y2="290"/>
    <!-- Extra depth: N4 has a child -->
    <line class="edge" x1="50"  y1="290" x2="110" y2="370"/>
    <line class="edge" x1="170" y1="290" x2="110" y2="370"/>
    <!-- Cross-links between subtrees -->
    <line class="cross" x1="170" y1="290" x2="210" y2="290"/>
    <line class="cross" x1="340" y1="290" x2="380" y2="290"/>
  </g>

  <g id="layer_nodes">
    <circle class="node-root" cx="275" cy="50"  r="14"/>
    <circle class="node-mid"  cx="110" cy="160" r="14"/>
    <circle class="node-mid"  cx="275" cy="160" r="14"/>
    <circle class="node-mid"  cx="440" cy="160" r="14"/>
    <circle class="node-leaf" cx="50"  cy="290" r="14"/>
    <circle class="node-leaf" cx="170" cy="290" r="14"/>
    <circle class="node-leaf" cx="210" cy="290" r="14"/>
    <circle class="node-leaf" cx="275" cy="290" r="14"/>
    <circle class="node-leaf" cx="340" cy="290" r="14"/>
    <circle class="node-leaf" cx="380" cy="290" r="14"/>
    <circle class="node-leaf" cx="500" cy="290" r="14"/>
    <circle class="node-leaf" cx="110" cy="370" r="14"/>
  </g>

  <g id="layer_labels">
    <text class="label" x="265" y="54">R</text>
    <text class="label" x="100" y="164">N1</text>
    <text class="label" x="265" y="164">N2</text>
    <text class="label" x="430" y="164">N3</text>
    <text class="label" x="40"  y="294">N4</text>
    <text class="label" x="160" y="294">N5</text>
    <text class="label" x="200" y="294">N6</text>
    <text class="label" x="265" y="294">N7</text>
    <text class="label" x="330" y="294">N8</text>
    <text class="label" x="370" y="294">N9</text>
    <text class="label" x="490" y="294">N10</text>
    <text class="label" x="100" y="374">N11</text>
  </g>
</svg>"""


def _make_synthetic_audio(duration_sec: float = 3.0, sr: int = 22050, seed: int = 42):
    """Generate a multi-frequency synthetic audio signal."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, duration_sec, int(duration_sec * sr))
    signal = (
        1.0  * np.sin(2 * np.pi * 220 * t + 0.0) +
        0.7  * np.sin(2 * np.pi * 330 * t + 0.3) +
        0.5  * np.sin(2 * np.pi * 440 * t + 0.7) +
        0.3  * np.sin(2 * np.pi * 660 * t + 1.1) +
        0.15 * np.sin(2 * np.pi * 880 * t + 1.8) +
        0.05 * rng.standard_normal(len(t))
    )
    return signal / np.abs(signal).max(), sr


# (generate_random_adjoint_operators removed — Exp7 uses real Type B data
#  from SAMPLE_TEXTS_B + SAMPLE_SVG_B instead of random operators)


# ============================================================
#  APPROXIMATE CLOSURE MINING  (identical to Exp6)
# ============================================================

def mine_closure_cycles(
    ops: np.ndarray,
    n_cycles: int = 50,
    max_closure_hol: float | None = None,
    seed: int = 42,
) -> list:
    """Mine approximate algebraic closure cycles (auto 25th-percentile threshold)."""
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
                residual = ops[k].T @ product
                hol = np.linalg.norm(residual - I, 'fro')
                cycles.append((i, j, k, hol))

    cycles.sort(key=lambda x: x[3])

    if max_closure_hol is None:
        all_hols = [h for _, _, _, h in cycles]
        max_closure_hol = float(np.percentile(all_hols, 25)) if all_hols else 2.5

    good = [(i, j, k) for i, j, k, h in cycles if h < max_closure_hol]
    if not good:
        good = [(i, j, k) for i, j, k, _ in cycles[:max(n_cycles // 4, 1)]]

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
    """Build consistent/inconsistent cycle tensors from mined closures."""
    rng = np.random.default_rng(seed)
    n = ops.shape[0]
    I = np.eye(ops.shape[1])

    if len(closure_cycles) >= n_consistent:
        idxs = rng.choice(len(closure_cycles), n_consistent, replace=False)
        con = [closure_cycles[i] for i in idxs]
    else:
        con = closure_cycles * (n_consistent // max(len(closure_cycles), 1) + 1)
        con = con[:n_consistent]

    closure_set = set(map(tuple, closure_cycles))
    incon = []
    attempts = 0
    while len(incon) < n_inconsistent and attempts < 50_000:
        tri = tuple(rng.choice(n, 3, replace=False))
        if tri not in closure_set:
            incon.append(list(tri))
        attempts += 1

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
#  WILSON LOOP CURVATURE  (NEW in Exp7)
# ============================================================

def compute_pairwise_curvature(ops: torch.Tensor) -> torch.Tensor:
    """
    Compute the Wilson loop curvature matrix F[r1, r2] for all operator pairs.

    F[r1, r2] = ||Q_{r1}^T @ Q_{r2}^T @ Q_{r1} @ Q_{r2} - I||_F

    This is the Frobenius norm of the group commutator minus identity.
    For abelianized operators: Q_{r1} Q_{r2} ≈ Q_{r2} Q_{r1} → F → 0.
    For generic operators in a non-abelian group: F > 0.

    Args:
        ops: (n_rel, n_gen, n_gen) tensor of orthogonal transport operators.

    Returns:
        F: (n_rel, n_rel) curvature matrix (symmetric, zero diagonal).
    """
    n, ng, _ = ops.shape
    I = torch.eye(ng, dtype=ops.dtype, device=ops.device)
    ops_T = ops.transpose(-1, -2)  # (n, ng, ng)
    F = torch.zeros(n, n, dtype=ops.dtype, device=ops.device)

    for r1 in range(n):
        # Q_{r1}^T @ Q_{r2}^T @ Q_{r1} @ Q_{r2}  for all r2 simultaneously
        # step1: Q_{r1}^T @ Q_{r2}^T  shape (n, ng, ng)
        step1 = ops_T[r1].unsqueeze(0) @ ops_T          # (n, ng, ng)
        # step2: step1 @ Q_{r1}
        step2 = step1 @ ops[r1].unsqueeze(0)             # (n, ng, ng)
        # step3: step2 @ Q_{r2}
        step3 = step2 @ ops                               # (n, ng, ng)
        diff = step3 - I.unsqueeze(0)                     # (n, ng, ng)
        F[r1] = diff.norm(dim=(-2, -1))                   # (n,)

    return F


def per_relation_curvature(
    ops: torch.Tensor,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute per-relation mean curvature C_r = mean_{r2 != r1} F[r1, r2].

    Returns:
        C: (n_rel,) mean curvature per relation (numpy).
        F_full: (n_rel, n_rel) full curvature matrix (numpy).
    """
    with torch.no_grad():
        F = compute_pairwise_curvature(ops)
    n = ops.shape[0]
    mask = ~torch.eye(n, dtype=torch.bool, device=ops.device)
    denom = mask.float().sum(dim=1).clamp(min=1.0)
    C = (F * mask.float()).sum(dim=1) / denom
    return C.cpu().numpy(), F.cpu().numpy()


def analyze_curvature(
    model: 'AlgKGModel',
    type_labels: np.ndarray,
    n_M: int,
) -> dict:
    """
    Compute Wilson loop curvature metrics and AUC for M vs B classification.

    Returns dict with:
      curve_C      : (n_rel,) per-relation mean curvature
      curve_M_mean : mean C_r for Type M operators
      curve_B_mean : mean C_r for Type B operators
      curve_auc    : ROC-AUC of (M vs B) by C_r
                     (hypothesis: C_r^M < C_r^B → AUC < 0.5 if label=0 for M)
                     We report max(auc, 1-auc) for symmetry.
      curve_auc_raw: raw ROC-AUC score
      F_matrix     : (n_rel, n_rel) full curvature matrix
      closure_M    : T_closure tension for Type M operators
      closure_B    : T_closure tension for Type B operators
    """
    ops = model.get_relation_matrices()  # (n_rel, n_gen, n_gen) — nn.Parameter
    ops_t = ops.detach().to(dtype=torch.float64, device=DEVICE)

    C, F_matrix = per_relation_curvature(ops_t)

    # Closure tension per block
    ops_M_t = ops_t[:n_M]
    ops_B_t = ops_t[n_M:]
    with torch.no_grad():
        cl_M = float(closure_tension(ops_M_t).item()) if n_M >= 2 else float('nan')
        cl_B = float(closure_tension(ops_B_t).item()) if ops_B_t.shape[0] >= 2 else float('nan')

    try:
        curve_auc_raw = roc_auc_score(type_labels, C)
        # Normalize so AUC > 0.5 means "M has lower curvature than B"
        curve_auc = max(curve_auc_raw, 1.0 - curve_auc_raw)
    except ValueError:
        curve_auc_raw = float('nan')
        curve_auc = float('nan')

    return {
        'curve_C': C,
        'curve_M_mean': float(C[:n_M].mean()),
        'curve_B_mean': float(C[n_M:].mean()),
        'curve_auc': curve_auc,
        'curve_auc_raw': curve_auc_raw,
        'F_matrix': F_matrix,
        'closure_M': cl_M,
        'closure_B': cl_B,
    }


# ============================================================
#  CONFIG
# ============================================================

D = 3
ANGLE_SCALE = 0.5
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

MIN_GAP = 0.0
PHASE3_LOG_EVERY = 50
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
    print("  Exp7 — Wilson Loop Curvature of the KG Connection")
    print("=" * 70)
    print(f"  Device : {DEVICE}")
    print(f"  Seed   : {seed}")
    print(f"  d      : {D}  (n_gen = {D*D-1})")
    print(f"  Hypothesis: M operators form flatter gauge connection than B")

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

    if _gap > 0:
        before = enc.n_relations
        enc._relations = [r for r in enc._relations if r.spectral_gap >= _gap]
        removed = before - enc.n_relations
        if removed:
            print(f"  [gap filter] Removed {removed} operators with gap < {_gap} "
                  f"({enc.n_relations} remaining)")

    if enc.n_relations < 2:
        raise RuntimeError("Not enough multimodal operators extracted.")

    Q_M_arr, relations_M = enc.to_kg_operators()
    N_M = Q_M_arr.shape[0]
    print(f"\n  Total multimodal operators : {N_M}")
    print(f"  Source modalities : {set(r.source_modality for r in relations_M)}")
    print(f"  Relation names    :")
    for rel in relations_M:
        print(f"    [{rel.source_modality:8s}]  {rel.name:30s}  "
              f"gap={rel.spectral_gap:.4f}  edges={rel.n_edges_in_subgraph}")

    # ========================================================
    #  PHASE 0B: TYPE B OPERATORS (DIFFERENT REAL DOMAIN)
    # ========================================================

    print(f"\n--- Extracting Type B operators (different real domain) ---")
    enc_B = MultimodalEncoder(d=D, angle_scale=ANGLE_SCALE, seed=seed + 1000)

    for i, txt in enumerate(SAMPLE_TEXTS_B):
        enc_B.encode_text(txt, use_spacy=False, window=4, min_edges_per_type=2)
        print(f"  text_B[{i}]: +{enc_B.n_relations} relations so far")

    n_after_text_B = enc_B.n_relations
    enc_B.encode_svg(SAMPLE_SVG_B, use_geometric=True, min_edges_per_type=1)
    print(f"  svg_B:      +{enc_B.n_relations - n_after_text_B} new relations "
          f"({enc_B.n_relations} total)")

    n_after_svg_B = enc_B.n_relations
    signal_B, sr_B = _make_synthetic_audio(duration_sec=3.0, sr=SR, seed=seed + 1000)
    enc_B.encode_audio(signal_B, sr_B, n_fft=512, hop_length=256,
                       n_freq_bands=16, max_frames=64,
                       use_chroma=True, min_edges_per_type=3)
    print(f"  audio_B:    +{enc_B.n_relations - n_after_svg_B} new relations "
          f"({enc_B.n_relations} total)")

    if _gap > 0:
        enc_B._relations = [r for r in enc_B._relations if r.spectral_gap >= _gap]

    if enc_B.n_relations < 2:
        raise RuntimeError("Not enough Type B operators extracted.")

    Q_B_arr, relations_B = enc_B.to_kg_operators()
    N_B = Q_B_arr.shape[0]
    print(f"\n  Total Type B operators : {N_B}")
    print(f"  Source modalities B : {set(r.source_modality for r in relations_B)}")
    print(f"  Relation names B  :")
    for rel in relations_B:
        print(f"    [{rel.source_modality:8s}]  {rel.name:30s}  "
              f"gap={rel.spectral_gap:.4f}  edges={rel.n_edges_in_subgraph}")

    I_np = np.eye(n_gen)
    dist_M = np.mean([np.linalg.norm(Q_M_arr[r] - I_np, 'fro') for r in range(N_M)])
    dist_B = np.mean([np.linalg.norm(Q_B_arr[r] - I_np, 'fro') for r in range(N_B)])
    print(f"  ||Q - I|| mean: Type M={dist_M:.4f}, Type B={dist_B:.4f}")

    Q_all = np.concatenate([Q_M_arr, Q_B_arr], axis=0)
    N_REL = Q_all.shape[0]
    type_labels = np.array([0] * N_M + [1] * N_B)
    print(f"  Total relations: {N_REL}  (M={N_M}, B={N_B})")
    print(f"  Comparison: physics-domain ops vs non-physics-domain ops (real vs real)")

    # ========================================================
    #  PHASE 0B2: GROUND-TRUTH CURVATURE ANALYSIS (before training)
    # ========================================================

    print("\n--- Ground-truth Wilson loop curvature analysis ---")
    Q_all_t = torch.tensor(Q_all, dtype=torch.float64, device=DEVICE)
    C_gt, F_gt = per_relation_curvature(Q_all_t)
    print(f"  GT curvature  Type M: mean={C_gt[:N_M].mean():.4f}, "
          f"std={C_gt[:N_M].std():.4f}")
    print(f"  GT curvature  Type B: mean={C_gt[N_M:].mean():.4f}, "
          f"std={C_gt[N_M:].std():.4f}")
    try:
        gt_curve_auc_raw = roc_auc_score(type_labels, C_gt)
        gt_curve_auc = max(gt_curve_auc_raw, 1.0 - gt_curve_auc_raw)
    except ValueError:
        gt_curve_auc = float('nan')
    print(f"  GT curvature AUC (M vs B): {gt_curve_auc:.4f}")

    # Ground-truth closure tension
    Q_M_t = torch.tensor(Q_M_arr, dtype=torch.float64, device=DEVICE)
    Q_B_t = torch.tensor(Q_B_arr, dtype=torch.float64, device=DEVICE)
    with torch.no_grad():
        gt_cl_M = float(closure_tension(Q_M_t).item()) if N_M >= 2 else float('nan')
        gt_cl_B = float(closure_tension(Q_B_t).item()) if Q_B_t.shape[0] >= 2 else float('nan')
    print(f"  GT T_closure  Type M: {gt_cl_M:.6f}")
    print(f"  GT T_closure  Type B: {gt_cl_B:.6f}")

    # Ground-truth Killing for comparison
    with torch.no_grad():
        K_gt = killing_metric(Q_all_t).cpu().numpy()
    K_gt_diag = np.diag(K_gt)
    try:
        gt_k_auc = roc_auc_score(type_labels, K_gt_diag)
    except ValueError:
        gt_k_auc = float('nan')
    print(f"  GT K_rr AUC (M vs B):     {gt_k_auc:.4f}  (expected ~0.27, Ad-invariant)")

    # ========================================================
    #  PHASE 0C: MINE CLOSURE CYCLES
    # ========================================================

    print("\n--- Mining approximate closure cycles from Type M operators ---")
    raw_closures = mine_closure_cycles(
        Q_M_arr,
        n_cycles=300,
        max_closure_hol=None,
        seed=seed,
    )
    print(f"  Raw mined closure candidates: {len(raw_closures)}")

    hols_mined = []
    for (i, j, k) in raw_closures[:10]:
        H = Q_M_arr[k] @ Q_M_arr[j] @ Q_M_arr[i]
        hol = np.linalg.norm(H - I_np, 'fro')
        hols_mined.append(hol)
    if hols_mined:
        print(f"  Best 10 closure holonomies: min={min(hols_mined):.4f}, "
              f"max={max(hols_mined):.4f}, mean={np.mean(hols_mined):.4f}")

    rng0 = np.random.default_rng(seed)
    rand_hols = []
    for _ in range(50):
        i, j, k = rng0.choice(N_M, 3, replace=False)
        H = Q_M_arr[k] @ Q_M_arr[j] @ Q_M_arr[i]
        rand_hols.append(np.linalg.norm(H - I_np, 'fro'))
    print(f"  Random triple holonomies (baseline): mean={np.mean(rand_hols):.4f}")

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
    #  PHASE 0D: SYNTHETIC KG GENERATION
    # ========================================================

    print("\n--- Generating synthetic KG ---")
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

    print(f"\n  Phase 1 holonomy evaluation:")
    hol_p1 = evaluate_holonomy(
        model, cycle_data, con_paths, incon_paths,
        gt_con_hol, gt_incon_hol, label="Phase1"
    )

    # ========================================================
    #  PHASE 2: KILLING + CURVATURE ANALYSIS (after Phase 1)
    # ========================================================

    print(f"\n{'='*70}")
    print(f"  PHASE 2: Killing + Curvature Analysis (after Phase 1)")
    print(f"{'='*70}")

    analysis_p1 = analyze_killing_spectrum(model, type_labels, N_M, N_B)
    curve_p1 = analyze_curvature(model, type_labels, N_M)

    print(f"\n  Phase 1 Killing classification:")
    print(f"    K_rr ROC-AUC:     {analysis_p1['k_diag_auc']:.4f}")
    print(f"    CommTotal ROC-AUC: {analysis_p1['comm_total_auc']:.4f}")

    print(f"\n  Phase 1 Wilson loop curvature classification:")
    print(f"    C_r^M mean: {curve_p1['curve_M_mean']:.4f}")
    print(f"    C_r^B mean: {curve_p1['curve_B_mean']:.4f}")
    print(f"    Curvature AUC:    {curve_p1['curve_auc']:.4f}  (normalized, > 0.5 = separable)")
    print(f"    Curvature AUC raw:{curve_p1['curve_auc_raw']:.4f}")
    print(f"    T_closure M: {curve_p1['closure_M']:.6f}")
    print(f"    T_closure B: {curve_p1['closure_B']:.6f}")

    state_after_p1 = copy.deepcopy(model.state_dict())
    adam_after_p1 = copy.deepcopy(optimizer.state_dict())

    # ========================================================
    #  PHASE 3: SELECTIVE CYCLE LOSS + CURVATURE TRACKING
    # ========================================================

    print(f"\n{'='*70}")
    print(f"  PHASE 3: Selective Cycle Loss ({PHASE3_EPOCHS} ep)")
    print(f"  Tracking K_rr, CommTotal, Wilson curvature, T_closure every {PHASE3_LOG_EVERY} ep")
    print(f"{'='*70}")

    model.load_state_dict(copy.deepcopy(state_after_p1))
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    optimizer.load_state_dict(adam_after_p1)

    p3_trace = []   # per-chunk metrics
    hist_p3 = {'val_acc': [], 'tk': []}
    n_chunks = max(1, PHASE3_EPOCHS // PHASE3_LOG_EVERY)
    remainder = PHASE3_EPOCHS % PHASE3_LOG_EVERY
    analysis_p3 = None
    curve_p3 = None

    for chunk_idx in range(n_chunks + (1 if remainder > 0 else 0)):
        chunk_ep = PHASE3_LOG_EVERY if chunk_idx < n_chunks else remainder
        if chunk_ep == 0:
            continue
        ep_end = chunk_idx * PHASE3_LOG_EVERY + chunk_ep

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

        analysis_p3 = analyze_killing_spectrum(model, type_labels, N_M, N_B)
        curve_p3 = analyze_curvature(model, type_labels, N_M)

        p3_trace.append({
            'epoch': ep_end,
            'k_auc': analysis_p3['k_diag_auc'],
            'comm_auc': analysis_p3['comm_total_auc'],
            'curve_M': curve_p3['curve_M_mean'],
            'curve_B': curve_p3['curve_B_mean'],
            'curve_auc': curve_p3['curve_auc'],
            'curve_auc_raw': curve_p3['curve_auc_raw'],
            'closure_M': curve_p3['closure_M'],
            'closure_B': curve_p3['closure_B'],
        })

    print(f"\n  Phase 3 holonomy evaluation:")
    hol_p3 = evaluate_holonomy(
        model, cycle_data, con_paths, incon_paths,
        gt_con_hol, gt_incon_hol, label="Phase3"
    )

    # ========================================================
    #  SUMMARY
    # ========================================================

    t_total = time.time() - t0

    print(f"\n{'='*70}")
    print(f"  SUMMARY — Exp7: Wilson Loop Curvature of the KG Connection")
    print(f"{'='*70}")

    print(f"\n  Type M operators (physics domain): {N_M}")
    print(f"    Sources M: "
          f"text={sum(1 for r in relations_M if r.source_modality=='text')}, "
          f"svg={sum(1 for r in relations_M if r.source_modality=='svg')}, "
          f"audio={sum(1 for r in relations_M if r.source_modality=='audio')}")
    print(f"  Type B operators (non-physics domain): {N_B}")
    print(f"    Sources B: "
          f"text={sum(1 for r in relations_B if r.source_modality=='text')}, "
          f"svg={sum(1 for r in relations_B if r.source_modality=='svg')}, "
          f"audio={sum(1 for r in relations_B if r.source_modality=='audio')}")

    print(f"\n  Ground-truth metrics (before training):")
    print(f"    GT K_rr AUC:        {gt_k_auc:.4f}  (Ad-invariant → expected ~0.27)")
    print(f"    GT curvature AUC:   {gt_curve_auc:.4f}  (pre-training baseline)")
    print(f"    GT T_closure M:     {gt_cl_M:.6f}")
    print(f"    GT T_closure B:     {gt_cl_B:.6f}")

    print(f"\n  Phase 1 (triple-only, {PHASE1_EPOCHS}ep):")
    print(f"    Val acc:        {hist_p1['val_acc'][-1]:.4f}")
    print(f"    T_K final:      {hist_p1['tk'][-1]:.4f}")
    print(f"    path_AUC:       {hol_p1['path_auc']:.4f}")
    print(f"    K_rr AUC:       {analysis_p1['k_diag_auc']:.4f}")
    print(f"    CommTotal AUC:  {analysis_p1['comm_total_auc']:.4f}")
    print(f"    Curvature AUC:  {curve_p1['curve_auc']:.4f}  (norm)  "
          f"{curve_p1['curve_auc_raw']:.4f}  (raw)")
    print(f"    T_closure M:    {curve_p1['closure_M']:.6f}")
    print(f"    T_closure B:    {curve_p1['closure_B']:.6f}")

    print(f"\n  Phase 3 (+ cycle loss, {PHASE3_EPOCHS}ep):")
    print(f"    Val acc:        {hist_p3['val_acc'][-1]:.4f}")
    print(f"    T_K final:      {hist_p3['tk'][-1]:.4f}")
    print(f"    path_AUC:       {hol_p3['path_auc']:.4f}")
    print(f"    K_rr AUC:       {analysis_p3['k_diag_auc']:.4f}")
    print(f"    CommTotal AUC:  {analysis_p3['comm_total_auc']:.4f}")
    print(f"    Curvature AUC:  {curve_p3['curve_auc']:.4f}  (norm)  "
          f"{curve_p3['curve_auc_raw']:.4f}  (raw)")
    print(f"    T_closure M:    {curve_p3['closure_M']:.6f}")
    print(f"    T_closure B:    {curve_p3['closure_B']:.6f}")
    print(f"    C_r^M mean:     {curve_p3['curve_M_mean']:.4f}")
    print(f"    C_r^B mean:     {curve_p3['curve_B_mean']:.4f}")

    print(f"\n  Phase 3 co-evolution trace (every {PHASE3_LOG_EVERY} ep):")
    hdr = f"  {'Epoch':>6s}  {'K_rr':>8s}  {'CommAUC':>8s}  "
    hdr += f"{'CurveAUC':>9s}  {'C_r^M':>8s}  {'C_r^B':>8s}  "
    hdr += f"{'T_cl_M':>9s}  {'T_cl_B':>9s}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for e in p3_trace:
        print(f"  {e['epoch']:6d}  {e['k_auc']:8.4f}  {e['comm_auc']:8.4f}  "
              f"{e['curve_auc']:9.4f}  {e['curve_M']:8.4f}  {e['curve_B']:8.4f}  "
              f"{e['closure_M']:9.6f}  {e['closure_B']:9.6f}")

    print(f"\n  KEY RESULT — Curvature domain separation:")
    delta_curve = curve_p3['curve_B_mean'] - curve_p3['curve_M_mean']
    lower_domain = "physics (M)" if delta_curve > 0 else "non-physics (B)"
    print(f"    Delta C_r (B - M) at end of Phase 3: {delta_curve:+.4f}")
    print(f"    Lower-curvature domain: {lower_domain}")
    if curve_p3['curve_auc'] > 0.65:
        verdict = (f"SEPARABLE  (AUC={curve_p3['curve_auc']:.4f} > 0.65: "
                   f"cycle loss creates domain-specific curvature, "
                   f"{lower_domain} domain abelianizes more)")
    elif curve_p3['curve_auc'] > 0.55:
        verdict = (f"WEAK SIGNAL  (AUC={curve_p3['curve_auc']:.4f}: "
                   f"marginal curvature difference, {lower_domain} domain slightly lower)")
    else:
        verdict = "NOT SEPARABLE  (AUC ≤ 0.55: curvature does not distinguish domains)"
    print(f"    Verdict: {verdict}")
    print(f"\n    Curvature AUC (P3):  {curve_p3['curve_auc']:.4f}")
    print(f"    CommTotal AUC (P3):  {analysis_p3['comm_total_auc']:.4f}  (reference)")
    print(f"    K_rr AUC (P3):       {analysis_p3['k_diag_auc']:.4f}  (expected unstable)")

    print(f"\n  Closure tension co-evolution:")
    if len(p3_trace) >= 2:
        delta_cl_M = curve_p3['closure_M'] - curve_p1['closure_M']
        delta_cl_B = curve_p3['closure_B'] - curve_p1['closure_B']
        print(f"    T_closure M: {curve_p1['closure_M']:.6f} → {curve_p3['closure_M']:.6f}  "
              f"(Δ={delta_cl_M:+.6f})")
        print(f"    T_closure B: {curve_p1['closure_B']:.6f} → {curve_p3['closure_B']:.6f}  "
              f"(Δ={delta_cl_B:+.6f})")
        if delta_cl_M < 0 and abs(delta_cl_M) > abs(delta_cl_B):
            print(f"    Cycle loss reduced T_closure(M) more than T_closure(B) — as expected")

    print(f"\n  Total time: {t_total:.1f}s")

    return {
        'n_multimodal': N_M,
        'n_b': N_B,
        'n_mined_closures': len(raw_closures),
        'gt_k_auc': gt_k_auc,
        'gt_curve_auc': gt_curve_auc,
        'gt_cl_M': gt_cl_M,
        'gt_cl_B': gt_cl_B,
        'hol_sep_gt': hol_sep_gt,
        'p1_val_acc': hist_p1['val_acc'][-1],
        'p1_tk': hist_p1['tk'][-1],
        'p1_path_auc': hol_p1['path_auc'],
        'p1_k_auc': analysis_p1['k_diag_auc'],
        'p1_comm_auc': analysis_p1['comm_total_auc'],
        'p1_curve_auc': curve_p1['curve_auc'],
        'p1_closure_M': curve_p1['closure_M'],
        'p1_closure_B': curve_p1['closure_B'],
        'p3_val_acc': hist_p3['val_acc'][-1],
        'p3_path_auc': hol_p3['path_auc'],
        'p3_k_auc': analysis_p3['k_diag_auc'],
        'p3_comm_auc': analysis_p3['comm_total_auc'],
        'p3_curve_auc': curve_p3['curve_auc'],
        'p3_curve_auc_raw': curve_p3['curve_auc_raw'],
        'p3_curve_M': curve_p3['curve_M_mean'],
        'p3_curve_B': curve_p3['curve_B_mean'],
        'p3_closure_M': curve_p3['closure_M'],
        'p3_closure_B': curve_p3['closure_B'],
        'p3_trace': p3_trace,
        'time_s': t_total,
    }


# ============================================================
#  ENTRY POINT
# ============================================================

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Exp7: Wilson Loop Curvature of the KG Connection")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--seeds", type=int, nargs="+", default=None,
                        help="Multiple seeds (e.g. --seeds 42 123 777)")
    parser.add_argument("--quick", action="store_true",
                        help="Short run for smoke-testing (200+100 epochs)")
    parser.add_argument("--gap-filter", type=float, default=None, dest="gap_filter")
    parser.add_argument("--angle-scale", type=float, default=None, dest="angle_scale")
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
            'gt_k_auc', 'gt_curve_auc',
            'p1_k_auc', 'p1_comm_auc', 'p1_curve_auc',
            'p1_path_auc',
            'p3_k_auc', 'p3_comm_auc', 'p3_curve_auc',
            'p3_path_auc',
            'p3_curve_M', 'p3_curve_B',
            'p3_closure_M', 'p3_closure_B',
            'n_multimodal', 'n_b',
        ]
        print(f"  {'Metric':<32s}  {'Mean':>8s}  {'Std':>8s}  {'Min':>8s}  {'Max':>8s}")
        print(f"  {'-'*32}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")
        for key in keys:
            vals = [r[key] for r in all_results if not np.isnan(r.get(key, float('nan')))]
            if not vals:
                print(f"  {key:<32s}  {'nan':>8s}")
                continue
            print(f"  {key:<32s}  {np.mean(vals):8.4f}  {np.std(vals):8.4f}  "
                  f"{min(vals):8.4f}  {max(vals):8.4f}")

        p3_curve_vals = [f'{r["p3_curve_auc"]:.4f}' for r in all_results]
        p3_comm_vals  = [f'{r["p3_comm_auc"]:.4f}'  for r in all_results]
        p3_k_vals     = [f'{r["p3_k_auc"]:.4f}'     for r in all_results]
        print(f"\n  Per-seed curvature AUC P3: {p3_curve_vals}")
        print(f"  Per-seed CommTotal AUC P3: {p3_comm_vals}")
        print(f"  Per-seed K_rr AUC P3:      {p3_k_vals}")

        # Aggregate verdict — direction determined by data
        mean_curve_auc = np.mean([r['p3_curve_auc'] for r in all_results])
        mean_delta = np.mean([r['p3_curve_B'] - r['p3_curve_M'] for r in all_results])
        lower_domain = "physics (M)" if mean_delta > 0 else "non-physics (B)"
        if mean_curve_auc > 0.65:
            verdict = f"SEPARABLE — {lower_domain} domain abelianizes more"
        elif mean_curve_auc > 0.55:
            verdict = f"WEAK SIGNAL — {lower_domain} domain marginally lower curvature"
        else:
            verdict = "NOT SEPARABLE — curvature does not distinguish domains"
        print(f"\n  AGGREGATE VERDICT: {verdict}  "
              f"(mean curvature AUC = {mean_curve_auc:.4f}, "
              f"mean delta C_r (B-M) = {mean_delta:+.4f})")
