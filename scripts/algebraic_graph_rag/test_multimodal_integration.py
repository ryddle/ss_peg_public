"""
Integration smoke test for the multimodal encoder pipeline.
============================================================
Verifies that all four modality encoders and the unified MultimodalEncoder:

  1. Produce AlgebraicRelation objects with valid shapes.
  2. Produce O(n_gen) matrices (Q Q^T ≈ I, error < 0.01).
  3. Are compatible with the AlgKGModel from Exp4.
  4. Produce the correct stacked array shape for the Killing analysis.

Run from the project root:
  python algebraic_graph_rag/test_multimodal_integration.py
"""

import sys
import os
import traceback

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'training_ai_test'))


PASS_MARK = "[PASS]"
FAIL_MARK = "[FAIL]"
SKIP_MARK = "[SKIP]"


def check_orthogonality(ops: np.ndarray, tol: float = 0.01) -> float:
    """Return max |Q Q^T - I| over all operators."""
    n = ops.shape[0]
    n_gen = ops.shape[1]
    max_err = 0.0
    for i in range(n):
        Q = ops[i]
        err = np.max(np.abs(Q @ Q.T - np.eye(n_gen)))
        max_err = max(max_err, err)
    return max_err


# ============================================================
#  TEXT
# ============================================================

def test_text():
    from text_encoder import text_to_operators       # noqa: PLC0415

    sample = (
        "Quarks carry colour charge mediated by gluons. "
        "The strong force confines quarks inside hadrons. "
        "Protons and neutrons are examples of baryons composed of three quarks. "
        "Gluons are the force carriers of the strong nuclear force."
    )
    ops = text_to_operators(sample, d=3, angle_scale=0.25, use_spacy=False)
    assert len(ops) > 0, "No operators from text"
    arr = np.stack([o.operator for o in ops])
    assert arr.shape[1:] == (8, 8), f"Wrong shape: {arr.shape}"
    err = check_orthogonality(arr)
    assert err < 0.01, f"Orthogonality error too high: {err:.4f}"
    return len(ops), err


# ============================================================
#  AUDIO
# ============================================================

def test_audio():
    from audio_encoder import audio_to_operators     # noqa: PLC0415

    rng = np.random.default_rng(42)
    sr = 22050
    # Synthetic audio: 3 seconds of multi-frequency signal
    t = np.linspace(0, 3.0, 3 * sr)
    signal = (np.sin(2 * np.pi * 440 * t) +
              0.5 * np.sin(2 * np.pi * 880 * t) +
              0.25 * np.sin(2 * np.pi * 1320 * t) +
              0.1 * rng.standard_normal(len(t)))

    ops = audio_to_operators(signal, sample_rate=sr, d=3, angle_scale=0.25)
    assert len(ops) > 0, "No operators from audio"
    arr = np.stack([o.operator for o in ops])
    assert arr.shape[1:] == (8, 8), f"Wrong shape: {arr.shape}"
    err = check_orthogonality(arr)
    assert err < 0.01, f"Orthogonality error too high: {err:.4f}"
    return len(ops), err


# ============================================================
#  IMAGE
# ============================================================

def test_image():
    from image_encoder import image_to_operators     # noqa: PLC0415

    rng = np.random.default_rng(42)
    # Synthetic 64×64 image with gradient structure
    img = np.zeros((64, 64), dtype=np.float64)
    img[:32, :] = rng.uniform(0.1, 0.4, (32, 64))
    img[32:, :] = rng.uniform(0.6, 0.9, (32, 64))
    # Add a few rectangles
    img[10:20, 10:30] = 0.8
    img[40:55, 35:60] = 0.15

    ops = image_to_operators(img, d=3, angle_scale=0.25,
                             use_regions=False,   # skip skimage if not installed
                             use_patches=True,
                             max_dim=64,
                             patch_size=8,
                             min_edges_per_type=3)
    assert len(ops) > 0, "No operators from image"
    arr = np.stack([o.operator for o in ops])
    assert arr.shape[1:] == (8, 8), f"Wrong shape: {arr.shape}"
    err = check_orthogonality(arr)
    assert err < 0.01, f"Orthogonality error too high: {err:.4f}"
    return len(ops), err


# ============================================================
#  SVG
# ============================================================

SAMPLE_SVG = """<svg xmlns="http://www.w3.org/2000/svg" width="200" height="200">
  <defs>
    <style>.red { fill: red; } .blue { fill: blue; }</style>
  </defs>
  <g id="group1" transform="translate(10, 10)">
    <rect id="r1" class="red" x="0" y="0" width="50" height="50"/>
    <circle id="c1" class="blue" cx="80" cy="25" r="20"/>
    <rect id="r2" class="red" x="120" y="10" width="40" height="30"/>
  </g>
  <g id="group2" transform="translate(10, 100)">
    <rect id="r3" x="0" y="0" width="60" height="40" fill="url(#r1)"/>
    <circle id="c2" cx="90" cy="20" r="15" class="blue"/>
    <line id="l1" x1="0" y1="50" x2="180" y2="50" stroke="black"/>
  </g>
</svg>"""


def test_svg():
    from svg_encoder import svg_to_operators         # noqa: PLC0415

    ops = svg_to_operators(SAMPLE_SVG, d=3, angle_scale=0.25,
                           use_geometric=True, min_edges_per_type=1)
    assert len(ops) > 0, "No operators from SVG"
    arr = np.stack([o.operator for o in ops])
    assert arr.shape[1:] == (8, 8), f"Wrong shape: {arr.shape}"
    err = check_orthogonality(arr)
    assert err < 0.01, f"Orthogonality error too high: {err:.4f}"
    return len(ops), err


# ============================================================
#  UNIFIED ENCODER
# ============================================================

def test_unified():
    from relational_encoder import MultimodalEncoder  # noqa: PLC0415

    rng = np.random.default_rng(0)
    t = np.linspace(0, 1.0, 22050)
    signal = np.sin(2 * np.pi * 440 * t) + 0.1 * rng.standard_normal(len(t))

    enc = MultimodalEncoder(d=3, angle_scale=0.25)
    enc.encode_text(
        "Energy and momentum are conserved in elastic collisions.",
        use_spacy=False,
    )
    enc.encode_audio(signal, sample_rate=22050, use_chroma=True)
    enc.encode_svg(SAMPLE_SVG, use_geometric=True)

    arr, relations = enc.to_kg_operators()

    assert arr.ndim == 3, f"Expected 3D array, got {arr.ndim}D"
    assert arr.shape[1:] == (8, 8), f"Wrong shape: {arr.shape}"
    err = check_orthogonality(arr)
    assert err < 0.01, f"Orthogonality error too high: {err:.4f}"
    assert len(relations) == arr.shape[0]
    return arr.shape[0], err


# ============================================================
#  COMPATIBILITY WITH AlgKGModel
# ============================================================

def test_algkg_compat():
    """Verify that the produced operators tensor can be used inside AlgKGModel."""
    import torch                                       # noqa: PLC0415
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'training_ai_test'))

    # Import the model from exp4
    exp4_path = os.path.join(os.path.dirname(__file__), 'exp4_killing_classifier.py')
    import importlib.util                              # noqa: PLC0415
    spec = importlib.util.spec_from_file_location("exp4", exp4_path)
    exp4 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(exp4)

    from relational_encoder import MultimodalEncoder   # noqa: PLC0415

    enc = MultimodalEncoder(d=3, angle_scale=0.25)
    enc.encode_text(
        "Symmetry breaking gives mass to gauge bosons via the Higgs mechanism.",
        use_spacy=False,
    )
    enc.encode_svg(SAMPLE_SVG)

    arr, _ = enc.to_kg_operators()
    K, n_gen, _ = arr.shape

    # Build minimal AlgKGModel
    n_entities = 20
    model = exp4.AlgKGModel(n_entities=n_entities, n_relations=K, n_gen=n_gen)

    # Inject encoded operators (model uses float64 internally)
    device = next(model.parameters()).device
    ops_tensor = torch.tensor(arr, dtype=torch.float64).to(device)
    model.relation_ops.data = ops_tensor

    # Forward pass with dummy triples (format: [head, relation, tail])
    head_ids = torch.randint(0, n_entities, (10,)).to(device)
    rel_ids  = torch.randint(0, K,          (10,)).to(device)
    tail_ids = torch.randint(0, n_entities, (10,)).to(device)
    triples  = torch.stack([head_ids, rel_ids, tail_ids], dim=1)
    scores = model(triples)

    assert scores.shape[0] == 10, f"Unexpected score shape: {scores.shape}"
    assert not torch.isnan(scores).any(), "NaN in scores"
    return K, n_gen


# ============================================================
#  RUNNER
# ============================================================

def run_all():
    tests = [
        ("text_encoder",    test_text),
        ("audio_encoder",   test_audio),
        ("image_encoder",   test_image),
        ("svg_encoder",     test_svg),
        ("unified_encoder", test_unified),
        ("algkg_compat",    test_algkg_compat),
    ]

    results = {}
    for name, fn in tests:
        try:
            ret = fn()
            results[name] = ("pass", ret)
            print(f"{PASS_MARK} {name:25s}  → {ret}")
        except AssertionError as e:
            results[name] = ("fail", str(e))
            print(f"{FAIL_MARK} {name:25s}  → AssertionError: {e}")
        except ImportError as e:
            results[name] = ("skip", str(e))
            print(f"{SKIP_MARK} {name:25s}  → ImportError (optional dep): {e}")
        except Exception:
            tb = traceback.format_exc()
            results[name] = ("error", tb)
            print(f"{FAIL_MARK} {name:25s}  → EXCEPTION:\n{tb}")

    n_pass = sum(1 for v in results.values() if v[0] == "pass")
    n_skip = sum(1 for v in results.values() if v[0] == "skip")
    n_fail = len(results) - n_pass - n_skip

    print(f"\n{'=' * 50}")
    print(f"Results: {n_pass} passed, {n_skip} skipped, {n_fail} failed")
    if n_fail == 0:
        print("All required tests passed.")
    return n_fail == 0


if __name__ == "__main__":
    success = run_all()
    sys.exit(0 if success else 1)
