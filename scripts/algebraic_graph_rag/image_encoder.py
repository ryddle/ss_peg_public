"""
Image → Relational Graph → Algebraic Operators
================================================
Converts a raster image (numpy array or PIL Image) into a RelationalGraph
and then into SU(d) adjoint operators for the Killing classifier.

Three graph scales are extracted independently (multi-scale → more operators):

  Scale 1 — Geometric  (pixel/edge level)
    Edge type: "boundary"   — edges between pixels on opposite sides of a gradient
    Edge type: "contour"    — edges along detected contours (Canny)

  Scale 2 — Regional  (superpixel level, ~100 regions)
    Edge type: "region_adj"    — region adjacency graph (shared boundary)
    Edge type: "color_similar" — regions with similar mean colour

  Scale 3 — Semantic  (patch / embedding level)
    Edge type: "patch_sim"    — cosine similarity between 8×8 patch embeddings
    Edge type: "spatial_near" — spatially adjacent patches

Dependencies (all optional):
  - PIL / Pillow   : loading images
  - scikit-image   : Canny, SLIC superpixels, region props
  (Without these, falls back to raw numpy gradient edges.)
"""

import sys
import os
from typing import List, Optional, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from topology_to_algebra import RelationalGraph, AlgebraicRelation, graph_to_operators


# ============================================================
#  IMAGE LOADING / PREPROCESSING
# ============================================================

def _load_image(img) -> np.ndarray:
    """
    Accept numpy array (H, W) or (H, W, C), PIL Image, or file path.
    Returns (H, W) float64 grayscale in [0, 1].
    """
    # File path
    if isinstance(img, str):
        try:
            from PIL import Image as PILImage     # noqa: PLC0415
            img = np.array(PILImage.open(img).convert("L"), dtype=np.float64) / 255.0
            return img
        except ImportError:
            raise RuntimeError("PIL required for loading image files: pip install Pillow")

    # PIL Image
    try:
        from PIL import Image as PILImage         # noqa: PLC0415
        if isinstance(img, PILImage.Image):
            img = np.array(img.convert("L"), dtype=np.float64) / 255.0
            return img
    except ImportError:
        pass

    # Numpy array
    arr = np.asarray(img, dtype=np.float64)
    if arr.ndim == 3:
        # RGB/RGBA → grayscale (luminance weights)
        if arr.shape[2] == 4:
            arr = arr[:, :, :3]
        arr = 0.2126 * arr[:, :, 0] + 0.7152 * arr[:, :, 1] + 0.0722 * arr[:, :, 2]
    if arr.max() > 1.0:
        arr = arr / 255.0
    return arr


def _resize_if_needed(gray: np.ndarray, max_dim: int = 256) -> np.ndarray:
    """Downsample image if larger than max_dim on any side."""
    h, w = gray.shape
    if max(h, w) <= max_dim:
        return gray
    scale = max_dim / max(h, w)
    new_h, new_w = int(h * scale), int(w * scale)
    # Simple block average downsampling
    from skimage.transform import resize as sk_resize         # noqa: PLC0415
    return sk_resize(gray, (new_h, new_w), anti_aliasing=True)


# ============================================================
#  SCALE 1: GEOMETRIC GRAPH
# ============================================================

def _build_geometric_graph(gray: np.ndarray) -> List[Tuple[int, int, str, float]]:
    """
    Build boundary and contour edges from pixel gradient structure.

    Nodes: downsampled pixel grid (at most 64×64 = 4096 nodes)
    Edges:
      "boundary"  — pairs of 4-connected pixels with large gradient
      "contour"   — Canny edge pixels connected to their 8-neighbours
    """
    # Downsample to ≤ 64×64 for tractability
    h, w = gray.shape
    if max(h, w) > 64:
        scale = 64 / max(h, w)
        new_h, new_w = max(1, int(h * scale)), max(1, int(w * scale))
        try:
            from skimage.transform import resize as sk_resize
            gray = sk_resize(gray, (new_h, new_w), anti_aliasing=True)
        except ImportError:
            # Crude stride-based resize
            sh, sw = max(1, h // new_h), max(1, w // new_w)
            gray = gray[::sh, ::sw][:new_h, :new_w]
        h, w = gray.shape

    n_nodes = h * w
    edges = []

    def idx(r, c):
        return r * w + c

    # Horizontal gradient → boundary edges (left-right pairs)
    grad_h = np.abs(np.diff(gray, axis=1))
    for r in range(h):
        for c in range(w - 1):
            weight = float(grad_h[r, c])
            if weight > 0.05:
                edges.append((idx(r, c), idx(r, c + 1), "boundary", weight))

    # Vertical gradient → boundary edges (top-bottom pairs)
    grad_v = np.abs(np.diff(gray, axis=0))
    for r in range(h - 1):
        for c in range(w):
            weight = float(grad_v[r, c])
            if weight > 0.05:
                edges.append((idx(r, c), idx(r + 1, c), "boundary", weight))

    # Canny contour edges (optional)
    try:
        from skimage.feature import canny                     # noqa: PLC0415
        edges_mask = canny(gray, sigma=1.0)
        contour_pixels = list(zip(*np.where(edges_mask)))
        for (r, c) in contour_pixels:
            for dr, dc in [(-1, 0), (0, -1), (1, 0), (0, 1),
                           (-1, -1), (-1, 1), (1, -1), (1, 1)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < h and 0 <= nc < w and edges_mask[nr, nc]:
                    edges.append((idx(r, c), idx(nr, nc), "contour", 1.0))
    except ImportError:
        pass

    nodes = [f"px_{r}_{c}" for r in range(h) for c in range(w)]
    return nodes, edges, n_nodes


# ============================================================
#  SCALE 2: REGION ADJACENCY GRAPH
# ============================================================

def _build_region_graph(
    gray: np.ndarray,
    n_segments: int = 100,
) -> Tuple[List[str], List[Tuple[int, int, str, float]]]:
    """
    Build region adjacency graph using SLIC superpixels (requires scikit-image).

    Edge types:
      "region_adj"    — regions sharing a boundary
      "color_similar" — regions with |mean_color_diff| < 0.1
    """
    try:
        from skimage.segmentation import slic                  # noqa: PLC0415
        from skimage.measure import regionprops                # noqa: PLC0415
        from skimage.future.graph import rag_mean_color        # noqa: PLC0415
    except ImportError:
        return [], []

    # Superpixel segmentation
    segments = slic(gray, n_segments=n_segments, compactness=10, channel_axis=None)
    n_regions = segments.max() + 1
    nodes = [f"region_{i}" for i in range(n_regions)]

    # Mean colour per region
    mean_color = np.zeros(n_regions)
    for i in range(n_regions):
        mask = segments == i
        if mask.any():
            mean_color[i] = gray[mask].mean()

    edges = []

    # Adjacency: find touching pairs by scanning horizontal/vertical borders
    seen = set()
    for r in range(segments.shape[0]):
        for c in range(segments.shape[1] - 1):
            a, b = int(segments[r, c]), int(segments[r, c + 1])
            if a != b:
                key = (min(a, b), max(a, b))
                if key not in seen:
                    seen.add(key)
                    edges.append((a, b, "region_adj", 1.0))
    for r in range(segments.shape[0] - 1):
        for c in range(segments.shape[1]):
            a, b = int(segments[r, c]), int(segments[r + 1, c])
            if a != b:
                key = (min(a, b), max(a, b))
                if key not in seen:
                    seen.add(key)
                    edges.append((a, b, "region_adj", 1.0))

    # Colour similarity edges
    for i in range(n_regions):
        for j in range(i + 1, n_regions):
            diff = abs(float(mean_color[i]) - float(mean_color[j]))
            if diff < 0.1:
                edges.append((i, j, "color_similar", 1.0 - diff / 0.1))

    return nodes, edges


# ============================================================
#  SCALE 3: PATCH SIMILARITY GRAPH
# ============================================================

def _build_patch_graph(
    gray: np.ndarray,
    patch_size: int = 8,
    sim_threshold: float = 0.85,
    max_pairs: int = 3000,
) -> Tuple[List[str], List[Tuple[int, int, str, float]]]:
    """
    Build a patch-level graph.

    Each patch_size×patch_size block is a node.
    Edge types:
      "spatial_near" — patches that are immediate neighbours on the grid
      "patch_sim"    — patches with cosine similarity > sim_threshold
    """
    h, w = gray.shape
    ph = h // patch_size
    pw = w // patch_size
    if ph < 2 or pw < 2:
        return [], []

    n_patches = ph * pw

    # Extract patch feature vectors (flattened, normalised)
    patches = np.zeros((n_patches, patch_size * patch_size), dtype=np.float64)
    for i in range(ph):
        for j in range(pw):
            patch = gray[i * patch_size:(i + 1) * patch_size,
                         j * patch_size:(j + 1) * patch_size].ravel()
            norm = np.linalg.norm(patch)
            patches[i * pw + j] = patch / norm if norm > 1e-8 else patch

    nodes = [f"patch_{i // pw}_{i % pw}" for i in range(n_patches)]
    edges = []

    def pidx(i, j):
        return i * pw + j

    # Spatial neighbour edges
    for i in range(ph):
        for j in range(pw):
            if j + 1 < pw:
                edges.append((pidx(i, j), pidx(i, j + 1), "spatial_near", 1.0))
            if i + 1 < ph:
                edges.append((pidx(i, j), pidx(i + 1, j), "spatial_near", 1.0))

    # Cosine similarity edges (among non-neighbouring patches)
    sim = patches @ patches.T
    pairs = []
    for i in range(n_patches):
        for j in range(i + 1, n_patches):
            s = float(sim[i, j])
            if s > sim_threshold:
                pairs.append((i, j, s))

    pairs.sort(key=lambda x: -x[2])
    for i, j, s in pairs[:max_pairs]:
        edges.append((i, j, "patch_sim", s))

    return nodes, edges


# ============================================================
#  PUBLIC API
# ============================================================

def image_to_graph(
    img,
    use_regions: bool = True,
    use_patches: bool = True,
    max_dim: int = 256,
    patch_size: int = 8,
    n_segments: int = 64,
) -> RelationalGraph:
    """
    Convert an image to a multi-scale RelationalGraph.

    Parameters
    ----------
    img          : numpy array (H,W) or (H,W,C), PIL Image, or file path
    use_regions  : extract region adjacency graph (requires scikit-image)
    use_patches  : extract patch similarity graph
    max_dim      : maximum image dimension before downsampling
    patch_size   : patch side length for scale 3
    n_segments   : number of superpixel segments
    """
    gray = _load_image(img)
    gray = _resize_if_needed(gray, max_dim=max_dim)

    all_nodes: List[str] = []
    all_edges: List[Tuple[int, int, str, float]] = []
    node_offset = 0

    # Scale 1: geometric
    geo_nodes, geo_edges, n_geo = _build_geometric_graph(gray)
    all_nodes.extend(geo_nodes)
    all_edges.extend(geo_edges)
    node_offset += n_geo

    # Scale 2: region adjacency
    if use_regions:
        reg_nodes, reg_edges = _build_region_graph(gray, n_segments=n_segments)
        if reg_nodes:
            # Offset node indices to avoid clashes with scale 1
            reg_edges_off = [
                (s + node_offset, t + node_offset, et, w)
                for s, t, et, w in reg_edges
            ]
            all_nodes.extend(reg_nodes)
            all_edges.extend(reg_edges_off)
            node_offset += len(reg_nodes)

    # Scale 3: patch similarity
    if use_patches:
        pat_nodes, pat_edges = _build_patch_graph(gray, patch_size=patch_size)
        if pat_nodes:
            pat_edges_off = [
                (s + node_offset, t + node_offset, et, w)
                for s, t, et, w in pat_edges
            ]
            all_nodes.extend(pat_nodes)
            all_edges.extend(pat_edges_off)

    return RelationalGraph(nodes=all_nodes, edges=all_edges)


def image_to_operators(
    img,
    d: int = 3,
    angle_scale: float = 0.25,
    use_regions: bool = True,
    use_patches: bool = True,
    max_dim: int = 256,
    patch_size: int = 8,
    n_segments: int = 64,
    min_edges_per_type: int = 5,
    seed: int = 42,
) -> List[AlgebraicRelation]:
    """
    Full pipeline: image → operators Q_r ∈ O(n_gen).
    """
    graph = image_to_graph(
        img,
        use_regions=use_regions,
        use_patches=use_patches,
        max_dim=max_dim,
        patch_size=patch_size,
        n_segments=n_segments,
    )

    if graph.n_edges == 0:
        return []

    return graph_to_operators(
        graph,
        d=d,
        default_angle_scale=angle_scale,
        modality="image",
        min_edges_per_type=min_edges_per_type,
        seed=seed,
    )
