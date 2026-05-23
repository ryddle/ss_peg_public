"""
SVG → Relational Graph → Algebraic Operators
=============================================
Converts an SVG document (string or file) into a RelationalGraph and
then into SU(d) adjoint operators for the Killing classifier.

SVG has a natural graph structure (DOM tree) plus geometric relations.
Four relation types are extracted:

  Structural (DOM-level):
    "dom_child"    : parent → child element
    "dom_sibling"  : consecutive sibling elements

  Reference-based:
    "attr_ref"     : elements linked via xlink:href, fill="url(#id)", etc.
    "style_class"  : elements sharing the same CSS class

  Geometric (bounding-box level, best-effort):
    "bbox_overlap" : element bounding boxes that overlap in 2-D
    "bbox_near"    : element centres within a proximity threshold

Edge weights encode visual/structural salience:
  - transform magnitude (Frobenius norm of matrix part) for geometric edges
  - 1 / depth for DOM edges (higher = shallower in tree)
  - 1.0 for reference/class edges

Only stdlib + numpy — no external SVG renderer required.
"""

import re
import sys
import os
import xml.etree.ElementTree as ET
from typing import List, Tuple, Optional, Dict, Any
from dataclasses import dataclass, field

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from topology_to_algebra import RelationalGraph, AlgebraicRelation, graph_to_operators


# ============================================================
#  SVG TRANSFORM PARSER
# ============================================================

def _parse_transform(transform_str: str) -> np.ndarray:
    """
    Parse an SVG transform attribute and return the 3×3 homogeneous matrix.

    Supported: matrix(), translate(), rotate(), scale().
    Returns identity if parse fails.
    """
    M = np.eye(3)
    if not transform_str:
        return M

    # Split multiple transforms and apply left-to-right
    token_re = re.compile(
        r'(matrix|translate|rotate|scale)\s*\(([^)]*)\)',
        re.IGNORECASE,
    )
    for m in token_re.finditer(transform_str):
        func = m.group(1).lower()
        args = [float(x) for x in re.split(r'[\s,]+', m.group(2).strip()) if x]
        try:
            if func == "matrix" and len(args) == 6:
                a, b, c, d, e, f = args
                T = np.array([[a, c, e], [b, d, f], [0, 0, 1]])
            elif func == "translate":
                tx = args[0]
                ty = args[1] if len(args) > 1 else 0.0
                T = np.array([[1, 0, tx], [0, 1, ty], [0, 0, 1]])
            elif func == "rotate":
                deg = args[0]
                rad = np.deg2rad(deg)
                cx = args[1] if len(args) > 1 else 0.0
                cy = args[2] if len(args) > 2 else 0.0
                cos_a, sin_a = np.cos(rad), np.sin(rad)
                T = np.array([
                    [cos_a, -sin_a, cx * (1 - cos_a) + cy * sin_a],
                    [sin_a,  cos_a, cy * (1 - cos_a) - cx * sin_a],
                    [0, 0, 1],
                ])
            elif func == "scale":
                sx = args[0]
                sy = args[1] if len(args) > 1 else sx
                T = np.array([[sx, 0, 0], [0, sy, 0], [0, 0, 1]])
            else:
                continue
            M = M @ T
        except Exception:
            continue

    return M


def _transform_magnitude(M: np.ndarray) -> float:
    """Frobenius norm of the 2×2 linear part of a 3×3 homogeneous matrix."""
    return float(np.linalg.norm(M[:2, :2]))


# ============================================================
#  BOUNDING BOX EXTRACTION (best-effort)
# ============================================================

def _bbox_from_element(elem: ET.Element) -> Optional[Tuple[float, float, float, float]]:
    """
    Attempt to extract (x, y, width, height) from common SVG attributes.
    Returns None if not possible.
    """
    tag = elem.tag.split("}")[-1].lower() if "}" in elem.tag else elem.tag.lower()
    attrib = elem.attrib

    def fget(key, default=None):
        val = attrib.get(key, "")
        try:
            return float(re.sub(r"[a-z%]+$", "", val))
        except ValueError:
            return default

    if tag in ("rect", "image", "foreignobject", "use"):
        x, y = fget("x", 0.0), fget("y", 0.0)
        w, h = fget("width"), fget("height")
        if w is not None and h is not None:
            return x, y, w, h

    if tag == "circle":
        cx, cy = fget("cx", 0.0), fget("cy", 0.0)
        r = fget("r")
        if r is not None:
            return cx - r, cy - r, 2 * r, 2 * r

    if tag == "ellipse":
        cx, cy = fget("cx", 0.0), fget("cy", 0.0)
        rx, ry = fget("rx"), fget("ry")
        if rx is not None and ry is not None:
            return cx - rx, cy - ry, 2 * rx, 2 * ry

    if tag == "line":
        x1, y1 = fget("x1", 0.0), fget("y1", 0.0)
        x2, y2 = fget("x2", 0.0), fget("y2", 0.0)
        return (min(x1, x2), min(y1, y2),
                abs(x2 - x1) or 1.0, abs(y2 - y1) or 1.0)

    return None


def _boxes_overlap(b1, b2) -> bool:
    if b1 is None or b2 is None:
        return False
    x1, y1, w1, h1 = b1
    x2, y2, w2, h2 = b2
    return not (x1 + w1 < x2 or x2 + w2 < x1 or y1 + h1 < y2 or y2 + h2 < y1)


def _box_centre(bbox) -> Optional[Tuple[float, float]]:
    if bbox is None:
        return None
    x, y, w, h = bbox
    return x + w / 2, y + h / 2


# ============================================================
#  DOM WALKER
# ============================================================

@dataclass
class _SVGElement:
    idx: int
    tag: str
    attrib: Dict[str, str]
    depth: int
    parent_idx: Optional[int]
    transform: np.ndarray
    bbox: Optional[Tuple[float, float, float, float]]


def _walk_svg_tree(root: ET.Element) -> List[_SVGElement]:
    """Flatten SVG DOM tree into a list of _SVGElement records."""
    elements: List[_SVGElement] = []
    idx_counter = [0]

    def _strip_ns(tag: str) -> str:
        return tag.split("}")[-1] if "}" in tag else tag

    def _recurse(elem: ET.Element, depth: int, parent_idx: Optional[int]):
        # Skip non-visual / meta elements
        tag = _strip_ns(elem.tag).lower()
        if tag in ("defs", "style", "script", "title", "desc", "metadata"):
            return

        my_idx = idx_counter[0]
        idx_counter[0] += 1

        transform_str = elem.attrib.get("transform", "")
        M = _parse_transform(transform_str)
        bbox = _bbox_from_element(elem)

        elements.append(_SVGElement(
            idx=my_idx,
            tag=tag,
            attrib=dict(elem.attrib),
            depth=depth,
            parent_idx=parent_idx,
            transform=M,
            bbox=bbox,
        ))

        children = list(elem)
        for i, child in enumerate(children):
            _recurse(child, depth + 1, my_idx)

    for child in root:
        _recurse(child, 0, None)

    return elements


# ============================================================
#  EDGE BUILDERS
# ============================================================

def _dom_edges(elements: List[_SVGElement]) -> List[Tuple[int, int, str, float]]:
    edges = []
    # Child edges
    for el in elements:
        if el.parent_idx is not None:
            depth_w = 1.0 / (el.depth + 1)
            edges.append((el.parent_idx, el.idx, "dom_child", depth_w))

    # Sibling edges (consecutive children of the same parent)
    by_parent: Dict[Optional[int], List[_SVGElement]] = {}
    for el in elements:
        by_parent.setdefault(el.parent_idx, []).append(el)

    for parent_idx, children in by_parent.items():
        if len(children) < 2:
            continue
        for i in range(len(children) - 1):
            a, b = children[i], children[i + 1]
            depth_w = 1.0 / (a.depth + 1)
            edges.append((a.idx, b.idx, "dom_sibling", depth_w))

    return edges


def _reference_edges(elements: List[_SVGElement]) -> List[Tuple[int, int, str, float]]:
    """
    href / xlink:href and fill='url(#id)' reference edges.
    """
    id_to_idx: Dict[str, int] = {}
    for el in elements:
        eid = el.attrib.get("id", "")
        if eid:
            id_to_idx[eid] = el.idx

    edges = []
    href_re = re.compile(r"url\(#([^)]+)\)|^#(.+)$")

    for el in elements:
        # xlink:href / href
        href = (el.attrib.get("{http://www.w3.org/1999/xlink}href", "") or
                el.attrib.get("href", ""))
        if href.startswith("#"):
            target_id = href[1:]
            if target_id in id_to_idx:
                edges.append((el.idx, id_to_idx[target_id], "attr_ref", 1.0))

        # fill / stroke URL references
        for attr_key in ("fill", "stroke", "clip-path", "mask"):
            val = el.attrib.get(attr_key, "")
            m = href_re.match(val)
            if m:
                target_id = m.group(1) or m.group(2)
                if target_id in id_to_idx:
                    edges.append((el.idx, id_to_idx[target_id], "attr_ref", 0.8))

    return edges


def _class_edges(elements: List[_SVGElement]) -> List[Tuple[int, int, str, float]]:
    """
    Style-class edges: elements sharing the same CSS class.
    Each class group forms a clique (bounded by max_per_class).
    """
    class_groups: Dict[str, List[int]] = {}
    for el in elements:
        cls_str = el.attrib.get("class", "")
        for cls in cls_str.split():
            cls = cls.strip()
            if cls:
                class_groups.setdefault(cls, []).append(el.idx)

    edges = []
    max_per_class = 20
    for cls, group in class_groups.items():
        group = group[:max_per_class]
        for i in range(len(group) - 1):
            edges.append((group[i], group[i + 1], "style_class", 1.0))

    return edges


def _geometric_edges(
    elements: List[_SVGElement],
    proximity_threshold: Optional[float] = None,
) -> List[Tuple[int, int, str, float]]:
    """
    Geometric overlap and proximity edges between elements with bounding boxes.
    """
    # Auto threshold: 10% of the average bounding box diagonal
    diagonals = []
    for el in elements:
        if el.bbox:
            _, _, w, h = el.bbox
            diagonals.append(np.sqrt(w ** 2 + h ** 2))
    if diagonals:
        avg_diag = float(np.mean(diagonals))
        if proximity_threshold is None:
            proximity_threshold = avg_diag * 0.15
    else:
        proximity_threshold = proximity_threshold or 50.0

    with_bbox = [el for el in elements if el.bbox is not None]
    edges = []
    n = len(with_bbox)

    for i in range(n):
        for j in range(i + 1, n):
            a, b = with_bbox[i], with_bbox[j]

            # Overlap
            if _boxes_overlap(a.bbox, b.bbox):
                # Weight = transform magnitude of the more transformed element
                tm = max(_transform_magnitude(a.transform),
                         _transform_magnitude(b.transform))
                w = np.clip(tm / 4.0, 0.1, 1.0)
                edges.append((a.idx, b.idx, "bbox_overlap", float(w)))
                continue

            # Proximity
            ca, cb = _box_centre(a.bbox), _box_centre(b.bbox)
            if ca is not None and cb is not None:
                dist = np.sqrt((ca[0] - cb[0]) ** 2 + (ca[1] - cb[1]) ** 2)
                if dist < proximity_threshold:
                    w = float(1.0 - dist / proximity_threshold)
                    edges.append((a.idx, b.idx, "bbox_near", w))

    return edges


# ============================================================
#  PUBLIC API
# ============================================================

def svg_to_graph(
    svg: str,
    use_geometric: bool = True,
) -> RelationalGraph:
    """
    Convert an SVG document string to a RelationalGraph.

    Parameters
    ----------
    svg           : SVG content as string (or file path ending in .svg)
    use_geometric : include bounding-box overlap / proximity edges
    """
    if svg.endswith(".svg"):
        with open(svg, "r", encoding="utf-8") as f:
            svg = f.read()

    # Strip DOCTYPE and XML declaration for cleaner parsing
    svg = re.sub(r"<!DOCTYPE[^>]*>", "", svg)
    svg = re.sub(r"<\?xml[^?]*\?>", "", svg)

    try:
        root = ET.fromstring(svg.strip())
    except ET.ParseError:
        # Try wrapping in a dummy SVG element
        try:
            root = ET.fromstring(f"<svg>{svg}</svg>")
        except ET.ParseError:
            return RelationalGraph(nodes=[], edges=[])

    elements = _walk_svg_tree(root)
    if not elements:
        return RelationalGraph(nodes=[], edges=[])

    nodes = [f"{el.tag}_{el.idx}" for el in elements]
    all_edges: List[Tuple[int, int, str, float]] = []

    all_edges.extend(_dom_edges(elements))
    all_edges.extend(_reference_edges(elements))
    all_edges.extend(_class_edges(elements))
    if use_geometric:
        all_edges.extend(_geometric_edges(elements))

    return RelationalGraph(nodes=nodes, edges=all_edges)


def svg_to_operators(
    svg: str,
    d: int = 3,
    angle_scale: float = 0.25,
    use_geometric: bool = True,
    min_edges_per_type: int = 2,
    seed: int = 42,
) -> List[AlgebraicRelation]:
    """
    Full pipeline: SVG document → operators Q_r ∈ O(n_gen).
    """
    graph = svg_to_graph(svg, use_geometric=use_geometric)

    if graph.n_edges == 0:
        return []

    return graph_to_operators(
        graph,
        d=d,
        default_angle_scale=angle_scale,
        modality="svg",
        min_edges_per_type=min_edges_per_type,
        seed=seed,
    )
