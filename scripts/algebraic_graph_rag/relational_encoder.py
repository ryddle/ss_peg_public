"""
Unified Multimodal Encoder
===========================
Single entry point for converting any supported data modality into
SU(d) adjoint operators Q_r ∈ O(n_gen) for the Killing classifier pipeline.

Supported modalities:
  text    → text_encoder.text_to_operators
  image   → image_encoder.image_to_operators
  audio   → audio_encoder.audio_to_operators
  svg     → svg_encoder.svg_to_operators

Usage example
-------------
>>> from relational_encoder import MultimodalEncoder
>>> enc = MultimodalEncoder(d=3, angle_scale=0.25)
>>> enc.encode_text("Quarks carry colour charge mediated by gluons.")
>>> enc.encode_svg(open("diagram.svg").read())
>>> ops_array, relations = enc.to_kg_operators()
>>> # ops_array.shape == (K, 8, 8)  for d=3
"""

import sys
import os
from typing import List, Optional, Dict, Any

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from topology_to_algebra import AlgebraicRelation, operators_to_array


# ============================================================
#  LAZY IMPORTS (import only if the modality is actually used)
# ============================================================

def _text_to_operators(text: str, **kw) -> List[AlgebraicRelation]:
    from text_encoder import text_to_operators       # noqa: PLC0415
    return text_to_operators(text, **kw)


def _image_to_operators(img, **kw) -> List[AlgebraicRelation]:
    from image_encoder import image_to_operators     # noqa: PLC0415
    return image_to_operators(img, **kw)


def _audio_to_operators(signal, sample_rate: int, **kw) -> List[AlgebraicRelation]:
    from audio_encoder import audio_to_operators     # noqa: PLC0415
    return audio_to_operators(signal, sample_rate, **kw)


def _svg_to_operators(svg: str, **kw) -> List[AlgebraicRelation]:
    from svg_encoder import svg_to_operators         # noqa: PLC0415
    return svg_to_operators(svg, **kw)


# ============================================================
#  MULTIMODAL ENCODER CLASS
# ============================================================

class MultimodalEncoder:
    """
    Accumulates operators from multiple data modalities.

    Each call to encode_* appends relations to an internal list.
    Call to_kg_operators() to obtain the stacked array.
    """

    def __init__(
        self,
        d: int = 3,
        angle_scale: float = 0.25,
        seed: int = 42,
    ):
        self.d = d
        self.angle_scale = angle_scale
        self.seed = seed
        self._relations: List[AlgebraicRelation] = []

    # ----------------------------------------------------------
    #  Per-modality encoders
    # ----------------------------------------------------------

    def encode_text(
        self,
        text: str,
        window: int = 4,
        use_spacy: bool = True,
        use_semantic: bool = False,
        min_edges_per_type: int = 2,
    ) -> "MultimodalEncoder":
        """Encode text and add relations to the pool."""
        ops = _text_to_operators(
            text,
            d=self.d,
            angle_scale=self.angle_scale,
            window=window,
            use_spacy=use_spacy,
            use_semantic=use_semantic,
            min_edges_per_type=min_edges_per_type,
            seed=self.seed,
        )
        self._relations.extend(ops)
        return self

    def encode_image(
        self,
        img,
        use_regions: bool = True,
        use_patches: bool = True,
        max_dim: int = 256,
        patch_size: int = 8,
        n_segments: int = 64,
        min_edges_per_type: int = 5,
    ) -> "MultimodalEncoder":
        """Encode image and add relations to the pool."""
        ops = _image_to_operators(
            img,
            d=self.d,
            angle_scale=self.angle_scale,
            use_regions=use_regions,
            use_patches=use_patches,
            max_dim=max_dim,
            patch_size=patch_size,
            n_segments=n_segments,
            min_edges_per_type=min_edges_per_type,
            seed=self.seed,
        )
        self._relations.extend(ops)
        return self

    def encode_audio(
        self,
        signal: np.ndarray,
        sample_rate: int,
        n_fft: int = 512,
        hop_length: int = 256,
        n_freq_bands: int = 16,
        max_frames: int = 128,
        use_chroma: bool = True,
        min_edges_per_type: int = 3,
    ) -> "MultimodalEncoder":
        """Encode audio signal and add relations to the pool."""
        ops = _audio_to_operators(
            signal,
            sample_rate,
            d=self.d,
            angle_scale=self.angle_scale,
            n_fft=n_fft,
            hop_length=hop_length,
            n_freq_bands=n_freq_bands,
            max_frames=max_frames,
            use_chroma=use_chroma,
            min_edges_per_type=min_edges_per_type,
            seed=self.seed,
        )
        self._relations.extend(ops)
        return self

    def encode_svg(
        self,
        svg: str,
        use_geometric: bool = True,
        min_edges_per_type: int = 2,
    ) -> "MultimodalEncoder":
        """Encode SVG document and add relations to the pool."""
        ops = _svg_to_operators(
            svg,
            d=self.d,
            angle_scale=self.angle_scale,
            use_geometric=use_geometric,
            min_edges_per_type=min_edges_per_type,
            seed=self.seed,
        )
        self._relations.extend(ops)
        return self

    def encode_graph(
        self,
        nodes: List[str],
        edges: List,
        modality: str = "graph",
        min_edges_per_type: int = 1,
    ) -> "MultimodalEncoder":
        """
        Encode an already-constructed RelationalGraph (or raw edge list).

        edges: list of (src_idx, dst_idx, edge_type, weight) tuples.
        """
        from topology_to_algebra import RelationalGraph, graph_to_operators  # noqa: PLC0415
        graph = RelationalGraph(nodes=nodes, edges=edges)
        ops = graph_to_operators(
            graph,
            d=self.d,
            default_angle_scale=self.angle_scale,
            modality=modality,
            min_edges_per_type=min_edges_per_type,
            seed=self.seed,
        )
        self._relations.extend(ops)
        return self

    # ----------------------------------------------------------
    #  Convenience batch encoder
    # ----------------------------------------------------------

    def encode_all(
        self,
        text: Optional[str] = None,
        image=None,
        audio: Optional[np.ndarray] = None,
        sample_rate: int = 22050,
        svg: Optional[str] = None,
        **kwargs,
    ) -> "MultimodalEncoder":
        """Encode all provided modalities in one call."""
        if text is not None:
            self.encode_text(text, **{k: v for k, v in kwargs.items()
                                       if k in ("window", "use_spacy",
                                                "use_semantic", "min_edges_per_type")})
        if image is not None:
            self.encode_image(image, **{k: v for k, v in kwargs.items()
                                         if k in ("use_regions", "use_patches",
                                                  "max_dim", "patch_size",
                                                  "n_segments")})
        if audio is not None:
            self.encode_audio(audio, sample_rate, **{k: v for k, v in kwargs.items()
                                                      if k in ("n_fft", "hop_length",
                                                               "n_freq_bands",
                                                               "max_frames",
                                                               "use_chroma")})
        if svg is not None:
            self.encode_svg(svg, **{k: v for k, v in kwargs.items()
                                     if k in ("use_geometric",)})
        return self

    # ----------------------------------------------------------
    #  Output
    # ----------------------------------------------------------

    @property
    def relations(self) -> List[AlgebraicRelation]:
        """All collected AlgebraicRelation objects."""
        return list(self._relations)

    @property
    def n_relations(self) -> int:
        return len(self._relations)

    def to_kg_operators(self) -> tuple:
        """
        Return (ops_array, relations).

        ops_array : (K, n_gen, n_gen) float64 numpy array
        relations : List[AlgebraicRelation] with metadata
        """
        if not self._relations:
            raise ValueError("No operators encoded yet.")
        arr = operators_to_array(self._relations)
        return arr, list(self._relations)

    def clear(self) -> "MultimodalEncoder":
        """Reset the operator pool."""
        self._relations.clear()
        return self

    def summary(self) -> str:
        """Human-readable summary of collected relations."""
        if not self._relations:
            return "MultimodalEncoder: no relations encoded."
        lines = [f"MultimodalEncoder — {self.n_relations} operators (d={self.d}):"]
        for rel in self._relations:
            lines.append(
                f"  [{rel.source_modality:10s}] {rel.name:30s} "
                f"  gap={rel.spectral_gap:.4f}  edges={rel.n_edges_in_subgraph}"
            )
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (f"MultimodalEncoder(d={self.d}, angle_scale={self.angle_scale}, "
                f"n_relations={self.n_relations})")


# ============================================================
#  FUNCTIONAL API (one-shot, stateless)
# ============================================================

def encode_multimodal(
    text: Optional[str] = None,
    image=None,
    audio: Optional[np.ndarray] = None,
    sample_rate: int = 22050,
    svg: Optional[str] = None,
    d: int = 3,
    angle_scale: float = 0.25,
    seed: int = 42,
    **kwargs,
) -> tuple:
    """
    One-shot functional interface.

    Returns (ops_array, relations) directly.
    """
    enc = MultimodalEncoder(d=d, angle_scale=angle_scale, seed=seed)
    enc.encode_all(text=text, image=image, audio=audio,
                   sample_rate=sample_rate, svg=svg, **kwargs)
    return enc.to_kg_operators()
