"""
Text → Relational Graph → Algebraic Operators
===============================================
Converts a text document (or list of texts) into a RelationalGraph and
then into SU(d) adjoint operators for the Killing classifier.

Relation types extracted:
  - syntactic_dep    : syntactic dependency edges from a dependency parse
  - cooccurrence     : word/token co-occurrence within a window
  - coreference      : coreference chains (entity mentions of the same referent)
  - semantic_sim     : high cosine-similarity token pairs (semantic proximity)
  - sentence_seq     : sentence-level sequential ordering
  - entity_type      : edges between tokens of the same named-entity type

Dependencies (all optional — graceful degradation if missing):
  - spacy              : for dependency parsing and NER
  - sentence_transformers : for semantic similarity edges
  (Without these, falls back to pure co-occurrence + sequential edges.)
"""

import re
import sys
import os
from collections import defaultdict
from typing import List, Optional, Tuple, Dict, Any

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from topology_to_algebra import RelationalGraph, AlgebraicRelation, graph_to_operators


# ============================================================
#  TOKEN / SENTENCE UTILITIES
# ============================================================

def _simple_tokenize(text: str) -> List[str]:
    """Minimal tokenizer: lowercase, split on non-alphanumeric."""
    tokens = re.findall(r"[a-z0-9']+", text.lower())
    stopwords = {
        "the", "a", "an", "is", "was", "are", "were", "be", "been",
        "have", "has", "had", "do", "does", "did", "will", "would",
        "can", "could", "should", "may", "might", "of", "in", "on",
        "at", "to", "for", "with", "by", "from", "and", "or", "but",
        "not", "it", "its", "this", "that", "these", "those",
    }
    return [t for t in tokens if t not in stopwords and len(t) > 1]


def _split_sentences(text: str) -> List[str]:
    """Split text into sentences on '.', '!', '?'."""
    parts = re.split(r'(?<=[.!?])\s+', text.strip())
    return [p for p in parts if len(p) > 0]


# ============================================================
#  GRAPH BUILDERS
# ============================================================

def _build_cooccurrence_graph(
    tokens: List[str],
    window: int = 4,
    vocab: Optional[Dict[str, int]] = None,
) -> Tuple[Dict[str, int], List[Tuple[int, int, str, float]]]:
    """
    Co-occurrence graph: edges between token pairs within a sliding window.
    Edge weight = normalised co-occurrence count.
    """
    if vocab is None:
        unique = list(dict.fromkeys(tokens))      # preserves order, dedup
        vocab = {t: i for i, t in enumerate(unique)}

    counts: Dict[Tuple[int, int], float] = defaultdict(float)
    for i, t in enumerate(tokens):
        if t not in vocab:
            continue
        for j in range(i + 1, min(i + window + 1, len(tokens))):
            u = tokens[j]
            if u not in vocab:
                continue
            a, b = vocab[t], vocab[u]
            if a != b:
                key = (min(a, b), max(a, b))
                counts[key] += 1.0

    if not counts:
        return vocab, []

    max_count = max(counts.values())
    edges = [
        (a, b, "cooccurrence", cnt / max_count)
        for (a, b), cnt in counts.items()
    ]
    return vocab, edges


def _build_sequential_graph(
    sentences: List[str],
    vocab: Dict[str, int],
) -> List[Tuple[int, int, str, float]]:
    """
    Sentence-level sequential edges: consecutive sentences share an edge.
    Weight = 1 / sentence_index (earlier pairs are more fundamental).
    """
    edges = []
    sent_tokens = [_simple_tokenize(s) for s in sentences]
    for i in range(len(sent_tokens) - 1):
        a_tokens = [t for t in sent_tokens[i] if t in vocab]
        b_tokens = [t for t in sent_tokens[i + 1] if t in vocab]
        # Connect last token of sentence i to first token of sentence i+1
        if a_tokens and b_tokens:
            a_node = vocab[a_tokens[-1]]
            b_node = vocab[b_tokens[0]]
            w = 1.0 / (i + 1)
            edges.append((a_node, b_node, "sentence_seq", w))
    return edges


def _build_spacy_graph(
    text: str,
    vocab: Dict[str, int],
) -> List[Tuple[int, int, str, float]]:
    """
    Dependency parse edges (requires spaCy).
    Edge types: syntactic_dep, coreference, entity_type.
    Gracefully returns [] if spaCy is not available.
    """
    try:
        import spacy                                   # noqa: PLC0415
        try:
            nlp = spacy.load("en_core_web_sm")
        except OSError:
            return []
        doc = nlp(text[:100_000])                      # cap at 100k chars
    except ImportError:
        return []

    edges = []

    # Dependency edges
    for token in doc:
        if token.dep_ == "ROOT":
            continue
        parent = token.head
        a_text = token.text.lower()
        b_text = parent.text.lower()
        if a_text in vocab and b_text in vocab and a_text != b_text:
            a, b = vocab[a_text], vocab[b_text]
            edges.append((a, b, "syntactic_dep", 1.0))

    # Named-entity type edges: tokens of the same entity type
    ent_groups: Dict[str, List[str]] = defaultdict(list)
    for ent in doc.ents:
        for tok in ent:
            t = tok.text.lower()
            if t in vocab:
                ent_groups[ent.label_].append(t)

    for label, group in ent_groups.items():
        unique_group = list(dict.fromkeys(group))
        for i in range(len(unique_group) - 1):
            a_text = unique_group[i]
            b_text = unique_group[i + 1]
            if a_text in vocab and b_text in vocab:
                a, b = vocab[a_text], vocab[b_text]
                edges.append((a, b, "entity_type", 0.8))

    return edges


def _build_semantic_sim_graph(
    unique_tokens: List[str],
    vocab: Dict[str, int],
    threshold: float = 0.75,
    max_pairs: int = 2000,
) -> List[Tuple[int, int, str, float]]:
    """
    Semantic similarity edges via sentence-transformers (optional).
    Edges between token pairs with cosine similarity > threshold.
    Gracefully returns [] if not available.
    """
    try:
        from sentence_transformers import SentenceTransformer         # noqa: PLC0415
    except ImportError:
        return []

    # Only embed tokens present in vocabulary (cap at 512 tokens)
    tokens_to_embed = [t for t in unique_tokens if t in vocab][:512]
    if len(tokens_to_embed) < 2:
        return []

    model = SentenceTransformer("all-MiniLM-L6-v2")
    embeddings = model.encode(tokens_to_embed, normalize_embeddings=True)

    # Cosine similarity matrix
    sim = embeddings @ embeddings.T

    edges = []
    n = len(tokens_to_embed)
    pairs = [(i, j, float(sim[i, j])) for i in range(n) for j in range(i + 1, n)]
    pairs.sort(key=lambda x: -x[2])

    for i, j, s in pairs[:max_pairs]:
        if s < threshold:
            break
        a = vocab[tokens_to_embed[i]]
        b = vocab[tokens_to_embed[j]]
        edges.append((a, b, "semantic_sim", s))

    return edges


# ============================================================
#  PUBLIC API
# ============================================================

def text_to_graph(
    text: str,
    window: int = 4,
    use_spacy: bool = True,
    use_semantic: bool = False,
    semantic_threshold: float = 0.75,
) -> RelationalGraph:
    """
    Convert a text document to a RelationalGraph.

    Parameters
    ----------
    text               : input document (any length)
    window             : co-occurrence window size
    use_spacy          : try to load spaCy for dep-parse and NER edges
    use_semantic       : try to use sentence-transformers for sim edges
    semantic_threshold : minimum cosine similarity for semantic edges
    """
    tokens = _simple_tokenize(text)
    sentences = _split_sentences(text)

    if not tokens:
        return RelationalGraph(nodes=[], edges=[])

    # Vocabulary
    unique_tokens = list(dict.fromkeys(tokens))
    vocab = {t: i for i, t in enumerate(unique_tokens)}

    all_edges: List[Tuple[int, int, str, float]] = []

    # 1. Co-occurrence
    _, cooc_edges = _build_cooccurrence_graph(tokens, window=window, vocab=vocab)
    all_edges.extend(cooc_edges)

    # 2. Sequential sentence
    seq_edges = _build_sequential_graph(sentences, vocab)
    all_edges.extend(seq_edges)

    # 3. Syntactic + NER (optional)
    if use_spacy:
        spacy_edges = _build_spacy_graph(text, vocab)
        all_edges.extend(spacy_edges)

    # 4. Semantic similarity (optional, slow)
    if use_semantic:
        sem_edges = _build_semantic_sim_graph(
            unique_tokens, vocab, threshold=semantic_threshold
        )
        all_edges.extend(sem_edges)

    return RelationalGraph(nodes=unique_tokens, edges=all_edges)


def text_to_operators(
    text: str,
    d: int = 3,
    angle_scale: float = 0.25,
    window: int = 4,
    use_spacy: bool = True,
    use_semantic: bool = False,
    semantic_threshold: float = 0.75,
    min_edges_per_type: int = 2,
    seed: int = 42,
) -> List[AlgebraicRelation]:
    """
    Full pipeline: text → operators Q_r ∈ O(n_gen).

    Returns a list of AlgebraicRelation, one per detected relation type.
    Typically produces 2–6 operators from plain text, more with spaCy.
    """
    graph = text_to_graph(
        text,
        window=window,
        use_spacy=use_spacy,
        use_semantic=use_semantic,
        semantic_threshold=semantic_threshold,
    )

    if graph.n_edges == 0:
        return []

    return graph_to_operators(
        graph,
        d=d,
        default_angle_scale=angle_scale,
        modality="text",
        min_edges_per_type=min_edges_per_type,
        seed=seed,
    )


def texts_to_operators(
    texts: List[str],
    d: int = 3,
    angle_scale: float = 0.25,
    window: int = 4,
    use_spacy: bool = True,
    use_semantic: bool = False,
    merge_mode: str = "concat",
    seed: int = 42,
) -> List[AlgebraicRelation]:
    """
    Convert multiple text documents to operators.

    merge_mode:
      "concat"  : concatenate all texts, build one big graph (default)
      "per_doc" : build one graph per document, concatenate operator lists
    """
    if merge_mode == "concat":
        merged = "\n\n".join(texts)
        return text_to_operators(merged, d=d, angle_scale=angle_scale,
                                 window=window, use_spacy=use_spacy,
                                 use_semantic=use_semantic, seed=seed)
    else:
        all_ops: List[AlgebraicRelation] = []
        for i, t in enumerate(texts):
            ops = text_to_operators(t, d=d, angle_scale=angle_scale,
                                    window=window, use_spacy=use_spacy,
                                    use_semantic=use_semantic, seed=seed + i)
            all_ops.extend(ops)
        return all_ops
