"""
Audio → Relational Graph → Algebraic Operators
================================================
Converts an audio signal (1-D numpy array + sample rate) into a
RelationalGraph and then into SU(d) adjoint operators for the Killing
classifier.

Three graph layers are extracted:

  Layer 1 — Frequency (STFT magnitude, n_freq_bands nodes)
    "harmonic"    : frequency bands whose power ratio approximates 2:1
    "spectral_nb" : neighbouring frequency bands

  Layer 2 — Temporal (frame-level, n_frames nodes)
    "temporal_seq" : consecutive frames (short-term dynamics)
    "onset"        : frames around energy onsets connected to their neighbours

  Layer 3 — Chroma / tonal (12 chroma bins, pure numpy FFT)
    "chroma_nb"    : neighbouring chroma bins
    "chroma_harm"  : chroma bins separated by a perfect 5th (7 semitones)

Dependencies:
  - numpy (required)
  - librosa (optional) — if available, uses librosa's STFT and onset detection;
    otherwise falls back to pure numpy FFT.
"""

import sys
import os
from typing import List, Tuple, Optional

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from topology_to_algebra import RelationalGraph, AlgebraicRelation, graph_to_operators


# ============================================================
#  STFT / SPECTROGRAM (pure numpy fallback)
# ============================================================

def _stft_numpy(
    signal: np.ndarray,
    n_fft: int = 512,
    hop_length: int = 256,
) -> np.ndarray:
    """
    Short-Time Fourier Transform using numpy.

    Returns magnitude spectrogram of shape (n_fft//2+1, n_frames).
    """
    n_frames = 1 + (len(signal) - n_fft) // hop_length
    if n_frames <= 0:
        # Signal shorter than one frame: pad and take a single frame
        signal = np.pad(signal, (0, n_fft - len(signal)))
        n_frames = 1

    window = np.hanning(n_fft).astype(np.float64)
    n_bins = n_fft // 2 + 1
    S = np.zeros((n_bins, max(1, n_frames)), dtype=np.float64)

    for t in range(n_frames):
        start = t * hop_length
        frame = signal[start:start + n_fft]
        if len(frame) < n_fft:
            frame = np.pad(frame, (0, n_fft - len(frame)))
        fft = np.fft.rfft(frame * window)
        S[:, t] = np.abs(fft)

    return S


def _get_spectrogram(
    signal: np.ndarray,
    sample_rate: int,
    n_fft: int = 512,
    hop_length: int = 256,
) -> np.ndarray:
    """
    Return magnitude spectrogram. Uses librosa if available, else numpy.
    Shape: (n_freq_bins, n_frames)
    """
    try:
        import librosa                              # noqa: PLC0415
        S = np.abs(librosa.stft(signal.astype(np.float32), n_fft=n_fft,
                                hop_length=hop_length))
        return S.astype(np.float64)
    except ImportError:
        return _stft_numpy(signal.astype(np.float64), n_fft=n_fft,
                           hop_length=hop_length)


# ============================================================
#  LAYER 1: FREQUENCY GRAPH
# ============================================================

def _build_frequency_graph(
    S: np.ndarray,
    n_freq_bands: int = 16,
) -> Tuple[List[str], List[Tuple[int, int, str, float]]]:
    """
    Aggregate spectrogram into n_freq_bands mel-like bands.
    Build harmonic and neighbourhood edges between bands.
    """
    n_bins, n_frames = S.shape
    bins_per_band = max(1, n_bins // n_freq_bands)

    # Band mean power
    band_power = np.zeros(n_freq_bands)
    for k in range(n_freq_bands):
        start = k * bins_per_band
        end = min(start + bins_per_band, n_bins)
        band_power[k] = S[start:end, :].mean()

    max_power = band_power.max() if band_power.max() > 0 else 1.0
    band_power_norm = band_power / max_power

    nodes = [f"freq_{k}" for k in range(n_freq_bands)]
    edges = []

    # Spectral neighbour edges
    for k in range(n_freq_bands - 1):
        w = 0.5 * (band_power_norm[k] + band_power_norm[k + 1])
        edges.append((k, k + 1, "spectral_nb", float(w)))

    # Harmonic edges: approximate integer frequency ratios
    # Bands k and j are "harmonic" if j ≈ 2k (one octave up)
    for k in range(1, n_freq_bands):
        j = 2 * k
        if j < n_freq_bands:
            w = 0.5 * (band_power_norm[k] + band_power_norm[j])
            edges.append((k, j, "harmonic", float(w)))

    return nodes, edges


# ============================================================
#  LAYER 2: TEMPORAL GRAPH
# ============================================================

def _build_temporal_graph(
    S: np.ndarray,
    max_frames: int = 128,
) -> Tuple[List[str], List[Tuple[int, int, str, float]]]:
    """
    Build temporal graph from spectral flux between adjacent frames.

    Nodes: time frames (downsampled to max_frames if needed).
    Edges:
      "temporal_seq" : every consecutive pair
      "onset"        : frames where spectral flux exceeds a threshold
    """
    n_bins, n_frames = S.shape

    # Downsample frames
    if n_frames > max_frames:
        step = n_frames // max_frames
        S = S[:, ::step][:, :max_frames]
        n_frames = S.shape[1]

    nodes = [f"frame_{t}" for t in range(n_frames)]
    edges = []

    # Spectral flux: L2 difference between consecutive frames
    flux = np.zeros(n_frames)
    for t in range(1, n_frames):
        diff = S[:, t] - S[:, t - 1]
        flux[t] = float(np.sqrt((diff ** 2).sum()))

    max_flux = flux.max() if flux.max() > 0 else 1.0
    flux_norm = flux / max_flux

    # Sequential edges
    for t in range(n_frames - 1):
        w = float(flux_norm[t + 1])
        edges.append((t, t + 1, "temporal_seq", max(0.05, w)))

    # Onset edges: frames with flux > 0.5 connected to neighbours
    onset_threshold = 0.5
    for t in range(1, n_frames - 1):
        if flux_norm[t] > onset_threshold:
            edges.append((t - 1, t + 1, "onset", float(flux_norm[t])))

    return nodes, edges


# ============================================================
#  LAYER 3: CHROMA / TONAL GRAPH
# ============================================================

def _chroma_from_spectrogram(
    S: np.ndarray,
    sample_rate: int,
    n_fft: int = 512,
) -> np.ndarray:
    """
    Compute chroma (12-bin pitch class) vector from magnitude spectrogram.

    Uses frequency-to-midi mapping and wraps into 12 pitch classes.
    Returns (12,) chroma vector normalised to [0,1].
    """
    n_bins = S.shape[0]
    chroma = np.zeros(12)

    for k in range(1, n_bins):
        freq_hz = k * sample_rate / n_fft
        # Map frequency to MIDI note
        if freq_hz <= 0:
            continue
        midi = 69.0 + 12.0 * np.log2(freq_hz / 440.0)
        pitch_class = int(round(midi)) % 12
        chroma[pitch_class] += S[k, :].mean()

    max_c = chroma.max()
    return chroma / max_c if max_c > 0 else chroma


def _build_chroma_graph(
    S: np.ndarray,
    sample_rate: int,
    n_fft: int = 512,
) -> Tuple[List[str], List[Tuple[int, int, str, float]]]:
    """
    12-note chroma graph.
    Edges:
      "chroma_nb"   : adjacent semitones
      "chroma_harm" : perfect-5th interval (7 semitones)
    """
    chroma = _chroma_from_spectrogram(S, sample_rate, n_fft)
    note_names = ["C", "C#", "D", "D#", "E", "F",
                  "F#", "G", "G#", "A", "A#", "B"]
    nodes = note_names[:]
    edges = []

    # Adjacent semitones
    for i in range(12):
        j = (i + 1) % 12
        w = 0.5 * (float(chroma[i]) + float(chroma[j]))
        edges.append((i, j, "chroma_nb", max(0.05, w)))

    # Perfect 5th (7 semitones = most consonant interval)
    for i in range(12):
        j = (i + 7) % 12
        w = 0.5 * (float(chroma[i]) + float(chroma[j]))
        edges.append((i, j, "chroma_harm", max(0.05, w)))

    return nodes, edges


# ============================================================
#  PUBLIC API
# ============================================================

def audio_to_graph(
    signal: np.ndarray,
    sample_rate: int,
    n_fft: int = 512,
    hop_length: int = 256,
    n_freq_bands: int = 16,
    max_frames: int = 128,
    use_chroma: bool = True,
) -> RelationalGraph:
    """
    Convert an audio signal to a multi-layer RelationalGraph.

    Parameters
    ----------
    signal      : 1-D float numpy array (mono waveform)
    sample_rate : samples per second (e.g. 22050)
    n_fft       : FFT window size
    hop_length  : hop size between frames
    n_freq_bands: number of frequency band nodes
    max_frames  : maximum temporal frame nodes
    use_chroma  : include chroma/tonal graph layer
    """
    signal = np.asarray(signal, dtype=np.float64).ravel()
    # Normalise amplitude
    max_amp = np.abs(signal).max()
    if max_amp > 0:
        signal = signal / max_amp

    S = _get_spectrogram(signal, sample_rate, n_fft=n_fft, hop_length=hop_length)

    all_nodes: List[str] = []
    all_edges: List[Tuple[int, int, str, float]] = []
    node_offset = 0

    # Layer 1: frequency
    freq_nodes, freq_edges = _build_frequency_graph(S, n_freq_bands=n_freq_bands)
    all_nodes.extend(freq_nodes)
    all_edges.extend(freq_edges)
    node_offset += len(freq_nodes)

    # Layer 2: temporal
    temp_nodes, temp_edges = _build_temporal_graph(S, max_frames=max_frames)
    temp_edges_off = [
        (s + node_offset, t + node_offset, et, w)
        for s, t, et, w in temp_edges
    ]
    all_nodes.extend(temp_nodes)
    all_edges.extend(temp_edges_off)
    node_offset += len(temp_nodes)

    # Layer 3: chroma
    if use_chroma:
        chr_nodes, chr_edges = _build_chroma_graph(S, sample_rate, n_fft=n_fft)
        chr_edges_off = [
            (s + node_offset, t + node_offset, et, w)
            for s, t, et, w in chr_edges
        ]
        all_nodes.extend(chr_nodes)
        all_edges.extend(chr_edges_off)

    return RelationalGraph(nodes=all_nodes, edges=all_edges)


def audio_to_operators(
    signal: np.ndarray,
    sample_rate: int,
    d: int = 3,
    angle_scale: float = 0.25,
    n_fft: int = 512,
    hop_length: int = 256,
    n_freq_bands: int = 16,
    max_frames: int = 128,
    use_chroma: bool = True,
    min_edges_per_type: int = 3,
    seed: int = 42,
) -> List[AlgebraicRelation]:
    """
    Full pipeline: audio signal → operators Q_r ∈ O(n_gen).
    """
    graph = audio_to_graph(
        signal,
        sample_rate,
        n_fft=n_fft,
        hop_length=hop_length,
        n_freq_bands=n_freq_bands,
        max_frames=max_frames,
        use_chroma=use_chroma,
    )

    if graph.n_edges == 0:
        return []

    return graph_to_operators(
        graph,
        d=d,
        default_angle_scale=angle_scale,
        modality="audio",
        min_edges_per_type=min_edges_per_type,
        seed=seed,
    )
