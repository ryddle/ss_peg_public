#!/usr/bin/env python3
"""
P8: Ramanujan Analysis of Non-Compact Lie Algebras
===================================================

P1-P7 tested only compact (and abelian) Lie algebras.  Physics requires
non-compact groups for spacetime symmetries:

  - so(3,1)    Lorentz group                   dim 6
  - so(4,1)    de Sitter                        dim 10
  - so(3,2)    anti-de Sitter / conformal 2+1   dim 10
  - so(4,2)    conformal group in 3+1           dim 15
  - sl(2,R)    simplest non-compact semisimple   dim 3
  - su(2,1)    pseudo-unitary                    dim 8
  - sl(2,C)_R  Lorentz algebra (real form)       dim 6

Key question:  The adjoint commutation graph depends ONLY on the structure
constants |[T_a, T_b]| != 0, which are the SAME for a real Lie algebra and
its compact real form (they share the complexification).  But the signs
in the commutation relations differ (Killing form signature).

Hypothesis: The adjoint graph of so(3,1) is ISOMORPHIC to that of so(4)
(= su(2) + su(2)), because the non-zero pattern of [T_a, T_b] is identical.

This script tests that hypothesis and extends the analysis to all
physically relevant non-compact algebras.

Author: SS-PEG project
Date:   June 2025
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(__file__))
from ramanujan_extended import (
    roots_A, roots_B, roots_C, roots_D,
    roots_G2, roots_F4, roots_E6,
    adjacency_from_roots, spectral_analysis
)

HASH_DIM_LIMIT = 200


def commutator(A, B):
    return A @ B - B @ A


def adjacency_from_generators(generators, threshold=1e-8):
    """Build adjacency from explicit matrix generators."""
    n = len(generators)
    adj = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            C = commutator(generators[i], generators[j])
            if np.linalg.norm(C) > threshold:
                adj[i, j] = 1
                adj[j, i] = 1
    return adj


def killing_form(generators):
    """Compute the Killing form K_{ab} = Tr(ad_a . ad_b)."""
    n = len(generators)
    # Build ad matrices
    ad = np.zeros((n, n, n), dtype=np.complex128)
    for a in range(n):
        for b in range(n):
            C = commutator(generators[a], generators[b])
            # Expand C in the generator basis
            for c in range(n):
                # Use trace inner product
                ad[a, b, c] = np.trace(C @ generators[c].conj().T).real
    # Killing: K_{ab} = sum_c ad[a,c,d] * ad[b,d,c] (trace of ad_a . ad_b)
    K = np.zeros((n, n), dtype=np.float64)
    for a in range(n):
        for b in range(n):
            K[a, b] = np.sum(ad[a] * ad[b].T).real
    return K


def killing_signature(K):
    """Return (n_pos, n_neg, n_zero) eigenvalue counts of the Killing form."""
    eigs = np.linalg.eigvalsh(K)
    n_pos = np.sum(eigs > 1e-10)
    n_neg = np.sum(eigs < -1e-10)
    n_zero = len(eigs) - n_pos - n_neg
    return int(n_pos), int(n_neg), int(n_zero)


# ─── Non-compact generator factories ─────────────────────────────────

def generators_so(n):
    """Compact so(n): antisymmetric real n x n matrices."""
    gens = []
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n), dtype=np.complex128)
            E[i, j] = 1
            E[j, i] = -1
            gens.append(E)
    return gens


def generators_so_pq(p, q):
    """
    Non-compact so(p,q): preserves eta = diag(+1,...,+1,-1,...,-1) [p pluses, q minuses].
    Generators X satisfy X^T eta + eta X = 0.

    Basis:
      J_{ij} (i<j, both in 'space' block 0..p-1): X[i,j]=-X[j,i]=1  (compact, rotation)
      J_{ab} (a<b, both in 'time' block p..p+q-1): same form           (compact, rotation)
      K_{ia} (i in space, a in time): X[i,a]=X[a,i]=1                  (non-compact, boost)
    """
    n = p + q
    gens = []
    # Compact: within first p indices
    for i in range(p):
        for j in range(i + 1, p):
            E = np.zeros((n, n), dtype=np.complex128)
            E[i, j] = 1; E[j, i] = -1
            gens.append(E)
    # Compact: within last q indices
    for a in range(p, n):
        for b in range(a + 1, n):
            E = np.zeros((n, n), dtype=np.complex128)
            E[a, b] = 1; E[b, a] = -1
            gens.append(E)
    # Non-compact boosts: mixed indices
    for i in range(p):
        for a in range(p, n):
            E = np.zeros((n, n), dtype=np.complex128)
            E[i, a] = 1; E[a, i] = 1   # symmetric (boost)
            gens.append(E)
    return gens


def generators_su_pq(p, q):
    """
    Non-compact su(p,q): preserves eta = diag(+1,...,-1,...).
    Generators X satisfy X^dag eta + eta X = 0, Tr(X) = 0.
    dim = (p+q)^2 - 1.
    """
    n = p + q
    eta = np.diag([1.0] * p + [-1.0] * q).astype(np.complex128)
    gens = []

    # Off-diagonal: for i,j with eta[i,i] = eta[j,j] (compact block)
    #   symmetric and antisymmetric hermitian
    # For i,j with eta[i,i] != eta[j,j] (non-compact block)
    #   different sign structure
    for i in range(n):
        for j in range(i + 1, n):
            if eta[i, i] == eta[j, j]:
                # Compact pair — same as su(n)
                E1 = np.zeros((n, n), dtype=np.complex128)
                E1[i, j] = 1; E1[j, i] = 1
                gens.append(E1 / 2)
                E2 = np.zeros((n, n), dtype=np.complex128)
                E2[i, j] = -1j; E2[j, i] = 1j
                gens.append(E2 / 2)
            else:
                # Non-compact pair — boost generators
                E1 = np.zeros((n, n), dtype=np.complex128)
                E1[i, j] = 1; E1[j, i] = -1
                gens.append(E1 / 2)
                E2 = np.zeros((n, n), dtype=np.complex128)
                E2[i, j] = 1j; E2[j, i] = 1j
                gens.append(E2 / 2)

    # Diagonal (Cartan): same as su(n)
    for k in range(1, n):
        D = np.zeros((n, n), dtype=np.complex128)
        for i in range(k):
            D[i, i] = 1
        D[k, k] = -k
        D /= np.sqrt(2 * k * (k + 1))
        gens.append(D)

    return gens


def generators_sl_n_R(n):
    """
    sl(n, R): real traceless n x n matrices.
    dim = n^2 - 1.
    Basis: E_{ij} (off-diagonal) + diagonal traceless.
    """
    gens = []
    # Off-diagonal
    for i in range(n):
        for j in range(n):
            if i != j:
                E = np.zeros((n, n), dtype=np.complex128)
                E[i, j] = 1
                gens.append(E)
    # Diagonal traceless
    for k in range(1, n):
        D = np.zeros((n, n), dtype=np.complex128)
        for i in range(k):
            D[i, i] = 1
        D[k, k] = -k
        D /= np.sqrt(k * (k + 1))
        gens.append(D)
    return gens


def generators_su(n):
    """Compact su(n): traceless anti-hermitian or hermitian convention."""
    gens = []
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n), dtype=np.complex128)
            E[i, j] = 1; E[j, i] = 1
            gens.append(E / 2)
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n), dtype=np.complex128)
            E[i, j] = -1j; E[j, i] = 1j
            gens.append(E / 2)
    for k in range(1, n):
        D = np.zeros((n, n), dtype=np.complex128)
        for i in range(k):
            D[i, i] = 1
        D[k, k] = -k
        D /= np.sqrt(2 * k * (k + 1))
        gens.append(D)
    return gens


# ─── Catalog ─────────────────────────────────────────────────────────

def build_catalog():
    """Build ordered catalog of compact/non-compact pairs."""
    cat = OrderedDict()

    # ── Compact references (from root systems) ──
    cat['so(3)'] = {
        'type': 'compact', 'method': 'generators',
        'gen_fn': lambda: generators_so(3),
        'dim': 3, 'complexification': 'A1', 'physics': 'Rotation group'}
    cat['so(4)'] = {
        'type': 'compact', 'method': 'generators',
        'gen_fn': lambda: generators_so(4),
        'dim': 6, 'complexification': 'D2=A1+A1', 'physics': '4D rotation'}
    cat['su(2)'] = {
        'type': 'compact', 'method': 'roots',
        'root_fn': lambda: roots_A(1),
        'dim': 3, 'complexification': 'A1', 'physics': 'Isospin'}
    cat['su(3)'] = {
        'type': 'compact', 'method': 'roots',
        'root_fn': lambda: roots_A(2),
        'dim': 8, 'complexification': 'A2', 'physics': 'Color SU(3)'}
    cat['so(5)'] = {
        'type': 'compact', 'method': 'roots',
        'root_fn': lambda: roots_B(2),
        'dim': 10, 'complexification': 'B2', 'physics': 'Compact B2'}
    cat['so(6)'] = {
        'type': 'compact', 'method': 'roots',
        'root_fn': lambda: roots_D(3),
        'dim': 15, 'complexification': 'D3=A3', 'physics': 'Compact D3'}
    cat['su(4)'] = {
        'type': 'compact', 'method': 'roots',
        'root_fn': lambda: roots_A(3),
        'dim': 15, 'complexification': 'A3', 'physics': 'Compact A3'}
    cat['so(10)_cpt'] = {
        'type': 'compact', 'method': 'roots',
        'root_fn': lambda: roots_D(5),
        'dim': 45, 'complexification': 'D5', 'physics': 'GUT SO(10)'}
    cat['G2'] = {
        'type': 'compact', 'method': 'roots',
        'root_fn': roots_G2,
        'dim': 14, 'complexification': 'G2', 'physics': 'Exceptional'}
    cat['F4'] = {
        'type': 'compact', 'method': 'roots',
        'root_fn': roots_F4,
        'dim': 52, 'complexification': 'F4', 'physics': 'Exceptional'}
    cat['E6'] = {
        'type': 'compact', 'method': 'roots',
        'root_fn': roots_E6,
        'dim': 78, 'complexification': 'E6', 'physics': 'GUT E6'}

    # ── Non-compact: Lorentz and spacetime ──
    cat['so(3,1)'] = {
        'type': 'non-compact', 'method': 'generators',
        'gen_fn': lambda: generators_so_pq(3, 1),
        'dim': 6, 'complexification': 'D2=A1+A1',
        'physics': 'Lorentz group',
        'compact_partner': 'so(4)'}
    cat['so(2,1)'] = {
        'type': 'non-compact', 'method': 'generators',
        'gen_fn': lambda: generators_so_pq(2, 1),
        'dim': 3, 'complexification': 'A1',
        'physics': '2+1 Lorentz',
        'compact_partner': 'so(3)'}
    cat['so(4,1)'] = {
        'type': 'non-compact', 'method': 'generators',
        'gen_fn': lambda: generators_so_pq(4, 1),
        'dim': 10, 'complexification': 'B2',
        'physics': 'de Sitter',
        'compact_partner': 'so(5)'}
    cat['so(3,2)'] = {
        'type': 'non-compact', 'method': 'generators',
        'gen_fn': lambda: generators_so_pq(3, 2),
        'dim': 10, 'complexification': 'B2',
        'physics': 'Anti-de Sitter / Conformal 2+1',
        'compact_partner': 'so(5)'}
    cat['so(4,2)'] = {
        'type': 'non-compact', 'method': 'generators',
        'gen_fn': lambda: generators_so_pq(4, 2),
        'dim': 15, 'complexification': 'D3=A3',
        'physics': 'Conformal 3+1',
        'compact_partner': 'so(6)'}
    cat['so(5,5)'] = {
        'type': 'non-compact', 'method': 'generators',
        'gen_fn': lambda: generators_so_pq(5, 5),
        'dim': 45, 'complexification': 'D5',
        'physics': 'Split D5',
        'compact_partner': 'so(10)_cpt'}
    cat['so(9,1)'] = {
        'type': 'non-compact', 'method': 'generators',
        'gen_fn': lambda: generators_so_pq(9, 1),
        'dim': 45, 'complexification': 'D5',
        'physics': '10D Lorentz',
        'compact_partner': 'so(10)_cpt'}

    # ── Non-compact: sl and su variants ──
    cat['sl(2,R)'] = {
        'type': 'non-compact', 'method': 'generators',
        'gen_fn': lambda: generators_sl_n_R(2),
        'dim': 3, 'complexification': 'A1',
        'physics': 'Simplest non-compact',
        'compact_partner': 'su(2)'}
    cat['sl(3,R)'] = {
        'type': 'non-compact', 'method': 'generators',
        'gen_fn': lambda: generators_sl_n_R(3),
        'dim': 8, 'complexification': 'A2',
        'physics': 'Split A2',
        'compact_partner': 'su(3)'}
    cat['su(2,1)'] = {
        'type': 'non-compact', 'method': 'generators',
        'gen_fn': lambda: generators_su_pq(2, 1),
        'dim': 8, 'complexification': 'A2',
        'physics': 'Pseudo-unitary',
        'compact_partner': 'su(3)'}
    cat['sl(4,R)'] = {
        'type': 'non-compact', 'method': 'generators',
        'gen_fn': lambda: generators_sl_n_R(4),
        'dim': 15, 'complexification': 'A3',
        'physics': 'Split A3',
        'compact_partner': 'su(4)'}
    cat['su(3,1)'] = {
        'type': 'non-compact', 'method': 'generators',
        'gen_fn': lambda: generators_su_pq(3, 1),
        'dim': 15, 'complexification': 'A3',
        'physics': 'Pseudo-unitary (3,1)',
        'compact_partner': 'su(4)'}

    return cat


# ─── Hashimoto D4 ────────────────────────────────────────────────────

def compute_d4_hashimoto(adj):
    """Compute D4 (Hashimoto) bound for a single adjacency matrix."""
    n = adj.shape[0]
    if n > HASH_DIM_LIMIT:
        return None
    edges = []
    for i in range(n):
        for j in range(n):
            if adj[i, j] > 0.5:
                edges.append((i, j))
    m = len(edges)
    if m == 0:
        return 0.0
    T_e = np.zeros((m, m), dtype=np.float64)
    edge_index = {e: idx for idx, e in enumerate(edges)}
    for idx_e, (i, j) in enumerate(edges):
        for k in range(n):
            if k != i and adj[j, k] > 0.5:
                jk = edge_index.get((j, k))
                if jk is not None:
                    T_e[jk, idx_e] = 1.0
    hash_eigs = np.linalg.eigvals(T_e)
    return float(np.max(np.abs(hash_eigs)))


# ─── 6-definition analysis ───────────────────────────────────────────

def compute_6defs(adj, rho_hash=None):
    """Compute Ramanujan under 6 definitions (reused logic from P6)."""
    n = adj.shape[0]
    eigs = np.sort(np.linalg.eigvalsh(adj))[::-1]
    degrees = adj.sum(axis=1)

    d_mean = float(np.mean(degrees))
    d_max = float(np.max(degrees))
    is_regular = np.allclose(degrees, degrees[0])

    nontrivial = eigs[1:]
    if n > 2 and is_regular and abs(eigs[-1] + eigs[0]) < 0.01:
        nontrivial = eigs[1:-1]
    max_nt = float(np.max(np.abs(nontrivial))) if len(nontrivial) > 0 else 0.0

    nonzero_d = degrees[degrees > 0]
    d_harm = float(len(nonzero_d) / np.sum(1.0 / nonzero_d)) if len(nonzero_d) > 0 else 0.0
    dm1 = degrees - 1
    dm1_pos = dm1[dm1 > 0]
    d_geom_m1 = float(np.exp(np.mean(np.log(dm1_pos)))) if len(dm1_pos) == n else max(d_mean - 1, 0)
    sorted_deg = np.sort(degrees)[::-1]
    d_2nd = float(sorted_deg[1]) if n > 1 else d_mean

    defs = {}
    for label, q_eff in [('D1', d_mean - 1), ('D2', d_max - 1),
                          ('D3', d_geom_m1), ('D5', d_harm - 1),
                          ('D6', d_2nd - 1)]:
        bound = 2 * np.sqrt(max(q_eff, 0))
        is_ram = bool(max_nt <= bound + 1e-10)
        ratio = max_nt / bound if bound > 0 else float('inf')
        defs[label] = {'q': q_eff, 'bound': bound, 'is_ram': is_ram, 'ratio': ratio}

    if rho_hash is not None:
        bound4 = 2 * np.sqrt(max(rho_hash, 0))
        is_ram4 = bool(max_nt <= bound4 + 1e-10)
        ratio4 = max_nt / bound4 if bound4 > 0 else float('inf')
        defs['D4'] = {'q': rho_hash, 'bound': bound4, 'is_ram': is_ram4, 'ratio': ratio4}
    else:
        defs['D4'] = None

    return {
        'max_nt': max_nt, 'd_mean': d_mean, 'd_max': d_max,
        'is_regular': is_regular, 'eigs': eigs, 'degrees': degrees,
        'definitions': defs,
    }


# ─── Graph isomorphism test (via sorted spectrum) ────────────────────

def graphs_isospectral(adj1, adj2):
    """Check if two graphs have the same spectrum (necessary for isomorphism)."""
    if adj1.shape != adj2.shape:
        return False
    eigs1 = np.sort(np.linalg.eigvalsh(adj1))
    eigs2 = np.sort(np.linalg.eigvalsh(adj2))
    return np.allclose(eigs1, eigs2, atol=1e-8)


def graphs_isomorphic_degree_seq(adj1, adj2):
    """Check degree sequence match (necessary for isomorphism)."""
    if adj1.shape != adj2.shape:
        return False
    d1 = np.sort(adj1.sum(axis=1))
    d2 = np.sort(adj2.sum(axis=1))
    return np.allclose(d1, d2)


# ─── Visualization ───────────────────────────────────────────────────

def create_figure(results, pairs, output_path):
    """6-panel figure for P8."""
    fig, axes = plt.subplots(2, 3, figsize=(19, 12))
    fig.suptitle('P8: Ramanujan Analysis of Non-Compact Lie Algebras',
                 fontsize=14, fontweight='bold')

    names = list(results.keys())
    compact = [n for n in names if results[n]['type'] == 'compact']
    noncompact = [n for n in names if results[n]['type'] == 'non-compact']

    # ── Panel 1: D1 ratio comparison (compact vs non-compact) ──
    ax = axes[0, 0]
    pair_names = list(pairs.keys())
    x = np.arange(len(pair_names))
    c_ratios = [pairs[p]['compact_ratio'] for p in pair_names]
    nc_ratios = [pairs[p]['noncompact_ratio'] for p in pair_names]
    w = 0.35
    ax.bar(x - w / 2, c_ratios, w, label='Compact', color='#1565C0')
    ax.bar(x + w / 2, nc_ratios, w, label='Non-compact', color='#EF5350')
    ax.axhline(1.0, color='black', ls='--', lw=1, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(pair_names, rotation=45, ha='right', fontsize=7)
    ax.set_ylabel('D1 ratio')
    ax.set_title('Compact vs Non-Compact: D1 Ratio')
    ax.legend(fontsize=8)

    # ── Panel 2: Isomorphism test results ──
    ax = axes[0, 1]
    iso_labels = pair_names
    iso_spec = [1 if pairs[p]['isospectral'] else 0 for p in pair_names]
    iso_deg = [1 if pairs[p]['same_degree_seq'] else 0 for p in pair_names]
    mat = np.array([iso_deg, iso_spec]).T
    from matplotlib.colors import ListedColormap
    cmap2 = ListedColormap(['#EF5350', '#66BB6A'])
    ax.imshow(mat, aspect='auto', cmap=cmap2, vmin=0, vmax=1,
              interpolation='nearest')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Same degree\nsequence', 'Isospectral'], fontsize=9)
    ax.set_yticks(range(len(iso_labels)))
    ax.set_yticklabels(iso_labels, fontsize=7)
    ax.set_title('Graph Isomorphism Tests')
    for i in range(len(iso_labels)):
        for j in range(2):
            ax.text(j, i, 'YES' if mat[i, j] else 'NO',
                    ha='center', va='center', fontsize=8, fontweight='bold')

    # ── Panel 3: Killing form signature ──
    ax = axes[0, 2]
    all_names = compact + noncompact
    sigs = []
    for n in all_names:
        r = results[n]
        p, neg, z = r.get('killing_sig', (0, 0, 0))
        sigs.append((p, neg))
    pos_vals = [s[0] for s in sigs]
    neg_vals = [s[1] for s in sigs]
    y = np.arange(len(all_names))
    ax.barh(y, pos_vals, height=0.4, label='+', color='#66BB6A', align='center')
    ax.barh(y + 0.4, [-v for v in neg_vals], height=0.4, label='-', color='#EF5350',
            align='center')
    ax.set_yticks(y + 0.2)
    ax.set_yticklabels(all_names, fontsize=6)
    ax.set_xlabel('# eigenvalues')
    ax.set_title('Killing Form Signature (+/-)')
    ax.legend(fontsize=8)

    # ── Panel 4: Ramanujan heatmap (6 defs) ──
    ax = axes[1, 0]
    def_keys = ['D1', 'D2', 'D3', 'D4', 'D5', 'D6']
    mat_ram = np.full((len(all_names), 6), 0.5)
    for i, n in enumerate(all_names):
        for j, dk in enumerate(def_keys):
            d = results[n]['defs']['definitions'].get(dk)
            if d is None:
                mat_ram[i, j] = 0.5
            elif d['is_ram']:
                mat_ram[i, j] = 1.0
            else:
                mat_ram[i, j] = 0.0
    cmap3 = ListedColormap(['#EF5350', '#BDBDBD', '#66BB6A'])
    ax.imshow(mat_ram, aspect='auto', cmap=cmap3, vmin=0, vmax=1,
              interpolation='nearest')
    ax.set_xticks(range(6))
    ax.set_xticklabels(def_keys, fontsize=8)
    ax.set_yticks(range(len(all_names)))
    ax.set_yticklabels(all_names, fontsize=6)
    ax.set_title('Ramanujan Classification (6 defs)')
    # Horizontal line separating compact / non-compact
    ax.axhline(len(compact) - 0.5, color='yellow', lw=2, ls='--')

    # ── Panel 5: Dimension vs ratio scatter ──
    ax = axes[1, 1]
    for alg_type, color, marker in [('compact', '#1565C0', 'o'),
                                     ('non-compact', '#EF5350', 's')]:
        mask = [n for n in names if results[n]['type'] == alg_type]
        dims = [results[n]['dim'] for n in mask]
        rats = [results[n]['defs']['definitions']['D1']['ratio'] for n in mask]
        ax.scatter(dims, rats, c=color, marker=marker, s=80,
                   label=alg_type, edgecolors='black', zorder=3)
        for n, d, r in zip(mask, dims, rats):
            ax.annotate(n, (d, r), fontsize=5.5, alpha=0.7,
                        xytext=(3, 3), textcoords='offset points')
    ax.axhline(1.0, color='red', ls='--', lw=1.5)
    ax.set_xlabel('Dimension')
    ax.set_ylabel('D1 ratio')
    ax.set_title('Dimension vs Ramanujan Ratio')
    ax.legend(fontsize=8)

    # ── Panel 6: Spectral comparison for Lorentz vs so(4) ──
    ax = axes[1, 2]
    if 'so(4)' in results and 'so(3,1)' in results:
        eigs_c = results['so(4)']['defs']['eigs']
        eigs_nc = results['so(3,1)']['defs']['eigs']
        ax.stem(range(len(eigs_c)), eigs_c, linefmt='b-', markerfmt='bo',
                basefmt='b-', label='so(4) compact')
        ax.stem(range(len(eigs_nc)), eigs_nc + 0.05, linefmt='r-',
                markerfmt='rs', basefmt='r-', label='so(3,1) non-compact')
        ax.set_xlabel('Eigenvalue index')
        ax.set_ylabel('Eigenvalue')
        ax.set_title('Spectrum: so(4) vs so(3,1)')
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, 'No Lorentz pair', transform=ax.transAxes,
                ha='center')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved: {output_path}")


# ─── Main analysis ───────────────────────────────────────────────────

def run_analysis():
    catalog = build_catalog()
    results = OrderedDict()

    print("=" * 90)
    print("P8: RAMANUJAN ANALYSIS OF NON-COMPACT LIE ALGEBRAS")
    print("=" * 90)

    # ── Phase 1: Build all adjacency graphs ──
    print("\n--- Phase 1: Adjacency Graph Construction ---")
    hdr = (f"{'Name':15s} {'Type':12s} {'Dim':>4s} {'d_mean':>6s} "
           f"|lam_nt|{'':<2s} {'bound':>6s} {'ratio':>6s} {'D1':>4s} "
           f"{'Cmplx':>12s}")
    print(hdr)
    print("-" * len(hdr))

    for name, info in catalog.items():
        if info['method'] == 'generators':
            gens = info['gen_fn']()
            adj = adjacency_from_generators(gens)
            K = killing_form(gens)
            k_sig = killing_signature(K)
        else:
            rank, roots, simple = info['root_fn']()
            adj = adjacency_from_roots(rank, roots, simple)
            k_sig = (0, info['dim'], 0)  # Compact: all negative

        dim_actual = adj.shape[0]
        sa = spectral_analysis(adj)
        sa['is_ramanujan'] = bool(sa['is_ramanujan'])
        rho_hash = compute_d4_hashimoto(adj)
        defs = compute_6defs(adj, rho_hash)

        results[name] = {
            'adj': adj, 'dim': dim_actual, 'type': info['type'],
            'complexification': info['complexification'],
            'physics': info['physics'],
            'spec': sa, 'defs': defs, 'killing_sig': k_sig,
            'rho_hash': rho_hash,
        }

        d1 = defs['definitions']['D1']
        tag = 'RAM' if d1['is_ram'] else 'non-Ram'
        ct = 'C' if info['type'] == 'compact' else 'NC'
        print(f"{name:15s} {ct:12s} {dim_actual:4d} {defs['d_mean']:6.2f} "
              f"{defs['max_nt']:8.4f} {d1['bound']:6.4f} {d1['ratio']:6.4f} "
              f"{tag:>4s} {info['complexification']:>12s}")

    # ── Phase 2: Compact/Non-compact pair comparison ──
    print("\n--- Phase 2: Compact vs Non-Compact Pair Comparison ---")
    pairs = OrderedDict()
    for name, info in catalog.items():
        if info['type'] == 'non-compact' and 'compact_partner' in info:
            partner = info['compact_partner']
            if partner in results:
                adj_c = results[partner]['adj']
                adj_nc = results[name]['adj']
                iso_spec = graphs_isospectral(adj_c, adj_nc)
                iso_deg = graphs_isomorphic_degree_seq(adj_c, adj_nc)
                c_d1 = results[partner]['defs']['definitions']['D1']
                nc_d1 = results[name]['defs']['definitions']['D1']
                label = f"{name} vs {partner}"
                pairs[label] = {
                    'compact': partner, 'noncompact': name,
                    'isospectral': iso_spec,
                    'same_degree_seq': iso_deg,
                    'compact_ratio': c_d1['ratio'],
                    'noncompact_ratio': nc_d1['ratio'],
                    'compact_ram': c_d1['is_ram'],
                    'noncompact_ram': nc_d1['is_ram'],
                }

    print(f"\n{'Pair':45s} {'IsoSpec':>8s} {'DegSeq':>8s} "
          f"{'C ratio':>8s} {'NC ratio':>9s} {'Same?':>6s}")
    print("-" * 90)
    n_identical = 0
    for label, p in pairs.items():
        same = 'YES' if p['isospectral'] else 'NO'
        deg_same = 'YES' if p['same_degree_seq'] else 'NO'
        ratio_match = 'YES' if abs(p['compact_ratio'] - p['noncompact_ratio']) < 1e-6 else 'NO'
        print(f"{label:45s} {same:>8s} {deg_same:>8s} "
              f"{p['compact_ratio']:8.4f} {p['noncompact_ratio']:9.4f} {ratio_match:>6s}")
        if p['isospectral']:
            n_identical += 1

    print(f"\nIsospectral pairs: {n_identical}/{len(pairs)}")

    # ── Phase 3: Killing form analysis ──
    print("\n--- Phase 3: Killing Form Signatures ---")
    print(f"{'Name':15s} {'Type':5s} {'Dim':>4s} "
          f"{'(+)':>4s} {'(-)':>4s} {'(0)':>4s} {'Definite?':>10s}")
    print("-" * 55)
    for name, res in results.items():
        ct = 'C' if res['type'] == 'compact' else 'NC'
        p, neg, z = res['killing_sig']
        if p == 0 and z == 0:
            definite = 'neg-def'
        elif neg == 0 and z == 0:
            definite = 'pos-def'
        elif z == res['dim']:
            definite = 'zero'
        else:
            definite = 'INDEFINITE'
        print(f"{name:15s} {ct:5s} {res['dim']:4d} "
              f"{p:4d} {neg:4d} {z:4d} {definite:>10s}")

    # ── Phase 4: Full 6-definition table ──
    print("\n--- Phase 4: Full 6-Definition Table ---")
    def_keys = ['D1', 'D2', 'D3', 'D4', 'D5', 'D6']
    print(f"{'Name':15s} {'Type':5s} ", end='')
    for dk in def_keys:
        print(f"  {dk:>4s}", end='')
    print(f"  {'#Ram':>5s}")
    print("-" * 65)
    for name, res in results.items():
        ct = 'C' if res['type'] == 'compact' else 'NC'
        print(f"{name:15s} {ct:5s} ", end='')
        n_ram = 0
        n_tested = 0
        for dk in def_keys:
            d = res['defs']['definitions'].get(dk)
            if d is None:
                print(f"  {'N/A':>4s}", end='')
            elif d['is_ram']:
                print(f"  {'Y':>4s}", end='')
                n_ram += 1
                n_tested += 1
            else:
                print(f"  {'N':>4s}", end='')
                n_tested += 1
        print(f"  {n_ram:2d}/{n_tested}")

    # ── Phase 5: Ramanujan preservation under real forms ──
    print("\n--- Phase 5: Does Ramanujan Depend on the Real Form? ---")
    form_groups = {}
    for name, res in results.items():
        cx = res['complexification']
        if cx not in form_groups:
            form_groups[cx] = []
        form_groups[cx].append(name)

    for cx, members in sorted(form_groups.items()):
        if len(members) < 2:
            continue
        d1_results = {}
        for m in members:
            d1 = results[m]['defs']['definitions']['D1']
            d1_results[m] = (d1['is_ram'], d1['ratio'])
        all_same = len(set(v[0] for v in d1_results.values())) == 1
        print(f"\n  Complexification {cx}:")
        for m, (ram, ratio) in d1_results.items():
            tag = 'RAM' if ram else 'non-Ram'
            ct = 'C' if results[m]['type'] == 'compact' else 'NC'
            print(f"    {m:15s} [{ct}] D1 ratio={ratio:.4f}  {tag}")
        verdict = "SAME classification" if all_same else "DIFFERENT classification"
        print(f"    => {verdict}")

    # ── Phase 6: SM spacetime symmetry ──
    print("\n--- Phase 6: Standard Model Full Symmetry ---")
    print("  SM internal:  su(3) + su(2) + u(1)   [from P7: D1 ratio = 0.657]")
    if 'so(3,1)' in results:
        lor = results['so(3,1)']['defs']['definitions']['D1']
        print(f"  SM spacetime: so(3,1)                 D1 ratio = {lor['ratio']:.4f}")
        print(f"  so(3,1) Ramanujan? {'YES' if lor['is_ram'] else 'NO'}")

    print("\n  Physically relevant non-compact spacetime algebras:")
    for name in ['so(3,1)', 'so(4,1)', 'so(3,2)', 'so(4,2)']:
        if name in results:
            d1 = results[name]['defs']['definitions']['D1']
            tag = 'RAM' if d1['is_ram'] else 'non-Ram'
            print(f"    {name:12s} ({results[name]['physics']:30s}) "
                  f"ratio={d1['ratio']:.4f}  {tag}")

    # ── Summary ──
    print("\n" + "=" * 90)
    print("SUMMARY: P8 KEY FINDINGS")
    print("=" * 90)

    print(f"\n1. Isospectral pairs (compact=non-compact graph): "
          f"{n_identical}/{len(pairs)}")

    # Count Ramanujan by type
    compact_names = [n for n in results if results[n]['type'] == 'compact']
    nc_names = [n for n in results if results[n]['type'] == 'non-compact']
    c_ram = sum(1 for n in compact_names
                if results[n]['defs']['definitions']['D1']['is_ram'])
    nc_ram = sum(1 for n in nc_names
                 if results[n]['defs']['definitions']['D1']['is_ram'])
    print(f"\n2. Compact Ramanujan (D1): {c_ram}/{len(compact_names)}")
    print(f"   Non-compact Ramanujan (D1): {nc_ram}/{len(nc_names)}")

    # Form independence
    form_indep = all(
        len(set(bool(results[m]['defs']['definitions']['D1']['is_ram'])
                for m in members)) <= 1
        for members in form_groups.values() if len(members) >= 2
    )
    print(f"\n3. Ramanujan classification independent of real form? "
          f"{'YES' if form_indep else 'NO'}")

    # Lorentz
    if 'so(3,1)' in results:
        lor = results['so(3,1)']['defs']['definitions']['D1']
        print(f"\n4. Lorentz algebra so(3,1): D1 ratio = {lor['ratio']:.4f}, "
              f"Ramanujan = {'YES' if lor['is_ram'] else 'NO'}")

    out = os.path.join(os.path.dirname(__file__), 'ramanujan_noncompact.png')
    create_figure(results, pairs, out)


if __name__ == '__main__':
    run_analysis()
