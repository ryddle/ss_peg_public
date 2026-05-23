"""Deeper analysis of E6 generator structure."""
import sys
sys.path.insert(0, 'verification')
import numpy as np
from octonion_jordan import build_e6_generators, build_f4_generators
from shared_utils_gpu import killing_metric_np, commutator_np, check_in_span_np

print("=== E6 STRUCTURE ANALYSIS ===")
gens = build_e6_generators()

# Decompose each generator into antisym + sym parts
antisym_norms = []
sym_norms = []
traces = []
for g in gens:
    A = 0.5 * (g - g.T)  # antisymmetric
    S = 0.5 * (g + g.T)  # symmetric
    antisym_norms.append(np.linalg.norm(A))
    sym_norms.append(np.linalg.norm(S))
    traces.append(abs(np.trace(g)))

print(f"  Antisymmetric component norms: min={min(antisym_norms):.4f}, max={max(antisym_norms):.4f}, mean={np.mean(antisym_norms):.4f}")
print(f"  Symmetric component norms:     min={min(sym_norms):.4f}, max={max(sym_norms):.4f}, mean={np.mean(sym_norms):.4f}")
print(f"  Trace norms:                   max={max(traces):.2e}")

# How many are purely antisymmetric? Purely symmetric?
n_antisym = sum(1 for s in sym_norms if s < 1e-8)
n_sym = sum(1 for a in antisym_norms if a < 1e-8)
n_mixed = len(gens) - n_antisym - n_sym
print(f"  Purely antisymmetric: {n_antisym}")
print(f"  Purely symmetric:    {n_sym}")
print(f"  Mixed:               {n_mixed}")

# The structure of e6 in the 27-dim rep:
# e6 = f4 + J3(O)_0 as vector spaces
# where J3(O)_0 = traceless Jordan elements (26-dim)
# f4 generators: antisymmetric (derivations)
# J3(O)_0 generators: symmetric (Jordan multiplications L_a)
# Total: 52 + 26 = 78 ✓

# Compute Killing using structure_constants_general_np for proper treatment
from shared_utils_gpu import structure_constants_general_np

print("\n  Computing structure constants (slow for 78 gens)...")
import time
t0 = time.time()
c = structure_constants_general_np(gens)
t1 = time.time()
print(f"  Structure constants computed in {t1-t0:.1f}s")

# Killing from structure constants: K_ab = c_{acd} c_{bdc}
K_struct = np.einsum('acd,bdc->ab', c, c)
K_diag = np.diag(K_struct).real
K_offdiag = np.abs(K_struct - np.diag(K_diag)).max()
print(f"  K from struct consts:")
print(f"    K_diag: mean={np.mean(K_diag):.6f}, std={np.std(K_diag):.6f}")
print(f"    K_offdiag max: {K_offdiag:.2e}")
K_eigs = np.linalg.eigvalsh(K_struct.real)
print(f"    K eigenvalues: min={K_eigs[0]:.4f} max={K_eigs[-1]:.4f}")
# Is it negative definite? (should be for compact real form)
print(f"    Negative definite: {all(K_eigs < 0)}")

# Compare with Tr-based Killing
K_tr = killing_metric_np(gens)
print(f"\n  K from Tr(ad_a ad_b):")
print(f"    K_diag: mean={np.mean(np.diag(K_tr).real):.6f}")
print(f"    Diff from struct const K: {np.max(np.abs(K_tr - K_struct)):.2e}")

# kappa and h_dual with the proper approach
kappa = np.mean([np.trace(g @ g).real for g in gens])
k_ratio = abs(np.mean(K_diag)) / abs(kappa)**2 if abs(kappa) > 1e-12 else 0.0
print(f"\n  kappa (mean Tr(T^2)): {kappa:.6f}")
print(f"  |K_diag_mean|/kappa^2: {k_ratio:.4f}")
print(f"  For E6: h∨=12, so I_R = 12/{k_ratio:.4f} = {12.0/k_ratio:.4f}")
