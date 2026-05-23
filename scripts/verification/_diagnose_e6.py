"""Diagnose E6 k_ratio and f4 subalgebra embedding."""
import sys, time
sys.path.insert(0, 'verification')
import numpy as np
from octonion_jordan import build_e6_generators, build_f4_generators
from shared_utils_gpu import killing_metric_np, orthonormalize_killing_np, commutator_np, check_in_span_np

print("=== E6 DIAGNOSIS ===")
gens_e6 = build_e6_generators()
n = len(gens_e6)

# Killing-orthonormalize
print("\n  Killing-orthonormalizing E6...")
gens_e6_k = orthonormalize_killing_np(gens_e6)
K = killing_metric_np(gens_e6_k)
K_diag = np.diag(K).real
K_offdiag = np.abs(K - np.diag(K_diag)).max()
K_diag_mean = np.mean(K_diag)
K_diag_std = np.std(K_diag)
print(f"  K_diag: {K_diag_mean:.6f} +/- {K_diag_std:.6f}")
print(f"  K_offdiag: {K_offdiag:.2e}")
K_eigs = np.linalg.eigvalsh(K.real)
print(f"  K eigenvalues: min={K_eigs[0]:.4f} max={K_eigs[-1]:.4f}")

kappa = np.mean([np.trace(g @ g).real for g in gens_e6_k])
k_ratio = abs(K_diag_mean) / abs(kappa)**2
print(f"  kappa = {kappa:.6f}")
print(f"  |K|/kappa^2 = {k_ratio:.4f}")

# Check closure after K-ortho
closes = 0
total = n * (n-1) // 2
for i in range(n):
    for j in range(i+1, n):
        C = commutator_np(gens_e6_k[i], gens_e6_k[j])
        norm_C = np.linalg.norm(C)
        if norm_C < 1e-10 or check_in_span_np(C, gens_e6_k, tol=1e-4):
            closes += 1
print(f"  Closure: {closes}/{total} ({closes/total:.1%})")

# Check f4 in e6 using projection residual instead of check_in_span
print("\n  F4 subalgebra embedding check (projection method)...")
gens_f4 = build_f4_generators()

# Build E6 span matrix for projection
E6_basis = np.array([g.flatten() for g in gens_e6])  # (78, 729)
# For each F4 generator, project onto E6 span and check residual
f4_in_e6 = 0
max_residual = 0.0
for i, f4_g in enumerate(gens_f4):
    v = f4_g.flatten()
    # Project: coeffs = E6_basis @ v / (E6_basis @ E6_basis^T) — use lstsq
    coeffs, _, _, _ = np.linalg.lstsq(E6_basis.T, v, rcond=None)
    proj = E6_basis.T @ coeffs
    residual = np.linalg.norm(v - proj) / np.linalg.norm(v)
    max_residual = max(max_residual, residual)
    if residual < 1e-8:
        f4_in_e6 += 1
    elif i < 5:
        print(f"    F4 gen {i}: residual = {residual:.2e}")

print(f"  f4 generators in span(e6): {f4_in_e6}/{len(gens_f4)}")
print(f"  Max relative residual: {max_residual:.2e}")

# Also verify: f4 closure within e6 span
# Take f4 gens projected onto e6, check if their commutators still close
print("\n  Verifying f4 commutator closure within e6 span...")
f4_projected = []
for f4_g in gens_f4:
    v = f4_g.flatten()
    coeffs, _, _, _ = np.linalg.lstsq(E6_basis.T, v, rcond=None)
    proj = (E6_basis.T @ coeffs).reshape(27, 27)
    f4_projected.append(proj)

closes_f4 = 0
total_f4 = len(f4_projected) * (len(f4_projected)-1) // 2
for i in range(len(f4_projected)):
    for j in range(i+1, len(f4_projected)):
        C = commutator_np(f4_projected[i], f4_projected[j])
        norm_C = np.linalg.norm(C)
        if norm_C < 1e-10 or check_in_span_np(C, f4_projected, tol=1e-4):
            closes_f4 += 1
print(f"  f4 projected closure: {closes_f4}/{total_f4} ({closes_f4/total_f4:.1%})")
