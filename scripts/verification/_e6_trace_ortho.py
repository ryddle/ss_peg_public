"""
Construct compact E6 properly using Tits construction with CORRECT bracket structure.

e6(-78) = f4 ⊕ p26 where:
  - f4 = Der(J₃(𝕆)) acts on J₃(𝕆) by derivation: D(X)
  - p26 = J₃(𝕆)₀ (traceless) acts by: i*L_X where L_X(Y) = X∘Y - (1/3)Tr(X∘Y)I

The bracket structure:
  [D₁, D₂] = D₁D₂ - D₂D₁                     (derivation commutator in f4)
  [D, L_X] = L_{D(X)}                            (f4 acts on p by derivation)  
  [L_X, L_Y] = D_{X,Y} + L_{[X,Y]_J}           (Lie algebra of p)
      where D_{X,Y} is a specific derivation and [X,Y]_J is a Jordan-Lie bracket

Actually for COMPACT form, the generators of p must be ANTI-HERMITIAN w.r.t.
the Killing form. Since we're in a REAL representation, "anti-hermitian" means
antisymmetric w.r.t. some inner product.

Let's use a different approach entirely: construct the ADJOINT representation of e6
from the structure constants of the 78-dim algebra, which are well-defined by:
  [f4, f4] ⊂ f4    (the f4 structure constants)
  [f4, p] ⊂ p      (derivation action on Jordan elements)
  [p, p] ⊂ f4 + p  (Jordan-related bracket)

For [p, p]: [L_X, L_Y] = [L_X, L_Y]_f4 + [L_X, L_Y]_p
The f4 component: 2(L_{X∘Y} L_Z - L_{Y∘Z} L_X + L_Z L_{X∘Y} - ...) = 2[L_X, L_Y]
which involves the derivation D(X,Y)(Z) = [[L_X, L_Y], L_Z] projected onto f4.

Actually there's a much simpler approach: since we already know the generators
in the 27-dim rep close at 100%, we can compute the structure constants
PROPERLY using the representation trace form, which just requires orthogonalizing
the generators w.r.t. Tr(T_a @ T_b).

The key insight: structure_constants_general_np assumes 
  [T_a, T_b] = Σ_c c_{abc} T_c
  c_{abc} = Tr([T_a,T_b] T_c) / Tr(T_c²)
  
This requires Tr(T_a T_b) ~ δ_{ab} to work. Let's orthogonalize w.r.t. trace form.
"""
import sys, time
sys.path.insert(0, 'verification')
import numpy as np
from octonion_jordan import build_e6_generators, build_f4_generators
from shared_utils_gpu import (
    killing_metric_np, commutator_np, check_in_span_np,
    structure_constants_general_np
)

print("=== E6: TRACE-ORTHOGONALIZED BASIS ===")
gens_raw = build_e6_generators()
n = len(gens_raw)

# Compute trace metric: G[a,b] = Tr(T_a @ T_b)
print("\n  Computing trace metric Tr(T_a T_b)...")
G = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        G[i,j] = np.trace(gens_raw[i] @ gens_raw[j]).real
        G[j,i] = G[i,j]
G_eigs = np.linalg.eigvalsh(G)
print(f"  Trace metric eigenvalues: min={G_eigs[0]:.6f} max={G_eigs[-1]:.6f}")
print(f"  Negative: {sum(e < -1e-8 for e in G_eigs)}")
print(f"  Zero: {sum(abs(e) < 1e-8 for e in G_eigs)}")
print(f"  Positive: {sum(e > 1e-8 for e in G_eigs)}")

# If G is indefinite, generators have mixed signature in trace form.
# This is expected for e6(split) where f4-part has Tr(D²) < 0 (antisym)
# and p-part has Tr(L²) > 0 (symmetric, since Tr(L_X²) = Tr(X∘X) > 0).

# Let's orthogonalize: diagonalize G and normalize to |eigenvalue| = 1
evals, evecs = np.linalg.eigh(G)
print(f"\n  Diagonalizing trace metric...")
gens_ortho = []
signs = []
for a in range(n):
    g = sum(evecs[b, a] * gens_raw[b] for b in range(n))
    if abs(evals[a]) > 1e-10:
        g = g / np.sqrt(abs(evals[a]))
        gens_ortho.append(g)
        signs.append(1 if evals[a] > 0 else -1)

print(f"  Got {len(gens_ortho)} orthonormalized generators")
print(f"  Signs: {sum(s>0 for s in signs)} positive, {sum(s<0 for s in signs)} negative")

# Verify trace-orthogonality
G_new = np.zeros((len(gens_ortho), len(gens_ortho)))
for i in range(len(gens_ortho)):
    for j in range(i, len(gens_ortho)):
        G_new[i,j] = np.trace(gens_ortho[i] @ gens_ortho[j]).real
        G_new[j,i] = G_new[i,j]
G_new_offdiag = np.abs(G_new - np.diag(np.diag(G_new))).max()
print(f"  Off-diagonal of new trace metric: {G_new_offdiag:.2e}")
print(f"  Diagonal values (first 5): {[f'{G_new[i,i]:.4f}' for i in range(5)]}")

# Now compute structure constants — should work since Tr(T_a T_b) ∝ δ_{ab}
print("\n  Computing structure constants...")
t0 = time.time()
c = structure_constants_general_np(gens_ortho)
t1 = time.time()
print(f"  Computed in {t1-t0:.1f}s")

# Jacobi identity check
max_jacobi = 0.0
for i in range(min(15, n)):
    for j in range(i+1, min(15, n)):
        for k in range(j+1, min(15, n)):
            for m in range(n):
                val = 0.0
                for l in range(n):
                    val += c[i,j,l]*c[l,k,m] + c[j,k,l]*c[l,i,m] + c[k,i,l]*c[l,j,m]
                max_jacobi = max(max_jacobi, abs(val))
print(f"  Max Jacobi error (sample): {max_jacobi:.2e}")

# Closure check
closes = 0
total = n*(n-1)//2
for i in range(n):
    for j in range(i+1, n):
        C = commutator_np(gens_ortho[i], gens_ortho[j])
        norm_C = np.linalg.norm(C)
        if norm_C < 1e-10 or check_in_span_np(C, gens_ortho, tol=1e-4):
            closes += 1
print(f"  Closure: {closes}/{total}")

# Killing form from structure constants (TRUE adjoint Killing)
K_adj = np.einsum('iac,jca->ij', c, c)
K_adj_diag = np.diag(K_adj).real
K_adj_eigs = np.linalg.eigvalsh(K_adj.real)
K_adj_offdiag = np.abs(K_adj - np.diag(K_adj_diag)).max()
print(f"\n  TRUE Killing (adjoint):")
print(f"    Diag: mean={np.mean(K_adj_diag):.4f}, std={np.std(K_adj_diag):.4f}")
print(f"    Off-diag max: {K_adj_offdiag:.2e}")
print(f"    Eigenvalues: min={K_adj_eigs[0]:.4f} max={K_adj_eigs[-1]:.4f}")

# B_V from killing_metric_np
K_bv = killing_metric_np(gens_ortho)
K_bv_diag = np.diag(K_bv).real
K_bv_eigs = np.linalg.eigvalsh(K_bv.real)
K_bv_offdiag = np.abs(K_bv - np.diag(K_bv_diag)).max()
print(f"\n  B_V = sum_a Tr([T_i,T_a][T_j,T_a]):")
print(f"    Diag: mean={np.mean(K_bv_diag):.4f}, std={np.std(K_bv_diag):.4f}")
print(f"    Off-diag max: {K_bv_offdiag:.2e}")
print(f"    Eigenvalues: min={K_bv_eigs[0]:.4f} max={K_bv_eigs[-1]:.4f}")

# The relationship between B_V and K_adj:
# For a simple algebra, B_V = I_V * K_adj (proportional)
# Let's check the ratio
if abs(np.mean(K_adj_diag)) > 1e-6:
    ratio = np.mean(K_bv_diag) / np.mean(K_adj_diag)
    print(f"\n  B_V / K_adj ratio (diagonal means): {ratio:.6f}")

# h∨ calculation
kappa = np.mean([np.trace(g @ g).real for g in gens_ortho])
print(f"\n  kappa (mean Tr(T²)): {kappa:.6f}")

# Using B_V approach:
k_ratio_bv = abs(np.mean(K_bv_diag)) / abs(kappa)**2
print(f"  |B_V_diag| / kappa²: {k_ratio_bv:.4f}")

# For compact generators (f4-type, antisym): Tr(T²) < 0
# For non-compact generators (p-type, sym): Tr(T²) > 0
# The mean kappa is a mix.

# Let's separate and analyze by sign:
neg_idx = [i for i in range(n) if signs[i] < 0]
pos_idx = [i for i in range(n) if signs[i] > 0]
print(f"\n  Negative-kappa generators ({len(neg_idx)}):")
if neg_idx:
    kappa_neg = np.mean([np.trace(gens_ortho[i] @ gens_ortho[i]).real for i in neg_idx])
    K_neg = np.mean([K_bv_diag[i] for i in neg_idx])
    kr_neg = abs(K_neg) / abs(kappa_neg)**2
    print(f"    kappa: {kappa_neg:.6f}, |K_diag|/kappa²: {kr_neg:.4f}")

print(f"  Positive-kappa generators ({len(pos_idx)}):")
if pos_idx:
    kappa_pos = np.mean([np.trace(gens_ortho[i] @ gens_ortho[i]).real for i in pos_idx])
    K_pos = np.mean([K_bv_diag[i] for i in pos_idx])
    kr_pos = abs(K_pos) / abs(kappa_pos)**2
    print(f"    kappa: {kappa_pos:.6f}, |K_diag|/kappa²: {kr_pos:.4f}")
