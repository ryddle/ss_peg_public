"""
KEY INSIGHT: The compact form e6(-78) preserves a POSITIVE DEFINITE Hermitian form on J3(O).
The trace form Tr(X•Y) IS positive definite on J3(O), and f4 = Der(J3(O)) automatically 
preserves it. BUT the representation on R^27 uses a SPECIFIC basis for J3(O), and in that 
basis, the f4 generators may not be literally antisymmetric (T + T^T = 0) unless the basis 
is orthonormal w.r.t. the trace form.

Our basis for J3(O):
  e_0 = diag(1,0,0), e_1 = diag(0,1,0), e_2 = diag(0,0,1)   (diagonal ξ₁,ξ₂,ξ₃)
  e_{3..10} = off-diagonal x₁ components (octonionic)
  e_{11..18} = off-diagonal x₂ components
  e_{19..26} = off-diagonal x₃ components

The trace form on J3(O) is: Tr(X•Y) where X•Y is Jordan product.
For diagonal elements: Tr(e_i • e_j) = Tr(e_i*e_j+e_j*e_i)/2 for i,j < 3 → this is ½δ_{ij}*2
Actually Tr(e_i) for diag elements: Tr(diag(1,0,0)) = 1.

Let me compute the trace bilinear form explicitly, then change basis to orthonormalize.
"""
import sys, time
sys.path.insert(0, 'verification')
import numpy as np
from octonion_jordan import (
    build_f4_generators, build_e6_generators, cubic_norm_tensor,
    jordan_product_matrix, _orthonormalize_frobenius,
    verify_derivation, verify_cubic_preservation
)
from shared_utils_gpu import (
    killing_metric_np, commutator_np, check_in_span_np,
    structure_constants_general_np
)

D = 27
J = jordan_product_matrix()

# Trace bilinear form: Tr(X•Y) = Σ_c J[a,b,c] * δ_c (where δ_c picks out trace elements)
# Actually, Tr is the linear form that picks out the trace of the matrix.
# For J3(O) element with coords v[0..26]:
#   The matrix has diagonal entries ξ₁=v[0], ξ₂=v[1], ξ₃=v[2]
#   So Tr(X) = v[0] + v[1] + v[2]
# The trace bilinear form is: (X,Y) = Tr(X•Y)
# Tr(X•Y) = Σ_c T_c * (X•Y)_c where T = trace vector = (1,1,1,0,...,0)
#          = Σ_{a,b,c} T_c * J[a,b,c] * x_a * y_b

T = np.zeros(D)
T[0] = T[1] = T[2] = 1.0

# Bilinear form: B[a,b] = Σ_c T[c] * J[a,b,c]
B = np.einsum('c,abc->ab', T, J)
print("Trace bilinear form B[a,b] = Tr(e_a • e_b):")
print(f"  Shape: {B.shape}")
print(f"  Eigenvalues: {np.sort(np.linalg.eigvalsh(B))[:5]} ... {np.sort(np.linalg.eigvalsh(B))[-5:]}")
print(f"  Is symmetric: {np.allclose(B, B.T)}")
print(f"  Is positive definite: {all(np.linalg.eigvalsh(B) > 0)}")
print(f"  Diagonal: {np.diag(B)[:6]}...")

# Change basis to orthonormalize w.r.t. B
# B = P^T @ P (Cholesky if PD)
# New basis: e'_a = Σ_b P^{-1}_{ba} e_b → B' = I
evals, evecs = np.linalg.eigh(B)
print(f"\n  B eigenvalues: {evals}")
if all(evals > 0):
    sqrt_evals = np.sqrt(evals)
    # P = diag(sqrt_evals) @ evecs^T → B = P^T P
    # Change basis: v' = P @ v
    # Generator in new basis: T' = P @ T @ P^{-1}
    P = np.diag(sqrt_evals) @ evecs.T
    P_inv = evecs @ np.diag(1.0/sqrt_evals)
    
    print("\nChanging basis to trace-orthonormal frame...")
    
    # Transform f4 generators
    gens_f4_raw = build_f4_generators()
    gens_f4_new = [P @ g @ P_inv for g in gens_f4_raw]
    
    # Check antisymmetry in new basis
    asym_errors = [np.abs(g + g.T).max() for g in gens_f4_new]
    print(f"\n=== F4 in trace-orthonormal basis ===")
    print(f"  Max antisym error across all gens: {max(asym_errors):.2e}")
    print(f"  All antisymmetric: {all(e < 1e-10 for e in asym_errors)}")
    
    # Transform e6 generators
    gens_e6_raw = build_e6_generators()
    gens_e6_new = [P @ g @ P_inv for g in gens_e6_raw]
    
    asym_errors_e6 = [np.abs(g + g.T).max() for g in gens_e6_new]
    sym_norms = [np.linalg.norm(g + g.T) for g in gens_e6_new]
    asym_norms = [np.linalg.norm(g - g.T) for g in gens_e6_new]
    
    print(f"\n=== E6 in trace-orthonormal basis ===")
    print(f"  Antisym errors: min={min(asym_errors_e6):.2e}, max={max(asym_errors_e6):.2e}")
    n_antisym = sum(1 for e in asym_errors_e6 if e < 1e-10)
    n_sym = sum(1 for n in sym_norms if n > 1e-10)
    print(f"  Purely antisymmetric: {n_antisym}")
    print(f"  Have symmetric component: {n_sym}")
    
    # Decompose each into antisym + sym
    antisym_parts = [(g - g.T)/2 for g in gens_e6_new]
    sym_parts = [(g + g.T)/2 for g in gens_e6_new]
    asym_norms2 = [np.linalg.norm(a) for a in antisym_parts]
    sym_norms2 = [np.linalg.norm(s) for s in sym_parts]
    
    print(f"\n  Antisym part norms: min={min(asym_norms2):.4f}, max={max(asym_norms2):.4f}")
    print(f"  Sym part norms: min={min(sym_norms2):.4f}, max={max(sym_norms2):.4f}")
    
    # Now if we project onto purely antisymmetric subspace of the null space:
    # Parameterize antisymmetric matrices and intersect with cubic preservation
    # This should give the COMPACT form
    
    # Transform cubic norm tensor to new basis
    # N'_{abc} = Σ_{ijk} P^{-1}_{ia} P^{-1}_{jb} P^{-1}_{kc} N_{ijk}
    # = Σ_{ijk} (P_inv^T)_{ai} N_{ijk} (P_inv^T)_{bj} (P_inv^T)_{ck}
    N_raw = cubic_norm_tensor()
    N_new = np.einsum('ia,jb,kc,ijk->abc', P_inv, P_inv, P_inv, N_raw)
    
    print(f"\n  Transformed cubic tensor. Checking symmetry...")
    print(f"  |N - N_perm|_max: {np.abs(N_new - np.transpose(N_new, (1,0,2))).max():.2e}")
    
    # Now build antisym + cubic-preserving generators in the NEW basis
    n_params = D*(D-1)//2  # 351
    
    def param_index(a, b):
        return a * D - a * (a + 1) // 2 + (b - a - 1)
    
    # Count constraints
    triples = []
    for a in range(D):
        for b in range(a, D):
            for c in range(b, D):
                triples.append((a,b,c))
    n_constraints = len(triples)
    
    print(f"\n  Building constraint matrix (antisym + cubic in new basis)...")
    print(f"  {n_constraints} constraints × {n_params} params")
    
    constraint = np.zeros((n_constraints, n_params))
    for p_idx in range(n_params):
        i = 0
        while param_index(i, i+1) + (D - i - 1) <= p_idx:
            i += 1
        j = p_idx - (i * D - i * (i + 1) // 2) + i + 1
        
        for row_idx, (a, b, c) in enumerate(triples):
            val = 0.0
            if a == i: val += N_new[j,b,c]
            if a == j: val -= N_new[i,b,c]
            if b == i: val += N_new[a,j,c]
            if b == j: val -= N_new[a,i,c]
            if c == i: val += N_new[a,b,j]
            if c == j: val -= N_new[a,b,i]
            constraint[row_idx, p_idx] = val
    
    print("  Computing SVD...")
    U, S, Vt = np.linalg.svd(constraint, full_matrices=True)
    tol = 1e-10
    rank = np.sum(S > tol)
    null_dim = n_params - rank
    print(f"  Rank: {rank}, Null dim: {null_dim}")
    
    # Build generators
    null_vecs = Vt[rank:]
    gens_compact = []
    for k in range(null_dim):
        L = np.zeros((D, D))
        for p_idx in range(n_params):
            i = 0
            while param_index(i, i+1) + (D - i - 1) <= p_idx:
                i += 1
            j = p_idx - (i * D - i * (i + 1) // 2) + i + 1
            L[i, j] = null_vecs[k, p_idx]
            L[j, i] = -null_vecs[k, p_idx]
        gens_compact.append(L)
    
    # Orthonormalize
    gens_compact = _orthonormalize_frobenius(gens_compact)
    print(f"  Orthonormalized: {len(gens_compact)} generators")
    
    # Verify antisymmetry
    max_asym = max(np.abs(g + g.T).max() for g in gens_compact)
    print(f"  Antisymmetry error: {max_asym:.2e}")
    
    # All kappas should be negative
    kappas = [np.trace(g @ g).real for g in gens_compact]
    print(f"  kappa: min={min(kappas):.6f}, max={max(kappas):.6f}")
    
    # Closure
    closes = 0
    total = null_dim*(null_dim-1)//2
    for i in range(null_dim):
        for j in range(i+1, null_dim):
            C = commutator_np(gens_compact[i], gens_compact[j])
            if np.linalg.norm(C) < 1e-10 or check_in_span_np(C, gens_compact, tol=1e-4):
                closes += 1
    print(f"\n  Closure: {closes}/{total}")
    
    # Killing
    K_bv = killing_metric_np(gens_compact)
    K_bv_diag = np.diag(K_bv).real
    K_bv_eigs = np.linalg.eigvalsh(K_bv.real)
    print(f"  B_V diag: mean={np.mean(K_bv_diag):.6f}")
    print(f"  B_V eigenvalues: min={K_bv_eigs[0]:.6f}, max={K_bv_eigs[-1]:.6f}")
    print(f"  Negative definite: {all(e < -1e-8 for e in K_bv_eigs)}")
    
    k_ratio = abs(np.mean(K_bv_diag)) / np.mean(kappas)**2
    print(f"  |B_V|/kappa² = {k_ratio:.4f}")
    
    # Structure constants + Jacobi
    print("\n  Structure constants...")
    c_sc = structure_constants_general_np(gens_compact)
    max_jacobi = 0.0
    for i2 in range(min(15, null_dim)):
        for j2 in range(i2+1, min(15, null_dim)):
            for k2 in range(j2+1, min(15, null_dim)):
                for m in range(null_dim):
                    val = 0.0
                    for l in range(null_dim):
                        val += c_sc[i2,j2,l]*c_sc[l,k2,m] + c_sc[j2,k2,l]*c_sc[l,i2,m] + c_sc[k2,i2,l]*c_sc[l,j2,m]
                    max_jacobi = max(max_jacobi, abs(val))
    print(f"  Jacobi error: {max_jacobi:.2e}")
    
    K_adj = np.einsum('iac,jca->ij', c_sc, c_sc)
    K_adj_eigs = np.linalg.eigvalsh(K_adj.real)
    print(f"  K_adj eigenvalues: min={K_adj_eigs[0]:.6f}, max={K_adj_eigs[-1]:.6f}")
    print(f"  K_adj negative definite: {all(e < -1e-8 for e in K_adj_eigs)}")
    
else:
    print("B is NOT positive definite — cannot use Cholesky")

print("\n=== DONE ===")
