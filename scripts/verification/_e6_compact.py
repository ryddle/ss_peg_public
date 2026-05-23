"""
E6 compact form: the issue is that cubic preservation alone gives e6(6) (split).
To get e6(-78) (compact), we need generators that are ANTISYMMETRIC in the 27-dim rep.

Let's check: for F4, all generators ARE antisymmetric (derivations of J3(O) preserve the 
trace form Tr(X•Y), and derivations are antisymmetric in that metric).

For E6, we need: L preserves N(X) AND L is antisymmetric w.r.t. trace form:
  Tr(L(X)•Y) + Tr(X•L(Y)) = 0  ⟺  L^T + L = 0

The 52 negative-kappa generators already ARE antisymmetric.
The 26 positive-kappa generators are SYMMETRIC → they are in the split complement.
Removing the symmetric generators should give COMPACT e6(-78) if we can find 26 more
antisymmetric generators that preserve the cubic form.

Alternative: Instead of cubic norm preservation alone, add the ANTISYMMETRY constraint.
Antisym + cubic ⟹ compact e6(-78).

Let's try building e6 as: solve for L in gl(27) such that:
  1) L + L^T = 0  (antisymmetric)
  2) N(L·X, Y, Z) + N(X, L·Y, Z) + N(X, Y, L·Z) = 0  (cubic preservation)
"""
import sys, time
sys.path.insert(0, 'verification')
import numpy as np
from octonion_jordan import (
    cubic_norm_tensor, jordan_product_matrix, 
    _orthonormalize_frobenius
)
from shared_utils_gpu import (
    killing_metric_np, commutator_np, check_in_span_np,
    structure_constants_general_np
)

D = 27  # dimension of J3(O)

# Build cubic norm tensor
print("Building cubic norm tensor...")
N = cubic_norm_tensor()
print(f"  N shape: {N.shape}")

# Build ANTISYMMETRIC + CUBIC PRESERVING constraint
# Variable: L is a 27×27 antisymmetric matrix
# Antisymmetric: L[a,b] = -L[b,a], so we have D*(D-1)/2 = 351 free parameters
# Parameterize by upper triangle: L[a,b] for a < b

n_params = D*(D-1)//2  # 351
print(f"\nParameterizing antisymmetric 27×27 matrices: {n_params} free params")

# Cubic preservation constraint: Σ_d L[a,d] N[d,b,c] + L[b,d] N[a,d,c] + L[c,d] N[a,b,d] = 0
# For each (a,b,c) with a ≤ b ≤ c (symmetry of N), this gives one constraint.
# N is fully symmetric.

# Map (a,b) with a<b to parameter index
def param_index(a, b):
    """For a < b, return index in the 351-vector."""
    return a * D - a * (a + 1) // 2 + (b - a - 1)

# Number of unique (a,b,c) triples with a ≤ b ≤ c
n_constraints = 0
for a in range(D):
    for b in range(a, D):
        for c in range(b, D):
            n_constraints += 1
print(f"Number of cubic constraints: {n_constraints}")

# Build constraint matrix: each row corresponds to (a,b,c), each column to param (i,j) with i<j
print("Building antisym + cubic constraint matrix...")
t0 = time.time()

# For efficiency, build the full constraint as: for each L-parameter (i,j) with i<j,
# L[i,j] = +1, L[j,i] = -1 (antisymmetric basis element E_{ij})
# Then the cubic constraint for (a,b,c) from E_{ij} is:
#   (δ_{a,i} - δ_{a,j}...) → 
# Let me do it directly.

# Actually, for each basis element E_{ij} (antisym), L = E_ij - E_ji where i<j:
# L[i,j] = 1, L[j,i] = -1, rest 0
# Cubic: Σ_d L[a,d] N[d,b,c] + Σ_d L[b,d] N[a,d,c] + Σ_d L[c,d] N[a,b,d]
# = δ_{d=j}(δ_{a=i}) N[j,b,c] - δ_{d=i}(δ_{a=j}) N[i,b,c]
# Wait, L[a,d] = δ_{a,i}δ_{d,j} - δ_{a,j}δ_{d,i}
# So Σ_d L[a,d] N[d,b,c] = δ_{a,i} N[j,b,c] - δ_{a,j} N[i,b,c]

# The full constraint for E_{ij} on (a,b,c) is:
# δ_{a,i} N[j,b,c] - δ_{a,j} N[i,b,c] + δ_{b,i} N[a,j,c] - δ_{b,j} N[a,i,c] + δ_{c,i} N[a,b,j] - δ_{c,j} N[a,b,i]

# Since N is symmetric, N[a,b,i] = N[i,a,b] etc.

# Build constraint matrix
rows = []
row_idx = 0
for a in range(D):
    for b in range(a, D):
        for c in range(b, D):
            row = np.zeros(n_params)
            for p in range(n_params):
                # Decode parameter index back to (i,j) with i<j
                pass  # This is slow, let me use a different approach
            rows.append(row)
            row_idx += 1

# Actually this approach is too slow. Let me build it column by column instead.
constraint = np.zeros((n_constraints, n_params))

for p_idx in range(n_params):
    # Decode (i,j) from param_index
    # Find i,j such that param_index(i,j) = p_idx
    i = 0
    while param_index(i, i+1) + (D - i - 1) <= p_idx:
        i += 1
    j = p_idx - (i * D - i * (i + 1) // 2) + i + 1
    
    # For E_{ij} basis element, compute cubic constraint on each (a,b,c)
    row_idx = 0
    for a in range(D):
        for b in range(a, D):
            for c in range(b, D):
                val = 0.0
                # δ_{a,i} N[j,b,c] - δ_{a,j} N[i,b,c]
                if a == i: val += N[j,b,c]
                if a == j: val -= N[i,b,c]
                # δ_{b,i} N[a,j,c] - δ_{b,j} N[a,i,c]
                if b == i: val += N[a,j,c]
                if b == j: val -= N[a,i,c]
                # δ_{c,i} N[a,b,j] - δ_{c,j} N[a,b,i]
                if c == i: val += N[a,b,j]
                if c == j: val -= N[a,b,i]
                constraint[row_idx, p_idx] = val
                row_idx += 1

t1 = time.time()
print(f"  Built in {t1-t0:.1f}s, shape: {constraint.shape}")

# SVD to find null space
print("Computing SVD...")
t0 = time.time()
U, S, Vt = np.linalg.svd(constraint, full_matrices=True)
t1 = time.time()
print(f"  SVD in {t1-t0:.1f}s")

tol = 1e-10
rank = np.sum(S > tol)
null_dim = n_params - rank
print(f"  Rank: {rank}, Null dim: {null_dim}")
print(f"  Expected: 78 for e6 (if antisym+cubic gives all of e6)")
print(f"  Actually expected: f4 has dim 52, but f4 ⊂ antisym∩cubic")
print(f"  So we expect at least 52 (f4) if not more")

# Extract null space vectors and build generators
null_vecs = Vt[rank:]  # shape (null_dim, n_params)
print(f"\n  Building {null_dim} antisymmetric generators...")
gens = []
for k in range(null_dim):
    L = np.zeros((D, D))
    for p_idx in range(n_params):
        # Decode (i,j)
        i = 0
        while param_index(i, i+1) + (D - i - 1) <= p_idx:
            i += 1
        j = p_idx - (i * D - i * (i + 1) // 2) + i + 1
        L[i, j] = null_vecs[k, p_idx]
        L[j, i] = -null_vecs[k, p_idx]
    gens.append(L)

# Verify antisymmetry
max_asym = max(np.abs(g + g.T).max() for g in gens)
print(f"  Max antisymmetry error: {max_asym:.2e}")

# Verify trace-orthonormality
G = np.zeros((null_dim, null_dim))
for i in range(null_dim):
    for j in range(i, null_dim):
        G[i,j] = np.trace(gens[i] @ gens[j]).real
        G[j,i] = G[i,j]
print(f"  Trace metric diagonal: {np.diag(G)[:5]}")

# Orthonormalize w.r.t. Frobenius (= trace for antisymmetric since Tr(A^T B) = -Tr(AB))
gens_ortho = _orthonormalize_frobenius(gens)
print(f"  Orthonormalized: {len(gens_ortho)} generators")

# All should have Tr(T²) < 0 (antisymmetric)
kappas = [np.trace(g @ g).real for g in gens_ortho]
print(f"  kappa values: min={min(kappas):.6f}, max={max(kappas):.6f}")
print(f"  All negative: {all(k < -1e-8 for k in kappas)}")

# Closure check
closes = 0
total = null_dim*(null_dim-1)//2
for i in range(null_dim):
    for j in range(i+1, null_dim):
        C = commutator_np(gens_ortho[i], gens_ortho[j])
        if np.linalg.norm(C) < 1e-10 or check_in_span_np(C, gens_ortho, tol=1e-4):
            closes += 1
print(f"\n  Closure: {closes}/{total}")

# Killing metric
print("\n  Computing Killing metrics...")
K_bv = killing_metric_np(gens_ortho)
K_bv_diag = np.diag(K_bv).real
K_bv_eigs = np.linalg.eigvalsh(K_bv.real)
K_bv_offdiag = np.abs(K_bv - np.diag(K_bv_diag)).max()
print(f"  B_V diag: mean={np.mean(K_bv_diag):.6f}, std={np.std(K_bv_diag):.6f}")
print(f"  B_V off-diag max: {K_bv_offdiag:.2e}")
print(f"  B_V eigenvalues: min={K_bv_eigs[0]:.6f}, max={K_bv_eigs[-1]:.6f}")
print(f"  B_V negative definite: {all(e < -1e-8 for e in K_bv_eigs)}")

# h∨ from B_V
kappa = np.mean(kappas)
k_ratio = abs(np.mean(K_bv_diag)) / kappa**2
I_R = 3  # Dynkin index for 27-dim rep of E6
h_dual = I_R * k_ratio
print(f"\n  kappa = {kappa:.6f}")
print(f"  |B_V_diag|/kappa² = {k_ratio:.4f}")
print(f"  I_R = {I_R}")
print(f"  h∨ = I_R × k_ratio = {h_dual:.4f}  (expected: 12)")

# Structure constants and adjoint Killing
print("\n  Computing structure constants...")
c = structure_constants_general_np(gens_ortho)
max_jacobi = 0.0
for i in range(min(15, null_dim)):
    for j in range(i+1, min(15, null_dim)):
        for k in range(j+1, min(15, null_dim)):
            for m in range(null_dim):
                val = 0.0
                for l in range(null_dim):
                    val += c[i,j,l]*c[l,k,m] + c[j,k,l]*c[l,i,m] + c[k,i,l]*c[l,j,m]
                max_jacobi = max(max_jacobi, abs(val))
print(f"  Jacobi error: {max_jacobi:.2e}")

K_adj = np.einsum('iac,jca->ij', c, c)
K_adj_diag = np.diag(K_adj).real
K_adj_eigs = np.linalg.eigvalsh(K_adj.real)
print(f"  K_adj diag: mean={np.mean(K_adj_diag):.6f}")
print(f"  K_adj eigenvalues: min={K_adj_eigs[0]:.6f}, max={K_adj_eigs[-1]:.6f}")
print(f"  K_adj negative definite: {all(e < -1e-8 for e in K_adj_eigs)}")

print("\n=== DONE ===")
