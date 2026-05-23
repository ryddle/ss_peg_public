"""Fix E6 construction: add antisymmetry constraint for compact real form."""
import sys, time
sys.path.insert(0, 'verification')
import numpy as np
from octonion_jordan import cubic_norm_tensor, _orthonormalize_frobenius
from shared_utils_gpu import killing_metric_np, commutator_np, check_in_span_np, structure_constants_general_np

print("=== E6 COMPACT FORM: Cubic + Antisymmetry ===")

# Method: intersect cubic-preservation null space with antisymmetric matrices
# so(27) has dim 351. Cubic preservation reduces 729 -> 78.
# Antisymmetry: D + D^T = 0, picks 351 out of 729.
# Intersection should give compact e6(-78) = 78 generators (all antisymmetric).

# Actually no: e6(-78) has dimension 78 but only 36 generators are antisymmetric
# (maximal compact subalgebra = f4, dim=52... wait, f4 has 52 gens all antisymmetric
# and they are derivations. E6 compact is 78-dim with ALL Killing < 0.)

# Let me reconsider. The compact form e6(-78) acts on the COMPLEX 27-dim space.
# In the REAL 27-dim J3(O) representation, e6 generators split into:
#   - f4 part (52-dim): antisymmetric w.r.t. trace form
#   - p part (26-dim): symmetric w.r.t. trace form
# Together they form the compact e6(-78) with Cartan decomposition e6 = f4 + p.

# The Killing form on the compact form should still be negative definite.
# The problem is that our trace-based Killing computation is wrong for mixed
# antisym+sym generators.

# Let me try: restrict to antisymmetric part of cubic null space.
print("  Computing cubic norm tensor N_{abc}...")
N = cubic_norm_tensor()

print("  Building constraint: cubic preservation + antisymmetry...")
dim = 27
n_vars = dim * (dim - 1) // 2  # 351 parameters for antisymmetric matrix

# Parametrize L as antisymmetric: L[i,j] = a[idx(i,j)] for i<j, L[j,i] = -a[idx(i,j)]
def idx(i, j):
    """Index for (i,j) with i<j in upper triangle."""
    return i * dim - i * (i + 1) // 2 + (j - i - 1)

rows = []
for a in range(dim):
    for b in range(a, dim):
        for c in range(b, dim):
            row = np.zeros(n_vars)
            for d in range(dim):
                # L[d,a] * N[d,b,c] + L[d,b] * N[a,d,c] + L[d,c] * N[a,b,d]
                # L[d,a]: if d<a, L[d,a]=param[idx(d,a)]; if d>a, L[d,a]=-param[idx(a,d)]; if d==a, L[d,a]=0
                def add_L(row, d, x, coeff):
                    """Add coeff * L[d,x] to row."""
                    if d == x:
                        return
                    if d < x:
                        row[idx(d, x)] += coeff
                    else:
                        row[idx(x, d)] -= coeff
                
                if abs(N[d, b, c]) > 1e-15:
                    add_L(row, d, a, N[d, b, c])
                if abs(N[a, d, c]) > 1e-15:
                    add_L(row, d, b, N[a, d, c])
                if abs(N[a, b, d]) > 1e-15:
                    add_L(row, d, c, N[a, b, d])
            if np.any(np.abs(row) > 1e-15):
                rows.append(row)

M = np.array(rows)
print(f"  Constraint matrix: {M.shape[0]} x {M.shape[1]}")

_, S, Vt = np.linalg.svd(M, full_matrices=True)
rank = np.sum(S > 1e-8)
null_dim = n_vars - rank
print(f"  SVD: rank={rank}, null space dim={null_dim}")
print(f"  Expected: 52 (= f4, since these are antisym + cubic-preserving)")

# These should be exactly f4 generators!
# Because: antisymmetric + cubic-preserving = trace-form-preserving + cubic-preserving
# = derivations of J3(O) = f4.

# So the compact e6 CANNOT be obtained just by intersecting with so(27).
# For the full compact e6, we need a DIFFERENT approach.

# CORRECT APPROACH: Use the Tits construction or the explicit Cartan decomposition.
#
# e6(-78) = f4 + p, where p ≅ J3(O)_0 (traceless Jordan elements)
# as f4-modules.
#
# p is NOT antisymmetric in the 27-dim rep. Instead, each X in J3(O)_0
# acts on J3(O) via L_X(Y) = X ∘ Y - (1/3)Tr(X∘Y)I
# (traceless Jordan multiplication).
#
# The bracket [p, p] ⊂ f4 is: [L_X, L_Y] = [X,Y]_f4 = L_{X∘Y - Y∘X} + D_{X,Y}
# Actually, [L_X, L_Y](Z) = ... this gives the commutator structure.
#
# Let me construct e6 explicitly:
# Basis = {D_1,...,D_52} (f4 derivations) ∪ {L_X1,...,L_X26} (traceless L_X)

print("\n=== CONSTRUCTING E6 VIA TITS: f4 + L_X ===")

# Get f4 generators (derivations)
from octonion_jordan import build_f4_generators, jordan_product_matrix, j3_basis
f4_gens = build_f4_generators()
print(f"  f4 generators: {len(f4_gens)}")

# Get Jordan product structure constants
J = jordan_product_matrix()

# Build J3(O)_0 basis: traceless elements of J3(O)
# Standard basis e_0, e_1, ..., e_26 of J3(O)
# Trace of e_i: Tr(X) = X_0 + X_1 + X_2 (the three diagonal entries)
# Traceless subspace: {X : X_0 + X_1 + X_2 = 0}

basis_27 = j3_basis()  # standard basis
# Traceless basis: 26-dimensional
# Use: e_0 - e_1, e_0 - e_2, e_3, e_4, ..., e_26
traceless_basis = []
# First two: diagonal traceless combinations
v1 = basis_27[0] - basis_27[1]  # xi1 - xi2
v2 = basis_27[0] - basis_27[2]  # xi1 - xi3
traceless_basis.append(v1)
traceless_basis.append(v2)
# The off-diagonal ones are already traceless
for i in range(3, 27):
    traceless_basis.append(basis_27[i])
print(f"  Traceless J3(O)_0 basis: {len(traceless_basis)} vectors")

# Build L_X operators: L_X(Y) = X ∘ Y - (1/3)Tr(X∘Y)·I
# In matrix form: (L_X)_{ij} = J[X_idx, j, i] - (1/3) Σ_k δ_{k∈{0,1,2}} J[X_idx, j, k] δ_i?
# Actually, let's compute it directly.

# I = basis with trace 3: I = e_0 + e_1 + e_2
I_j3 = np.zeros(27)
I_j3[0] = I_j3[1] = I_j3[2] = 1.0

def build_L_X(X, J):
    """Build the 27x27 matrix for L_X(Y) = X∘Y - (1/3)Tr(X∘Y)·I."""
    dim = 27
    L = np.zeros((dim, dim))
    for j in range(dim):
        # X ∘ e_j = Σ_i J[X, j, i] e_i (where X is expressed in basis)
        # But X is a vector, so X ∘ e_j = Σ_a X_a · (e_a ∘ e_j)_i = Σ_a,i X_a J[a,j,i]
        Xej = np.zeros(dim)
        for a in range(dim):
            if abs(X[a]) > 1e-15:
                Xej += X[a] * J[a, j, :]
        # Tr(X∘e_j) = Xej[0] + Xej[1] + Xej[2]
        tr_Xej = Xej[0] + Xej[1] + Xej[2]
        # L_X(e_j) = X∘e_j - (1/3)tr·I
        result = Xej - (tr_Xej / 3.0) * I_j3
        L[:, j] = result
    return L

p_gens = []
for X in traceless_basis:
    L = build_L_X(X, J)
    p_gens.append(L)
print(f"  L_X generators: {len(p_gens)}")

# Full e6 basis: f4 ∪ p
e6_gens = list(f4_gens) + p_gens
print(f"  Total e6 generators: {len(e6_gens)}")

# Check: are L_X traceless? (They should be since X is traceless and L is traceless projection.)
max_trace_p = max(abs(np.trace(g)) for g in p_gens)
print(f"  Max trace of L_X: {max_trace_p:.2e}")

# Verify closure
print("\n  Checking closure...")
closes = 0
total = len(e6_gens) * (len(e6_gens)-1) // 2
for i in range(len(e6_gens)):
    for j in range(i+1, len(e6_gens)):
        C = commutator_np(e6_gens[i], e6_gens[j])
        norm_C = np.linalg.norm(C)
        if norm_C < 1e-10:
            closes += 1
        elif check_in_span_np(C, e6_gens, tol=1e-4):
            closes += 1
print(f"  Closure: {closes}/{total} ({closes/total:.1%})")

# If closure works, compute Killing
if closes/total > 0.5:
    K = killing_metric_np(e6_gens)
    K_diag = np.diag(K).real
    print(f"  K_diag mean: {np.mean(K_diag):.6f} +/- {np.std(K_diag):.6f}")
    K_eigs = np.linalg.eigvalsh(K.real)
    print(f"  K eigenvalues: min={K_eigs[0]:.4f} max={K_eigs[-1]:.4f}")
    
    kappa = np.mean([np.trace(g @ g).real for g in e6_gens])
    k_ratio = abs(np.mean(K_diag)) / abs(kappa)**2
    print(f"  kappa: {kappa:.6f}")
    print(f"  |K|/kappa^2: {k_ratio:.4f}")
