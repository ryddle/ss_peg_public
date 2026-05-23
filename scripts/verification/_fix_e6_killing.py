"""Fix E6 by Killing-diagonalizing within the null space."""
import sys, time
sys.path.insert(0, 'verification')
import numpy as np
from octonion_jordan import build_e6_generators, build_f4_generators, cubic_norm_tensor, verify_cubic_preservation
from shared_utils_gpu import killing_metric_np, commutator_np, check_in_span_np

print("=== E6 WITH KILLING-DIAGONALIZED BASIS ===")

# Build raw generators (Frobenius-orthonormal)
gens_raw = build_e6_generators()
n = len(gens_raw)

# Compute Killing metric in this basis
print("\n  Computing Killing metric in Frobenius basis...")
t0 = time.time()
K_raw = killing_metric_np(gens_raw)
t1 = time.time()
print(f"  Killing computed in {t1-t0:.1f}s")
print(f"  K eigenvalues: {np.sort(np.linalg.eigvalsh(K_raw.real))[:5]}...{np.sort(np.linalg.eigvalsh(K_raw.real))[-5:]}")

# Diagonalize K within the null space
# K_raw = V D V^T where D is diagonal
evals, evecs = np.linalg.eigh(K_raw.real)
print(f"  Eigenvalues: min={evals[0]:.4f}, max={evals[-1]:.4f}")
print(f"  Positive eigenvalues: {sum(e > 0.01 for e in evals)}")
print(f"  Near-zero eigenvalues: {sum(abs(e) < 0.01 for e in evals)}")
print(f"  Negative eigenvalues: {sum(e < -0.01 for e in evals)}")

# The issue: E6 (compact real form e6(-78)) should have NEGATIVE DEFINITE Killing.
# If we get positive eigenvalues, either:
# a) We have the wrong real form (split e6 instead of compact)
# b) The Killing computation is wrong
# c) The generators are from the complexification

# Let's check: for the COMPACT e6, all generators should satisfy K(X,X) < 0.
# This means [X, [X, *]] has all negative eigenvalues when X is non-zero.
# In the 27-dim rep, the compact form acts on REAL 27-space.

# Actually, the cubic norm preservation condition ALONE gives the
# COMPLEXIFIED e6 null space (real dimension 78 for real matrices).
# The compact real form e6(-78) is the subset that additionally
# preserves a positive-definite hermitian form on the 27-dim space.
# But since we work over R, the constraint is: generators must preserve
# the trace bilinear form Tr(X o Y).

# Let's try a different approach: compute the Killing form properly
# using the ad-representation within this 78-dim algebra.

# The Killing form K(X,Y) = Tr(ad_X ad_Y) where ad acts on the 78-dim algebra itself.
# Our killing_metric_np computes Tr(ad_X ad_Y) but acting on the 27-dim representation,
# which gives the weight-space trace B_27(X,Y) = Tr_27(rho(X) rho(Y)).
# These are proportional only when the Killing is well-defined.

# Actually killing_metric_np computes:
# K(X,Y) = sum_a,b Tr([X,T_a][Y,T_b]) * delta_ab component = sum_a Tr(ad_X(T_a) @ ad_Y(T_a))
# Wait no, let me re-read the code...

# killing_metric_np: K[i,j] = sum_a Tr([T_i, T_a] @ [T_j, T_a])
# This is the Killing form of the ADJOINT representation = true Killing.
# But actually this is K[i,j] = Tr(ad_i ad_j) only if the ad acts in the 
# representation space, not the algebra itself.

# Hmm. Let me think more carefully.
# For a Lie algebra g with representation rho: g -> gl(V),
# the Killing form is K(X,Y) = Tr_g(ad_X ad_Y)
# where ad_X(Z) = [X,Z] acts on g.
# 
# What killing_metric_np computes is:
# B_rho(X,Y) = Tr_V(rho(X) rho(Y)) -- trace in the representation V.
# 
# These are DIFFERENT. B_rho = l(rho) * K where l(rho) is the Dynkin index.
# But wait, killing_metric_np computes:
# K[i,j] = sum_a Tr([T_i, T_a] @ [T_j, T_a])
# This uses T_a as basis of g and computes [T_i, T_a] in the rep space.
# 
# Actually: [T_i, T_a] = sum_b c_{iab} T_b (expansion in algebra basis)
# So Tr([T_i, T_a] @ [T_j, T_a]) = sum_{b,c} c_{iab} c_{jac} Tr(T_b T_c)
# This is NOT the Killing form unless Tr(T_b T_c) = delta_{bc}.
# It's K[i,j] = sum_{a,b,c} c_{iab} c_{jac} G_{bc} where G = Tr(T_b T_c).

# For Frobenius-orthonormal generators (our case), G_{bc} ≈ delta_{bc} * kappa,
# but not exactly since mixed antisym+sym generators.

# The TRUE Killing form is:
# K_true[i,j] = sum_{a,b} c_{iab} c_{jba} (trace of ad_i @ ad_j as matrices on g)

# Let me compute both and compare.
print("\n  Computing TRUE Killing form from structure constants...")
# First get structure constants
from shared_utils_gpu import structure_constants_general_np
c = structure_constants_general_np(gens_raw)

# True Killing: K[i,j] = sum_{a,b} c[i,a,b] * c[j,b,a]
# (This is Tr on g-space, using the structure constants)
K_true = np.einsum('iab,jba->ij', c, c)
evals_true = np.linalg.eigvalsh(K_true.real)
print(f"  True Killing eigenvalues: min={evals_true[0]:.4f} max={evals_true[-1]:.4f}")
print(f"  Negative definite: {all(evals_true < -1e-6)}")
print(f"  K_true diag mean: {np.mean(np.diag(K_true).real):.4f}")
print(f"  K_true diag std:  {np.std(np.diag(K_true).real):.4f}")

# Diagonalize TRUE Killing
evals_t, evecs_t = np.linalg.eigh(K_true.real)

# If negative definite, use |K| gens
if all(evals_t < -1e-6):
    print("  Good: compact form e6(-78) confirmed!")
    # New basis: T'_a = sum_b V_{ba} T_b, normalized so K(T'_a, T'_a) = -1
    new_gens = []
    for a in range(n):
        g = sum(evecs_t[b, a] * gens_raw[b] for b in range(n))
        # Normalize: K(g,g) = evals_t[a], want K = -1
        g = g / np.sqrt(abs(evals_t[a]))
        new_gens.append(g)
    
    # Verify K is now -1 * I
    K_new = killing_metric_np(new_gens)
    # Wait, this is NOT the true Killing. Need to recompute from struct consts.
    c_new = structure_constants_general_np(new_gens)
    K_new_true = np.einsum('iab,jba->ij', c_new, c_new)
    K_new_diag = np.diag(K_new_true).real
    K_new_offdiag = np.abs(K_new_true - np.diag(K_new_diag)).max()
    print(f"\n  After Killing-normalization:")
    print(f"    K_true diag mean: {np.mean(K_new_diag):.6f} +/- {np.std(K_new_diag):.6f}")
    print(f"    K_true offdiag max: {K_new_offdiag:.2e}")
    
    # Now compute kappa and k_ratio
    kappa = np.mean([np.trace(g @ g).real for g in new_gens])
    k_ratio = abs(np.mean(K_new_diag)) / abs(kappa)**2
    print(f"    kappa: {kappa:.6f}")
    print(f"    |K_true|/kappa^2: {k_ratio:.4f}")
    print(f"    h∨ = 12, I_R = 12/{k_ratio:.4f} = {12.0/k_ratio:.4f}")
    
    # Closure
    closes = 0
    total = n*(n-1)//2
    for i in range(n):
        for j in range(i+1, n):
            C = commutator_np(new_gens[i], new_gens[j])
            norm_C = np.linalg.norm(C)
            if norm_C < 1e-10 or check_in_span_np(C, new_gens, tol=1e-4):
                closes += 1
    print(f"    Closure: {closes}/{total} ({closes/total:.1%})")
else:
    print("  WARNING: NOT negative definite - may not be compact real form")
    # Show distribution
    print(f"  Eigenvalue distribution:")
    for e in sorted(evals_t):
        sign = "+" if e > 0 else "-"
        print(f"    {sign} {abs(e):.4f}")
