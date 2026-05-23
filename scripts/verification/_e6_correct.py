"""
Correct approach: Use the 78 cubic-null-space generators (which DO close at 100%)
and compute the TRUE Killing form via structure constants, then analyze.

The Killing form must be computed correctly:
  K(T_a, T_b) = c_{acd} c_{bdc} (sum over c,d)
where c_{abc} are structure constants in the 78-dim algebra basis.

For a SIMPLE Lie algebra, K should be proportional to identity in the Cartan-Weyl basis.
In our SVD basis it won't be diagonal, but it WILL be non-degenerate.
"""
import sys, time
sys.path.insert(0, 'verification')
import numpy as np
from octonion_jordan import build_e6_generators, _orthonormalize_frobenius
from shared_utils_gpu import (
    killing_metric_np, commutator_np, check_in_span_np,
    structure_constants_general_np, orthonormalize_killing_np
)

print("=== E6 CORRECT ANALYSIS ===")
gens = build_e6_generators()
n = len(gens)

# We know closure is 100%. The issue is the Killing metric computation.
# 
# killing_metric_np computes: K[i,j] = sum_a Tr([T_i, T_a] @ [T_j, T_a])
# This equals: sum_a Tr(ad_i(T_a) @ ad_j(T_a))
# where ad_i(T_a) = [T_i, T_a] is computed in the REPRESENTATION space (27-dim).
# 
# This is actually B_V(T_i, T_j) = Tr_V(rho(T_i) rho(T_j)) but evaluated differently.
# Specifically:
# killing_metric_np gives: sum_a Tr_V([T_i, T_a][T_j, T_a])
#   = sum_{a,p,q} ([T_i, T_a])_{pq} ([T_j, T_a])_{qp}
#   = sum_{a,p,q} (T_i T_a - T_a T_i)_{pq} (T_j T_a - T_a T_j)_{qp}
#
# This is the "Killing-like" form using the representation.
# For the TRUE Killing form we need ad acting on the algebra (78-dim), not V (27-dim).
#
# TRUE Killing:
# First compute structure constants c: [T_i, T_j] = sum_k c_{ijk} T_k
# Then K_true[i,j] = sum_{c,d} c_{icd} c_{jdc}

# The B_V form (what killing_metric_np gives) is proportional to K_true
# by the Dynkin index: B_V = l_V * K_true
# But only if the representation is irreducible.

# For E6 in the 27-dim representation: this IS the fundamental irrep, 
# and l_27 (Dynkin index) = 3 (standard convention with long roots squared = 2).

# So: B_V = 3 * K_true, or K_true = B_V / 3
# And our formula h∨ = I_R * |K_diag| / κ² uses B_V directly with I_R normalization.

# Let me verify: for the raw Frobenius-orthonormal generators, what is
# B_V (killing_metric_np output)?

print("\n  Computing B_V = Tr_27([T_i,T_a][T_j,T_a])...")
t0 = time.time()
B_V = killing_metric_np(gens)
t1 = time.time()
print(f"  Computed in {t1-t0:.1f}s")

B_V_diag = np.diag(B_V).real
B_V_eigs = np.linalg.eigvalsh(B_V.real)
print(f"  B_V diag: mean={np.mean(B_V_diag):.6f}, std={np.std(B_V_diag):.6f}")
print(f"  B_V eigenvalues: min={B_V_eigs[0]:.4f} max={B_V_eigs[-1]:.4f}")
print(f"  Negative eigenvalues: {sum(e < -0.01 for e in B_V_eigs)}")
print(f"  Near-zero: {sum(abs(e) < 0.01 for e in B_V_eigs)}")
print(f"  Positive: {sum(e > 0.01 for e in B_V_eigs)}")

# Diagonalize B_V to get Killing-adapted basis
evals, evecs = np.linalg.eigh(B_V.real)
print(f"\n  Eigenvalue distribution of B_V:")
neg = [e for e in evals if e < -0.01]
zero = [e for e in evals if abs(e) < 0.01]
pos = [e for e in evals if e > 0.01]
print(f"  Negative: {len(neg)}, range [{min(neg):.4f}, {max(neg):.4f}]" if neg else "  No negative")
print(f"  Zero: {len(zero)}" if zero else "  No zero")
print(f"  Positive: {len(pos)}, range [{min(pos):.4f}, {max(pos):.4f}]" if pos else "  No positive")

# If we have both positive and negative: we have a non-compact real form.
# For e6, the split form e6(6) has signature (42, 36) on the Killing.
# The compact form e6(-78) has signature (0, 78).
# Our cubic-only constraint gives e6(6) (split), not e6(-78) (compact).

# For the framework verification, we have 3 options:
# A) Use the split form as-is (not physically relevant for gauge theory)
# B) Add the compactness constraint (antisymmetry w.r.t. an appropriate bilinear form)  
# C) Compute adjoint Killing correctly and extract h∨

# Let's try option C: compute the TRUE Killing via struct constants
# and check if it gives h∨ = 12 for ANY real form of e6.

# Actually wait. killing_metric_np computes:
# K[i,j] = sum_a Tr( [T_i, T_a] @ [T_j, T_a] )
# If representations are REAL, this is:
# = sum_a Tr(T_i T_a T_j T_a - T_i T_a^2 T_j - T_a T_i T_j T_a + T_a T_i T_a T_j)
# Hmm, this is NOT the same as B_V(T_i, T_j) = Tr(T_i T_j).
#
# Actually killing_metric_np IS computing the "trace form in the ad-representation":
# By definition, for the ad-representation on V = span{T_a}:
# The representation rho: g -> End(V) maps T_i -> ad_{T_i}
# where ad_{T_i}(T_a) = [T_i, T_a]
#
# Then: Tr_V(rho(T_i) rho(T_j)) = sum_a <[T_i, [T_j, T_a]], T_a>_trace
# where <A,B> = Tr(A B) (trace form on End(R^27)).
#
# This equals: sum_a Tr_27( [T_i, [T_j, T_a]] @ T_a )
#
# But killing_metric_np computes: sum_a Tr_27( [T_i, T_a] @ [T_j, T_a] )
# These differ by a sign from cyclicity:
# Tr([T_i, T_a] [T_j, T_a]) = Tr(T_i T_a T_j T_a) - Tr(T_a T_i T_j T_a)
# while
# Tr([T_i, [T_j, T_a]] T_a) = Tr(T_i T_j T_a^2) - Tr(T_i T_a T_j T_a)
# These ARE the same when you sum over a (by cyclicity of trace)? No...

# Let me just compute structure constants directly and get the true adjoint Killing.
print("\n  Computing structure constants [T_i, T_j] = c_{ijk} T_k ...")
t0 = time.time()
c = structure_constants_general_np(gens)
t1 = time.time()
print(f"  Computed in {t1-t0:.1f}s")

# Verify Jacobi identity (sanity check)
print("  Checking Jacobi identity...")
max_jacobi = 0.0
for i in range(min(10, n)):
    for j in range(i+1, min(10, n)):
        for k in range(j+1, min(10, n)):
            val = 0.0
            for m in range(n):
                val += c[i,j,m]*c[m,k,0] + c[j,k,m]*c[m,i,0] + c[k,i,m]*c[m,j,0]
            max_jacobi = max(max_jacobi, abs(val))
print(f"  Max Jacobi error (sample): {max_jacobi:.2e}")

# True adjoint Killing
K_adj = np.einsum('iac,jca->ij', c, c)
K_adj_diag = np.diag(K_adj).real
K_adj_eigs = np.linalg.eigvalsh(K_adj.real)
print(f"\n  Adjoint Killing K_adj[i,j] = c_iac c_jca:")
print(f"    Diag mean: {np.mean(K_adj_diag):.4f}, std: {np.std(K_adj_diag):.4f}")
print(f"    Eigenvalues: min={K_adj_eigs[0]:.4f} max={K_adj_eigs[-1]:.4f}")
print(f"    Negative: {sum(e < -0.01 for e in K_adj_eigs)}")
print(f"    Zero: {sum(abs(e) < 0.01 for e in K_adj_eigs)}")
print(f"    Positive: {sum(e > 0.01 for e in K_adj_eigs)}")

# For a SIMPLE algebra, K_adj should be non-degenerate (no zero eigenvalues).
# If it has both signs, it's the split form. If all negative, compact.

# Try computing h∨ using the relation between Killing and trace form
# For irreducible rep V of dimension d_V:
# Tr_V(T_i T_j) = I_V / dim(g) * K_adj(T_i, T_j)  ... no, that's not right.
# 
# The correct relation: Tr_V(T_i T_j) = I_V * g_{ij} where g is some 
# fixed bilinear form... the Dynkin index definition varies by convention.

# Let me just use the DIRECT approach that worked for all other algebras:
# h∨ = I_R × |K_diag| / κ²
# where K is what killing_metric_np returns and I_R = Dynkin index

# For E6 fundamental (27-dim): I_27 = 3 (with normalization θ²=2)
# This means: Tr_27(T_a T_b) = 3 * (2/ψ²) * K(T_a, T_b) where ψ = long root length

# The formula h∨ = I_R * |K| / κ² was derived assuming K is the representation-space 
# Killing (what killing_metric_np gives). Let's compute this for proper basis.

# First, diagonalize B_V and normalize to uniform eigenvalues
# The issue: B_V has mixed signs because e6(6) generators have indefinite B_V.
# That's because some generators are "compact" (B < 0) and some "non-compact" (B > 0).

# For our framework verification: the h∨ formula works when generators are
# Killing-orthonormalized. For non-compact forms, the Killing has indefinite signature.
# The h∨ is still well-defined as a property of the COMPLEXIFIED algebra.

# Let's try diagonalizing and using |eigenvalue| for normalization:
print("\n  Diagonalizing B_V and renormalizing...")
new_gens = []
for a in range(n):
    g = sum(evecs[b, a] * gens[b] for b in range(n))
    if abs(evals[a]) > 1e-8:
        g = g / np.sqrt(abs(evals[a]))
    new_gens.append(g)

# Now B_V should have eigenvalues ±1
B_new = killing_metric_np(new_gens)
B_new_diag = np.diag(B_new).real
print(f"  B_V diag after normalization: mean={np.mean(B_new_diag):.6f}, std={np.std(B_new_diag):.6f}")
print(f"  B_V diag values (first 10): {[f'{d:.4f}' for d in B_new_diag[:10]]}")

# Check closure still holds
closes = 0
total = n*(n-1)//2
for i in range(n):
    for j in range(i+1, n):
        C = commutator_np(new_gens[i], new_gens[j])
        norm_C = np.linalg.norm(C)
        if norm_C < 1e-10 or check_in_span_np(C, new_gens, tol=1e-4):
            closes += 1
print(f"  Closure: {closes}/{total} ({closes/total:.1%})")

# Compute kappa and h∨
kappa = np.mean([np.trace(g @ g).real for g in new_gens])
print(f"  kappa (mean Tr(T^2)): {kappa:.6f}")
neg_K = [d for d in B_new_diag if d < 0]
pos_K = [d for d in B_new_diag if d > 0]
print(f"  Negative K generators: {len(neg_K)}, mean={np.mean(neg_K):.4f}" if neg_K else "  No negative")
print(f"  Positive K generators: {len(pos_K)}, mean={np.mean(pos_K):.4f}" if pos_K else "  No positive")

# For the compact direction: |K|/κ² should give h∨/I_R
# Use the MEAN of |K_diag|
mean_abs_K = np.mean(np.abs(B_new_diag))
k_ratio_abs = mean_abs_K / abs(kappa)**2
print(f"  mean|K_diag|/kappa^2: {k_ratio_abs:.4f}")
print(f"  h∨ = I_R * k_ratio. For I_R=3: h∨ = {3*k_ratio_abs:.4f} (theo: 12)")

# For split form: the right formula uses the absolute value of the Killing
# which gives the SAME h∨ as the compact form.
