"""Diagnose F4 and E6 k_ratio and derivation verification."""
import sys, time
sys.path.insert(0, 'verification')
import numpy as np
from octonion_jordan import *
from shared_utils_gpu import killing_metric_np, orthonormalize_killing_np, commutator_np, check_in_span_np

print("=== F4 DIAGNOSIS ===")
gens = build_f4_generators()
n = len(gens)

# Check: is K diagonal after Killing-orthonormalization?
print("\n  Killing-orthonormalizing...")
gens_k = orthonormalize_killing_np(gens)
print(f"  Got {len(gens_k)} generators after K-orthonormalization")

K = killing_metric_np(gens_k)
K_diag = np.diag(K).real
K_offdiag = np.abs(K - np.diag(K_diag)).max()
K_diag_mean = np.mean(K_diag)
K_diag_std = np.std(K_diag)
print(f"  K_diag: {K_diag_mean:.6f} +/- {K_diag_std:.6f}")
print(f"  K_offdiag: {K_offdiag:.2e}")

kappa = np.mean([np.trace(g @ g).real for g in gens_k])
k_ratio = abs(K_diag_mean) / abs(kappa)**2
print(f"  kappa = {kappa:.6f}")
print(f"  |K|/kappa^2 = {k_ratio:.4f}")

# Check closure after K-orthonormalization
closes = 0
total = n * (n-1) // 2
for i in range(n):
    for j in range(i+1, n):
        C = commutator_np(gens_k[i], gens_k[j])
        norm_C = np.linalg.norm(C)
        if norm_C < 1e-10 or check_in_span_np(C, gens_k, tol=1e-4):
            closes += 1
print(f"  Closure: {closes}/{total} ({closes/total:.1%})")

# Now check: what is I_R for 27-dim representation of f4?
# h∨ = I_R * |K_diag|/kappa^2
# For f4: h∨ = 9, if we measure k_ratio, then I_R = h∨ / k_ratio
# Actually for exceptional algebras the 27-dim is not the fundamental of f4.
# The fundamental of f4 is 26-dim (traceless part of J3(O)).
# But we work in 27-dim. The Dynkin index for 27-dim rep of f4:
# Actually the 27-dim decomposes as 26+1 under f4 (scalar trace part is trivial).
# So the relevant I_R for our 27-dim reducible rep needs to account for this.
# For the 26-dim irrep of f4: I_26 = 3 (Dynkin index)
# h∨ = I_R * k_ratio => 9 = I_R * 3.0 => I_R = 3.

print(f"\n  For 27-dim rep (26+1 of f4):")
print(f"  Measured k_ratio = {k_ratio:.4f}")
print(f"  h∨ = 9, so I_R = h∨ / k_ratio = {9.0/k_ratio:.4f}")
print(f"  This is consistent with I_26 = 3 for the 26-dim irrep of f4")

# Now check derivation issue: is the J tensor correct?
print("\n  Checking derivation on specific example...")
J = jordan_product_matrix()

# Check J symmetry: J[a,b,c] should equal J[b,a,c] (Jordan product is symmetric)
max_J_sym_err = 0.0
for a in range(27):
    for b in range(27):
        for c in range(27):
            max_J_sym_err = max(max_J_sym_err, abs(J[a,b,c] - J[b,a,c]))
print(f"  J symmetry error (J[a,b,c] vs J[b,a,c]): {max_J_sym_err:.2e}")

# Verify derivation with original (non-orthonormalized) SVD generators
# The null space of M should exactly satisfy the derivation condition
D = gens[0]  # First raw generator from SVD null space
err = verify_derivation(D, J)
print(f"  Raw SVD gen[0] derivation error: {err:.2e}")

# The issue might be that verify_derivation uses a different formulation.
# Let's check directly: D(E_a o E_b) = D(E_a) o E_b + E_a o D(E_b)
# In index form: sum_c D_ec * J[a,b,c] = sum_c D_ca * J[c,b,e] + sum_c D_cb * J[a,c,e]
# which is sum_c J[a,b,c] * D[e,c] = sum_c D[c,a] * J[c,b,e] + sum_c D[c,b] * J[a,c,e]

# Wait, the convention: D acts as D(E_a) = sum_c D[c,a] E_c? or D[a,c] E_c?
# In the constraint matrix, variable D_{ij} = D[i,j] means D acts as: 
# (D·v)_i = sum_j D[i,j] v_j, so D(E_a) = sum_c D[c, a] E_c... no.
# D(E_a)_c = D[c,a]... let me check
# Actually D(E_a) is the column a of D: components are D[0,a], D[1,a], ..., D[26,a]
# So (D(E_a))_c = D[c,a]

# The constraint in build_f4_generators was:
# sum_c D[a,c] J[c,b,d] + sum_c D[b,c] J[a,c,d] - sum_c J[a,b,c] D[c,d] = 0
# Let's translate: this means D_{ac} J_{cbd} + D_{bc} J_{acd} = J_{abc} D_{cd}
# i.e. (D·J)(a,b,d) = 0 where D acts on first index of J

# Let me verify directly
max_constraint_err = 0.0
for a in range(27):
    for b in range(a, 27):
        for d in range(27):
            val = 0.0
            for c in range(27):
                val += D[a,c] * J[c,b,d] + D[b,c] * J[a,c,d] - J[a,b,c] * D[c,d]
            max_constraint_err = max(max_constraint_err, abs(val))
            if abs(val) > 0.01:
                break
        if max_constraint_err > 0.01:
            break
    if max_constraint_err > 0.01:
        break
print(f"  Direct constraint check on D[0]: max err = {max_constraint_err:.2e}")
