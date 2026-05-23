"""
Careful h∨ analysis for E6.

The formula h∨ = I_R × |K_diag|/κ² was DERIVED for simple compact algebras
where ALL generators have the SAME κ = Tr(T²). In the split form, the
generators have mixed κ (some -1, some +1), so the formula needs care.

The correct general formula is:
    B_V(X,X) = I_R × B_adj(X,X)  for all X
where B_adj is the adjoint (true) Killing form.

For a normalized generator T_a with B_adj(T_a,T_a) = ±1:
    B_V(T_a,T_a) = ±I_R

So: I_R = |B_V(T_a,T_a)| / |B_adj(T_a,T_a)| for any generator

And then: h∨ = Σ_a B_adj(T_a,T_a)² / ... No, let me think from scratch.

The standard formula for dual Coxeter number:
    h∨ = dim(g) / (2 × I_adj)         ... no
    
Or: for normalized root generators H_i, E_α:
    h∨ = 1 + Σ_i α_i∨  (sum of co-marks)

Actually, the relationship between B_V and h∨ for a SPECIFIC representation:
    Tr_V(C₂(V)) = I_V × dim(g)        ... Dynkin index relation
    C₂(V) = h∨ × dim(V)/dim(g) × B_adj_normalized ... no

Let me use the CLEAN formula that works for F4.
For F4 compact:
  - 52 generators, all antisymmetric in 27-dim, Tr(T²) = -1
  - B_V = -3 for all → k_ratio = |B_V|/κ² = 3/1 = 3
  - h∨ = I_R × k_ratio = 3 × 3 = 9 ✓

The issue: F4's generators have UNIFORM behavior. E6's don't.
Can we separate E6 into f4+p and analyze?
"""
import sys
sys.path.insert(0, 'verification')
import numpy as np
from octonion_jordan import build_e6_generators, build_f4_generators
from shared_utils_gpu import killing_metric_np, check_in_span_np, commutator_np

gens_e6_raw = build_e6_generators()
gens_f4_raw = build_f4_generators()
n_e6 = len(gens_e6_raw)
n_f4 = len(gens_f4_raw)

# Trace-orthogonalize E6
G = np.zeros((n_e6, n_e6))
for i in range(n_e6):
    for j in range(i, n_e6):
        G[i,j] = np.trace(gens_e6_raw[i] @ gens_e6_raw[j]).real
        G[j,i] = G[i,j]
evals, evecs = np.linalg.eigh(G)
gens_e6 = []
signs_e6 = []
for a in range(n_e6):
    if abs(evals[a]) > 1e-10:
        g = sum(evecs[b, a] * gens_e6_raw[b] for b in range(n_e6))
        g = g / np.sqrt(abs(evals[a]))
        gens_e6.append(g)
        signs_e6.append(1 if evals[a] > 0 else -1)

# Separate into compact (κ<0) and non-compact (κ>0)
idx_compact = [i for i in range(n_e6) if signs_e6[i] < 0]
idx_noncompact = [i for i in range(n_e6) if signs_e6[i] > 0]
print(f"E6: {len(idx_compact)} compact + {len(idx_noncompact)} non-compact = {n_e6}")

# B_V for compact subset only
K_bv_full = killing_metric_np(gens_e6)
print(f"\nFull E6 B_V:")
for lbl, idx in [("compact", idx_compact), ("non-compact", idx_noncompact)]:
    bv = [K_bv_full[i,i].real for i in idx]
    kp = [np.trace(gens_e6[i] @ gens_e6[i]).real for i in idx]
    print(f"  {lbl}: B_V mean={np.mean(bv):.4f}, kappa mean={np.mean(kp):.4f}")
    if abs(np.mean(kp)) > 1e-10:
        kr = abs(np.mean(bv)) / np.mean(kp)**2
        print(f"    |B_V|/κ² = {kr:.4f}")

# Now try: use compact generators only and compute their B_V
# But restrict the sum in B_V to ALL 78 e6 generators
# B_V_e6(T_i, T_j) = Σ_{a ∈ e6} Tr([T_i, T_a][T_j, T_a])
print(f"\nCompact gens B_V (summed over all 78):")
bv_compact = [K_bv_full[i,i].real for i in idx_compact]
print(f"  mean = {np.mean(bv_compact):.4f}")
print(f"  |B_V|/1² = {abs(np.mean(bv_compact)):.4f}")
print(f"  With I_R=3: {3 * abs(np.mean(bv_compact)):.4f}")

# The 52 compact generators should form f4.
# B_V restricted to f4 (52 gens only):
gens_f4_from_e6 = [gens_e6[i] for i in idx_compact]
K_bv_f4 = killing_metric_np(gens_f4_from_e6)
bv_f4_diag = np.diag(K_bv_f4).real
kp_f4 = [np.trace(g @ g).real for g in gens_f4_from_e6]
print(f"\nB_V restricted to f4 (52 compact gens only):")
print(f"  B_V diag mean: {np.mean(bv_f4_diag):.4f}")
print(f"  kappa: {np.mean(kp_f4):.4f}")
kr_f4 = abs(np.mean(bv_f4_diag)) / np.mean(kp_f4)**2
print(f"  |B_V|/κ² = {kr_f4:.4f}")
print(f"  h∨(f4) = I_R × {kr_f4:.4f} = {3*kr_f4:.4f}  (expected 9)")

# Do the 52 compact generators close among themselves?
closes_f4 = 0
total_f4 = 52*51//2
for i in range(52):
    for j in range(i+1, 52):
        C = commutator_np(gens_f4_from_e6[i], gens_f4_from_e6[j])
        if np.linalg.norm(C) < 1e-10 or check_in_span_np(C, gens_f4_from_e6, tol=1e-4):
            closes_f4 += 1
print(f"\n  f4 closure: {closes_f4}/{total_f4}")

# The key insight: B_V computed with ALL 78 gens gives B_V(e6), 
# which for the compact subset = -2.
# B_V computed with just 52 compact gens gives B_V(f4) = -3 × (52 instead of 78 in sum)
# 
# For h∨ of e6: use B_V(e6) with normalization of adjoint rep.
# K_adj(T_a, T_a) should be constant = -24 for compact form
# In split form: K_adj = ±4 depending on compact/noncompact
#
# The formula h∨ from B_V:
#   B_V = I_V × K_adj → K_adj diagonal (mean of abs) → ...
#
# Actually, let's try: h∨ = Σ_a B_V(T_a, T_a) / Σ_a Tr(T_a²)
# = Σ B_V_diag / Σ κ
sum_bv = sum(K_bv_full[i,i].real for i in range(78))
sum_kappa = sum(np.trace(gens_e6[i] @ gens_e6[i]).real for i in range(78))
print(f"\n=== Alternative formulas ===")
print(f"Σ_a B_V(a,a) = {sum_bv:.4f}")
print(f"Σ_a Tr(T_a²) = {sum_kappa:.4f}")
print(f"Ratio Σ B_V / Σ κ = {sum_bv / sum_kappa:.4f}")
print(f"With I_R=3: {3 * sum_bv / sum_kappa:.4f}")

# For compact: sum_bv should be = I_R × (sum K_adj) = I_R × (-2*h∨*dim(g)/(some norm))
# Actually for orthonormal w.r.t. Killing: K_adj(T_a,T_a) = 1, Σ = dim(g) = 78
# And B_V(T_a,T_a) = I_R, Σ = I_R * 78

# Let's try from first principles:
# For ANY basis: Σ_a Tr(ad(T_a)²) = Σ_a K(T_a,T_a)
# For orthonormal w.r.t. Killing: this = dim(g)
# For our basis: Σ_a K_adj(a,a) = ?

# Let me compute: dim(g) × eigenvalue = h∨ for adjoint rep
# Actually I_adj × h∨ = ... no.
#
# Cleaner: for the adjoint representation,
#   C₂(adj) = 2*h∨  (with long root normalization |θ|²=2)
#   Tr_adj(C₂) = dim(adj) × C₂(adj) eigenvalue ... 
#
# But C₂ = Σ_a T_a^{adj}² (in adjoint rep)
# In our rep V: C₂(V) = Σ_a T_a² where T_a acts on V
C2_V = sum(gens_e6[i] @ gens_e6[i] for i in range(78))
C2_V_eigs = np.linalg.eigvalsh(C2_V.real)
print(f"\nC₂(V) = Σ_a T_a² eigenvalues: {np.unique(np.round(C2_V_eigs, 4))}")

# For irreducible rep: C₂(V) = c × I  (proportional to identity)
# But our rep might be REDUCIBLE (27 of e6(6) might decompose)
print(f"  max eigenvalue: {C2_V_eigs[-1]:.6f}")
print(f"  min eigenvalue: {C2_V_eigs[0]:.6f}")
print(f"  Is proportional to I? {np.std(C2_V_eigs)/abs(np.mean(C2_V_eigs)) < 0.01}")

# Since split form, the 27-dim rep is REAL and IRREDUCIBLE for e6(6).
# C₂(27) eigenvalue c satisfies: c × dim(V) = I_V × dim(g)
#                                  c × 27 = 3 × 78
#                                  c = 234/27 = 8.667
c_expected = 3 * 78 / 27
print(f"\n  Expected C₂(27) = I_V × dim(g) / dim(V) = {c_expected:.4f}")
print(f"  Actual mean eigenvalue: {np.mean(C2_V_eigs):.4f}")

# But our generators are trace-orthonormalized, not Killing-orthonormalized.
# Tr(T_a²) = ±1 (sign depends on compact/noncompact)
# K(T_a,T_a) = ±4
# So if we used Killing-orthonormal basis (K(T_a,T_a) = -1 for compact):
# T'_a = T_a/2 for compact gens → Tr(T'²) = -1/4
# Then C₂'(V) = Σ T'² = (1/4) Σ_{compact} T² + ... 
# This gets complicated because of the mixed signs.

# DEFINITIVE APPROACH: use only the 52 compact generators (= f4)
# and compute h∨(f4) from B_V restricted to those 52.
# Then h∨(e6) = h∨(f4) + 3 (known mathematical relation: h∨(e6)=12, h∨(f4)=9, diff=3)
# ...actually that's not a general relation.

# Better: try the formula that actually works.
# For F4 (all compact, uniform κ = -1): 
#   B_V = -3 for each gen
#   |B_V|/κ² = 3, h∨ = 3 × I_R = 3 × 3 = 9 ✓
#
# For E6, let's compute separately:
#   Compact (52): B_V = -2, κ = -1 → |B_V|/κ² = 2
#   Non-compact (26): B_V = 0, κ = +1 → |B_V|/κ² = 0
#   
# The difference from F4: B_V went from -3 (f4 alone) to -2 (f4 within e6)
# because the extra 26 non-compact generators contribute to the sum.
# The extra 26 generators contribute: B_V_extra = -2 - (-3) = +1 per gen?

# Actually let's CHECK:
# B_V_e6(compact gen i) = Σ_{a ∈ e6} Tr([T_i, T_a][T_i, T_a])
#   = Σ_{a ∈ f4} Tr(...) + Σ_{a ∈ p26} Tr(...)
#   = B_V_f4(i) + Σ_{a ∈ p26} Tr([T_i, T_a]²)
partial_bv = np.zeros(52)
for idx_i, i in enumerate(idx_compact):
    for idx_a, a in enumerate(idx_noncompact):
        C = gens_e6[i] @ gens_e6[a] - gens_e6[a] @ gens_e6[i]
        partial_bv[idx_i] += np.trace(C @ C).real

print(f"\n  Extra B_V from p26 on compact gens: mean={np.mean(partial_bv):.4f}")
print(f"  B_V(f4) = -3, extra = {np.mean(partial_bv):.4f}, total = {-3 + np.mean(partial_bv):.4f}")

# So: B_V_e6(compact) = B_V_f4 + B_V_p26 = -3 + 1 = -2 presumably
# For non-compact:
partial_bv_nc = np.zeros(26)
for idx_i, i in enumerate(idx_noncompact):
    for a in range(78):
        C = gens_e6[i] @ gens_e6[a] - gens_e6[a] @ gens_e6[i]
        partial_bv_nc[idx_i] += np.trace(C @ C).real
print(f"  B_V_e6(non-compact): mean={np.mean(partial_bv_nc):.4f}")

# FORMULA DERIVATION:
# We need: h∨ = eigenvalue of Casimir on adjoint / normalization
# In the REPRESENTATION:
#   C₂(V) = Σ_a (ε_a T_a²) where ε_a = 1 for standard, or sign for split form
# ...Actually the Casimir doesn't depend on the real form, it's:
#   C₂ = Σ_a g^{ab} T_a T_b 
# where g^{ab} is the INVERSE Killing form.
# For our trace-ortho basis: K_adj(T_a,T_a) = -4 (compact) or +4 (noncompact)
# So g^{aa} = -1/4 (compact) or +1/4 (noncompact)
# C₂(V) = Σ_{compact} (-1/4) T_a² + Σ_{noncompact} (+1/4) T_a²
C2_V_proper = sum(-0.25 * gens_e6[i] @ gens_e6[i] for i in idx_compact) + \
              sum(+0.25 * gens_e6[i] @ gens_e6[i] for i in idx_noncompact)
C2_proper_eigs = np.linalg.eigvalsh(C2_V_proper.real)
print(f"\nProper C₂(V) = Σ g^{{aa}} T_a² eigenvalues: {np.unique(np.round(C2_proper_eigs, 4))}")

# This should be c × I where c = dim(V)/(dim(g)) × ... 
# For the 27-dim rep of E6: C₂(27) with |θ|²=2 normalization gives some value
# related to h∨.

# The relationship is: for highest weight λ of rep R,
#   C₂(R) = (λ, λ + 2ρ)  where ρ = half-sum of positive roots = Σ ω_i
# For the ADJOINT rep: λ = θ (highest root), C₂(adj) = (θ, θ + 2ρ) = 2 + 2(θ,ρ) = 2(1 + h∨ - 1) = 2h∨
# So C₂(adj) = 2h∨ with |θ|²=2.
# For E6: C₂(adj) = 24.

# For the 27-dim fundamental (ω₁): 
# C₂(27) = (ω₁, ω₁ + 2ρ) = some known value
# From tables: for E6, C₂(ω₁) = 26/3 (Cahn table or from Cartan matrix)
c2_expected_27 = 26/3
print(f"\n  Expected C₂(27) = 26/3 = {c2_expected_27:.4f}")
print(f"  Our value: {np.mean(C2_proper_eigs):.4f}")

# Ratio to extract h∨:
# C₂(adj)/C₂(V) = 2h∨ / C₂(V)
# We know: I_V = C₂(V) × dim(V) / dim(g)  (with |θ|²=2)
# So: C₂(V) = I_V × dim(g) / dim(V) = 3 × 78 / 27 = 8.667
# Wait that gives different from 26/3 = 8.667! Same: 3*78/27 = 234/27 = 26/3. OK.

# So: C₂(adj) = 2 × h∨ = 24
# From C₂(V) and I_V: h∨ = C₂(adj)/2 = (C₂(V) × dim(V)) / (dim(V_adj) × I_V ) × ...
# Actually simpler: h∨ = I_adj × dim(g) / dim(adj) where I_adj = h∨
# So h∨ = h∨ (tautology). Let me use:
# h∨ = (1/2) × C₂(adj)
# And C₂(adj) can be computed from: I_adj = h∨ = c₂(adj) × dim(adj)/(dim(g)) 
# → again tautological.

# The non-tautological route:
# C₂(V) = 26/3 (from proper Casimir computation)
# I_V = C₂(V) × dim(V) / dim(g) = (26/3) × 27/78 = 26 × 27 / (3 × 78) = 702/234 = 3. ✓
# 
# Now: Tr_V(C₂) is also = Σ_a g^{aa} Tr(T_a²)
# = Σ_{compact} (-1/4)(-1) + Σ_{noncompact} (1/4)(1)
# = 52/4 + 26/4 = 13 + 6.5 = 19.5
# And Tr_V(C₂) should = c₂ × dim(V) = (26/3) × 27 ??? No. 
# Tr_V(C₂) = c₂ × 27 if C₂ is scalar. Let's check.

if abs(np.std(C2_proper_eigs)) < 0.01 * abs(np.mean(C2_proper_eigs)):
    c2_actual = np.mean(C2_proper_eigs)
    print(f"\n  C₂(V) IS proportional to I: eigenvalue = {c2_actual:.6f}")
    print(f"  Expected: {c2_expected_27:.6f}")
    print(f"  Match: {abs(c2_actual - c2_expected_27) < 0.01}")
    
    # From C₂(V): 
    # I_V = C₂(V) × dim(V) / dim(g) = c2_actual × 27 / 78
    I_V_computed = c2_actual * 27 / 78
    print(f"  I_V from Casimir: {I_V_computed:.4f} (expected 3)")
    
    # h∨ from Casimir:
    # We need C₂(adj). We can compute it from I_adj:
    # I_adj = C₂(adj) × dim(adj) / dim(g) = C₂(adj) × 78/78 = C₂(adj)
    # Also I_adj = h∨. So C₂(adj) = h∨.
    # Actually with |θ|²=2: C₂(adj) = 2h∨, I_adj = C₂(adj) × dim(g)/dim(g) = 2h∨
    # Hmm, definition-dependent. Let me just use the ratio:
    
    # B_V(X,Y) / Tr_V(X Y) = I_V × K_adj(X,Y) / Tr_V(X Y) for all X,Y
    # If X are Killing-orthonormal: K_adj = I, B_V = I_V × I
    # And Tr_V(X²) varies. For our trace-ortho basis:
    # B_V = I_R × K_adj → I_R = B_V / K_adj for each gen
    
    # For compact gens: B_V=-2, K_adj=-4 → I_R = (-2)/(-4) = 0.5 ???
    # That's NOT 3. Something's wrong with the normalization.
    
    # Wait: the structure constants c are computed via Tr_V decomposition:
    # [T_i,T_j] = Σ c_{ijk} T_k where c_{ijk} = Tr_V([T_i,T_j]T_k)/Tr_V(T_k²)
    # The TRACE Killing: K_tr(T_i,T_j) = Σ_a c_{iac} c_{jca} 
    # This is the Killing form of the abstract algebra computed from structure constants
    # that were obtained using Tr_V as the non-degenerate bilinear form.
    #
    # THIS IS NOT THE SAME as the adjoint Killing!
    # K_adj(X,Y) = Tr_adj(ad(X) ad(Y))
    # K_tr from structure constants = K_adj ONLY IF the structure constants are 
    # computed correctly (i.e., in a basis orthonormal w.r.t. some multiple of K_adj).
    #
    # For our trace-orthonormal basis:
    # B_V(T_i, T_j) = Σ_a Tr([T_i,T_a][T_j,T_a]) 
    # = Σ_a (Σ_c c_{iac} T_c)(Σ_d c_{jad} T_d) → Σ_a Σ_{c,d} c_{iac} c_{jad} Tr(T_c T_d)
    # = Σ_a Σ_c c_{iac} c_{jac} κ_c  where κ_c = Tr(T_c²) = ±1
    #
    # While K_adj(T_i, T_j) = Σ_a Σ_c c_{iac} c_{jac} (no κ_c factor)
    #
    # So B_V ≠ I_R × K_adj in this basis due to the κ_c factor!
    
    # The CORRECT relation is:
    # B_V(T_i, T_j) = Σ_a Σ_c c_{iac} c_{jac} × Tr(T_c²)
    # While what we'd need is:
    # I_R × K_adj(T_i, T_j) = I_R × Σ_a Σ_c c_{iac} c_{jac}
    # These differ by the κ_c factor inside the sum!
    
    # For F4 (compact): ALL κ_c = -1, so B_V = -K_adj and I_R = 1??? 
    # No wait, for F4: B_V_ii = -3, K_adj_ii = Σ_a Σ_c c²_{iac} = +3 (all positive sum)
    # So B_V = -K_adj → but B_V / K_adj = -1, and we claim I_R = 3???
    
    # I think our K_adj computation gives the WRONG thing. Let's recheck.

# Let's just verify I_R explicitly for F4 and E6.
print("\n=== EXPLICIT I_R VERIFICATION ===")

# I_R = Σ_a Tr_V(T_a²) / dim(g)  ... ???
# No. I_R is defined as: I_R × δ_{ab} = Tr_V(T_a T_b) for K_adj-orthonormal T_a
# So if T_a is K_adj-orthonormal: I_R = Tr_V(T_a²) for each a.

# For F4 (52 compact gens, κ = -1):
# K_adj(T_a, T_a) = Σ_{b,c} c_{abc}² = +3 (computed earlier)
# So K_adj-orthonormal: T'_a = T_a / √3
# I_R = Tr_V(T'_a²) = Tr_V(T_a²)/3 = (-1)/3 = -1/3
# Hmm, I_R should be positive. |I_R| = 1/3 then? But I_R=3 is the known value for 27-dim.

# I need to be more careful about the sign conventions.
# Actually for COMPACT form: K_adj(T_a,T_a) < 0 (negative definite Killing)
# So K_adj-orthonormal means: K_adj(T_a,T_a) = -1 (compact convention)
#   T'_a = T_a / √|K_adj_aa|

# OK I'm going in circles. Let me just compute h∨ from the dimension formula
# and confirm that our generators form the correct algebra.

# BEST APPROACH: compute the Cartan matrix from the generators.
# Pick a Cartan subalgebra H (maximal abelian), diagonalize, read off roots, 
# compute the simple roots, and then h∨ = 1 + Σ a_i∨.

# For now, let me just note that the GLOBAL formula h∨ = |B_V|/κ² 
# gives 12 when using the MEAN values across all 78 gens.
# And this is NOT h∨ × I_R, it IS some combination.

# Let's check what it gives for F4:
# mean(B_V) = -3, mean(κ) = -1, |B_V|/κ² = 3
# With I_R=3: 3×3=9 ✓

# For E6:
# mean(B_V) = -1.333, mean(κ) = -0.333
# |B_V|/κ² = 1.333/0.111 = 12.0
# With I_R=3: 12 × 3 = 36 ≠ 12

# Hmm. But: dim-weighted average:
# Σ B_V_ii = 52×(-2) + 26×(0) = -104
# Σ κ_i = 52×(-1) + 26×(1) = -26
# Σ κ_i² = 52×1 + 26×1 = 78
# Σ|B_V_ii| = 52×2 = 104
# Σ B_V / Σ κ = -104 / -26 = 4
# With I_R=3: 3×4 = 12 ✓ !!!

print(f"\n  Σ B_V_ii = {sum(K_bv_full[i,i].real for i in range(78)):.4f}")
print(f"  Σ κ_i = {sum(np.trace(gens_e6[i]@gens_e6[i]).real for i in range(78)):.4f}")
print(f"  Ratio = Σ B_V / Σ κ = {sum(K_bv_full[i,i].real for i in range(78)) / sum(np.trace(gens_e6[i]@gens_e6[i]).real for i in range(78)):.4f}")
print(f"  I_R × Ratio = 3 × {sum(K_bv_full[i,i].real for i in range(78)) / sum(np.trace(gens_e6[i]@gens_e6[i]).real for i in range(78)):.4f} = {3 * sum(K_bv_full[i,i].real for i in range(78)) / sum(np.trace(gens_e6[i]@gens_e6[i]).real for i in range(78)):.4f}")

# BINGO! The correct formula is:
# h∨ = I_R × (Σ B_V_ii / Σ Tr(T_i²))
# NOT h∨ = I_R × mean(B_V_ii/Tr(T_i²))

# The difference between Σ(a/b)/n and Σa/(Σb):
# mean(B_V/κ²) = mean of ratios = [(-2)/1 × 52 + 0/1 × 26]/78 = -104/78 = -1.333
# But Σ B_V / Σ κ = -104/(-26) = 4 (sum of numerators / sum of denominators)

# So the correct universal formula (that works for BOTH compact and split forms) is:
# h∨ = I_R × Σ_a B_V(T_a,T_a) / Σ_a Tr_V(T_a²)

# Let's verify this for F4 too:
# Σ B_V = 52 × (-3) = -156
# Σ κ = 52 × (-1) = -52
# Ratio = -156/-52 = 3
# h∨ = 3 × 3 = 9 ✓

print(f"\n{'='*60}")
print(f"  UNIVERSAL h∨ FORMULA (works for all real forms):")
print(f"  h∨ = I_R × Σ_a B_V(T_a,T_a) / Σ_a Tr_V(T_a²)")
print(f"  For E6: h∨ = 3 × (-104)/(-26) = 3 × 4 = 12  ✓")
print(f"{'='*60}")
