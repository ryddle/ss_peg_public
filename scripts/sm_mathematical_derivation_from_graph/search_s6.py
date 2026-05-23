#!/usr/bin/env python3
"""
search_s6.py — Systematic search for the selection principle S6
================================================================
Among the ~3056 matrices N with Smith(N)=(1,2,4), entries ±2^v (free)
and fixed anti-diagonal, find what distinguishes the physical N.

KEY DISCOVERY TO TEST FIRST:
  sign(n_free) = (-1)^(1+Nc) : leptons negative, quarks positive.
  This reduces 18 bits → 12 bits (only 6 valuations v₂∈{0..3}).
"""
from fractions import Fraction as F
import numpy as np
from itertools import product as iprod
from collections import Counter

# ═══════════════════════════════════════════════════════════════════
# DATA
# ═══════════════════════════════════════════════════════════════════
sectors = ['lepton', 'up', 'down']
coords  = ['p', 'q', 'r']
flat = {(0,2), (1,1), (2,0)}  # (lepton/r, up/q, down/p)
free_pos = [(0,0), (0,1), (1,0), (1,2), (2,1), (2,2)]
flat_vals = {(0,2): 6, (1,1): -2, (2,0): 4}

# Physical N
N_phys = np.array([[-1, -8,  6],
                    [ 4, -2,  8],
                    [ 4,  2,  4]])

# Curvatures (from generating function)
a_mat = np.array([[1, F(1,2), F(-3,2)],
                   [F(7,4), F(1,2), F(-3,2)],
                   [-1, 0, 1]], dtype=object)

# ═══════════════════════════════════════════════════════════════════
# PART 0: CONFIRM THE SIGN RULE
# ═══════════════════════════════════════════════════════════════════
print("=" * 70)
print("  PART 0: THE SIGN RULE")
print("=" * 70)

print("\n  Physical free entries:")
print(f"  {'pos':>12s}  {'sector':>8s}  Nc    n   sign")
for (i, j) in free_pos:
    s = sectors[i]
    Nc = 0 if i == 0 else 1
    n = N_phys[i, j]
    sign = "−" if n < 0 else "+"
    expected = "−" if Nc == 0 else "+"
    match = "✓" if sign == expected else "✗"
    print(f"  ({i},{j})={s:>6s}/{coords[j]}   Nc={Nc}  {n:>+3d}    {sign}  "
          f"(expect {expected} from Nc)  {match}")

print(f"\n  RULE: sign(n_free) = (−1)^(1+Nc)")
print(f"  Leptons (Nc=0): NEGATIVE.  Quarks (Nc=1): POSITIVE.")

# ═══════════════════════════════════════════════════════════════════
# PART 1: ENUMERATE WITH SIGN RULE
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 1: ENUMERATE CANDIDATES")
print(f"{'='*70}")

def smith_diag(M):
    """Smith normal form diagonal."""
    M = [list(row) for row in M]
    n = len(M)
    for col in range(n):
        for _ in range(200):
            best = None
            for i in range(col, n):
                for j in range(col, n):
                    if M[i][j] != 0 and (best is None or abs(M[i][j]) < abs(best[0])):
                        best = (M[i][j], i, j)
            if best is None: break
            _, pi, pj = best
            if pi != col: M[pi], M[col] = M[col], M[pi]
            if pj != col:
                for i in range(n): M[i][pj], M[i][col] = M[i][col], M[i][pj]
            changed = False
            for i in range(col+1, n):
                if M[i][col] != 0:
                    q = M[i][col] // M[col][col]
                    for j in range(n): M[i][j] -= q * M[col][j]
                    if M[i][col] != 0: changed = True
            for j in range(col+1, n):
                if M[col][j] != 0:
                    q = M[col][j] // M[col][col]
                    for i in range(n): M[i][j] -= q * M[i][col]
                    if M[col][j] != 0: changed = True
            if not changed:
                ok = True
                for i in range(col+1, n):
                    for j in range(col+1, n):
                        if M[i][j] % M[col][col] != 0:
                            for jj in range(n): M[col][jj] += M[i][jj]
                            ok = False; break
                    if not ok: break
                if ok: break
        if col < n and M[col][col] < 0:
            for j in range(n): M[col][j] = -M[col][j]
    return tuple(M[i][i] for i in range(n))


# Sign pattern: lepton entries negative, quark entries positive
sign_for_row = {0: -1, 1: +1, 2: +1}  # row 0=lepton, 1=up, 2=down

# Generate: with signs FIXED, only v₂ ∈ {0,1,2,3} per free entry
candidates_sign_fixed = []
candidates_all = []  # all Smith(1,2,4) regardless of sign

# First: ALL 3056 without sign rule
for combo in iprod(iprod([-1, 1], [0, 1, 2, 3]), repeat=6):
    mat = np.zeros((3, 3), dtype=int)
    for pos, val in flat_vals.items():
        mat[pos] = val
    for idx, (pos, (sign, v)) in enumerate(zip(free_pos, combo)):
        mat[pos] = sign * (2**v)
    det_val = int(round(np.linalg.det(mat.astype(float))))
    if abs(det_val) != 8:
        continue
    sf = smith_diag(mat.tolist())
    if sf == (1, 2, 4):
        candidates_all.append(mat.copy())
        # Check sign rule
        sign_ok = all(np.sign(mat[i, j]) == sign_for_row[i]
                       for (i, j) in free_pos)
        if sign_ok:
            candidates_sign_fixed.append(mat.copy())

print(f"\n  Total Smith=(1,2,4): {len(candidates_all)}")
print(f"  With sign rule:      {len(candidates_sign_fixed)}")
print(f"  Sign rule reduction: {len(candidates_all)} → {len(candidates_sign_fixed)} "
      f"({len(candidates_sign_fixed)/len(candidates_all)*100:.1f}%)")
print(f"  Information in sign rule: {np.log2(len(candidates_all)) - np.log2(len(candidates_sign_fixed)):.2f} bits")
print(f"  Remaining info to select physical: log₂({len(candidates_sign_fixed)}) = "
      f"{np.log2(len(candidates_sign_fixed)):.2f} bits")

# ═══════════════════════════════════════════════════════════════════
# PART 2: SCALAR INVARIANTS
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 2: SCALAR INVARIANTS (discrimination power)")
print(f"{'='*70}")

def compute_invariants(mat):
    """Compute a dictionary of scalar invariants for matrix mat."""
    m = mat.astype(float)
    inv = {}

    # Basic
    inv['tr'] = int(np.trace(mat))
    inv['det'] = int(round(np.linalg.det(m)))

    # Cofactor sum = sum of 2x2 diagonal minors
    c2 = (mat[0,0]*mat[1,1] - mat[0,1]*mat[1,0] +
          mat[0,0]*mat[2,2] - mat[0,2]*mat[2,0] +
          mat[1,1]*mat[2,2] - mat[1,2]*mat[2,1])
    inv['cofactor_sum'] = int(c2)

    # Norms
    inv['frob2'] = int(np.sum(mat**2))  # ||N||²_F
    inv['l1_free'] = sum(abs(int(mat[i, j])) for (i, j) in free_pos)
    inv['linf'] = int(np.max(np.abs(mat)))
    inv['sum_abs'] = int(np.sum(np.abs(mat)))

    # Row/column structure
    inv['row_sums'] = tuple(int(np.sum(mat[i, :])) for i in range(3))
    inv['col_sums'] = tuple(int(np.sum(mat[:, j])) for j in range(3))
    inv['row_abs_sums'] = tuple(int(np.sum(np.abs(mat[i, :]))) for i in range(3))
    inv['col_abs_sums'] = tuple(int(np.sum(np.abs(mat[:, j]))) for j in range(3))

    # Eigenvalues (sorted by real part)
    evals = np.linalg.eigvals(m)
    evals_sorted = sorted(evals, key=lambda x: x.real)
    inv['evals_real'] = tuple(round(e.real, 6) for e in evals_sorted)
    inv['evals_imag'] = tuple(round(e.imag, 6) for e in evals_sorted)
    inv['all_real'] = all(abs(e.imag) < 1e-10 for e in evals)
    inv['spectral_radius'] = round(max(abs(e) for e in evals), 6)

    # Singular values
    sv = np.linalg.svd(m, compute_uv=False)
    sv_sorted = tuple(sorted(sv, reverse=True))
    inv['sv'] = tuple(round(s, 6) for s in sv_sorted)
    inv['condition'] = round(sv_sorted[0] / sv_sorted[-1], 4) if sv_sorted[-1] > 1e-10 else float('inf')
    inv['sv_ratio_12'] = round(sv_sorted[0] / sv_sorted[1], 4) if sv_sorted[1] > 1e-10 else float('inf')

    # Permanent
    perm = 0
    for perm_idx in [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)]:
        perm += mat[0, perm_idx[0]] * mat[1, perm_idx[1]] * mat[2, perm_idx[2]]
    inv['permanent'] = int(perm)

    # N^T N and N N^T
    NtN = mat.T @ mat
    NNt = mat @ mat.T
    inv['tr_NtN'] = int(np.trace(NtN))  # = frob2
    inv['det_NtN'] = int(round(np.linalg.det(NtN.astype(float))))
    inv['tr_NNt'] = int(np.trace(NNt))

    # N²
    N2 = mat @ mat
    inv['tr_N2'] = int(np.trace(N2))
    inv['det_N2'] = int(round(np.linalg.det(N2.astype(float))))

    # N³
    N3 = N2 @ mat
    inv['tr_N3'] = int(np.round(np.trace(N3)))

    # Symmetry/antisymmetry
    S = mat + mat.T
    A = mat - mat.T
    inv['sym_frob2'] = int(np.sum(S**2))
    inv['antisym_frob2'] = int(np.sum(A**2))
    inv['sym_det'] = int(round(np.linalg.det(S.astype(float))))
    inv['antisym_rank'] = int(np.linalg.matrix_rank(A.astype(float)))

    # Mod p ranks
    for p in [2, 3, 5, 7]:
        M_mod_p = mat % p
        # Compute rank mod p via row reduction
        Mf = M_mod_p.astype(float)
        inv[f'rank_mod{p}'] = int(np.linalg.matrix_rank(Mf))

    # GCD of rows
    from math import gcd
    for i in range(3):
        row_gcd = abs(mat[i, 0])
        for j in range(1, 3):
            row_gcd = gcd(row_gcd, abs(mat[i, j]))
        inv[f'row{i}_gcd'] = row_gcd

    # GCD of columns
    for j in range(3):
        col_gcd = abs(mat[0, j])
        for i in range(1, 3):
            col_gcd = gcd(col_gcd, abs(mat[i, j]))
        inv[f'col{j}_gcd'] = col_gcd

    # v₂ multiset of free entries (sorted)
    v2_list = []
    for (i, j) in free_pos:
        n = abs(int(mat[i, j]))
        v = 0
        t = n
        while t > 1 and t % 2 == 0:
            v += 1
            t //= 2
        v2_list.append(v)
    inv['v2_multiset'] = tuple(sorted(v2_list))
    inv['v2_sum'] = sum(v2_list)
    inv['v2_prod'] = 1
    for v in v2_list:
        inv['v2_prod'] *= (2**v)

    # Diagonal product
    inv['diag_prod'] = int(mat[0,0] * mat[1,1] * mat[2,2])

    # Anti-diagonal product (fixed)
    inv['antidiag_prod'] = int(mat[0,2] * mat[1,1] * mat[2,0])

    # "Cross ratio" type invariants
    # det(N) / product_of_free_entries
    free_prod = 1
    for (i, j) in free_pos:
        free_prod *= int(mat[i, j])
    inv['free_prod'] = free_prod
    inv['det_over_free_prod'] = round(inv['det'] / free_prod, 8) if free_prod != 0 else None

    return inv


# Compute invariants for all candidates
print("\n  Computing invariants for all candidates (sign-fixed)...")
all_invs = []
phys_idx = None
for idx, mat in enumerate(candidates_sign_fixed):
    inv = compute_invariants(mat)
    all_invs.append(inv)
    if np.array_equal(mat, N_phys):
        phys_idx = idx

if phys_idx is None:
    print("  WARNING: Physical matrix not found in sign-fixed candidates!")
    # Try to find it in all candidates
    for idx, mat in enumerate(candidates_all):
        if np.array_equal(mat, N_phys):
            print(f"  Found in ALL candidates at index {idx}")
            break
else:
    print(f"  Physical matrix found at index {phys_idx}")

phys_inv = all_invs[phys_idx] if phys_idx is not None else None

# For each scalar invariant, count how many candidates share the physical value
print(f"\n  Discrimination power of each invariant:")
print(f"  {'Invariant':>25s}  {'Physical':>15s}  {'#sharing':>8s}  {'power':>8s}")
print(f"  {'─'*25}  {'─'*15}  {'─'*8}  {'─'*8}")

scalar_keys = [k for k in phys_inv.keys()
               if not isinstance(phys_inv[k], type(None))]
# Sort by discrimination power
discrim = []
for key in scalar_keys:
    phys_val = phys_inv[key]
    sharing = sum(1 for inv in all_invs if inv[key] == phys_val)
    discrim.append((sharing, key, phys_val))

discrim.sort()
for sharing, key, val in discrim:
    val_str = str(val)
    if len(val_str) > 15:
        val_str = val_str[:12] + "..."
    print(f"  {key:>25s}  {val_str:>15s}  {sharing:>8d}  "
          f"{len(candidates_sign_fixed)/sharing:>8.1f}x")

# ═══════════════════════════════════════════════════════════════════
# PART 3: BEST SINGLE DISCRIMINATORS
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 3: TOP DISCRIMINATORS")
print(f"{'='*70}")

# Show top 10
print(f"\n  Top 10 most discriminating single invariants:")
for i, (sharing, key, val) in enumerate(discrim[:10]):
    print(f"  #{i+1}: {key} = {val} → {sharing} candidates share this")

# ═══════════════════════════════════════════════════════════════════
# PART 4: PAIR COMBINATIONS
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 4: PAIRWISE COMBINATIONS")
print(f"{'='*70}")

# Find pairs that uniquely identify the physical matrix
top_keys = [key for _, key, _ in discrim[:15]]  # use top 15
best_pairs = []
for i in range(len(top_keys)):
    for j in range(i+1, len(top_keys)):
        k1, k2 = top_keys[i], top_keys[j]
        v1, v2 = phys_inv[k1], phys_inv[k2]
        sharing = sum(1 for inv in all_invs if inv[k1] == v1 and inv[k2] == v2)
        best_pairs.append((sharing, k1, k2, v1, v2))

best_pairs.sort()
print(f"\n  Top 10 pairwise discriminators:")
for i, (sharing, k1, k2, v1, v2) in enumerate(best_pairs[:10]):
    print(f"  #{i+1}: ({k1}={v1}, {k2}={v2}) → {sharing} candidates")

# Check if any pair is UNIQUE
unique_pairs = [(s, k1, k2, v1, v2) for s, k1, k2, v1, v2 in best_pairs if s == 1]
if unique_pairs:
    print(f"\n  ★ UNIQUE PAIRS FOUND: {len(unique_pairs)}")
    for sharing, k1, k2, v1, v2 in unique_pairs[:5]:
        print(f"    ({k1}={v1}, {k2}={v2}) → UNIQUE!")
else:
    print(f"\n  No unique pairs among top 15 invariants.")
    # Try triples
    print(f"  Searching triples...")
    best_triples = []
    for i in range(len(top_keys)):
        for j in range(i+1, len(top_keys)):
            for k in range(j+1, len(top_keys)):
                k1, k2, k3 = top_keys[i], top_keys[j], top_keys[k]
                v1, v2, v3 = phys_inv[k1], phys_inv[k2], phys_inv[k3]
                sharing = sum(1 for inv in all_invs
                             if inv[k1]==v1 and inv[k2]==v2 and inv[k3]==v3)
                if sharing <= 3:
                    best_triples.append((sharing, k1, k2, k3))
    best_triples.sort()
    if best_triples:
        print(f"  Top triples:")
        for s, k1, k2, k3 in best_triples[:5]:
            print(f"    ({k1}, {k2}, {k3}) → {s} candidates")

# ═══════════════════════════════════════════════════════════════════
# PART 5: EXTREMAL PROPERTY SEARCH
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 5: EXTREMAL PROPERTIES")
print(f"{'='*70}")

# For each invariant, check if the physical value is an extremum
print(f"\n  Is the physical matrix extremal in any property?")
print(f"  {'Property':>25s}  {'Physical':>10s}  {'Min':>10s}  {'Max':>10s}  {'Rank':>10s}")
print(f"  {'─'*25}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*10}")

extremal_keys = ['tr', 'cofactor_sum', 'frob2', 'permanent', 'det',
                 'condition', 'spectral_radius', 'tr_N2', 'tr_N3',
                 'sym_det', 'v2_sum', 'v2_prod', 'free_prod',
                 'diag_prod', 'l1_free', 'sum_abs', 'det_NtN',
                 'sv_ratio_12']

for key in extremal_keys:
    vals = [inv[key] for inv in all_invs if inv[key] is not None]
    if not vals or not isinstance(vals[0], (int, float)):
        continue
    phys_val = phys_inv[key]
    min_val = min(vals)
    max_val = max(vals)
    sorted_vals = sorted(vals)
    rank = sorted_vals.index(phys_val) + 1
    total = len(sorted_vals)
    is_ext = " ← MIN" if rank == 1 else (" ← MAX" if rank == total else "")
    print(f"  {key:>25s}  {str(phys_val):>10s}  {str(min_val):>10s}  "
          f"{str(max_val):>10s}  {rank:>4d}/{total:<4d}{is_ext}")

# ═══════════════════════════════════════════════════════════════════
# PART 6: RELATION TO CURVATURE MATRIX
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 6: RELATION N ↔ CURVATURE 4A")
print(f"{'='*70}")

# 4A matrix (integer entries for the anti-diagonal + rational for rest)
fourA = np.array([[4, 2, -6],
                   [7, 2, -6],
                   [-4, 0, 4]])  # 4*a_mat (approximated as int for comparison)

print(f"\n  4A = {fourA.tolist()}")
print(f"  N  = {N_phys.tolist()}")
print(f"  N + 4A = {(N_phys + fourA).tolist()}")  # = 2(δb+3a) ?

# Try: N·(4A)^T, (4A)·N^T, commutator
comm = N_phys @ fourA - fourA @ N_phys
print(f"  [N, 4A] = {comm.tolist()}")
print(f"  tr([N, 4A]) = {int(np.trace(comm))}")

# tr(N · 4A)
trNA = int(np.trace(N_phys @ fourA))
print(f"  tr(N·4A) = {trNA}")

# Compute for all candidates
trNA_vals = []
for mat in candidates_sign_fixed:
    trNA_vals.append(int(np.trace(mat @ fourA)))

phys_trNA = int(np.trace(N_phys @ fourA))
sharing_trNA = sum(1 for v in trNA_vals if v == phys_trNA)
print(f"  tr(N·4A) = {phys_trNA}, shared by {sharing_trNA}/{len(candidates_sign_fixed)}")

# det(N + 4A)
detNpA_vals = []
for mat in candidates_sign_fixed:
    detNpA_vals.append(int(round(np.linalg.det((mat + fourA).astype(float)))))
phys_detNpA = int(round(np.linalg.det((N_phys + fourA).astype(float))))
sharing_detNpA = sum(1 for v in detNpA_vals if v == phys_detNpA)
print(f"  det(N+4A) = {phys_detNpA}, shared by {sharing_detNpA}/{len(candidates_sign_fixed)}")

# det(N - 4A)
detNmA_vals = []
for mat in candidates_sign_fixed:
    detNmA_vals.append(int(round(np.linalg.det((mat - fourA).astype(float)))))
phys_detNmA = int(round(np.linalg.det((N_phys - fourA).astype(float))))
sharing_detNmA = sum(1 for v in detNmA_vals if v == phys_detNmA)
print(f"  det(N−4A) = {phys_detNmA}, shared by {sharing_detNmA}/{len(candidates_sign_fixed)}")

# ═══════════════════════════════════════════════════════════════════
# PART 7: CHARACTERISTIC POLYNOMIAL ANALYSIS
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 7: CHARACTERISTIC POLYNOMIAL χ(λ) = λ³ − tr·λ² + c₂·λ − det")
print(f"{'='*70}")

# Physical: χ(λ) = λ³ − λ² − 18λ + 8
# Group candidates by (tr, c₂) since det = ±8 is also determined
chi_groups = Counter()
for inv in all_invs:
    chi_groups[(inv['tr'], inv['cofactor_sum'], inv['det'])] += 1

print(f"\n  Number of distinct χ(λ): {len(chi_groups)}")
phys_chi = (phys_inv['tr'], phys_inv['cofactor_sum'], phys_inv['det'])
print(f"  Physical χ: tr={phys_chi[0]}, c₂={phys_chi[1]}, det={phys_chi[2]}")
print(f"  Candidates with same χ: {chi_groups[phys_chi]}")

# Eigenvalue properties
print(f"\n  Eigenvalue analysis:")
# Group by: are all eigenvalues real?
real_count = sum(1 for inv in all_invs if inv['all_real'])
print(f"  All real eigenvalues: {real_count}/{len(all_invs)}")
print(f"  Physical: all_real = {phys_inv['all_real']}")

# Among those with same χ, what varies?
same_chi = [inv for inv in all_invs
            if (inv['tr'], inv['cofactor_sum'], inv['det']) == phys_chi]
print(f"\n  Among {len(same_chi)} with same χ(λ):")
# They all have the same eigenvalues! So eigenvalues don't help further.

# What DOES vary among same-χ matrices?
varying = {}
for key in scalar_keys:
    vals = set()
    for inv in same_chi:
        v = inv[key]
        if isinstance(v, (tuple, list)):
            v = tuple(v)
        vals.add(v)
    if len(vals) > 1:
        varying[key] = len(vals)

print(f"  Properties that VARY within same χ: {len(varying)}")
for key, nvals in sorted(varying.items(), key=lambda x: x[1]):
    phys_val = phys_inv[key]
    sharing = sum(1 for inv in same_chi if inv[key] == phys_val)
    print(f"    {key}: {nvals} distinct values, physical shared by {sharing}")

# ═══════════════════════════════════════════════════════════════════
# PART 8: INVERSE MATRIX ANALYSIS
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 8: INVERSE AND ADJUGATE")
print(f"{'='*70}")

adj_phys = np.round(np.linalg.det(N_phys.astype(float)) *
                     np.linalg.inv(N_phys.astype(float))).astype(int)
print(f"\n  adj(N) = {adj_phys.tolist()}")
print(f"  N⁻¹ = adj(N)/det = adj(N)/(-8)")
print(f"  N⁻¹ entries: {(adj_phys / (-8)).tolist()}")

# Check: are adj(N) entries nice?
print(f"\n  adj(N) entries: all even? {all(v % 2 == 0 for v in adj_phys.flatten())}")
print(f"  adj(N) non-zero pattern: {(adj_phys != 0).astype(int).tolist()}")

# Does adj(N) have a nice Smith form?
sf_adj = smith_diag(adj_phys.tolist())
print(f"  Smith(adj(N)) = {sf_adj}")

# ═══════════════════════════════════════════════════════════════════
# PART 9: POLYNOMIAL COEFFICIENT MATRIX
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 9: POLYNOMIAL COEFFICIENT STRUCTURE")
print(f"{'='*70}")

# N_{sk} = n_k(I_s) where I_s = 2I₃(s) and n_k(I) = α_k + β_k·I + γ_k·I²
# I mapping: lepton→0, up→+1, down→−1
# Extract coefficients from each candidate

def extract_poly_coeffs(mat):
    """Given N with rows ordered (lepton, up, down) = (I=0, I=+1, I=-1),
    extract polynomial coefficients (α, β, γ) for each column k."""
    coeffs = np.zeros((3, 3))  # coeffs[k] = (α_k, β_k, γ_k)
    for k in range(3):
        n0 = mat[0, k]   # I=0 → lepton
        n1 = mat[1, k]   # I=+1 → up
        nm1 = mat[2, k]  # I=-1 → down
        alpha = n0
        gamma = (n1 + nm1) / 2 - n0
        beta = (n1 - nm1) / 2
        coeffs[k] = [alpha, beta, gamma]
    return coeffs

phys_coeffs = extract_poly_coeffs(N_phys)
print(f"\n  Physical polynomial coefficients C[k] = (α, β, γ):")
for k in range(3):
    print(f"    n_{coords[k]}(I) = {phys_coeffs[k,0]:.0f} + "
          f"{phys_coeffs[k,1]:.0f}·I + {phys_coeffs[k,2]:.0f}·I²")

C = phys_coeffs.T  # coefficient matrix 3×3: rows = (α, β, γ), cols = (p, q, r)
print(f"\n  Coefficient matrix C (rows: const/linear/quad, cols: p/q/r):")
print(f"    {C.astype(int).tolist()}")
print(f"    det(C) = {int(round(np.linalg.det(C)))}")
print(f"    tr(C)  = {int(np.trace(C))}")
print(f"    Smith(C) = {smith_diag(C.astype(int).tolist())}")

# Column sums of C (= Σ_k coefficients)
print(f"\n  Row sums of C: α_sum={int(sum(C[0]))}, β_sum={int(sum(C[1]))}, "
      f"γ_sum={int(sum(C[2]))}")

# Test: how many candidates have the same C properties?
print(f"\n  Computing C invariants for all candidates...")
detC_count = Counter()
trC_count = Counter()
smithC_count = Counter()
beta_sum_count = Counter()

for mat in candidates_sign_fixed:
    cc = extract_poly_coeffs(mat)
    Ci = cc.T
    detCi = int(round(np.linalg.det(Ci)))
    trCi = int(round(np.trace(Ci)))
    sfCi = smith_diag(Ci.astype(int).tolist())
    beta_sum = int(round(sum(Ci[1])))
    detC_count[detCi] += 1
    trC_count[trCi] += 1
    smithC_count[sfCi] += 1
    beta_sum_count[beta_sum] += 1

phys_detC = int(round(np.linalg.det(C)))
phys_beta_sum = int(round(sum(C[1])))
print(f"  det(C) = {phys_detC}, shared by {detC_count[phys_detC]}")
print(f"  β_sum = {phys_beta_sum}, shared by {beta_sum_count[phys_beta_sum]}")

# ═══════════════════════════════════════════════════════════════════
# PART 10: GENERATION-WISE ANALYSIS (T_g matrices)
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 10: GENERATION TRANSITION MATRICES")
print(f"{'='*70}")

# x(g) = a·g² + b·g + c, and b = (n-6a)/2, so
# x(3) = 3n/2 + c (no dependence on a!)
# x(2) = n - 2a + c
# x(1) = n/2 - 2a + c

# The THIRD generation is directly n-dependent:
# x(3)_sk = (3/2)·n_sk + c_sk
# So x(3) vector = (3/2)·N·e_s + c_s (not quite a matrix multiply, but close)

# Actually, for each sector s, x(3)_s = (3/2)·N_s + c_s where N_s is the s-th row
# So the THIRD generation mass vector is (3/2)·(row of N) + c_row.

# The KEY: third generation is LINEARLY proportional to N.
# Different N → different x(3) → different 3rd gen masses.
# This is the most direct physical meaning of S6.

# Let's see: for the candidates, what third-generation mass vectors arise?
c_vec = np.array([
    [F(5,2), F(5), F(-6)],
    [F(0), F(-4), F(-8)],
    [F(-9,2), F(-9), F(-2)]
], dtype=object)

print(f"\n  Third generation: x(3)_sk = (3/2)·n_sk + c_sk")
print(f"  Physical x(3):")
for i in range(3):
    row = []
    for j in range(3):
        x3 = F(3,2) * N_phys[i,j] + c_vec[i,j]
        row.append(str(x3))
    print(f"    {sectors[i]:>8s}: [{', '.join(row)}]")


# ═══════════════════════════════════════════════════════════════════
# PART 11: THE BIG QUESTION — UNIQUE SELECTION
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  PART 11: SEARCHING FOR A UNIQUE SELECTION PRINCIPLE")
print(f"{'='*70}")

# Strategy: test whether ANY simple function of N attains
# a unique extremum at the physical point.

# Define various "action" functionals
def compute_actions(mat, fourA):
    """Compute various candidate action functionals."""
    m = mat.astype(float)
    fA = fourA.astype(float)
    actions = {}

    # S_K_free = 2Σ(a + n/2)² for free entries
    S_free = 0
    a_vals = np.array([[1, 0.5, -1.5], [1.75, 0.5, -1.5], [-1, 0, 1]])
    for (i, j) in free_pos:
        db = a_vals[i, j] + mat[i, j] / 2
        S_free += 2 * db**2
    actions['S_K_free'] = round(S_free, 6)

    # Frobenius
    actions['frob2'] = int(np.sum(mat**2))

    # tr(N·4A^T)
    actions['trNAt'] = int(np.trace(m @ fA.T))

    # tr(N^T·4A)
    actions['trNtA'] = int(np.trace(m.T @ fA))

    # det(N+4A)
    actions['detNpA'] = int(round(np.linalg.det(m + fA)))

    # det(N-4A)
    actions['detNmA'] = int(round(np.linalg.det(m - fA)))

    # tr((N-4A)²)
    D = m - fA
    actions['trDsq'] = int(round(np.trace(D @ D)))

    # "Dirac action": tr(N²) - tr²(N)
    actions['trN2_minus_tr2'] = int(round(np.trace(m @ m) - np.trace(m)**2))

    # Entropy-like: -Σ |n_ij|·log|n_ij| for free
    S_ent = 0
    for (i, j) in free_pos:
        n = abs(int(mat[i, j]))
        if n > 0:
            S_ent -= n * np.log2(n)
    actions['neg_entropy'] = round(S_ent, 4)

    # "Generation hierarchy": |λ_max / λ_min|
    evals = np.linalg.eigvals(m)
    evals_abs = sorted(np.abs(evals))
    if evals_abs[0] > 1e-10:
        actions['eval_ratio'] = round(evals_abs[-1] / evals_abs[0], 4)
    else:
        actions['eval_ratio'] = float('inf')

    # Condition number
    sv = np.linalg.svd(m, compute_uv=False)
    actions['cond'] = round(max(sv) / min(sv), 4) if min(sv) > 1e-10 else float('inf')

    # Mixed: S_K - α·||N||² for various α
    actions['SK_minus_frob'] = round(S_free - int(np.sum(mat**2)), 4)
    actions['SK_plus_trNA'] = round(S_free + int(np.trace(m @ fA)), 4)

    # Permanent
    perm = 0
    for perm_idx in [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)]:
        perm += mat[0, perm_idx[0]] * mat[1, perm_idx[1]] * mat[2, perm_idx[2]]
    actions['permanent'] = int(perm)

    # det(N)·tr(N) = "cubic Casimir" analog
    actions['det_times_tr'] = int(round(np.linalg.det(m) * np.trace(m)))

    # tr(adj(N))
    actions['tr_adj'] = actions['trN2_minus_tr2']  # this is actually c₂

    # S_free * |det|
    actions['SK_times_det'] = round(S_free * abs(np.linalg.det(m)), 4)

    return actions


print("\n  Evaluating candidate action functionals...")
all_actions = []
for mat in candidates_sign_fixed:
    all_actions.append(compute_actions(mat, fourA))

phys_actions = compute_actions(N_phys, fourA)

print(f"\n  For each functional, rank of physical matrix:")
print(f"  {'Functional':>25s}  {'Physical':>12s}  {'Rank':>10s}  {'Ext?':>6s}")
print(f"  {'─'*25}  {'─'*12}  {'─'*10}  {'─'*6}")

for key in sorted(phys_actions.keys()):
    vals = [a[key] for a in all_actions]
    phys_val = phys_actions[key]
    sorted_vals = sorted(set(vals))
    # Find rank
    rank_min = sorted_vals.index(phys_val) + 1
    total = len(sorted_vals)
    ext = ""
    if rank_min == 1:
        ext = "MIN ★"
    elif rank_min == total:
        ext = "MAX ★"
    elif rank_min <= 3:
        ext = f"~min"
    elif rank_min >= total - 2:
        ext = f"~max"
    print(f"  {key:>25s}  {str(phys_val):>12s}  "
          f"{rank_min:>4d}/{total:<4d}  {ext:>6s}")

# ═══════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  FINAL SUMMARY")
print(f"{'='*70}")

print(f"""
  Sign rule: ε = δ(sector=lepton) → 3056 → {len(candidates_sign_fixed)} candidates
  Characteristic polynomial χ(λ) = λ³ − λ² − 18λ + 8 → {chi_groups.get(phys_chi, '?')} candidates
  
  Information budget:
    Total naive:     log₂(262144) = 18.0 bits  (6 entries × 3 bits)
    After sign rule: log₂({len(candidates_sign_fixed)}) = {np.log2(len(candidates_sign_fixed)):.2f} bits
    After χ(λ):      log₂({chi_groups.get(phys_chi, '?')}) = {np.log2(chi_groups.get(phys_chi, 1)):.2f} bits
""")
