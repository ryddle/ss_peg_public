#!/usr/bin/env python3
"""
complete_action.py
==================
THE COMPLETE ACTION — can we reconstruct ALL 27 coordinates 
(= 9 fermion masses) from ONLY the rules discovered?

RULES:
  R1. Generating function F(2I₃, Nc, g) gives a_k and c_k 
      (curvatures and initial positions) — 18 of 27 params
  R2. Flatness: one (s,k) per sector has b = -5a  (3 constraints)
      → 3 of 9 b_k fixed → 24 of 27 known
  R3. Kinetic optimum: b* = -4a defines the reference
      δb = b + 4a is the free deviation
  R4. 2-adic quantization: n = 2(δb - a) = ±2^v  for free entries
      where v ∈ {0,1,2,3}
  R5. The specific (sign, v) assignment per (sector, coord)

QUESTION: Do R1-R4 constrain enough, or do we need R5 as INPUT?
If R5 is needed, how much information does it carry?
"""
from fractions import Fraction as F
import numpy as np

# ═══════════════════════════════════════════════════════════════════
# KNOWN FROM GRAPH TOPOLOGY (R1: generating function)
# ═══════════════════════════════════════════════════════════════════
sectors = ['lepton', 'up', 'down']
coords  = ['p', 'q', 'r']
flat = {('lepton','r'), ('up','q'), ('down','p')}

# These come from the generating function F(2I₃, Nc, g)
a_data = {
    ('lepton','p'): F(1),    ('lepton','q'): F(1,2),  ('lepton','r'): F(-3,2),
    ('up','p'):     F(7,4),  ('up','q'):     F(1,2),   ('up','r'):     F(-3,2),
    ('down','p'):   F(-1),   ('down','q'):   F(0),     ('down','r'):   F(1),
}
c_data = {
    ('lepton','p'): F(5,2),  ('lepton','q'): F(5),     ('lepton','r'): F(-6),
    ('up','p'):     F(0),    ('up','q'):     F(-4),     ('up','r'):     F(-8),
    ('down','p'):   F(-9,2), ('down','q'):   F(-9),     ('down','r'):   F(-2),
}
# Physical b values (what we want to DERIVE)
b_physical = {
    ('lepton','p'): F(-7,2), ('lepton','q'): F(-11,2), ('lepton','r'): F(15,2),
    ('up','p'):     F(-13,4),('up','q'):     F(-5,2),  ('up','r'):     F(17,2),
    ('down','p'):   F(5),    ('down','q'):   F(1),     ('down','r'):   F(-1),
}

print("=" * 70)
print("  THE COMPLETE ACTION: FROM RULES TO MASSES")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════
# STEP 1: Apply R2 (flatness) → fixes 3 b_k
# ═══════════════════════════════════════════════════════════════════
print("\n§1. FLATNESS CONSTRAINT (R2): b_flat = -5a")
print("-" * 70)

b_derived = {}
for s in sectors:
    for k in coords:
        if (s,k) in flat:
            b_derived[(s,k)] = -5 * a_data[(s,k)]
            check = b_derived[(s,k)] == b_physical[(s,k)]
            print(f"  {s:>6s}/{k} [FLAT]: b = -5×({a_data[(s,k)]}) = {b_derived[(s,k)]}  "
                  f"(physical: {b_physical[(s,k)]})  {'✓' if check else '✗'}")

# ═══════════════════════════════════════════════════════════════════
# STEP 2: Apply R3+R4 (2-adic quantization) → restricts 6 b_k
# ═══════════════════════════════════════════════════════════════════
print("\n§2. 2-ADIC QUANTIZATION (R4): n = ±2^v, δb = a + n/2")
print("-" * 70)

# The 6 free entries need (sign, v₂) = 6 × 3 = 18 bits.
# But first: what values of (sign, v) are POSSIBLE given
# the lattice constraint δb ∈ a + Z/2 (i.e., n ∈ Z)?

# n = ±2^v must be an integer. It always is. So R4 gives:
# δb = a ± 2^v / 2 = a ± 2^(v-1)
# b = δb - 4a = -3a ± 2^(v-1)

# The CANDIDATE set for free entries:
# v ∈ {0,1,2,3}, sign ∈ {+,-} → 8 candidates per entry
# 8^6 = 262144 total combinations

print("\n  For each free (s,k), the physical (sign, v₂):")
free_entries = []
for s in sectors:
    for k in coords:
        if (s,k) not in flat:
            a = a_data[(s,k)]
            b_phys = b_physical[(s,k)]
            db = b_phys + 4*a  # δb
            n = 2*(db - a)     # n = 2(δb - a)
            n_int = int(n)
            sign = "+" if n_int > 0 else "−"
            v = 0
            temp = abs(n_int)
            while temp > 1 and temp % 2 == 0:
                v += 1
                temp //= 2
            free_entries.append((s, k, n_int, sign, v))
            print(f"  {s:>6s}/{k}: n = {n_int:>+3d} = {sign}2^{v}  "
                  f"→ δb = {db} → b = {b_phys}")

# ═══════════════════════════════════════════════════════════════════
# STEP 3: What is R5? The specific assignment matrix
# ═══════════════════════════════════════════════════════════════════
print("\n§3. THE ASSIGNMENT MATRIX (R5)")
print("-" * 70)

# v₂ matrix for free entries:
# lepton: (0, 3, --)
# up:     (2, --, 3)
# down:   (--, 1, 2)
# sign matrix for free entries:
# lepton: (-, -, --)
# up:     (+, --, +)
# down:   (--, +, +)

print("\n  Valuation v₂ (free only):")
print(f"  {'':>8s}     p    q    r")
for s in sectors:
    row = ""
    for k in coords:
        if (s,k) in flat:
            row += "   -- "
        else:
            n = int(2*(b_physical[(s,k)] + 4*a_data[(s,k)] - a_data[(s,k)]))
            v = 0; t = abs(n)
            while t > 1 and t % 2 == 0: v += 1; t //= 2
            row += f"   {v:>2d} "
    print(f"  {s:>6s}  {row}")

print("\n  Sign ε (free only, 1=negative):")
print(f"  {'':>8s}     p    q    r")
for s in sectors:
    row = ""
    for k in coords:
        if (s,k) in flat:
            row += "   -- "
        else:
            n = int(2*(b_physical[(s,k)] + 4*a_data[(s,k)] - a_data[(s,k)]))
            e = 1 if n < 0 else 0
            row += f"   {e:>2d} "
    print(f"  {s:>6s}  {row}")

# ═══════════════════════════════════════════════════════════════════
# STEP 4: How many bits does R5 carry?
# ═══════════════════════════════════════════════════════════════════
print("\n§4. INFORMATION CONTENT OF R5")
print("-" * 70)

# Naive: 6 entries × (1 bit sign + 2 bits v₂) = 18 bits
# But are there constraints that reduce this?

# v₂ multiset is {0, 1, 2, 2, 3, 3}
# Sign pattern: 2 negative (both in lepton), 4 positive

# Count: how many distinct N matrices have Smith form (1,2,4)?
# That is the real question.

print("\n  Naive information: 6 × (1 sign + 2 bits for v₂∈{0..3}) = 18 bits")
print("  But v₂ only takes values {0,1,2,3} (2 bits) → confirmed 18 bits")
print()

# Let's check: given Smith form (1,2,4) and det = ±8,
# how many integer 3×3 matrices N exist with:
#   - entries ±2^v (v∈{0..3}) off the anti-diagonal
#   - entries -4a on the anti-diagonal (fixed)
#   - Smith form = (1,2,4)?

print("  CONSTRAINT: Smith(N) = (1,2,4) → |det| = 8")
print("  Testing: how many (sign,v) combinations give |det| = 8...")

count_det8 = 0
count_smith124 = 0
total = 0

# Free positions: (0,0), (0,1), (1,0), (1,2), (2,1), (2,2)
# Flat positions: (0,2)=6, (1,1)=-2, (2,0)=4
free_pos = [(0,0), (0,1), (1,0), (1,2), (2,1), (2,2)]
flat_vals = {(0,2): 6, (1,1): -2, (2,0): 4}

from itertools import product as iprod

for combo in iprod(iprod([-1,1], [0,1,2,3]), repeat=6):
    total += 1
    mat = np.zeros((3,3), dtype=int)
    # Fill flat
    for pos, val in flat_vals.items():
        mat[pos] = val
    # Fill free
    for idx, (pos, (sign, v)) in enumerate(zip(free_pos, combo)):
        mat[pos] = sign * (2**v)
    
    det_val = int(round(np.linalg.det(mat.astype(float))))
    if abs(det_val) == 8:
        count_det8 += 1

print(f"  Total (sign,v) combinations: {total}")
print(f"  Combinations with |det| = 8: {count_det8}")
print(f"  Fraction: {count_det8}/{total} = {count_det8/total:.4f}")

# Now check which of those have Smith = (1,2,4)
print(f"\n  Checking Smith form for the {count_det8} candidates...")

def smith_diag(M):
    """Compute Smith normal form diagonal of integer matrix."""
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

smith_counts = {}
physical_found = False
physical_combo = None

for combo in iprod(iprod([-1,1], [0,1,2,3]), repeat=6):
    mat = np.zeros((3,3), dtype=int)
    for pos, val in flat_vals.items():
        mat[pos] = val
    for idx, (pos, (sign, v)) in enumerate(zip(free_pos, combo)):
        mat[pos] = sign * (2**v)
    
    det_val = int(round(np.linalg.det(mat.astype(float))))
    if abs(det_val) == 8:
        sf = smith_diag(mat.tolist())
        smith_counts[sf] = smith_counts.get(sf, 0) + 1
        
        # Check if this is the physical matrix
        if np.array_equal(mat, np.array([[-1,-8,6],[4,-2,8],[4,2,4]])):
            physical_found = True
            physical_combo = combo

print(f"\n  Smith form distribution among |det|=8 matrices:")
for sf, cnt in sorted(smith_counts.items()):
    marker = " ← PHYSICAL" if sf == (1,2,4) else ""
    print(f"    Smith = {sf}: {cnt} matrices{marker}")

total_124 = smith_counts.get((1,2,4), 0)

print(f"\n  Physical matrix found: {physical_found}")
print(f"\n  RESULT: {total_124} out of {total} have Smith = (1,2,4)")
print(f"  Information to select physical: log₂({total_124}) = {np.log2(total_124):.2f} bits")
print(f"  (vs naive 18 bits = log₂({total}) = {np.log2(total):.2f})")

# ═══════════════════════════════════════════════════════════════════
# STEP 5: THE COMPLETE ACTION — formalize it
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  §5. THE COMPLETE ACTION S")
print(f"{'='*70}")

def x_at(s, k, g, b_dict):
    return a_data[(s,k)]*g**2 + b_dict[(s,k)]*g + c_data[(s,k)]

# Compute S_K
S_K = F(0)
for s in sectors:
    for k in coords:
        for g in [1, 2]:
            dx = x_at(s, k, g+1, b_physical) - x_at(s, k, g, b_physical)
            S_K += dx**2

print(f"\n  S_K (kinetic) = Σ(Δx)² = {S_K} = {float(S_K):.4f}")

# Compute V_lattice at each generation
V_latt = F(0)
for s in sectors:
    for k in coords:
        for g in [1, 2, 3]:
            x = x_at(s, k, g, b_physical)
            # Distance to nearest half-integer: d = x - round(2x)/2
            x2 = 2*x
            d = x2 - round(float(x2))
            V_latt += F(d)**2 / 4  # d²(x, Z/2) = (d/2)²

print(f"  V_lattice = Σ d²(x(g), Z/2) = {V_latt} = {float(V_latt):.6f}")
print(f"  (Should be 0 since all x(g) ∈ Z/2)")

# Decompose S_K
S_curv = F(0)
S_vel = F(0)
for s in sectors:
    for k in coords:
        a = a_data[(s,k)]
        db = b_physical[(s,k)] + 4*a
        S_curv += a**2
        S_vel += db**2

S_curv = F(8,3) * S_curv  # Wait, let me recompute properly

# S_K = Σ_{s,k} [(3a+b)² + (5a+b)²]
# = Σ [(3a+b)² + (5a+b)²]
# Let δb = b+4a, then b = δb-4a
# 3a+b = 3a + δb - 4a = δb - a
# 5a+b = 5a + δb - 4a = δb + a
# So S_K = Σ [(δb-a)² + (δb+a)²] = Σ [2δb² + 2a²] = 2Σa² + 2Σδb²

S_2a2 = F(0)
S_2db2 = F(0)
for s in sectors:
    for k in coords:
        a = a_data[(s,k)]
        db = b_physical[(s,k)] + 4*a
        S_2a2 += 2 * a**2
        S_2db2 += 2 * db**2

print(f"\n  S_K decomposition:")
print(f"    S_K = 2Σa² + 2Σδb² = {S_2a2} + {S_2db2} = {S_2a2 + S_2db2}")
print(f"    Check: {S_2a2 + S_2db2} == {S_K}: {S_2a2 + S_2db2 == S_K}")

# Separate flat and free contributions to S_vel
S_vel_flat = F(0)
S_vel_free = F(0)
for s in sectors:
    for k in coords:
        db = b_physical[(s,k)] + 4*a_data[(s,k)]
        if (s,k) in flat:
            S_vel_flat += 2 * db**2
        else:
            S_vel_free += 2 * db**2

print(f"\n    2Σδb² breakdown:")
print(f"      Flat (forced):  {S_vel_flat} = {float(S_vel_flat):.4f}")
print(f"      Free (chosen):  {S_vel_free} = {float(S_vel_free):.4f}")
print(f"      Total:          {S_vel_flat + S_vel_free}")

# ═══════════════════════════════════════════════════════════════════
# STEP 6: CAN WE DERIVE R5 FROM AN ACTION?
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  §6. DOES MINIMIZING S SELECT THE PHYSICAL POINT?")
print(f"{'='*70}")

# The free part of S is 2Σ(δb_free)² = 2Σ(a + n/2)² where n = ±2^v
# Physical values:
print(f"\n  Free δb values and their S_vel contributions:")
print(f"  {'key':>12s}  {'a':>6s}  {'n':>4s}  {'δb':>6s}  {'2δb²':>8s}")
total_free_S = F(0)
for s in sectors:
    for k in coords:
        if (s,k) not in flat:
            a = a_data[(s,k)]
            db = b_physical[(s,k)] + 4*a
            n = int(2*(db - a))
            contrib = 2 * db**2
            total_free_S += contrib
            print(f"  {s:>6s}/{k}:  {str(a):>6s}  {n:>+4d}  {str(db):>6s}  {str(contrib):>8s}")
print(f"  {'Total':>12s}  {'':>6s}  {'':>4s}  {'':>6s}  {str(total_free_S):>8s}")

# Now: what is the MINIMUM of 2Σ(a + n/2)² over all n ∈ ±2^v, v∈{0..3}?
print(f"\n  For each free entry, find n that minimizes 2(a + n/2)²:")
min_total = F(0)
min_assignment = []
for s in sectors:
    for k in coords:
        if (s,k) not in flat:
            a = a_data[(s,k)]
            best_n = None
            best_val = None
            for sign in [-1, 1]:
                for v in range(4):
                    n = sign * (2**v)
                    db = a + F(n, 2)
                    val = 2 * db**2
                    if best_val is None or val < best_val:
                        best_val = val
                        best_n = n
            min_total += best_val
            db_min = a + F(best_n, 2)
            
            # Physical
            db_phys = b_physical[(s,k)] + 4*a
            n_phys = int(2*(db_phys - a))
            val_phys = 2 * db_phys**2
            
            match = "✓" if best_n == n_phys else "✗"
            print(f"  {s:>6s}/{k}: min at n={best_n:>+3d} (2δb²={str(best_val):>8s})  "
                  f"phys n={n_phys:>+3d} (2δb²={str(val_phys):>8s})  {match}")
            min_assignment.append((s, k, best_n, n_phys))

print(f"\n  Min Σ(free) = {min_total} = {float(min_total):.4f}")
print(f"  Physical Σ(free) = {total_free_S} = {float(total_free_S):.4f}")
print(f"  Excess = {total_free_S - min_total} = {float(total_free_S - min_total):.4f}")

matches = sum(1 for _,_,bn,pn in min_assignment if bn == pn)
print(f"\n  Entries matching minimum: {matches}/6")

# ═══════════════════════════════════════════════════════════════════
# STEP 7: WHAT ABOUT det(N) = ±8?
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  §7. MINIMIZING S SUBJECT TO Smith(N) = (1,2,4)")
print(f"{'='*70}")

# Among the matrices with Smith = (1,2,4), which minimizes S_vel_free?
print(f"\n  Scanning all {total_124} matrices with Smith=(1,2,4)...")
print(f"  Finding the one with minimum Σ(δb_free)²...")

min_S_smith = None
min_mat_smith = None
phys_S = float(total_free_S)
phys_rank = None
all_S_vals = []

for combo in iprod(iprod([-1,1], [0,1,2,3]), repeat=6):
    mat = np.zeros((3,3), dtype=int)
    for pos, val in flat_vals.items():
        mat[pos] = val
    for idx, (pos, (sign, v)) in enumerate(zip(free_pos, combo)):
        mat[pos] = sign * (2**v)
    
    det_val = int(round(np.linalg.det(mat.astype(float))))
    if abs(det_val) != 8:
        continue
    
    sf = smith_diag(mat.tolist())
    if sf != (1, 2, 4):
        continue
    
    # Compute S_vel_free for this assignment
    S_free = F(0)
    valid = True
    for idx, pos in enumerate(free_pos):
        i, j = pos
        s = sectors[i]
        k = coords[j]
        n = mat[i,j]
        a = a_data[(s,k)]
        db = a + F(n, 2)
        S_free += 2 * db**2
    
    all_S_vals.append((float(S_free), mat.copy()))
    
    if min_S_smith is None or S_free < min_S_smith:
        min_S_smith = S_free
        min_mat_smith = mat.copy()

all_S_vals.sort(key=lambda x: x[0])

print(f"\n  Minimum S_vel_free with Smith=(1,2,4): {float(min_S_smith):.4f}")
print(f"  Physical S_vel_free:                    {phys_S:.4f}")

# Find rank of physical
for rank, (sv, mat) in enumerate(all_S_vals):
    if np.array_equal(mat, np.array([[-1,-8,6],[4,-2,8],[4,2,4]])):
        phys_rank = rank + 1
        break

print(f"  Physical matrix rank: #{phys_rank} out of {len(all_S_vals)}")
print(f"  Is physical = minimum? {'YES ✓' if phys_rank == 1 else 'NO ✗'}")

# Show top 5
print(f"\n  Top 5 by S_vel_free:")
for i, (sv, mat) in enumerate(all_S_vals[:5]):
    marker = " ← PHYSICAL" if np.array_equal(mat, np.array([[-1,-8,6],[4,-2,8],[4,2,4]])) else ""
    print(f"    #{i+1}: S = {sv:>8.4f}  N = {mat.tolist()}{marker}")

# Show the physical one if not in top 5
if phys_rank > 5:
    sv, mat = all_S_vals[phys_rank-1]
    print(f"    ...")
    print(f"    #{phys_rank}: S = {sv:>8.4f}  N = {mat.tolist()} ← PHYSICAL")

# ═══════════════════════════════════════════════════════════════════
# STEP 8: TRY OTHER ACTIONS
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  §8. ALTERNATIVE ACTIONS")
print(f"{'='*70}")

# Try: minimize TOTAL S_K (including curvature)
# Or: minimize |det|, or tr(N^T N), etc.

# S_total = 2Σa² + 2Σδb² (a are fixed, so S_total ∝ S_free + const)
# So S_K minimum ⟺ S_free minimum. Same ranking.

# Try: minimize Frobenius norm ||N||²
print(f"\n  Alternative 1: minimize ||N||² = Σ n²")
norms = []
for combo in iprod(iprod([-1,1], [0,1,2,3]), repeat=6):
    mat = np.zeros((3,3), dtype=int)
    for pos, val in flat_vals.items():
        mat[pos] = val
    for idx, (pos, (sign, v)) in enumerate(zip(free_pos, combo)):
        mat[pos] = sign * (2**v)
    det_val = int(round(np.linalg.det(mat.astype(float))))
    if abs(det_val) != 8:
        continue
    sf = smith_diag(mat.tolist())
    if sf != (1, 2, 4):
        continue
    frob = int(np.sum(mat**2))
    norms.append((frob, mat.copy()))

norms.sort(key=lambda x: x[0])
phys_N = np.array([[-1,-8,6],[4,-2,8],[4,2,4]])
phys_frob = int(np.sum(phys_N**2))

for rank, (f, mat) in enumerate(norms):
    if np.array_equal(mat, phys_N):
        print(f"  Physical ||N||² = {phys_frob}, rank #{rank+1}/{len(norms)}")
        break
print(f"  Min ||N||² = {norms[0][0]}, Max = {norms[-1][0]}")

# Try: minimize |tr(N)|
print(f"\n  Alternative 2: tr(N) = {int(np.trace(phys_N))}")
traces = [(abs(int(np.trace(m))), int(np.trace(m)), m.tolist()) for _, m in norms]
traces.sort(key=lambda x: (x[0], x[1]))
for i, (at, t, m) in enumerate(traces[:3]):
    marker = " ← PHYS" if m == phys_N.tolist() else ""
    print(f"    #{i+1}: tr = {t}{marker}")

# Try: maximize |det(N)| ... all have |det|=8, so tied

# Try: some entropy/info measure
print(f"\n  Alternative 3: minimize Σ|n|(free entries only)")
l1_norms = []
for sv, mat in all_S_vals:
    l1 = sum(abs(int(mat[pos])) for pos in free_pos)
    l1_norms.append((l1, sv, mat))
l1_norms.sort(key=lambda x: (x[0], x[1]))
for rank, (l1, sv, mat) in enumerate(l1_norms):
    if np.array_equal(mat, phys_N):
        print(f"  Physical Σ|n| = {l1}, rank #{rank+1}/{len(l1_norms)}")
        break
print(f"  Min Σ|n| = {l1_norms[0][0]}")

# ═══════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  FINAL: THE COMPLETE ACTION")
print(f"{'='*70}")

print(f"""
  The complete specification of the mass spectrum requires:

  STRUCTURAL (from graph topology, 0 free parameters):
  ─────────────────────────────────────────────────────
  S1. Generating function F(2I₃, Nc, g)     → a_k, c_k  (18 params)
  S2. Flatness constraint (anti-diagonal)    → 3 of 9 b_k (b=-5a)
  S3. Lattice quantization: x(g) ∈ Z/2      → g ∈ {{1,2,3}}

  ACTION (variational, reduces freedom):
  ─────────────────────────────────────────────────────
  S4. 2-adic quantization: n ∈ {{±2^v : v=0..3}}     → 8 candidates/entry
  S5. Smith constraint: Smith(N) = (1,2,4)             → ~{total_124} matrices

  SELECTION (not yet derived from a principle):
  ─────────────────────────────────────────────────────
  S6. The specific (sign, v₂) assignment     → selects 1 from ~{total_124}
      Carries ~log₂({total_124}) ≈ {np.log2(total_124):.1f} bits of information.

  STATUS:
  • S1-S3 are PROVEN (from the graph).
  • S4 is OBSERVED (empirical: all free n are ±2^v).
  • S5 is OBSERVED (Smith form computed).
  • S6 is NOT derived — it is the remaining free information.

  The question "do we have the complete action?" becomes:
  Is there a principle that selects the physical matrix 
  from the ~{total_124} Smith=(1,2,4) candidates?
  
  Minimizing S_K does NOT work (physical is #{phys_rank}/{len(all_S_vals)}).
  Minimizing ||N||² does NOT work either.
  
  The ~{np.log2(total_124):.1f} bits in S6 are the IRREDUCIBLE free information 
  of the mass spectrum within this framework.
""")
