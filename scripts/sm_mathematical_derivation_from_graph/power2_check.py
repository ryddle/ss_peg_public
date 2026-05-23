#!/usr/bin/env python3
"""Check if the lattice displacement integers n are all powers of 2."""
from fractions import Fraction as F
import math

a = {
    ('lepton','p'): F(1),    ('lepton','q'): F(1,2),  ('lepton','r'): F(-3,2),
    ('up','p'):     F(7,4),  ('up','q'):     F(1,2),   ('up','r'):     F(-3,2),
    ('down','p'):   F(-1),   ('down','q'):   F(0),     ('down','r'):   F(1),
}
b = {
    ('lepton','p'): F(-7,2), ('lepton','q'): F(-11,2), ('lepton','r'): F(15,2),
    ('up','p'):     F(-13,4),('up','q'):     F(-5,2),  ('up','r'):     F(17,2),
    ('down','p'):   F(5),    ('down','q'):   F(1),     ('down','r'):   F(-1),
}
flat = {('lepton','r'), ('up','q'), ('down','p')}
sectors = ['lepton', 'up', 'down']
coords = ['p', 'q', 'r']

# Compute n = 2(δb - a) for each entry
print("="*70)
print("  n = 2(δb − a) = number of half-steps from coset origin")
print("="*70)
print()

n_matrix = {}
for s in sectors:
    for k in coords:
        db = b[(s,k)] + 4*a[(s,k)]
        n = int(2*(db - a[(s,k)]))
        n_matrix[(s,k)] = n

# Display n matrix
print("  n matrix (sector × coord):")
print(f"  {'':>8s}    p     q     r")
for s in sectors:
    row = ""
    for k in coords:
        n = n_matrix[(s,k)]
        tag = "*" if (s,k) in flat else " "
        row += f"  {n:>+3d}{tag} "
    print(f"  {s:>8s}  {row}")
print("  (* = flat)")
print()

# Check powers of 2
print("  Power-of-2 check:")
print(f"  {'key':>12s}  {'n':>5s}  {'|n|':>5s}  {'2^k':>6s}")
all_pow2 = True
for s in sectors:
    for k in coords:
        n = n_matrix[(s,k)]
        abs_n = abs(n)
        is_pow2 = abs_n > 0 and (abs_n & (abs_n - 1)) == 0
        if is_pow2:
            exp = int(math.log2(abs_n))
            p2str = f"2^{exp}"
        else:
            p2str = "NO"
            all_pow2 = False
        tag = " [FLAT]" if (s,k) in flat else ""
        print(f"  {s:>6s}/{k}:  {n:>+5d}  {abs_n:>5d}  {p2str:>6s}{tag}")

print()
print(f"  ALL |n| are powers of 2? {all_pow2}")
print()

# Exponent matrix
print("  Exponent matrix k where |n| = 2^k:")
print(f"  {'':>8s}    p     q     r    Σrow")
for s in sectors:
    row_sum = 0
    row = ""
    for k in coords:
        n = n_matrix[(s,k)]
        abs_n = abs(n)
        exp = int(math.log2(abs_n)) if (abs_n > 0 and (abs_n & (abs_n-1)) == 0) else -1
        tag = "*" if (s,k) in flat else " "
        row += f"  {exp:>3d}{tag} "
        row_sum += exp
    print(f"  {s:>8s}  {row}  {row_sum}")

col_sums = ""
for k in coords:
    col = sum(int(math.log2(abs(n_matrix[(s,k)]))) for s in sectors
              if abs(n_matrix[(s,k)]) > 0 and (abs(n_matrix[(s,k)]) & (abs(n_matrix[(s,k)])-1)) == 0)
    col_sums += f"  {col:>3d}  "
print(f"  {'Σcol':>8s}  {col_sums}")
print()

# Sign matrix
print("  Sign matrix:")
print(f"  {'':>8s}    p     q     r")
for s in sectors:
    row = ""
    for k in coords:
        n = n_matrix[(s,k)]
        sign = "+" if n > 0 else "−"
        tag = "*" if (s,k) in flat else " "
        row += f"    {sign}{tag} "
    print(f"  {s:>8s}  {row}")
print()

# Products
print("  Product of n per row:")
for s in sectors:
    prod = 1
    for k in coords:
        prod *= n_matrix[(s,k)]
    print(f"    {s:>8s}: {prod:>+8d}  = {'+' if prod > 0 else '-'}2^{int(math.log2(abs(prod)))}" 
          if abs(prod) > 0 and (abs(prod) & (abs(prod)-1)) == 0
          else f"    {s:>8s}: {prod:>+8d}")

print("  Product of n per column:")
for k in coords:
    prod = 1
    for s in sectors:
        prod *= n_matrix[(s,k)]
    print(f"    {k}: {prod:>+8d}  = {'+' if prod > 0 else '-'}2^{int(math.log2(abs(prod)))}"
          if abs(prod) > 0 and (abs(prod) & (abs(prod)-1)) == 0
          else f"    {k}: {prod:>+8d}")

prod_all = 1
for s in sectors:
    for k in coords:
        prod_all *= n_matrix[(s,k)]
print(f"\n  Product of ALL 9: {prod_all}")
if abs(prod_all) > 0 and (abs(prod_all) & (abs(prod_all)-1)) == 0:
    print(f"  = {'+'if prod_all>0 else '-'}2^{int(math.log2(abs(prod_all)))}")
print()

# Sum of exponents
total_exp = sum(int(math.log2(abs(n_matrix[(s,k)]))) for s in sectors for k in coords
                if abs(n_matrix[(s,k)]) > 0 and (abs(n_matrix[(s,k)]) & (abs(n_matrix[(s,k)])-1)) == 0)
print(f"  Sum of ALL exponents: {total_exp}")
free_exp = sum(int(math.log2(abs(n_matrix[(s,k)]))) for s in sectors for k in coords
               if (s,k) not in flat and abs(n_matrix[(s,k)]) > 0 and (abs(n_matrix[(s,k)]) & (abs(n_matrix[(s,k)])-1)) == 0)
flat_exp = sum(int(math.log2(abs(n_matrix[(s,k)]))) for s in sectors for k in coords
               if (s,k) in flat and abs(n_matrix[(s,k)]) > 0 and (abs(n_matrix[(s,k)]) & (abs(n_matrix[(s,k)])-1)) == 0)
print(f"  Free exponents sum:   {free_exp}")
print(f"  Flat exponents sum:   {flat_exp}")
print()

# Connection to quantum numbers
print("="*70)
print("  CONNECTION TO QUANTUM NUMBERS")
print("="*70)
print()
print("  Recall: the 3 coords are (p, q, r) = coefficients of (ln2, lnR3, ln3)")
print("  p ↔ 2I3 (isospin), q ↔ Nc (color), r ↔ generation structure")
print()

# Physical quantum numbers
I3 = {'lepton': {1: F(0), 2: F(0), 3: F(0)},   # e,μ,τ: I3=0 (right-handed)
      'up':     {1: F(1,2), 2: F(1,2), 3: F(1,2)},
      'down':   {1: F(-1,2), 2: F(-1,2), 3: F(-1,2)}}
Nc = {'lepton': 1, 'up': 3, 'down': 3}
Q  = {'lepton': -1, 'up': F(2,3), 'down': F(-1,3)}

print("  Sector properties and n pattern:")
for s in sectors:
    ns = [n_matrix[(s,k)] for k in coords]
    exps = [int(math.log2(abs(n))) if abs(n) > 0 and (abs(n) & (abs(n)-1))==0 else -1 for n in ns]
    print(f"    {s:>8s}: Q={float(Q[s]):>+5.2f}  Nc={Nc[s]}  "
          f"n=({ns[0]:>+3d},{ns[1]:>+3d},{ns[2]:>+3d})  "
          f"exp=({exps[0]},{exps[1]},{exps[2]})")
print()

# Is there a formula? n_k = f(Q, Nc, I3, k)?
print("  Testing: is exp(|n|) = |charge numerator| or similar?")
print(f"    lepton: |Q|=1,  Nc=1, exponents = (0, 3, 2)")
print(f"    up:     |Q|=2/3,Nc=3, exponents = (2, 1, 3)")  
print(f"    down:   |Q|=1/3,Nc=3, exponents = (1, 0, 2)")
print()

# Note the flat coordinate positions
print("  Flat coords and their exponents:")
for s,k in sorted(flat):
    n = n_matrix[(s,k)]
    exp = int(math.log2(abs(n)))
    print(f"    {s:>6s}/{k}: n={n:>+3d}, exp={exp}")
print()

# Check: do the free exponents satisfy some constraint?
print("  Free n values (sorted by |n|):")
free_ns = [(s, k, n_matrix[(s,k)]) for s in sectors for k in coords if (s,k) not in flat]
free_ns.sort(key=lambda x: abs(x[2]))
for s, k, n in free_ns:
    exp = int(math.log2(abs(n)))
    print(f"    {s:>6s}/{k}: n={n:>+3d} = {'+' if n > 0 else '-'}2^{exp}")

print()

# Check if n relates to δb directly through 2I3, Nc, g
# Recall generating function: x_k(g) = f(2I3, Nc, g) along each coord
# The free δb values are the "initial velocities"
# n = 2(δb - a) is the lattice displacement

# Let's check ratios
print("  Ratios between sectors (same coord):")
for k in coords:
    ns = [n_matrix[(s,k)] for s in sectors]
    flat_in_col = any((s,k) in flat for s in sectors)
    if not flat_in_col:
        print(f"    {k}: lepton:up:down = {ns[0]}:{ns[1]}:{ns[2]}")
    else:
        print(f"    {k}: {ns[0]}:{ns[1]}:{ns[2]} (contains flat)")
