"""Quick check: is max_nt = h - 1 (Coxeter - 1) for simply-laced algebras?"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from ramanujan_extended import (roots_A, roots_B, roots_C, roots_D,
    roots_G2, roots_F4, roots_E6, roots_E7, roots_E8,
    adjacency_from_roots, spectral_analysis)

cox_e = {'G2': 6, 'F4': 12, 'E6': 12, 'E7': 18, 'E8': 30}

print("Checking max_nt vs h-1 (Coxeter - 1):")
print(f"{'Alg':>10s}  {'h':>3s}  {'h-1':>4s}  {'max_nt':>10s}  {'diff':>8s}  exact?  lacing")
print("-" * 65)

for n in range(1, 13):
    r, roots, simple = roots_A(n)
    A = adjacency_from_roots(r, roots, simple)
    s = spectral_analysis(A)
    h = n + 1
    d = s['max_nontrivial'] - (h - 1)
    ex = "YES" if abs(d) < 0.01 else "no"
    print(f"  su({n+1:2d})    {h:3d}  {h-1:4d}  {s['max_nontrivial']:10.4f}  {d:+8.4f}  {ex:>5s}   SL")

for n in range(2, 9):
    r, roots, simple = roots_B(n)
    A = adjacency_from_roots(r, roots, simple)
    s = spectral_analysis(A)
    h = 2 * n
    d = s['max_nontrivial'] - (h - 1)
    ex = "YES" if abs(d) < 0.01 else "no"
    print(f"  so({2*n+1:2d})   {h:3d}  {h-1:4d}  {s['max_nontrivial']:10.4f}  {d:+8.4f}  {ex:>5s}   NSL")

for n in range(3, 10):
    r, roots, simple = roots_D(n)
    A = adjacency_from_roots(r, roots, simple)
    s = spectral_analysis(A)
    h = 2 * (n - 1)
    d = s['max_nontrivial'] - (h - 1)
    ex = "YES" if abs(d) < 0.01 else "no"
    print(f"  so({2*n:2d})   {h:3d}  {h-1:4d}  {s['max_nontrivial']:10.4f}  {d:+8.4f}  {ex:>5s}   SL")

for name, fn, rnk in [('G2', roots_G2, 2), ('F4', roots_F4, 4),
                       ('E6', roots_E6, 6), ('E7', roots_E7, 7), ('E8', roots_E8, 8)]:
    r, roots, simple = fn()
    A = adjacency_from_roots(r, roots, simple)
    s = spectral_analysis(A)
    h = cox_e[name]
    d = s['max_nontrivial'] - (h - 1)
    ex = "YES" if abs(d) < 0.01 else "no"
    lac = "SL" if name.startswith('E') else "NSL"
    print(f"  {name:>6s}   {h:3d}  {h-1:4d}  {s['max_nontrivial']:10.4f}  {d:+8.4f}  {ex:>5s}   {lac}")
