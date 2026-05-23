"""
Reconcile the two R values for su(3):
  - R_CW = 0.552 (Cartan-Weyl root system graph)
  - R_adj = 0.436 (Physics/Gell-Mann adjoint graph)
Build both graphs explicitly, verify spectra, find exact algebraic forms.
"""
import numpy as np
from itertools import combinations

np.set_printoptions(precision=8, linewidth=120)

sep = '=' * 70

# ==============================================================
# GRAPH 1: Cartan-Weyl root system graph for A_2 = su(3)
# ==============================================================

print(sep)
print('  GRAPH 1: Cartan-Weyl root system graph for SU(3)')
print(sep)

rank = 2
# A_2 roots in R^3: e_i - e_j for i != j
roots = []
for i in range(3):
    for j in range(3):
        if i != j:
            r = [0, 0, 0]
            r[i] = 1
            r[j] = -1
            roots.append(tuple(r))

simple = [(1, -1, 0), (0, 1, -1)]

print(f'\n  Roots ({len(roots)}):')
for i, r in enumerate(roots):
    print(f'    alpha_{i} = {r}')

# Build adjacency
dim = rank + len(roots)  # 2 + 6 = 8
roots_arr = [np.array(r, dtype=float) for r in roots]
simple_arr = [np.array(s, dtype=float) for s in simple]
root_set = set(roots)

A_cw = np.zeros((dim, dim))

# H_i -- E_alpha: iff alpha . simple_i != 0
for i in range(rank):
    for j in range(len(roots)):
        dot = np.dot(roots_arr[j], simple_arr[i])
        if abs(dot) > 1e-10:
            A_cw[i, rank + j] = 1
            A_cw[rank + j, i] = 1

# E_alpha -- E_beta: iff alpha+beta in roots OR beta = -alpha
for j in range(len(roots)):
    for k in range(j + 1, len(roots)):
        s = tuple(int(roots_arr[j][l] + roots_arr[k][l]) for l in range(3))
        if np.linalg.norm(roots_arr[j] + roots_arr[k]) < 1e-10:
            A_cw[rank + j, rank + k] = 1
            A_cw[rank + k, rank + j] = 1
        elif s in root_set:
            A_cw[rank + j, rank + k] = 1
            A_cw[rank + k, rank + j] = 1

degs_cw = A_cw.sum(axis=1).astype(int)
edges_cw = int(A_cw.sum() / 2)
eigs_cw = np.sort(np.linalg.eigvalsh(A_cw))[::-1]

print(f'\n  Vertices: {dim}')
print(f'  Edges: {edges_cw}')
print(f'  Degree sequence: {degs_cw}')
print(f'  d_mean = {degs_cw.mean():.4f}')
print(f'  q = d_mean - 1 = {degs_cw.mean() - 1:.4f}')
print(f'\n  Eigenvalues: {eigs_cw}')
print(f'  lambda_1 (Perron) = {eigs_cw[0]:.8f}')
print(f'  lambda_2 = {eigs_cw[1]:.8f}')
print(f'  lambda_min = {eigs_cw[-1]:.8f}')

lam_nt_cw = max(abs(eigs_cw[1]), abs(eigs_cw[-1]))
q_cw = degs_cw.mean() - 1
R_cw = lam_nt_cw / (2 * np.sqrt(q_cw))
print(f'\n  max|lambda_nt| = {lam_nt_cw:.8f}')
print(f'  2*sqrt(q) = {2*np.sqrt(q_cw):.8f}')
print(f'  R(CW) = {R_cw:.6f}')

# ==============================================================
# GRAPH 2: Physics (Gell-Mann) adjoint graph
# ==============================================================

print(f'\n{sep}')
print('  GRAPH 2: Physics (Gell-Mann) adjoint graph for SU(3)')
print(sep)

triplets = [
    (0,1,2), (0,3,6), (0,4,5), (1,3,5),
    (1,4,6), (2,3,4), (2,5,6), (3,4,7), (5,6,7),
]
A_phys = np.zeros((8, 8))
for (a, b, c) in triplets:
    A_phys[a, b] = A_phys[b, a] = 1
    A_phys[a, c] = A_phys[c, a] = 1
    A_phys[b, c] = A_phys[c, b] = 1

degs_phys = A_phys.sum(axis=1).astype(int)
edges_phys = int(A_phys.sum() / 2)
eigs_phys = np.sort(np.linalg.eigvalsh(A_phys))[::-1]

print(f'\n  Vertices: 8')
print(f'  Edges: {edges_phys}')
print(f'  Degree sequence: {degs_phys}')
print(f'  d_mean = {degs_phys.mean():.4f}')

lam_nt_phys = max(abs(eigs_phys[1]), abs(eigs_phys[-1]))
q_phys = degs_phys.mean() - 1
R_phys = lam_nt_phys / (2 * np.sqrt(q_phys))
print(f'\n  Eigenvalues: {eigs_phys}')
print(f'  max|lambda_nt| = {lam_nt_phys:.8f}')
print(f'  R(phys) = {R_phys:.6f}')

# ==============================================================
# EXACT ALGEBRAIC FORM OF R(CW)
# ==============================================================

print(f'\n{sep}')
print('  EXACT ALGEBRAIC ANALYSIS of Cartan-Weyl graph')
print(sep)

print(f'\n  Degrees: {degs_cw} -> [6,6,5,5,5,5,5,5]')
print(f'  d_mean = 42/8 = 21/4 = {21/4}')
print(f'  q = 21/4 - 1 = 17/4')
print(f'  2*sqrt(q) = 2*sqrt(17/4) = sqrt(17) = {np.sqrt(17):.8f}')

# Check if eigenvalues have nice forms
print(f'\n  Testing algebraic forms for eigenvalues:')
# The characteristic polynomial of A_cw
# If degrees are [6,6,5,5,5,5,5,5], the graph is almost 5-regular
# except 2 vertices of degree 6

# Let me check specific values:
for candidate in [5, -3, 2, -2, 1, -1, 3, 4, -4, 0]:
    for m in range(1, 8):
        # Check if lambda = candidate is an eigenvalue
        pass

# Instead, compute characteristic polynomial numerically
from numpy.polynomial import polynomial as P
coeffs = np.poly(eigs_cw)
print(f'  Characteristic polynomial coefficients (highest to lowest):')
print(f'  {np.round(coeffs, 2)}')

# Try to identify integer eigenvalues
print(f'\n  Checking if eigenvalues are algebraic:')
for i, e in enumerate(eigs_cw):
    # Check nearby integers and simple fractions
    close_int = round(e)
    if abs(e - close_int) < 1e-6:
        print(f'    eig[{i}] = {e:.8f} = {close_int} (exact integer)')
    else:
        # Try quadratic irrationals: (a + sqrt(b))/c
        found = False
        for a in range(-10, 11):
            for b in range(1, 50):
                for c in [1, 2, 3, 4]:
                    for sign in [1, -1]:
                        val = (a + sign * np.sqrt(b)) / c
                        if abs(e - val) < 1e-6:
                            sstr = '+' if sign > 0 else '-'
                            print(f'    eig[{i}] = {e:.8f} = ({a}{sstr}sqrt({b}))/{c}')
                            found = True
                            break
                    if found: break
                if found: break
            if found: break
        if not found:
            print(f'    eig[{i}] = {e:.8f} (no simple algebraic form found)')

# Compute R_CW with exact q
q_exact = 17/4
R_cw_exact = lam_nt_cw / (2 * np.sqrt(q_exact))
print(f'\n  R(CW) with q = 17/4: {R_cw_exact:.8f}')

# Key formula: alpha = (1/2) * R_CW / (4*pi*3)
pi = np.pi
alpha_cw = 0.5 * R_cw / (4 * pi * 3)
print(f'\n  Formula: alpha = (1/2) * R(CW) / (4*pi*3)')
print(f'  alpha = {alpha_cw:.8f}')
print(f'  1/alpha = {1/alpha_cw:.4f}')
print(f'  Error vs 137.036: {abs(1/alpha_cw - 137.036)/137.036*100:.3f}%')

# ==============================================================
# COMPARE ALL (Cartan-Weyl vs Physics) 
# ==============================================================

print(f'\n{sep}')
print('  COMPARISON: Both graphs for SU(3)')
print(sep)

print(f'\n  {"Property":<25s} {"Cartan-Weyl":>15s} {"Physics":>15s}')
print(f'  {"-"*25} {"-"*15} {"-"*15}')
print(f'  {"Vertices":<25s} {dim:>15d} {8:>15d}')
print(f'  {"Edges":<25s} {edges_cw:>15d} {edges_phys:>15d}')
print(f'  {"Degrees":<25s} {"[6,6,5x6]":>15s} {"[6x3,7x4,4]":>15s}')
print(f'  {"d_mean":<25s} {degs_cw.mean():>15.4f} {degs_phys.mean():>15.4f}')
print(f'  {"q = d_mean - 1":<25s} {q_cw:>15.4f} {q_phys:>15.4f}')
print(f'  {"lambda_PF":<25s} {eigs_cw[0]:>15.6f} {eigs_phys[0]:>15.6f}')
print(f'  {"max|lambda_nt|":<25s} {lam_nt_cw:>15.6f} {lam_nt_phys:>15.6f}')
print(f'  {"R":<25s} {R_cw:>15.6f} {R_phys:>15.6f}')
print(f'  {"alpha = R/(24*pi)":<25s} {R_cw/(24*pi):>15.8f} {R_phys/(24*pi):>15.8f}')
print(f'  {"1/alpha":<25s} {24*pi/R_cw:>15.2f} {24*pi/R_phys:>15.2f}')
print(f'  {"Error vs 137.036":<25s} {abs(24*pi/R_cw-137.036)/137.036*100:>14.3f}% {abs(24*pi/R_phys-137.036)/137.036*100:>14.3f}%')

# ==============================================================
# IS THE CW EIGENVALUE EXACT? CHECK VIA SYMMETRY
# ==============================================================

print(f'\n{sep}')
print('  SYMMETRY ANALYSIS of Cartan-Weyl graph')
print(sep)

# The CW graph for A_2: vertices = {H1, H2, E1, E2, E3, E4, E5, E6}
# Labelss for printout
labels = ['H1', 'H2'] + [f'E({r})' for r in roots]
print('\n  Adjacency matrix:')
for i in range(8):
    row = '  '.join([str(int(A_cw[i, j])) for j in range(8)])
    print(f'    {labels[i]:>12s}: {row}')

# Check graph automorphisms: S3 acts on roots by permuting e1,e2,e3
# Also Z2 acts by -alpha
# The H-vertices are a clique of isolated vertices connected to all E
# Since H1, H2 connect to ALL 6 E-vertices, they're interchangeable

# Compute Laplacian
D_cw = np.diag(degs_cw.astype(float))
L_cw = D_cw - A_cw
eigs_L = np.sort(np.linalg.eigvalsh(L_cw))
print(f'\n  Laplacian eigenvalues: {np.round(eigs_L, 6)}')

# For CW graph: check if it's the complete bipartite K_{2,6} + extra E-E edges
# K_{2,6} would have edges_HE = 12, which we have
# Plus E-E edges = 9
# Total = 21 ✓

# The E-E subgraph: 6 vertices, 9 edges
A_EE = A_cw[2:, 2:]
eigs_EE = np.sort(np.linalg.eigvalsh(A_EE))[::-1]
degs_EE = A_EE.sum(axis=1).astype(int)
print(f'\n  E-E subgraph:')
print(f'    Vertices: 6, Edges: {int(A_EE.sum()/2)}')
print(f'    Degrees: {degs_EE}')
print(f'    Eigenvalues: {np.round(eigs_EE, 6)}')

# Check: is E-E subgraph = 3*K_2 + 3*P_3?
# Each root connects to its negative (3 pairs) and 2 neighbors via addition
# This is actually the OCTAHEDRON graph minus 3 edges? 
# No - octahedron has 6 vertices, 12 edges. We have 9 edges.

# The complement of E-E: 6 vertices, 15-9 = 6 edges
A_EE_comp = 1 - A_EE - np.eye(6)
print(f'\n  E-E complement: {int(A_EE_comp.sum()/2)} edges')
print(f'    Complement degrees: {A_EE_comp.sum(axis=1).astype(int)}')

# ==============================================================
# KEY: Can we derive alpha exactly from the CW spectrum?
# ==============================================================

print(f'\n{sep}')
print('  EXACT FORMULA DERIVATION (Cartan-Weyl)')
print(sep)

# We need the exact value of max|lambda_nt| for CW graph
# From eigenvalues: check if it's a root of a polynomial with small coefficients

# The eigenvalues are approximately:
# [5.276, 2.276, 0, 0, -1, -1, -2.276, -3.276]
# Pattern: looks like {5+x, 2+x, 0, 0, -1, -1, -2-x, -3-x} with x ≈ 0.276

# Test: is x = (sqrt(17)-4)/2?
x_test = (np.sqrt(17) - 4) / 2
print(f'  Test x = (sqrt(17)-4)/2 = {x_test:.8f}')
print(f'  eig[0] = {eigs_cw[0]:.8f}, 5+x = {5+x_test:.8f}')
print(f'  eig[1] = {eigs_cw[1]:.8f}, 2+x = {2+x_test:.8f}')

# Try: eigenvalues = (a + sqrt(b)) / c
# eig[0] ≈ 5.276: try (10+sqrt(b))/2
for b in range(1, 200):
    val = (10 + np.sqrt(b)) / 2
    if abs(val - eigs_cw[0]) < 1e-4:
        print(f'\n  eig[0] = (10+sqrt({b}))/2 = {val:.8f} ✓')
        break

# Actually, eig[0]+eig[-1] ≈ 5.276 + (-3.276) = 2.0
# eig[1]+eig[-2] ≈ 2.276 + (-2.276) = 0.0
# eig[2]+eig[3] = 0
# eig[4]+eig[5] = -2
# Sum of all eigs = tr(A) = 0 ✓ (2+0-2+0 = 0)

print(f'\n  Pairing analysis:')
print(f'    eig[0]+eig[7] = {eigs_cw[0]+eigs_cw[7]:.8f}') 
print(f'    eig[1]+eig[6] = {eigs_cw[1]+eigs_cw[6]:.8f}')
print(f'    eig[2]+eig[5] = {eigs_cw[2]+eigs_cw[5]:.8f}')
print(f'    eig[3]+eig[4] = {eigs_cw[3]+eigs_cw[4]:.8f}')

# Compute det(A_cw)
det_cw = np.linalg.det(A_cw)
print(f'\n  det(A_cw) = {det_cw:.4f}')
print(f'  Product of eigenvalues = {np.prod(eigs_cw):.4f}')

# trace of A^k
for k in range(1, 6):
    print(f'  tr(A^{k}) = {np.trace(np.linalg.matrix_power(A_cw, k)):.4f}')

# ==============================================================
# SEARCH FOR EXACT FORMULA COMBINING BOTH GRAPHS
# ==============================================================

print(f'\n{sep}')
print('  SYSTEMATIC FORMULA SEARCH using both graphs')
print(sep)

alpha_exp = 1/137.036
pi = np.pi

# Available exact quantities:
quantities = {
    'R_CW': R_cw,
    'R_phys': R_phys, 
    'lam_nt_CW': lam_nt_cw,
    'lam_nt_phys': lam_nt_phys,
    'q_CW': q_cw,
    'q_phys': q_phys,
    'edges_CW': 21,
    'edges_phys': 25,
    'dim': 8,
    'rank': 2,
    'hv': 3,
    'h': 3,
    '1/2': 0.5,  # R(SU(2)) exactly
    'pi': pi,
    'sqrt_17': np.sqrt(17),
    'sqrt_21': np.sqrt(21),
}

# Search: a * b / (c * d) for all combinations
best = []
names = list(quantities.keys())
vals = [quantities[n] for n in names]
N = len(names)

for i1 in range(N):
    for i2 in range(N):
        for i3 in range(N):
            for i4 in range(N):
                if len(set([i1,i2,i3,i4])) < 2: continue
                num = vals[i1] * vals[i2]
                den = vals[i3] * vals[i4]
                if abs(den) < 1e-10: continue
                cand = num / den
                if 0 < cand < 0.1:
                    err = abs(cand - alpha_exp) / alpha_exp * 100
                    if err < 0.5:
                        best.append((err, names[i1], names[i2], names[i3], names[i4], cand))

best.sort()
if best:
    print(f'\n  Formulas within 0.5% of 1/137.036:')
    seen = set()
    for err, a, b, c, d, val in best[:30]:
        key = tuple(sorted([a,b])) + tuple(sorted([c,d]))
        if key in seen: continue
        seen.add(key)
        print(f'    {a}*{b} / ({c}*{d}) = {val:.7f} -> 1/a = {1/val:.2f} (err {err:.3f}%)')
else:
    print('  No exact formulas found within 0.5%!')

# ==============================================================
# THE KEY QUESTION: Is R_CW/sqrt(17) rational?
# ==============================================================

print(f'\n{sep}')
print('  ALGEBRAIC FORM: Is max|lambda_nt(CW)| a quadratic irrational?')
print(sep)

# We know q_CW = 17/4, so R = lam_nt / sqrt(17)
ratio = lam_nt_cw / np.sqrt(17)
print(f'  max|lam_nt| / sqrt(17) = {ratio:.10f}')
print(f'  Test: 0.55201... = ?')

# Try p/q for small p, q
for p in range(1, 50):
    for q in range(1, 50):
        if abs(p/q - lam_nt_cw) < 1e-5:
            print(f'  max|lam_nt| ≈ {p}/{q} = {p/q:.10f}')

# Try (p + sqrt(r)) / q
for p in range(-20, 20):
    for r in range(1, 100):
        for q in [1, 2, 3, 4, 6, 8, 12]:
            for s in [1, -1]:
                val = (p + s*np.sqrt(r)) / q
                if abs(val - lam_nt_cw) < 1e-6:
                    sign = '+' if s > 0 else '-'
                    print(f'  max|lam_nt| = ({p}{sign}sqrt({r}))/{q} = {val:.10f}')

# Minimal polynomial approach
# If lam^4 + a*lam^3 + b*lam^2 + c*lam + d = 0, find a,b,c,d
print(f'\n  Minimal polynomial of max|lam_nt|:')
# Try degrees 1-4
lam = lam_nt_cw
for deg in [1, 2, 3, 4]:
    # Solve: sum_{k=0}^{deg} c_k * lam^k = 0 with c_deg = 1
    powers = [lam**k for k in range(deg)]
    target = -lam**deg
    if deg == 1:
        c = [target]
    else:
        c = np.linalg.lstsq(np.array(powers).reshape(1, -1).T.reshape(-1, deg), 
                             np.array([target]), rcond=None)[0] if deg > 1 else None
        # Simpler approach
        A_mat = np.zeros((1, deg))
        for k in range(deg):
            A_mat[0, k] = lam**k
        # Need more equations - use multiple eigenvalues
    # Use companion matrix approach: compute x^deg residual
    residual = lam**deg
    for k in range(deg):
        # Need to solve for coefficients
        pass

# Brute force: check if lam satisfies p(x)=0 for small integer polynomial
print(f'  lam = {lam:.10f}')
for a in range(-10, 11):
    for b in range(-20, 21):
        for c in range(-20, 21):
            val = lam**2 + a*lam + b + c/lam if abs(lam) > 1e-10 else float('inf')
            if abs(lam**2 + a*lam + b) < 0.01:
                # Check more precisely
                residual = lam**2 + a*lam + b
                if abs(residual) < 1e-5:
                    print(f'  lam^2 + {a}*lam + {b} = {residual:.2e} -> lam = ({-a}+-sqrt({a**2-4*b}))/2')

# Check cubic
for a in range(-10, 11):
    for b in range(-20, 21):
        for c in range(-50, 51):
            residual = lam**3 + a*lam**2 + b*lam + c
            if abs(residual) < 1e-4:
                print(f'  lam^3 + {a}*lam^2 + {b}*lam + {c} = {residual:.2e}')
