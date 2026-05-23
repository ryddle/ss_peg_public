"""
SPECTRAL ANALYSIS OF THE FULL SU(2)⊕SU(3) CARTAN-WEYL GRAPH
=============================================================

Goal: Does the full CW graph (11 vertices) encode the mass spectrum
building blocks {ln 2, ln R3, ln 3} directly in its spectrum?

Two versions:
  A. Disconnected: su(2) CW + su(3) CW (two components)
  B. Connected: add bifundamental edges between Cartan subalgebras
"""
import numpy as np
from scipy import sparse

# ============================================================
# 1. BUILD SU(2) CW GRAPH (3 vertices)
# ============================================================

# su(2): rank=1, dim=3
# Vertices: H1 (Cartan), E_α, E_-α (roots)
# Edges: [H1, E_α]≠0, [H1, E_-α]≠0
# E-E: α+(-α)=0, so E_α and E_-α are connected
# Total: 3 vertices, 2 edges → K3 (complete graph)

A_su2 = np.array([
    [0, 1, 1],  # H1
    [1, 0, 1],  # E_α
    [1, 1, 0],  # E_-α
])

eigs_su2 = np.sort(np.linalg.eigvalsh(A_su2))[::-1]
degs_su2 = A_su2.sum(axis=1)

print("=" * 70)
print("  SU(2) CW GRAPH (3 vertices, 2 edges)")
print("=" * 70)
print(f"  Vertices: H1, E_α, E_-α")
print(f"  Degree sequence: {degs_su2.astype(int)}")
print(f"  Eigenvalues: {np.round(eigs_su2, 8)}")
print(f"  d_max = {degs_su2.max()}")
print(f"  d_mean = {degs_su2.mean():.4f}")
print(f"  R = λ1/d_max = {eigs_su2[0]/degs_su2.max():.6f}")
print(f"  Expected: R = 1/2 = 0.5")
print()

# ============================================================
# 2. BUILD SU(3) CW GRAPH (8 vertices)
# ============================================================

# su(3): rank=2, dim=8
# Vertices: H1, H2 (Cartan), E_α1, E_α2, E_α3, E_-α1, E_-α2, E_-α3 (roots)
# Simple roots: α1, α2, α3 = α1+α2
# Root system A2: 6 roots at 60° intervals

# Adjacency: H_i connected to E_α if α(H_i) ≠ 0 (i.e., α is a root)
# E_α connected to E_β if α+β is a root or β = -α

A_su3 = np.zeros((8, 8))

# Labeling:
# 0,1: H1, H2 (Cartan)
# 2: E_α1 (simple root 1)
# 3: E_α2 (simple root 2)
# 4: E_α3 (composite root α1+α2)
# 5: E_-α1
# 6: E_-α2
# 7: E_-α3

# H-E edges: every root pairs with every Cartan
for i in range(2):       # H_i
    for j in range(2, 8):  # E_α
        A_su3[i, j] = A_su3[j, i] = 1

# E-E edges: roots connected if α+β is a root or β=-α
# α1 + α2 = α3 → E_α1-E_α2 connected
# α1 + (-α1) = 0 → E_α1-E_-α1 connected
# α2 + (-α2) = 0 → E_α2-E_-α2 connected
# α3 + (-α3) = 0 → E_α3-E_-α3 connected
# α1 + α3 = 2α1+α2 (not a root) → no edge
# α2 + α3 = α1+2α2 (not a root) → no edge
# -α1 + -α2 = -(α1+α2) = -α3 → E_-α1-E_-α2 connected
# -α1 + -α3 = -(2α1+α2) (not a root) → no edge
# -α2 + -α3 = -(α1+2α2) (not a root) → no edge

ee_edges = [
    (2, 3),   # α1 + α2 = α3
    (2, 5),   # α1 + (-α1) = 0
    (3, 6),   # α2 + (-α2) = 0
    (4, 7),   # α3 + (-α3) = 0
    (5, 6),   # -α1 + -α2 = -α3
]

for i, j in ee_edges:
    A_su3[i, j] = A_su3[j, i] = 1

eigs_su3 = np.sort(np.linalg.eigvalsh(A_su3))[::-1]
degs_su3 = A_su3.sum(axis=1)

print("=" * 70)
print("  SU(3) CW GRAPH (8 vertices, 21 edges)")
print("=" * 70)
print(f"  Vertices: H1, H2, E_α1, E_α2, E_α3, E_-α1, E_-α2, E_-α3")
print(f"  Degree sequence: {degs_su3.astype(int)}")
print(f"  Eigenvalues: {np.round(eigs_su3, 8)}")
print(f"  d_max = {degs_su3.max()}")
print(f"  d_mean = {degs_su3.mean():.4f}")
print(f"  q = d_mean - 1 = {degs_su3.mean() - 1:.4f}")
print(f"  2√q = {2*np.sqrt(degs_su3.mean() - 1):.6f}")
print(f"  max|λ_nt| = max(|{eigs_su3[1]:.6f}|, |{eigs_su3[-1]:.6f}|) = {max(abs(eigs_su3[1]), abs(eigs_su3[-1])):.6f}")

lam_nt = max(abs(eigs_su3[1]), abs(eigs_su3[-1]))
q = degs_su3.mean() - 1
R_su3_D1 = lam_nt / (2*np.sqrt(q))
print(f"  R (D1, mean-deg) = {R_su3_D1:.6f}")
print(f"  Expected: (√57-3)/(2√17) = {(np.sqrt(57)-3)/(2*np.sqrt(17)):.6f}")

# Characteristic polynomial factorization
print(f"\n  Characteristic polynomial factors:")
print(f"    λ^3 (from Cartan nullity)")
print(f"    (λ - 1) (from E-E subgraph)")
print(f"    (λ + 2)^2 (from E-E subgraph)")
print(f"    λ^2 - 3λ - 12 (extremal pair)")
print(f"    → roots: (3±√57)/2")

print()

# ============================================================
# 3. FULL DISCONNECTED GRAPH (11 vertices)
# ============================================================

print("=" * 70)
print("  FULL GRAPH: su(2) ⊕ su(3) — DISCONNECTED (11 vertices)")
print("=" * 70)

A_full_disc = np.block([
    [A_su2, np.zeros((3, 8))],
    [np.zeros((8, 3)), A_su3]
])

eigs_full_disc = np.sort(np.linalg.eigvalsh(A_full_disc))[::-1]
degs_full = A_full_disc.sum(axis=1)

print(f"\n  Full spectrum (11 eigenvalues):")
for i, e in enumerate(eigs_full_disc):
    component = "su(2)" if i < 3 else "su(3)"
    print(f"    λ_{i+1:2d} = {e:12.8f}  ({component})")

print(f"\n  Degree sequence: {degs_full.astype(int)}")
print(f"  d_max = {degs_full.max()}")
print(f"  d_mean = {degs_full.mean():.4f}")
print(f"  q = {degs_full.mean() - 1:.4f}")

# Ramanujan ratio for disconnected graph
lam_nt_full = max(abs(eigs_full_disc[1]), abs(eigs_full_disc[-1]))
q_full = degs_full.mean() - 1
R_full_D1 = lam_nt_full / (2*np.sqrt(q_full))
print(f"\n  R (D1, mean-deg) = {R_full_D1:.6f}")

# Individual R values
R_su2 = eigs_su2[0] / degs_su2.max()  # = 1/2
print(f"  R_su2 = {R_su2:.6f} (exact: 0.5)")
print(f"  R_su3 = {R_su3_D1:.6f}")
print(f"  R_su2 × R_su3 = {R_su2 * R_su3_D1:.6f}")

# ============================================================
# 4. CONNECTED GRAPH: add bifundamental edges
# ============================================================

print("\n" + "=" * 70)
print("  FULL GRAPH: su(2) ⊕ su(3) — CONNECTED (11 vertices)")
print("  Adding bifundamental edges between Cartan subalgebras")
print("=" * 70)

# The bifundamental (2,3) of SU(2)×SU(3) means:
# - The SU(2) Cartan H1 acts on the SU(3) fundamental (3 states)
# - The SU(3) Cartan H1, H2 act on the SU(2) fundamental (2 states)
# Add edges between Cartan generators representing this shared action

# Add 1 edge between H1(su2) and H1(su3)
# Add 1 edge between H1(su2) and H2(su3)
# This represents the bifundamental structure

A_connected = A_full_disc.copy()
# Edge between H1(su2, vertex 0) and H1(su3, vertex 3)
A_connected[0, 3] = A_connected[3, 0] = 1
# Edge between H1(su2, vertex 0) and H2(su3, vertex 4)
A_connected[0, 4] = A_connected[4, 0] = 1

eigs_connected = np.sort(np.linalg.eigvalsh(A_connected))[::-1]
degs_connected = A_connected.sum(axis=1)

print(f"\n  New degree sequence:")
for i, d in enumerate(degs_connected):
    if i < 3:
        name = ["H1(su2)", "E_α(su2)", "E_-α(su2)"][i]
    else:
        name = ["H1(su3)", "H2(su3)", "E_α1", "E_α2", "E_α3", "E_-α1", "E_-α2", "E_-α3"][i-3]
    print(f"    {name:>12s}: {int(d)}")

print(f"\n  Full spectrum (11 eigenvalues):")
for i, e in enumerate(eigs_connected):
    print(f"    λ_{i+1:2d} = {e:12.8f}")

print(f"\n  d_max = {degs_connected.max()}")
print(f"  d_mean = {degs_connected.mean():.4f}")
print(f"  q = {degs_connected.mean() - 1:.4f}")

lam_nt_conn = max(abs(eigs_connected[1]), abs(eigs_connected[-1]))
q_conn = degs_connected.mean() - 1
R_conn_D1 = lam_nt_conn / (2*np.sqrt(q_conn))
print(f"\n  R (D1, mean-deg) = {R_conn_D1:.6f}")

# Characteristic polynomial
char_coeffs = np.poly(eigs_connected)
print(f"\n  Characteristic polynomial coefficients:")
for i, c in enumerate(char_coeffs):
    power = len(char_coeffs) - 1 - i
    if abs(c) > 1e-10:
        sign = "+" if c > 0 else "-"
        print(f"    {sign} {abs(c):10.6f}·λ^{power}")

# ============================================================
# 5. ANALYZE: do eigenvalues encode {ln 2, ln R3, ln 3}?
# ============================================================

print("\n" + "=" * 70)
print("  ANALYSIS: Do eigenvalues encode the mass building blocks?")
print("=" * 70)

# The building blocks in log-space:
ln2 = np.log(2)
lnR3 = np.log((np.sqrt(57) - 3) / (2 * np.sqrt(17)))
ln3 = np.log(3)

print(f"\n  Target basis (log-space):")
print(f"    ln 2  = {ln2:.8f}")
print(f"    ln R3 = {lnR3:.8f}")
print(f"    ln 3  = {ln3:.8f}")

print(f"\n  Eigenvalues of full graph (disconnected):")
for i, e in enumerate(eigs_full_disc):
    print(f"    λ_{i+1:2d} = {e:.8f}")

print(f"\n  Eigenvalues of connected graph:")
for i, e in enumerate(eigs_connected):
    print(f"    λ_{i+1:2d} = {e:.8f}")

# Check: do any eigenvalues match the building blocks directly?
print(f"\n  Direct matches (eigenvalue ≈ building block):")
all_eigs = np.concatenate([eigs_full_disc, eigs_connected])
targets = [ln2, lnR3, ln3, -ln2, -lnR3, -ln3, np.log(R_su2), np.log(R_su3_D1), np.log(R_full_D1)]
target_names = ["ln 2", "ln R3", "ln 3", "-ln 2", "-ln R3", "-ln 3", "ln R_su2", "ln R_su3", "ln R_full"]

for t, name in zip(targets, target_names):
    matches = np.abs(all_eigs - t) < 0.01
    if np.any(matches):
        print(f"    {name:>12s} = {t:.8f}: MATCH at λ_{np.argmin(np.abs(all_eigs - t)) + 1}")
    else:
        best_idx = np.argmin(np.abs(all_eigs - t))
        print(f"    {name:>12s} = {t:.8f}: closest = λ_{best_idx+1} = {all_eigs[best_idx]:.8f} (diff = {abs(all_eigs[best_idx] - t):.6f})")

# ============================================================
# 6. TRY: spectral ratios and combinations
# ============================================================

print("\n" + "=" * 70)
print("  ANALYSIS: Spectral ratios and combinations")
print("=" * 70)

# Check ratios of eigenvalues
print(f"\n  Ratios of eigenvalues (full disconnected graph):")
for i in range(len(eigs_full_disc)):
    for j in range(i+1, len(eigs_full_disc)):
        ratio = eigs_full_disc[i] / eigs_full_disc[j]
        for target, name in [(ln2, "ln 2"), (lnR3, "ln R3"), (ln3, "ln 3"),
                             (np.log(R_su2), "ln R_su2"), (np.log(R_su3_D1), "ln R_su3"),
                             (1, "1"), (2, "2"), (3, "3"), (1/2, "1/2"), (1/3, "1/3")]:
            if abs(ratio - target) < 0.01:
                print(f"    λ_{i+1}/λ_{j+1} = {ratio:.6f} ≈ {name} = {target:.6f}")

# Check: does the spectral gap encode anything?
print(f"\n  Spectral gaps:")
print(f"    su(2): λ1 - λ2 = {eigs_su2[0] - eigs_su2[1]:.6f}")
print(f"    su(3): λ1 - λ2 = {eigs_su3[0] - eigs_su3[1]:.6f}")
print(f"    full:  λ1 - λ2 = {eigs_full_disc[0] - eigs_full_disc[1]:.6f}")

# Check: does the nullity encode anything?
nullity_full = np.sum(np.abs(eigs_full_disc) < 1e-10)
print(f"\n  Nullity (zero eigenvalues): {int(nullity_full)}")
print(f"  This equals rank(su(2)) + rank(su(3)) = 1 + 2 = 3")

# ============================================================
# 7. TRY: different Ramanujan definitions on the full graph
# ============================================================

print("\n" + "=" * 70)
print("  ANALYSIS: Ramanujan ratios for the full graph")
print("=" * 70)

d_max = degs_full.max()
d_mean = degs_full.mean()
d_geom = np.exp(np.mean(np.log(degs_full[degs_full > 0])))

print(f"\n  Degree statistics:")
print(f"    d_max = {d_max}")
print(f"    d_mean = {d_mean:.4f}")
print(f"    d_geom = {d_geom:.4f}")

# D1: mean degree
q_D1 = d_mean - 1
R_D1 = lam_nt_full / (2*np.sqrt(q_D1))
print(f"\n  D1 (mean degree):")
print(f"    q = {q_D1:.4f}")
print(f"    2√q = {2*np.sqrt(q_D1):.6f}")
print(f"    R = {R_D1:.6f}")

# D2: max degree
q_D2 = d_max - 1
R_D2 = lam_nt_full / (2*np.sqrt(q_D2))
print(f"\n  D2 (max degree):")
print(f"    q = {q_D2}")
print(f"    R = {R_D2:.6f}")

# D3: geometric mean degree
q_D3 = d_geom - 1
R_D3 = lam_nt_full / (2*np.sqrt(q_D3))
print(f"\n  D3 (geometric mean):")
print(f"    q = {q_D3:.4f}")
print(f"    R = {R_D3:.6f}")

# Check: which R is closest to the individual R2 and R3?
print(f"\n  Comparison:")
print(f"    R_su2 = {R_su2:.6f}")
print(f"    R_su3 = {R_su3_D1:.6f}")
print(f"    R_full_D1 = {R_D1:.6f}")
print(f"    R_full_D2 = {R_D2:.6f}")
print(f"    R_full_D3 = {R_D3:.6f}")
print(f"    R_EE = {1/np.sqrt(2):.6f}")

# ============================================================
# 8. TRY: spectral invariants of the connected graph
# ============================================================

print("\n" + "=" * 70)
print("  ANALYSIS: Spectral invariants of CONNECTED graph")
print("=" * 70)

d_max_c = degs_connected.max()
d_mean_c = degs_connected.mean()
d_geom_c = np.exp(np.mean(np.log(degs_connected[degs_connected > 0])))

print(f"\n  Degree statistics (connected):")
print(f"    d_max = {d_max_c}")
print(f"    d_mean = {d_mean_c:.4f}")
print(f"    d_geom = {d_geom_c:.4f}")

lam_nt_c = max(abs(eigs_connected[1]), abs(eigs_connected[-1]))
print(f"\n  Nontrivial eigenvalues:")
print(f"    λ2 = {eigs_connected[1]:.8f}")
print(f"    λ11 = {eigs_connected[-1]:.8f}")
print(f"    max|λ_nt| = {lam_nt_c:.8f}")

q_D1_c = d_mean_c - 1
R_D1_c = lam_nt_c / (2*np.sqrt(q_D1_c))
print(f"\n  D1 (mean degree):")
print(f"    q = {q_D1_c:.4f}")
print(f"    R = {R_D1_c:.6f}")

# Check: does the connected graph have a Ramanujan ratio closer to R_su3?
print(f"\n  Comparison:")
print(f"    R_su3 = {R_su3_D1:.6f}")
print(f"    R_full_connected_D1 = {R_D1_c:.6f}")
print(f"    Difference: {abs(R_D1_c - R_su3_D1):.6f}")

# ============================================================
# 9. TRY: eigenvalue patterns and the 3D basis
# ============================================================

print("\n" + "=" * 70)
print("  ANALYSIS: Eigenvalue patterns vs. 3D basis")
print("=" * 70)

print(f"\n  The 3D basis is {ln2:.6f}, {lnR3:.6f}, {ln3:.6f}")
print(f"  These are logarithms of: 2, R3={R_su3_D1:.6f}, 3")
print(f"  The eigenvalues are NOT logarithms — they are adjacency eigenvalues.")
print(f"  So the question is: can we extract {ln2}, {lnR3}, {ln3} from the spectrum?")

# Try: log of eigenvalue ratios
print(f"\n  Log of eigenvalue ratios:")
for i in range(len(eigs_full_disc)):
    for j in range(i+1, len(eigs_full_disc)):
        log_ratio = np.log(abs(eigs_full_disc[i] / eigs_full_disc[j]))
        for target, name in [(ln2, "ln 2"), (lnR3, "ln R3"), (ln3, "ln 3")]:
            if abs(log_ratio - target) < 0.05:
                print(f"    log|λ_{i+1}/λ_{j+1}| = {log_ratio:.6f} ≈ {name} = {target:.6f}")

# Try: eigenvalue differences
print(f"\n  Eigenvalue differences:")
for i in range(len(eigs_full_disc)):
    for j in range(i+1, len(eigs_full_disc)):
        diff = abs(eigs_full_disc[i] - eigs_full_disc[j])
        for target, name in [(ln2, "ln 2"), (lnR3, "ln R3"), (ln3, "ln 3"), (1, "1"), (2, "2"), (3, "3")]:
            if abs(diff - target) < 0.05:
                print(f"    |λ_{i+1} - λ_{j+1}| = {diff:.6f} ≈ {name} = {target:.6f}")

# ============================================================
# 10. SUMMARY
# ============================================================

print("\n" + "=" * 70)
print("  SUMMARY")
print("=" * 70)

print(f"""
  Key findings:

  1. The disconnected graph spectrum is just the UNION of the two
     component spectra. No new information appears.

  2. The individual R values (R_su2 = 1/2, R_su3 = {R_su3_D1:.6f})
     are computed FROM the component graphs, not from the full graph.

  3. The 3D basis {ln2:.6f}, {lnR3:.6f}, {ln3:.6f} does NOT appear
     directly as eigenvalues or simple combinations of eigenvalues.

  4. Adding bifundamental edges changes the spectrum but does not
     obviously produce the building blocks.

  CONCLUSION: The building blocks (ln 2, ln R3, ln 3) are NOT
  spectral invariants of the full CW graph in a direct sense.
  They are derived quantities from the COMPONENT graphs.

  This means the mass spectrum is encoded in the SPECTRAL PROPERTIES
  of the individual factors (su(2) and su(3)), not in the spectrum
  of their direct sum. The graph topology matters, but at the level
  of the individual gauge group factors, not their combination.

  NEXT STEPS:
  - Try: are the building blocks related to the CHARACTERISTIC
    POLYNOMIAL coefficients of the component graphs?
  - Try: are they related to the GRAPH LAPLACIAN spectrum (not
    adjacency)?
  - Try: is the 3D basis related to the EIGENVECTOR structure
    (not just eigenvalues)?
""")
