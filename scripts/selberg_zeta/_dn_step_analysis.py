"""Analyze the D_n step anomaly visible in spectral evolution plots."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from ramanujan_extended import roots_A, roots_B, roots_C, roots_D
from ramanujan_extended import adjacency_from_roots, spectral_analysis

print("=== D_n family detailed ===")
print(f"{'n':>3} {'h':>4} {'dim':>5} {'deg':>8} {'max_nt':>8} {'ratio':>7} {'gap':>7} {'ram':>5}")
print("-" * 60)
dn_data = {}
for n in range(3, 10):
    h = 2 * (n - 1)
    rank, roots, simple = roots_D(n)
    adj = adjacency_from_roots(rank, roots, simple)
    sa = spectral_analysis(adj)
    deg = sa['degree'] if sa['is_regular'] else np.mean(np.sum(adj, axis=1))
    dn_data[n] = sa
    print(f"{n:3d} {h:4d} {adj.shape[0]:5d} {deg:8.2f} {sa['max_nontrivial']:8.3f} "
          f"{sa['ratio']:7.4f} {sa['spectral_gap']:7.3f} {bool(sa['is_ramanujan']):>5}")

print()
print("=== A_n for comparison (A3 ~ D3) ===")
print(f"{'n':>3} {'h':>4} {'dim':>5} {'deg':>8} {'max_nt':>8} {'ratio':>7} {'gap':>7} {'ram':>5}")
print("-" * 60)
for n in [2, 3, 4, 5]:
    h = n + 1
    rank, roots, simple = roots_A(n)
    adj = adjacency_from_roots(rank, roots, simple)
    sa = spectral_analysis(adj)
    deg = sa['degree'] if sa['is_regular'] else np.mean(np.sum(adj, axis=1))
    print(f"{n:3d} {h:4d} {adj.shape[0]:5d} {deg:8.2f} {sa['max_nontrivial']:8.3f} "
          f"{sa['ratio']:7.4f} {sa['spectral_gap']:7.3f} {bool(sa['is_ramanujan']):>5}")

print()
print("=== D3 vs A3 spectra comparison ===")
rank_d, roots_d3, sim_d = roots_D(3)
adj_d = adjacency_from_roots(rank_d, roots_d3, sim_d)
eigs_d = sorted(np.linalg.eigvalsh(adj_d.astype(float)))

rank_a, roots_a3, sim_a = roots_A(3)
adj_a = adjacency_from_roots(rank_a, roots_a3, sim_a)
eigs_a = sorted(np.linalg.eigvalsh(adj_a.astype(float)))

print(f"D3 eigenvalues: {[round(x, 3) for x in eigs_d]}")
print(f"A3 eigenvalues: {[round(x, 3) for x in eigs_a]}")
print(f"Spectra identical: {np.allclose(eigs_d, eigs_a)}")

print()
print("=== D4 triality check: eigenvalue multiplicities ===")
rank4, roots4, sim4 = roots_D(4)
adj4 = adjacency_from_roots(rank4, roots4, sim4)
eigs4 = sorted(np.linalg.eigvalsh(adj4.astype(float)))
unique4, counts4 = np.unique(np.round(eigs4, 6), return_counts=True)
print(f"D4 = so(8), dim={adj4.shape[0]}, unique eigenvalues={len(unique4)}")
print("Eigenvalue multiplicities (triality => divisible by 3):")
n_div3 = 0
for v, c in zip(unique4, counts4):
    div3 = " <-- div by 3" if c % 3 == 0 else ""
    if c % 3 == 0:
        n_div3 += 1
    print(f"  {v:8.3f}  mult={c:2d}{div3}")
print(f"Fraction with mult divisible by 3: {n_div3}/{len(unique4)}")

print()
print("=== D5 check: eigenvalue multiplicities ===")
rank5, roots5, sim5 = roots_D(5)
adj5 = adjacency_from_roots(rank5, roots5, sim5)
eigs5 = sorted(np.linalg.eigvalsh(adj5.astype(float)))
unique5, counts5 = np.unique(np.round(eigs5, 6), return_counts=True)
print(f"D5 = so(10), dim={adj5.shape[0]}, unique eigenvalues={len(unique5)}")
for v, c in zip(unique5, counts5):
    div3 = " <-- div by 3" if c % 3 == 0 else ""
    print(f"  {v:8.3f}  mult={c:2d}{div3}")

# Quantify the step: compute delta(ratio) between consecutive ranks
print()
print("=== Step quantification: delta(ratio) between consecutive ranks ===")
print(f"{'Family':>8} {'n':>3}->{'n+1':>3} {'h':>4}->{'h+1':>4} {'drho':>8}")
print("-" * 45)
# D_n
for n in range(3, 9):
    h1 = 2 * (n - 1)
    h2 = 2 * n
    r1 = dn_data[n]['ratio']
    r2 = dn_data[n + 1]['ratio']
    delta = r2 - r1
    flag = " <<<" if abs(delta) > 0.10 else ""
    print(f"{'Dn':>8} {n:3d}->{n+1:3d} {h1:4d}->{h2:4d} {delta:8.4f}{flag}")

# Compare with B_n and C_n steps
print()
print("=== B_n steps for comparison ===")
for n in range(2, 9):
    h = 2 * n
    rank, roots, simple = roots_B(n)
    adj = adjacency_from_roots(rank, roots, simple)
    sa1 = spectral_analysis(adj)
    h2 = 2 * (n + 1)
    rank2, roots2, simple2 = roots_B(n + 1)
    adj2 = adjacency_from_roots(rank2, roots2, simple2)
    sa2 = spectral_analysis(adj2)
    delta = sa2['ratio'] - sa1['ratio']
    flag = " <<<" if abs(delta) > 0.10 else ""
    print(f"{'Bn':>8} {n:3d}->{n+1:3d} {h:4d}->{h2:4d} {delta:8.4f}{flag}")

# Degree sequence analysis for D3, D4, D5
print()
print("=== Degree sequences D3, D4, D5 ===")
for n in [3, 4, 5]:
    rank, roots, simple = roots_D(n)
    adj = adjacency_from_roots(rank, roots, simple)
    degs = np.sum(adj, axis=1)
    unique_d, counts_d = np.unique(degs, return_counts=True)
    cv = np.std(degs) / np.mean(degs) if np.mean(degs) > 0 else 0
    print(f"D{n}: dim={len(degs)}, mean_deg={np.mean(degs):.2f}, CV={cv:.4f}")
    for d, c in zip(unique_d, counts_d):
        print(f"  degree {d:3d}: count={c:2d}")
