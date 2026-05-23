"""
Quick diagnostic: does su(2) ⊂ su(3) crystallise with evolve_generators_generic?
Tests n_gen=3, dim=3 with different killing_target values.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from shared_utils_gpu import (
    evolve_generators_generic,
    killing_metric_np,
    check_in_span_np,
    commutator_np,
)


def classify(gens):
    generators = list(gens)
    n = len(generators)
    K = killing_metric_np(generators)
    K_diag = float(np.mean(np.diag(K)))
    K_offdiag = float(np.mean(np.abs(K - np.diag(np.diag(K)))))

    total_pairs = n * (n - 1) // 2
    closes = 0
    comm_norms = []
    for i in range(n):
        for j in range(i + 1, n):
            C = commutator_np(generators[i], generators[j])
            cn = np.linalg.norm(C)
            comm_norms.append(cn)
            if cn < 1e-6:
                closes += 1
            elif check_in_span_np(C, generators, tol=0.1):
                closes += 1
    closure = closes / total_pairs if total_pairs > 0 else 0.0
    return K_diag, K_offdiag, closure, comm_norms


print("=" * 60)
print("  DIAGNOSTIC: su(2) crystallisation in SU(3) param")
print("=" * 60)

for kt in [-1.0, -2.0, -4.0, -8.0]:
    print(f"\n  killing_target = {kt}")
    results = evolve_generators_generic(
        n_gen=3, dim=3,
        steps=20_000, lr=0.001,
        alpha=5.0, beta=1.0,
        killing_target=kt,
        batch_size=4,
        seeds=[1, 2, 3, 4],
        verbose=False,
    )
    for i, (gens_np, hist) in enumerate(results):
        kd, ko, cl, cn = classify(gens_np)
        print(
            f"    seed={i+1}: K_diag={kd:+.4f}, K_off={ko:.4f}, "
            f"closure={cl:.2f}, ||comm||={[f'{c:.4f}' for c in cn]}, "
            f"T_total={hist['T_total'][-1]:.6f}"
        )
