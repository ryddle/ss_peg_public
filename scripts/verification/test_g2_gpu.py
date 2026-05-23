"""
Test B.1: G₂ exceptional Lie algebra — research attempt.

G₂ = aut(O) is the smallest exceptional Lie algebra:
  - 14 generators, rank 2
  - Fundamental representation: dim=7
  - h∨ = 4, I_fund = 1 → expected |K|/κ² = 4

Two approaches:
  1. BLIND: SO(7) parameterization with n_gen=14 — optimizer must discover g₂
  2. CONSTRUCTIVE: Build g₂ generators from octonion 3-form, verify with framework

Approach 1 fails (0% closure) — g₂ occupies a 14D submanifold of 21D SO(7)
that gradient descent cannot find blindly.

Approach 2 succeeds — constructing g₂ from φ_{ijk} and verifying emergence.

Usage:
    python test_g2_gpu.py                  # constructive + brief blind test
    python test_g2_gpu.py --blind-only     # only the blind optimization attempt
    python test_g2_gpu.py --seeds 5        # more seeds for blind test
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils_gpu import (
    params_to_antisymmetric,
    antisymmetric_to_params,
    random_params_so,
    evolve_generators_generic,
    killing_metric_np,
    structure_constants_general_np,
    check_in_span_np,
    commutator_np,
    get_device_info,
    DTYPE,
)


def analyze_g2_result(generators, seed_idx, verbose=True):
    """Full analysis of a G₂ candidate."""
    n_gen = len(generators)
    dim = generators[0].shape[0]

    # Antisymmetry: T^T = -T
    antisym_err = max(np.linalg.norm(g + g.T) for g in generators)
    # Reality: Im(T) = 0
    reality_err = max(np.max(np.abs(g.imag)) for g in generators)

    # Killing metric
    K = killing_metric_np(generators)
    K_diag = np.diag(K).real
    K_diag_mean = np.mean(K_diag)
    K_diag_std = np.std(K_diag)
    K_offdiag = np.abs(K - np.diag(K_diag)).max()

    # Killing eigenvalues
    K_eigs = np.linalg.eigvalsh(K.real)

    # Closure
    closes_count = 0
    total_pairs = n_gen * (n_gen - 1) // 2
    max_err = 0.0
    for i in range(n_gen):
        for j in range(i + 1, n_gen):
            C = commutator_np(generators[i], generators[j])
            norm_C = np.linalg.norm(C)
            if norm_C < 1e-10:
                closes_count += 1
            elif check_in_span_np(C, generators, tol=1e-4):
                closes_count += 1
            else:
                max_err = max(max_err, norm_C)
    closure_frac = closes_count / total_pairs if total_pairs > 0 else 1.0

    # Generator norm
    kappa = np.mean([np.trace(g @ g).real for g in generators])
    kappa_abs = abs(kappa)
    k_ratio = abs(K_diag_mean) / kappa_abs**2 if kappa_abs > 1e-12 else 0.0

    # G₂ invariants: h∨ = 4, I_fund = 1
    h_dual_theo = 4
    I_R = 1.0
    h_dual_measured = I_R * k_ratio

    # Rank estimation: count near-zero eigenvalues of ad(h₁), ad(h₂)
    # (For a quick check, count how many generators have [T_i, T_j] ≈ 0
    #  for all other generators — these form the Cartan subalgebra)
    diag_generators = []
    for i in range(n_gen):
        n_commuting = sum(
            np.linalg.norm(commutator_np(generators[i], generators[j])) < 1e-4
            for j in range(n_gen) if j != i
        )
        diag_generators.append((i, n_commuting))
    max_commuting = max(c for _, c in diag_generators) if diag_generators else 0

    result = {
        "seed": seed_idx,
        "n_gen": n_gen,
        "dim": dim,
        "antisym_err": antisym_err,
        "reality_err": reality_err,
        "K_diag_mean": K_diag_mean,
        "K_diag_std": K_diag_std,
        "K_offdiag": K_offdiag,
        "K_eigs": K_eigs,
        "closure_frac": closure_frac,
        "closes_count": closes_count,
        "total_pairs": total_pairs,
        "kappa": kappa,
        "k_ratio": k_ratio,
        "h_dual_measured": h_dual_measured,
        "h_dual_theo": h_dual_theo,
        "max_commuting": max_commuting,
    }

    if verbose:
        tag = "PASS" if closure_frac > 0.99 else "PARTIAL" if closure_frac > 0.5 else "FAIL"
        print(f"\n  seed={seed_idx} [{tag}]:")
        print(f"    Antisymmetry err: {antisym_err:.2e}")
        print(f"    Reality err:      {reality_err:.2e}")
        print(f"    K_diag:           {K_diag_mean:.6f} +/- {K_diag_std:.6f}")
        print(f"    K_offdiag max:    {K_offdiag:.2e}")
        print(f"    K eigenvalues:    min={K_eigs[0]:.4f} max={K_eigs[-1]:.4f}")
        print(f"    Closure:          {closes_count}/{total_pairs} ({closure_frac:.1%})")
        print(f"    kappa=Tr(T^2):    {kappa:.6f}")
        print(f"    |K|/kappa^2:      {k_ratio:.4f}")
        print(f"    h_dual measured:  {h_dual_measured:.4f} (theo: {h_dual_theo})")
        print(f"    Max commuting:    {max_commuting} (rank indicator)")

    return result


def run_g2_test(n_seeds=3, steps=None, lr=None):
    """Run G₂ optimization: 14 generators in SO(7) parameterization."""
    n_gen = 14
    dim = 7

    if steps is None:
        steps = 40000
    if lr is None:
        lr = 0.0001

    print("=" * 78)
    print("  G₂ EXCEPTIONAL ALGEBRA TEST")
    print("=" * 78)
    print(f"  Device:     {get_device_info()}")
    print(f"  Approach:   14 generators in SO(7) parameterization")
    print(f"  Expected:   g₂ ⊂ so(7), h∨ = 4, |K|/κ² = 4")
    print(f"  Steps:      {steps}, LR: {lr}, Seeds: {n_seeds}")
    print()

    t0 = time.time()
    raw = evolve_generators_generic(
        n_gen=n_gen,
        dim=dim,
        steps=steps,
        lr=lr,
        alpha=5.0,
        beta=1.0,
        killing_target=-1.0,
        batch_size=n_seeds,
        seeds=list(range(n_seeds)),
        verbose=True,
        convergence_threshold=1e-4,
        param_to_gen_fn=params_to_antisymmetric,
        gen_to_param_fn=antisymmetric_to_params,
        random_param_fn=lambda bs, ng, device, seeds: random_params_so(bs, ng, dim, device, seeds),
        param_to_gen_kwargs={"dim": dim},
    )
    gpu_time = time.time() - t0

    results = []
    for seed_idx, (gens_np, hist) in enumerate(raw):
        gen_list = list(gens_np)
        r = analyze_g2_result(gen_list, seed_idx)
        r["gpu_time"] = gpu_time / n_seeds
        results.append(r)

    any_pass = any(r["closure_frac"] > 0.99 for r in results)
    avg_closure = np.mean([r["closure_frac"] for r in results])

    print(f"\n  SUMMARY G₂:")
    print(f"    GPU time:     {gpu_time:.1f}s ({gpu_time/n_seeds:.1f}s/seed)")
    print(f"    Avg closure:  {avg_closure:.1%}")
    print(f"    Any PASS:     {any_pass}")
    if any_pass:
        best = max(results, key=lambda r: r["closure_frac"])
        print(f"    Best seed:    {best['seed']}")
        print(f"    h∨ measured:  {best['h_dual_measured']:.4f} (expected: 4)")
        print(f"    |K|/κ²:       {best['k_ratio']:.4f} (expected: 4)")
    else:
        print("    Note: G₂ is a challenging optimization target.")
        print("    Consider: more steps, different lr, or larger batch_size.")

    return results


# ===================================================================
# CONSTRUCTIVE APPROACH: Build G₂ from octonion 3-form
# ===================================================================
def build_g2_generators():
    """Construct the 14 generators of g₂ ⊂ so(7) from octonion structure.

    G₂ = {A ∈ so(7) : A preserves the octonion cross product}
    Constraint: Σ_a A_{ia}φ_{ajk} + A_{ja}φ_{iak} + A_{ka}φ_{ija} = 0

    Returns: list of 14 orthonormalized 7×7 real antisymmetric matrices.
    """
    # Fano plane triples: φ_{ijk} = +1 for these (0-indexed)
    # Convention: e_a · e_b = φ_{abc} e_c for imaginary octonions
    fano_triples = [
        (0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6),
        (4, 5, 0), (5, 6, 1), (6, 0, 2)
    ]

    phi = np.zeros((7, 7, 7))
    for (a, b, c) in fano_triples:
        phi[a, b, c] = 1
        phi[b, c, a] = 1
        phi[c, a, b] = 1
        phi[a, c, b] = -1
        phi[c, b, a] = -1
        phi[b, a, c] = -1

    # so(7) basis: E_{ij} for i < j
    basis = []
    for i in range(7):
        for j in range(i + 1, 7):
            E = np.zeros((7, 7))
            E[i, j] = 1
            E[j, i] = -1
            basis.append(E)
    n_basis = len(basis)  # 21

    # Constraint matrix: M · coeffs = 0 gives g₂
    ijk_list = [(i, j, k) for i in range(7) for j in range(i+1, 7) for k in range(j+1, 7)]
    n_constraints = len(ijk_list)  # C(7,3) = 35

    M = np.zeros((n_constraints, n_basis))
    for row, (i, j, k) in enumerate(ijk_list):
        for col, E in enumerate(basis):
            val = 0.0
            for a in range(7):
                val += E[i, a] * phi[a, j, k]
                val += E[j, a] * phi[i, a, k]
                val += E[k, a] * phi[i, j, a]
            M[row, col] = val

    # Null space of M = g₂ coefficients in so(7) basis
    _, S, Vt = np.linalg.svd(M)
    rank = np.sum(S > 1e-10)
    null_dim = n_basis - rank

    if null_dim != 14:
        print(f"  WARNING: null space dimension = {null_dim} (expected 14)")

    g2_coefficients = Vt[rank:]  # (14, 21) — null space rows
    generators = []
    for coeffs in g2_coefficients:
        G = sum(c * E for c, E in zip(coeffs, basis))
        generators.append(G)

    # Orthonormalize w.r.t. Frobenius inner product
    ortho = []
    for g in generators:
        v = g.copy()
        for u in ortho:
            proj = np.sum(v * u) / np.sum(u * u)
            v = v - proj * u
        norm = np.sqrt(np.sum(v * v))
        if norm > 1e-12:
            ortho.append(v / norm)
    generators = ortho

    print(f"  Constraint matrix: rank={rank}, null space={null_dim}")
    print(f"  Constructed {len(generators)} orthonormal g₂ generators")
    return generators


def run_constructive_g2():
    """Build g₂ generators and analyze them."""
    print("=" * 78)
    print("  G₂ CONSTRUCTIVE TEST — From Octonion 3-Form")
    print("=" * 78)

    generators = build_g2_generators()
    results = analyze_g2_result(generators, seed_idx="constructive")

    # Additional: verify antisymmetry and trace
    for i, g in enumerate(generators):
        assert np.allclose(g, -g.T), f"Generator {i} not antisymmetric!"
        assert np.allclose(g.imag, 0), f"Generator {i} not real!"

    # Verify: are these REALLY in g₂? Check φ-preservation for a few triples
    fano_triples = [
        (0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6),
        (4, 5, 0), (5, 6, 1), (6, 0, 2)
    ]
    phi = np.zeros((7, 7, 7))
    for (a, b, c) in fano_triples:
        phi[a, b, c] = 1; phi[b, c, a] = 1; phi[c, a, b] = 1
        phi[a, c, b] = -1; phi[c, b, a] = -1; phi[b, a, c] = -1

    max_phi_err = 0.0
    for g in generators:
        for i in range(7):
            for j in range(i+1, 7):
                for k in range(j+1, 7):
                    val = sum(g[i, a] * phi[a, j, k] +
                              g[j, a] * phi[i, a, k] +
                              g[k, a] * phi[i, j, a] for a in range(7))
                    max_phi_err = max(max_phi_err, abs(val))
    print(f"  φ-preservation error: {max_phi_err:.2e}")

    # Summary
    is_pass = results["closure_frac"] > 0.99
    print(f"\n  CONSTRUCTIVE G₂ VERDICT: {'PASS' if is_pass else 'FAIL'}")
    if is_pass:
        print(f"    14 generators, 100% closure, h∨ = {results['h_dual_measured']:.4f} (theo: 4)")
        print(f"    |K|/κ² = {results['k_ratio']:.4f}")
        print(f"    This IS the exceptional Lie algebra G₂!")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--steps", type=int, default=40000)
    parser.add_argument("--lr", type=float, default=0.0001)
    parser.add_argument("--blind-only", action="store_true")
    parser.add_argument("--no-blind", action="store_true")
    args = parser.parse_args()

    if not args.blind_only:
        run_constructive_g2()
        print()

    if not args.no_blind:
        run_g2_test(n_seeds=args.seeds, steps=args.steps, lr=args.lr)
