"""
Test A.4: Sp(2n) Lie algebra emergence from Killing metric optimization.

Tests:
  - Sp(2) = USp(2): 3 generators in dim=2 (≅ SU(2), cross-check)
  - Sp(4) = USp(4): 10 generators in dim=4 (≅ SO(5), cross-check)
  - Sp(6) = USp(6): 21 generators in dim=6

Usage:
    python test_sp_n_gpu.py                     # run Sp(2,4,6) sweep
    python test_sp_n_gpu.py --groups 1          # Sp(2) only (half_dim=1)
    python test_sp_n_gpu.py --groups 1,2        # Sp(2) and Sp(4)
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils_gpu import (
    params_to_symplectic,
    symplectic_to_params,
    random_params_sp,
    evolve_generators_generic,
    killing_metric_np,
    structure_constants_real_np,
    compute_casimir_np,
    check_in_span_np,
    commutator_np,
    get_device_info,
    _build_omega,
    DTYPE,
)


def n_generators_sp(half_dim):
    """Number of generators for Sp(2n) = USp(2n): n(2n+1)."""
    n = half_dim
    return n * (2 * n + 1)


def analyze_sp_result(generators, half_dim, seed_idx, verbose=True):
    """Full analysis of optimized Sp(2n) generators."""
    n = half_dim
    dim = 2 * n
    n_gen = len(generators)
    expected = n_generators_sp(n)

    # Hermiticity: T† = T
    herm_err = max(np.linalg.norm(g - g.conj().T) for g in generators)

    # Symplectic constraint: X^T·Ω + Ω·X = 0 where X = -iT
    Omega = np.zeros((dim, dim), dtype=complex)
    Omega[:n, n:] = np.eye(n)
    Omega[n:, :n] = -np.eye(n)
    sp_errs = []
    for g in generators:
        X = -1j * g
        err = np.abs(X.T @ Omega + Omega @ X).max()
        sp_errs.append(err)
    sp_err = max(sp_errs)

    # Killing metric
    K = killing_metric_np(generators)
    K_diag = np.diag(K)
    K_diag_mean = np.mean(K_diag)
    K_diag_std = np.std(K_diag)
    K_offdiag = np.abs(K - np.diag(K_diag)).max()

    # Closure
    closes_count = 0
    total_pairs = n_gen * (n_gen - 1) // 2
    for i in range(n_gen):
        for j in range(i + 1, n_gen):
            C = commutator_np(generators[i], generators[j])
            if np.linalg.norm(C) < 1e-10:
                closes_count += 1
            elif check_in_span_np(C, generators, tol=1e-4):
                closes_count += 1
    closure_frac = closes_count / total_pairs if total_pairs > 0 else 1.0

    # Structure constants (hermitian convention: [T_a,T_b] = i f_{abc} T_c)
    f = structure_constants_real_np(generators)
    f_rms = np.sqrt(np.mean(f**2))

    # Jacobi (sample)
    jacobi_max = 0.0
    lim = min(n_gen, 6)
    for a in range(lim):
        for b in range(a + 1, lim):
            for c in range(b + 1, lim):
                j_val = sum(
                    f[a, b, d] * f[d, c, :].sum() +
                    f[b, c, d] * f[d, a, :].sum() +
                    f[c, a, d] * f[d, b, :].sum()
                    for d in range(n_gen)
                )
                jacobi_max = max(jacobi_max, abs(j_val))

    # Casimir
    C_mat = compute_casimir_np(generators)
    C_diag = np.diag(C_mat).real
    C_lambda = np.mean(C_diag)
    C_spread = np.std(C_diag)

    # kappa
    kappa = np.mean([np.trace(g @ g).real for g in generators])

    # |K|/κ² — for Sp(2n), h∨ = n+1
    kappa_abs = abs(kappa)
    k_ratio = abs(K_diag_mean) / kappa_abs**2 if kappa_abs > 1e-12 else 0.0
    h_dual_theoretical = n + 1

    result = {
        "half_dim": n,
        "dim": dim,
        "n_gen": n_gen,
        "n_gen_expected": expected,
        "seed": seed_idx,
        "herm_err": herm_err,
        "sp_err": sp_err,
        "K_diag_mean": K_diag_mean,
        "K_diag_std": K_diag_std,
        "K_offdiag": K_offdiag,
        "closure_frac": closure_frac,
        "closes_count": closes_count,
        "total_pairs": total_pairs,
        "f_rms": f_rms,
        "jacobi_max": jacobi_max,
        "C_lambda": C_lambda,
        "C_spread": C_spread,
        "kappa": kappa,
        "k_ratio": k_ratio,
        "h_dual_theoretical": h_dual_theoretical,
    }

    if verbose:
        tag = "PASS" if closure_frac > 0.99 else "PARTIAL" if closure_frac > 0.5 else "FAIL"
        print(f"\n  seed={seed_idx} [{tag}]:")
        print(f"    Hermiticity err:  {herm_err:.2e}")
        print(f"    Symplectic err:   {sp_err:.2e}")
        print(f"    K_diag:           {K_diag_mean:.6f} +/- {K_diag_std:.6f}")
        print(f"    K_offdiag max:    {K_offdiag:.2e}")
        print(f"    Closure:          {closes_count}/{total_pairs} ({closure_frac:.1%})")
        print(f"    f_rms:            {f_rms:.6f}")
        print(f"    Jacobi max:       {jacobi_max:.2e}")
        print(f"    C_lambda:         {C_lambda:.6f} (spread: {C_spread:.2e})")
        print(f"    kappa=Tr(T^2):    {kappa:.6f}")
        print(f"    |K|/kappa^2:      {k_ratio:.4f} (h_dual_theo: {h_dual_theoretical})")

    return result


def run_test_sp(half_dim, n_seeds=3, steps=None, lr=None):
    """Run Sp(2n) optimizer and analyze."""
    n = half_dim
    dim = 2 * n
    n_gen = n_generators_sp(n)

    if steps is None:
        steps = max(20000, 5000 * dim)
    if lr is None:
        lr = max(0.0001, 0.001 / dim)

    print(f"\n{'=' * 78}")
    print(f"  Sp({dim}) = USp({dim}): {n_gen} generators in dim={dim}")
    print(f"  steps={steps}, lr={lr}, seeds={n_seeds}")
    print(f"{'=' * 78}")

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
        param_to_gen_fn=params_to_symplectic,
        gen_to_param_fn=lambda gens: symplectic_to_params(gens, half_dim=n),
        random_param_fn=lambda bs, ng, device, seeds: random_params_sp(bs, ng, n, device, seeds),
        param_to_gen_kwargs={"half_dim": n},
    )
    gpu_time = time.time() - t0

    results = []
    for seed_idx, (gens_np, hist) in enumerate(raw):
        gen_list = list(gens_np)
        r = analyze_sp_result(gen_list, n, seed_idx)
        r["gpu_time"] = gpu_time / n_seeds
        results.append(r)

    any_pass = any(r["closure_frac"] > 0.99 for r in results)
    avg_closure = np.mean([r["closure_frac"] for r in results])
    avg_K = np.mean([r["K_diag_mean"] for r in results])

    print(f"\n  SUMMARY Sp({dim}):")
    print(f"    GPU time:          {gpu_time:.1f}s ({gpu_time/n_seeds:.1f}s/seed)")
    print(f"    Avg closure:       {avg_closure:.1%}")
    print(f"    Avg K_diag:        {avg_K:.6f}")
    print(f"    Verdict:           {'PASS' if any_pass else 'FAIL'}")

    return {
        "group": f"Sp({dim})",
        "half_dim": n,
        "dim": dim,
        "n_gen": n_gen,
        "results": results,
        "gpu_time": gpu_time,
        "verdict": "PASS" if any_pass else "FAIL",
    }


def run_test_sp_sweep(half_dims=None, n_seeds=3):
    """Run Sp(2n) sweep."""
    if half_dims is None:
        half_dims = [1, 2, 3]

    print("=" * 78)
    print("  Sp(2n) SWEEP — Killing Metric Optimization")
    print("=" * 78)
    print(f"  Device: {get_device_info()}")
    print(f"  Groups: {', '.join(f'Sp({2*n})' for n in half_dims)}")
    print(f"  Seeds per group: {n_seeds}")

    all_results = []
    t_total = time.time()

    for n in half_dims:
        r = run_test_sp(n, n_seeds=n_seeds)
        all_results.append(r)

    t_total = time.time() - t_total

    # Summary
    print(f"\n\n{'=' * 78}")
    print(f"  COMPARATIVE SUMMARY — Sp(2n) sweep")
    print(f"{'=' * 78}")
    print(f"\n  {'Group':8s} {'n_gen':>6s} {'dim':>5s} {'Closure':>9s} {'K_diag':>10s} "
          f"{'|K|/k^2':>10s} {'h_theo':>10s} {'C_lam':>10s} {'Time':>8s} {'Verdict':>8s}")
    print("  " + "─" * 90)
    for r in all_results:
        avg = lambda key: np.mean([s[key] for s in r["results"]])
        print(f"  {r['group']:8s} {r['n_gen']:>6d} {r['dim']:>5d} {avg('closure_frac'):>8.1%} "
              f"{avg('K_diag_mean'):>10.6f} {avg('k_ratio'):>10.4f} "
              f"{int(avg('h_dual_theoretical')):>10d} {avg('C_lambda'):>10.6f} "
              f"{r['gpu_time']:>7.1f}s {r['verdict']:>8s}")
    print("  " + "─" * 90)
    print(f"\n  Total wall clock: {t_total:.1f}s")

    # Cross-checks
    sp2 = next((r for r in all_results if r["half_dim"] == 1), None)
    if sp2:
        avg_cl = np.mean([s["closure_frac"] for s in sp2["results"]])
        print(f"\n  Cross-check Sp(2) ~ SU(2) (3 gen, dim=2):")
        print(f"    Sp(2) closure = {avg_cl:.1%}")
        print(f"    SU(2) closure = 100% (from prior runs)")

    sp4 = next((r for r in all_results if r["half_dim"] == 2), None)
    if sp4:
        avg_cl = np.mean([s["closure_frac"] for s in sp4["results"]])
        avg_K = np.mean([s["K_diag_mean"] for s in sp4["results"]])
        print(f"\n  Cross-check Sp(4) ~ SO(5) (10 gen, different dim):")
        print(f"    Sp(4) closure = {avg_cl:.1%} (dim=4)")
        print(f"    SO(5) closure = 100% (dim=5, from prior runs)")

    return all_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--groups", default="1,2,3",
                        help="Comma-separated half-dims: Sp(2n) for n=...")
    parser.add_argument("--seeds", type=int, default=3)
    args = parser.parse_args()

    half_dims = [int(x) for x in args.groups.split(",")]
    run_test_sp_sweep(half_dims=half_dims, n_seeds=args.seeds)
