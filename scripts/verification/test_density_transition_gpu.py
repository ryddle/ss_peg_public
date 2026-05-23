"""
Test C.3: Density transition — algebra closure vs generator density.

Fix dim=d and sweep n_gen from 1 to d²-1 (SU) or d(d-1)/2 (SO),
measuring what fraction of commutator pairs close into the span.

This tests the central prediction: closure is a sharp transition
at n_gen = dim(algebra), not a gradual convergence.

Usage:
    python test_density_transition_gpu.py                    # default: SU dim=4
    python test_density_transition_gpu.py --group su --dim 3
    python test_density_transition_gpu.py --group so --dim 5
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils_gpu import (
    params_to_generators,
    generators_to_params,
    random_params,
    params_to_antisymmetric,
    antisymmetric_to_params,
    random_params_so,
    evolve_generators_generic,
    killing_metric_np,
    check_in_span_np,
    commutator_np,
    get_device_info,
)


def measure_closure(generators):
    """Measure the fraction of commutator pairs that close in the span."""
    n_gen = len(generators)
    if n_gen < 2:
        return 1.0, 0, 0
    closes = 0
    total = n_gen * (n_gen - 1) // 2
    for i in range(n_gen):
        for j in range(i + 1, n_gen):
            C = commutator_np(generators[i], generators[j])
            if np.linalg.norm(C) < 1e-10:
                closes += 1
            elif check_in_span_np(C, generators, tol=1e-3):
                closes += 1
    return closes / total if total > 0 else 1.0, closes, total


def run_density_transition(group_type="su", dim=4, n_seeds=3):
    """Sweep n_gen at fixed dim and measure closure fraction."""

    if group_type == "su":
        n_gen_full = dim**2 - 1
        label_fn = lambda d: f"SU({d})"
        param_to_gen = params_to_generators
        gen_to_param = generators_to_params
        rand_param = lambda bs, ng, device, seeds: random_params(bs, ng, dim, device, seeds)
        gen_kwargs = {"dim": dim}
    elif group_type == "so":
        n_gen_full = dim * (dim - 1) // 2
        label_fn = lambda d: f"SO({d})"
        param_to_gen = params_to_antisymmetric
        gen_to_param = antisymmetric_to_params
        rand_param = lambda bs, ng, device, seeds: random_params_so(bs, ng, dim, device, seeds)
        gen_kwargs = {"dim": dim}
    else:
        raise ValueError(f"Unknown group type: {group_type}")

    # Sweep points: sparse below transition, dense near and above
    sweep_points = sorted(set(
        list(range(1, min(4, n_gen_full + 1))) +
        list(range(max(1, n_gen_full - 3), n_gen_full + 1)) +
        [n_gen_full // 2, n_gen_full // 4, 3 * n_gen_full // 4]
    ))
    sweep_points = [p for p in sweep_points if 1 <= p <= n_gen_full]

    print("=" * 78)
    print(f"  DENSITY TRANSITION — {label_fn(dim)} (dim={dim})")
    print("=" * 78)
    print(f"  Device: {get_device_info()}")
    print(f"  Full algebra: n_gen = {n_gen_full}")
    print(f"  Sweep: {len(sweep_points)} points, n_gen = {sweep_points}")
    print(f"  Seeds per point: {n_seeds}")

    results = []
    t_total = time.time()

    for n_gen in sweep_points:
        steps = max(15000, 3000 * dim)
        lr = max(0.0001, 0.001 / dim)

        print(f"\n  n_gen = {n_gen}/{n_gen_full} ...", end=" ", flush=True)

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
            verbose=False,
            convergence_threshold=1e-4,
            param_to_gen_fn=param_to_gen,
            gen_to_param_fn=gen_to_param,
            random_param_fn=rand_param,
            param_to_gen_kwargs=gen_kwargs,
        )
        gpu_time = time.time() - t0

        seed_closures = []
        seed_K_diag = []
        for seed_idx, (gens_np, hist) in enumerate(raw):
            gen_list = list(gens_np)
            frac, closes, total = measure_closure(gen_list)
            seed_closures.append(frac)

            K = killing_metric_np(gen_list)
            K_diag_mean = np.mean(np.diag(K))
            seed_K_diag.append(K_diag_mean)

        avg_closure = np.mean(seed_closures)
        avg_K = np.mean(seed_K_diag)
        tag = "FULL" if avg_closure > 0.99 else "PARTIAL" if avg_closure > 0.5 else "SPARSE"

        print(f"closure={avg_closure:.1%} K={avg_K:.4f} [{tag}] ({gpu_time:.1f}s)")

        results.append({
            "n_gen": n_gen,
            "n_gen_full": n_gen_full,
            "density": n_gen / n_gen_full,
            "avg_closure": avg_closure,
            "avg_K_diag": avg_K,
            "gpu_time": gpu_time,
            "tag": tag,
        })

    t_total = time.time() - t_total

    # Summary table
    print(f"\n\n{'=' * 78}")
    print(f"  DENSITY TRANSITION TABLE — {label_fn(dim)} (dim={dim})")
    print(f"{'=' * 78}")

    print(f"\n  {'n_gen':>6s} {'density':>9s} {'closure':>9s} {'K_diag':>10s} {'status':>9s}")
    print("  " + "─" * 50)
    for r in results:
        print(f"  {r['n_gen']:>6d} {r['density']:>8.1%} {r['avg_closure']:>8.1%} "
              f"{r['avg_K_diag']:>10.4f} {r['tag']:>9s}")
    print("  " + "─" * 50)

    # Identify transition point (skip n_gen=1 which is trivially closed)
    for i, r in enumerate(results):
        if r["n_gen"] > 1 and r["avg_closure"] > 0.99:
            print(f"\n  Transition: closure reaches 100% at n_gen={r['n_gen']} "
                  f"(density={r['density']:.1%} of full algebra)")
            if r['n_gen'] == r['n_gen_full']:
                print(f"  => SHARP TRANSITION: algebra requires ALL {r['n_gen_full']} generators")
            break
    else:
        print(f"\n  No full closure achieved in sweep range.")

    print(f"\n  Total wall clock: {t_total:.1f}s")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--group", default="su", choices=["su", "so"])
    parser.add_argument("--dim", type=int, default=4)
    parser.add_argument("--seeds", type=int, default=3)
    args = parser.parse_args()

    run_density_transition(group_type=args.group, dim=args.dim, n_seeds=args.seeds)
