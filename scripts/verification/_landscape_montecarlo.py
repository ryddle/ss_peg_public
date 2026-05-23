"""
Monte Carlo Landscape: probability measure on the variational theory space.

Runs N independent direct+anneal optimizations for 52 generators in so(27)
and classifies each outcome into:
  - F4           : simple, closed, h_dual ~ 9
  - SIMPLE_OTHER : simple, closed, h_dual != 9
  - CLOSED_NON_SIMPLE : closed, K not proportional to I
  - PARTIAL      : 50-99% closure
  - NOT_CLOSED   : < 50% closure

Then measures:
  1. Frequency of each category  (landscape measure)
  2. Energy (T_total) per category (depth of minima)
  3. Stability of F4 solutions under perturbation (basin width)

Batch elements share a single Adam optimizer but have independent gradients
(each T_total[b] depends only on params[b]), so they evolve independently.
This allows GPU-parallel Monte Carlo.

Usage:
  python _landscape_montecarlo.py                          # 32 seeds, default
  python _landscape_montecarlo.py --n-seeds 64             # more statistics
  python _landscape_montecarlo.py --n-seeds 16 --quick     # fast test run
"""

import sys, os, time, argparse, json
import numpy as np
import torch

sys.path.insert(0, os.path.dirname(__file__))

from shared_utils_gpu import (
    batched_killing_metric,
    batched_non_commutativity_tension,
    batched_norm_tension,
    params_to_antisymmetric,
    antisymmetric_to_params,
    random_params_so,
    killing_metric_np,
    DTYPE, FDTYPE, DEVICE,
)
from test_closure_driven import batched_closure_tension

N_GEN = 52
DIM = 27
ALPHA = 5.0
BETA = 1.0


# ============================================================
# PHASE 1: OPTIMIZER
# ============================================================

def run_batch(seeds, s1_steps, s2_steps, s1_gamma, s2_gamma,
              s1_lr, s2_lr, verbose=True):
    """Run direct+anneal for a batch of independent seeds.

    Returns:
        all_gens: (B, N_GEN, DIM, DIM) numpy complex128
        tensions: dict with per-element (B,) numpy arrays
    """
    B = len(seeds)
    device = DEVICE

    params = random_params_so(B, N_GEN, DIM, device, seeds)

    # Normalize initial generators to ||T||_F = 1
    with torch.no_grad():
        gens = params_to_antisymmetric(params, N_GEN, dim=DIM)
        norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
        gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
        params = antisymmetric_to_params(gens)
    params = params.clone().detach().requires_grad_(True)

    K_target = -1.0 * torch.eye(N_GEN, dtype=FDTYPE, device=device)

    for stage, (steps, lr, gamma) in enumerate([
        (s1_steps, s1_lr, s1_gamma),
        (s2_steps, s2_lr, s2_gamma),
    ], 1):
        if verbose:
            print(f"  Stage {stage}: gamma={gamma}, lr={lr}, {steps} steps")
        optimizer = torch.optim.Adam([params], lr=lr)
        for step in range(steps):
            optimizer.zero_grad()
            generators = params_to_antisymmetric(params, N_GEN, dim=DIM)
            T_nc = batched_non_commutativity_tension(generators)
            K = batched_killing_metric(generators)
            T_k = ((K - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
            T_norm = batched_norm_tension(generators)
            T_cl = batched_closure_tension(generators)
            T_total = T_nc + ALPHA * T_k + BETA * T_norm + gamma * T_cl
            loss = T_total.sum()
            loss.backward()
            optimizer.step()

            if verbose and step % 2000 == 0:
                with torch.no_grad():
                    tcl_vals = T_cl.cpu().numpy()
                    print(f"    step {step:5d}: T_cl avg={tcl_vals.mean():.4e} "
                          f"min={tcl_vals.min():.4e} max={tcl_vals.max():.4e}",
                          flush=True)

    # Compute final per-element tensions
    with torch.no_grad():
        generators = params_to_antisymmetric(params, N_GEN, dim=DIM)
        T_nc = batched_non_commutativity_tension(generators).cpu().numpy()
        K = batched_killing_metric(generators)
        T_k = ((K - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1)).cpu().numpy()
        T_norm = batched_norm_tension(generators).cpu().numpy()
        T_cl = batched_closure_tension(generators).cpu().numpy()

        final_gamma = s2_gamma
        T_total = T_nc + ALPHA * T_k + BETA * T_norm + final_gamma * T_cl
        T_phys = T_nc + ALPHA * T_k + BETA * T_norm  # without closure penalty

        all_gens = generators.cpu().numpy()

    return all_gens, {
        'T_nc': T_nc, 'T_k': T_k, 'T_norm': T_norm,
        'T_cl': T_cl, 'T_total': T_total, 'T_phys': T_phys,
    }


# ============================================================
# PHASE 2: CLASSIFICATION
# ============================================================

def classify_solution(gens_np):
    """Classify a set of 52 generators in so(27).

    Uses vectorized lstsq (all 1326 pairs solved in one LAPACK call).

    Returns: dict with category and all diagnostic metrics.
    """
    gens_np = np.asarray(gens_np)
    n = len(gens_np)
    d = gens_np[0].shape[0]

    # Vectorized closure analysis
    idx_i, idx_j = np.triu_indices(n, k=1)
    G_i = gens_np[idx_i]  # (n_pairs, d, d)
    G_j = gens_np[idx_j]
    C = G_i @ G_j - G_j @ G_i  # (n_pairs, d, d)
    C_flat = C.reshape(-1, d * d).T  # (d^2, n_pairs)

    B_mat = np.array([g.flatten() for g in gens_np]).T  # (d^2, n)

    # Solve all 1326 least-squares problems in one call
    coeffs, _, _, _ = np.linalg.lstsq(B_mat, C_flat, rcond=None)
    proj = B_mat @ coeffs  # (d^2, n_pairs)
    residual_norms = np.linalg.norm(C_flat - proj, axis=0)  # (n_pairs,)

    closure_frac = float(np.mean(residual_norms < 1e-4))
    max_res = float(residual_norms.max())

    # Killing metric
    K = killing_metric_np(list(gens_np))
    K_diag = np.diag(K).real
    K_mean = float(K_diag.mean())
    K_std = float(K_diag.std())

    # Simplicity: ||K/mean - I|| / ||I||
    if abs(K_mean) > 1e-12:
        K_normed = K / K_mean
        simplicity = float(
            np.linalg.norm(K_normed - np.eye(n)) / np.linalg.norm(np.eye(n))
        )
    else:
        simplicity = 999.0

    # h_dual: normalize generators to ||T||_F = 1, then compute
    gens_r = [g / np.linalg.norm(g) for g in gens_np]
    K_r = killing_metric_np(gens_r)
    Tr_sq = np.mean([np.trace(g @ g).real for g in gens_r])
    h_raw = float(abs(np.diag(K_r).real.mean() / Tr_sq)) if abs(Tr_sq) > 1e-12 else 0.0
    h_dual_ir3 = 3.0 * h_raw  # h^v = I_R * h_raw, I_R = 3 for 27-dim rep

    # Classification
    is_closed = closure_frac >= 0.99
    is_simple = simplicity < 0.05
    is_f4 = is_closed and is_simple and 7.0 < h_dual_ir3 < 11.0

    if is_f4:
        category = "F4"
    elif is_closed and is_simple:
        category = "SIMPLE_OTHER"
    elif is_closed:
        category = "CLOSED_NON_SIMPLE"
    elif closure_frac >= 0.5:
        category = "PARTIAL"
    else:
        category = "NOT_CLOSED"

    return {
        'category': category,
        'closure_frac': closure_frac,
        'max_residual': max_res,
        'K_mean': K_mean,
        'K_std': K_std,
        'simplicity': simplicity,
        'h_raw': h_raw,
        'h_dual_ir3': h_dual_ir3,
    }


# ============================================================
# PHASE 3: STABILITY TEST
# ============================================================

def stability_test(f4_gens_np, epsilons, recovery_steps=5000,
                   gamma=200.0, lr=0.0001, verbose=True):
    """Perturb F4 generators and test recovery.

    f4_gens_np: (N_GEN, DIM, DIM) numpy complex128
    epsilons: list of noise scales (relative to param std)
    """
    device = DEVICE

    # Convert generators -> params
    gens_torch = torch.tensor(
        np.array(f4_gens_np), dtype=DTYPE, device=device
    ).unsqueeze(0)  # (1, N_GEN, DIM, DIM)
    params_f4 = antisymmetric_to_params(gens_torch)  # (1, N_GEN, n_params)
    param_std = params_f4.std().item()

    K_target = -1.0 * torch.eye(N_GEN, dtype=FDTYPE, device=device)

    results = []
    rng = torch.Generator(device=device).manual_seed(12345)

    for eps in epsilons:
        if verbose:
            print(f"\n  eps = {eps:.0e}: ", end="", flush=True)

        sigma = eps * param_std
        noise = sigma * torch.randn(
            params_f4.shape, dtype=FDTYPE, device=device, generator=rng
        )
        params = (params_f4 + noise).clone().detach().requires_grad_(True)

        optimizer = torch.optim.Adam([params], lr=lr)
        for step in range(recovery_steps):
            optimizer.zero_grad()
            generators = params_to_antisymmetric(params, N_GEN, dim=DIM)
            T_nc = batched_non_commutativity_tension(generators)
            K = batched_killing_metric(generators)
            T_k = ((K - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
            T_norm = batched_norm_tension(generators)
            T_cl = batched_closure_tension(generators)
            T_total = T_nc + ALPHA * T_k + BETA * T_norm + gamma * T_cl
            loss = T_total.sum()
            loss.backward()
            optimizer.step()

            if verbose and step % 1000 == 0:
                print(f"T_cl={T_cl.item():.2e} ", end="", flush=True)

        # Classify recovered solution
        with torch.no_grad():
            final_gens = params_to_antisymmetric(params, N_GEN, dim=DIM)
            final_np = final_gens[0].cpu().numpy()

        gens_list = [final_np[k] for k in range(N_GEN)]
        info = classify_solution(gens_list)
        info['epsilon'] = eps
        results.append(info)

        if verbose:
            print(f"\n    -> {info['category']}, "
                  f"closure={info['closure_frac']:.1%}, "
                  f"simplicity={info['simplicity']:.4f}, "
                  f"h_dual(IR=3)={info['h_dual_ir3']:.2f}")

    return results


# ============================================================
# REPORTING
# ============================================================

def print_report(all_results, actual_seeds, elapsed):
    """Print formatted Monte Carlo results."""

    print(f"\n{'=' * 80}")
    print("  RESULTS TABLE")
    print(f"{'=' * 80}")
    header = (f"{'Seed':>4s}  {'Category':20s}  {'Close%':>7s}  {'MaxRes':>9s}  "
              f"{'K_std':>9s}  {'Simple':>8s}  {'h_v(3)':>7s}  {'T_total':>8s}  {'T_phys':>8s}")
    print(header)
    print("-" * len(header))
    for r in all_results:
        print(f"{r['seed']:4d}  {r['category']:20s}  {r['closure_frac']:6.1%}  "
              f"{r['max_residual']:9.2e}  {r['K_std']:9.2e}  "
              f"{r['simplicity']:8.4f}  {r['h_dual_ir3']:7.2f}  "
              f"{r['T_total']:8.1f}  {r['T_phys']:8.1f}")

    # Category summary
    categories = {}
    for r in all_results:
        categories.setdefault(r['category'], []).append(r)

    print(f"\n{'=' * 80}")
    print("  CATEGORY FREQUENCIES (landscape measure)")
    print(f"{'=' * 80}")
    for cat in ['F4', 'SIMPLE_OTHER', 'CLOSED_NON_SIMPLE', 'PARTIAL', 'NOT_CLOSED']:
        if cat in categories:
            count = len(categories[cat])
            pct = count / actual_seeds * 100
            print(f"  {cat:20s}: {count:3d}/{actual_seeds} = {pct:5.1f}%")
        else:
            print(f"  {cat:20s}:   0/{actual_seeds} =   0.0%")

    # Energy comparison
    print(f"\n{'=' * 80}")
    print("  ENERGY COMPARISON")
    print(f"{'=' * 80}")
    print(f"  {'Category':20s}  {'T_total':>12s}  {'T_phys':>12s}  "
          f"{'T_cl':>12s}  {'T_k':>12s}")
    print(f"  {'-' * 20}  {'-' * 12}  {'-' * 12}  {'-' * 12}  {'-' * 12}")
    for cat in ['F4', 'SIMPLE_OTHER', 'CLOSED_NON_SIMPLE', 'PARTIAL', 'NOT_CLOSED']:
        if cat in categories:
            rs = categories[cat]
            tt = [r['T_total'] for r in rs]
            tp = [r['T_phys'] for r in rs]
            tc = [r['T_cl'] for r in rs]
            tk = [r['T_k'] for r in rs]
            n = len(rs)
            tt_str = f"{np.mean(tt):.2f}" + (f" +/- {np.std(tt):.2f}" if n > 1 else "")
            tp_str = f"{np.mean(tp):.2f}" + (f" +/- {np.std(tp):.2f}" if n > 1 else "")
            tc_str = f"{np.mean(tc):.2e}" + (f"" if n <= 1 else "")
            tk_str = f"{np.mean(tk):.2e}" + (f"" if n <= 1 else "")
            print(f"  {cat:20s}  {tt_str:>12s}  {tp_str:>12s}  "
                  f"{tc_str:>12s}  {tk_str:>12s}")

    print(f"\n  Total elapsed: {elapsed:.0f}s ({elapsed / 60:.1f} min)")

    return categories


def print_stability_report(stab_results):
    """Print formatted stability results."""
    print(f"\n  {'eps':>7s}  {'Category':15s}  {'Closure':>8s}  "
          f"{'Simple':>8s}  {'h_v(3)':>7s}  {'Recovered':>9s}")
    print(f"  {'-' * 7}  {'-' * 15}  {'-' * 8}  {'-' * 8}  {'-' * 7}  {'-' * 9}")
    for sr in stab_results:
        recovered = "YES" if sr['category'] == 'F4' else "no"
        print(f"  {sr['epsilon']:7.0e}  {sr['category']:15s}  "
              f"{sr['closure_frac']:7.1%}  {sr['simplicity']:8.4f}  "
              f"{sr['h_dual_ir3']:7.2f}  {recovered:>9s}")


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Monte Carlo landscape study for F4 in so(27)')
    parser.add_argument('--n-seeds', type=int, default=32,
                        help='Total number of independent seeds')
    parser.add_argument('--batch-size', type=int, default=16,
                        help='Seeds per GPU batch')
    parser.add_argument('--s1-steps', type=int, default=8000,
                        help='Stage 1 optimizer steps')
    parser.add_argument('--s2-steps', type=int, default=3000,
                        help='Stage 2 optimizer steps')
    parser.add_argument('--s1-gamma', type=float, default=50.0)
    parser.add_argument('--s2-gamma', type=float, default=200.0)
    parser.add_argument('--s1-lr', type=float, default=0.001)
    parser.add_argument('--s2-lr', type=float, default=0.0001)
    parser.add_argument('--stability-steps', type=int, default=5000,
                        help='Recovery steps per perturbation level')
    parser.add_argument('--no-stability', action='store_true',
                        help='Skip stability test')
    parser.add_argument('--save', type=str, default='landscape_results.npz',
                        help='Output file for results')
    parser.add_argument('--save-gens', action='store_true',
                        help='Also save generator tensors to npz')
    parser.add_argument('--seeds', type=int, nargs='+', default=None,
                        help='Explicit seed list (overrides --n-seeds)')
    parser.add_argument('--quick', action='store_true',
                        help='Quick test: 2k+1k steps')
    args = parser.parse_args()

    if args.quick:
        args.s1_steps = 2000
        args.s2_steps = 1000
        args.stability_steps = 1000

    # Seed list: explicit or auto-generated
    if args.seeds is not None:
        all_seeds = args.seeds
        batch_size = min(args.batch_size, len(all_seeds))
    else:
        n_seeds = args.n_seeds
        batch_size = args.batch_size
        n_batches_tmp = (n_seeds + batch_size - 1) // batch_size
        all_seeds = list(range(n_batches_tmp * batch_size))

    # Partition into batches
    seed_batches = [all_seeds[i:i+batch_size]
                    for i in range(0, len(all_seeds), batch_size)]
    n_batches = len(seed_batches)
    actual_seeds = len(all_seeds)

    print("=" * 80)
    print("  MONTE CARLO LANDSCAPE: F4 in so(27)")
    print("=" * 80)
    print(f"  Seeds:   {actual_seeds} ({n_batches} batch(es), bs<={batch_size})")
    if args.seeds is not None:
        print(f"  Explicit seeds: {args.seeds}")
    print(f"  Stage 1: {args.s1_steps} steps, gamma={args.s1_gamma}, lr={args.s1_lr}")
    print(f"  Stage 2: {args.s2_steps} steps, gamma={args.s2_gamma}, lr={args.s2_lr}")
    print(f"  Density: {N_GEN}/{DIM * (DIM - 1) // 2} = "
          f"{N_GEN / (DIM * (DIM - 1) // 2) * 100:.1f}%")
    print(f"  Device:  {DEVICE}")
    print()

    all_results = []
    all_gens_list = []
    t0 = time.time()

    for batch_idx in range(n_batches):
        seeds = seed_batches[batch_idx]
        print(f"--- Batch {batch_idx + 1}/{n_batches}: seeds {seeds} ---")

        t_batch = time.time()
        gens_batch, tensions = run_batch(
            seeds, args.s1_steps, args.s2_steps,
            args.s1_gamma, args.s2_gamma, args.s1_lr, args.s2_lr,
            verbose=True,
        )
        dt_opt = time.time() - t_batch
        print(f"  Optimization: {dt_opt:.0f}s ({dt_opt / 60:.1f} min)")

        # Classify each seed
        print("  Classifying...", flush=True)
        t_class = time.time()
        for b in range(len(seeds)):
            seed = seeds[b]
            gens_list = [gens_batch[b, k] for k in range(N_GEN)]
            info = classify_solution(gens_list)
            info['seed'] = seed
            info['T_nc'] = float(tensions['T_nc'][b])
            info['T_k'] = float(tensions['T_k'][b])
            info['T_norm'] = float(tensions['T_norm'][b])
            info['T_cl'] = float(tensions['T_cl'][b])
            info['T_total'] = float(tensions['T_total'][b])
            info['T_phys'] = float(tensions['T_phys'][b])
            all_results.append(info)
            all_gens_list.append(gens_batch[b])
            print(f"    Seed {seed:3d}: {info['category']:20s} "
                  f"closure={info['closure_frac']:.1%} "
                  f"simple={info['simplicity']:.4f} "
                  f"h_v(3)={info['h_dual_ir3']:.2f} "
                  f"T={info['T_total']:.1f}")
        dt_class = time.time() - t_class
        print(f"  Classification: {dt_class:.1f}s")

        # Save intermediate checkpoint
        save_path = os.path.join(os.path.dirname(__file__), args.save)
        save_kwargs = dict(
            results=json.dumps(all_results),
            n_seeds_done=len(all_results),
            n_seeds_total=actual_seeds,
            config=json.dumps(vars(args)),
        )
        if args.save_gens:
            save_kwargs['generators'] = np.array(all_gens_list)
        np.savez(save_path, **save_kwargs)
        print(f"  Checkpoint saved ({len(all_results)}/{actual_seeds} seeds)")
        print()

    elapsed = time.time() - t0

    # ============================================================
    # REPORT
    # ============================================================
    categories = print_report(all_results, actual_seeds, elapsed)

    # ============================================================
    # STABILITY TEST
    # ============================================================
    if not args.no_stability and 'F4' in categories:
        print(f"\n{'=' * 80}")
        print("  STABILITY TEST")
        print(f"{'=' * 80}")

        # Pick best F4 (lowest simplicity = closest to K proportional to I)
        f4_results = categories['F4']
        best = min(f4_results, key=lambda r: r['simplicity'])
        best_idx = all_results.index(best)
        best_gens = all_gens_list[best_idx]

        print(f"  Source: seed {best['seed']} "
              f"(simplicity={best['simplicity']:.6f}, "
              f"h_v(3)={best['h_dual_ir3']:.2f})")
        print(f"  Recovery: {args.stability_steps} steps, "
              f"gamma={args.s2_gamma}, lr={args.s2_lr}")

        epsilons = [1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0]
        gens_list = [best_gens[k] for k in range(N_GEN)]
        stab_results = stability_test(
            gens_list, epsilons,
            recovery_steps=args.stability_steps,
            gamma=args.s2_gamma, lr=args.s2_lr,
        )

        print_stability_report(stab_results)

        # Save stability results too
        stab_kwargs = dict(
            results=json.dumps(all_results),
            stability=json.dumps(stab_results),
            n_seeds_done=actual_seeds,
            n_seeds_total=actual_seeds,
            config=json.dumps(vars(args)),
        )
        if args.save_gens:
            stab_kwargs['generators'] = np.array(all_gens_list)
        np.savez(save_path, **stab_kwargs)

    elif not args.no_stability:
        print(f"\n  No F4 found in {actual_seeds} seeds -- skipping stability test.")

    total_elapsed = time.time() - t0
    print(f"\n{'=' * 80}")
    print(f"  DONE. Total: {total_elapsed:.0f}s ({total_elapsed / 60:.1f} min)")
    print(f"  Results: {save_path}")
    print(f"{'=' * 80}")


if __name__ == '__main__':
    main()
