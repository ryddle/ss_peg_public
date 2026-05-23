"""
128-seed Monte Carlo landscape for F₄ in so(27) — Cloud-adapted.

Extension of §11f (32-seed study) to 128 seeds for robust statistics.
Adapted for T4 GPUs on Lightning.ai / Google Colab:
  - Smaller default batch_size (8) for T4 VRAM safety
  - Auto-resume from checkpoint on re-run (--resume)
  - ETA reporting per batch
  - Graceful Ctrl+C handling (saves checkpoint before exit)

Estimated time on T4: ~17-21 hours total.

Usage:
  # Full 128-seed run (recommended for cloud)
  python _landscape_128seeds_cloud.py

  # Resume after interruption
  python _landscape_128seeds_cloud.py --resume

  # Quick local test (2 seeds, reduced steps)
  python _landscape_128seeds_cloud.py --quick

  # Custom configuration
  python _landscape_128seeds_cloud.py --n-seeds 64 --batch-size 4
"""

import sys, os, time, argparse, json, signal
import numpy as np
import torch

sys.path.insert(0, os.path.dirname(__file__))

from _landscape_montecarlo import (
    run_batch, classify_solution, stability_test,
    print_report, print_stability_report,
    N_GEN, DIM, ALPHA, BETA,
)
from shared_utils_gpu import DEVICE


# ============================================================
# CHECKPOINT I/O
# ============================================================

def load_checkpoint(path):
    """Load checkpoint and return (results_list, gens_list, n_done)."""
    if not os.path.exists(path):
        return [], [], 0
    data = np.load(path, allow_pickle=True)
    results = json.loads(str(data['results']))
    n_done = int(data['n_seeds_done'])
    gens = list(data['generators']) if 'generators' in data else []
    print(f"  Resumed from checkpoint: {n_done} seeds completed")
    return results, gens, n_done


def save_checkpoint(path, results, gens_list, n_done, n_total, config,
                    save_gens=True, stability=None):
    """Save checkpoint with all accumulated results."""
    kw = dict(
        results=json.dumps(results),
        n_seeds_done=n_done,
        n_seeds_total=n_total,
        config=json.dumps(config),
    )
    if save_gens and gens_list:
        kw['generators'] = np.array(gens_list)
    if stability is not None:
        kw['stability'] = json.dumps(stability)
    np.savez(path, **kw)


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='128-seed Monte Carlo landscape for F4 in so(27) [Cloud]')
    parser.add_argument('--n-seeds', type=int, default=128,
                        help='Total number of independent seeds (default: 128)')
    parser.add_argument('--batch-size', type=int, default=8,
                        help='Seeds per GPU batch (default: 8, safe for T4)')
    parser.add_argument('--s1-steps', type=int, default=8000,
                        help='Stage 1 optimizer steps')
    parser.add_argument('--s2-steps', type=int, default=3000,
                        help='Stage 2 optimizer steps')
    parser.add_argument('--s1-gamma', type=float, default=50.0)
    parser.add_argument('--s2-gamma', type=float, default=200.0)
    parser.add_argument('--s1-lr', type=float, default=0.001)
    parser.add_argument('--s2-lr', type=float, default=0.0001)
    parser.add_argument('--stability-steps', type=int, default=5000)
    parser.add_argument('--no-stability', action='store_true')
    parser.add_argument('--save', type=str,
                        default='landscape_128seeds_results.npz')
    parser.add_argument('--resume', action='store_true',
                        help='Resume from existing checkpoint')
    parser.add_argument('--quick', action='store_true',
                        help='Quick test: 2 seeds, 2k+1k steps')
    args = parser.parse_args()

    if args.quick:
        args.n_seeds = 2
        args.s1_steps = 2000
        args.s2_steps = 1000
        args.stability_steps = 1000
        args.batch_size = 2

    save_path = os.path.join(os.path.dirname(__file__), args.save)

    # Build seed list: deterministic sequence [0, 1, ..., n_seeds-1]
    n_seeds = args.n_seeds
    batch_size = args.batch_size
    n_batches = (n_seeds + batch_size - 1) // batch_size
    all_seeds = list(range(n_batches * batch_size))[:n_seeds]

    # Partition into batches
    seed_batches = [all_seeds[i:i + batch_size]
                    for i in range(0, len(all_seeds), batch_size)]

    # Resume from checkpoint if requested
    all_results = []
    all_gens_list = []
    n_done = 0
    if args.resume:
        all_results, all_gens_list, n_done = load_checkpoint(save_path)

    # Find which batch to start from
    start_batch = 0
    for i, batch_seeds in enumerate(seed_batches):
        if all(any(r['seed'] == s for r in all_results) for s in batch_seeds):
            start_batch = i + 1
        else:
            break

    remaining_batches = len(seed_batches) - start_batch

    print("=" * 80)
    print("  128-SEED MONTE CARLO LANDSCAPE: F4 in so(27)")
    print("=" * 80)
    print(f"  Seeds:     {n_seeds} ({len(seed_batches)} batches, bs={batch_size})")
    print(f"  Completed: {n_done} seeds ({start_batch} batches)")
    print(f"  Remaining: {n_seeds - n_done} seeds ({remaining_batches} batches)")
    print(f"  Stage 1:   {args.s1_steps} steps, gamma={args.s1_gamma}, lr={args.s1_lr}")
    print(f"  Stage 2:   {args.s2_steps} steps, gamma={args.s2_gamma}, lr={args.s2_lr}")
    print(f"  Density:   {N_GEN}/{DIM * (DIM - 1) // 2} = "
          f"{N_GEN / (DIM * (DIM - 1) // 2) * 100:.1f}%")
    print(f"  Device:    {DEVICE}")
    if torch.cuda.is_available():
        print(f"  GPU:       {torch.cuda.get_device_name()}")
        print(f"  VRAM:      {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    print()

    if remaining_batches == 0:
        print("  All seeds already completed. Use without --resume to re-run.")
        # Still print report from loaded data
        elapsed = 0
        print_report(all_results, n_seeds, elapsed)
        return

    # Graceful interrupt handler
    interrupted = [False]

    def handler(sig, frame):
        if interrupted[0]:
            print("\n  Force quit (checkpoint already saved).")
            sys.exit(1)
        interrupted[0] = True
        print("\n  Interrupt received — saving checkpoint after current batch...")

    signal.signal(signal.SIGINT, handler)

    # ============================================================
    # OPTIMIZATION LOOP
    # ============================================================
    t0 = time.time()
    batch_times = []

    for batch_idx in range(start_batch, len(seed_batches)):
        if interrupted[0]:
            break

        seeds = seed_batches[batch_idx]
        batch_num = batch_idx + 1
        total_batches = len(seed_batches)
        print(f"--- Batch {batch_num}/{total_batches}: seeds {seeds} ---")

        t_batch = time.time()
        gens_batch, tensions = run_batch(
            seeds, args.s1_steps, args.s2_steps,
            args.s1_gamma, args.s2_gamma, args.s1_lr, args.s2_lr,
            verbose=True,
        )
        dt_batch = time.time() - t_batch
        batch_times.append(dt_batch)

        # Classify each seed
        print("  Classifying...", flush=True)
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

        n_done = len(all_results)

        # Save checkpoint
        save_checkpoint(save_path, all_results, all_gens_list,
                        n_done, n_seeds, vars(args), save_gens=True)
        print(f"  Checkpoint saved ({n_done}/{n_seeds} seeds)")

        # ETA
        avg_batch_time = np.mean(batch_times)
        remaining = total_batches - batch_num
        eta_seconds = avg_batch_time * remaining
        eta_h = eta_seconds / 3600
        print(f"  Batch time: {dt_batch:.0f}s | "
              f"Avg: {avg_batch_time:.0f}s | "
              f"ETA: {eta_h:.1f}h ({remaining} batches left)")
        print()

    elapsed = time.time() - t0

    if interrupted[0]:
        print(f"\n  Interrupted after {n_done}/{n_seeds} seeds. "
              f"Re-run with --resume to continue.")
        return

    # ============================================================
    # REPORT
    # ============================================================
    categories = print_report(all_results, n_seeds, elapsed)

    # ============================================================
    # STABILITY TEST
    # ============================================================
    if not args.no_stability and 'F4' in categories:
        print(f"\n{'=' * 80}")
        print("  STABILITY TEST")
        print(f"{'=' * 80}")

        f4_results = categories['F4']
        best = min(f4_results, key=lambda r: r['simplicity'])
        best_idx = all_results.index(best)
        best_gens = all_gens_list[best_idx]

        print(f"  Source: seed {best['seed']} "
              f"(simplicity={best['simplicity']:.6f}, "
              f"h_v(3)={best['h_dual_ir3']:.2f})")

        epsilons = [1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0]
        gens_list_best = [best_gens[k] for k in range(N_GEN)]
        stab_results = stability_test(
            gens_list_best, epsilons,
            recovery_steps=args.stability_steps,
            gamma=args.s2_gamma, lr=args.s2_lr,
        )

        print_stability_report(stab_results)

        save_checkpoint(save_path, all_results, all_gens_list,
                        n_seeds, n_seeds, vars(args),
                        save_gens=True, stability=stab_results)

    elif not args.no_stability:
        print(f"\n  No F4 in {n_seeds} seeds — skipping stability test.")

    total_elapsed = time.time() - t0
    print(f"\n{'=' * 80}")
    print(f"  DONE. Total: {total_elapsed:.0f}s ({total_elapsed / 3600:.1f}h)")


if __name__ == "__main__":
    main()
