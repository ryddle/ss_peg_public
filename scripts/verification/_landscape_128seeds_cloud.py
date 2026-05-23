"""
128-seed Monte Carlo landscape for F₄ in so(27) — Resumable edition.

Extension of §11f (32-seed study) to 128 seeds for robust statistics.
Fine-grained checkpointing every N steps allows interruption and resume
from the exact point where it stopped — no wasted computation.

Checkpoint format (landscape_*_checkpoint.pt):
  - params: current optimizer tensor (B, n_gen, n_params)
  - optimizer_state: Adam state_dict (includes momentum buffers)
  - batch_idx: which batch is in progress
  - stage: 1 or 2 (optimizer stage)
  - step: current step within stage
  - completed_results: JSON list of finished seed results
  - completed_gens: list of numpy arrays for finished seeds
  - config: JSON of all hyperparameters

Usage:
  python _landscape_128seeds_cloud.py                               # full 128
  python _landscape_128seeds_cloud.py --n-seeds 16 --batch-size 4   # 16-seed pilot
  python _landscape_128seeds_cloud.py --resume                      # continue
  python _landscape_128seeds_cloud.py --quick                       # 2-seed test
"""

import sys, os, time, argparse, json, signal, io
import numpy as np
import torch

sys.path.insert(0, os.path.dirname(__file__))

from _landscape_montecarlo import (
    classify_solution, stability_test,
    print_report, print_stability_report,
    N_GEN, DIM, ALPHA, BETA,
)
from shared_utils_gpu import (
    batched_killing_metric,
    batched_non_commutativity_tension,
    batched_norm_tension,
    params_to_antisymmetric,
    antisymmetric_to_params,
    random_params_so,
    DTYPE, FDTYPE, DEVICE,
)
from test_closure_driven import batched_closure_tension


# ============================================================
# CHECKPOINT I/O (step-level)
# ============================================================

def save_step_checkpoint(path, params, optimizer, batch_idx, stage, step,
                         completed_results, completed_gens, config):
    """Save fine-grained checkpoint including optimizer state."""
    data = {
        'params': params.detach().cpu(),
        'optimizer_state': optimizer.state_dict(),
        'batch_idx': batch_idx,
        'stage': stage,
        'step': step,
        'completed_results': json.dumps(completed_results),
        'completed_gens': [g if isinstance(g, np.ndarray) else g
                           for g in completed_gens],
        'config': json.dumps(config),
    }

    # Atomic write: write to temp then rename
    tmp_path = path + '.tmp'
    torch.save(data, tmp_path)
    if os.path.exists(path):
        os.replace(tmp_path, path)
    else:
        os.rename(tmp_path, path)


def load_step_checkpoint(path, device):
    """Load step-level checkpoint. Returns dict or None."""
    if not os.path.exists(path):
        return None
    data = torch.load(path, map_location='cpu', weights_only=False)
    data['params'] = data['params'].to(device)
    data['completed_results'] = json.loads(data['completed_results'])
    data['config'] = json.loads(data['config'])
    return data


def save_final_results(path, all_results, all_gens_list, n_seeds, config,
                       stability=None):
    """Save final .npz results (backward-compatible format)."""
    kw = dict(
        results=json.dumps(all_results),
        n_seeds_done=len(all_results),
        n_seeds_total=n_seeds,
        config=json.dumps(config),
    )
    if all_gens_list:
        kw['generators'] = np.array(all_gens_list)
    if stability is not None:
        kw['stability'] = json.dumps(stability)
    np.savez(path, **kw)


# ============================================================
# INLINE BATCH OPTIMIZER (with step-level checkpointing)
# ============================================================

def run_batch_resumable(seeds, s1_steps, s2_steps, s1_gamma, s2_gamma,
                        s1_lr, s2_lr, batch_idx,
                        completed_results, completed_gens,
                        config, checkpoint_path,
                        checkpoint_every=1000,
                        resume_stage=None, resume_step=None,
                        resume_params=None, resume_optimizer_state=None,
                        interrupted=None, verbose=True):
    """Run optimization for a batch with intra-batch checkpointing.

    Returns: (gens_batch_np, tensions_dict) or (None, None) on interruption.
    """
    B = len(seeds)
    device = DEVICE

    # --- Initialize or restore params ---
    if resume_params is not None:
        params = resume_params.clone().detach().requires_grad_(True)
    else:
        params = random_params_so(B, N_GEN, DIM, device, seeds)
        # Normalize initial generators
        with torch.no_grad():
            gens = params_to_antisymmetric(params, N_GEN, dim=DIM)
            norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
            gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
            params = antisymmetric_to_params(gens)
        params = params.clone().detach().requires_grad_(True)

    K_target = -1.0 * torch.eye(N_GEN, dtype=FDTYPE, device=device)

    stages = [
        (1, s1_steps, s1_lr, s1_gamma),
        (2, s2_steps, s2_lr, s2_gamma),
    ]

    for stage_num, steps, lr, gamma in stages:
        # Skip already-completed stages on resume
        if resume_stage is not None and stage_num < resume_stage:
            continue

        start_step = 0
        if resume_stage is not None and stage_num == resume_stage and resume_step is not None:
            start_step = resume_step

        if verbose:
            print(f"  Stage {stage_num}: gamma={gamma}, lr={lr}, {steps} steps"
                  + (f" (resuming from step {start_step})" if start_step > 0 else ""))

        optimizer = torch.optim.Adam([params], lr=lr)

        # Restore optimizer state if resuming in this stage
        if (resume_optimizer_state is not None
                and resume_stage is not None
                and stage_num == resume_stage):
            try:
                optimizer.load_state_dict(resume_optimizer_state)
            except Exception as e:
                print(f"  WARNING: Could not restore optimizer state: {e}")
                print(f"  Continuing with fresh Adam (momentum lost)")
            resume_optimizer_state = None  # only use once

        for step in range(start_step, steps):
            if interrupted is not None and interrupted[0]:
                # Save checkpoint before breaking
                save_step_checkpoint(
                    checkpoint_path, params, optimizer,
                    batch_idx, stage_num, step,
                    completed_results, completed_gens, config)
                print(f"\n  Checkpoint saved mid-batch: batch={batch_idx}, "
                      f"stage={stage_num}, step={step}")
                return None, None  # signal interruption

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

            # Periodic checkpoint
            if (step + 1) % checkpoint_every == 0:
                save_step_checkpoint(
                    checkpoint_path, params, optimizer,
                    batch_idx, stage_num, step + 1,
                    completed_results, completed_gens, config)
                if verbose and (step + 1) % 2000 != 0:  # avoid double-print
                    print(f"    [checkpoint at stage {stage_num}, step {step + 1}]",
                          flush=True)

        # Stage done — clear resume markers
        resume_stage = None
        resume_step = None

    # --- Compute final tensions ---
    with torch.no_grad():
        generators = params_to_antisymmetric(params, N_GEN, dim=DIM)
        T_nc = batched_non_commutativity_tension(generators).cpu().numpy()
        K = batched_killing_metric(generators)
        T_k = ((K - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1)).cpu().numpy()
        T_norm = batched_norm_tension(generators).cpu().numpy()
        T_cl = batched_closure_tension(generators).cpu().numpy()

        final_gamma = s2_gamma
        T_total = T_nc + ALPHA * T_k + BETA * T_norm + final_gamma * T_cl
        T_phys = T_nc + ALPHA * T_k + BETA * T_norm
        all_gens = generators.cpu().numpy()

    return all_gens, {
        'T_nc': T_nc, 'T_k': T_k, 'T_norm': T_norm,
        'T_cl': T_cl, 'T_total': T_total, 'T_phys': T_phys,
    }


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='128-seed Monte Carlo landscape for F4 in so(27) [Resumable]')
    parser.add_argument('--n-seeds', type=int, default=128,
                        help='Total number of independent seeds (default: 128)')
    parser.add_argument('--batch-size', type=int, default=8,
                        help='Seeds per GPU batch (default: 8)')
    parser.add_argument('--s1-steps', type=int, default=8000,
                        help='Stage 1 optimizer steps')
    parser.add_argument('--s2-steps', type=int, default=3000,
                        help='Stage 2 optimizer steps')
    parser.add_argument('--s1-gamma', type=float, default=50.0)
    parser.add_argument('--s2-gamma', type=float, default=200.0)
    parser.add_argument('--s1-lr', type=float, default=0.001)
    parser.add_argument('--s2-lr', type=float, default=0.0001)
    parser.add_argument('--checkpoint-every', type=int, default=1000,
                        help='Save checkpoint every N steps (default: 1000)')
    parser.add_argument('--stability-steps', type=int, default=5000)
    parser.add_argument('--no-stability', action='store_true')
    parser.add_argument('--save', type=str,
                        default='landscape_128seeds_results.npz')
    parser.add_argument('--resume', action='store_true',
                        help='Resume from checkpoint (auto-detected)')
    parser.add_argument('--quick', action='store_true',
                        help='Quick test: 2 seeds, 500+250 steps')
    args = parser.parse_args()

    if args.quick:
        args.n_seeds = 2
        args.s1_steps = 500
        args.s2_steps = 250
        args.stability_steps = 500
        args.batch_size = 2
        args.checkpoint_every = 250

    save_dir = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(save_dir, args.save)
    checkpoint_path = os.path.join(save_dir,
                                   args.save.replace('.npz', '_checkpoint.pt'))

    # Build seed list
    n_seeds = args.n_seeds
    batch_size = args.batch_size
    n_batches = (n_seeds + batch_size - 1) // batch_size
    all_seeds_list = list(range(n_batches * batch_size))[:n_seeds]
    seed_batches = [all_seeds_list[i:i + batch_size]
                    for i in range(0, len(all_seeds_list), batch_size)]

    config = vars(args)

    # --- Resume logic ---
    all_results = []
    all_gens_list = []
    resume_batch_idx = None
    resume_stage = None
    resume_step = None
    resume_params = None
    resume_optimizer_state = None

    if args.resume:
        # Try step-level checkpoint first (highest priority)
        ckpt = load_step_checkpoint(checkpoint_path, DEVICE)
        if ckpt is not None:
            all_results = ckpt['completed_results']
            all_gens_list = list(ckpt['completed_gens']) if ckpt['completed_gens'] else []
            resume_batch_idx = ckpt['batch_idx']
            resume_stage = ckpt['stage']
            resume_step = ckpt['step']
            resume_params = ckpt['params']
            resume_optimizer_state = ckpt['optimizer_state']
            print(f"  Resumed: batch={resume_batch_idx}, stage={resume_stage}, "
                  f"step={resume_step}, completed={len(all_results)} seeds")
        elif os.path.exists(save_path):
            # Fall back to final .npz (batch-level)
            data = np.load(save_path, allow_pickle=True)
            all_results = json.loads(str(data['results']))
            all_gens_list = list(data['generators']) if 'generators' in data else []
            print(f"  Resumed from .npz: {len(all_results)} seeds completed")

    # Find start batch
    n_done = len(all_results)
    start_batch = 0
    if resume_batch_idx is not None:
        start_batch = resume_batch_idx  # resume mid-batch
    else:
        for i, batch_seeds in enumerate(seed_batches):
            if all(any(r['seed'] == s for r in all_results) for s in batch_seeds):
                start_batch = i + 1
            else:
                break

    remaining_batches = len(seed_batches) - start_batch

    print("=" * 80)
    print("  128-SEED MONTE CARLO LANDSCAPE: F4 in so(27) [Resumable]")
    print("=" * 80)
    print(f"  Seeds:       {n_seeds} ({len(seed_batches)} batches, bs={batch_size})")
    print(f"  Completed:   {n_done} seeds ({start_batch} batches)")
    print(f"  Remaining:   {n_seeds - n_done} seeds ({remaining_batches} batches)")
    print(f"  Stage 1:     {args.s1_steps} steps, gamma={args.s1_gamma}, lr={args.s1_lr}")
    print(f"  Stage 2:     {args.s2_steps} steps, gamma={args.s2_gamma}, lr={args.s2_lr}")
    print(f"  Checkpoint:  every {args.checkpoint_every} steps")
    print(f"  Density:     {N_GEN}/{DIM * (DIM - 1) // 2} = "
          f"{N_GEN / (DIM * (DIM - 1) // 2) * 100:.1f}%")
    print(f"  Device:      {DEVICE}")
    if torch.cuda.is_available():
        print(f"  GPU:         {torch.cuda.get_device_name()}")
        print(f"  VRAM:        {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    print()

    if remaining_batches == 0:
        print("  All seeds already completed. Use without --resume to re-run.")
        elapsed = 0
        print_report(all_results, n_seeds, elapsed)
        return

    # Graceful interrupt handler
    interrupted = [False]

    def handler(sig, frame):
        if interrupted[0]:
            print("\n  Force quit.")
            sys.exit(1)
        interrupted[0] = True
        print("\n  Interrupt received — will save checkpoint at next opportunity...")

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

        # Determine resume state for this batch
        r_stage = resume_stage if batch_idx == resume_batch_idx else None
        r_step = resume_step if batch_idx == resume_batch_idx else None
        r_params = resume_params if batch_idx == resume_batch_idx else None
        r_optim = resume_optimizer_state if batch_idx == resume_batch_idx else None

        gens_batch, tensions = run_batch_resumable(
            seeds, args.s1_steps, args.s2_steps,
            args.s1_gamma, args.s2_gamma, args.s1_lr, args.s2_lr,
            batch_idx=batch_idx,
            completed_results=all_results,
            completed_gens=all_gens_list,
            config=config,
            checkpoint_path=checkpoint_path,
            checkpoint_every=args.checkpoint_every,
            resume_stage=r_stage,
            resume_step=r_step,
            resume_params=r_params,
            resume_optimizer_state=r_optim,
            interrupted=interrupted,
            verbose=True,
        )

        # Clear resume state after first batch
        resume_batch_idx = None
        resume_stage = None
        resume_step = None
        resume_params = None
        resume_optimizer_state = None

        if gens_batch is None:
            # Interrupted mid-batch — checkpoint already saved
            break

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

        # Save batch-level results (.npz)
        save_final_results(save_path, all_results, all_gens_list, n_seeds, config)
        print(f"  Results saved ({n_done}/{n_seeds} seeds)")

        # Remove step-level checkpoint (batch is done)
        if os.path.exists(checkpoint_path):
            os.remove(checkpoint_path)

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
        print(f"\n  Interrupted after {len(all_results)}/{n_seeds} seeds. "
              f"Re-run with --resume to continue.")
        return

    # ============================================================
    # REPORT
    # ============================================================
    categories = print_report(all_results, n_seeds, elapsed)

    # Clean up step checkpoint
    if os.path.exists(checkpoint_path):
        os.remove(checkpoint_path)

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
              f"(simplicity={best['simplicity']:.6f})")

        epsilons = [1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1]
        stab_results = stability_test(
            list(best_gens), epsilons,
            recovery_steps=args.stability_steps,
            verbose=True,
        )
        print_stability_report(stab_results)

        # Save with stability
        save_final_results(save_path, all_results, all_gens_list,
                           n_seeds, config, stability=stab_results)

    print(f"\n{'=' * 80}")
    print(f"  DONE. Total: {elapsed:.0f}s ({elapsed / 3600:.1f}h)")


if __name__ == "__main__":
    main()
