"""
Bifurcation Perturbation Test: measure F4 basin solid angle at the saddle.

Strategy — STACKED SINGLE-BATCH:
  1. Evolve seed 256 (bs=1) to step 1000 — still in the deterministic regime
  2. Create ALL replicas across ALL sigma levels in ONE batch
  3. Run ONE optimization for 7000 S1 + 3000 S2 steps
  4. Classify all outcomes and group by sigma → P(F4 | σ)

The diagnostic test showed:
  - Steps 0–1000: bs=1 and bs=3 trajectories are BIT-IDENTICAL (Δp ~ 1e-17)
  - Steps 1000–1500: chaotic amplification (Δp jumps from 1e-17 to 1e-4)
  - After step 1500: macroscopic divergence → different basins

This test probes the geometry of the F4 basin boundary at the critical region.
"""

import sys, os, time, argparse
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
import json

N_GEN = 52
DIM   = 27
ALPHA = 5.0
BETA  = 1.0


def evolve_to_checkpoint(seed, steps, gamma, lr):
    """Evolve seed for `steps` iterations, return params tensor."""
    params = random_params_so(1, N_GEN, DIM, DEVICE, [seed])

    with torch.no_grad():
        gens = params_to_antisymmetric(params, N_GEN, dim=DIM)
        norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
        gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
        params = antisymmetric_to_params(gens)
    params = params.clone().detach().requires_grad_(True)

    K_target = -1.0 * torch.eye(N_GEN, dtype=FDTYPE, device=DEVICE)
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

        if step % 200 == 0:
            with torch.no_grad():
                tcl = T_cl[0].item()
                tt = T_total[0].item()
            print(f"    step {step:5d}: T_total={tt:.4e} T_cl={tcl:.4e}",
                  flush=True)

    return params.detach().clone()


def create_stacked_batch(params_checkpoint, sigmas, n_replicas):
    """Create one big batch with n_replicas per sigma level.

    Returns:
        params: (B_total, N_GEN, n_params) tensor
        sigma_labels: list of sigma values, one per batch element
    """
    param_scale = params_checkpoint.abs().mean().item()
    all_params = []
    sigma_labels = []

    for sigma in sigmas:
        base = params_checkpoint.expand(n_replicas, -1, -1).clone()
        if sigma > 0:
            torch.manual_seed(hash(f"bifurcation_{sigma}") % (2**31))
            noise = torch.randn_like(base) * sigma * param_scale
            base = base + noise
        all_params.append(base)
        sigma_labels.extend([sigma] * n_replicas)

    stacked = torch.cat(all_params, dim=0)   # (B_total, N_GEN, n_params)
    return stacked, sigma_labels


def run_stacked_optimization(params_all, s1_remaining, s2_steps,
                             s1_gamma, s2_gamma, s1_lr, s2_lr,
                             max_bs=16, verbose=True):
    """Run optimization in chunks of max_bs to avoid OOM.

    Returns final generators for all elements.
    """
    B_total = params_all.shape[0]
    n_chunks = (B_total + max_bs - 1) // max_bs
    K_target = -1.0 * torch.eye(N_GEN, dtype=FDTYPE, device=DEVICE)

    all_gens = []

    for chunk_idx in range(n_chunks):
        i0 = chunk_idx * max_bs
        i1 = min(i0 + max_bs, B_total)
        chunk_size = i1 - i0

        if verbose:
            print(f"\n  Chunk {chunk_idx+1}/{n_chunks} "
                  f"(elements {i0}-{i1-1}, bs={chunk_size})", flush=True)

        params = params_all[i0:i1].clone().requires_grad_(True)

        for stage, (steps, lr, gamma) in enumerate([
            (s1_remaining, s1_lr, s1_gamma),
            (s2_steps, s2_lr, s2_gamma),
        ], 1):
            if verbose:
                print(f"    Stage {stage}: γ={gamma}, lr={lr}, {steps} steps",
                      flush=True)
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
                        tcl = T_cl.cpu().numpy()
                        print(f"      step {step:5d}: T_cl avg={tcl.mean():.4e} "
                              f"min={tcl.min():.4e} max={tcl.max():.4e}",
                              flush=True)

        with torch.no_grad():
            generators = params_to_antisymmetric(params, N_GEN, dim=DIM)
            all_gens.append(generators.cpu().numpy())

    return np.concatenate(all_gens, axis=0)   # (B_total, N_GEN, DIM, DIM)


def classify_solution(gens_np):
    """Classify generators (same logic as _landscape_montecarlo.py)."""
    gens_np = np.asarray(gens_np)
    n = len(gens_np)
    d = gens_np[0].shape[0]

    idx_i, idx_j = np.triu_indices(n, k=1)
    G_i = gens_np[idx_i]
    G_j = gens_np[idx_j]
    C = G_i @ G_j - G_j @ G_i
    C_flat = C.reshape(-1, d * d).T

    B_mat = np.array([g.flatten() for g in gens_np]).T
    coeffs, _, _, _ = np.linalg.lstsq(B_mat, C_flat, rcond=None)
    proj = B_mat @ coeffs
    residual_norms = np.linalg.norm(C_flat - proj, axis=0)

    closure_frac = float(np.mean(residual_norms < 1e-4))
    max_res = float(residual_norms.max())

    K = killing_metric_np(list(gens_np))
    K_diag = np.diag(K).real
    K_mean = float(K_diag.mean())

    if abs(K_mean) > 1e-12:
        K_normed = K / K_mean
        simplicity = float(
            np.linalg.norm(K_normed - np.eye(n)) / np.linalg.norm(np.eye(n))
        )
    else:
        simplicity = 999.0

    gens_r = [g / np.linalg.norm(g) for g in gens_np]
    K_r = killing_metric_np(gens_r)
    Tr_sq = np.mean([np.trace(g @ g).real for g in gens_r])
    h_raw = float(abs(np.diag(K_r).real.mean() / Tr_sq)) if abs(Tr_sq) > 1e-12 else 0.0
    h_dual_ir3 = 3.0 * h_raw

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
        'simplicity': simplicity,
        'h_dual_ir3': h_dual_ir3,
        'max_residual': max_res,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=256)
    parser.add_argument('--checkpoint-step', type=int, default=1000,
                        help='Step at which to checkpoint (before bifurcation)')
    parser.add_argument('--n-replicas', type=int, default=4,
                        help='Number of perturbed replicas per noise level')
    parser.add_argument('--max-bs', type=int, default=16,
                        help='Max batch size per GPU chunk')
    parser.add_argument('--sigmas', type=float, nargs='+',
                        default=[0.0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
                        help='Noise levels (relative to param scale)')
    parser.add_argument('--save', type=str, default=None)
    args = parser.parse_args()

    S1_STEPS = 8000
    S2_STEPS = 3000
    S1_GAMMA = 50.0
    S2_GAMMA = 200.0
    S1_LR = 0.001
    S2_LR = 0.0001

    s1_remaining = S1_STEPS - args.checkpoint_step
    n_sigmas = len(args.sigmas)
    B_total = n_sigmas * args.n_replicas
    n_chunks = (B_total + args.max_bs - 1) // args.max_bs

    print("=" * 70)
    print("  BIFURCATION PERTURBATION TEST (stacked batch)")
    print("=" * 70)
    print(f"  Seed:            {args.seed}")
    print(f"  Checkpoint step: {args.checkpoint_step}")
    print(f"  Remaining S1:    {s1_remaining} steps")
    print(f"  S2:              {S2_STEPS} steps")
    print(f"  Sigma levels:    {args.sigmas}")
    print(f"  Replicas/sigma:  {args.n_replicas}")
    print(f"  Total replicas:  {B_total} ({n_chunks} chunk(s) of ≤{args.max_bs})")
    print(f"  Device:          {DEVICE}")
    print()

    # Phase 1: Evolve to checkpoint
    print(f"[1] Evolving seed {args.seed} to step {args.checkpoint_step}...")
    t0 = time.time()
    params_ckpt = evolve_to_checkpoint(
        args.seed, args.checkpoint_step, S1_GAMMA, S1_LR
    )
    t_ckpt = time.time() - t0
    print(f"    Checkpoint reached in {t_ckpt:.1f}s")
    print(f"    Params shape: {params_ckpt.shape}")
    print(f"    Params |mean|: {params_ckpt.abs().mean().item():.6e}")
    print(f"    Params std:    {params_ckpt.std().item():.6e}")

    # Phase 2: Create stacked batch
    print(f"\n[2] Creating stacked batch ({B_total} replicas)...")
    params_all, sigma_labels = create_stacked_batch(
        params_ckpt, args.sigmas, args.n_replicas
    )
    print(f"    Stacked params shape: {params_all.shape}")

    # Phase 3: Run optimization (chunked if needed)
    print(f"\n[3] Running optimization ({n_chunks} chunk(s))...")
    t_opt = time.time()
    all_gens = run_stacked_optimization(
        params_all, s1_remaining, S2_STEPS,
        S1_GAMMA, S2_GAMMA, S1_LR, S2_LR,
        max_bs=args.max_bs, verbose=True,
    )
    t_opt_elapsed = time.time() - t_opt
    print(f"\n  Optimization total: {t_opt_elapsed:.0f}s ({t_opt_elapsed/60:.1f} min)")

    # Phase 4: Classify all
    print(f"\n[4] Classifying {B_total} solutions...")
    t_class = time.time()
    all_results = []
    for i in range(B_total):
        result = classify_solution(all_gens[i])
        result['sigma'] = sigma_labels[i]
        all_results.append(result)
        cat = result['category']
        s = result['simplicity']
        h = result['h_dual_ir3']
        cl = result['closure_frac']
        print(f"    [{i:3d}] σ={sigma_labels[i]:.1e}  {cat:25s}  "
              f"closure={100*cl:.1f}%  simplicity={s:.4f}  h∨={h:.2f}",
              flush=True)
    t_class_elapsed = time.time() - t_class
    print(f"  Classification: {t_class_elapsed:.0f}s")

    # Phase 5: Group by sigma and report
    print(f"\n{'='*70}")
    print(f"  RESULTS BY SIGMA LEVEL")
    print(f"{'='*70}")

    summary = {}
    for sigma in args.sigmas:
        group = [r for r in all_results if r['sigma'] == sigma]
        cat_counts = {}
        for r in group:
            c = r['category']
            cat_counts[c] = cat_counts.get(c, 0) + 1

        n_f4 = cat_counts.get('F4', 0)
        p_f4 = n_f4 / len(group) if group else 0

        print(f"\n  σ = {sigma:.1e}  (N={len(group)}):")
        for cat in ['F4', 'SIMPLE_OTHER', 'CLOSED_NON_SIMPLE', 'PARTIAL', 'NOT_CLOSED']:
            count = cat_counts.get(cat, 0)
            if count > 0:
                pct = 100 * count / len(group)
                print(f"    {cat:25s}: {count:3d}/{len(group)} = {pct:5.1f}%")

        simps = [r['simplicity'] for r in group if r['closure_frac'] >= 0.99]
        if simps:
            print(f"    simplicity (closed): mean={np.mean(simps):.4f} "
                  f"min={np.min(simps):.4f} max={np.max(simps):.4f}")

        summary[sigma] = {
            'n': len(group), 'n_f4': n_f4, 'p_f4': p_f4,
            'categories': cat_counts,
        }

    # Summary table
    print(f"\n{'='*70}")
    print(f"  SUMMARY: P(F4) vs σ")
    print(f"{'='*70}")
    print(f"  {'σ':>12s}  {'P(F4)':>8s}  {'F4':>4s}  {'NonSimp':>8s}  "
          f"{'Partial':>8s}  {'NotCl':>6s}  {'N':>4s}")
    print(f"  {'-'*12}  {'-'*8}  {'-'*4}  {'-'*8}  {'-'*8}  {'-'*6}  {'-'*4}")
    for sigma in args.sigmas:
        s = summary[sigma]
        ns = s['categories'].get('CLOSED_NON_SIMPLE', 0)
        pa = s['categories'].get('PARTIAL', 0)
        nc = s['categories'].get('NOT_CLOSED', 0)
        print(f"  {sigma:12.1e}  {s['p_f4']:8.1%}  {s['n_f4']:4d}  "
              f"{ns:8d}  {pa:8d}  {nc:6d}  {s['n']:4d}")

    total_time = time.time() - t0
    print(f"\n  Total elapsed: {total_time:.0f}s ({total_time/60:.1f} min)")

    # Save
    if args.save:
        save_path = os.path.join(os.path.dirname(__file__), args.save)
        save_data = {
            'seed': args.seed,
            'checkpoint_step': args.checkpoint_step,
            'sigmas': args.sigmas,
            'n_replicas': args.n_replicas,
            'total_replicas': B_total,
            'results': all_results,
            'summary': {str(k): v for k, v in summary.items()},
            'optimization_time_s': t_opt_elapsed,
            'total_time_s': total_time,
        }
        np.savez_compressed(save_path,
                            results_json=json.dumps(save_data, default=str))
        print(f"  Saved to {save_path}")

    print(f"\n{'='*70}")
    print(f"  DONE")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
