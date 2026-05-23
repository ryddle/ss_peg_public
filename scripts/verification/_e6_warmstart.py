"""
E₆ warm-start: skip Phase A, load pre-converged F₄ from 128-seed landscape.

The 128-seed landscape campaign produced 68 CLOSED seeds with 100% closure.
These are authentic F₄ sub-algebras in so(27). We use them directly to skip
Phase A (~11k steps, ~55 min per trial) and go straight to Phase B+C.

Strategy:
  Phase B: Add 26 complement generators on top of frozen/semi-frozen F₄.
  Phase C: Joint annealing of all 78 generators.

Key improvements over _e6_f4_scaffold.py:
  - Phase A eliminated (warm-start from landscape)
  - Multiple F₄ bases from different simplicity bands
  - Gradual F₄ unfreezing schedule in Phase B
  - Multiple complement seeds per F₄ base
  - Fine-grained checkpointing

Usage:
  python _e6_warmstart.py                           # Default: top 5 F4, 2 comp seeds each
  python _e6_warmstart.py --n-f4 10 --n-comp 3      # Top 10 F4, 3 complement seeds
  python _e6_warmstart.py --resume                   # Resume from checkpoint
  python _e6_warmstart.py --quick                    # Quick test (1 trial)
  python _e6_warmstart.py --f4-seeds 62 35 97        # Specific F4 seeds
  python _e6_warmstart.py --unfreeze-schedule 0.0 0.01 0.05  # Gradual unfreeze
"""

import sys, os, time, argparse, json, signal
import numpy as np
import torch

sys.path.insert(0, os.path.dirname(__file__))

from shared_utils_gpu import (
    batched_killing_metric,
    batched_non_commutativity_tension,
    batched_norm_tension,
    params_to_antisymmetric,
    antisymmetric_to_params,
    DTYPE, FDTYPE, DEVICE,
)
from test_closure_driven import batched_closure_tension

N_GEN_F4 = 52
N_GEN_E6 = 78
N_GEN_COMPLEMENT = N_GEN_E6 - N_GEN_F4  # 26
DIM = 27
ALPHA = 5.0
BETA = 1.0


# ============================================================
# F₄ SELECTION
# ============================================================

def select_f4_seeds(landscape_path, n_f4=5, specific_seeds=None, diverse=True):
    """Select the best F₄ seeds from the 128-seed landscape.

    Args:
        landscape_path: path to landscape_128seeds_results.npz
        n_f4: how many F₄ bases to pick
        specific_seeds: list of seed indices to use (overrides n_f4)
        diverse: if True, pick seeds from different simplicity bands

    Returns:
        list of (seed_idx, seed_number, generators_52x27x27, info_dict)
    """
    data = np.load(landscape_path, allow_pickle=True)
    results = json.loads(str(data['results']))
    generators = data['generators']  # (128, 52, 27, 27)

    if specific_seeds is not None:
        selected = []
        for s in specific_seeds:
            for idx, r in enumerate(results):
                if r['seed'] == s:
                    cat = r.get('category', '')
                    if 'CLOSED' not in cat:
                        print(f"  WARNING: seed {s} is {cat}, not CLOSED")
                    selected.append((idx, r['seed'], generators[idx], r))
                    break
            else:
                print(f"  WARNING: seed {s} not found in landscape")
        return selected

    # Filter to CLOSED seeds only
    closed = [(idx, r) for idx, r in enumerate(results)
              if r.get('category', '') in ('CLOSED', 'CLOSED_NON_SIMPLE')]
    closed.sort(key=lambda x: x[1]['T_cl'])

    if not diverse:
        # Simple: just pick top N by T_cl
        selected = [(idx, r['seed'], generators[idx], r)
                     for idx, r in closed[:n_f4]]
        return selected

    # Diverse selection: pick from different simplicity bands
    band_centers = {
        'B1': 0.08, 'B2': 0.14, 'B3': 0.20,
        'B4': 0.25, 'B5': 0.29,
    }

    def get_band(s_val):
        best_band, best_dist = 'B2', 999.0
        for name, center in band_centers.items():
            d = abs(s_val - center)
            if d < best_dist:
                best_dist = d
                best_band = name
        return best_band

    # Group by band
    by_band = {}
    for idx, r in closed:
        band = get_band(r['simplicity'])
        by_band.setdefault(band, []).append((idx, r))

    # Round-robin across bands, picking best from each
    selected = []
    band_order = sorted(by_band.keys())
    band_ptrs = {b: 0 for b in band_order}

    while len(selected) < n_f4:
        added = False
        for band in band_order:
            if len(selected) >= n_f4:
                break
            ptr = band_ptrs[band]
            if ptr < len(by_band[band]):
                idx, r = by_band[band][ptr]
                selected.append((idx, r['seed'], generators[idx], r))
                band_ptrs[band] = ptr + 1
                added = True
        if not added:
            break

    return selected


# ============================================================
# CHECKPOINT I/O
# ============================================================

def save_checkpoint(path, trial_idx, phase, stage, step,
                    params_dict, optimizer, completed_trials,
                    f4_gens=None, config=None):
    data = {
        'trial_idx': trial_idx,
        'phase': phase,
        'stage': stage,
        'step': step,
        'optimizer_state': optimizer.state_dict(),
        'completed_trials': json.dumps(completed_trials),
        'config': json.dumps(config) if config else '{}',
    }
    for k, v in params_dict.items():
        data[f'param_{k}'] = v.detach().cpu()
    if f4_gens is not None:
        data['f4_gens'] = (np.array(f4_gens) if not isinstance(f4_gens, np.ndarray)
                           else f4_gens)

    tmp_path = path + '.tmp'
    torch.save(data, tmp_path)
    if os.path.exists(path):
        os.replace(tmp_path, path)
    else:
        os.rename(tmp_path, path)


def load_checkpoint(path, device):
    if not os.path.exists(path):
        return None
    data = torch.load(path, map_location='cpu', weights_only=False)
    data['completed_trials'] = json.loads(data['completed_trials'])
    data['config'] = json.loads(data['config'])
    for k in list(data.keys()):
        if k.startswith('param_') and isinstance(data[k], torch.Tensor):
            data[k] = data[k].to(device)
    return data


# ============================================================
# PHASE B: SCAFFOLD — GROW COMPLEMENT
# ============================================================

def phase_b_grow(f4_gens_np, complement_seed=999,
                 stage1_steps=20000, stage1_gamma=50.0, stage1_lr=0.001,
                 stage2_steps=20000, stage2_gamma=200.0, stage2_lr=0.0001,
                 unfreeze_schedule=None,
                 trial_idx=0, completed_trials=None,
                 config=None, checkpoint_path=None,
                 checkpoint_every=1000,
                 resume_stage=None, resume_step=None,
                 resume_f4_params=None, resume_complement_params=None,
                 resume_optimizer_state=None,
                 interrupted=None, verbose=True):
    """Grow 26 complement generators around F₄ scaffold.

    Args:
        unfreeze_schedule: list of f4_lr_factor values per segment.
            E.g. [0.0, 0.01, 0.05] splits each stage into 3 equal
            segments with progressively unfreezing F₄.
            Default: [0.0] (fully frozen).
    """
    device = DEVICE
    if completed_trials is None:
        completed_trials = []
    if unfreeze_schedule is None:
        unfreeze_schedule = [0.0]

    # Convert F₄ generators to parameter space
    if resume_f4_params is not None:
        f4_params = resume_f4_params.clone().detach().requires_grad_(True)
    else:
        f4_torch = torch.tensor(
            np.array(f4_gens_np), dtype=DTYPE, device=device
        ).unsqueeze(0)
        f4_params = antisymmetric_to_params(f4_torch).squeeze(0)
        f4_params = f4_params.clone().detach().requires_grad_(True)

    # Initialize complement generators
    if resume_complement_params is not None:
        complement_params = resume_complement_params.clone().detach().requires_grad_(True)
    else:
        n_p = DIM * (DIM - 1) // 2
        np.random.seed(complement_seed)
        complement_init = np.random.randn(N_GEN_COMPLEMENT, n_p) * 0.1
        complement_params = torch.tensor(
            complement_init, dtype=FDTYPE, device=device
        ).requires_grad_(True)

    K_target = -1.0 * torch.eye(N_GEN_E6, dtype=FDTYPE, device=device)

    stages = [
        (1, stage1_steps, stage1_lr, stage1_gamma),
        (2, stage2_steps, stage2_lr, stage2_gamma),
    ]

    for stage_num, total_steps, base_lr, gamma in stages:
        if resume_stage is not None and stage_num < resume_stage:
            continue

        start_step = 0
        if resume_stage is not None and stage_num == resume_stage and resume_step is not None:
            start_step = resume_step

        # Compute unfreeze segment boundaries
        n_segments = len(unfreeze_schedule)
        segment_len = total_steps // n_segments

        if verbose:
            sched_str = " -> ".join(f"{f:.3f}" for f in unfreeze_schedule)
            print(f"\n  Phase B Stage {stage_num}: gamma={gamma}, "
                  f"comp_lr={base_lr}, {total_steps} steps")
            print(f"    F4 unfreeze schedule: [{sched_str}] "
                  f"(segments of {segment_len} steps)")
            if start_step > 0:
                print(f"    (resuming from step {start_step})")

        # Initial optimizer — will be rebuilt at segment boundaries
        current_segment = start_step // segment_len if segment_len > 0 else 0
        current_segment = min(current_segment, n_segments - 1)
        f4_lr_factor = unfreeze_schedule[current_segment]

        param_groups = [
            {'params': [complement_params], 'lr': base_lr},
            {'params': [f4_params], 'lr': base_lr * f4_lr_factor},
        ]
        optimizer = torch.optim.Adam(param_groups)

        if (resume_optimizer_state is not None
                and resume_stage is not None
                and stage_num == resume_stage):
            try:
                optimizer.load_state_dict(resume_optimizer_state)
            except Exception:
                pass
            resume_optimizer_state = None

        for step in range(start_step, total_steps):
            if interrupted is not None and interrupted[0]:
                save_checkpoint(
                    checkpoint_path, trial_idx, 'B', stage_num, step,
                    {'f4_params': f4_params, 'complement_params': complement_params},
                    optimizer, completed_trials,
                    f4_gens=f4_gens_np, config=config)
                print(f"\n    Checkpoint saved: phase=B, stage={stage_num}, step={step}")
                return None, None

            # Check if we crossed a segment boundary
            seg = min(step // segment_len, n_segments - 1) if segment_len > 0 else 0
            if seg != current_segment:
                current_segment = seg
                f4_lr_factor = unfreeze_schedule[current_segment]
                # Rebuild optimizer with new F₄ lr
                param_groups = [
                    {'params': [complement_params], 'lr': base_lr},
                    {'params': [f4_params], 'lr': base_lr * f4_lr_factor},
                ]
                optimizer = torch.optim.Adam(param_groups)
                if verbose:
                    print(f"    [segment {seg}] f4_lr_factor -> {f4_lr_factor:.4f}",
                          flush=True)

            optimizer.zero_grad()
            all_params = torch.cat([f4_params, complement_params], dim=0)
            all_params_batch = all_params.unsqueeze(0)
            generators = params_to_antisymmetric(all_params_batch, N_GEN_E6, dim=DIM)

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
                    tcl = T_cl.item()
                    kd_std = torch.diagonal(K[0]).std().item()
                    gens_all = params_to_antisymmetric(
                        all_params_batch, N_GEN_E6, dim=DIM)
                    comp_norms = (gens_all[0, N_GEN_F4:].abs() ** 2
                                  ).sum(dim=(-2, -1)).sqrt().mean().item()
                    print(f"    step {step:5d}: T_cl={tcl:.4e} "
                          f"K_diag_std={kd_std:.4e} "
                          f"comp_norm={comp_norms:.4f} "
                          f"[f4_lr_f={f4_lr_factor:.4f}]",
                          flush=True)

            if checkpoint_path and (step + 1) % checkpoint_every == 0:
                save_checkpoint(
                    checkpoint_path, trial_idx, 'B', stage_num, step + 1,
                    {'f4_params': f4_params, 'complement_params': complement_params},
                    optimizer, completed_trials,
                    f4_gens=f4_gens_np, config=config)
                if verbose and (step + 1) % 2000 != 0:
                    print(f"    [checkpoint at stage {stage_num}, step {step + 1}]",
                          flush=True)

        resume_stage = None
        resume_step = None

    # Extract final generators
    with torch.no_grad():
        all_params = torch.cat([f4_params, complement_params], dim=0)
        final_gens = params_to_antisymmetric(
            all_params.unsqueeze(0), N_GEN_E6, dim=DIM)
        gens_78_np = final_gens[0].cpu().numpy()

    return gens_78_np, {'completed': True}


# ============================================================
# PHASE C: JOINT ANNEALING
# ============================================================

def phase_c_anneal(gens_78_np, steps=15000, gamma=200.0, lr=5e-5,
                   trial_idx=0, completed_trials=None,
                   f4_gens_np=None, config=None,
                   checkpoint_path=None, checkpoint_every=1000,
                   resume_step=None, resume_params=None,
                   resume_optimizer_state=None,
                   interrupted=None, verbose=True):
    device = DEVICE
    if completed_trials is None:
        completed_trials = []

    if verbose:
        print(f"\n  Phase C: Joint annealing ({steps} steps, lr={lr}, gamma={gamma})")

    if resume_params is not None:
        params_78 = resume_params.clone().detach().requires_grad_(True)
    else:
        gens_torch = torch.tensor(
            np.array(gens_78_np), dtype=DTYPE, device=device
        ).unsqueeze(0)
        params_78 = antisymmetric_to_params(gens_torch)
        params_78 = params_78.clone().detach().requires_grad_(True)

    K_target = -1.0 * torch.eye(N_GEN_E6, dtype=FDTYPE, device=device)

    start_step = 0 if resume_step is None else resume_step
    if verbose and start_step > 0:
        print(f"    (resuming from step {start_step})")

    optimizer = torch.optim.Adam([params_78], lr=lr)
    if resume_optimizer_state is not None:
        try:
            optimizer.load_state_dict(resume_optimizer_state)
        except Exception:
            pass

    for step in range(start_step, steps):
        if interrupted is not None and interrupted[0]:
            save_checkpoint(
                checkpoint_path, trial_idx, 'C', 1, step,
                {'params_78': params_78}, optimizer, completed_trials,
                f4_gens=f4_gens_np, config=config)
            print(f"\n    Checkpoint saved: phase=C, step={step}")
            return None

        optimizer.zero_grad()
        generators = params_to_antisymmetric(params_78, N_GEN_E6, dim=DIM)
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
                tcl = T_cl.item()
                kd_std = torch.diagonal(K[0]).std().item()
                print(f"    step {step:5d}: T_cl={tcl:.4e} K_diag_std={kd_std:.4e}",
                      flush=True)

        if checkpoint_path and (step + 1) % checkpoint_every == 0:
            save_checkpoint(
                checkpoint_path, trial_idx, 'C', 1, step + 1,
                {'params_78': params_78}, optimizer, completed_trials,
                f4_gens=f4_gens_np, config=config)
            if verbose and (step + 1) % 2000 != 0:
                print(f"    [checkpoint at step {step + 1}]", flush=True)

    with torch.no_grad():
        final_gens = params_to_antisymmetric(params_78, N_GEN_E6, dim=DIM)
        return final_gens[0].cpu().numpy()


# ============================================================
# ANALYSIS
# ============================================================

def analyze_e6(gens_78_np, name=""):
    """Detailed E₆ analysis."""
    from shared_utils_gpu import killing_metric_np, commutator_np

    gens_78_np = np.asarray(gens_78_np)
    n = len(gens_78_np)
    total_pairs = n * (n - 1) // 2

    # Vectorized closure
    idx_i, idx_j = np.triu_indices(n, k=1)
    G_i = gens_78_np[idx_i]
    G_j = gens_78_np[idx_j]
    C = G_i @ G_j - G_j @ G_i
    C_flat = C.reshape(-1, DIM ** 2).T
    B_mat = np.array([g.flatten() for g in gens_78_np]).T
    coeffs, _, _, _ = np.linalg.lstsq(B_mat, C_flat, rcond=None)
    proj = B_mat @ coeffs
    residuals = np.linalg.norm(C_flat - proj, axis=0)

    closure_frac = float(np.mean(residuals < 1e-4))
    max_res = float(residuals.max())

    # Killing metric
    K = killing_metric_np(list(gens_78_np))
    K_diag = np.diag(K).real
    K_mean = float(K_diag.mean())
    K_std = float(K_diag.std())

    # Simplicity
    if abs(K_mean) > 1e-12:
        K_normed = K / K_mean
        simplicity = float(
            np.linalg.norm(K_normed - np.eye(n)) / np.linalg.norm(np.eye(n)))
    else:
        simplicity = 999.0

    # h_dual
    gens_r = [g / np.linalg.norm(g) for g in gens_78_np]
    K_r = killing_metric_np(gens_r)
    Tr_sq = np.mean([np.trace(g @ g).real for g in gens_r])
    h_raw = (float(abs(np.diag(K_r).real.mean() / Tr_sq))
             if abs(Tr_sq) > 1e-12 else 0.0)
    h_dual_ir3 = 3.0 * h_raw

    # Rank estimation
    np.random.seed(42)
    coeffs_h = np.random.randn(n)
    H = sum(c * g for c, g in zip(coeffs_h, gens_78_np))
    ad_H = np.array([commutator_np(H, g).flatten() for g in gens_78_np]).T
    _, S, _ = np.linalg.svd(ad_H)
    rank_estimate = int(np.sum(S < 1e-6 * S[0]))

    print(f"\n{'=' * 60}")
    print(f"  {name}")
    print(f"{'=' * 60}")
    print(f"  Generators: {n}")
    print(f"  Closure: {closure_frac:.1%} ({int(closure_frac * total_pairs)}/{total_pairs})")
    print(f"  Max residual: {max_res:.4e}")
    print(f"  K_diag: mean={K_mean:.6f}  std={K_std:.6e}")
    print(f"  Simplicity: {simplicity:.6e}")
    print(f"  h_v (I_R=1): {h_raw:.4f}")
    print(f"  h_v (I_R=3): {h_dual_ir3:.4f}  (E6 expected: 36.0, F4: 9.0)")
    print(f"  Estimated rank: {rank_estimate}  (E6: 6, F4: 4)")

    is_closed = closure_frac >= 0.99
    is_simple = simplicity < 0.05
    is_e6 = is_closed and is_simple and 30.0 < h_dual_ir3 < 42.0 and rank_estimate >= 5
    is_f4 = is_closed and is_simple and 7.0 < h_dual_ir3 < 11.0

    if is_e6:
        verdict = "E6 (PASS)"
    elif is_f4:
        verdict = "F4 embedded in 78 gens"
    elif is_closed and is_simple:
        verdict = f"SIMPLE (h_v={h_dual_ir3:.1f}, rank~{rank_estimate})"
    elif is_closed:
        verdict = "CLOSED_NON_SIMPLE"
    else:
        verdict = f"INCOMPLETE (closure={closure_frac:.1%})"

    print(f"  VERDICT: {verdict}")

    return {
        'closure_frac': float(closure_frac),
        'max_residual': float(max_res),
        'K_mean': float(K_mean),
        'K_std': float(K_std),
        'simplicity': float(simplicity),
        'h_raw': float(h_raw),
        'h_dual_ir3': float(h_dual_ir3),
        'rank_estimate': int(rank_estimate),
        'verdict': verdict,
    }


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='E6 warm-start from pre-converged F4 landscape')

    # F₄ selection
    parser.add_argument('--landscape', type=str,
                        default='landscape_128seeds_results.npz',
                        help='Path to landscape .npz with converged F4')
    parser.add_argument('--n-f4', type=int, default=5,
                        help='Number of F4 bases to use (default: 5)')
    parser.add_argument('--f4-seeds', type=int, nargs='+', default=None,
                        help='Specific F4 seed numbers (overrides --n-f4)')
    parser.add_argument('--diverse', action='store_true', default=True,
                        help='Pick F4 from diverse bands (default: True)')
    parser.add_argument('--no-diverse', dest='diverse', action='store_false')

    # Complement
    parser.add_argument('--n-comp', type=int, default=2,
                        help='Complement seeds per F4 base (default: 2)')
    parser.add_argument('--comp-seed-base', type=int, default=5000,
                        help='Base for complement seed generation')

    # Phase B
    parser.add_argument('--b-s1-steps', type=int, default=20000)
    parser.add_argument('--b-s2-steps', type=int, default=20000)
    parser.add_argument('--b-s1-gamma', type=float, default=50.0)
    parser.add_argument('--b-s2-gamma', type=float, default=200.0)
    parser.add_argument('--b-s1-lr', type=float, default=0.001)
    parser.add_argument('--b-s2-lr', type=float, default=0.0001)
    parser.add_argument('--unfreeze-schedule', type=float, nargs='+',
                        default=[0.0, 0.0, 0.01],
                        help='F4 lr factors per segment (default: 0.0 0.0 0.01)')

    # Phase C
    parser.add_argument('--c-steps', type=int, default=15000)
    parser.add_argument('--c-lr', type=float, default=5e-5)
    parser.add_argument('--c-gamma', type=float, default=200.0)

    # General
    parser.add_argument('--checkpoint-every', type=int, default=1000)
    parser.add_argument('--save', type=str, default='e6_warmstart_results.npz')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--quick', action='store_true')

    args = parser.parse_args()

    if args.quick:
        args.n_f4 = 1
        args.n_comp = 1
        args.b_s1_steps = 3000
        args.b_s2_steps = 2000
        args.c_steps = 2000
        args.checkpoint_every = 500

    save_dir = os.path.dirname(os.path.abspath(__file__))
    landscape_path = os.path.join(save_dir, args.landscape)
    save_path = os.path.join(save_dir, args.save)
    checkpoint_path = os.path.join(
        save_dir, args.save.replace('.npz', '_checkpoint.pt'))

    # ---- Select F₄ bases ----
    print("=" * 80)
    print("  E6 WARM-START FROM PRE-CONVERGED F4")
    print("=" * 80)

    f4_selections = select_f4_seeds(
        landscape_path,
        n_f4=args.n_f4,
        specific_seeds=args.f4_seeds,
        diverse=args.diverse,
    )

    # Generate complement seeds
    comp_seeds_base = args.comp_seed_base
    trials = []
    for sel_idx, (arr_idx, seed_num, gens, info) in enumerate(f4_selections):
        for ci in range(args.n_comp):
            comp_seed = comp_seeds_base + sel_idx * 100 + ci
            trials.append({
                'f4_seed': seed_num,
                'f4_idx': arr_idx,
                'f4_gens': gens,
                'f4_info': info,
                'complement_seed': comp_seed,
            })

    n_trials = len(trials)
    config = vars(args)
    config['f4_seeds_used'] = [t['f4_seed'] for t in trials]
    config['comp_seeds_used'] = [t['complement_seed'] for t in trials]

    print(f"  F4 bases: {len(f4_selections)} "
          f"(seeds: {[s[1] for s in f4_selections]})")
    print(f"  Complement seeds per base: {args.n_comp}")
    print(f"  Total trials: {n_trials}")
    print(f"  Phase B: {args.b_s1_steps}+{args.b_s2_steps} steps")
    sched_str = " -> ".join(f"{f:.3f}" for f in args.unfreeze_schedule)
    print(f"  Unfreeze schedule: [{sched_str}]")
    print(f"  Phase C: {args.c_steps} steps")
    print(f"  Checkpoint: every {args.checkpoint_every} steps")
    print(f"  Device: {DEVICE}")
    if torch.cuda.is_available():
        print(f"  GPU: {torch.cuda.get_device_name()}")

    print()
    print("  F4 bases selected:")
    for arr_idx, seed_num, gens, info in f4_selections:
        print(f"    seed {seed_num:3d}: T_cl={info['T_cl']:.2e}  "
              f"S={info['simplicity']:.4f}  "
              f"closure={info['closure_frac']:.1%}  "
              f"K_std={info['K_std']:.4f}  "
              f"h_v(3)={info['h_dual_ir3']:.2f}")
    print()

    # ---- Resume logic ----
    all_results = []
    all_gens_list = []
    resume_trial_idx = None
    resume_phase = None
    resume_stage = None
    resume_step = None
    resume_ckpt = None

    if args.resume:
        ckpt = load_checkpoint(checkpoint_path, DEVICE)
        if ckpt is not None:
            all_results = ckpt['completed_trials']
            resume_trial_idx = ckpt['trial_idx']
            resume_phase = ckpt['phase']
            resume_stage = ckpt['stage']
            resume_step = ckpt['step']
            resume_ckpt = ckpt
            print(f"  Resumed: trial={resume_trial_idx}, "
                  f"phase={resume_phase}, stage={resume_stage}, "
                  f"step={resume_step}, completed={len(all_results)}")
        elif os.path.exists(save_path):
            data = np.load(save_path, allow_pickle=True)
            all_results = json.loads(str(data['results']))
            print(f"  Resumed from .npz: {len(all_results)} trials completed")

    start_trial = len(all_results)
    if resume_trial_idx is not None:
        start_trial = resume_trial_idx

    # Graceful interrupt handler
    interrupted = [False]

    def handler(sig, frame):
        if interrupted[0]:
            print("\n  Force quit.")
            sys.exit(1)
        interrupted[0] = True
        print("\n  Interrupt received - saving checkpoint...")

    signal.signal(signal.SIGINT, handler)

    # ============================================================
    # TRIAL LOOP
    # ============================================================
    t0_global = time.time()

    for trial_idx in range(start_trial, n_trials):
        if interrupted[0]:
            break

        trial = trials[trial_idx]
        f4_seed = trial['f4_seed']
        comp_seed = trial['complement_seed']
        f4_gens = trial['f4_gens']

        print(f"\n{'#' * 80}")
        print(f"  TRIAL {trial_idx + 1}/{n_trials}: "
              f"F4 seed={f4_seed}, complement={comp_seed}")
        print(f"  F4 quality: T_cl={trial['f4_info']['T_cl']:.2e}, "
              f"S={trial['f4_info']['simplicity']:.4f}")
        print(f"{'#' * 80}")

        t0_trial = time.time()

        r_phase = (resume_phase if trial_idx == resume_trial_idx else None)

        # ---- Phase B ----
        gens_78 = None
        if not interrupted[0]:
            r_b_stage = resume_stage if r_phase == 'B' else None
            r_b_step = resume_step if r_phase == 'B' else None
            r_b_f4 = (resume_ckpt.get('param_f4_params')
                      if r_phase == 'B' and resume_ckpt else None)
            r_b_comp = (resume_ckpt.get('param_complement_params')
                        if r_phase == 'B' and resume_ckpt else None)
            r_b_optim = (resume_ckpt.get('optimizer_state')
                         if r_phase == 'B' and resume_ckpt else None)

            gens_78, history = phase_b_grow(
                f4_gens,
                complement_seed=comp_seed,
                stage1_steps=args.b_s1_steps,
                stage1_gamma=args.b_s1_gamma,
                stage1_lr=args.b_s1_lr,
                stage2_steps=args.b_s2_steps,
                stage2_gamma=args.b_s2_gamma,
                stage2_lr=args.b_s2_lr,
                unfreeze_schedule=args.unfreeze_schedule,
                trial_idx=trial_idx,
                completed_trials=all_results,
                config=config,
                checkpoint_path=checkpoint_path,
                checkpoint_every=args.checkpoint_every,
                resume_stage=r_b_stage,
                resume_step=r_b_step,
                resume_f4_params=r_b_f4,
                resume_complement_params=r_b_comp,
                resume_optimizer_state=r_b_optim,
                interrupted=interrupted,
                verbose=True,
            )
            if gens_78 is None:
                break

        # ---- Phase C ----
        if not interrupted[0] and args.c_steps > 0:
            r_c_step = resume_step if r_phase == 'C' else None
            r_c_params = (resume_ckpt.get('param_params_78')
                          if r_phase == 'C' and resume_ckpt else None)
            r_c_optim = (resume_ckpt.get('optimizer_state')
                         if r_phase == 'C' and resume_ckpt else None)

            gens_78 = phase_c_anneal(
                gens_78,
                steps=args.c_steps,
                gamma=args.c_gamma,
                lr=args.c_lr,
                trial_idx=trial_idx,
                completed_trials=all_results,
                f4_gens_np=f4_gens,
                config=config,
                checkpoint_path=checkpoint_path,
                checkpoint_every=args.checkpoint_every,
                resume_step=r_c_step,
                resume_params=r_c_params,
                resume_optimizer_state=r_c_optim,
                interrupted=interrupted,
                verbose=True,
            )
            if gens_78 is None:
                break

        if interrupted[0]:
            break

        # ---- Analysis ----
        result = analyze_e6(
            gens_78,
            name=f"Trial {trial_idx + 1} (F4={f4_seed}, comp={comp_seed})")
        result['f4_seed'] = f4_seed
        result['complement_seed'] = comp_seed
        result['f4_band'] = trial['f4_info']['simplicity']
        result['f4_T_cl'] = trial['f4_info']['T_cl']
        all_results.append(result)
        all_gens_list.append(np.array(gens_78))

        dt_trial = time.time() - t0_trial
        print(f"\n  Trial time: {dt_trial:.0f}s ({dt_trial / 60:.1f} min)")

        # Save
        save_data = {
            'results': json.dumps(all_results),
            'n_trials_done': len(all_results),
            'config': json.dumps(config),
        }
        if all_gens_list:
            save_data['generators'] = np.array(
                [g if isinstance(g, np.ndarray) else np.array(g)
                 for g in all_gens_list])
        np.savez(save_path, **save_data)
        print(f"  Results saved ({len(all_results)}/{n_trials} trials)")

        if os.path.exists(checkpoint_path):
            os.remove(checkpoint_path)

        resume_trial_idx = None
        resume_phase = None
        resume_stage = None
        resume_step = None
        resume_ckpt = None

    # ============================================================
    # SUMMARY
    # ============================================================
    elapsed = time.time() - t0_global

    if interrupted[0]:
        print(f"\n  Interrupted after {len(all_results)}/{n_trials} trials. "
              f"Re-run with --resume to continue.")
        return

    if all_results:
        print(f"\n{'=' * 80}")
        print("  SUMMARY")
        print(f"{'=' * 80}")
        print(f"  {'F4':>4s}  {'Comp':>5s}  {'Closure':>8s}  "
              f"{'Simple':>8s}  {'h_v(3)':>8s}  {'Rank':>4s}  Verdict")
        print(f"  {'-' * 4}  {'-' * 5}  {'-' * 8}  "
              f"{'-' * 8}  {'-' * 8}  {'-' * 4}  {'-' * 20}")
        for r in all_results:
            print(f"  {r['f4_seed']:4d}  {r['complement_seed']:5d}  "
                  f"{r['closure_frac']:7.1%}  "
                  f"{r['simplicity']:8.4f}  {r['h_dual_ir3']:8.2f}  "
                  f"{r['rank_estimate']:4d}  {r['verdict']}")

        # Best result
        best = max(all_results, key=lambda r: r['closure_frac'])
        print(f"\n  Best: F4={best['f4_seed']}, "
              f"comp={best['complement_seed']}, "
              f"closure={best['closure_frac']:.1%}, "
              f"h_v(3)={best['h_dual_ir3']:.2f}")

    print(f"\n  Total elapsed: {elapsed:.0f}s ({elapsed / 60:.1f} min)")


if __name__ == '__main__':
    main()
