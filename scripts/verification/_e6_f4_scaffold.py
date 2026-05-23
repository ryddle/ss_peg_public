"""
E₆ blind discovery via F₄ scaffolding — Resumable edition.

Strategy (motivated by §11e failure):
  E₆ direct-from-random fails because the optimizer converges to an F₄ subalgebra
  but cannot "jump" to the full 78-dim algebra. The scaffolding approach:

  Phase A: Converge F₄ (52 gens in so(27)) — proven to work reliably.
  Phase B: Add 26 complement generators. Optimize with two param groups:
           - F₄ generators: frozen (lr=0) or very low lr (1e-6)
           - Complement generators: normal lr (1e-3 → 1e-4)
           The complement must complete a 78-dim simple, closed algebra.
  Phase C: Unfreeze all 78 generators at low lr for final annealing.

Fine-grained checkpointing every N steps allows interruption and resume
from the exact point where it stopped — no wasted computation.

Checkpoint format (e6_*_checkpoint.pt):
  - trial_idx, phase (A/B/C), stage (1/2), step
  - params (+ optimizer state)
  - completed_trials results
  - f4_gens (for phase B/C restore)

Usage:
  python _e6_f4_scaffold.py                    # Full run with 3 seeds
  python _e6_f4_scaffold.py --quick            # Quick local test
  python _e6_f4_scaffold.py --resume           # Resume from checkpoint
  python _e6_f4_scaffold.py --f4-checkpoint landscape_128seeds_results.npz --f4-seed 0
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
    random_params_so,
    killing_metric_np,
    commutator_np,
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
# STEP-LEVEL CHECKPOINT I/O
# ============================================================

def save_step_checkpoint(path, trial_idx, phase, stage, step,
                         params_dict, optimizer, completed_trials,
                         f4_gens=None, config=None):
    """Save fine-grained checkpoint including optimizer state.

    Args:
        params_dict: dict of tensors. For Phase A: {'params': T}.
                     For Phase B: {'f4_params': T, 'complement_params': T}.
                     For Phase C: {'params_78': T}.
    """
    data = {
        'trial_idx': trial_idx,
        'phase': phase,
        'stage': stage,
        'step': step,
        'optimizer_state': optimizer.state_dict(),
        'completed_trials': json.dumps(completed_trials),
        'config': json.dumps(config) if config else '{}',
    }
    # Store each param tensor
    for k, v in params_dict.items():
        data[f'param_{k}'] = v.detach().cpu()

    # Store F4 generators for Phase B/C resume
    if f4_gens is not None:
        data['f4_gens'] = np.array(f4_gens) if not isinstance(f4_gens, np.ndarray) else f4_gens

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
    data['completed_trials'] = json.loads(data['completed_trials'])
    data['config'] = json.loads(data['config'])
    # Move param tensors to device
    for k in list(data.keys()):
        if k.startswith('param_') and isinstance(data[k], torch.Tensor):
            data[k] = data[k].to(device)
    return data


# ============================================================
# PHASE A: CONVERGE F₄ (inline with checkpointing)
# ============================================================

def phase_a_converge_f4(seed, s1_steps=8000, s2_steps=3000,
                        s1_gamma=50.0, s2_gamma=200.0,
                        s1_lr=0.001, s2_lr=0.0001,
                        trial_idx=0, completed_trials=None,
                        config=None, checkpoint_path=None,
                        checkpoint_every=1000,
                        resume_stage=None, resume_step=None,
                        resume_params=None, resume_optimizer_state=None,
                        interrupted=None, verbose=True):
    """Converge a single F₄ solution (52 gens in so(27)) with checkpointing.

    Returns: (gens_np, info_dict) or (None, None) on interrupt.
    """
    device = DEVICE
    if completed_trials is None:
        completed_trials = []

    if verbose:
        print(f"\n  Phase A: Converging F4 (seed={seed})")

    # Initialize or restore params
    if resume_params is not None:
        params = resume_params.clone().detach().requires_grad_(True)
    else:
        params = random_params_so(1, N_GEN_F4, DIM, device, [seed])
        with torch.no_grad():
            gens = params_to_antisymmetric(params, N_GEN_F4, dim=DIM)
            norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
            gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
            params = antisymmetric_to_params(gens)
        params = params.clone().detach().requires_grad_(True)

    K_target = -1.0 * torch.eye(N_GEN_F4, dtype=FDTYPE, device=device)

    stages = [
        (1, s1_steps, s1_lr, s1_gamma),
        (2, s2_steps, s2_lr, s2_gamma),
    ]

    for stage_num, steps, lr, gamma in stages:
        if resume_stage is not None and stage_num < resume_stage:
            continue

        start_step = 0
        if resume_stage is not None and stage_num == resume_stage and resume_step is not None:
            start_step = resume_step

        if verbose:
            print(f"  Stage {stage_num}: gamma={gamma}, lr={lr}, {steps} steps"
                  + (f" (resuming from step {start_step})" if start_step > 0 else ""))

        optimizer = torch.optim.Adam([params], lr=lr)
        if (resume_optimizer_state is not None
                and resume_stage is not None
                and stage_num == resume_stage):
            try:
                optimizer.load_state_dict(resume_optimizer_state)
            except Exception as e:
                print(f"  WARNING: Could not restore optimizer state: {e}")
            resume_optimizer_state = None

        for step in range(start_step, steps):
            if interrupted is not None and interrupted[0]:
                save_step_checkpoint(
                    checkpoint_path, trial_idx, 'A', stage_num, step,
                    {'params': params}, optimizer, completed_trials,
                    config=config)
                print(f"\n  Checkpoint saved: phase=A, stage={stage_num}, step={step}")
                return None, None

            optimizer.zero_grad()
            generators = params_to_antisymmetric(params, N_GEN_F4, dim=DIM)
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
                save_step_checkpoint(
                    checkpoint_path, trial_idx, 'A', stage_num, step + 1,
                    {'params': params}, optimizer, completed_trials,
                    config=config)
                if verbose and (step + 1) % 2000 != 0:
                    print(f"    [checkpoint at stage {stage_num}, step {step + 1}]",
                          flush=True)

        resume_stage = None
        resume_step = None

    # Classify result
    with torch.no_grad():
        final_gens = params_to_antisymmetric(params, N_GEN_F4, dim=DIM)
        gens_np = final_gens[0].cpu().numpy()

    from _landscape_montecarlo import classify_solution
    info = classify_solution(list(gens_np))

    if verbose:
        print(f"  Phase A result: {info['category']}, "
              f"closure={info['closure_frac']:.1%}, "
              f"h_v(3)={info['h_dual_ir3']:.2f}")

    if info['category'] != 'F4':
        print(f"  WARNING: Phase A did not converge to F4 "
              f"(got {info['category']}). Proceeding anyway.")

    return gens_np, info


def load_f4_from_checkpoint(checkpoint_path, seed):
    """Load a converged F₄ solution from an existing landscape checkpoint."""
    data = np.load(checkpoint_path, allow_pickle=True)
    results = json.loads(str(data['results']))
    generators = data['generators']

    for idx, r in enumerate(results):
        if r['seed'] == seed and r['category'] == 'F4':
            gens_np = generators[idx]
            print(f"  Loaded F4 from checkpoint: seed={seed}, "
                  f"h_v(3)={r['h_dual_ir3']:.2f}")
            return gens_np, r

    raise ValueError(f"No F4 solution found for seed={seed} in {checkpoint_path}")


# ============================================================
# PHASE B: SCAFFOLD — GROW COMPLEMENT (with checkpointing)
# ============================================================

def phase_b_scaffold(f4_gens_np, complement_seed=999,
                     stage1_steps=15000, stage1_gamma=50.0, stage1_lr=0.001,
                     stage2_steps=15000, stage2_gamma=200.0, stage2_lr=0.0001,
                     f4_lr_factor=0.0,
                     trial_idx=0, completed_trials=None,
                     config=None, checkpoint_path=None,
                     checkpoint_every=1000,
                     resume_stage=None, resume_step=None,
                     resume_f4_params=None, resume_complement_params=None,
                     resume_optimizer_state=None,
                     interrupted=None, verbose=True):
    """Grow 26 complement generators around frozen F₄ scaffold, with checkpointing.

    Returns: (gens_78_np, history) or (None, None) on interrupt.
    """
    device = DEVICE
    if completed_trials is None:
        completed_trials = []

    # Convert F₄ generators to parameter space
    if resume_f4_params is not None:
        f4_params = resume_f4_params.clone().detach().requires_grad_(True)
    else:
        f4_torch = torch.tensor(
            np.array(f4_gens_np), dtype=DTYPE, device=device
        ).unsqueeze(0)
        f4_params = antisymmetric_to_params(f4_torch).squeeze(0)
        f4_params = f4_params.clone().detach().requires_grad_(True)

    # Initialize or restore complement generators
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
    history = {"T_cl": [], "K_diag_std": [], "stage": []}

    stages = [
        (1, stage1_steps, stage1_lr, stage1_gamma),
        (2, stage2_steps, stage2_lr, stage2_gamma),
    ]

    for stage_num, steps, lr, gamma in stages:
        if resume_stage is not None and stage_num < resume_stage:
            continue

        start_step = 0
        if resume_stage is not None and stage_num == resume_stage and resume_step is not None:
            start_step = resume_step

        f4_lr_effective = lr * f4_lr_factor
        if verbose:
            print(f"\n  Phase B Stage {stage_num}: gamma={gamma}, "
                  f"complement_lr={lr}, f4_lr={f4_lr_effective:.1e}, "
                  f"{steps} steps"
                  + (f" (resuming from step {start_step})" if start_step > 0 else ""))

        param_groups = [
            {'params': [complement_params], 'lr': lr},
            {'params': [f4_params], 'lr': f4_lr_effective},
        ]
        optimizer = torch.optim.Adam(param_groups)

        if (resume_optimizer_state is not None
                and resume_stage is not None
                and stage_num == resume_stage):
            try:
                optimizer.load_state_dict(resume_optimizer_state)
            except Exception as e:
                print(f"  WARNING: Could not restore optimizer state: {e}")
            resume_optimizer_state = None

        for step in range(start_step, steps):
            if interrupted is not None and interrupted[0]:
                save_step_checkpoint(
                    checkpoint_path, trial_idx, 'B', stage_num, step,
                    {'f4_params': f4_params, 'complement_params': complement_params},
                    optimizer, completed_trials,
                    f4_gens=f4_gens_np, config=config)
                print(f"\n  Checkpoint saved: phase=B, stage={stage_num}, step={step}")
                return None, None

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
                          f"comp_norm={comp_norms:.4f}",
                          flush=True)
                    history["T_cl"].append(tcl)
                    history["K_diag_std"].append(kd_std)
                    history["stage"].append(stage_num)

            if checkpoint_path and (step + 1) % checkpoint_every == 0:
                save_step_checkpoint(
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

    return gens_78_np, history


# ============================================================
# PHASE C: FULL ANNEALING (with checkpointing)
# ============================================================

def phase_c_anneal(gens_78_np, steps=10000, gamma=200.0, lr=5e-5,
                   trial_idx=0, completed_trials=None,
                   f4_gens_np=None, config=None,
                   checkpoint_path=None, checkpoint_every=1000,
                   resume_stage=None, resume_step=None,
                   resume_params=None, resume_optimizer_state=None,
                   interrupted=None, verbose=True):
    """Final joint annealing of all 78 generators at low lr, with checkpointing.

    Returns: gens_78_np or None on interrupt.
    """
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

    start_step = 0
    if resume_step is not None:
        start_step = resume_step
        if verbose:
            print(f"  (resuming from step {start_step})")

    optimizer = torch.optim.Adam([params_78], lr=lr)
    if resume_optimizer_state is not None:
        try:
            optimizer.load_state_dict(resume_optimizer_state)
        except Exception as e:
            print(f"  WARNING: Could not restore optimizer state: {e}")

    for step in range(start_step, steps):
        if interrupted is not None and interrupted[0]:
            save_step_checkpoint(
                checkpoint_path, trial_idx, 'C', 1, step,
                {'params_78': params_78}, optimizer, completed_trials,
                f4_gens=f4_gens_np, config=config)
            print(f"\n  Checkpoint saved: phase=C, step={step}")
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
            save_step_checkpoint(
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
    """Detailed E₆ analysis — closure, Killing, h∨, rank."""
    gens_78_np = np.asarray(gens_78_np)
    n = len(gens_78_np)
    total_pairs = n * (n - 1) // 2

    # Vectorized closure
    idx_i, idx_j = np.triu_indices(n, k=1)
    G_i = gens_78_np[idx_i]
    G_j = gens_78_np[idx_j]
    C = G_i @ G_j - G_j @ G_i
    C_flat = C.reshape(-1, gens_78_np[0].shape[0] ** 2).T
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
    h_raw = float(abs(np.diag(K_r).real.mean() / Tr_sq)) if abs(Tr_sq) > 1e-12 else 0.0
    h_dual_ir3 = 3.0 * h_raw

    # Rank estimation via Cartan subalgebra
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
    print(f"  h∨ (I_R=1): {h_raw:.4f}")
    print(f"  h∨ (I_R=3): {h_dual_ir3:.4f}  (E6 expected: 36.0, F4: 9.0)")
    print(f"  Estimated rank: {rank_estimate}  (E6: 6, F4: 4)")

    is_closed = closure_frac >= 0.99
    is_simple = simplicity < 0.05
    is_e6 = is_closed and is_simple and 30.0 < h_dual_ir3 < 42.0 and rank_estimate >= 5
    is_f4 = is_closed and is_simple and 7.0 < h_dual_ir3 < 11.0

    if is_e6:
        verdict = "E6 (PASS)"
    elif is_f4:
        verdict = "F4 embedded in 78 gens (partial success)"
    elif is_closed and is_simple:
        verdict = f"SIMPLE (h∨={h_dual_ir3:.1f}, rank≈{rank_estimate})"
    elif is_closed:
        verdict = "CLOSED_NON_SIMPLE"
    else:
        verdict = f"INCOMPLETE (closure={closure_frac:.1%})"

    print(f"  VERDICT: {verdict}")

    return {
        'closure_frac': closure_frac,
        'max_residual': max_res,
        'K_mean': K_mean,
        'K_std': K_std,
        'simplicity': simplicity,
        'h_raw': h_raw,
        'h_dual_ir3': h_dual_ir3,
        'rank_estimate': rank_estimate,
        'verdict': verdict,
    }


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='E6 blind discovery via F4 scaffolding [Resumable]')
    parser.add_argument('--seeds', type=int, nargs='+', default=[42, 137, 256],
                        help='Seeds for F4 Phase A (default: 42 137 256)')
    parser.add_argument('--complement-seeds', type=int, nargs='+', default=None,
                        help='Seeds for complement init (default: seed+1000)')
    # Phase A
    parser.add_argument('--a-s1-steps', type=int, default=8000)
    parser.add_argument('--a-s2-steps', type=int, default=3000)
    # Phase B
    parser.add_argument('--b-s1-steps', type=int, default=15000)
    parser.add_argument('--b-s2-steps', type=int, default=15000)
    parser.add_argument('--b-s1-gamma', type=float, default=50.0)
    parser.add_argument('--b-s2-gamma', type=float, default=200.0)
    parser.add_argument('--b-s1-lr', type=float, default=0.001)
    parser.add_argument('--b-s2-lr', type=float, default=0.0001)
    parser.add_argument('--f4-lr-factor', type=float, default=0.0)
    # Phase C
    parser.add_argument('--c-steps', type=int, default=10000)
    parser.add_argument('--c-lr', type=float, default=5e-5)
    parser.add_argument('--c-gamma', type=float, default=200.0)
    # Checkpointing
    parser.add_argument('--checkpoint-every', type=int, default=1000,
                        help='Save checkpoint every N steps (default: 1000)')
    # External F₄
    parser.add_argument('--f4-checkpoint', type=str, default=None)
    parser.add_argument('--f4-seed', type=int, default=None)
    # General
    parser.add_argument('--save', type=str, default='e6_scaffold_results.npz')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--quick', action='store_true')
    args = parser.parse_args()

    if args.quick:
        args.seeds = [42]
        args.a_s1_steps = 2000
        args.a_s2_steps = 1000
        args.b_s1_steps = 3000
        args.b_s2_steps = 2000
        args.c_steps = 2000
        args.checkpoint_every = 500

    if args.complement_seeds is None:
        args.complement_seeds = [s + 1000 for s in args.seeds]

    save_dir = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(save_dir, args.save)
    checkpoint_path = os.path.join(save_dir,
                                   args.save.replace('.npz', '_checkpoint.pt'))

    config = vars(args)

    print("=" * 80)
    print("  E₆ BLIND DISCOVERY VIA F₄ SCAFFOLDING [Resumable]")
    print("=" * 80)
    print(f"  Seeds:          {args.seeds}")
    print(f"  Complement:     {args.complement_seeds}")
    print(f"  Phase A (F4):   {args.a_s1_steps}+{args.a_s2_steps} steps")
    print(f"  Phase B (grow): {args.b_s1_steps}+{args.b_s2_steps} steps, "
          f"f4_lr_factor={args.f4_lr_factor}")
    print(f"  Phase C (ann):  {args.c_steps} steps" if args.c_steps > 0
          else "  Phase C: SKIPPED")
    print(f"  Checkpoint:     every {args.checkpoint_every} steps")
    print(f"  Device:         {DEVICE}")
    if torch.cuda.is_available():
        print(f"  GPU:            {torch.cuda.get_device_name()}")
    print()

    # --- Resume logic ---
    all_results = []
    all_gens_list = []
    resume_trial_idx = None
    resume_phase = None
    resume_stage = None
    resume_step = None
    resume_ckpt = None  # full checkpoint data

    if args.resume:
        ckpt = load_step_checkpoint(checkpoint_path, DEVICE)
        if ckpt is not None:
            all_results = ckpt['completed_trials']
            resume_trial_idx = ckpt['trial_idx']
            resume_phase = ckpt['phase']
            resume_stage = ckpt['stage']
            resume_step = ckpt['step']
            resume_ckpt = ckpt
            print(f"  Resumed: trial={resume_trial_idx}, phase={resume_phase}, "
                  f"stage={resume_stage}, step={resume_step}, "
                  f"completed_trials={len(all_results)}")
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
        print("\n  Interrupt received — will save checkpoint at next opportunity...")

    signal.signal(signal.SIGINT, handler)

    # ============================================================
    # TRIAL LOOP
    # ============================================================
    t0_global = time.time()

    for trial_idx in range(start_trial, len(args.seeds)):
        if interrupted[0]:
            break

        seed = args.seeds[trial_idx]
        comp_seed = args.complement_seeds[trial_idx]

        print(f"\n{'#' * 80}")
        print(f"  TRIAL {trial_idx + 1}/{len(args.seeds)}: "
              f"seed={seed}, complement_seed={comp_seed}")
        print(f"{'#' * 80}")

        t0_trial = time.time()

        # Determine what to resume within this trial
        r_phase = resume_phase if trial_idx == resume_trial_idx else None
        f4_gens = None

        # ---- Phase A: Get F₄ ----
        if r_phase is not None and r_phase > 'A':
            # Phase A already done, load F4 from checkpoint
            if resume_ckpt is not None and 'f4_gens' in resume_ckpt:
                f4_gens = resume_ckpt['f4_gens']
                print(f"  Phase A: loaded from checkpoint")
            else:
                print(f"  WARNING: Phase A not in checkpoint, re-running")
                r_phase = 'A'

        if f4_gens is None:
            if args.f4_checkpoint is not None:
                f4_seed = args.f4_seed if args.f4_seed is not None else seed
                f4_gens, f4_info = load_f4_from_checkpoint(
                    os.path.join(save_dir, args.f4_checkpoint), f4_seed)
            else:
                # Resume params for Phase A
                r_a_stage = resume_stage if r_phase == 'A' else None
                r_a_step = resume_step if r_phase == 'A' else None
                r_a_params = (resume_ckpt.get('param_params')
                              if r_phase == 'A' and resume_ckpt else None)
                r_a_optim = (resume_ckpt.get('optimizer_state')
                             if r_phase == 'A' and resume_ckpt else None)

                f4_gens, f4_info = phase_a_converge_f4(
                    seed,
                    s1_steps=args.a_s1_steps,
                    s2_steps=args.a_s2_steps,
                    trial_idx=trial_idx,
                    completed_trials=all_results,
                    config=config,
                    checkpoint_path=checkpoint_path,
                    checkpoint_every=args.checkpoint_every,
                    resume_stage=r_a_stage,
                    resume_step=r_a_step,
                    resume_params=r_a_params,
                    resume_optimizer_state=r_a_optim,
                    interrupted=interrupted,
                    verbose=True,
                )

                if f4_gens is None:
                    break  # interrupted

            if f4_gens is not None and not isinstance(f4_info, dict):
                f4_info = {}
            print(f"  Phase A complete: {f4_info.get('category', '?')}")

        # ---- Phase B: Grow complement ----
        if not interrupted[0]:
            print(f"\n  Phase B: Growing 26 complement generators...")

            r_b_stage = resume_stage if r_phase == 'B' else None
            r_b_step = resume_step if r_phase == 'B' else None
            r_b_f4 = (resume_ckpt.get('param_f4_params')
                       if r_phase == 'B' and resume_ckpt else None)
            r_b_comp = (resume_ckpt.get('param_complement_params')
                        if r_phase == 'B' and resume_ckpt else None)
            r_b_optim = (resume_ckpt.get('optimizer_state')
                         if r_phase == 'B' and resume_ckpt else None)

            gens_78, history = phase_b_scaffold(
                f4_gens,
                complement_seed=comp_seed,
                stage1_steps=args.b_s1_steps,
                stage1_gamma=args.b_s1_gamma,
                stage1_lr=args.b_s1_lr,
                stage2_steps=args.b_s2_steps,
                stage2_gamma=args.b_s2_gamma,
                stage2_lr=args.b_s2_lr,
                f4_lr_factor=args.f4_lr_factor,
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
                break  # interrupted

        # ---- Phase C: Joint annealing ----
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
                break  # interrupted

        if interrupted[0]:
            break

        # ---- Analysis ----
        gens_list = list(gens_78)
        result = analyze_e6(gens_list,
                            name=f"Trial {trial_idx + 1} "
                                 f"(seed={seed}, comp={comp_seed})")
        result['seed'] = seed
        result['complement_seed'] = comp_seed
        result['f4_category'] = f4_info.get('category', '?') if 'f4_info' in dir() else '?'
        all_results.append(result)
        all_gens_list.append(np.array(gens_list))

        dt_trial = time.time() - t0_trial
        print(f"\n  Trial time: {dt_trial:.0f}s ({dt_trial / 60:.1f} min)")

        # Save batch-level results
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
        print(f"  Results saved ({len(all_results)}/{len(args.seeds)} trials)")

        # Remove step checkpoint (trial complete)
        if os.path.exists(checkpoint_path):
            os.remove(checkpoint_path)

        # Clear resume state
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
        print(f"\n  Interrupted after {len(all_results)}/{len(args.seeds)} trials. "
              f"Re-run with --resume to continue.")
        return

    if all_results:
        print(f"\n{'=' * 80}")
        print("  SUMMARY")
        print(f"{'=' * 80}")
        print(f"  {'Seed':>4s}  {'Comp':>4s}  {'Closure':>8s}  "
              f"{'Simple':>8s}  {'h_v(3)':>8s}  {'Rank':>4s}  {'Verdict'}")
        print(f"  {'-' * 4}  {'-' * 4}  {'-' * 8}  "
              f"{'-' * 8}  {'-' * 8}  {'-' * 4}  {'-' * 20}")
        for r in all_results:
            print(f"  {r['seed']:4d}  {r['complement_seed']:4d}  "
                  f"{r['closure_frac']:7.1%}  "
                  f"{r['simplicity']:8.4f}  {r['h_dual_ir3']:8.2f}  "
                  f"{r['rank_estimate']:4d}  {r['verdict']}")

    # Clean up
    if os.path.exists(checkpoint_path):
        os.remove(checkpoint_path)

    print(f"\n  Total elapsed: {elapsed:.0f}s ({elapsed / 60:.1f} min)")


if __name__ == "__main__":
    main()
