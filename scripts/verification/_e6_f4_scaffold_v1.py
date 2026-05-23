"""
E₆ blind discovery via F₄ scaffolding.

Strategy (motivated by §11e failure):
  E₆ direct-from-random fails because the optimizer converges to an F₄ subalgebra
  but cannot "jump" to the full 78-dim algebra. The scaffolding approach:

  Phase A: Converge F₄ (52 gens in so(27)) — proven to work reliably.
  Phase B: Add 26 complement generators. Optimize with two param groups:
           - F₄ generators: frozen (lr=0) or very low lr (1e-6)
           - Complement generators: normal lr (1e-3 → 1e-4)
           The complement must complete a 78-dim simple, closed algebra.
  Phase C: Unfreeze all 78 generators at low lr for final annealing.

  Expected signatures if E₆ emerges:
    - 78 generators, all commutators closed
    - Killing metric ∝ I₇₈  (simple)
    - h∨ = 12 (with I_R = 3 for 27-dim rep: h∨(I_R=3) = 36)
    - Cartan subalgebra: rank 6

Usage:
  # Full run with 3 seeds
  python _e6_f4_scaffold.py

  # Quick local test
  python _e6_f4_scaffold.py --quick

  # Resume after interruption
  python _e6_f4_scaffold.py --resume

  # Use existing F₄ checkpoint (skip Phase A)
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
# PHASE A: CONVERGE F₄
# ============================================================

def phase_a_converge_f4(seed, s1_steps=8000, s2_steps=3000,
                        s1_gamma=50.0, s2_gamma=200.0,
                        s1_lr=0.001, s2_lr=0.0001,
                        verbose=True):
    """Converge a single F₄ solution (52 gens in so(27)).

    Returns: (gens_np, params_torch) — numpy (52, 27, 27) and torch params.
    """
    from _f4_direct_anneal import run_direct_anneal

    if verbose:
        print(f"\n  Phase A: Converging F4 (seed={seed})")

    all_gens, history = run_direct_anneal(
        n_gen=N_GEN_F4, dim=DIM,
        stage1_steps=s1_steps, stage1_lr=s1_lr, stage1_gamma=s1_gamma,
        stage2_steps=s2_steps, stage2_lr=s2_lr, stage2_gamma=s2_gamma,
        batch_size=1,
        seeds=[seed],
        verbose=verbose,
    )

    # all_gens shape: (1, 52, 27, 27)
    gens_np = all_gens[0]  # (52, 27, 27)

    # Verify F₄ convergence
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
    generators = data['generators']  # (n_seeds, 52, 27, 27)

    # Find the index for the requested seed
    for idx, r in enumerate(results):
        if r['seed'] == seed and r['category'] == 'F4':
            gens_np = generators[idx]  # (52, 27, 27)
            print(f"  Loaded F4 from checkpoint: seed={seed}, "
                  f"h_v(3)={r['h_dual_ir3']:.2f}")
            return gens_np, r

    raise ValueError(f"No F4 solution found for seed={seed} in {checkpoint_path}")


# ============================================================
# PHASE B: SCAFFOLD — GROW COMPLEMENT
# ============================================================

def phase_b_scaffold(f4_gens_np, complement_seed=999,
                     stage1_steps=15000, stage1_gamma=50.0, stage1_lr=0.001,
                     stage2_steps=15000, stage2_gamma=200.0, stage2_lr=0.0001,
                     f4_lr_factor=0.0,
                     verbose=True):
    """Grow 26 complement generators around frozen F₄ scaffold.

    Args:
        f4_gens_np: (52, 27, 27) converged F₄ generators
        complement_seed: RNG seed for initializing the 26 complement generators
        f4_lr_factor: Learning rate multiplier for F₄ params.
            0.0 = fully frozen, 0.01 = very slow adaptation, 1.0 = equal to complement
        stage1_steps/stage2_steps: Optimizer steps per stage
        stage1_gamma/stage2_gamma: Closure penalty weight per stage

    Returns: (gens_78_np, history)
        gens_78_np: (78, 27, 27) numpy complex128
    """
    device = DEVICE

    # Convert F₄ generators to parameter space
    f4_torch = torch.tensor(
        np.array(f4_gens_np), dtype=DTYPE, device=device
    ).unsqueeze(0)  # (1, 52, 27, 27)
    f4_params = antisymmetric_to_params(f4_torch)  # (1, 52, n_p)
    f4_params = f4_params.squeeze(0)  # (52, n_p)

    # Initialize 26 complement generators randomly
    n_p = DIM * (DIM - 1) // 2
    np.random.seed(complement_seed)
    complement_init = np.random.randn(N_GEN_COMPLEMENT, n_p) * 0.1
    complement_params = torch.tensor(
        complement_init, dtype=FDTYPE, device=device
    )

    # Separate parameter groups
    f4_params = f4_params.clone().detach().requires_grad_(True)
    complement_params = complement_params.clone().detach().requires_grad_(True)

    K_target = -1.0 * torch.eye(N_GEN_E6, dtype=FDTYPE, device=device)

    history = {"T_cl": [], "K_diag_std": [], "stage": []}

    for stage, (steps, lr, gamma) in enumerate([
        (stage1_steps, stage1_lr, stage1_gamma),
        (stage2_steps, stage2_lr, stage2_gamma),
    ], 1):
        if verbose:
            f4_lr_effective = lr * f4_lr_factor
            print(f"\n  Phase B Stage {stage}: gamma={gamma}, "
                  f"complement_lr={lr}, f4_lr={f4_lr_effective:.1e}, "
                  f"{steps} steps")

        param_groups = [
            {'params': [complement_params], 'lr': lr},
            {'params': [f4_params], 'lr': lr * f4_lr_factor},
        ]
        optimizer = torch.optim.Adam(param_groups)

        for step in range(steps):
            optimizer.zero_grad()

            # Assemble full 78-gen parameter tensor: (1, 78, n_p)
            all_params = torch.cat([f4_params, complement_params], dim=0)
            all_params = all_params.unsqueeze(0)  # (1, 78, n_p)

            generators = params_to_antisymmetric(all_params, N_GEN_E6, dim=DIM)

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
                    kd = torch.diagonal(K[0]).mean().item()
                    kd_std = torch.diagonal(K[0]).std().item()
                    # Complement norm (are they growing?)
                    gens_all = params_to_antisymmetric(
                        all_params, N_GEN_E6, dim=DIM)
                    comp_norms = (gens_all[0, N_GEN_F4:].abs() ** 2
                                  ).sum(dim=(-2, -1)).sqrt().mean().item()
                    print(f"    step {step:5d}: T_cl={tcl:.4e} "
                          f"K_diag={kd:.4f}±{kd_std:.4e} "
                          f"comp_norm={comp_norms:.4f}")
                    history["T_cl"].append(tcl)
                    history["K_diag_std"].append(kd_std)
                    history["stage"].append(stage)

    # Extract final generators
    with torch.no_grad():
        all_params = torch.cat([f4_params, complement_params], dim=0)
        all_params = all_params.unsqueeze(0)
        final_gens = params_to_antisymmetric(all_params, N_GEN_E6, dim=DIM)
        gens_78_np = final_gens[0].cpu().numpy()

    return gens_78_np, history


# ============================================================
# PHASE C: FULL ANNEALING (OPTIONAL)
# ============================================================

def phase_c_anneal(gens_78_np, steps=10000, gamma=200.0, lr=5e-5,
                   verbose=True):
    """Final joint annealing of all 78 generators at low lr.

    This allows small adjustments to break any artifacts from the
    frozen/unfrozen boundary.
    """
    device = DEVICE

    if verbose:
        print(f"\n  Phase C: Joint annealing ({steps} steps, lr={lr}, gamma={gamma})")

    gens_torch = torch.tensor(
        np.array(gens_78_np), dtype=DTYPE, device=device
    ).unsqueeze(0)  # (1, 78, 27, 27)
    params_78 = antisymmetric_to_params(gens_torch)
    params_78 = params_78.clone().detach().requires_grad_(True)

    K_target = -1.0 * torch.eye(N_GEN_E6, dtype=FDTYPE, device=device)
    optimizer = torch.optim.Adam([params_78], lr=lr)

    for step in range(steps):
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
                print(f"    step {step:5d}: T_cl={tcl:.4e} K_diag_std={kd_std:.4e}")

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
    # Diagonalize ad_H for random H in the algebra to find maximal torus
    structure_constants = coeffs  # (78, n_pairs) — structure constants

    # Simple rank estimate: count near-zero eigenvalues of ad_H
    # for a generic element H
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

    # Verdict
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
# CHECKPOINT
# ============================================================

def save_checkpoint(path, phase, seed, data):
    """Save checkpoint with phase marker."""
    kw = dict(phase=phase, seed=seed)
    kw.update(data)
    np.savez(path, **kw)
    print(f"  Checkpoint saved: phase={phase}, seed={seed}")


def load_checkpoint(path):
    """Load checkpoint. Returns (phase, seed, data_dict) or None."""
    if not os.path.exists(path):
        return None
    data = np.load(path, allow_pickle=True)
    phase = str(data['phase'])
    seed = int(data['seed'])
    return phase, seed, data


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='E6 blind discovery via F4 scaffolding')
    parser.add_argument('--seeds', type=int, nargs='+', default=[42, 137, 256],
                        help='Seeds for F4 Phase A (default: 42 137 256)')
    parser.add_argument('--complement-seeds', type=int, nargs='+', default=None,
                        help='Seeds for complement init (default: seed+1000)')
    # Phase A
    parser.add_argument('--a-s1-steps', type=int, default=8000)
    parser.add_argument('--a-s2-steps', type=int, default=3000)
    # Phase B
    parser.add_argument('--b-s1-steps', type=int, default=15000,
                        help='Phase B Stage 1 steps (complement growth)')
    parser.add_argument('--b-s2-steps', type=int, default=15000,
                        help='Phase B Stage 2 steps (complement annealing)')
    parser.add_argument('--b-s1-gamma', type=float, default=50.0)
    parser.add_argument('--b-s2-gamma', type=float, default=200.0)
    parser.add_argument('--b-s1-lr', type=float, default=0.001)
    parser.add_argument('--b-s2-lr', type=float, default=0.0001)
    parser.add_argument('--f4-lr-factor', type=float, default=0.0,
                        help='LR multiplier for F4 params in Phase B. '
                             '0.0=frozen, 0.01=slow drift (default: 0.0)')
    # Phase C
    parser.add_argument('--c-steps', type=int, default=10000,
                        help='Phase C annealing steps (0 to skip)')
    parser.add_argument('--c-lr', type=float, default=5e-5)
    parser.add_argument('--c-gamma', type=float, default=200.0)
    # External F₄ checkpoint
    parser.add_argument('--f4-checkpoint', type=str, default=None,
                        help='Load F4 from existing landscape checkpoint (.npz)')
    parser.add_argument('--f4-seed', type=int, default=None,
                        help='Which seed to load from F4 checkpoint')
    # General
    parser.add_argument('--save', type=str, default='e6_scaffold_results.npz')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--quick', action='store_true',
                        help='Quick test: reduced steps')
    args = parser.parse_args()

    if args.quick:
        args.a_s1_steps = 2000
        args.a_s2_steps = 1000
        args.b_s1_steps = 3000
        args.b_s2_steps = 2000
        args.c_steps = 2000

    if args.complement_seeds is None:
        args.complement_seeds = [s + 1000 for s in args.seeds]

    save_dir = os.path.dirname(__file__)
    save_path = os.path.join(save_dir, args.save)

    print("=" * 80)
    print("  E₆ BLIND DISCOVERY VIA F₄ SCAFFOLDING")
    print("=" * 80)
    print(f"  Seeds:          {args.seeds}")
    print(f"  Complement:     {args.complement_seeds}")
    print(f"  Phase A (F4):   {args.a_s1_steps}+{args.a_s2_steps} steps")
    print(f"  Phase B (grow): {args.b_s1_steps}+{args.b_s2_steps} steps, "
          f"f4_lr_factor={args.f4_lr_factor}")
    print(f"  Phase C (ann):  {args.c_steps} steps" if args.c_steps > 0
          else "  Phase C: SKIPPED")
    print(f"  Device:         {DEVICE}")
    if torch.cuda.is_available():
        print(f"  GPU:            {torch.cuda.get_device_name()}")
    print()

    all_results = []
    t0_global = time.time()

    for trial_idx, (seed, comp_seed) in enumerate(
            zip(args.seeds, args.complement_seeds)):
        print(f"\n{'#' * 80}")
        print(f"  TRIAL {trial_idx + 1}/{len(args.seeds)}: "
              f"seed={seed}, complement_seed={comp_seed}")
        print(f"{'#' * 80}")

        t0_trial = time.time()

        # ---- Phase A: Get F₄ ----
        if args.f4_checkpoint is not None:
            f4_seed = args.f4_seed if args.f4_seed is not None else seed
            f4_gens, f4_info = load_f4_from_checkpoint(
                os.path.join(save_dir, args.f4_checkpoint), f4_seed)
        else:
            f4_gens, f4_info = phase_a_converge_f4(
                seed,
                s1_steps=args.a_s1_steps,
                s2_steps=args.a_s2_steps,
                verbose=True,
            )

        print(f"  Phase A complete: {f4_info.get('category', '?')}")

        # ---- Phase B: Grow complement ----
        print(f"\n  Phase B: Growing 26 complement generators...")
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
            verbose=True,
        )

        # ---- Phase C: Joint annealing ----
        if args.c_steps > 0:
            gens_78 = phase_c_anneal(
                gens_78,
                steps=args.c_steps,
                gamma=args.c_gamma,
                lr=args.c_lr,
                verbose=True,
            )

        # ---- Analysis ----
        gens_list = list(gens_78)
        result = analyze_e6(gens_list,
                            name=f"Trial {trial_idx + 1} "
                                 f"(seed={seed}, comp={comp_seed})")
        result['seed'] = seed
        result['complement_seed'] = comp_seed
        result['f4_category'] = f4_info.get('category', '?')
        all_results.append(result)

        dt_trial = time.time() - t0_trial
        print(f"\n  Trial time: {dt_trial:.0f}s ({dt_trial / 60:.1f} min)")

        # Save checkpoint after each trial
        save_data = {
            'results': json.dumps(all_results),
            'generators': np.array(gens_list),
            'history_T_cl': np.array(history['T_cl']),
        }
        save_checkpoint(save_path, f'done_trial_{trial_idx + 1}',
                        seed, save_data)

    # ============================================================
    # SUMMARY
    # ============================================================
    elapsed = time.time() - t0_global
    print(f"\n{'=' * 80}")
    print("  SUMMARY")
    print(f"{'=' * 80}")
    print(f"  {'Seed':>4s}  {'Comp':>4s}  {'F4':>4s}  {'Closure':>8s}  "
          f"{'Simple':>8s}  {'h_v(3)':>8s}  {'Rank':>4s}  {'Verdict'}")
    print(f"  {'-' * 4}  {'-' * 4}  {'-' * 4}  {'-' * 8}  "
          f"{'-' * 8}  {'-' * 8}  {'-' * 4}  {'-' * 20}")
    for r in all_results:
        print(f"  {r['seed']:4d}  {r['complement_seed']:4d}  "
              f"{r['f4_category']:>4s}  {r['closure_frac']:7.1%}  "
              f"{r['simplicity']:8.4f}  {r['h_dual_ir3']:8.2f}  "
              f"{r['rank_estimate']:4d}  {r['verdict']}")

    print(f"\n  Total elapsed: {elapsed:.0f}s ({elapsed / 60:.1f} min)")


if __name__ == "__main__":
    main()
