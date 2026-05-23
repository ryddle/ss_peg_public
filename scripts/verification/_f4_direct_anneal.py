"""
F₄ blind discovery: direct + learning-rate annealing.

Lessons from previous runs:
  - Direct gamma=50, 30k:   K∝I perfect (std=4e-6), T_cl=0.0002, max_res=4e-4
  - Curriculum gamma→200:    seed 42 found CLOSED non-simple (not F₄)

Strategy: 
  Stage 1: gamma=50, lr=0.001, 30k steps  → establishes K∝I + partial closure
  Stage 2: gamma=200, lr=0.0001, 30k steps → anneals into exact closure

The direct approach maintains Killing from step 0 (avoiding curriculum trap),
while stage 2 adds aggressive closure plus fine lr for breaking through
the residual barrier.
"""
import sys, os, time
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


def run_direct_anneal(
    n_gen=52, dim=27,
    stage1_steps=30000, stage1_lr=0.001, stage1_gamma=50.0,
    stage2_steps=30000, stage2_lr=0.0001, stage2_gamma=200.0,
    alpha=5.0, beta=1.0,
    killing_target=-1.0,
    batch_size=1,
    seeds=None,
    verbose=True,
):
    """Direct optimizer with lr annealing, no curriculum trap."""
    device = DEVICE

    params = random_params_so(batch_size, n_gen, dim, device, seeds)

    # Normalize initial generators
    with torch.no_grad():
        gens = params_to_antisymmetric(params, n_gen, dim=dim)
        norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
        gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
        params = antisymmetric_to_params(gens)
    params = params.clone().detach().requires_grad_(True)

    K_target = killing_target * torch.eye(n_gen, dtype=FDTYPE, device=device)

    history = {"T_cl": [], "K_diag": [], "stage": []}

    for stage, (steps, lr, gamma) in enumerate([
        (stage1_steps, stage1_lr, stage1_gamma),
        (stage2_steps, stage2_lr, stage2_gamma),
    ], 1):
        if verbose:
            print(f"\n=== Stage {stage}: gamma={gamma}, lr={lr}, {steps} steps ===")

        optimizer = torch.optim.Adam([params], lr=lr)

        for step in range(steps):
            optimizer.zero_grad()
            generators = params_to_antisymmetric(params, n_gen, dim=dim)
            T_nc = batched_non_commutativity_tension(generators)
            K = batched_killing_metric(generators)
            T_k = ((K - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
            T_norm = batched_norm_tension(generators)
            T_cl = batched_closure_tension(generators)
            T_total = T_nc + alpha * T_k + beta * T_norm + gamma * T_cl
            loss = T_total.sum()
            loss.backward()
            optimizer.step()

            if verbose and step % 2000 == 0:
                with torch.no_grad():
                    kd = torch.stack([torch.diagonal(K[b]).mean()
                                      for b in range(batch_size)]).mean().item()
                    kd_std = torch.stack([torch.diagonal(K[b]).std()
                                          for b in range(batch_size)]).mean().item()
                    print(f"  Step {step:5d}: T={T_total.mean().item():.4f} "
                          f"T_cl={T_cl.mean().item():.6f} "
                          f"K_diag={kd:.4f}±{kd_std:.4e}")
                    history["T_cl"].append(T_cl.mean().item())
                    history["K_diag"].append(kd)
                    history["stage"].append(stage)

    # Final extraction
    with torch.no_grad():
        final_gens = params_to_antisymmetric(params, n_gen, dim=dim)
        final_np = final_gens.cpu().numpy()

    return final_np, history


def analyze_detailed(gens_np, name=""):
    """Detailed analysis with per-pair residuals."""
    n = len(gens_np)
    total_pairs = n * (n - 1) // 2
    B_mat = np.array([g.flatten() for g in gens_np]).T
    residuals = []
    for i in range(n):
        for j in range(i + 1, n):
            C = commutator_np(gens_np[i], gens_np[j])
            norm_c = np.linalg.norm(C)
            if norm_c < 1e-10:
                residuals.append(0.0)
            else:
                coeffs, _, _, _ = np.linalg.lstsq(B_mat, C.flatten(), rcond=None)
                res = np.linalg.norm(C.flatten() - B_mat @ coeffs)
                residuals.append(res)

    residuals = np.array(residuals)
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")
    print(f"  Residuals ({total_pairs} pairs):")
    print(f"    Max:    {residuals.max():.6e}")
    print(f"    Mean:   {residuals.mean():.6e}")
    print(f"    Median: {np.median(residuals):.6e}")
    print(f"    P99:    {np.percentile(residuals, 99):.6e}")

    for tol in [1e-2, 5e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5]:
        closes = int(np.sum(residuals < tol))
        pct = closes / total_pairs
        mark = " <<<" if pct == 1.0 else ""
        print(f"    tol={tol:.0e}: {closes}/{total_pairs} = {pct:.1%}{mark}")

    K = killing_metric_np(list(gens_np))
    K_diag = np.diag(K).real
    K_offdiag = np.abs(K - np.diag(K_diag)).mean()
    print(f"  K_diag: mean={K_diag.mean():.6f}  std={K_diag.std():.6e}")
    print(f"  K_offdiag: {K_offdiag:.6e}")

    # Simplicity
    K_normed = K / K_diag.mean() if abs(K_diag.mean()) > 1e-12 else K
    dev = np.linalg.norm(K_normed - np.eye(n)) / np.linalg.norm(np.eye(n))
    print(f"  ||K/tr - I||/||I|| = {dev:.6e}")

    # h_dual
    gens_r = [g / np.linalg.norm(g) for g in gens_np]
    K_r = killing_metric_np(gens_r)
    K_r_diag = np.diag(K_r).real
    Tr_sq = np.mean([np.trace(g @ g).real for g in gens_r])
    h_dual = abs(K_r_diag.mean() / Tr_sq) if abs(Tr_sq) > 1e-12 else 0
    print(f"  h_dual (I_R=1): {h_dual:.4f}")
    print(f"  h_dual (I_R=3): {3*h_dual:.4f} (F4: 9.0)")

    return residuals


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, nargs="+", default=[42, 137, 256])
    parser.add_argument("--s1-gamma", type=float, default=50.0)
    parser.add_argument("--s2-gamma", type=float, default=200.0)
    parser.add_argument("--s1-steps", type=int, default=30000)
    parser.add_argument("--s2-steps", type=int, default=30000)
    parser.add_argument("--s2-lr", type=float, default=0.0001)
    args = parser.parse_args()

    batch_size = len(args.seeds)
    print(f"F4 direct+anneal: {batch_size} seeds")
    print(f"  Stage 1: gamma={args.s1_gamma}, lr=0.001, {args.s1_steps} steps")
    print(f"  Stage 2: gamma={args.s2_gamma}, lr={args.s2_lr}, {args.s2_steps} steps")
    print(f"  Seeds: {args.seeds}")

    t0 = time.time()
    all_gens, history = run_direct_anneal(
        stage1_steps=args.s1_steps,
        stage1_lr=0.001,
        stage1_gamma=args.s1_gamma,
        stage2_steps=args.s2_steps,
        stage2_lr=args.s2_lr,
        stage2_gamma=args.s2_gamma,
        batch_size=batch_size,
        seeds=args.seeds,
    )
    elapsed = time.time() - t0
    print(f"\nTotal: {elapsed:.1f}s ({elapsed/60:.1f} min)")

    for b in range(batch_size):
        gens_b = [all_gens[b, k] for k in range(52)]
        analyze_detailed(gens_b, name=f"Seed {args.seeds[b]}")
