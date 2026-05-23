"""
F₄ blind discovery with curriculum gamma scheduling.

Strategy:
  Phase 1: gamma = 0, lr = 0.003, 10k steps → fast Killing convergence
  Phase 2: gamma ramps 0 → gamma_max, lr = lr2, 50k steps → gradual closure

Previous results:
  - gamma=5,  30k: T_closure 80.9→0.45, max_res=3.9e-2 (FAR)
  - gamma=50, 30k: T_closure 80.9→0.0002, max_res=4.0e-4 (VERY CLOSE, tol=1e-4)
  
Hypothesis: curriculum ramp avoids early wasted gradient in wrong direction,
and the gradual increase lets the algebra "crystallize" incrementally.
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
    check_in_span_np,
    DTYPE, FDTYPE, DEVICE,
)
from test_closure_driven import batched_closure_tension


def run_curriculum(
    n_gen=52, dim=27,
    phase1_steps=10000, phase1_lr=0.003,
    phase2_steps=50000, phase2_lr=0.0005,
    gamma_max=200.0,
    alpha=5.0, beta=1.0,
    killing_target=-1.0,
    batch_size=3,
    seeds=None,
    verbose=True,
):
    """Two-phase optimizer: Killing first, then ramp closure."""
    device = DEVICE
    total_steps = phase1_steps + phase2_steps

    params = random_params_so(batch_size, n_gen, dim, device, seeds)

    # Normalize initial generators
    with torch.no_grad():
        gens = params_to_antisymmetric(params, n_gen, dim=dim)
        norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
        gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
        params = antisymmetric_to_params(gens)
    params = params.clone().detach().requires_grad_(True)

    K_target = killing_target * torch.eye(n_gen, dtype=FDTYPE, device=device)

    history = {
        "T_total": [], "T_nc": [], "T_k": [], "T_closure": [],
        "K_diag_mean": [], "gamma": [],
    }

    # Phase 1: gamma = 0
    if verbose:
        print(f"\n=== Phase 1: gamma=0, lr={phase1_lr}, {phase1_steps} steps ===")
    optimizer = torch.optim.Adam([params], lr=phase1_lr)
    for step in range(phase1_steps):
        optimizer.zero_grad()
        generators = params_to_antisymmetric(params, n_gen, dim=dim)
        T_nc = batched_non_commutativity_tension(generators)
        K = batched_killing_metric(generators)
        T_k = ((K - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
        T_norm = batched_norm_tension(generators)
        T_total = T_nc + alpha * T_k + beta * T_norm
        loss = T_total.sum()
        loss.backward()
        optimizer.step()

        if verbose and step % 2000 == 0:
            with torch.no_grad():
                T_cl = batched_closure_tension(generators)
                kd = torch.stack([torch.diagonal(K[b]).mean()
                                  for b in range(batch_size)]).mean().item()
                print(f"  Step {step:5d}: T={T_total.mean().item():.4f} "
                      f"T_cl={T_cl.mean().item():.4f} K_diag={kd:.4f} gamma=0")
                history["gamma"].append(0.0)
                history["T_total"].append(T_total.mean().item())
                history["T_nc"].append(T_nc.mean().item())
                history["T_k"].append(T_k.mean().item())
                history["T_closure"].append(T_cl.mean().item())
                history["K_diag_mean"].append(kd)

    # Phase 2: gamma ramps 0 → gamma_max
    if verbose:
        print(f"\n=== Phase 2: gamma 0→{gamma_max}, lr={phase2_lr}, {phase2_steps} steps ===")
    optimizer2 = torch.optim.Adam([params], lr=phase2_lr)

    for step in range(phase2_steps):
        frac = step / max(phase2_steps - 1, 1)
        gamma = gamma_max * frac

        optimizer2.zero_grad()
        generators = params_to_antisymmetric(params, n_gen, dim=dim)
        T_nc = batched_non_commutativity_tension(generators)
        K = batched_killing_metric(generators)
        T_k = ((K - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
        T_norm = batched_norm_tension(generators)
        T_cl = batched_closure_tension(generators)
        T_total = T_nc + alpha * T_k + beta * T_norm + gamma * T_cl
        loss = T_total.sum()
        loss.backward()
        optimizer2.step()

        if verbose and step % 2000 == 0:
            with torch.no_grad():
                kd = torch.stack([torch.diagonal(K[b]).mean()
                                  for b in range(batch_size)]).mean().item()
                print(f"  Step {step:5d}: T={T_total.mean().item():.4f} "
                      f"T_cl={T_cl.mean().item():.4f} K_diag={kd:.4f} "
                      f"gamma={gamma:.1f}")
                history["gamma"].append(gamma)
                history["T_total"].append(T_total.mean().item())
                history["T_nc"].append(T_nc.mean().item())
                history["T_k"].append(T_k.mean().item())
                history["T_closure"].append(T_cl.mean().item())
                history["K_diag_mean"].append(kd)

    # Extract final generators
    with torch.no_grad():
        final_gens = params_to_antisymmetric(params, n_gen, dim=dim)
        final_np = final_gens.cpu().numpy()

    return final_np, history


def analyze_with_tolerances(gens_np, name=""):
    """Analyze closure at multiple tolerances + Killing structure."""
    n = len(gens_np)
    total_pairs = n * (n - 1) // 2

    # Compute all residuals
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
    print(f"  Residual stats ({total_pairs} pairs):")
    print(f"    Max:        {residuals.max():.6e}")
    print(f"    Mean:       {residuals.mean():.6e}")
    print(f"    Median:     {np.median(residuals):.6e}")
    print(f"    P99:        {np.percentile(residuals, 99):.6e}")
    print(f"    P99.9:      {np.percentile(residuals, 99.9):.6e}")

    for tol in [1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5]:
        closes = int(np.sum(residuals < tol))
        pct = closes / total_pairs
        print(f"    tol={tol:.0e}: {closes}/{total_pairs} = {pct:.1%}"
              + (" <<<" if pct == 1.0 else ""))

    # Killing metric
    K = killing_metric_np(list(gens_np))
    K_diag = np.diag(K).real
    K_offdiag = np.abs(K - np.diag(K_diag)).mean()
    print(f"  K_diag: mean={K_diag.mean():.6f}  std={K_diag.std():.6e}")
    print(f"  K_offdiag: {K_offdiag:.6e}")

    # Renormalize and check h_dual
    norms = np.array([np.linalg.norm(g) for g in gens_np])
    gens_r = [g / np.linalg.norm(g) for g in gens_np]
    K_r = killing_metric_np(gens_r)
    K_r_diag = np.diag(K_r).real
    Tr_sq = np.mean([np.trace(g @ g).real for g in gens_r])
    h_dual = abs(K_r_diag.mean() / Tr_sq) if abs(Tr_sq) > 1e-12 else 0
    print(f"  h_dual (I_R=1, renorm): {h_dual:.4f}")
    print(f"  h_dual (I_R=3, renorm): {3*h_dual:.4f}  (F4 expected: 9.0)")

    # Simplicity check
    K_normed = K / K_diag.mean() if abs(K_diag.mean()) > 1e-12 else K
    dev = np.linalg.norm(K_normed - np.eye(n)) / np.linalg.norm(np.eye(n))
    print(f"  ||K/tr - I|| / ||I|| = {dev:.6e}  (simple if << 1)")

    return residuals


# ================================================================
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, nargs="+", default=[42, 137, 256])
    parser.add_argument("--gamma-max", type=float, default=200.0)
    parser.add_argument("--p1-steps", type=int, default=10000)
    parser.add_argument("--p2-steps", type=int, default=50000)
    parser.add_argument("--p2-lr", type=float, default=0.0005)
    args = parser.parse_args()

    batch_size = len(args.seeds)
    print(f"F4 curriculum: {batch_size} seeds, Phase1={args.p1_steps}+Phase2={args.p2_steps} steps")
    print(f"  gamma ramp: 0 → {args.gamma_max}, Phase2 lr={args.p2_lr}")
    print(f"  Seeds: {args.seeds}")

    t0 = time.time()
    all_gens, history = run_curriculum(
        n_gen=52, dim=27,
        phase1_steps=args.p1_steps,
        phase2_steps=args.p2_steps,
        phase1_lr=0.003,
        phase2_lr=args.p2_lr,
        gamma_max=args.gamma_max,
        batch_size=batch_size,
        seeds=args.seeds,
    )
    elapsed = time.time() - t0
    print(f"\nTotal elapsed: {elapsed:.1f}s ({elapsed/60:.1f} min)")

    # Analyze each seed
    for b in range(batch_size):
        gens_b = [all_gens[b, k] for k in range(52)]
        analyze_with_tolerances(gens_b, name=f"Seed {args.seeds[b]}")
