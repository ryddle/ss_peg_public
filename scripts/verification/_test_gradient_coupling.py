"""
Diagnostic test: Compare gradients for seed 256 at bs=1 vs bs=3.

If gradients are bit-identical, Adam trajectories must be identical,
and the batch_size effect is impossible → some other mechanism is at play.

If gradients differ at bit-level, floating-point non-determinism
in CUDA kernels explains the different convergence outcomes.
"""

import torch
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

from shared_utils_gpu import (
    DEVICE, FDTYPE, DTYPE,
    params_to_antisymmetric,
    antisymmetric_to_params,
    batched_killing_metric,
    batched_non_commutativity_tension,
    batched_norm_tension,
    random_params_so,
)
from test_closure_driven import batched_closure_tension

N_GEN = 52
DIM   = 27
ALPHA = 1.0
BETA  = 0.1


def init_params(seeds):
    """Initialize and normalize params for given seeds (same as _landscape_montecarlo.py)."""
    B = len(seeds)
    params = random_params_so(B, N_GEN, DIM, DEVICE, seeds)
    with torch.no_grad():
        gens = params_to_antisymmetric(params, N_GEN, dim=DIM)
        norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
        gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
        params = antisymmetric_to_params(gens)
    return params.clone().detach().requires_grad_(True)


def forward_and_grad(params, gamma=50.0):
    """Run one forward pass and return loss and gradients."""
    K_target = -1.0 * torch.eye(N_GEN, dtype=FDTYPE, device=DEVICE)
    generators = params_to_antisymmetric(params, N_GEN, dim=DIM)
    T_nc = batched_non_commutativity_tension(generators)
    K = batched_killing_metric(generators)
    T_k = ((K - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
    T_norm = batched_norm_tension(generators)
    T_cl = batched_closure_tension(generators)
    T_total = T_nc + ALPHA * T_k + BETA * T_norm + gamma * T_cl
    loss = T_total.sum()
    loss.backward()
    return T_total.detach(), params.grad.detach().clone()


def main():
    print("=" * 70)
    print("  GRADIENT COUPLING DIAGNOSTIC")
    print("=" * 70)

    # --- Test 1: Init params comparison ---
    print("\n[1] Comparing initial params for seed 256...")
    p1 = init_params([256])        # bs=1
    p3 = init_params([42, 137, 256])  # bs=3, seed 256 at index 2

    # Extract seed 256's params
    p256_alone = p1[0].detach().cpu()
    p256_batch = p3[2].detach().cpu()

    diff_init = (p256_alone - p256_batch).abs().max().item()
    print(f"  Max |diff| in initial params: {diff_init:.2e}")
    if diff_init == 0.0:
        print("  ✓ Initial params are BIT-IDENTICAL")
    else:
        print("  ✗ Initial params DIFFER — seed generation depends on batch context!")
        return

    # --- Test 2: Single forward+backward gradient comparison ---
    print("\n[2] Comparing gradients after one forward+backward pass...")

    # bs=1
    p1 = init_params([256])
    T1, g1 = forward_and_grad(p1)

    # bs=3 — seed 256 at index 2
    p3 = init_params([42, 137, 256])
    T3, g3 = forward_and_grad(p3)

    g256_alone = g1[0].cpu()
    g256_batch = g3[2].cpu()

    T256_alone = T1[0].cpu().item()
    T256_batch = T3[2].cpu().item()

    diff_T = abs(T256_alone - T256_batch)
    diff_grad_max = (g256_alone - g256_batch).abs().max().item()
    diff_grad_mean = (g256_alone - g256_batch).abs().mean().item()
    diff_grad_rel = (g256_alone - g256_batch).abs() / (g256_alone.abs() + 1e-30)
    diff_grad_rel_max = diff_grad_rel.max().item()

    print(f"  T_total[256] alone : {T256_alone:.10e}")
    print(f"  T_total[256] batch : {T256_batch:.10e}")
    print(f"  |diff T|           : {diff_T:.2e}")
    print()
    print(f"  grad256 alone shape: {g256_alone.shape}")
    print(f"  grad256 batch shape: {g256_batch.shape}")
    print(f"  Max |grad diff|    : {diff_grad_max:.2e}")
    print(f"  Mean |grad diff|   : {diff_grad_mean:.2e}")
    print(f"  Max relative diff  : {diff_grad_rel_max:.2e}")

    if diff_grad_max == 0.0:
        print("\n  ✓ Gradients are BIT-IDENTICAL")
        print("  → Adam trajectories MUST be identical.")
        print("  → batch_size effect is IMPOSSIBLE through this mechanism.")
        print("  → Some OTHER mechanism must explain the difference.")
    elif diff_grad_max < 1e-10:
        print(f"\n  ~ Gradients differ at epsilon level (max={diff_grad_max:.2e})")
        print("  → Floating-point non-determinism in CUDA kernels.")
        print("  → Over 11000 steps, these tiny diffs compound → different basins.")
    else:
        print(f"\n  ✗ Gradients differ SIGNIFICANTLY (max={diff_grad_max:.2e})")
        print("  → There IS computational coupling through batch operations.")

    # --- Test 3: Different partner seeds ---
    print("\n[3] Comparing gradients with different partners...")
    p3b = init_params([100, 200, 256])  # different partners
    T3b, g3b = forward_and_grad(p3b)

    g256_other = g3b[2].cpu()
    T256_other = T3b[2].cpu().item()

    diff_partners = (g256_batch - g256_other).abs().max().item()
    print(f"  T_total[256] {'{42,137}':>10s}: {T256_batch:.10e}")
    print(f"  T_total[256] {'{100,200}':>10s}: {T256_other:.10e}")
    print(f"  Max |grad diff| between partners: {diff_partners:.2e}")

    if diff_partners == 0.0:
        print("  ✓ Gradients identical regardless of partners → no partner coupling")
    else:
        print(f"  ~ Partners affect gradient at {diff_partners:.2e} level")

    # --- Test 4: Extended divergence ---
    TOTAL_STEPS = 8000
    LOG_STEPS = set([0,1,2,5,10,20,50,100,200,500,1000,1500,2000,
                     3000,4000,5000,6000,7000,7999])
    print(f"\n[4] Extended divergence test ({TOTAL_STEPS} Adam steps, Stage 1)...")
    p_solo = init_params([256])
    p_trio = init_params([42, 137, 256])

    opt_solo = torch.optim.Adam([p_solo], lr=0.001)
    opt_trio = torch.optim.Adam([p_trio], lr=0.001)

    K_target = -1.0 * torch.eye(N_GEN, dtype=FDTYPE, device=DEVICE)

    for step in range(TOTAL_STEPS):
        # Solo
        opt_solo.zero_grad()
        g_s = params_to_antisymmetric(p_solo, N_GEN, dim=DIM)
        T_nc_s = batched_non_commutativity_tension(g_s)
        K_s = batched_killing_metric(g_s)
        T_k_s = ((K_s - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
        T_norm_s = batched_norm_tension(g_s)
        T_cl_s = batched_closure_tension(g_s)
        T_s = T_nc_s + ALPHA * T_k_s + BETA * T_norm_s + 50.0 * T_cl_s
        T_s.sum().backward()
        opt_solo.step()

        # Trio
        opt_trio.zero_grad()
        g_t = params_to_antisymmetric(p_trio, N_GEN, dim=DIM)
        T_nc_t = batched_non_commutativity_tension(g_t)
        K_t = batched_killing_metric(g_t)
        T_k_t = ((K_t - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
        T_norm_t = batched_norm_tension(g_t)
        T_cl_t = batched_closure_tension(g_t)
        T_t = T_nc_t + ALPHA * T_k_t + BETA * T_norm_t + 50.0 * T_cl_t
        T_t.sum().backward()
        opt_trio.step()

        if step in LOG_STEPS:
            with torch.no_grad():
                diff = (p_solo[0] - p_trio[2]).abs().max().item()
                g_s_eval = params_to_antisymmetric(p_solo, N_GEN, dim=DIM)
                g_t_eval = params_to_antisymmetric(p_trio, N_GEN, dim=DIM)
                Tcl_s = batched_closure_tension(g_s_eval)[0].item()
                Tcl_t = batched_closure_tension(g_t_eval)[2].item()
            print(f"  step {step:5d}: max|Δp| = {diff:.2e}  "
                  f"T_solo={T_s[0].item():.4e}  T_trio={T_t[2].item():.4e}  "
                  f"Tcl_solo={Tcl_s:.4e}  Tcl_trio={Tcl_t:.4e}  "
                  f"ΔTcl={abs(Tcl_s - Tcl_t):.2e}", flush=True)

    # --- Stage 2 ---
    print(f"\n[5] Stage 2 ({3000} steps, lr=0.0001, gamma=200)...")
    opt_solo2 = torch.optim.Adam([p_solo], lr=0.0001)
    opt_trio2 = torch.optim.Adam([p_trio], lr=0.0001)

    S2_STEPS = 3000
    S2_LOG = set([0,100,500,1000,1500,2000,2500,2999])

    for step in range(S2_STEPS):
        opt_solo2.zero_grad()
        g_s = params_to_antisymmetric(p_solo, N_GEN, dim=DIM)
        T_nc_s = batched_non_commutativity_tension(g_s)
        K_s = batched_killing_metric(g_s)
        T_k_s = ((K_s - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
        T_norm_s = batched_norm_tension(g_s)
        T_cl_s = batched_closure_tension(g_s)
        T_s = T_nc_s + ALPHA * T_k_s + BETA * T_norm_s + 200.0 * T_cl_s
        T_s.sum().backward()
        opt_solo2.step()

        opt_trio2.zero_grad()
        g_t = params_to_antisymmetric(p_trio, N_GEN, dim=DIM)
        T_nc_t = batched_non_commutativity_tension(g_t)
        K_t = batched_killing_metric(g_t)
        T_k_t = ((K_t - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
        T_norm_t = batched_norm_tension(g_t)
        T_cl_t = batched_closure_tension(g_t)
        T_t = T_nc_t + ALPHA * T_k_t + BETA * T_norm_t + 200.0 * T_cl_t
        T_t.sum().backward()
        opt_trio2.step()

        if step in S2_LOG:
            with torch.no_grad():
                diff = (p_solo[0] - p_trio[2]).abs().max().item()
                g_s_eval = params_to_antisymmetric(p_solo, N_GEN, dim=DIM)
                g_t_eval = params_to_antisymmetric(p_trio, N_GEN, dim=DIM)
                Tcl_s = batched_closure_tension(g_s_eval)[0].item()
                Tcl_t = batched_closure_tension(g_t_eval)[2].item()
            print(f"  s2 step {step:5d}: max|Δp| = {diff:.2e}  "
                  f"T_solo={T_s[0].item():.4e}  T_trio={T_t[2].item():.4e}  "
                  f"Tcl_solo={Tcl_s:.4e}  Tcl_trio={Tcl_t:.4e}  "
                  f"ΔTcl={abs(Tcl_s - Tcl_t):.2e}", flush=True)

    # Final classification comparison
    print("\n[6] Final state comparison...")
    with torch.no_grad():
        g_s_final = params_to_antisymmetric(p_solo, N_GEN, dim=DIM)
        g_t_final = params_to_antisymmetric(p_trio, N_GEN, dim=DIM)
        Tcl_s_final = batched_closure_tension(g_s_final)[0].item()
        Tcl_t_final = batched_closure_tension(g_t_final)[2].item()
        Tnc_s = batched_non_commutativity_tension(g_s_final)[0].item()
        Tnc_t = batched_non_commutativity_tension(g_t_final)[2].item()
        K_s = batched_killing_metric(g_s_final)
        K_t = batched_killing_metric(g_t_final)
        Tk_s = ((K_s[0] - K_target) ** 2).sum().item()
        Tk_t = ((K_t[2] - K_target) ** 2).sum().item()

        param_diff = (p_solo[0] - p_trio[2]).abs().max().item()

    print(f"  max|Δparams|  : {param_diff:.2e}")
    print(f"  T_cl  solo    : {Tcl_s_final:.6e}")
    print(f"  T_cl  trio    : {Tcl_t_final:.6e}")
    print(f"  T_nc  solo    : {Tnc_s:.6e}")
    print(f"  T_nc  trio    : {Tnc_t:.6e}")
    print(f"  T_k   solo    : {Tk_s:.6e}")
    print(f"  T_k   trio    : {Tk_t:.6e}")

    print("\n" + "=" * 70)
    print("  DIAGNOSTIC COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
