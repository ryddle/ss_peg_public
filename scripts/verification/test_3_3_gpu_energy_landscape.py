"""
Test 3.3 GPU: Energy landscape scan — batched on GPU.
Question: Is SU(2) at the global minimum of the energy landscape?

GPU version: Evaluates 2000+ random configurations in parallel batches on GPU.

EXPECTED: SU(2) configurations cluster at the global minimum of T_Killing.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import torch
from shared_utils_gpu import (
    params_to_generators,
    batched_killing_metric,
    batched_non_commutativity_tension,
    batched_norm_tension,
    pauli_matrices_torch,
    generators_to_params,
    evolve_generators_gpu,
    get_device_info,
    DEVICE, DTYPE, FDTYPE,
)


def check_is_su2_structure_np(generators, tol=0.1):
    """Heuristic: closure + K proportional to I (numpy)."""
    if len(generators) != 3:
        return False
    for i, j in [(0, 1), (1, 2), (2, 0)]:
        comm = generators[i] @ generators[j] - generators[j] @ generators[i]
        B = np.array([b.flatten() for b in generators]).T
        v = comm.flatten()
        coeffs, _, _, _ = np.linalg.lstsq(B, v, rcond=None)
        if np.linalg.norm(B @ coeffs - v) > tol:
            return False
    from shared_utils_gpu import killing_metric_np
    K = killing_metric_np(generators)
    K_offdiag = np.max(np.abs(K - np.diag(np.diag(K))))
    return K_offdiag <= tol


def run_test_landscape_gpu(n_samples=2000, n_optimized=100):
    print("=" * 70)
    print("TEST 3.3 GPU: ENERGY LANDSCAPE SCAN")
    print("=" * 70)
    print(f"  Device: {get_device_info()}")
    print()

    device = DEVICE
    n_gen, dim = 3, 2

    # ---- Phase 1: Random samples (batched on GPU) ----
    print(f"  Phase 1: Evaluating {n_samples} random configurations on GPU...")
    t0 = time.time()

    # Generate random parameters
    torch.manual_seed(0)
    n_params = dim * dim - 1
    params = torch.randn(n_samples, n_gen, n_params, dtype=FDTYPE, device=device)

    # Convert to generators and normalize
    gens = params_to_generators(params, n_gen, dim)
    norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
    gens = gens / norms.unsqueeze(-1).unsqueeze(-1)

    # Evaluate tensions
    T_nc = batched_non_commutativity_tension(gens).cpu().numpy()
    K = batched_killing_metric(gens)
    K_target = -1.0 * torch.eye(n_gen, dtype=FDTYPE, device=device)
    T_k = ((K - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1)).cpu().numpy()
    T_total = T_nc + 5.0 * T_k

    random_time = time.time() - t0
    print(f"    Done in {random_time:.2f}s")

    random_results = []
    for i in range(n_samples):
        random_results.append({
            "T_nc": T_nc[i],
            "T_k": T_k[i],
            "T_total": T_total[i],
        })

    # ---- Phase 2: Optimized samples (batched GPU optimizer) ----
    print(f"\n  Phase 2: Running {n_optimized} optimized seeds on GPU...")
    t0 = time.time()

    opt_results_raw, opt_hist = evolve_generators_gpu(
        n_gen=3, dim=2, steps=15000,
        lr=0.001, alpha=5.0, beta=1.0,
        killing_target=-1.0,
        batch_size=n_optimized,
        seeds=list(range(n_optimized)),
        verbose=True,
    )

    # Evaluate final tensions for optimized
    opt_gens_torch = torch.tensor(opt_results_raw, dtype=DTYPE, device=device)
    T_nc_opt = batched_non_commutativity_tension(opt_gens_torch).cpu().numpy()
    K_opt = batched_killing_metric(opt_gens_torch)
    T_k_opt = ((K_opt - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1)).cpu().numpy()

    opt_time = time.time() - t0
    print(f"    Optimizer done in {opt_time:.1f}s")

    opt_results = []
    for i in range(n_optimized):
        opt_results.append({
            "T_nc": T_nc_opt[i],
            "T_k": T_k_opt[i],
        })

    # ---- Phase 3: Reference (Pauli/2) ----
    pauli = pauli_matrices_torch(device)
    pauli_half = torch.stack([s / 2 for s in pauli]).unsqueeze(0)  # (1, 3, 2, 2)
    T_nc_pauli = batched_non_commutativity_tension(pauli_half).item()
    K_pauli = batched_killing_metric(pauli_half)
    T_k_pauli = ((K_pauli - K_target.unsqueeze(0)) ** 2).sum().item()

    # ---- Analysis ----
    print("\n--- LANDSCAPE ANALYSIS ---")
    print(f"  Random samples ({n_samples}):")
    print(f"    T_nc range:  [{min(T_nc):.2f}, {max(T_nc):.2f}]")
    print(f"    T_k range:   [{min(T_k):.2f}, {max(T_k):.2f}]")

    print(f"\n  Optimized ({n_optimized}):")
    opt_tnc = [r["T_nc"] for r in opt_results]
    opt_tk = [r["T_k"] for r in opt_results]
    print(f"    T_nc range:  [{min(opt_tnc):.4f}, {max(opt_tnc):.4f}]")
    print(f"    T_k range:   [{min(opt_tk):.6f}, {max(opt_tk):.6f}]")

    print(f"\n  Pauli/2 (exact SU(2)):")
    print(f"    T_nc = {T_nc_pauli:.4f}")
    print(f"    T_k  = {T_k_pauli:.6f}")

    # ---- Save plot ----
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(10, 7))

        ax.scatter(T_nc, T_k, alpha=0.15, s=8, c='gray', label=f'Random ({n_samples})')
        ax.scatter(opt_tnc, opt_tk, alpha=0.6, s=25, c='blue', label=f'Optimized ({n_optimized})')
        ax.scatter([T_nc_pauli], [T_k_pauli], s=200, c='red', marker='*',
                   zorder=10, label='Pauli/2 (exact)')

        ax.set_xlabel('T_nc (non-commutativity tension)')
        ax.set_ylabel('T_Killing (Killing metric tension)')
        ax.set_title('Energy Landscape: Random vs Optimized vs Exact SU(2) [GPU]')
        ax.legend()
        ax.set_yscale('log')

        plot_path = os.path.join(os.path.dirname(__file__), "test_3_3_gpu_energy_landscape.png")
        fig.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"\n  Plot saved: {plot_path}")
    except ImportError:
        print("\n  (matplotlib not available — skipping plot)")

    # ---- Verdict ----
    print("\n--- VERDICT ---")
    optimized_near_pauli = sum(1 for r in opt_results if r["T_k"] < 0.5)
    if T_k_pauli < min(T_k) and optimized_near_pauli > n_optimized * 0.8:
        print(f"  PASS: Pauli/2 sits at global minimum of T_Killing.")
        print(f"  {optimized_near_pauli}/{n_optimized} optimized configs near Pauli basin.")
    else:
        print(f"  PARTIAL: {optimized_near_pauli}/{n_optimized} optimized near Pauli basin.")

    total_time = random_time + opt_time
    print(f"\n  [GPU total: {total_time:.1f}s ({n_samples} random + {n_optimized} optimized)]")
    return random_results, opt_results


if __name__ == "__main__":
    run_test_landscape_gpu(n_samples=2000, n_optimized=100)
