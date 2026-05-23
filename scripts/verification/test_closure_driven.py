"""
Closure-driven blind search for exceptional Lie algebras.

KEY IDEA: The existing optimizer uses only Killing metric + non-commutativity
tension. This is sufficient for classical algebras (SU, SO, Sp) but fails
for exceptionals (G₂, F₄, E₆) because the optimizer converges to generic
SO(n) minima that have the right Killing metric but NO closure.

This test adds a DIFFERENTIABLE CLOSURE TENSION to the loss function:

  T_closure = Σ_{i<j} ||[T_i, T_j] - proj_span([T_i, T_j])||²

which explicitly penalizes commutators that leak outside the algebra span.

  T_total = T_nc + α·T_K + β·T_norm + γ·T_closure

Hypothesis: With closure penalty, the optimizer can discover G₂ (and
potentially F₄, E₆) from a blind start.

No existing files are modified. Self-contained optimizer with closure term.

Usage:
    python test_closure_driven.py                      # G₂ test (default)
    python test_closure_driven.py --algebra g2          # G₂ only
    python test_closure_driven.py --algebra g2 --compare # G₂ with/without closure
    python test_closure_driven.py --algebra su3         # SU(3) sanity check
    python test_closure_driven.py --algebra f4          # F₄ (long!)
    python test_closure_driven.py --gamma 0             # disable closure (baseline)
"""
import sys
import os
import time
import argparse

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
    params_to_generators,
    generators_to_params,
    random_params,
    killing_metric_np,
    commutator_np,
    check_in_span_np,
    structure_constants_general_np,
    get_device_info,
    DTYPE,
    FDTYPE,
    DEVICE,
)


# ============================================================
# DIFFERENTIABLE CLOSURE TENSION (GPU)
# ============================================================

def batched_closure_tension(generators):
    """Differentiable closure tension: penalizes commutators outside the span.

    generators: (B, n, d, d) complex tensor (requires grad)
    Returns: (B,) real tensor

    For each pair (i<j), computes [T_i, T_j] and projects onto span{T_k}
    using COMPLEX least squares. This correctly handles both:
      - Antisymmetric generators (SO/G₂): [T,T] ∈ span{T}  (real coefficients)
      - Hermitian generators (SU):        [T,T] = i·Σf T_k  (imaginary coefficients)

    Uses: proj_V(c) = V @ (V^H V)^{-1} @ V^H @ c  (complex inner product)

    Complexity: O(n² · d² + n³) — dominated by n² commutators.
    """
    B, n, d, _ = generators.shape

    # Flatten generators: (B, n, d²) complex
    G_flat = generators.reshape(B, n, d * d)

    # Complex Gram matrix: Gram[a,c] = <T_a, T_c> = Σ_k conj(T_a_k) * T_c_k
    Gram = torch.bmm(G_flat.conj(), G_flat.transpose(1, 2))  # (B, n, n)

    # Regularize for numerical stability
    reg = 1e-8 * torch.eye(n, dtype=FDTYPE, device=generators.device).unsqueeze(0)
    Gram_inv = torch.linalg.inv(Gram + reg.to(DTYPE))  # (B, n, n)

    # Commutators for upper triangle pairs
    idx_i, idx_j = torch.triu_indices(n, n, offset=1, device=generators.device)
    G_i = generators[:, idx_i]  # (B, n_pairs, d, d)
    G_j = generators[:, idx_j]
    C = G_i @ G_j - G_j @ G_i  # (B, n_pairs, d, d)
    C_flat = C.reshape(B, -1, d * d)  # (B, n_pairs, d²)

    # Inner products: v[a,p] = <T_a, C_p> = Σ_k conj(G_flat[a,k]) * C_flat[p,k]
    inner = torch.bmm(G_flat.conj(), C_flat.transpose(1, 2))  # (B, n, n_pairs)

    # Coefficients: c = Gram_inv @ v
    coeffs = torch.bmm(Gram_inv, inner)  # (B, n, n_pairs)

    # Projections: proj[k,p] = Σ_a coeffs[a,p] * G_flat[a,k]
    proj = torch.bmm(G_flat.transpose(1, 2), coeffs)  # (B, d², n_pairs)

    # Residuals
    residuals = C_flat.transpose(1, 2) - proj  # (B, d², n_pairs)

    # Total closure tension: Σ |residual|²
    T_closure = (residuals.abs() ** 2).sum(dim=1).sum(dim=1)  # (B,)

    return T_closure.real


# ============================================================
# CLOSURE-DRIVEN OPTIMIZER
# ============================================================

def evolve_with_closure(
    n_gen,
    dim,
    steps=10000,
    lr=0.003,
    alpha=5.0,
    beta=1.0,
    gamma=10.0,
    killing_target=-1.0,
    batch_size=1,
    seeds=None,
    verbose=True,
    convergence_threshold=1e-6,
    param_type="antisymmetric",
    snapshot_every=None,
):
    """GPU optimizer with differentiable closure tension.

    Loss = T_nc + α·T_K + β·T_norm + γ·T_closure

    param_type: "antisymmetric" for SO(n), "hermitian" for SU(n)
    snapshot_every: if set, record numpy generators at these intervals

    Returns list of (generators_np, history) per seed.
    """
    device = DEVICE

    if param_type == "antisymmetric":
        param_fn = params_to_antisymmetric
        gen_to_param_fn = antisymmetric_to_params
        rand_fn = lambda bs, ng, device, seeds: random_params_so(bs, ng, dim, device, seeds)
        ptg_kwargs = {"dim": dim}
    else:
        param_fn = params_to_generators
        gen_to_param_fn = generators_to_params
        rand_fn = lambda bs, ng, device, seeds: random_params(bs, ng, dim, device, seeds)
        ptg_kwargs = {"dim": dim}

    params = rand_fn(batch_size, n_gen, device=device, seeds=seeds)

    # Normalize initial generators
    with torch.no_grad():
        gens = param_fn(params, n_gen, **ptg_kwargs)
        norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
        gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
        params = gen_to_param_fn(gens)
    params = params.clone().detach().requires_grad_(True)

    optimizer = torch.optim.Adam([params], lr=lr)
    K_target_mat = killing_target * torch.eye(n_gen, dtype=FDTYPE, device=device)

    per_seed = [{"T_total": [], "T_nc": [], "T_k": [], "T_norm": [],
                 "T_closure": [], "K_diag": [], "K_offdiag": []}
                for _ in range(batch_size)]

    snapshots = {b: [] for b in range(batch_size)} if snapshot_every else None

    prev_totals = torch.zeros(batch_size, device=device)
    active = torch.ones(batch_size, dtype=torch.bool, device=device)

    for step in range(steps):
        optimizer.zero_grad()
        generators = param_fn(params, n_gen, **ptg_kwargs)

        # Compute all tension components
        T_nc = batched_non_commutativity_tension(generators)
        K = batched_killing_metric(generators)
        T_k = ((K - K_target_mat.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
        T_norm = batched_norm_tension(generators)
        T_closure = batched_closure_tension(generators)

        T_total = T_nc + alpha * T_k + beta * T_norm + gamma * T_closure

        loss = (T_total * active.float()).sum()
        loss.backward()
        optimizer.step()

        with torch.no_grad():
            for b in range(batch_size):
                if not active[b]:
                    continue
                per_seed[b]["T_total"].append(T_total[b].item())
                per_seed[b]["T_nc"].append(T_nc[b].item())
                per_seed[b]["T_k"].append(T_k[b].item())
                per_seed[b]["T_norm"].append(T_norm[b].item())
                per_seed[b]["T_closure"].append(T_closure[b].item())
                per_seed[b]["K_diag"].append(torch.diagonal(K[b]).mean().item())
                per_seed[b]["K_offdiag"].append(
                    (K[b] - torch.diag(torch.diagonal(K[b]))).abs().mean().item())

            # Snapshots
            if snapshot_every and step % snapshot_every == 0:
                gens_np = generators.cpu().numpy()
                for b in range(batch_size):
                    if active[b]:
                        snapshots[b].append((step, gens_np[b].copy()))

            # Convergence / plateau detection
            converged = T_total < convergence_threshold
            if step > 0:
                plateau = (T_total - prev_totals).abs() < 1e-14
                newly_done = (converged | plateau) & active
                active = active & ~newly_done
            prev_totals = T_total.clone()

            if verbose and step % 1000 == 0:
                n_active = active.sum().item()
                mean_t = T_total[active].mean().item() if n_active > 0 else 0
                mean_cl = T_closure[active].mean().item() if n_active > 0 else 0
                mean_kd = torch.stack([torch.diagonal(K[b]).mean() for b in range(batch_size)
                                       if active[b]]).mean().item() if n_active > 0 else 0
                print(f"  Step {step:5d}: active={n_active}/{batch_size} "
                      f"T={mean_t:.4f} T_cl={mean_cl:.4f} K_diag={mean_kd:.4f}")

            if not active.any():
                if verbose:
                    print(f"  All seeds done at step {step}")
                break

    # Final extraction
    with torch.no_grad():
        final_gens = param_fn(params, n_gen, **ptg_kwargs)
        final_np = final_gens.cpu().numpy()

    results = []
    for b in range(batch_size):
        r = (final_np[b], per_seed[b])
        if snapshots:
            # Append final snapshot
            snapshots[b].append((step, final_np[b].copy()))
        results.append(r)

    return results, snapshots


# ============================================================
# ANALYSIS (numpy, post-optimization)
# ============================================================

def analyze_result(generators, name="", verbose=True):
    """Full analysis of an emerged generator set."""
    n = len(generators)
    dim = generators[0].shape[0]

    # Killing metric
    K = killing_metric_np(generators)
    K_diag = np.diag(K).real
    K_diag_mean = np.mean(K_diag)
    K_diag_std = np.std(K_diag)
    K_offdiag = np.abs(K - np.diag(K_diag)).mean()

    # Closure
    closes = 0
    total_pairs = n * (n - 1) // 2
    max_residual = 0.0
    closure_energy = 0.0
    total_comm_energy = 0.0

    for i in range(n):
        for j in range(i + 1, n):
            C = commutator_np(generators[i], generators[j])
            norm_c = np.linalg.norm(C)
            total_comm_energy += norm_c ** 2
            if norm_c < 1e-10:
                closes += 1
            elif check_in_span_np(C, generators, tol=1e-4):
                closes += 1
            else:
                # Compute residual
                B_mat = np.array([g.flatten() for g in generators]).T
                coeffs, _, _, _ = np.linalg.lstsq(B_mat, C.flatten(), rcond=None)
                res = np.linalg.norm(C.flatten() - B_mat @ coeffs)
                max_residual = max(max_residual, res)
                closure_energy += res ** 2

    closure_frac = closes / total_pairs if total_pairs > 0 else 1.0
    closure_leak = closure_energy / total_comm_energy if total_comm_energy > 1e-20 else 0.0

    # h∨ estimation
    kappa_sq = np.mean([np.trace(g @ g).real for g in generators])
    sum_K = np.sum(K_diag)
    sum_Tr = sum(np.trace(g @ g).real for g in generators)
    h_dual = abs(sum_K / sum_Tr) if abs(sum_Tr) > 1e-12 else 0

    if verbose:
        print(f"\n  Analysis: {name}")
        print(f"    n_gen={n}, dim={dim}")
        print(f"    K_diag: mean={K_diag_mean:.6f} std={K_diag_std:.6f}")
        print(f"    K_offdiag: {K_offdiag:.6f}")
        print(f"    Closure: {closes}/{total_pairs} = {closure_frac:.1%}")
        print(f"    Closure leak: {closure_leak:.4e}")
        print(f"    Max residual: {max_residual:.4e}")
        print(f"    h_dual (I_R=1): {h_dual:.4f}")

    return {
        "name": name,
        "n_gen": n, "dim": dim,
        "K_diag_mean": K_diag_mean,
        "K_diag_std": K_diag_std,
        "K_offdiag": K_offdiag,
        "closure_count": closes,
        "closure_frac": closure_frac,
        "closure_leak": closure_leak,
        "max_residual": max_residual,
        "h_dual": h_dual,
    }


# ============================================================
# TEST CONFIGURATIONS
# ============================================================

ALGEBRAS = {
    "su2": {
        "n_gen": 3, "dim": 2, "param_type": "hermitian",
        "expected_h": 2.0, "steps": 3000, "lr": 0.005,
        "label": "SU(2) [3 gens, 2x2]",
    },
    "su3": {
        "n_gen": 8, "dim": 3, "param_type": "hermitian",
        "expected_h": 3.0, "steps": 5000, "lr": 0.003,
        "label": "SU(3) [8 gens, 3x3]",
    },
    "so3": {
        "n_gen": 3, "dim": 3, "param_type": "antisymmetric",
        "expected_h": 2.0, "steps": 3000, "lr": 0.005,
        "label": "SO(3) [3 gens, 3x3]",
    },
    "g2": {
        "n_gen": 14, "dim": 7, "param_type": "antisymmetric",
        "expected_h": 4.0, "steps": 15000, "lr": 0.003,
        "label": "G2 [14 gens, 7x7 in SO(7)]",
    },
    "f4": {
        "n_gen": 52, "dim": 27, "param_type": "antisymmetric",
        "expected_h": 9.0, "steps": 30000, "lr": 0.001,
        "label": "F4 [52 gens, 27x27 in SO(27)]",
    },
}


def run_comparison(algebra_key, seeds=3, gamma_values=None):
    """Run blind optimization with different gamma values and compare."""
    if gamma_values is None:
        gamma_values = [0.0, 5.0, 20.0, 50.0]

    cfg = ALGEBRAS[algebra_key]
    print(f"\n{'='*70}")
    print(f"  CLOSURE-DRIVEN COMPARISON: {cfg['label']}")
    print(f"  Seeds: {seeds}, Steps: {cfg['steps']}")
    print(f"  Gamma values: {gamma_values}")
    print(f"{'='*70}")

    all_results = {}

    for gamma in gamma_values:
        label = f"gamma={gamma}"
        print(f"\n--- {label} ---")
        t0 = time.time()

        results, _ = evolve_with_closure(
            n_gen=cfg["n_gen"],
            dim=cfg["dim"],
            steps=cfg["steps"],
            lr=cfg["lr"],
            alpha=5.0,
            beta=1.0,
            gamma=gamma,
            killing_target=-1.0,
            batch_size=seeds,
            seeds=list(range(42, 42 + seeds)),
            verbose=True,
            param_type=cfg["param_type"],
        )

        elapsed = time.time() - t0
        print(f"  Elapsed: {elapsed:.1f}s")

        analyses = []
        for b, (gens, hist) in enumerate(results):
            a = analyze_result(list(gens), name=f"{label} seed={b}")
            analyses.append(a)

        all_results[gamma] = analyses

    # Summary table
    print(f"\n{'='*80}")
    print(f"  COMPARISON SUMMARY: {cfg['label']}")
    print(f"  Expected h_dual = {cfg['expected_h']}")
    print(f"{'='*80}")
    print(f"  {'gamma':>8} {'seed':>4} {'K_diag':>10} {'closure':>10} "
          f"{'leak':>10} {'h_dual':>8} {'match':>6}")
    print(f"  {'-'*8} {'-'*4} {'-'*10} {'-'*10} {'-'*10} {'-'*8} {'-'*6}")

    for gamma, analyses in all_results.items():
        for a in analyses:
            match = "YES" if (a["closure_frac"] > 0.95
                              and abs(a["h_dual"] - cfg["expected_h"]) < 0.5) else "no"
            seed_idx = analyses.index(a)
            print(f"  {gamma:>8.1f} {seed_idx:>4} {a['K_diag_mean']:>10.4f} "
                  f"{a['closure_frac']:>9.1%} {a['closure_leak']:>10.2e} "
                  f"{a['h_dual']:>8.4f} {match:>6}")

    return all_results


def run_single(algebra_key, gamma=20.0, seeds=3, steps=None):
    """Run a single configuration."""
    cfg = ALGEBRAS[algebra_key]
    if steps:
        cfg = {**cfg, "steps": steps}

    print(f"\n{'='*70}")
    print(f"  CLOSURE-DRIVEN BLIND SEARCH: {cfg['label']}")
    print(f"  gamma={gamma}, seeds={seeds}, steps={cfg['steps']}")
    print(f"{'='*70}")

    t0 = time.time()

    results, snapshots = evolve_with_closure(
        n_gen=cfg["n_gen"],
        dim=cfg["dim"],
        steps=cfg["steps"],
        lr=cfg["lr"],
        alpha=5.0,
        beta=1.0,
        gamma=gamma,
        killing_target=-1.0,
        batch_size=seeds,
        seeds=list(range(42, 42 + seeds)),
        verbose=True,
        param_type=cfg["param_type"],
        snapshot_every=cfg["steps"] // 10,
    )

    elapsed = time.time() - t0
    print(f"\n  Total elapsed: {elapsed:.1f}s")

    for b, (gens, hist) in enumerate(results):
        a = analyze_result(list(gens), name=f"Seed {b}")

    # History plot summary
    print(f"\n  TENSION HISTORY (final values, seed 0):")
    h = results[0][1]
    for key in ["T_total", "T_nc", "T_k", "T_closure"]:
        if h[key]:
            print(f"    {key}: {h[key][0]:.4f} -> {h[key][-1]:.6f}")

    return results


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Closure-driven blind search")
    parser.add_argument("--algebra", default="g2",
                        choices=list(ALGEBRAS.keys()),
                        help="Which algebra to test")
    parser.add_argument("--gamma", type=float, default=20.0,
                        help="Closure penalty weight (0 = disabled)")
    parser.add_argument("--seeds", type=int, default=3,
                        help="Number of random seeds")
    parser.add_argument("--steps", type=int, default=None,
                        help="Override step count")
    parser.add_argument("--compare", action="store_true",
                        help="Run comparison across multiple gamma values")
    args = parser.parse_args()

    print("=" * 70)
    print("  CLOSURE-DRIVEN BLIND SEARCH FOR EXCEPTIONAL LIE ALGEBRAS")
    print("=" * 70)
    print(f"  Device: {get_device_info()}")
    print(f"  Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    if args.compare:
        run_comparison(args.algebra, seeds=args.seeds)
    else:
        run_single(args.algebra, gamma=args.gamma, seeds=args.seeds,
                   steps=args.steps)


if __name__ == "__main__":
    main()
