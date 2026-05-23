"""
Test D.1: Exploratory test — sl(2,R) non-compact Lie algebra.

sl(2,R) is the REAL form of A₁ with Killing form signature (+,+,−).
Compare with su(2), which is the COMPACT form with Killing signature (−,−,−).

Both have 3 generators in dim=2, but:
  - su(2): hermitian traceless generators → K negative-definite
  - sl(2,R): real traceless generators → K indefinite

Experiments:
  A. Real traceless param, target K=-I (compact target) — should fail or produce garbage
  B. Real traceless param, custom indefinite Killing target — should converge to sl(2,R)
  C. Hermitian traceless param, target K=-I (control) — same as SU(2), must converge

Usage:
    python test_sl2r_gpu.py
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import torch

from shared_utils_gpu import (
    evolve_generators_generic,
    killing_metric_np,
    structure_constants_general_np,
    check_in_span_np,
    commutator_np,
    get_device_info,
    DTYPE, FDTYPE, DEVICE,
    batched_killing_metric,
    batched_non_commutativity_tension,
    batched_norm_tension,
)


# ===================================================================
# Real traceless 2×2 parameterization for sl(2,R)
# ===================================================================
# A real traceless 2×2 matrix [[a, b], [c, -a]]: 3 real parameters.
# We embed in complex tensors (imaginary part = 0) for compatibility.

def params_to_sl2r(params, n_gen):
    """Real traceless 2×2 generators from 3 real parameters each.

    params: (B, n_gen, 3) complex tensor (only real part used)
    Returns: (B, n_gen, 2, 2) complex tensor (imaginary part = 0)
    """
    p = params.real
    B = p.shape[0]
    T = torch.zeros(B, n_gen, 2, 2, dtype=DTYPE, device=params.device)
    a = p[..., 0]  # (B, n_gen)
    b = p[..., 1]
    c = p[..., 2]
    T[..., 0, 0] = a
    T[..., 0, 1] = b
    T[..., 1, 0] = c
    T[..., 1, 1] = -a
    return T


def sl2r_to_params(generators):
    """Inverse of params_to_sl2r: extract (a, b, c) from generators."""
    B, n_gen = generators.shape[:2]
    p = torch.zeros(B, n_gen, 3, dtype=DTYPE, device=generators.device)
    p[..., 0] = generators[..., 0, 0].real
    p[..., 1] = generators[..., 0, 1].real
    p[..., 2] = generators[..., 1, 0].real
    return p


def random_params_sl2r(batch_size, n_gen, device, seeds):
    """Random parameters for sl(2,R): 3 real params per generator."""
    all_params = []
    for s in seeds:
        rng = np.random.default_rng(s)
        p = rng.standard_normal((n_gen, 3))
        all_params.append(p)
    params_np = np.stack(all_params)  # (B, n_gen, 3)
    return torch.tensor(params_np, dtype=DTYPE, device=device)


# ===================================================================
# Custom optimizer for indefinite Killing target
# ===================================================================
def evolve_sl2r_indefinite(
    n_gen=3,
    steps=20000,
    lr=0.001,
    alpha=5.0,
    beta=1.0,
    batch_size=3,
    seeds=None,
    verbose=True,
    convergence_threshold=1e-4,
):
    """Optimizer for sl(2,R) with non-uniform (indefinite) Killing target.

    Target: K_diag → sorted(eigenvalues) matches signature (+,+,-).
    Instead of K → k·I, we minimize:
        loss = T_nc + alpha * T_k_signature + beta * T_norm

    T_k_signature = sum over sorted pairs: (sorted_K_diag - sorted_target)²
    where sorted_target = (-λ, +λ, +λ) for λ=1.
    """
    device = DEVICE

    if seeds is None:
        seeds = list(range(batch_size))

    params = random_params_sl2r(batch_size, n_gen, device, seeds)
    with torch.no_grad():
        gens = params_to_sl2r(params, n_gen)
        norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
        gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
        params = sl2r_to_params(gens)
    params = params.clone().detach().requires_grad_(True)

    optimizer = torch.optim.Adam([params], lr=lr)

    # Indefinite target: one negative, two positive
    # sorted ascending: [-1, +1, +1]
    K_target_sorted = torch.tensor([-1.0, 1.0, 1.0], dtype=FDTYPE, device=device)

    prev_loss = torch.zeros(batch_size, device=device)
    active = torch.ones(batch_size, dtype=torch.bool, device=device)

    results_per_seed = [{"loss": [], "K_diag": []} for _ in range(batch_size)]

    for step in range(steps):
        optimizer.zero_grad()
        generators = params_to_sl2r(params, n_gen)  # (B, n_gen, 2, 2)

        # Closure tension
        T_nc = batched_non_commutativity_tension(generators)

        # Killing metric
        K = batched_killing_metric(generators)  # (B, n_gen, n_gen)
        K_diag = torch.diagonal(K, dim1=-2, dim2=-1)  # (B, n_gen)
        K_offdiag_mask = ~torch.eye(n_gen, dtype=torch.bool, device=device).unsqueeze(0)
        K_offdiag = K * K_offdiag_mask.float()

        # Signature-matched Killing loss:
        # Sort K_diag, compare with sorted target
        K_sorted, _ = K_diag.sort(dim=-1)  # (B, n_gen) ascending
        T_k_sig = ((K_sorted - K_target_sorted.unsqueeze(0)) ** 2).sum(dim=-1)  # (B,)
        T_k_off = (K_offdiag ** 2).sum(dim=(-2, -1))  # (B,) off-diagonal penalty

        # Norm tension
        T_norm = batched_norm_tension(generators)

        total = T_nc + alpha * (T_k_sig + T_k_off) + beta * T_norm
        loss = (total * active.float()).sum()
        loss.backward()
        optimizer.step()

        with torch.no_grad():
            for b in range(batch_size):
                if active[b]:
                    results_per_seed[b]["loss"].append(total[b].item())
                    results_per_seed[b]["K_diag"].append(K_diag[b].cpu().numpy().copy())

            converged = total < convergence_threshold
            if step > 0:
                plateau = (total - prev_loss).abs() < 1e-14
                newly_done = (converged | plateau) & active
                active = active & ~newly_done
            prev_loss = total.clone()

            if verbose and step % 2000 == 0:
                n_act = active.sum().item()
                mean_l = total[active].mean().item() if n_act > 0 else 0
                print(f"  Step {step}: active={n_act}/{batch_size} mean_loss={mean_l:.6f}")

            if not active.any():
                if verbose:
                    print(f"  All converged at step {step}")
                break

    with torch.no_grad():
        final = params_to_sl2r(params, n_gen).cpu().numpy()

    return [(final[b], results_per_seed[b]) for b in range(batch_size)]


# ===================================================================
# Analysis
# ===================================================================
def analyze_result(generators, label, verbose=True):
    """Analyze a set of 2×2 generators."""
    n_gen = len(generators)

    # Properties
    is_real = all(np.max(np.abs(g.imag)) < 1e-10 for g in generators)
    is_herm = all(np.linalg.norm(g - g.conj().T) < 1e-8 for g in generators)
    traces = [np.trace(g).real for g in generators]
    is_traceless = all(abs(t) < 1e-8 for t in traces)

    # Killing metric
    K = killing_metric_np(generators)
    K_diag = np.diag(K).real
    K_offdiag = np.abs(K - np.diag(K_diag)).max()

    # Signature: count positive, zero, negative eigenvalues
    eigs = np.linalg.eigvalsh(K.real)
    n_pos = np.sum(eigs > 1e-6)
    n_neg = np.sum(eigs < -1e-6)
    n_zero = n_gen - n_pos - n_neg

    # Closure
    closes = 0
    total_pairs = n_gen * (n_gen - 1) // 2
    for i in range(n_gen):
        for j in range(i + 1, n_gen):
            C = commutator_np(generators[i], generators[j])
            if np.linalg.norm(C) < 1e-10:
                closes += 1
            elif check_in_span_np(C, generators, tol=1e-4):
                closes += 1
    closure = closes / total_pairs if total_pairs > 0 else 1.0

    # Structure constants
    f = structure_constants_general_np(generators)
    f_rms = np.sqrt(np.mean(f**2))

    # Trace metric g_ab = Tr(T_a T_b)
    g_trace = np.zeros((n_gen, n_gen))
    for i in range(n_gen):
        for j in range(n_gen):
            g_trace[i, j] = np.trace(generators[i] @ generators[j]).real
    g_diag = np.diag(g_trace)

    result = {
        "label": label,
        "n_gen": n_gen,
        "is_real": is_real,
        "is_herm": is_herm,
        "is_traceless": is_traceless,
        "K_diag": K_diag,
        "K_offdiag": K_offdiag,
        "K_eigs": eigs,
        "signature": (int(n_pos), int(n_zero), int(n_neg)),
        "closure": closure,
        "f_rms": f_rms,
        "g_diag": g_diag,
    }

    if verbose:
        tag = "PASS" if closure > 0.99 else "FAIL"
        print(f"\n  [{tag}] {label}:")
        print(f"    Real: {is_real}, Hermitian: {is_herm}, Traceless: {is_traceless}")
        print(f"    K_diag:        {K_diag}")
        print(f"    K_offdiag max: {K_offdiag:.2e}")
        print(f"    K eigenvalues: {eigs}")
        print(f"    K signature:   ({n_pos}, {n_zero}, {n_neg})")
        print(f"    Closure:       {closes}/{total_pairs} ({closure:.1%})")
        print(f"    f_rms:         {f_rms:.6f}")
        print(f"    Tr(T_a^2):     {g_diag}")

    return result


def main():
    print("=" * 78)
    print("  sl(2,R) EXPLORATORY TEST — Compact vs Non-Compact Real Forms")
    print("=" * 78)
    print(f"  Device: {get_device_info()}")
    print(f"  Both algebras: 3 generators in dim=2")
    print()

    n_seeds = 3

    # ---------------------------------------------------------------
    # Experiment A: Real traceless + compact target K=-I
    # ---------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  EXPERIMENT A: Real traceless param + compact target (K = -I)")
    print("  Expected: should NOT converge (indefinite K can't equal -I)")
    print("=" * 78)

    t0 = time.time()
    results_a = evolve_generators_generic(
        n_gen=3, dim=2, steps=20000, lr=0.001,
        alpha=5.0, beta=1.0, killing_target=-1.0,
        batch_size=n_seeds, seeds=list(range(n_seeds)),
        verbose=True, convergence_threshold=1e-4,
        param_to_gen_fn=params_to_sl2r,
        gen_to_param_fn=sl2r_to_params,
        random_param_fn=random_params_sl2r,
        param_to_gen_kwargs={},
    )
    time_a = time.time() - t0

    for s, (gens, _) in enumerate(results_a):
        analyze_result(list(gens), f"Exp A seed={s}")

    # ---------------------------------------------------------------
    # Experiment B: Real traceless + indefinite target
    # ---------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  EXPERIMENT B: Real traceless param + indefinite target")
    print("  Target K_diag sorted: [-1, +1, +1]")
    print("  Expected: should converge to sl(2,R)")
    print("=" * 78)

    t0 = time.time()
    results_b = evolve_sl2r_indefinite(
        n_gen=3, steps=20000, lr=0.001,
        alpha=5.0, beta=1.0,
        batch_size=n_seeds, seeds=list(range(n_seeds)),
        verbose=True, convergence_threshold=1e-4,
    )
    time_b = time.time() - t0

    for s, (gens, _) in enumerate(results_b):
        analyze_result(list(gens), f"Exp B seed={s}")

    # ---------------------------------------------------------------
    # Experiment C: Hermitian traceless + compact target (SU(2) control)
    # ---------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  EXPERIMENT C: Hermitian traceless param + compact target (SU(2) control)")
    print("  Expected: 100% closure, K ∝ -I")
    print("=" * 78)

    t0 = time.time()
    results_c = evolve_generators_generic(
        n_gen=3, dim=2, steps=20000, lr=0.001,
        alpha=5.0, beta=1.0, killing_target=-1.0,
        batch_size=n_seeds, seeds=list(range(n_seeds)),
        verbose=True, convergence_threshold=1e-4,
    )
    time_c = time.time() - t0

    for s, (gens, _) in enumerate(results_c):
        analyze_result(list(gens), f"Exp C seed={s}")

    # ---------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------
    print("\n\n" + "=" * 78)
    print("  SUMMARY")
    print("=" * 78)
    print(f"\n  {'Experiment':35s} {'Closure':>10s} {'K signature':>15s} {'Time':>10s}")
    print("  " + "-" * 75)

    for label, results, t in [("A: sl(2,R) + compact target", results_a, time_a),
                                ("B: sl(2,R) + indefinite target", results_b, time_b),
                                ("C: su(2) control",               results_c, time_c)]:
        for s, (gens, _) in enumerate(results):
            r = analyze_result(list(gens), f"{label} s{s}", verbose=False)
            sig = r["signature"]
            print(f"  {label + f' s{s}':35s} {r['closure']:>9.1%} "
                  f"  ({sig[0]}+,{sig[1]}0,{sig[2]}-) {t/n_seeds:>9.1f}s")

    print("  " + "-" * 75)
    print("\n  Key question: Does Killing form signature discriminate compact vs non-compact?")
    print("  - Compact (su(2)): K signature (0+, 0 zero, 3−)")
    print("  - Non-compact (sl(2,R)): K signature (2+, 0 zero, 1−)")


if __name__ == "__main__":
    main()
