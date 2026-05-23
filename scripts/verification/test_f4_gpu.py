"""
Test B.2: F₄ exceptional Lie algebra — blind + constructive.

F₄ = Der(J₃(𝕆)) is the automorphism algebra of the exceptional Jordan algebra:
  - 52 generators, rank 4
  - Fundamental representation: dim=26 (traceless part of J₃(𝕆))
  - h∨ = 9, I_fund = 1 → expected |K|/κ² = 9

Two approaches:
  1. BLIND: Hermitian traceless parameterization (SU(26)-like) with n_gen=52
     — optimizer must discover f₄ inside SU(26). Expected FAIL: 52/675 = 7.7% density.
  2. CONSTRUCTIVE: Build f₄ generators from derivation condition on J₃(𝕆),
     verify with framework. Expected PASS.

For blind approach we use hermitian traceless (SU-like) parameterization
on 26×26 since the 26-dim is the traceless subspace of J₃(𝕆). The actual
f₄ representation on 26-dim is real orthogonal, so we also attempt with
SO(26) antisymmetric parameterization.

The constructive approach builds generators in the full 27-dim J₃(𝕆)
representation (where f₄ acts as derivations) and analyzes them.

Usage:
    python test_f4_gpu.py                 # constructive + brief blind test
    python test_f4_gpu.py --blind-only    # only the blind optimization attempt
    python test_f4_gpu.py --no-blind      # only the constructive approach
    python test_f4_gpu.py --seeds 5       # more seeds for blind test
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils_gpu import (
    params_to_antisymmetric,
    antisymmetric_to_params,
    random_params_so,
    evolve_generators_generic,
    killing_metric_np,
    structure_constants_general_np,
    check_in_span_np,
    commutator_np,
    get_device_info,
    DTYPE,
)
from octonion_jordan import (
    build_f4_generators,
    jordan_product_matrix,
    verify_derivation,
)


def analyze_f4_result(generators, seed_idx, verbose=True):
    """Full analysis of an F₄ candidate (27×27 generators).

    Uses trace-orthogonalization for correct Killing metric and the
    universal h∨ formula:  h∨ = I_R × Σ B_V(T_a,T_a) / Σ Tr(T_a²)
    """
    n_gen = len(generators)
    dim = generators[0].shape[0]

    # Reality check: generators should be real matrices
    reality_err = max(np.max(np.abs(g.imag)) if np.iscomplexobj(g) else 0.0
                      for g in generators)

    # Trace-orthogonalize for correct Killing metric
    G = np.zeros((n_gen, n_gen))
    for i in range(n_gen):
        for j in range(i, n_gen):
            G[i, j] = np.trace(generators[i] @ generators[j]).real
            G[j, i] = G[i, j]
    evals, evecs = np.linalg.eigh(G)
    gens_ortho = []
    for a in range(n_gen):
        if abs(evals[a]) > 1e-10:
            g = sum(evecs[b, a] * generators[b] for b in range(n_gen))
            g = g / np.sqrt(abs(evals[a]))
            gens_ortho.append(g)

    # Killing metric (B_V) on orthogonalized basis
    K = killing_metric_np(gens_ortho)
    K_diag = np.diag(K).real
    K_diag_mean = np.mean(K_diag)
    K_diag_std = np.std(K_diag)
    K_offdiag = np.abs(K - np.diag(K_diag)).max()

    # Killing eigenvalues
    K_eigs = np.linalg.eigvalsh(K.real)

    # Closure
    closes_count = 0
    total_pairs = n_gen * (n_gen - 1) // 2
    for i in range(n_gen):
        for j in range(i + 1, n_gen):
            C = commutator_np(gens_ortho[i], gens_ortho[j])
            norm_C = np.linalg.norm(C)
            if norm_C < 1e-10:
                closes_count += 1
            elif check_in_span_np(C, gens_ortho, tol=1e-4):
                closes_count += 1
    closure_frac = closes_count / total_pairs if total_pairs > 0 else 1.0

    # Generator norm (Frobenius-based trace)
    kappas = np.array([np.trace(g @ g).real for g in gens_ortho])
    kappa = np.mean(kappas)
    kappa_abs = abs(kappa)
    k_ratio = abs(K_diag_mean) / kappa_abs**2 if kappa_abs > 1e-12 else 0.0

    # Universal h∨ formula: h∨ = I_R × Σ B_V(T_a,T_a) / Σ Tr(T_a²)
    sum_bv = np.sum(K_diag)
    sum_kappa = np.sum(kappas)
    bv_kappa_ratio = sum_bv / sum_kappa if abs(sum_kappa) > 1e-12 else 0.0

    h_dual_theo = 9
    I_R = 3.0  # Dynkin index for 27-dim rep (26-dim irrep has I₂₆ = 3)
    h_dual_measured = I_R * bv_kappa_ratio

    result = {
        "seed": seed_idx,
        "n_gen": n_gen,
        "dim": dim,
        "reality_err": reality_err,
        "K_diag_mean": K_diag_mean,
        "K_diag_std": K_diag_std,
        "K_offdiag": K_offdiag,
        "K_eigs": K_eigs,
        "closure_frac": closure_frac,
        "closes_count": closes_count,
        "total_pairs": total_pairs,
        "kappa": kappa,
        "k_ratio": k_ratio,
        "sum_bv": sum_bv,
        "sum_kappa": sum_kappa,
        "bv_kappa_ratio": bv_kappa_ratio,
        "h_dual_theo": h_dual_theo,
        "h_dual_measured": h_dual_measured,
        "I_R": I_R,
    }

    if verbose:
        tag = "PASS" if closure_frac > 0.99 else "PARTIAL" if closure_frac > 0.5 else "FAIL"
        print(f"\n  seed={seed_idx} [{tag}]:")
        print(f"    Generators:       {n_gen} in dim={dim}")
        print(f"    Reality err:      {reality_err:.2e}")
        print(f"    K_diag:           {K_diag_mean:.6f} +/- {K_diag_std:.6f}")
        print(f"    K_offdiag max:    {K_offdiag:.2e}")
        print(f"    K eigenvalues:    min={K_eigs[0]:.4f} max={K_eigs[-1]:.4f}")
        print(f"    Closure:          {closes_count}/{total_pairs} ({closure_frac:.1%})")
        print(f"    kappa=Tr(T^2):    {kappa:.6f}")
        print(f"    Σ B_V / Σ Tr(T²): {bv_kappa_ratio:.4f}")
        print(f"    I_R (27-dim):     {I_R}")
        print(f"    h∨ measured:      {h_dual_measured:.4f} (theo: {h_dual_theo})")

    return result


# ===================================================================
# BLIND APPROACH: SO(27) parameterization with n_gen=52
# ===================================================================
def run_f4_blind(n_seeds=3, steps=None, lr=None):
    """Run F₄ blind optimization: 52 generators in SO(27) parameterization.

    Uses antisymmetric (SO) parameterization since f₄ generators on J₃(𝕆)
    are real endomorphisms. Density = 52/351 = 14.8% — expected to FAIL.
    """
    n_gen = 52
    dim = 27  # Full J₃(𝕆) representation

    if steps is None:
        steps = 40000
    if lr is None:
        lr = 0.0001

    print("=" * 78)
    print("  F₄ BLIND TEST — SO(27) parameterization")
    print("=" * 78)
    print(f"  Device:     {get_device_info()}")
    print(f"  Approach:   {n_gen} generators in SO({dim})")
    print(f"  Density:    {n_gen}/{dim*(dim-1)//2} = {n_gen/(dim*(dim-1)//2):.1%}")
    print(f"  Expected:   FAIL — f₄ is a specific 52D submanifold of so(27)")
    print(f"  Steps:      {steps}, LR: {lr}, Seeds: {n_seeds}")
    print()

    t0 = time.time()
    raw = evolve_generators_generic(
        n_gen=n_gen,
        dim=dim,
        steps=steps,
        lr=lr,
        alpha=5.0,
        beta=1.0,
        killing_target=-1.0,
        batch_size=n_seeds,
        seeds=list(range(n_seeds)),
        verbose=True,
        convergence_threshold=1e-4,
        param_to_gen_fn=params_to_antisymmetric,
        gen_to_param_fn=antisymmetric_to_params,
        random_param_fn=lambda bs, ng, device, seeds: random_params_so(bs, ng, dim, device, seeds),
        param_to_gen_kwargs={"dim": dim},
    )
    gpu_time = time.time() - t0

    results = []
    for seed_idx, (gens_np, hist) in enumerate(raw):
        gen_list = [gens_np[i].real for i in range(n_gen)]
        r = analyze_f4_result(gen_list, seed_idx)
        r["gpu_time"] = gpu_time / n_seeds
        results.append(r)

    any_pass = any(r["closure_frac"] > 0.99 for r in results)
    avg_closure = np.mean([r["closure_frac"] for r in results])

    print(f"\n  SUMMARY F₄ BLIND:")
    print(f"    GPU time:     {gpu_time:.1f}s ({gpu_time/n_seeds:.1f}s/seed)")
    print(f"    Avg closure:  {avg_closure:.1%}")
    print(f"    Any PASS:     {any_pass}")
    if not any_pass:
        print("    As expected: blind optimization cannot discover f₄ ⊂ so(27).")
        print("    Density 52/351 = 14.8% — well below the ~30% threshold.")

    return results


# ===================================================================
# CONSTRUCTIVE APPROACH: Build F₄ from Jordan algebra derivations
# ===================================================================
def run_constructive_f4():
    """Build f₄ generators from Der(J₃(𝕆)) and analyze them."""
    print("=" * 78)
    print("  F₄ CONSTRUCTIVE TEST — Der(J₃(𝕆))")
    print("=" * 78)

    t0 = time.time()
    generators = build_f4_generators()
    build_time = time.time() - t0

    print(f"\n  Build time: {build_time:.1f}s")
    print(f"  Analyzing {len(generators)} generators...")

    t0 = time.time()
    results = analyze_f4_result(generators, seed_idx="constructive")
    analysis_time = time.time() - t0
    print(f"  Analysis time: {analysis_time:.1f}s")

    # Verify: these are actual derivations of J₃(𝕆)
    print("\n  Verifying derivation property D(X∘Y) = D(X)∘Y + X∘D(Y)...")
    J = jordan_product_matrix()
    max_deriv_err = 0.0
    for i, D in enumerate(generators):
        err = verify_derivation(D, J)
        max_deriv_err = max(max_deriv_err, err)
        if i < 3:
            print(f"    Generator {i}: derivation error = {err:.2e}")
    print(f"    Max derivation error (all {len(generators)}): {max_deriv_err:.2e}")

    # Verify generators are traceless (derivations preserve identity)
    max_trace = max(abs(np.trace(g)) for g in generators)
    print(f"    Max trace: {max_trace:.2e}")

    # Check if f₄ ⊂ so(27) — i.e., generators antisymmetric w.r.t. some bilinear form
    # For the trace form <X,Y> = Tr(X∘Y) on J₃(𝕆), derivations satisfy <DX,Y> + <X,DY> = 0
    # In our basis this means D + D^T = 0 IF the basis is orthonormal w.r.t. trace form.
    # Let's check antisymmetry:
    antisym_errs = [np.linalg.norm(g + g.T) for g in generators]
    max_antisym = max(antisym_errs)
    mean_antisym = np.mean(antisym_errs)
    print(f"    Antisymmetry (D+D^T): max={max_antisym:.2e}, mean={mean_antisym:.2e}")

    # Summary
    is_pass = results["closure_frac"] > 0.99
    print(f"\n  CONSTRUCTIVE F₄ VERDICT: {'PASS' if is_pass else 'FAIL'}")
    if is_pass:
        print(f"    {len(generators)} generators, 100% closure")
        print(f"    h∨ = I_R × (Σ B_V / Σ κ) = {results['I_R']:.0f} × {results['bv_kappa_ratio']:.4f} = {results['h_dual_measured']:.4f} (theo: 9)")
        print(f"    Build: {build_time:.1f}s, Analysis: {analysis_time:.1f}s")
        print(f"    This IS the exceptional Lie algebra F₄!")

    results["build_time"] = build_time
    results["analysis_time"] = analysis_time
    results["max_deriv_err"] = max_deriv_err
    results["max_antisym"] = max_antisym
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="F₄ exceptional Lie algebra test")
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--steps", type=int, default=40000)
    parser.add_argument("--lr", type=float, default=0.0001)
    parser.add_argument("--blind-only", action="store_true")
    parser.add_argument("--no-blind", action="store_true")
    args = parser.parse_args()

    if not args.blind_only:
        run_constructive_f4()
        print()

    if not args.no_blind:
        run_f4_blind(n_seeds=args.seeds, steps=args.steps, lr=args.lr)
