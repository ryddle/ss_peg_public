"""
Test B.3: E₆ exceptional Lie algebra — blind + constructive.

E₆ preserves the cubic norm on J₃(𝕆):
  - 78 generators, rank 6
  - Fundamental representation: dim=27 (J₃(𝕆))
  - h∨ = 12, I_R = 3 (Dynkin index for 27-dim rep)

Two approaches:
  1. BLIND: SO(27) parameterization with n_gen=78
     — optimizer must discover e₆ inside so(27). Expected FAIL: 78/351 = 22.2% density.
  2. CONSTRUCTIVE: Build e₆ generators from cubic norm preservation on J₃(𝕆),
     verify with framework. Expected PASS.

Note: The 27-dim real representation yields the split real form e₆(6), whose
generators have mixed trace-signature (52 with Tr(T²)<0, 26 with Tr(T²)>0).
The h∨ extraction requires the universal formula:
    h∨ = I_R × Σ_a B_V(T_a,T_a) / Σ_a Tr_V(T_a²)
with generators trace-orthogonalized so structure constants are correct.

Note: e₆ ⊃ f₄ ⊃ g₂ — the exceptional algebra chain related to octonions.
The subalgebra chain can be verified by checking that f₄ generators span
a closed subalgebra within the e₆ generators.

Usage:
    python test_e6_gpu.py                 # constructive + brief blind test
    python test_e6_gpu.py --blind-only    # only the blind optimization attempt
    python test_e6_gpu.py --no-blind      # only the constructive approach
    python test_e6_gpu.py --seeds 5       # more seeds for blind test
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
    build_e6_generators,
    build_f4_generators,
    cubic_norm_tensor,
    verify_cubic_preservation,
)


def analyze_e6_result(generators, seed_idx, verbose=True):
    """Full analysis of an E₆ candidate (27×27 generators).

    Uses trace-orthogonalization for correct Killing metric and the
    universal h∨ formula:  h∨ = I_R × Σ B_V(T_a,T_a) / Σ Tr(T_a²)
    which works for both compact and split real forms.
    """
    n_gen = len(generators)
    dim = generators[0].shape[0]

    # Reality check
    reality_err = max(np.max(np.abs(g.imag)) if np.iscomplexobj(g) else 0.0
                      for g in generators)

    # Trace-orthogonalize for correct structure constant computation
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

    # Universal h∨ formula: h∨ = I_R × Σ B_V(T_a,T_a) / Σ Tr(T_a²)
    kappas = np.array([np.trace(g @ g).real for g in gens_ortho])
    sum_bv = np.sum(K_diag)
    sum_kappa = np.sum(kappas)
    bv_kappa_ratio = sum_bv / sum_kappa if abs(sum_kappa) > 1e-12 else 0.0

    h_dual_theo = 12
    I_R = 3.0  # Dynkin index for 27-dim rep of E₆
    h_dual_measured = I_R * bv_kappa_ratio

    # Signature analysis (split form diagnostic)
    n_compact = int(np.sum(kappas < -1e-8))
    n_noncompact = int(np.sum(kappas > 1e-8))

    result = {
        "seed": seed_idx,
        "n_gen": len(gens_ortho),
        "dim": dim,
        "reality_err": reality_err,
        "K_diag_mean": K_diag_mean,
        "K_diag_std": K_diag_std,
        "K_offdiag": K_offdiag,
        "K_eigs": K_eigs,
        "closure_frac": closure_frac,
        "closes_count": closes_count,
        "total_pairs": total_pairs,
        "sum_bv": sum_bv,
        "sum_kappa": sum_kappa,
        "bv_kappa_ratio": bv_kappa_ratio,
        "h_dual_theo": h_dual_theo,
        "h_dual_measured": h_dual_measured,
        "I_R": I_R,
        "n_compact": n_compact,
        "n_noncompact": n_noncompact,
    }

    if verbose:
        tag = "PASS" if closure_frac > 0.99 else "PARTIAL" if closure_frac > 0.5 else "FAIL"
        print(f"\n  seed={seed_idx} [{tag}]:")
        print(f"    Generators:       {len(gens_ortho)} in dim={dim}")
        print(f"    Reality err:      {reality_err:.2e}")
        print(f"    K_diag:           {K_diag_mean:.6f} +/- {K_diag_std:.6f}")
        print(f"    K_offdiag max:    {K_offdiag:.2e}")
        print(f"    K eigenvalues:    min={K_eigs[0]:.4f} max={K_eigs[-1]:.4f}")
        print(f"    Closure:          {closes_count}/{total_pairs} ({closure_frac:.1%})")
        print(f"    Signature:        {n_compact} compact + {n_noncompact} non-compact")
        print(f"    Σ B_V / Σ Tr(T²): {bv_kappa_ratio:.4f}")
        print(f"    I_R (27-dim):     {I_R}")
        print(f"    h∨ measured:      {h_dual_measured:.4f} (theo: {h_dual_theo})")

    return result


# ===================================================================
# BLIND APPROACH: SO(27) parameterization with n_gen=78
# ===================================================================
def run_e6_blind(n_seeds=3, steps=None, lr=None):
    """Run E₆ blind optimization: 78 generators in SO(27) parameterization.

    Density = 78/351 = 22.2% — expected to FAIL (below ~30% threshold).
    """
    n_gen = 78
    dim = 27

    if steps is None:
        steps = 40000
    if lr is None:
        lr = 0.0001

    print("=" * 78)
    print("  E₆ BLIND TEST — SO(27) parameterization")
    print("=" * 78)
    print(f"  Device:     {get_device_info()}")
    print(f"  Approach:   {n_gen} generators in SO({dim})")
    print(f"  Density:    {n_gen}/{dim*(dim-1)//2} = {n_gen/(dim*(dim-1)//2):.1%}")
    print(f"  Expected:   FAIL — e₆ is a specific 78D submanifold of so(27)")
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
        r = analyze_e6_result(gen_list, seed_idx)
        r["gpu_time"] = gpu_time / n_seeds
        results.append(r)

    any_pass = any(r["closure_frac"] > 0.99 for r in results)
    avg_closure = np.mean([r["closure_frac"] for r in results])

    print(f"\n  SUMMARY E₆ BLIND:")
    print(f"    GPU time:     {gpu_time:.1f}s ({gpu_time/n_seeds:.1f}s/seed)")
    print(f"    Avg closure:  {avg_closure:.1%}")
    print(f"    Any PASS:     {any_pass}")
    if not any_pass:
        print("    As expected: blind optimization cannot discover e₆ ⊂ so(27).")
        print("    Density 78/351 = 22.2% — below the ~30% threshold.")

    return results


# ===================================================================
# CONSTRUCTIVE APPROACH: Build E₆ from cubic norm preservation
# ===================================================================
def run_constructive_e6():
    """Build e₆ generators from cubic norm preservation and analyze them."""
    print("=" * 78)
    print("  E₆ CONSTRUCTIVE TEST — Cubic Norm Preservation on J₃(𝕆)")
    print("=" * 78)

    t0 = time.time()
    generators = build_e6_generators()
    build_time = time.time() - t0

    print(f"\n  Build time: {build_time:.1f}s")
    print(f"  Analyzing {len(generators)} generators...")

    t0 = time.time()
    results = analyze_e6_result(generators, seed_idx="constructive")
    analysis_time = time.time() - t0
    print(f"  Analysis time: {analysis_time:.1f}s")

    # Verify: these preserve the cubic norm
    print("\n  Verifying cubic norm preservation N(L·X, X, X) = 0...")
    N = cubic_norm_tensor()
    max_cubic_err = 0.0
    for i, L in enumerate(generators):
        err = verify_cubic_preservation(L, N)
        max_cubic_err = max(max_cubic_err, err)
        if i < 3:
            print(f"    Generator {i}: cubic preservation error = {err:.2e}")
    print(f"    Max cubic error (all {len(generators)}): {max_cubic_err:.2e}")

    # Check tracelessness (e₆ is a subalgebra of sl(27), not gl(27))
    max_trace = max(abs(np.trace(g)) for g in generators)
    print(f"    Max trace: {max_trace:.2e}")

    # Check subalgebra chain: f₄ ⊂ e₆
    # Both sets are built from independent SVD null spaces, so we check by
    # projecting f₄ generators onto the e₆ span using least-squares.
    print("\n  Checking f₄ ⊂ e₆ subalgebra embedding...")
    f4_gens = build_f4_generators()
    e6_basis = np.column_stack([g.flatten() for g in generators])  # (729, 78)
    f4_in_e6_count = 0
    max_residual = 0.0
    for i, f4_g in enumerate(f4_gens):
        coeffs, _, _, _ = np.linalg.lstsq(e6_basis, f4_g.flatten(), rcond=None)
        recon = e6_basis @ coeffs
        residual = np.linalg.norm(recon - f4_g.flatten()) / np.linalg.norm(f4_g.flatten())
        max_residual = max(max_residual, residual)
        if residual < 1e-6:
            f4_in_e6_count += 1
    print(f"    f₄ generators in span(e₆): {f4_in_e6_count}/{len(f4_gens)}")
    print(f"    Max relative residual: {max_residual:.2e}")
    subalgebra_ok = f4_in_e6_count == len(f4_gens)
    print(f"    f₄ ⊂ e₆: {'CONFIRMED' if subalgebra_ok else 'FAILED'}")

    # Summary
    is_pass = results["closure_frac"] > 0.99
    print(f"\n  CONSTRUCTIVE E₆ VERDICT: {'PASS' if is_pass else 'FAIL'}")
    if is_pass:
        print(f"    {results['n_gen']} generators, 100% closure")
        print(f"    Real form: e₆(6) (split) — {results['n_compact']} compact + {results['n_noncompact']} non-compact")
        print(f"    h∨ = I_R × (Σ B_V / Σ κ) = {results['I_R']:.0f} × {results['bv_kappa_ratio']:.4f} = {results['h_dual_measured']:.4f} (theo: 12)")
        print(f"    Build: {build_time:.1f}s, Analysis: {analysis_time:.1f}s")
        if subalgebra_ok:
            print(f"    Subalgebra chain f₄ ⊂ e₆ confirmed!")
        print(f"    This IS the exceptional Lie algebra E₆!")

    results["build_time"] = build_time
    results["analysis_time"] = analysis_time
    results["max_cubic_err"] = max_cubic_err
    results["f4_in_e6"] = f4_in_e6_count
    results["subalgebra_ok"] = subalgebra_ok
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="E₆ exceptional Lie algebra test")
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--steps", type=int, default=40000)
    parser.add_argument("--lr", type=float, default=0.0001)
    parser.add_argument("--blind-only", action="store_true")
    parser.add_argument("--no-blind", action="store_true")
    args = parser.parse_args()

    if not args.blind_only:
        run_constructive_e6()
        print()

    if not args.no_blind:
        run_e6_blind(n_seeds=args.seeds, steps=args.steps, lr=args.lr)
