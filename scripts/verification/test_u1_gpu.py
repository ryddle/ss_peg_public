"""
§13: U(1) Abelian Emergence Tests

U(1) is the gauge group of electromagnetism — and we completely forgot about it!
It's abelian, so its Killing metric is identically zero: K_ab = Tr(ad_a · ad_b) = 0
(because [T,T] = 0 for a single generator).

This script explores several questions:

  Exp A: n_gen = 1, dim = 2, hermitian traceless (SU param)
         → Single generator from su(2) basis. Trivially K=0, closure=trivial.
         → What does the optimizer do when K_target=-1 but K must be 0?

  Exp B: n_gen = d², dim = d, general hermitian (WITH trace)
         → u(d) = u(1) ⊕ su(d) should emerge
         → The optimizer has d²-1 "useful" generators and 1 extra
         → Does it discover the u(1) ⊕ su(d) decomposition?
         → Does the trace generator decouple (K=0)?

  Exp C: n_gen = d²-1+1 = d², dim = d, for d = 2,3
         → Same as B but we analyze the Killing spectrum carefully
         → Look for the zero eigenvalue corresponding to U(1)

  Exp D: n_gen = 1, dim = d, pure abelian limit
         → With K_target = 0 (correct for abelian)
         → Does it converge to a normalized generator?

Usage:
    python test_u1_gpu.py
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import torch

from shared_utils_gpu import (
    evolve_generators_generic,
    killing_metric_np,
    check_in_span_np,
    commutator_np,
    get_device_info,
    params_to_generators,
    generators_to_params,
    random_params,
    DTYPE, FDTYPE, DEVICE,
    batched_killing_metric,
    batched_non_commutativity_tension,
    batched_norm_tension,
)


# ===================================================================
# PARAMETERIZATION: General hermitian (WITH trace)
# ===================================================================
# d² real parameters per generator (one more than traceless)

def params_to_hermitian(params, n_gen, dim=2):
    """Convert real parameters to general hermitian matrices (with trace).

    params: (B, n_gen, d²) real tensor
    Returns: (B, n_gen, dim, dim) complex tensor

    Structure: d diagonal reals + d(d-1)/2 complex off-diagonal = d² params
    """
    B = params.shape[0]
    generators = torch.zeros(B, n_gen, dim, dim, dtype=DTYPE, device=params.device)

    for g in range(n_gen):
        p = params[:, g, :]  # (B, d²)
        idx = 0
        H = torch.zeros(B, dim, dim, dtype=DTYPE, device=params.device)

        # Full diagonal (d independent values — NO tracelessness constraint)
        for i in range(dim):
            H[:, i, i] = p[:, idx].to(DTYPE)
            idx += 1

        # Off-diagonal: real and imaginary parts (same as traceless version)
        for i in range(dim):
            for j in range(i + 1, dim):
                re = p[:, idx]
                im = p[:, idx + 1]
                H[:, i, j] = (re + 1j * im).to(DTYPE)
                H[:, j, i] = (re - 1j * im).to(DTYPE)
                idx += 2

        generators[:, g] = H

    return generators


def hermitian_to_params(generators):
    """Extract real parameters from general hermitian matrices.

    generators: (B, n_gen, dim, dim) complex tensor
    Returns: (B, n_gen, d²) real tensor
    """
    B, n_gen, dim, _ = generators.shape
    n_params = dim * dim
    params = torch.zeros(B, n_gen, n_params, dtype=FDTYPE, device=generators.device)

    for g in range(n_gen):
        H = generators[:, g]
        idx = 0
        # Full diagonal
        for i in range(dim):
            params[:, g, idx] = H[:, i, i].real
            idx += 1
        # Off-diagonal
        for i in range(dim):
            for j in range(i + 1, dim):
                params[:, g, idx] = H[:, i, j].real
                params[:, g, idx + 1] = H[:, i, j].imag
                idx += 2

    return params


def random_params_hermitian(batch_size, n_gen, dim, device=DEVICE, seeds=None):
    """Random parameters for general hermitian generators."""
    n_params = dim * dim
    if seeds is not None:
        params_list = []
        for s in seeds:
            np.random.seed(s)
            p = np.random.randn(1, n_gen, n_params)
            params_list.append(p)
        params_np = np.concatenate(params_list, axis=0)
        return torch.tensor(params_np, dtype=FDTYPE, device=device)
    return torch.randn(batch_size, n_gen, n_params, dtype=FDTYPE, device=device)


# ===================================================================
# ANALYSIS
# ===================================================================

def analyze_result(generators, label, verbose=True):
    """Full analysis of emerged generators."""
    n_gen = len(generators)
    dim = generators[0].shape[0]

    # Killing metric
    K = killing_metric_np(generators)
    K_diag = np.diag(K).real
    K_eigs = np.sort(np.linalg.eigvalsh(K.real))
    K_offdiag = np.abs(K - np.diag(K_diag)).max()

    # Signature
    n_pos = int(np.sum(K_eigs > 1e-6))
    n_neg = int(np.sum(K_eigs < -1e-6))
    n_zero = n_gen - n_pos - n_neg

    # Hermiticity check
    herm_err = max(np.linalg.norm(g - g.conj().T) for g in generators)

    # Trace analysis
    traces = [np.trace(g) for g in generators]
    trace_norms = [abs(t) for t in traces]

    # Frobenius norms
    frob_norms = [np.linalg.norm(g) for g in generators]

    # Closure check
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

    # Identify U(1) vs SU generators via trace
    u1_indices = [i for i, t in enumerate(trace_norms) if t > 0.01]
    su_indices = [i for i, t in enumerate(trace_norms) if t <= 0.01]

    if verbose:
        print(f"\n  [{label}]:")
        print(f"    Hermitian: {herm_err < 1e-8} (max err: {herm_err:.2e})")
        print(f"    K eigenvalues: {K_eigs}")
        print(f"    K_diag: {K_diag}")
        print(f"    K signature: ({n_pos}+, {n_zero}0, {n_neg}−)")
        print(f"    K_offdiag max: {K_offdiag:.2e}")
        print(f"    Closure: {closes}/{total_pairs} ({closure:.1%})")
        print(f"    Frobenius norms: {[f'{f:.4f}' for f in frob_norms]}")
        print(f"    Traces: {[f'{t:.4f}' for t in traces]}")
        print(f"    |Traces|: {[f'{t:.4f}' for t in trace_norms]}")
        if u1_indices:
            print(f"    U(1) generators (|Tr|>0.01): indices {u1_indices}")
        if su_indices:
            print(f"    SU generators (|Tr|≤0.01): indices {su_indices}")

    return {
        "K": K, "K_eigs": K_eigs, "K_diag": K_diag,
        "K_offdiag": K_offdiag, "signature": (n_pos, n_zero, n_neg),
        "closure": closure, "frob_norms": frob_norms,
        "traces": traces, "trace_norms": trace_norms,
        "u1_indices": u1_indices, "su_indices": su_indices,
    }


# ===================================================================
# EXPERIMENT A: Single generator, SU parameterization
# ===================================================================
def experiment_A():
    """n_gen=1, dim=2, hermitian traceless. K_target=-1.
    A single generator has no commutator partners → K=0 trivially.
    What does the optimizer do?
    """
    print("=" * 80)
    print("  EXPERIMENT A: n_gen=1, dim=2, hermitian traceless (SU param)")
    print("  K_target = -1, but K must be 0 (single generator)")
    print("=" * 80)

    n_gen, dim = 1, 2
    seeds = [42, 137, 256]

    results = evolve_generators_generic(
        n_gen=n_gen, dim=dim, steps=10000, lr=0.001,
        alpha=5.0, beta=1.0, killing_target=-1.0,
        batch_size=len(seeds), seeds=seeds,
        verbose=True, convergence_threshold=1e-8,
    )

    for i, (gens_np, hist) in enumerate(results):
        analyze_result(list(gens_np), f"Exp A seed={seeds[i]}")
        print(f"    Final T_total: {hist['T_total'][-1]:.6f}")
        print(f"    Final T_k: {hist['T_k'][-1]:.6f}")
        print(f"    Final T_nc: {hist['T_nc'][-1]:.6f}")

    return results


# ===================================================================
# EXPERIMENT B: n_gen=d², general hermitian → u(d) emergence
# ===================================================================
def experiment_B(dim=2):
    """n_gen=d², dim=d, general hermitian. Should give u(d) = u(1) ⊕ su(d).
    One generator should end up proportional to I (the U(1) part).
    """
    n_gen = dim * dim

    print("\n" + "=" * 80)
    print(f"  EXPERIMENT B: n_gen={n_gen}, dim={dim}, general hermitian")
    print(f"  Expected: u({dim}) = u(1) ⊕ su({dim})")
    print(f"  {n_gen-1} generators with K≠0 (su({dim}) part)")
    print(f"   1 generator with K=0 (u(1) part)")
    print("=" * 80)

    seeds = [42, 137, 256, 314, 999]

    results = evolve_generators_generic(
        n_gen=n_gen, dim=dim, steps=15000, lr=0.001,
        alpha=5.0, beta=1.0, killing_target=-1.0,
        batch_size=len(seeds), seeds=seeds,
        verbose=True, convergence_threshold=1e-6,
        param_to_gen_fn=lambda p, ng: params_to_hermitian(p, ng, dim),
        gen_to_param_fn=hermitian_to_params,
        random_param_fn=lambda bs, ng, device=DEVICE, seeds=None: random_params_hermitian(
            bs, ng, dim, device, seeds),
        param_to_gen_kwargs={},
    )

    for i, (gens_np, hist) in enumerate(results):
        r = analyze_result(list(gens_np), f"Exp B dim={dim} seed={seeds[i]}")

        # Check if U(1) generator is proportional to identity
        for idx in r["u1_indices"]:
            g = gens_np[idx]
            g_normalized = g / np.linalg.norm(g) if np.linalg.norm(g) > 1e-10 else g
            identity_proj = np.abs(np.trace(g_normalized)) / dim
            print(f"    Gen {idx}: |Tr(g/||g||)|/d = {identity_proj:.6f}"
                  f" (1.0 = pure U(1))")

        # Decomposition summary
        n_u1 = len(r["u1_indices"])
        n_su = len(r["su_indices"])
        K_su = r["K_eigs"][r["K_eigs"] < -1e-6] if n_su > 0 else []
        K_u1 = r["K_eigs"][np.abs(r["K_eigs"]) < 1e-6]
        print(f"    Decomposition: u(1)^{n_u1} ⊕ su({dim})^?")
        print(f"    K_su eigenvalues: {K_su}")
        print(f"    K_u1 eigenvalues: {K_u1}")

    return results


# ===================================================================
# EXPERIMENT C: n_gen=d², general hermitian, K_target=0
# ===================================================================
def experiment_C(dim=2):
    """n_gen=d², dim=d, general hermitian, K_target=0 (abelian target).
    What happens when we target the abelian Killing form?
    Does the optimizer find an abelian subalgebra?
    """
    n_gen = dim * dim

    print("\n" + "=" * 80)
    print(f"  EXPERIMENT C: n_gen={n_gen}, dim={dim}, general hermitian, K_target=0")
    print(f"  Targeting abelian algebra (all generators commute)")
    print("=" * 80)

    seeds = [42, 137, 256]

    results = evolve_generators_generic(
        n_gen=n_gen, dim=dim, steps=15000, lr=0.001,
        alpha=5.0, beta=1.0, killing_target=0.0,
        batch_size=len(seeds), seeds=seeds,
        verbose=True, convergence_threshold=1e-8,
        param_to_gen_fn=lambda p, ng: params_to_hermitian(p, ng, dim),
        gen_to_param_fn=hermitian_to_params,
        random_param_fn=lambda bs, ng, device=DEVICE, seeds=None: random_params_hermitian(
            bs, ng, dim, device, seeds),
        param_to_gen_kwargs={},
    )

    for i, (gens_np, hist) in enumerate(results):
        r = analyze_result(list(gens_np), f"Exp C dim={dim} seed={seeds[i]}")

        # Check commutativity
        max_comm = 0
        for a in range(n_gen):
            for b in range(a + 1, n_gen):
                comm = np.linalg.norm(commutator_np(gens_np[a], gens_np[b]))
                max_comm = max(max_comm, comm)
        print(f"    Max commutator norm: {max_comm:.2e}")
        print(f"    Abelian: {max_comm < 1e-6}")

        # Check if generators are simultaneously diagonalizable
        # (abelian hermitian matrices are simultaneously diagonalizable)
        all_eigs = []
        for g in gens_np:
            eigs = np.sort(np.linalg.eigvalsh(g.astype(np.complex128)))
            all_eigs.append(eigs)
        print(f"    Generator eigenvalues:")
        for j, eigs in enumerate(all_eigs):
            print(f"      Gen {j}: {eigs.real}")

    return results


# ===================================================================
# EXPERIMENT D: U(2) vs SU(2) direct comparison
# ===================================================================
def experiment_D():
    """Side-by-side: n_gen=3 (SU(2)) vs n_gen=4 (U(2)), both dim=2.
    The 4th generator should be the U(1) factor.
    """
    print("\n" + "=" * 80)
    print("  EXPERIMENT D: U(2) vs SU(2) side-by-side comparison")
    print("  SU(2): n_gen=3, dim=2, traceless hermitian")
    print("  U(2):  n_gen=4, dim=2, general hermitian")
    print("=" * 80)

    seeds = [42, 137, 256]

    # SU(2)
    print("\n  --- SU(2): n_gen=3 traceless ---")
    su2_results = evolve_generators_generic(
        n_gen=3, dim=2, steps=12000, lr=0.001,
        alpha=5.0, beta=1.0, killing_target=-1.0,
        batch_size=len(seeds), seeds=seeds,
        verbose=True, convergence_threshold=1e-6,
    )

    su2_best = min(su2_results, key=lambda r: r[1]["T_total"][-1])
    r_su2 = analyze_result(list(su2_best[0]), "SU(2) best")

    # U(2)
    print("\n  --- U(2): n_gen=4 general hermitian ---")
    dim = 2
    u2_results = evolve_generators_generic(
        n_gen=4, dim=2, steps=12000, lr=0.001,
        alpha=5.0, beta=1.0, killing_target=-1.0,
        batch_size=len(seeds), seeds=seeds,
        verbose=True, convergence_threshold=1e-6,
        param_to_gen_fn=lambda p, ng: params_to_hermitian(p, ng, dim),
        gen_to_param_fn=hermitian_to_params,
        random_param_fn=lambda bs, ng, device=DEVICE, seeds=None: random_params_hermitian(
            bs, ng, dim, device, seeds),
        param_to_gen_kwargs={},
    )

    u2_best = min(u2_results, key=lambda r: r[1]["T_total"][-1])
    r_u2 = analyze_result(list(u2_best[0]), "U(2) best")

    # Compare
    print("\n  --- Comparison ---")
    print(f"  SU(2) K_eigs: {r_su2['K_eigs']}")
    print(f"  U(2)  K_eigs: {r_u2['K_eigs']}")
    print(f"  SU(2) closure: {r_su2['closure']:.1%}")
    print(f"  U(2)  closure: {r_u2['closure']:.1%}")

    # Check if the extra generator in U(2) commutes with the SU(2) part
    u2_gens = list(u2_best[0])
    u1_idx = r_u2["u1_indices"]
    su_idx = r_u2["su_indices"]
    if u1_idx:
        print(f"\n  U(1) generator(s): {u1_idx}")
        for ui in u1_idx:
            max_cross = 0
            for si in su_idx:
                c = np.linalg.norm(commutator_np(u2_gens[ui], u2_gens[si]))
                max_cross = max(max_cross, c)
            print(f"    Gen {ui}: max |[U(1), SU]| = {max_cross:.2e}")
            # Is it ∝ I?
            g = u2_gens[ui]
            g_n = g / np.linalg.norm(g) if np.linalg.norm(g) > 1e-10 else g
            proj = abs(np.trace(g_n)) / dim
            print(f"    Gen {ui}: |Tr(g/||g||)|/d = {proj:.6f} (1.0 = pure identity)")

    return su2_results, u2_results


# ===================================================================
# EXPERIMENT E: U(3) = U(1) ⊕ SU(3)
# ===================================================================
def experiment_E():
    """n_gen=9, dim=3, general hermitian.
    Expected: u(3) = u(1) ⊕ su(3), with 8 su(3) generators + 1 u(1).
    """
    dim = 3
    n_gen = dim * dim  # 9

    print("\n" + "=" * 80)
    print(f"  EXPERIMENT E: n_gen={n_gen}, dim={dim}, general hermitian")
    print(f"  Expected: u(3) = u(1) ⊕ su(3)")
    print(f"  8 generators with K≈-0.962 (su(3) part)")
    print(f"  1 generator with K=0 (u(1) part)")
    print("=" * 80)

    seeds = [42, 137, 256, 314, 999]

    results = evolve_generators_generic(
        n_gen=n_gen, dim=dim, steps=15000, lr=0.001,
        alpha=5.0, beta=1.0, killing_target=-1.0,
        batch_size=len(seeds), seeds=seeds,
        verbose=True, convergence_threshold=1e-6,
        param_to_gen_fn=lambda p, ng: params_to_hermitian(p, ng, dim),
        gen_to_param_fn=hermitian_to_params,
        random_param_fn=lambda bs, ng, device=DEVICE, seeds=None: random_params_hermitian(
            bs, ng, dim, device, seeds),
        param_to_gen_kwargs={},
    )

    for i, (gens_np, hist) in enumerate(results):
        r = analyze_result(list(gens_np), f"Exp E seed={seeds[i]}")

        # Detailed K spectrum
        K_nonzero = r["K_eigs"][np.abs(r["K_eigs"]) > 1e-4]
        K_zero = r["K_eigs"][np.abs(r["K_eigs"]) <= 1e-4]
        print(f"    K nonzero eigenvalues ({len(K_nonzero)}): "
              f"mean={np.mean(K_nonzero):.6f}" if len(K_nonzero) > 0 else "")
        print(f"    K zero eigenvalues ({len(K_zero)}): {K_zero}")

        # Identify U(1) generator and check it's ∝ I
        for idx in r["u1_indices"]:
            g = gens_np[idx]
            g_n = g / np.linalg.norm(g) if np.linalg.norm(g) > 1e-10 else g
            proj = abs(np.trace(g_n)) / dim
            print(f"    Gen {idx}: |Tr(g/||g||)|/d = {proj:.6f}")

    return results


# ===================================================================
# MAIN
# ===================================================================
def main():
    print("=" * 80)
    print("  §13: U(1) ABELIAN EMERGENCE TESTS")
    print("=" * 80)
    print(f"  Device: {get_device_info()}")
    print()

    t0 = time.time()

    # Exp A: Single generator, SU param
    experiment_A()

    # Exp D: U(2) vs SU(2) comparison
    experiment_D()

    # Exp E: U(3) = U(1) ⊕ SU(3)
    experiment_E()

    # Exp C: Abelian target
    experiment_C(dim=2)

    dt = time.time() - t0
    print(f"\n\n{'='*80}")
    print(f"  Total time: {dt:.1f}s")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
