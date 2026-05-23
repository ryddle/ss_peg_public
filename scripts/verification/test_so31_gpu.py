"""
Test §12.4: SO(3,1) — Lorentz algebra emergence from the variational principle.

so(3,1) is the Lie algebra of the Lorentz group, with 6 generators in dim=4:
  - 3 rotation generators J₁, J₂, J₃ (compact, forming so(3) subalgebra)
  - 3 boost generators K₁, K₂, K₃ (non-compact)

Defining condition: T^T η + η T = 0, where η = diag(-1,+1,+1,+1).
This is a 6-parameter family of real 4×4 matrices per generator.

Killing form signature: (3+, 3−) — 3 compact (rotations) + 3 non-compact (boosts).
Compare with so(4): all 6 compact, Killing signature (0+, 6−).

Experiments:
  A. η-preserving param + indefinite Killing target [-1,-1,-1,+1,+1,+1]
     → should converge to so(3,1) Lorentz algebra
  B. Antisymmetric param + compact target K=-I (so(4) control)
     → should converge to so(4) compact
  C. η-preserving param + compact target K=-I
     → should FAIL (indefinite K can't match -I)

Usage:
    python test_so31_gpu.py
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
    params_to_antisymmetric,
    antisymmetric_to_params,
    random_params_so,
    DTYPE, FDTYPE, DEVICE,
    batched_killing_metric,
    batched_non_commutativity_tension,
    batched_norm_tension,
)


# ===================================================================
# Minkowski metric η = diag(-1,+1,+1,+1)
# ===================================================================
ETA = torch.tensor([[-1., 0., 0., 0.],
                     [ 0., 1., 0., 0.],
                     [ 0., 0., 1., 0.],
                     [ 0., 0., 0., 1.]], dtype=FDTYPE, device=DEVICE)

N_GEN = 6   # so(3,1) has 6 generators
DIM = 4     # 4×4 real matrices


# ===================================================================
# η-preserving parameterization for so(3,1)
# ===================================================================
# Condition: T^T η + η T = 0  (Lorentz algebra condition)
#
# For η = diag(-1,+1,+1,+1), writing T as 4×4 with entries T_{μν}:
#   η_{μμ} T_{μν} + η_{νν} T_{νμ} = 0  for all μ,ν (no sum)
#
# This means:
#   - If η_{μμ} = η_{νν} (same sign): T_{μν} = -T_{νμ} (antisymmetric)
#     → spatial-spatial pairs (i,j) for i,j ∈ {1,2,3}
#   - If η_{μμ} ≠ η_{νν} (different sign): T_{μν} = +T_{νμ} (symmetric)
#     → time-space pairs (0,i) for i ∈ {1,2,3}
#
# Independent parameters: 6 = 3 (rotations: T_{12}, T_{13}, T_{23}) + 3 (boosts: T_{01}=T_{10}, T_{02}=T_{20}, T_{03}=T_{30})
# Plus: T_{00} = 0 (diagonal must vanish by the η-condition with μ=ν)

def params_to_so31(params, n_gen):
    """Convert 6 real parameters per generator to η-preserving 4×4 real matrices.

    The parameterization enforces T^T η + η T = 0 structurally.

    params: (B, n_gen, 6) real/complex tensor (only real part used)
    Returns: (B, n_gen, 4, 4) complex tensor (purely real)

    Parameter mapping:
      p[0] = T_{01} = T_{10}  (boost x)    → symmetric because η₀₀ ≠ η₁₁
      p[1] = T_{02} = T_{20}  (boost y)    → symmetric
      p[2] = T_{03} = T_{30}  (boost z)    → symmetric
      p[3] = T_{12}; T_{21} = -T_{12}  (rotation z) → antisymmetric
      p[4] = T_{13}; T_{31} = -T_{13}  (rotation y) → antisymmetric
      p[5] = T_{23}; T_{32} = -T_{23}  (rotation x) → antisymmetric
    """
    p = params.real if params.is_complex() else params
    B = p.shape[0]
    T = torch.zeros(B, n_gen, 4, 4, dtype=DTYPE, device=params.device)

    # Boost parameters (symmetric in time-space block)
    b1 = p[..., 0]  # T_{01}
    b2 = p[..., 1]  # T_{02}
    b3 = p[..., 2]  # T_{03}

    # Rotation parameters (antisymmetric in spatial block)
    r3 = p[..., 3]  # T_{12} (rotation around z)
    r2 = p[..., 4]  # T_{13} (rotation around y)
    r1 = p[..., 5]  # T_{23} (rotation around x)

    # Time-space block: symmetric (boosts)
    T[..., 0, 1] = b1;  T[..., 1, 0] = b1
    T[..., 0, 2] = b2;  T[..., 2, 0] = b2
    T[..., 0, 3] = b3;  T[..., 3, 0] = b3

    # Spatial block: antisymmetric (rotations)
    T[..., 1, 2] = r3;  T[..., 2, 1] = -r3
    T[..., 1, 3] = r2;  T[..., 3, 1] = -r2
    T[..., 2, 3] = r1;  T[..., 3, 2] = -r1

    return T


def so31_to_params(generators):
    """Inverse: extract 6 real parameters from η-preserving generators."""
    B, n_gen = generators.shape[:2]
    p = torch.zeros(B, n_gen, 6, dtype=FDTYPE, device=generators.device)

    # Boosts
    p[..., 0] = generators[..., 0, 1].real  # T_{01}
    p[..., 1] = generators[..., 0, 2].real  # T_{02}
    p[..., 2] = generators[..., 0, 3].real  # T_{03}

    # Rotations
    p[..., 3] = generators[..., 1, 2].real  # T_{12}
    p[..., 4] = generators[..., 1, 3].real  # T_{13}
    p[..., 5] = generators[..., 2, 3].real  # T_{23}

    return p


def random_params_so31(batch_size, n_gen, device=DEVICE, seeds=None):
    """Random parameters for so(3,1): 6 real params per generator."""
    if seeds is not None:
        all_params = []
        for s in seeds:
            rng = np.random.default_rng(s)
            p = rng.standard_normal((n_gen, 6))
            all_params.append(p)
        params_np = np.stack(all_params)
        return torch.tensor(params_np, dtype=FDTYPE, device=device)
    return torch.randn(batch_size, n_gen, 6, dtype=FDTYPE, device=device)


# ===================================================================
# Verify the parameterization constraint T^T η + η T = 0
# ===================================================================
def verify_eta_constraint(generators_np, tol=1e-8):
    """Check that all generators satisfy the Lorentz algebra condition."""
    eta = np.diag([-1., 1., 1., 1.])
    violations = []
    for i, T in enumerate(generators_np):
        residual = T.real.T @ eta + eta @ T.real
        err = np.max(np.abs(residual))
        violations.append(err)
        if err > tol:
            return False, violations
    return True, violations


# ===================================================================
# Custom optimizer for indefinite Killing target (SO(3,1))
# ===================================================================
def evolve_so31_indefinite(
    n_gen=6,
    steps=30000,
    lr=0.001,
    alpha=5.0,
    beta=1.0,
    batch_size=5,
    seeds=None,
    verbose=True,
    convergence_threshold=1e-4,
):
    """Optimizer for so(3,1) with indefinite Killing target.

    Target: K eigenvalues sorted → [-1,-1,-1,+1,+1,+1]
    (3 compact directions + 3 non-compact directions)

    Loss = T_nc + alpha * (T_k_sig + T_k_off) + beta * T_norm
    """
    device = DEVICE

    if seeds is None:
        seeds = list(range(batch_size))

    params = random_params_so31(batch_size, n_gen, device, seeds)

    # Normalize
    with torch.no_grad():
        gens = params_to_so31(params, n_gen)
        norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
        gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
        params = so31_to_params(gens)
    params = params.clone().detach().requires_grad_(True)

    optimizer = torch.optim.Adam([params], lr=lr)

    # Indefinite target: 3 negative (compact/rotations), 3 positive (non-compact/boosts)
    # Sorted ascending: [-1, -1, -1, +1, +1, +1]
    K_target_sorted = torch.tensor([-1., -1., -1., 1., 1., 1.],
                                    dtype=FDTYPE, device=device)

    prev_loss = torch.zeros(batch_size, device=device)
    active = torch.ones(batch_size, dtype=torch.bool, device=device)

    results_per_seed = [{"loss": [], "T_nc": [], "T_k_sig": [], "T_k_off": [],
                          "K_eigs": []} for _ in range(batch_size)]

    for step in range(steps):
        optimizer.zero_grad()
        generators = params_to_so31(params, n_gen)

        # Non-commutativity tension (closure driver)
        T_nc = batched_non_commutativity_tension(generators)

        # Killing metric
        K = batched_killing_metric(generators)  # (B, 6, 6)

        # Eigenvalue-based Killing loss (signature matching)
        K_real = K.real if K.is_complex() else K
        K_eigs = torch.linalg.eigvalsh(K_real)  # (B, 6) ascending
        T_k_sig = ((K_eigs - K_target_sorted.unsqueeze(0)) ** 2).sum(dim=-1)

        # Off-diagonal penalty
        K_offdiag_mask = ~torch.eye(n_gen, dtype=torch.bool, device=device).unsqueeze(0)
        T_k_off = (K * K_offdiag_mask.float()).pow(2).sum(dim=(-2, -1))

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
                    results_per_seed[b]["T_nc"].append(T_nc[b].item())
                    results_per_seed[b]["T_k_sig"].append(T_k_sig[b].item())
                    results_per_seed[b]["T_k_off"].append(T_k_off[b].item())
                    results_per_seed[b]["K_eigs"].append(K_eigs[b].cpu().numpy().copy())

            converged = total < convergence_threshold
            if step > 0:
                plateau = (total - prev_loss).abs() < 1e-14
                newly_done = (converged | plateau) & active
                active = active & ~newly_done
            prev_loss = total.clone()

            if verbose and step % 2000 == 0:
                n_act = active.sum().item()
                mean_l = total[active].mean().item() if n_act > 0 else 0
                best_Tnc = T_nc.min().item()
                print(f"  Step {step}: active={n_act}/{batch_size} "
                      f"mean_loss={mean_l:.6f} best_Tnc={best_Tnc:.6f}")

            if not active.any():
                if verbose:
                    print(f"  All converged at step {step}")
                break

    with torch.no_grad():
        final = params_to_so31(params, n_gen).cpu().numpy()

    return [(final[b], results_per_seed[b]) for b in range(batch_size)]


# ===================================================================
# Analysis (extended from sl2r test)
# ===================================================================
def analyze_result(generators, label, verbose=True):
    """Analyze a set of 4×4 generators for Lorentz algebra properties."""
    n_gen = len(generators)

    # Basic properties
    is_real = all(np.max(np.abs(g.imag)) < 1e-10 for g in generators)
    traces = [np.trace(g).real for g in generators]
    is_traceless = all(abs(t) < 1e-8 for t in traces)

    # η-preservation check
    eta_ok, eta_errs = verify_eta_constraint(generators)

    # Killing metric
    K = killing_metric_np(generators)
    K_offdiag = np.abs(K - np.diag(np.diag(K))).max()

    # Eigenvalues and signature
    eigs = np.linalg.eigvalsh(K.real)
    n_pos = int(np.sum(eigs > 1e-6))
    n_neg = int(np.sum(eigs < -1e-6))
    n_zero = n_gen - n_pos - n_neg

    # Closure test
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

    # Identify rotation vs boost generators via Killing sign
    # Diagonal of K: negative → compact (rotation), positive → non-compact (boost)
    K_diag = np.diag(K).real

    result = {
        "label": label,
        "n_gen": n_gen,
        "is_real": is_real,
        "is_traceless": is_traceless,
        "eta_ok": eta_ok,
        "eta_errs": eta_errs,
        "K": K.real,
        "K_diag": K_diag,
        "K_offdiag": K_offdiag,
        "K_eigs": eigs,
        "signature": (n_pos, n_zero, n_neg),
        "closure": closure,
        "f_rms": f_rms,
    }

    if verbose:
        tag = "PASS" if closure > 0.99 else "FAIL"
        print(f"\n  [{tag}] {label}:")
        print(f"    Real: {is_real}, Traceless: {is_traceless}")
        print(f"    η-preserving: {eta_ok} (max err: {max(eta_errs):.2e})")
        print(f"    K eigenvalues: {np.sort(eigs)}")
        print(f"    K signature:   ({n_pos}+, {n_zero}0, {n_neg}−)")
        print(f"    K_offdiag max: {K_offdiag:.2e}")
        print(f"    Closure:       {closes}/{total_pairs} ({closure:.1%})")
        print(f"    f_rms:         {f_rms:.6f}")

        # Classify generators
        if closure > 0.99 and n_pos == 3 and n_neg == 3:
            print(f"\n    *** so(3,1) Lorentz algebra identified! ***")
            # Attempt to identify J vs K generators
            print(f"    K_diag = {K_diag}")
            rot_idx = [i for i in range(n_gen) if K_diag[i] < 0]
            boost_idx = [i for i in range(n_gen) if K_diag[i] > 0]
            print(f"    Rotation generators (K<0): indices {rot_idx}")
            print(f"    Boost generators    (K>0): indices {boost_idx}")

    return result


def compare_with_canonical(generators, verbose=True):
    """Compare emergent generators with canonical so(3,1) basis.

    Canonical basis (physicist convention):
      J₁ = rotation in (2,3) plane: T_{23}=-T_{32}=1
      J₂ = rotation in (1,3) plane: T_{13}=-T_{31}=1
      J₃ = rotation in (1,2) plane: T_{12}=-T_{21}=1
      K₁ = boost in x: T_{01}=T_{10}=1
      K₂ = boost in y: T_{02}=T_{20}=1
      K₃ = boost in z: T_{03}=T_{30}=1
    """
    # Build canonical generators
    J = [np.zeros((4, 4), dtype=complex) for _ in range(3)]
    K_boost = [np.zeros((4, 4), dtype=complex) for _ in range(3)]

    # J₁: rotation in (2,3) plane
    J[0][2, 3] = 1; J[0][3, 2] = -1
    # J₂: rotation in (1,3) plane
    J[1][1, 3] = 1; J[1][3, 1] = -1
    # J₃: rotation in (1,2) plane
    J[2][1, 2] = 1; J[2][2, 1] = -1

    # K₁: boost in x
    K_boost[0][0, 1] = 1; K_boost[0][1, 0] = 1
    # K₂: boost in y
    K_boost[1][0, 2] = 1; K_boost[1][2, 0] = 1
    # K₃: boost in z
    K_boost[2][0, 3] = 1; K_boost[2][3, 0] = 1

    canonical = J + K_boost  # [J₁, J₂, J₃, K₁, K₂, K₃]
    labels_can = ["J₁", "J₂", "J₃", "K₁", "K₂", "K₃"]

    if verbose:
        print("\n    Canonical so(3,1) comparison:")

    # Check if emergent generators span the same space
    can_flat = np.array([g.flatten() for g in canonical]).T  # (16, 6)
    span_ok = True
    for i, g in enumerate(generators):
        coeffs, _, _, _ = np.linalg.lstsq(can_flat, g.flatten(), rcond=None)
        residual = np.linalg.norm(can_flat @ coeffs - g.flatten())
        in_span = residual < 1e-3
        span_ok = span_ok and in_span
        if verbose:
            coeffs_str = ", ".join(f"{c:.3f}" for c in coeffs.real)
            print(f"      Gen {i}: [{coeffs_str}] residual={residual:.2e} {'✓' if in_span else '✗'}")

    if verbose:
        print(f"    All emergent gens in canonical span: {span_ok}")

    # Verify canonical commutation relations
    if verbose:
        print("\n    Canonical commutation check:")
        # [J_i, J_j] = ε_ijk J_k  (compact subalgebra so(3))
        for i in range(3):
            for j in range(i+1, 3):
                C = commutator_np(canonical[i], canonical[j])
                k = 3 - i - j  # the third index
                sign = 1 if (i, j, k) in [(0,1,2), (1,2,0), (2,0,1)] else -1
                expected = sign * canonical[k]
                err = np.linalg.norm(C - expected)
                print(f"      [J{i+1},J{j+1}] = {'+'if sign>0 else '-'}J{k+1}  err={err:.2e}")

        # [K_i, K_j] = -ε_ijk J_k  (boosts don't close into boosts)
        for i in range(3):
            for j in range(i+1, 3):
                C = commutator_np(canonical[3+i], canonical[3+j])
                k = 3 - i - j
                sign = 1 if (i, j, k) in [(0,1,2), (1,2,0), (2,0,1)] else -1
                expected = -sign * canonical[k]  # note the minus!
                err = np.linalg.norm(C - expected)
                print(f"      [K{i+1},K{j+1}] = {'−'if sign>0 else '+'}J{k+1}  err={err:.2e}")

        # [J_i, K_j] = ε_ijk K_k
        for i in range(3):
            for j in range(3):
                if i == j:
                    continue
                C = commutator_np(canonical[i], canonical[3+j])
                k = 3 - i - j
                if k < 0 or k > 2:
                    continue
                sign = 1 if (i, j, k) in [(0,1,2), (1,2,0), (2,0,1)] else -1
                expected = sign * canonical[3+k]
                err = np.linalg.norm(C - expected)
                print(f"      [J{i+1},K{j+1}] = {'+'if sign>0 else '-'}K{k+1}  err={err:.2e}")

    return span_ok


# ===================================================================
# so(4) helper for Experiment B
# ===================================================================
def random_params_so4(batch_size, n_gen, device=DEVICE, seeds=None):
    """Random parameters for so(4): antisymmetric 4×4 → 6 params per generator."""
    return random_params_so(batch_size, n_gen, DIM, device, seeds)


def params_to_so4(params, n_gen):
    """so(4) = antisymmetric 4×4 real matrices."""
    return params_to_antisymmetric(params, n_gen, DIM)


def so4_to_params(generators):
    """Extract parameters from antisymmetric generators."""
    return antisymmetric_to_params(generators)


# ===================================================================
# Main
# ===================================================================
def main():
    print("=" * 78)
    print("  §12.4: SO(3,1) LORENTZ ALGEBRA EMERGENCE TEST")
    print("=" * 78)
    print(f"  Device: {get_device_info()}")
    print(f"  so(3,1): 6 generators in dim=4 (real), K signature (3+, 3−)")
    print(f"  so(4):   6 generators in dim=4 (real), K signature (0+, 6−)")
    print()

    n_seeds = 5

    # ---------------------------------------------------------------
    # Experiment A: η-preserving param + indefinite target → so(3,1)
    # ---------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  EXPERIMENT A: η-preserving param + indefinite Killing target")
    print("  Target K_eigs sorted: [-1, -1, -1, +1, +1, +1]")
    print("  Expected: should converge to so(3,1) Lorentz algebra")
    print("=" * 78)

    t0 = time.time()
    results_a = evolve_so31_indefinite(
        n_gen=N_GEN, steps=30000, lr=0.001,
        alpha=5.0, beta=1.0,
        batch_size=n_seeds, seeds=list(range(n_seeds)),
        verbose=True, convergence_threshold=1e-4,
    )
    time_a = time.time() - t0

    best_a = None
    best_closure_a = -1
    for s, (gens, hist) in enumerate(results_a):
        r = analyze_result(list(gens), f"Exp A seed={s}")
        if r["closure"] > best_closure_a:
            best_closure_a = r["closure"]
            best_a = (list(gens), r)

    if best_a is not None and best_a[1]["closure"] > 0.99:
        compare_with_canonical(best_a[0])

    # ---------------------------------------------------------------
    # Experiment B: so(4) antisymmetric param + compact target (control)
    # ---------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  EXPERIMENT B: Antisymmetric param + compact target (so(4) control)")
    print("  Target K = -I (all 6 compact directions)")
    print("  Expected: should converge to so(4)")
    print("=" * 78)

    t0 = time.time()
    results_b = evolve_generators_generic(
        n_gen=N_GEN, dim=DIM, steps=20000, lr=0.001,
        alpha=5.0, beta=1.0, killing_target=-1.0,
        batch_size=n_seeds, seeds=list(range(n_seeds)),
        verbose=True, convergence_threshold=1e-4,
        param_to_gen_fn=params_to_so4,
        gen_to_param_fn=so4_to_params,
        random_param_fn=random_params_so4,
        param_to_gen_kwargs={},
    )
    time_b = time.time() - t0

    for s, (gens, _) in enumerate(results_b):
        r = analyze_result(list(gens), f"Exp B seed={s}")
        # Verify so(4) has all-negative K
        if r["closure"] > 0.99 and r["signature"] == (0, 0, 6):
            print(f"    ✓ so(4) compact algebra confirmed")

    # ---------------------------------------------------------------
    # Experiment C: η-preserving param + compact target → should fail
    # ---------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  EXPERIMENT C: η-preserving param + compact target K=-I")
    print("  Expected: FAIL (indefinite K cannot match negative-definite target)")
    print("=" * 78)

    t0 = time.time()
    results_c = evolve_generators_generic(
        n_gen=N_GEN, dim=DIM, steps=20000, lr=0.001,
        alpha=5.0, beta=1.0, killing_target=-1.0,
        batch_size=n_seeds, seeds=list(range(n_seeds)),
        verbose=True, convergence_threshold=1e-4,
        param_to_gen_fn=params_to_so31,
        gen_to_param_fn=so31_to_params,
        random_param_fn=random_params_so31,
        param_to_gen_kwargs={},
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
    print(f"\n  {'Experiment':40s} {'Closure':>10s} {'K signature':>15s} {'η-pres':>8s} {'Time':>8s}")
    print("  " + "-" * 85)

    for label, results, t in [
        ("A: so(3,1) η-pres + indef target", results_a, time_a),
        ("B: so(4) antisymm + compact target", results_b, time_b),
        ("C: so(3,1) η-pres + compact target", results_c, time_c),
    ]:
        for s, (gens, _) in enumerate(results):
            r = analyze_result(list(gens), f"{label} s{s}", verbose=False)
            sig = r["signature"]
            eta_str = "✓" if r["eta_ok"] else "✗"
            print(f"  {label + f' s{s}':40s} {r['closure']:>9.1%} "
                  f"  ({sig[0]}+,{sig[1]}0,{sig[2]}−) "
                  f"  {eta_str:>5s} {t/n_seeds:>7.1f}s")

    print("  " + "-" * 85)
    print()
    print("  Key predictions:")
    print("  - Exp A: Closure ≈ 100%, K signature (3+,0,3−) → so(3,1) ✓")
    print("  - Exp B: Closure ≈ 100%, K signature (0+,0,6−) → so(4) control ✓")
    print("  - Exp C: Closure < 100%, frustrated (wrong target for constraint) ✗")
    print()
    print("  Physics significance:")
    print("  - so(3,1) emergence from η-constraint + indefinite Killing target")
    print("    connects SS-PEG directly to relativistic spacetime symmetry.")
    print("  - The variational principle (T_nc, K, T_norm) with Minkowski η")
    print("    selects the Lorentz algebra as the unique non-compact 6D algebra in dim=4.")
    print("  - Killing signature discriminates: (3+,3−) → spacetime, (0+,6−) → internal.")


if __name__ == "__main__":
    main()
