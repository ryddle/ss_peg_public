"""
Test §14: SU(2,2) — Conformal algebra emergence from the variational principle.

su(2,2) ≅ so(4,2) — the conformal algebra of 3+1 Minkowski spacetime.

15 generators in dim=4 (complex traceless η-preserving):
  - 7 compact generators forming s(u(2)⊕u(2)) — maximal compact subalgebra
  - 8 non-compact generators (translations, special conformal, dilation)

Defining condition: T†η + ηT = 0, Tr(T) = 0, where η = diag(-1,-1,+1,+1).
This is a 15-parameter family per generator (complex 4×4 traceless η-anti-hermitian).

Killing form signature: (8+, 0, 7−) — 7 compact + 8 non-compact.
Compare with su(4): all 15 compact, Killing signature (0+, 0, 15−).

Physical significance:
  - SU(2,2) ≅ SO(4,2): conformal group of 3+1 spacetime
  - Contains SO(3,1) Lorentz as subalgebra
  - Contains the Poincaré algebra as subalgebra
  - Symmetry of massless field equations (Maxwell, Dirac for m=0)
  - Foundation of AdS₅/CFT₄ correspondence
  - 15 generators: 6 Lorentz (M_μν) + 4 translations (P_μ)
    + 4 special conformal (K_μ) + 1 dilation (D)

Experiments:
  A. η(2,2)-preserving param + indefinite Killing target [-1]*7 + [+1]*8
     → should converge to su(2,2) conformal algebra
  B. su(4) traceless hermitian + compact target K=-I
     → should converge to su(4) compact (control)
  C. η(2,2)-preserving param + compact target K=-I
     → should FAIL (indefinite K can't match negative-definite target)

Usage:
    python test_su22_gpu.py
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
    params_to_generators,
    generators_to_params,
    random_params,
    DTYPE, FDTYPE, DEVICE,
    batched_killing_metric,
    batched_non_commutativity_tension,
    batched_norm_tension,
)


# ===================================================================
# η metric for su(2,2): η = diag(-1,-1,+1,+1)
# ===================================================================
ETA = torch.tensor([[-1., 0., 0., 0.],
                     [ 0.,-1., 0., 0.],
                     [ 0., 0., 1., 0.],
                     [ 0., 0., 0., 1.]], dtype=FDTYPE, device=DEVICE)

N_GEN = 15   # su(2,2) has 15 generators
DIM = 4      # 4×4 complex matrices
N_PARAMS = 15  # 15 real parameters per generator


# ===================================================================
# η(2,2)-preserving parameterization for su(2,2)
# ===================================================================
# Condition: T†η + ηT = 0, Tr(T) = 0
# where η = diag(-1,-1,+1,+1)
#
# Component-wise: T̄_{ji} η_{jj} + η_{ii} T_{ij} = 0
#                 ⟹  T_{ij} = -(η_{jj}/η_{ii}) T̄_{ji}
#
# Same-sign blocks (εᵢ = εⱼ): T_{ij} = -T̄_{ji}  (anti-hermitian pairs)
#   → (0,1): ε₀=ε₁=-1       and  (2,3): ε₂=ε₃=+1
# Cross-sign blocks (εᵢ ≠ εⱼ): T_{ij} = +T̄_{ji}  (hermitian pairs)
#   → (0,2), (0,3), (1,2), (1,3)
# Diagonal: T_{ii} = -T̄_{ii} → purely imaginary
#           Tr(T)=0: a₁+a₂+a₃+a₄=0 → 3 free real params
#
# Parameter count: 3 (diag) + 2×2 (same-sign) + 4×2 (cross-sign) = 15 ✓

def params_to_su22(params, n_gen):
    """Convert 15 real parameters per generator to η(2,2)-preserving
    traceless 4×4 complex matrices.

    Enforces T†η + ηT = 0 and Tr(T) = 0 structurally.

    params: (B, n_gen, 15) real tensor
    Returns: (B, n_gen, 4, 4) complex tensor

    Parameter mapping:
      p[0-2]:   diagonal imaginary (a₁, a₂, a₃; a₄ = -(a₁+a₂+a₃))
      p[3-4]:   (0,1) re,im — same-sign, anti-hermitian
      p[5-6]:   (2,3) re,im — same-sign, anti-hermitian
      p[7-8]:   (0,2) re,im — cross-sign, hermitian
      p[9-10]:  (0,3) re,im — cross-sign, hermitian
      p[11-12]: (1,2) re,im — cross-sign, hermitian
      p[13-14]: (1,3) re,im — cross-sign, hermitian
    """
    p = params.real if params.is_complex() else params
    B = p.shape[0]
    T = torch.zeros(B, n_gen, 4, 4, dtype=DTYPE, device=params.device)
    zero = torch.zeros_like(p[..., 0])

    # --- Diagonal: purely imaginary, traceless ---
    a1, a2, a3 = p[..., 0], p[..., 1], p[..., 2]
    a4 = -(a1 + a2 + a3)
    T[..., 0, 0] = torch.complex(zero, a1)
    T[..., 1, 1] = torch.complex(zero, a2)
    T[..., 2, 2] = torch.complex(zero, a3)
    T[..., 3, 3] = torch.complex(zero, a4)

    # --- Same-sign pairs: anti-hermitian (T_{ij} = -T̄_{ji}) ---
    # Pair (0,1): ε₀ = ε₁ = -1
    re01, im01 = p[..., 3], p[..., 4]
    T[..., 0, 1] = torch.complex(re01, im01)
    T[..., 1, 0] = torch.complex(-re01, im01)

    # Pair (2,3): ε₂ = ε₃ = +1
    re23, im23 = p[..., 5], p[..., 6]
    T[..., 2, 3] = torch.complex(re23, im23)
    T[..., 3, 2] = torch.complex(-re23, im23)

    # --- Cross-sign pairs: hermitian (T_{ij} = T̄_{ji}) ---
    # Pair (0,2): ε₀=-1, ε₂=+1
    re02, im02 = p[..., 7], p[..., 8]
    T[..., 0, 2] = torch.complex(re02, im02)
    T[..., 2, 0] = torch.complex(re02, -im02)

    # Pair (0,3): ε₀=-1, ε₃=+1
    re03, im03 = p[..., 9], p[..., 10]
    T[..., 0, 3] = torch.complex(re03, im03)
    T[..., 3, 0] = torch.complex(re03, -im03)

    # Pair (1,2): ε₁=-1, ε₂=+1
    re12, im12 = p[..., 11], p[..., 12]
    T[..., 1, 2] = torch.complex(re12, im12)
    T[..., 2, 1] = torch.complex(re12, -im12)

    # Pair (1,3): ε₁=-1, ε₃=+1
    re13, im13 = p[..., 13], p[..., 14]
    T[..., 1, 3] = torch.complex(re13, im13)
    T[..., 3, 1] = torch.complex(re13, -im13)

    return T


def su22_to_params(generators):
    """Inverse: extract 15 real parameters from su(2,2) generators."""
    B, n_gen = generators.shape[:2]
    p = torch.zeros(B, n_gen, N_PARAMS, dtype=FDTYPE, device=generators.device)

    # Diagonal
    p[..., 0] = generators[..., 0, 0].imag
    p[..., 1] = generators[..., 1, 1].imag
    p[..., 2] = generators[..., 2, 2].imag

    # Same-sign pairs
    p[..., 3] = generators[..., 0, 1].real
    p[..., 4] = generators[..., 0, 1].imag
    p[..., 5] = generators[..., 2, 3].real
    p[..., 6] = generators[..., 2, 3].imag

    # Cross-sign pairs
    p[..., 7]  = generators[..., 0, 2].real
    p[..., 8]  = generators[..., 0, 2].imag
    p[..., 9]  = generators[..., 0, 3].real
    p[..., 10] = generators[..., 0, 3].imag
    p[..., 11] = generators[..., 1, 2].real
    p[..., 12] = generators[..., 1, 2].imag
    p[..., 13] = generators[..., 1, 3].real
    p[..., 14] = generators[..., 1, 3].imag

    return p


def random_params_su22(batch_size, n_gen, device=DEVICE, seeds=None):
    """Random parameters for su(2,2): 15 real params per generator."""
    if seeds is not None:
        all_params = []
        for s in seeds:
            rng = np.random.default_rng(s)
            p = rng.standard_normal((n_gen, N_PARAMS))
            all_params.append(p)
        params_np = np.stack(all_params)
        return torch.tensor(params_np, dtype=FDTYPE, device=device)
    return torch.randn(batch_size, n_gen, N_PARAMS, dtype=FDTYPE, device=device)


# ===================================================================
# Verify su(2,2) constraints: T†η + ηT = 0 and Tr(T) = 0
# ===================================================================
def verify_su22_constraints(generators_np, tol=1e-8):
    """Check that all generators satisfy su(2,2) conditions."""
    eta = np.diag([-1., -1., 1., 1.])
    violations_eta = []
    violations_trace = []
    for T in generators_np:
        residual = T.conj().T @ eta + eta @ T
        violations_eta.append(np.max(np.abs(residual)))
        violations_trace.append(abs(np.trace(T)))

    eta_ok = all(e < tol for e in violations_eta)
    trace_ok = all(e < tol for e in violations_trace)
    return eta_ok, trace_ok, violations_eta, violations_trace


# ===================================================================
# Custom optimizer for indefinite Killing target (SU(2,2))
# ===================================================================
def evolve_su22_indefinite(
    n_gen=15,
    steps=50000,
    lr=0.001,
    alpha=5.0,
    beta=1.0,
    batch_size=5,
    seeds=None,
    verbose=True,
    convergence_threshold=5e-4,
):
    """Optimizer for su(2,2) with indefinite Killing target.

    Target: K eigenvalues sorted → [-1]*7 + [+1]*8
    (7 compact directions + 8 non-compact directions)

    Loss = T_nc + alpha * (T_k_sig + T_k_off) + beta * T_norm
    """
    device = DEVICE

    if seeds is None:
        seeds = list(range(batch_size))

    params = random_params_su22(batch_size, n_gen, device, seeds)

    # Normalize
    with torch.no_grad():
        gens = params_to_su22(params, n_gen)
        norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
        gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
        params = su22_to_params(gens)
    params = params.clone().detach().requires_grad_(True)

    optimizer = torch.optim.Adam([params], lr=lr)

    # Indefinite target: 7 negative (compact), 8 positive (non-compact)
    K_target_sorted = torch.tensor(
        [-1.]*7 + [1.]*8,
        dtype=FDTYPE, device=device
    )

    prev_loss = torch.zeros(batch_size, device=device)
    active = torch.ones(batch_size, dtype=torch.bool, device=device)

    results_per_seed = [{"loss": [], "T_nc": [], "T_k_sig": [], "T_k_off": [],
                          "K_eigs": []} for _ in range(batch_size)]

    for step in range(steps):
        optimizer.zero_grad()
        generators = params_to_su22(params, n_gen)

        # Non-commutativity tension (closure driver)
        T_nc = batched_non_commutativity_tension(generators)

        # Killing metric
        K = batched_killing_metric(generators)  # (B, 15, 15)

        # Eigenvalue-based Killing loss (signature matching)
        K_real = K.real if K.is_complex() else K
        K_eigs = torch.linalg.eigvalsh(K_real)  # (B, 15) ascending
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

            if verbose and step % 5000 == 0:
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
        final = params_to_su22(params, n_gen).cpu().numpy()

    return [(final[b], results_per_seed[b]) for b in range(batch_size)]


# ===================================================================
# Analysis
# ===================================================================
def analyze_result(generators, label, verbose=True):
    """Analyze a set of 4×4 complex generators for su(2,2) properties."""
    n_gen = len(generators)

    # Basic properties
    traces = [np.trace(g) for g in generators]
    is_traceless = all(abs(t) < 1e-6 for t in traces)

    # η(2,2)-preservation check
    eta_ok, trace_ok, eta_errs, trace_errs = verify_su22_constraints(generators)

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

    K_diag = np.diag(K).real

    result = {
        "label": label,
        "n_gen": n_gen,
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
        print(f"    Traceless: {is_traceless}")
        print(f"    eta(2,2)-preserving: {eta_ok} (max err: {max(eta_errs):.2e})")
        print(f"    K eigenvalues: {np.sort(eigs)}")
        print(f"    K signature:   ({n_pos}+, {n_zero}0, {n_neg}-)")
        print(f"    K_offdiag max: {K_offdiag:.2e}")
        print(f"    Closure:       {closes}/{total_pairs} ({closure:.1%})")
        print(f"    f_rms:         {f_rms:.6f}")

        if closure > 0.99 and n_pos == 8 and n_neg == 7:
            print(f"\n    *** su(2,2) conformal algebra identified! ***")
            print(f"    (7 compact + 8 non-compact = conformal algebra of 3+1 spacetime)")
        elif closure > 0.99 and n_pos == 0 and n_neg == 15:
            print(f"\n    *** su(4) compact algebra identified ***")

    return result


def compare_with_canonical(generators, verbose=True):
    """Compare emergent generators with canonical su(2,2) basis.

    Canonical basis: 7 compact + 8 non-compact generators.

    Compact (s(u(2) + u(2))):
      C1-C3: i*sigma_k/2 in top-left 2x2 block (su(2))
      C4-C6: i*sigma_k/2 in bottom-right 2x2 block (su(2))
      C7:    (i/2)*diag(1,1,-1,-1)  (u(1) traceless)

    Non-compact (off-diagonal block):
      N1-N8: T = [[0,B],[B^dag,0]] for B running over 2x2 complex basis
    """
    canonical = []

    # --- Compact: top-left su(2) ---
    C1 = np.zeros((4,4), dtype=complex)
    C1[0,1] = 1j/2; C1[1,0] = 1j/2
    canonical.append(C1)

    C2 = np.zeros((4,4), dtype=complex)
    C2[0,1] = 0.5; C2[1,0] = -0.5
    canonical.append(C2)

    C3 = np.zeros((4,4), dtype=complex)
    C3[0,0] = 1j/2; C3[1,1] = -1j/2
    canonical.append(C3)

    # --- Compact: bottom-right su(2) ---
    C4 = np.zeros((4,4), dtype=complex)
    C4[2,3] = 1j/2; C4[3,2] = 1j/2
    canonical.append(C4)

    C5 = np.zeros((4,4), dtype=complex)
    C5[2,3] = 0.5; C5[3,2] = -0.5
    canonical.append(C5)

    C6 = np.zeros((4,4), dtype=complex)
    C6[2,2] = 1j/2; C6[3,3] = -1j/2
    canonical.append(C6)

    # --- Compact: u(1) traceless ---
    C7 = np.zeros((4,4), dtype=complex)
    C7[0,0] = 1j/2; C7[1,1] = 1j/2; C7[2,2] = -1j/2; C7[3,3] = -1j/2
    canonical.append(C7)

    # --- Non-compact: off-diagonal hermitian ---
    # B_{00} real
    N1 = np.zeros((4,4), dtype=complex)
    N1[0,2] = 1; N1[2,0] = 1
    canonical.append(N1)

    # B_{00} imag
    N2 = np.zeros((4,4), dtype=complex)
    N2[0,2] = 1j; N2[2,0] = -1j
    canonical.append(N2)

    # B_{01} real
    N3 = np.zeros((4,4), dtype=complex)
    N3[0,3] = 1; N3[3,0] = 1
    canonical.append(N3)

    # B_{01} imag
    N4 = np.zeros((4,4), dtype=complex)
    N4[0,3] = 1j; N4[3,0] = -1j
    canonical.append(N4)

    # B_{10} real
    N5 = np.zeros((4,4), dtype=complex)
    N5[1,2] = 1; N5[2,1] = 1
    canonical.append(N5)

    # B_{10} imag
    N6 = np.zeros((4,4), dtype=complex)
    N6[1,2] = 1j; N6[2,1] = -1j
    canonical.append(N6)

    # B_{11} real
    N7 = np.zeros((4,4), dtype=complex)
    N7[1,3] = 1; N7[3,1] = 1
    canonical.append(N7)

    # B_{11} imag
    N8 = np.zeros((4,4), dtype=complex)
    N8[1,3] = 1j; N8[3,1] = -1j
    canonical.append(N8)

    if verbose:
        print("\n    Canonical su(2,2) comparison (15 generators):")

    # Build real representation for lstsq span check
    can_complex = np.array([g.flatten() for g in canonical])  # (15, 16) complex
    can_flat = np.concatenate([can_complex.real, can_complex.imag], axis=-1).T  # (32, 15) real

    span_ok = True
    for i, g in enumerate(generators):
        g_flat = np.concatenate([g.flatten().real, g.flatten().imag])
        coeffs, _, _, _ = np.linalg.lstsq(can_flat, g_flat, rcond=None)
        recon = can_flat @ coeffs
        residual = np.linalg.norm(recon - g_flat)
        in_span = residual < 1e-3
        span_ok = span_ok and in_span
        if verbose:
            if i < 5 or not in_span:
                print(f"      Gen {i}: residual={residual:.2e} "
                      f"{'OK' if in_span else 'NOT IN SPAN'}")

    if verbose:
        if len(generators) > 5:
            print(f"      (showing first 5 of {len(generators)})")
        print(f"    All emergent gens in canonical su(2,2) span: {span_ok}")

    return span_ok


def analyze_compact_subalgebra(generators, verbose=True):
    """Identify and analyze the compact subalgebra from Killing eigenvectors.

    For su(2,2), the 7 compact directions (K<0) should form the
    s(u(2) + u(2)) subalgebra (dim 7).
    """
    n_gen = len(generators)
    K = killing_metric_np(generators).real

    evals, evecs = np.linalg.eigh(K)

    compact_idx = [i for i in range(n_gen) if evals[i] < -1e-6]
    noncompact_idx = [i for i in range(n_gen) if evals[i] > 1e-6]

    n_compact = len(compact_idx)
    n_noncompact = len(noncompact_idx)

    if verbose:
        print(f"\n    Compact/non-compact decomposition:")
        print(f"      Compact directions:     {n_compact} (expected: 7)")
        print(f"      Non-compact directions: {n_noncompact} (expected: 8)")

    if n_compact == 7:
        # Build compact generators as linear combinations via K eigenvectors
        compact_gens = []
        for k in compact_idx:
            coeffs = evecs[:, k]
            gen = sum(c * generators[i] for i, c in enumerate(coeffs))
            compact_gens.append(gen)

        # Check closure of compact subalgebra
        closes = 0
        total_p = n_compact * (n_compact - 1) // 2
        for i in range(n_compact):
            for j in range(i + 1, n_compact):
                C = commutator_np(compact_gens[i], compact_gens[j])
                if np.linalg.norm(C) < 1e-10:
                    closes += 1
                elif check_in_span_np(C, np.array(compact_gens), tol=1e-3):
                    closes += 1
        compact_closure = closes / total_p if total_p > 0 else 1.0

        if verbose:
            print(f"      Compact subalgebra closure: {closes}/{total_p} ({compact_closure:.1%})")
            if compact_closure > 0.9:
                print(f"      ==> s(u(2)+u(2)) compact subalgebra confirmed")

    return n_compact, n_noncompact


# ===================================================================
# Main
# ===================================================================
def main():
    print("=" * 78)
    print("  S14: SU(2,2) CONFORMAL ALGEBRA EMERGENCE TEST")
    print("=" * 78)
    print(f"  Device: {get_device_info()}")
    print(f"  su(2,2) = so(4,2): 15 generators in dim=4 (complex)")
    print(f"  K signature: (8+, 0, 7-) -- conformal algebra of 3+1 spacetime")
    print(f"  su(4) compact: K signature: (0+, 0, 15-)")
    print()

    n_seeds = 5

    # ---------------------------------------------------------------
    # Experiment A: eta(2,2)-preserving + indefinite target -> su(2,2)
    # ---------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  EXPERIMENT A: eta(2,2)-preserving param + indefinite Killing target")
    print("  Target K_eigs sorted: [-1]*7 + [+1]*8")
    print("  Expected: should converge to su(2,2) conformal algebra")
    print("=" * 78)

    t0 = time.time()
    results_a = evolve_su22_indefinite(
        n_gen=N_GEN, steps=50000, lr=0.001,
        alpha=5.0, beta=1.0,
        batch_size=n_seeds, seeds=list(range(n_seeds)),
        verbose=True, convergence_threshold=5e-4,
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
        analyze_compact_subalgebra(best_a[0])

    # ---------------------------------------------------------------
    # Experiment B: su(4) traceless hermitian + compact target (control)
    # ---------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  EXPERIMENT B: su(4) traceless hermitian + compact target K=-I")
    print("  Standard SU(N) parameterization, 15 generators in dim=4")
    print("  Expected: should converge to su(4) compact")
    print("=" * 78)

    t0 = time.time()
    results_b = evolve_generators_generic(
        n_gen=N_GEN, dim=DIM, steps=30000, lr=0.001,
        alpha=5.0, beta=1.0, killing_target=-1.0,
        batch_size=n_seeds, seeds=list(range(n_seeds)),
        verbose=True, convergence_threshold=1e-4,
    )
    time_b = time.time() - t0

    for s, (gens, _) in enumerate(results_b):
        r = analyze_result(list(gens), f"Exp B seed={s}")
        if r["closure"] > 0.99 and r["signature"][2] == 15:
            print(f"    ==> su(4) compact algebra confirmed")

    # ---------------------------------------------------------------
    # Experiment C: eta(2,2)-preserving + compact target -> frustrated
    # ---------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  EXPERIMENT C: eta(2,2)-preserving param + compact target K=-I")
    print("  Expected: FAIL (indefinite K can't match negative-definite target)")
    print("=" * 78)

    t0 = time.time()
    results_c = evolve_generators_generic(
        n_gen=N_GEN, dim=DIM, steps=25000, lr=0.001,
        alpha=5.0, beta=1.0, killing_target=-1.0,
        batch_size=n_seeds, seeds=list(range(n_seeds)),
        verbose=True, convergence_threshold=1e-4,
        param_to_gen_fn=params_to_su22,
        gen_to_param_fn=su22_to_params,
        random_param_fn=random_params_su22,
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
    print(f"\n  {'Experiment':45s} {'Closure':>10s} {'K signature':>16s} "
          f"{'eta-pres':>9s} {'Tr=0':>6s} {'Time':>8s}")
    print("  " + "-" * 98)

    for label, results, t in [
        ("A: su(2,2) eta-pres + indef target", results_a, time_a),
        ("B: su(4) standard + compact target", results_b, time_b),
        ("C: su(2,2) eta-pres + compact target", results_c, time_c),
    ]:
        for s, (gens, _) in enumerate(results):
            r = analyze_result(list(gens), f"{label} s{s}", verbose=False)
            sig = r["signature"]
            eta_str = "Y" if r["eta_ok"] else "-"
            trace_str = "Y" if r["is_traceless"] else "N"
            print(f"  {label + f' s{s}':45s} {r['closure']:>9.1%} "
                  f"  ({sig[0]:>2d}+,{sig[1]}0,{sig[2]:>2d}-) "
                  f"     {eta_str:>2s}    {trace_str:>2s}  {t/n_seeds:>7.1f}s")

    print("  " + "-" * 98)
    print()
    print("  Key predictions:")
    print("  - Exp A: Closure ~ 100%, K signature (8+,0,7-) -> su(2,2) PASS")
    print("  - Exp B: Closure ~ 100%, K signature (0+,0,15-) -> su(4) PASS")
    print("  - Exp C: Closure < 100%, frustrated (wrong target for su(2,2) param)")
    print()
    print("  Physics significance:")
    print("  - su(2,2) = so(4,2) = conformal algebra of 3+1 Minkowski spacetime")
    print("  - Contains so(3,1) Lorentz algebra -> connects to S8b results")
    print("  - The variational principle with eta(2,2) + indefinite Killing")
    print("    selects the conformal algebra as the unique 15D algebra in dim=4")
    print("  - Killing signature discriminates:")
    print("    (8+,7-) -> conformal/spacetime  vs  (0+,15-) -> internal/gauge")
    print()
    print("  Conformal algebra structure:")
    print("    su(2,2) = so(4,2) generators decompose as:")
    print("    +-- 6 Lorentz generators M_uv (rotations + boosts)")
    print("    +-- 4 translation generators P_u")
    print("    +-- 4 special conformal generators K_u")
    print("    +-- 1 dilation generator D")
    print()
    print("  SS-PEG insight:")
    print("    'Not only does spacetime symmetry emerge (so(3,1)),")
    print("     but the FULL conformal symmetry structure of massless")
    print("     field theories is selected by the variational principle.'")


if __name__ == "__main__":
    main()
