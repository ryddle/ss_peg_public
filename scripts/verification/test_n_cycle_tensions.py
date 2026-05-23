"""
Independent test: n-cycle tensions and higher-order correlations.

Measures how higher-order algebraic structures (ternary, quaternary, ...)
manifest in different generator sets:

Metrics computed:
  T^(2) = Σ_{i<j} ||[T_i, T_j]||²              (binary / 2-cycle)
  T^(3) = Σ_{i<j<k} ||J(T_i,T_j,T_k)||²        (ternary / 3-cycle / Jacobiator)
  Φ^(3) = Σ_{i<j<k} Tr([T_i,T_j] T_k)           (3-cocycle / Cartan form)
  T^(4) = Σ_{i<j<k<l} ||A4(T_i,T_j,T_k,T_l)||²  (quaternary / 4-cycle)

where:
  J(a,b,c)    = [a,[b,c]] + [b,[c,a]] + [c,[a,b]]          (Jacobiator)
  A4(a,b,c,d) = [[a,b],[c,d]] + [[b,c],[a,d]] + [[c,a],[b,d]]  (4-associator)

Theoretical predictions:
  - Exact Lie algebras: T^(3) = 0 (Jacobi identity), T^(4) = 0
  - Random generators: T^(3) >> 0, T^(4) >> 0
  - The 3-cocycle Φ^(3) is a topological invariant (nonzero for any Lie algebra
    with nontrivial 3rd cohomology; zero only for abelian algebras)

Test bank:
  1. Classical constructive: SU(2), SU(3), SO(3)
  2. Exceptional constructive: G₂, F₄, E₆
  3. Random (non-Lie) antisymmetric matrices: control
  4. Trajectory snapshots during blind optimization: SU(2) vs G₂

No existing files are modified. Imports only from shared_utils_gpu.py
and octonion_jordan.py.

Usage:
    python test_n_cycle_tensions.py                # full test suite
    python test_n_cycle_tensions.py --quick         # skip trajectory analysis
    python test_n_cycle_tensions.py --only-static   # only constructive + random
"""
import sys
import os
import time
import argparse
from itertools import combinations

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from shared_utils_gpu import (
    killing_metric_np,
    commutator_np,
    check_in_span_np,
    structure_constants_general_np,
    get_device_info,
    # GPU optimization (for trajectory snapshots)
    evolve_generators_generic,
    params_to_antisymmetric,
    antisymmetric_to_params,
    random_params_so,
    params_to_generators,
    generators_to_params,
    random_params,
    DTYPE,
    FDTYPE,
    DEVICE,
)

try:
    from octonion_jordan import build_f4_generators, build_e6_generators
    HAS_OCTONION = True
except ImportError:
    HAS_OCTONION = False
    print("[WARN] octonion_jordan.py not found — F₄/E₆ tests will be skipped")

import torch


# ============================================================
# n-CYCLE TENSION METRICS (numpy, post-hoc analysis)
# ============================================================

def tension_2(generators):
    """T^(2) = Σ_{i<j} ||[T_i, T_j]||_F²  (binary / non-commutativity)."""
    n = len(generators)
    t2 = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            C = commutator_np(generators[i], generators[j])
            t2 += np.linalg.norm(C) ** 2
    return t2


def jacobiator(a, b, c):
    """J(a,b,c) = [a,[b,c]] + [b,[c,a]] + [c,[a,b]]."""
    return (commutator_np(a, commutator_np(b, c))
            + commutator_np(b, commutator_np(c, a))
            + commutator_np(c, commutator_np(a, b)))


def tension_3(generators):
    """T^(3) = Σ_{i<j<k} ||J(T_i, T_j, T_k)||_F²  (ternary / Jacobiator).

    NOTE: For matrix algebras this is ALWAYS zero because the Jacobi
    identity holds automatically in associative algebras.
    Kept as a consistency check.
    """
    n = len(generators)
    t3 = 0.0
    for i, j, k in combinations(range(n), 3):
        J = jacobiator(generators[i], generators[j], generators[k])
        t3 += np.linalg.norm(J) ** 2
    return t3


def closure_tension(generators):
    """T_closure = Σ_{i<j} ||[T_i,T_j] - proj_span([T_i,T_j])||²_F.

    Measures how much commutator information LEAKS outside the algebra.
    For a closed Lie algebra: T_closure = 0.
    For random matrices: T_closure >> 0.
    This is the TRUE non-trivial "ternary" probe — it asks whether
    binary interactions produce outputs containable within the algebra.
    """
    n = len(generators)
    if n < 2:
        return 0.0, 0.0, []

    # Build basis matrix for projection: each generator flattened as a column
    B = np.array([g.flatten() for g in generators]).T  # (dim², n)

    residuals = []
    total_leak = 0.0
    total_comm_norm = 0.0

    for i in range(n):
        for j in range(i + 1, n):
            C = commutator_np(generators[i], generators[j])
            c_flat = C.flatten()
            norm_c = np.linalg.norm(c_flat)
            total_comm_norm += norm_c ** 2

            if norm_c < 1e-14:
                residuals.append(0.0)
                continue

            # Project onto span
            coeffs, _, _, _ = np.linalg.lstsq(B, c_flat, rcond=None)
            proj = B @ coeffs
            residual = np.linalg.norm(c_flat - proj)
            residuals.append(residual)
            total_leak += residual ** 2

    # Fractional leak: what fraction of commutator energy escapes the algebra
    frac_leak = total_leak / total_comm_norm if total_comm_norm > 1e-20 else 0.0

    return total_leak, frac_leak, residuals


def triple_structure_tensor(generators):
    """Compute the completely antisymmetric 3-tensor d_{ijk} = Tr({T_i,T_j}T_k).

    For SU(N), d_{ijk} is the symmetric structure constant (anomaly coefficient).
    For SO(N), the anticommutator {T_i,T_j} = T_iT_j + T_jT_i may have different
    structure. The norm ||d||² and its ratio to ||f||² (from the commutator)
    characterizes how much "ternary symmetric content" the algebra has vs
    "ternary antisymmetric content".
    """
    n = len(generators)
    d_norm_sq = 0.0
    f_norm_sq = 0.0

    for i in range(n):
        for j in range(i + 1, n):
            # Antisymmetric (commutator) — structure constants f
            C = commutator_np(generators[i], generators[j])
            for k in range(n):
                fk = np.trace(C @ generators[k])
                f_norm_sq += abs(fk) ** 2

            # Symmetric (anticommutator) — d-symbol
            A = generators[i] @ generators[j] + generators[j] @ generators[i]
            for k in range(n):
                dk = np.trace(A @ generators[k])
                d_norm_sq += abs(dk) ** 2

    return d_norm_sq, f_norm_sq


def cocycle_3(generators):
    """Φ^(3) = Σ_{i<j<k} Tr([T_i,T_j] T_k)  (3-cocycle / Cartan form).

    This is the totally antisymmetric 3-form on the algebra.
    For a simple Lie algebra it is proportional to f_{ijk} (structure constants).
    """
    n = len(generators)
    phi = 0.0
    for i, j, k in combinations(range(n), 3):
        C = commutator_np(generators[i], generators[j])
        phi += np.trace(C @ generators[k]).real
    return phi


def associator_4(a, b, c, d):
    """A4(a,b,c,d) = [[a,b],[c,d]] + [[b,c],[a,d]] + [[c,a],[b,d]].

    Cyclic 4-associator: measures failure of nested pairwise commutators
    to form a consistent 4-body interaction. For any Lie algebra this
    is determined by the Jacobi identity and equals [[a,b],[c,d]] + cyc,
    which generally is NOT zero even for Lie algebras (unlike the Jacobiator).
    """
    ab = commutator_np(a, b)
    bc = commutator_np(b, c)
    ca = commutator_np(c, a)
    return (commutator_np(ab, commutator_np(c, d))
            + commutator_np(bc, commutator_np(a, d))
            + commutator_np(ca, commutator_np(b, d)))


def tension_4(generators, max_tuples=5000):
    """T^(4) = Σ_{i<j<k<l} ||A4(T_i,T_j,T_k,T_l)||_F².

    For large algebras the number of 4-tuples grows as C(n,4), so we
    optionally sample a random subset for feasibility.
    """
    n = len(generators)
    all_tuples = list(combinations(range(n), 4))
    if len(all_tuples) > max_tuples:
        rng = np.random.default_rng(42)
        indices = rng.choice(len(all_tuples), size=max_tuples, replace=False)
        tuples = [all_tuples[idx] for idx in indices]
        sampled = True
    else:
        tuples = all_tuples
        sampled = False

    t4 = 0.0
    for i, j, k, l in tuples:
        A = associator_4(generators[i], generators[j],
                         generators[k], generators[l])
        t4 += np.linalg.norm(A) ** 2

    # Scale to estimate total if sampled
    if sampled:
        t4 *= len(all_tuples) / len(tuples)

    return t4, sampled


def normalized_tensions(generators, compute_t4=True):
    """Compute all tensions, normalized per tuple for comparability.

    Returns dict with raw and per-tuple normalized values.
    """
    n = len(generators)
    n2 = n * (n - 1) // 2
    n3 = n * (n - 1) * (n - 2) // 6
    n4 = n * (n - 1) * (n - 2) * (n - 3) // 24

    t2 = tension_2(generators)
    t3 = tension_3(generators)
    phi3 = cocycle_3(generators)

    # Closure tension: the real discriminator
    cl_total, cl_frac, cl_residuals = closure_tension(generators)

    # Triple structure: d-symbol vs f-symbol
    d_norm_sq, f_norm_sq = triple_structure_tensor(generators)
    df_ratio = d_norm_sq / f_norm_sq if f_norm_sq > 1e-20 else 0.0

    result = {
        "n_gen": n,
        "dim": generators[0].shape[0],
        "T2_raw": t2,
        "T3_raw": t3,
        "Phi3_raw": phi3,
        "T2_per_pair": t2 / n2 if n2 > 0 else 0,
        "T3_per_triple": t3 / n3 if n3 > 0 else 0,
        "Phi3_per_triple": phi3 / n3 if n3 > 0 else 0,
        "n_pairs": n2,
        "n_triples": n3,
        "n_quads": n4,
        # Closure tension
        "closure_leak": cl_total,
        "closure_frac": cl_frac,
        "closure_max_residual": max(cl_residuals) if cl_residuals else 0.0,
        # d/f structure
        "d_norm_sq": d_norm_sq,
        "f_norm_sq": f_norm_sq,
        "d_over_f": df_ratio,
    }

    # T4 only if feasible and requested
    if compute_t4 and n4 > 0 and n4 < 100000:
        t4, sampled = tension_4(generators)
        result["T4_raw"] = t4
        result["T4_per_quad"] = t4 / n4 if n4 > 0 else 0
        result["T4_sampled"] = sampled
    elif compute_t4 and n4 > 0:
        t4, sampled = tension_4(generators, max_tuples=5000)
        result["T4_raw"] = t4
        result["T4_per_quad"] = t4 / n4 if n4 > 0 else 0
        result["T4_sampled"] = sampled
    else:
        result["T4_raw"] = 0.0
        result["T4_per_quad"] = 0.0
        result["T4_sampled"] = False

    return result


# ============================================================
# CONSTRUCTIVE GENERATOR BUILDERS
# ============================================================

def build_su2_generators():
    """SU(2): 3 generators in 2×2 (Pauli/2)."""
    sigma = [
        np.array([[0, 1], [1, 0]], dtype=complex),
        np.array([[0, -1j], [1j, 0]], dtype=complex),
        np.array([[1, 0], [0, -1]], dtype=complex),
    ]
    return [s / 2.0 for s in sigma]


def build_su3_generators():
    """SU(3): 8 Gell-Mann matrices / 2."""
    lam = [np.zeros((3, 3), dtype=complex) for _ in range(8)]
    lam[0][0, 1] = lam[0][1, 0] = 1
    lam[1][0, 1] = -1j; lam[1][1, 0] = 1j
    lam[2][0, 0] = 1; lam[2][1, 1] = -1
    lam[3][0, 2] = lam[3][2, 0] = 1
    lam[4][0, 2] = -1j; lam[4][2, 0] = 1j
    lam[5][1, 2] = lam[5][2, 1] = 1
    lam[6][1, 2] = -1j; lam[6][2, 1] = 1j
    lam[7][0, 0] = 1; lam[7][1, 1] = 1; lam[7][2, 2] = -2
    lam[7] /= np.sqrt(3)
    return [l / 2.0 for l in lam]


def build_so3_generators():
    """SO(3): 3 generators in 3×3 real antisymmetric."""
    L1 = np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]], dtype=float)
    L2 = np.array([[0, 0, 1], [0, 0, 0], [-1, 0, 0]], dtype=float)
    L3 = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], dtype=float)
    return [L1, L2, L3]


def build_g2_generators():
    """G₂: 14 generators in 7×7 real antisymmetric (from octonion 3-form)."""
    fano_triples = [
        (0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6),
        (4, 5, 0), (5, 6, 1), (6, 0, 2)
    ]
    phi = np.zeros((7, 7, 7))
    for (a, b, c) in fano_triples:
        phi[a, b, c] = 1; phi[b, c, a] = 1; phi[c, a, b] = 1
        phi[a, c, b] = -1; phi[c, b, a] = -1; phi[b, a, c] = -1

    basis = []
    for i in range(7):
        for j in range(i + 1, 7):
            E = np.zeros((7, 7))
            E[i, j] = 1; E[j, i] = -1
            basis.append(E)

    ijk_list = [(i, j, k) for i in range(7)
                for j in range(i + 1, 7) for k in range(j + 1, 7)]

    M = np.zeros((len(ijk_list), len(basis)))
    for row, (i, j, k) in enumerate(ijk_list):
        for col, E in enumerate(basis):
            val = 0.0
            for a in range(7):
                val += E[i, a] * phi[a, j, k]
                val += E[j, a] * phi[i, a, k]
                val += E[k, a] * phi[i, j, a]
            M[row, col] = val

    _, S, Vt = np.linalg.svd(M)
    rank = np.sum(S > 1e-10)
    g2_coefficients = Vt[rank:]

    generators = []
    for coeffs in g2_coefficients:
        G = sum(c * E for c, E in zip(coeffs, basis))
        generators.append(G)

    # Gram-Schmidt orthonormalize (Frobenius)
    ortho = []
    for g in generators:
        v = g.copy()
        for u in ortho:
            v -= (np.sum(v * u) / np.sum(u * u)) * u
        norm = np.sqrt(np.sum(v * v))
        if norm > 1e-12:
            ortho.append(v / norm)
    return ortho


def build_random_antisymmetric(n_gen, dim, seed=0):
    """Random real antisymmetric matrices (NOT a Lie algebra — control)."""
    rng = np.random.default_rng(seed)
    generators = []
    for _ in range(n_gen):
        A = rng.standard_normal((dim, dim))
        A = (A - A.T) / 2.0
        A /= np.linalg.norm(A)
        generators.append(A)
    return generators


def build_random_hermitian(n_gen, dim, seed=0):
    """Random hermitian traceless matrices (NOT a Lie algebra — control)."""
    rng = np.random.default_rng(seed)
    generators = []
    for _ in range(n_gen):
        A = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
        H = (A + A.conj().T) / 2.0
        H -= np.trace(H) / dim * np.eye(dim)
        H /= np.linalg.norm(H)
        generators.append(H)
    return generators


# ============================================================
# G₂ 3-FORM COCYCLE
# ============================================================

def g2_3form_cocycle(generators):
    """Compute the G₂ octonionic 3-form invariant.

    For G₂ generators in the 7-dim rep, the 3-form φ_{ijk} from the
    Fano plane provides an additional invariant:

      Φ_G₂ = Σ_{(i,j,k)∈Fano} Tr([T_i, T_j] @ some_projection)

    More precisely, we compute how much the structure constants align
    with the octonionic 3-form, measuring the "ternary content":

      align = Σ_{a<b} Σ_c f_{abc} · φ_{c,idx_a,idx_b}

    This is only meaningful for 7×7 generators.
    """
    if generators[0].shape[0] != 7:
        return None

    fano_triples = [
        (0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6),
        (4, 5, 0), (5, 6, 1), (6, 0, 2)
    ]
    phi = np.zeros((7, 7, 7))
    for (a, b, c) in fano_triples:
        phi[a, b, c] = 1; phi[b, c, a] = 1; phi[c, a, b] = 1
        phi[a, c, b] = -1; phi[c, b, a] = -1; phi[b, a, c] = -1

    # Structure constants
    f = structure_constants_general_np(generators)
    n = len(generators)

    # Alignment: how much f_{abc} correlates with φ (only if we can
    # map generator indices to octonion indices — this is nontrivial
    # in general, so we use a representation-level approach instead)

    # Representation-level: for 7×7 generators, compute
    # Σ_{ijk in Fano triples} Σ_a [T_a]_{ij} [T_a]_{jk} [T_a]_{ki}
    # This measures how much each generator "sees" the ternary structure

    ternary_signal = 0.0
    for (i, j, k) in fano_triples:
        for a in range(n):
            T = generators[a]
            ternary_signal += (T[i, j] * T[j, k] * T[k, i]).real
            ternary_signal -= (T[i, k] * T[k, j] * T[j, i]).real

    return ternary_signal


# ============================================================
# TRAJECTORY ANALYSIS (GPU → numpy snapshots)
# ============================================================

def trajectory_n_cycle_analysis(
    algebra_name,
    n_gen,
    dim,
    steps=5000,
    lr=0.005,
    alpha=5.0,
    killing_target=-1.0,
    snapshot_every=500,
    param_type="hermitian",
    seed=42,
):
    """Run blind optimization and measure n-cycle tensions at intervals.

    Returns list of dicts, one per snapshot, with step number + all tensions.
    """
    print(f"\n{'='*60}")
    print(f"  Trajectory analysis: {algebra_name}")
    print(f"  n_gen={n_gen}, dim={dim}, steps={steps}, snapshot_every={snapshot_every}")
    print(f"{'='*60}")

    device = DEVICE

    if param_type == "antisymmetric":
        param_fn = params_to_antisymmetric
        gen_to_param_fn = antisymmetric_to_params
        rand_fn = lambda bs, ng, device, seeds: random_params_so(bs, ng, dim, device, seeds)
        ptg_kwargs = {"dim": dim}
    else:  # hermitian
        param_fn = params_to_generators
        gen_to_param_fn = generators_to_params
        rand_fn = lambda bs, ng, device, seeds: random_params(bs, ng, dim, device, seeds)
        ptg_kwargs = {"dim": dim}

    params = rand_fn(1, n_gen, device=device, seeds=[seed])
    with torch.no_grad():
        gens = param_fn(params, n_gen, **ptg_kwargs)
        norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt().clamp(min=1e-10)
        gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
        params = gen_to_param_fn(gens)
    params = params.clone().detach().requires_grad_(True)

    optimizer = torch.optim.Adam([params], lr=lr)
    K_target_mat = killing_target * torch.eye(n_gen, dtype=FDTYPE, device=device)

    snapshots = []

    for step in range(steps + 1):
        if step > 0:
            optimizer.zero_grad()
            generators = param_fn(params, n_gen, **ptg_kwargs)
            from shared_utils_gpu import (
                batched_non_commutativity_tension,
                batched_killing_metric,
                batched_norm_tension,
            )
            T_nc = batched_non_commutativity_tension(generators)
            K = batched_killing_metric(generators)
            T_k = ((K - K_target_mat.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
            T_norm = batched_norm_tension(generators)
            T_total = T_nc + alpha * T_k + T_norm
            T_total.sum().backward()
            optimizer.step()

        # Snapshot
        if step % snapshot_every == 0 or step == steps:
            with torch.no_grad():
                generators_t = param_fn(params, n_gen, **ptg_kwargs)
                gens_np = generators_t[0].cpu().numpy()

            # Compute n-cycle tensions on numpy
            metrics = normalized_tensions(list(gens_np), compute_t4=False)
            metrics["step"] = step

            # Killing diagonal
            K_np = killing_metric_np(list(gens_np))
            metrics["K_diag_mean"] = np.diag(K_np).mean().real

            snapshots.append(metrics)

            print(f"  Step {step:5d}: T2/pair={metrics['T2_per_pair']:.6f}  "
                  f"Cl.frac={metrics['closure_frac']:.2e}  "
                  f"Phi3/tri={metrics['Phi3_per_triple']:.4f}  "
                  f"d/f={metrics['d_over_f']:.4f}  "
                  f"close={metrics.get('closure_frac_check', metrics.get('closure_frac', 0)):.2e}  "
                  f"K_diag={metrics['K_diag_mean']:.4f}")

    return snapshots


# ============================================================
# DISPLAY
# ============================================================

def print_tensions(name, metrics):
    """Pretty-print tension metrics for one algebra."""
    print(f"\n  {name}:")
    print(f"    n_gen={metrics['n_gen']}, dim={metrics['dim']}")
    print(f"    T^(2) total={metrics['T2_raw']:.6f}  per pair={metrics['T2_per_pair']:.6f}")
    print(f"    T^(3) total={metrics['T3_raw']:.2e}  per triple={metrics['T3_per_triple']:.2e}  "
          f"(Jacobi check)")
    print(f"    Phi^(3) total={metrics['Phi3_raw']:.6f}  per triple={metrics['Phi3_per_triple']:.6f}")
    print(f"    Closure leak={metrics['closure_leak']:.6f}  "
          f"frac={metrics['closure_frac']:.2e}  "
          f"max_residual={metrics['closure_max_residual']:.2e}")
    print(f"    d/f ratio={metrics['d_over_f']:.6f}  "
          f"(||d||^2={metrics['d_norm_sq']:.4f}, ||f||^2={metrics['f_norm_sq']:.4f})")
    if "T4_raw" in metrics:
        s = " (sampled)" if metrics.get("T4_sampled") else ""
        print(f"    T^(4) total={metrics['T4_raw']:.2e}  per quad={metrics['T4_per_quad']:.2e}{s}")


def print_comparison_table(results):
    """Print comparison table of all algebras."""
    print(f"\n{'='*110}")
    print(f"  {'Algebra':<20} {'n':>3} {'d':>3} {'T2/pair':>10} "
          f"{'Cl.frac':>10} {'Phi3/tri':>10} {'d/f':>10} {'T4/quad':>10}")
    print(f"  {'-'*20} {'-'*3} {'-'*3} {'-'*10} "
          f"{'-'*10} {'-'*10} {'-'*10} {'-'*10}")
    for name, m in results:
        t4 = f"{m['T4_per_quad']:.2e}" if "T4_per_quad" in m else "N/A"
        print(f"  {name:<20} {m['n_gen']:>3} {m['dim']:>3} "
              f"{m['T2_per_pair']:>10.6f} {m['closure_frac']:>10.2e} "
              f"{m['Phi3_per_triple']:>10.6f} {m['d_over_f']:>10.6f} {t4:>10}")
    print(f"{'='*110}")


# ============================================================
# MAIN TEST
# ============================================================

def run_static_tests():
    """Compute n-cycle tensions for constructive + random generator sets."""
    print("\n" + "=" * 70)
    print("  PART 1: STATIC n-CYCLE TENSION ANALYSIS")
    print("  (Constructive algebras vs random controls)")
    print("=" * 70)

    results = []

    # --- Classical ---
    print("\n--- Classical Lie algebras ---")

    gens = build_su2_generators()
    m = normalized_tensions(gens)
    print_tensions("SU(2) [3 gens, 2×2]", m)
    results.append(("SU(2)", m))

    gens = build_su3_generators()
    m = normalized_tensions(gens)
    print_tensions("SU(3) [8 gens, 3×3]", m)
    results.append(("SU(3)", m))

    gens = build_so3_generators()
    m = normalized_tensions(gens)
    print_tensions("SO(3) [3 gens, 3×3]", m)
    results.append(("SO(3)", m))

    # --- Exceptional ---
    print("\n--- Exceptional Lie algebras ---")

    gens = build_g2_generators()
    m = normalized_tensions(gens)
    m["g2_3form"] = g2_3form_cocycle(gens)
    print_tensions("G₂ [14 gens, 7×7]", m)
    if m["g2_3form"] is not None:
        print(f"    G₂ octonionic 3-form signal: {m['g2_3form']:.6f}")
    results.append(("G₂", m))

    if HAS_OCTONION:
        gens_f4 = build_f4_generators()
        m = normalized_tensions(gens_f4)
        print_tensions("F₄ [52 gens, 27×27]", m)
        results.append(("F₄", m))

        gens_e6 = build_e6_generators()
        m = normalized_tensions(gens_e6)
        print_tensions("E₆ [78 gens, 27×27]", m)
        results.append(("E₆", m))

    # --- Random controls ---
    print("\n--- Random controls (NOT Lie algebras) ---")

    # Random antisymmetric in same dims as G₂
    gens = build_random_antisymmetric(14, 7, seed=0)
    m = normalized_tensions(gens)
    print_tensions("Random SO(7)-type [14, 7×7]", m)
    results.append(("Rand-SO(7)", m))

    # Random hermitian in same dims as SU(3)
    gens = build_random_hermitian(8, 3, seed=0)
    m = normalized_tensions(gens)
    print_tensions("Random SU(3)-type [8, 3×3]", m)
    results.append(("Rand-SU(3)", m))

    # Random antisymmetric in 27×27 (same as F₄/E₆ rep)
    gens = build_random_antisymmetric(14, 27, seed=0)
    m = normalized_tensions(gens)
    print_tensions("Random 27×27 [14 gens]", m)
    results.append(("Rand-27", m))

    # --- Comparison table ---
    print_comparison_table(results)

    # --- Key predictions check ---
    print("\n--- Prediction verification ---")

    # All Lie algebras should have T^(3) ≈ 0 (Jacobi, auto in associative algebras)
    print("\n  [1] Jacobi identity (T³ = 0 for ALL matrix sets):")
    for name, m in results:
        t3 = m["T3_per_triple"]
        status = "OK" if abs(t3) < 1e-10 else "UNEXPECTED"
        print(f"    {name}: T3/triple = {t3:.2e}  {status}")

    # Closure: Lie algebras should have closure_frac ≈ 0, random should have >> 0
    print("\n  [2] Closure (commutators stay in span — TRUE discriminator):")
    for name, m in results:
        cf = m["closure_frac"]
        if name.startswith("Rand"):
            status = "PASS (leaks)" if cf > 1e-6 else "FAIL (unexpectedly closed)"
        else:
            status = "PASS (closed)" if cf < 1e-6 else "FAIL (leaks!)"
        print(f"    {name}: closure_frac = {cf:.2e}  {status}")

    # d/f ratio: characterizes the algebra type
    print("\n  [3] d/f ratio (symmetric vs antisymmetric ternary content):")
    for name, m in results:
        print(f"    {name}: d/f = {m['d_over_f']:.6f}")

    # Cocycle Φ^(3): should be nonzero for non-abelian Lie algebras
    print("\n  [4] 3-cocycle Phi^(3) (Cartan form):")
    for name, m in results:
        if name.startswith("Rand"):
            continue
        phi = m["Phi3_per_triple"]
        status = "nonzero" if abs(phi) > 1e-10 else "(abelian or traceless)"
        print(f"    {name}: Phi3/triple = {phi:.6f}  {status}")

    return results


def run_trajectory_tests():
    """Trajectory analysis: watch n-cycle tensions during optimization."""
    print("\n" + "=" * 70)
    print("  PART 2: TRAJECTORY n-CYCLE TENSION ANALYSIS")
    print("  (How do n-cycle tensions evolve during blind optimization?)")
    print("=" * 70)

    all_trajectories = {}

    # SU(2) blind: should converge (classical, easy)
    snaps = trajectory_n_cycle_analysis(
        "SU(2) blind [3 gens, 2×2]",
        n_gen=3, dim=2, steps=3000, lr=0.005,
        alpha=5.0, killing_target=-1.0,
        snapshot_every=300, param_type="hermitian", seed=42,
    )
    all_trajectories["SU(2)_blind"] = snaps

    # SO(3) blind: should converge
    snaps = trajectory_n_cycle_analysis(
        "SO(3) blind [3 gens, 3×3]",
        n_gen=3, dim=3, steps=3000, lr=0.005,
        alpha=5.0, killing_target=-1.0,
        snapshot_every=300, param_type="antisymmetric", seed=42,
    )
    all_trajectories["SO(3)_blind"] = snaps

    # G₂ blind in SO(7): will NOT converge but trajectory is informative
    snaps = trajectory_n_cycle_analysis(
        "G₂ blind [14 gens, 7×7 SO(7)]",
        n_gen=14, dim=7, steps=5000, lr=0.003,
        alpha=5.0, killing_target=-1.0,
        snapshot_every=500, param_type="antisymmetric", seed=42,
    )
    all_trajectories["G2_blind"] = snaps

    # SU(3) blind: should converge (classical, medium)
    snaps = trajectory_n_cycle_analysis(
        "SU(3) blind [8 gens, 3×3]",
        n_gen=8, dim=3, steps=5000, lr=0.003,
        alpha=5.0, killing_target=-1.0,
        snapshot_every=500, param_type="hermitian", seed=42,
    )
    all_trajectories["SU(3)_blind"] = snaps

    # --- Compare final states ---
    print(f"\n{'='*90}")
    print(f"  TRAJECTORY FINAL STATE COMPARISON")
    print(f"{'='*90}")
    print(f"  {'Algebra':<25} {'Step':>5} {'T2/pair':>10} {'Cl.frac':>10} "
          f"{'Phi3/tri':>10} {'d/f':>10}")
    print(f"  {'-'*25} {'-'*5} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")

    for name, snaps in all_trajectories.items():
        final = snaps[-1]
        print(f"  {name:<25} {final['step']:>5} "
              f"{final['T2_per_pair']:>10.4f} {final['closure_frac']:>10.2e} "
              f"{final['Phi3_per_triple']:>10.4f} {final['d_over_f']:>10.4f}")

    # --- Evolution summary ---
    print(f"\n  TRAJECTORY EVOLUTION: Closure fraction decay profile")
    print(f"  {'-'*70}")
    for name, snaps in all_trajectories.items():
        cl_vals = [s["closure_frac"] for s in snaps]
        if len(cl_vals) >= 2 and cl_vals[0] > 1e-20:
            ratio = cl_vals[-1] / cl_vals[0] if cl_vals[0] > 1e-20 else 0
            print(f"  {name:<25}: Cl.frac initial={cl_vals[0]:.2e} -> "
                  f"final={cl_vals[-1]:.2e}  (ratio={ratio:.2e})")
        else:
            print(f"  {name:<25}: Cl.frac initial={cl_vals[0]:.2e} -> "
                  f"final={cl_vals[-1]:.2e}")

    return all_trajectories


# ============================================================
# ENTRY POINT
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="n-cycle tension analysis")
    parser.add_argument("--quick", action="store_true",
                        help="Skip trajectory analysis")
    parser.add_argument("--only-static", action="store_true",
                        help="Only constructive + random (no trajectories)")
    parser.add_argument("--only-trajectory", action="store_true",
                        help="Only trajectory analysis")
    args = parser.parse_args()

    print("=" * 70)
    print("  n-CYCLE TENSION TEST")
    print("  Independent measurement of higher-order algebraic correlations")
    print("=" * 70)
    print(f"  Device: {get_device_info()}")
    print(f"  Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    t0 = time.time()

    if not args.only_trajectory:
        run_static_tests()

    if not args.only_static and not args.quick:
        run_trajectory_tests()

    elapsed = time.time() - t0
    print(f"\n  Total elapsed: {elapsed:.1f}s")
    print("  Done.")


if __name__ == "__main__":
    main()
