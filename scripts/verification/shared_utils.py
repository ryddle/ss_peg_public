"""
Shared utilities for the systematic verification plan.
All test scripts import from here to avoid code duplication.
"""
import numpy as np


# =========================
# BASIC ALGEBRA
# =========================

def pauli_matrices():
    sigma_x = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    return [sigma_x, sigma_y, sigma_z]


def gell_mann_matrices():
    """The 8 Gell-Mann matrices (generators of su(3))."""
    l1 = np.array([[0,1,0],[1,0,0],[0,0,0]], dtype=np.complex128)
    l2 = np.array([[0,-1j,0],[1j,0,0],[0,0,0]], dtype=np.complex128)
    l3 = np.array([[1,0,0],[0,-1,0],[0,0,0]], dtype=np.complex128)
    l4 = np.array([[0,0,1],[0,0,0],[1,0,0]], dtype=np.complex128)
    l5 = np.array([[0,0,-1j],[0,0,0],[1j,0,0]], dtype=np.complex128)
    l6 = np.array([[0,0,0],[0,0,1],[0,1,0]], dtype=np.complex128)
    l7 = np.array([[0,0,0],[0,0,-1j],[0,1j,0]], dtype=np.complex128)
    l8 = np.array([[1,0,0],[0,1,0],[0,0,-2]], dtype=np.complex128) / np.sqrt(3)
    return [l1, l2, l3, l4, l5, l6, l7, l8]


def commutator(A, B):
    return A @ B - B @ A


def random_hermitian(dim=2):
    A = np.random.randn(dim, dim) + 1j * np.random.randn(dim, dim)
    return (A + A.conj().T) / 2


def random_traceless_hermitian(dim=2):
    H = random_hermitian(dim)
    H = H - np.trace(H) * np.eye(dim) / dim
    return H


def random_generators(n=3, dim=2):
    return [random_traceless_hermitian(dim) for _ in range(n)]


# =========================
# TENSIONS
# =========================

def non_commutativity_tension(generators):
    T = 0.0
    n = len(generators)
    for i in range(n):
        for j in range(i + 1, n):
            C = commutator(generators[i], generators[j])
            T += np.linalg.norm(C) ** 2
    return T


def adjoint(T, basis):
    return [commutator(T, B) for B in basis]


def killing_metric(generators):
    n = len(generators)
    K = np.zeros((n, n), dtype=np.complex128)
    for i in range(n):
        ad_i = adjoint(generators[i], generators)
        for j in range(n):
            ad_j = adjoint(generators[j], generators)
            total = 0
            for A, B in zip(ad_i, ad_j):
                total += np.trace(A @ B)
            K[i, j] = total
    return np.real(K)


def norm_tension(generators, target_norm=1.0):
    return sum((np.linalg.norm(G) - target_norm) ** 2 for G in generators)


def non_comm_constraint(generators, target=1.0):
    T = 0.0
    n = len(generators)
    for i in range(n):
        for j in range(i + 1, n):
            C = commutator(generators[i], generators[j])
            val = np.linalg.norm(C)
            T += (val - target) ** 2
    return T


# =========================
# SPAN CHECK
# =========================

def check_in_span(vec, basis, tol=1e-6):
    B = np.array([b.flatten() for b in basis]).T
    v = vec.flatten()
    coeffs, _, _, _ = np.linalg.lstsq(B, v, rcond=None)
    return np.linalg.norm(B @ coeffs - v) < tol


# =========================
# STRUCTURE CONSTANTS
# =========================

def structure_constants(generators):
    n = len(generators)
    f = np.zeros((n, n, n), dtype=np.complex128)
    for i in range(n):
        for j in range(n):
            C = commutator(generators[i], generators[j])
            for k in range(n):
                f[i, j, k] = np.trace(C @ generators[k])
    return f


def structure_constants_real(generators):
    """Extract real structure constants f_abc from [T_a, T_b] = i f_abc T_c.

    For hermitian generators, [T_a, T_b] is anti-hermitian, so
    Tr([T_a, T_b] · T_c) is purely imaginary. The real structure
    constants are: f_abc = -i * Tr([T_a, T_b] · T_c) / Tr(T_c^2).
    """
    n = len(generators)

    norms_sq = np.array([np.trace(generators[k] @ generators[k]).real for k in range(n)])

    f = np.zeros((n, n, n))
    for i in range(n):
        for j in range(n):
            C = commutator(generators[i], generators[j])
            for k in range(n):
                if norms_sq[k] > 1e-12:
                    # -i * Tr([T_a, T_b] T_c) extracts the real part
                    f[i, j, k] = (-1j * np.trace(C @ generators[k]) / norms_sq[k]).real
    return f


# =========================
# ORTHONORMALIZATION
# =========================

def orthonormalize_killing(generators):
    """Gram-Schmidt orthonormalization using Killing metric as inner product.

    Returns a new list of generators that are orthonormal w.r.t.
    K(T_a, T_b) = Tr(ad_a . ad_b), i.e. K_ab = -delta_ab after
    normalization. Falls back to Frobenius Gram-Schmidt if Killing
    is degenerate.
    """
    n = len(generators)
    ortho = []
    for i in range(n):
        v = generators[i].copy()
        for u in ortho:
            # inner product via Killing: k(u,v) = sum_c Tr([u,G_c] @ [v,G_c])
            ip = sum(np.trace(commutator(u, g) @ commutator(v, g)).real
                     for g in generators)
            norm_u = sum(np.trace(commutator(u, g) @ commutator(u, g)).real
                         for g in generators)
            if abs(norm_u) > 1e-12:
                v = v - (ip / norm_u) * u
        # Normalize
        norm_v = sum(np.trace(commutator(v, g) @ commutator(v, g)).real
                     for g in generators)
        if norm_v > 1e-12:
            v = v / np.sqrt(abs(norm_v))
        else:
            # Fallback: Frobenius norm
            fn = np.linalg.norm(v)
            if fn > 1e-12:
                v = v / fn
        ortho.append(v)
    return ortho


# =========================
# CASIMIR
# =========================

def compute_casimir(generators):
    return sum(G @ G for G in generators)


# =========================
# CORE OPTIMIZER
# =========================

def evolve_generators(
    n_gen=3,
    dim=2,
    steps=10000,
    lr=0.001,
    alpha=5.0,
    beta=1.0,
    killing_target=-1.0,
    initial_generators=None,
    verbose=True,
    convergence_threshold=1e-6,
):
    """
    Gradient-descent optimizer that evolves random hermitian matrices
    toward a Lie algebra structure.

    Returns (generators, history).
    """
    if initial_generators is not None:
        generators = [G.copy() for G in initial_generators]
        n_gen = len(generators)
        dim = generators[0].shape[0]
    else:
        generators = random_generators(n=n_gen, dim=dim)

    # Normalize initial
    for i in range(len(generators)):
        norm = np.linalg.norm(generators[i])
        if norm > 1e-10:
            generators[i] = generators[i] / norm

    K_target = killing_target * np.eye(n_gen)

    history = {
        "T_nc": [],
        "T_k": [],
        "T_norm": [],
        "T_total": [],
        "K_diag": [],
        "K_offdiag": [],
    }

    for step in range(steps):
        T_nc = non_commutativity_tension(generators)
        K = killing_metric(generators)
        T_k = np.linalg.norm(K - K_target) ** 2
        T_norm = sum((np.linalg.norm(G) - 1.0) ** 2 for G in generators)
        T_total = T_nc + alpha * T_k + beta * T_norm

        history["T_nc"].append(T_nc)
        history["T_k"].append(T_k)
        history["T_norm"].append(T_norm)
        history["T_total"].append(T_total)
        history["K_diag"].append(np.mean(np.diag(K)))
        history["K_offdiag"].append(np.mean(np.abs(K - np.diag(np.diag(K)))))

        # Numerical gradient
        grad = []
        epsilon = 1e-6

        for i, G in enumerate(generators):
            grad_G = np.zeros_like(G)
            for r in range(G.shape[0]):
                for c in range(G.shape[1]):
                    G_plus = G.copy()
                    G_plus[r, c] += epsilon
                    gens_plus = [g.copy() for g in generators]
                    gens_plus[i] = G_plus
                    K_plus = killing_metric(gens_plus)
                    T_plus = (
                        non_commutativity_tension(gens_plus)
                        + alpha * np.linalg.norm(K_plus - K_target) ** 2
                        + beta * sum((np.linalg.norm(Gp) - 1.0) ** 2 for Gp in gens_plus)
                    )

                    G_minus = G.copy()
                    G_minus[r, c] -= epsilon
                    gens_minus = [g.copy() for g in generators]
                    gens_minus[i] = G_minus
                    K_minus = killing_metric(gens_minus)
                    T_minus = (
                        non_commutativity_tension(gens_minus)
                        + alpha * np.linalg.norm(K_minus - K_target) ** 2
                        + beta * sum((np.linalg.norm(Gm) - 1.0) ** 2 for Gm in gens_minus)
                    )

                    grad_G[r, c] = (T_plus - T_minus) / (2 * epsilon)
            grad.append(grad_G)

        # Update
        for i in range(len(generators)):
            generators[i] = generators[i] - lr * grad[i]
            # Project to hermitian
            generators[i] = (generators[i] + generators[i].conj().T) / 2
            # Project to traceless
            generators[i] = generators[i] - np.trace(generators[i]) * np.eye(dim) / dim

        if verbose and step % 1000 == 0:
            print(
                f"Step {step}: T_total={T_total:.6f}  T_nc={T_nc:.6f}  T_k={T_k:.6f}  T_norm={T_norm:.6f}"
            )

        if T_total < convergence_threshold:
            if verbose:
                print(f"\nConverged at step {step}!")
            break

        # Early stopping: if T_total hasn't changed between consecutive steps
        if len(history["T_total"]) >= 2:
            prev = history["T_total"][-2]
            if abs(T_total - prev) < 1e-14:
                if verbose:
                    print(f"\nPlateau at step {step}: T_total={T_total:.10f}")
                break

    return generators, history


# =========================
# PARALLEL EXECUTION
# =========================

def _run_evolve_task(task):
    """Worker function for ProcessPoolExecutor."""
    import numpy as np
    seed = task.get('seed')
    if seed is not None:
        np.random.seed(seed)
    params = {k: v for k, v in task.items() if k != 'seed'}
    return evolve_generators(**params)


def evolve_parallel(tasks, max_workers=None):
    """Run multiple evolve_generators tasks in parallel.

    Each task is a dict with evolve_generators kwargs + optional 'seed' key.
    Returns list of (generators, history) tuples.
    """
    from concurrent.futures import ProcessPoolExecutor
    import os
    if max_workers is None:
        max_workers = min(len(tasks), os.cpu_count() or 4)

    # Limit numpy internal threading so process-level parallelism can use all cores
    thread_vars = ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS')
    old_env = {v: os.environ.get(v) for v in thread_vars}
    for v in thread_vars:
        os.environ[v] = '1'

    try:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(_run_evolve_task, tasks))
    finally:
        for v, val in old_env.items():
            if val is None:
                os.environ.pop(v, None)
            else:
                os.environ[v] = val

    return results
