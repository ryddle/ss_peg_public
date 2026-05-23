"""
GPU-accelerated utilities for the systematic verification plan.
Uses PyTorch with CUDA for autograd-based optimization and batched execution.

Key advantages over CPU version:
  1. Analytical gradients via torch.autograd (vs numerical finite differences)
  2. Batched execution of multiple seeds on a single GPU tensor
  3. Native GPU matrix operations for Killing metric, commutators, etc.
"""
import torch
import numpy as np


DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
DTYPE = torch.complex128   # Match CPU version precision
FDTYPE = torch.float64


def get_device_info():
    if torch.cuda.is_available():
        name = torch.cuda.get_device_name(0)
        mem = torch.cuda.get_device_properties(0).total_memory / 1e9
        return f"CUDA: {name} ({mem:.1f} GB)"
    return "CPU (no CUDA available)"


# =========================
# BASIC ALGEBRA (single)
# =========================

def pauli_matrices_torch(device=DEVICE):
    sigma_x = torch.tensor([[0, 1], [1, 0]], dtype=DTYPE, device=device)
    sigma_y = torch.tensor([[0, -1j], [1j, 0]], dtype=DTYPE, device=device)
    sigma_z = torch.tensor([[1, 0], [0, -1]], dtype=DTYPE, device=device)
    return [sigma_x, sigma_y, sigma_z]


def gell_mann_matrices_torch(device=DEVICE):
    l1 = torch.tensor([[0,1,0],[1,0,0],[0,0,0]], dtype=DTYPE, device=device)
    l2 = torch.tensor([[0,-1j,0],[1j,0,0],[0,0,0]], dtype=DTYPE, device=device)
    l3 = torch.tensor([[1,0,0],[0,-1,0],[0,0,0]], dtype=DTYPE, device=device)
    l4 = torch.tensor([[0,0,1],[0,0,0],[1,0,0]], dtype=DTYPE, device=device)
    l5 = torch.tensor([[0,0,-1j],[0,0,0],[1j,0,0]], dtype=DTYPE, device=device)
    l6 = torch.tensor([[0,0,0],[0,0,1],[0,1,0]], dtype=DTYPE, device=device)
    l7 = torch.tensor([[0,0,0],[0,0,-1j],[0,1j,0]], dtype=DTYPE, device=device)
    l8 = torch.tensor([[1,0,0],[0,1,0],[0,0,-2]], dtype=DTYPE, device=device) / np.sqrt(3)
    return [l1, l2, l3, l4, l5, l6, l7, l8]


# =========================
# BATCHED ALGEBRA
# =========================
# All batched functions operate on tensors of shape (B, n_gen, dim, dim)
# where B is the batch size.

def batched_commutator(A, B):
    """Commutator [A, B] = AB - BA. A, B: (..., dim, dim)"""
    return A @ B - B @ A


def batched_killing_metric(generators):
    """Compute Killing metric for a batch of generator sets (vectorized).

    generators: (B, n_gen, dim, dim)
    Returns: (B, n_gen, n_gen) real tensor

    Fully vectorized: computes all n_gen² commutators in one batched op,
    then contracts with a single einsum.  No Python loops over generators.
    """
    # G_i: (B, n, 1, d, d),  G_c: (B, 1, n, d, d)
    G_i = generators.unsqueeze(2)
    G_c = generators.unsqueeze(1)
    # ad[i,c] = [G_i, G_c] → (B, n, n, d, d)
    ad = G_i @ G_c - G_c @ G_i

    # K[i,j] = Σ_c Tr(ad[i,c] @ ad[j,c])
    #        = Σ_c Σ_k ad[i,c,k,l] * ad[j,c,l,k]
    K = torch.einsum('bickl,bjclk->bij', ad, ad).real
    return K


def batched_non_commutativity_tension(generators):
    """T_nc = sum_{i<j} ||[G_i, G_j]||^2 (vectorized).

    generators: (B, n_gen, dim, dim)
    Returns: (B,) real tensor
    """
    # Compute all commutators [G_i, G_j] at once
    G_i = generators.unsqueeze(2)  # (B, n, 1, d, d)
    G_j = generators.unsqueeze(1)  # (B, 1, n, d, d)
    C = G_i @ G_j - G_j @ G_i     # (B, n, n, d, d)

    # ||C||² per pair
    C_sq = (C.abs() ** 2).sum(dim=(-2, -1))  # (B, n, n)

    # Sum only upper triangle (i < j)
    idx = torch.triu_indices(generators.shape[1], generators.shape[1], offset=1,
                             device=generators.device)
    return C_sq[:, idx[0], idx[1]].sum(dim=-1)


def batched_norm_tension(generators, target_norm=1.0):
    """T_norm = sum_i (||G_i|| - target)^2.

    generators: (B, n_gen, dim, dim)
    Returns: (B,) real tensor
    """
    # Frobenius norm per generator: (B, n_gen)
    norms = (generators.abs() ** 2).sum(dim=(-2, -1)).sqrt()
    return ((norms - target_norm) ** 2).sum(dim=-1)


def batched_total_tension(generators, alpha=5.0, beta=1.0, killing_target=-1.0):
    """Compute total tension for a batch.

    generators: (B, n_gen, dim, dim) complex tensor (needs grad)
    Returns: (B,) real tensor
    """
    B, n, d, _ = generators.shape
    K_target = killing_target * torch.eye(n, dtype=FDTYPE, device=generators.device)

    T_nc = batched_non_commutativity_tension(generators)
    K = batched_killing_metric(generators)
    T_k = ((K - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
    T_norm = batched_norm_tension(generators)

    return T_nc + alpha * T_k + beta * T_norm


# Try to apply torch.compile for kernel fusion (PyTorch 2+)
try:
    batched_total_tension = torch.compile(batched_total_tension, mode='reduce-overhead')
    print("[shared_utils_gpu] torch.compile applied to batched_total_tension")
except Exception:
    print("[shared_utils_gpu] torch.compile not available, using eager mode")


# =========================
# PARAMETERIZATION
# =========================
# We parameterize hermitian traceless matrices via their real degrees of freedom.
# A d×d hermitian traceless matrix has d²-1 real parameters.
# This enables clean torch.autograd differentiation.

def params_to_generators(params, n_gen, dim):
    """Convert real parameters to hermitian traceless generators.

    params: (B, n_gen, dim*dim - 1) real tensor
    Returns: (B, n_gen, dim, dim) complex tensor
    """
    B = params.shape[0]
    generators = torch.zeros(B, n_gen, dim, dim, dtype=DTYPE, device=params.device)

    for g in range(n_gen):
        p = params[:, g, :]  # (B, dim*dim - 1)
        idx = 0
        H = torch.zeros(B, dim, dim, dtype=DTYPE, device=params.device)

        # Diagonal (d-1 independent values, last one is determined by tracelessness)
        for i in range(dim - 1):
            H[:, i, i] = p[:, idx].to(DTYPE)
            idx += 1
        H[:, dim-1, dim-1] = -H[:, :dim-1, :dim-1].diagonal(dim1=-2, dim2=-1).sum(dim=-1).to(DTYPE)

        # Off-diagonal: real and imaginary parts
        for i in range(dim):
            for j in range(i + 1, dim):
                re = p[:, idx]
                im = p[:, idx + 1]
                H[:, i, j] = (re + 1j * im).to(DTYPE)
                H[:, j, i] = (re - 1j * im).to(DTYPE)
                idx += 2

        generators[:, g] = H

    return generators


def generators_to_params(generators):
    """Extract real parameters from hermitian traceless generators.

    generators: (B, n_gen, dim, dim) complex tensor
    Returns: (B, n_gen, dim*dim - 1) real tensor
    """
    B, n_gen, dim, _ = generators.shape
    n_params = dim * dim - 1
    params = torch.zeros(B, n_gen, n_params, dtype=FDTYPE, device=generators.device)

    for g in range(n_gen):
        H = generators[:, g]  # (B, dim, dim)
        idx = 0
        for i in range(dim - 1):
            params[:, g, idx] = H[:, i, i].real
            idx += 1
        for i in range(dim):
            for j in range(i + 1, dim):
                params[:, g, idx] = H[:, i, j].real
                params[:, g, idx + 1] = H[:, i, j].imag
                idx += 2

    return params


def random_params(batch_size, n_gen, dim, device=DEVICE, seeds=None):
    """Generate random parameters for hermitian traceless generators.

    If seeds is a list of ints, each batch element uses a different numpy seed
    to ensure reproducibility matching the CPU version.
    """
    n_params = dim * dim - 1
    if seeds is not None:
        params_list = []
        for s in seeds:
            np.random.seed(s)
            p = np.random.randn(1, n_gen, n_params)
            params_list.append(p)
        params_np = np.concatenate(params_list, axis=0)
        params = torch.tensor(params_np, dtype=FDTYPE, device=device)
    else:
        params = torch.randn(batch_size, n_gen, n_params, dtype=FDTYPE, device=device)
    return params


# =========================
# CORE GPU OPTIMIZER
# =========================

def evolve_generators_gpu(
    n_gen=3,
    dim=2,
    steps=10000,
    lr=0.001,
    alpha=5.0,
    beta=1.0,
    killing_target=-1.0,
    batch_size=1,
    seeds=None,
    verbose=True,
    convergence_threshold=1e-6,
    initial_params=None,
):
    """
    GPU-accelerated optimizer using torch.autograd for analytical gradients.

    Can run multiple seeds simultaneously in a single batched tensor.
    Returns (final_generators_np, histories).
    """
    device = DEVICE

    if initial_params is not None:
        params = initial_params.clone().detach().to(device).requires_grad_(True)
        batch_size = params.shape[0]
    else:
        params = random_params(batch_size, n_gen, dim, device=device, seeds=seeds)
        # Normalize initial generators
        with torch.no_grad():
            gens = params_to_generators(params, n_gen, dim)
            norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt()  # (B, n_gen)
            norms = norms.clamp(min=1e-10)
            gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
            params = generators_to_params(gens)
        params = params.clone().detach().requires_grad_(True)

    optimizer = torch.optim.Adam([params], lr=lr)

    n_target = n_gen
    K_target_mat = killing_target * torch.eye(n_target, dtype=FDTYPE, device=device)

    history = {
        "T_total": [],
        "T_nc": [],
        "T_k": [],
        "K_diag": [],
        "K_offdiag": [],
    }

    prev_total = None

    for step in range(steps):
        optimizer.zero_grad()

        generators = params_to_generators(params, n_gen, dim)

        T_total = batched_total_tension(generators, alpha=alpha, beta=beta,
                                        killing_target=killing_target)

        loss = T_total.sum()
        loss.backward()
        optimizer.step()

        with torch.no_grad():
            # Record history (mean across batch)
            t_total_mean = T_total.mean().item()
            t_nc = batched_non_commutativity_tension(generators).mean().item()
            K = batched_killing_metric(generators)
            K_mean = K.mean(dim=0)
            t_k_mean = ((K - K_target_mat.unsqueeze(0)) ** 2).sum(dim=(-2, -1)).mean().item()
            k_diag = torch.diagonal(K_mean, dim1=-2, dim2=-1).mean().item()
            k_offdiag = (K_mean - torch.diag(torch.diagonal(K_mean))).abs().mean().item()

            history["T_total"].append(t_total_mean)
            history["T_nc"].append(t_nc)
            history["T_k"].append(t_k_mean)
            history["K_diag"].append(k_diag)
            history["K_offdiag"].append(k_offdiag)

            if verbose and step % 1000 == 0:
                print(f"  Step {step}: T_total={t_total_mean:.6f} T_nc={t_nc:.6f} "
                      f"T_k={t_k_mean:.6f} K_diag={k_diag:.4f}")

            if t_total_mean < convergence_threshold:
                if verbose:
                    print(f"  Converged at step {step}!")
                break

            # Early stopping: plateau detection
            if prev_total is not None and abs(t_total_mean - prev_total) < 1e-14:
                if verbose:
                    print(f"  Plateau at step {step}: T_total={t_total_mean:.10f}")
                break

            prev_total = t_total_mean

    # Extract final generators as numpy
    with torch.no_grad():
        final_gens = params_to_generators(params, n_gen, dim)
        final_np = final_gens.cpu().numpy()

    return final_np, history


def evolve_generators_gpu_individual(
    n_gen=3,
    dim=2,
    steps=10000,
    lr=0.001,
    alpha=5.0,
    beta=1.0,
    killing_target=-1.0,
    batch_size=1,
    seeds=None,
    verbose=True,
    convergence_threshold=1e-6,
):
    """
    GPU optimizer that returns per-seed results (individual histories).

    Runs all seeds in parallel on GPU but tracks per-seed metrics.
    Returns list of (generators_np, history) tuples (same format as CPU version).
    """
    device = DEVICE

    params = random_params(batch_size, n_gen, dim, device=device, seeds=seeds)
    # Normalize initial generators
    with torch.no_grad():
        gens = params_to_generators(params, n_gen, dim)
        norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt()
        norms = norms.clamp(min=1e-10)
        gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
        params = generators_to_params(gens)
    params = params.clone().detach().requires_grad_(True)

    optimizer = torch.optim.Adam([params], lr=lr)

    K_target_mat = killing_target * torch.eye(n_gen, dtype=FDTYPE, device=device)

    # Per-seed tracking
    per_seed = [{"T_total": [], "T_nc": [], "T_k": [], "T_norm": [],
                 "K_diag": [], "K_offdiag": []} for _ in range(batch_size)]

    prev_totals = torch.zeros(batch_size, device=device)
    active = torch.ones(batch_size, dtype=torch.bool, device=device)

    for step in range(steps):
        optimizer.zero_grad()

        generators = params_to_generators(params, n_gen, dim)

        T_total = batched_total_tension(generators, alpha=alpha, beta=beta,
                                        killing_target=killing_target)

        # Only backprop for active seeds
        loss = (T_total * active.float()).sum()
        loss.backward()
        optimizer.step()

        with torch.no_grad():
            T_nc = batched_non_commutativity_tension(generators)
            K = batched_killing_metric(generators)
            T_k = ((K - K_target_mat.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
            T_norm = batched_norm_tension(generators)

            for b in range(batch_size):
                if not active[b]:
                    continue
                per_seed[b]["T_total"].append(T_total[b].item())
                per_seed[b]["T_nc"].append(T_nc[b].item())
                per_seed[b]["T_k"].append(T_k[b].item())
                per_seed[b]["T_norm"].append(T_norm[b].item())
                per_seed[b]["K_diag"].append(
                    torch.diagonal(K[b]).mean().item())
                per_seed[b]["K_offdiag"].append(
                    (K[b] - torch.diag(torch.diagonal(K[b]))).abs().mean().item())

            # Check convergence/plateau per seed
            converged = T_total < convergence_threshold
            if step > 0:
                plateau = (T_total - prev_totals).abs() < 1e-14
                newly_done = (converged | plateau) & active
                active = active & ~newly_done

            prev_totals = T_total.clone()

            if verbose and step % 2000 == 0:
                n_active = active.sum().item()
                mean_t = T_total[active].mean().item() if n_active > 0 else 0
                print(f"  Step {step}: active={n_active}/{batch_size} "
                      f"mean_T={mean_t:.6f}")

            if not active.any():
                if verbose:
                    print(f"  All seeds converged/plateau at step {step}")
                break

    # Extract results
    with torch.no_grad():
        final_gens = params_to_generators(params, n_gen, dim)
        final_np = final_gens.cpu().numpy()

    results = []
    for b in range(batch_size):
        results.append((final_np[b], per_seed[b]))

    return results


# =========================
# CPU-COMPATIBLE ANALYSIS FUNCTIONS
# =========================
# These operate on numpy arrays (output of the GPU optimizer)
# to reuse the same verification logic as the CPU tests.

def killing_metric_np(generators):
    """Killing metric from numpy generators. Same as CPU shared_utils."""
    n = len(generators)
    K = np.zeros((n, n))
    for i in range(n):
        ad_i = [generators[i] @ g - g @ generators[i] for g in generators]
        for j in range(n):
            ad_j = [generators[j] @ g - g @ generators[j] for g in generators]
            K[i, j] = sum(np.trace(a @ b).real for a, b in zip(ad_i, ad_j))
    return K


def check_in_span_np(vec, basis, tol=1e-6):
    B = np.array([b.flatten() for b in basis]).T
    v = vec.flatten()
    coeffs, _, _, _ = np.linalg.lstsq(B, v, rcond=None)
    return np.linalg.norm(B @ coeffs - v) < tol


def commutator_np(A, B):
    return A @ B - B @ A


def structure_constants_real_np(generators):
    """Structure constants for hermitian generators: [T_a, T_b] = i·f_{abc}·T_c.

    Only valid when generators are hermitian (SU(N) convention).
    """
    n = len(generators)
    norms_sq = np.array([np.trace(generators[k] @ generators[k]).real for k in range(n)])
    f = np.zeros((n, n, n))
    for i in range(n):
        for j in range(n):
            C = generators[i] @ generators[j] - generators[j] @ generators[i]
            for k in range(n):
                if abs(norms_sq[k]) > 1e-12:
                    f[i, j, k] = (-1j * np.trace(C @ generators[k]) / norms_sq[k]).real
    return f


def structure_constants_general_np(generators):
    """Structure constants for any convention: [T_a, T_b] = c_{abc}·T_c.

    Works for both hermitian (SU) and anti-hermitian (SO) generators.
    For hermitian gens, c_{abc} = i·f_{abc} where f is the standard physics constant.
    """
    n = len(generators)
    norms_sq = np.array([np.trace(generators[k] @ generators[k]).real for k in range(n)])
    c = np.zeros((n, n, n))
    for i in range(n):
        for j in range(n):
            C = generators[i] @ generators[j] - generators[j] @ generators[i]
            for k in range(n):
                if abs(norms_sq[k]) > 1e-12:
                    c[i, j, k] = (np.trace(C @ generators[k]) / norms_sq[k]).real
    return c


def compute_casimir_np(generators):
    return sum(G @ G for G in generators)


def orthonormalize_killing_np(generators):
    n = len(generators)
    ortho = []
    for i in range(n):
        v = generators[i].copy()
        for u in ortho:
            ip = sum(np.trace(commutator_np(u, g) @ commutator_np(v, g)).real
                     for g in generators)
            norm_u = sum(np.trace(commutator_np(u, g) @ commutator_np(u, g)).real
                         for g in generators)
            if abs(norm_u) > 1e-12:
                v = v - (ip / norm_u) * u
        norm_v = sum(np.trace(commutator_np(v, g) @ commutator_np(v, g)).real
                     for g in generators)
        if norm_v > 1e-12:
            v = v / np.sqrt(abs(norm_v))
        else:
            fn = np.linalg.norm(v)
            if fn > 1e-12:
                v = v / fn
        ortho.append(v)
    return ortho


# =========================
# SO(n) PARAMETERIZATION
# =========================
# Real antisymmetric d×d matrices: A^T = -A, purely real.
# d(d-1)/2 parameters per generator (upper triangle entries).

def params_to_antisymmetric(params, n_gen, dim):
    """Convert real parameters to real antisymmetric generators (vectorized).

    params: (B, n_gen, dim*(dim-1)//2) real tensor
    Returns: (B, n_gen, dim, dim) complex tensor (purely real, A^T = -A)
    """
    B = params.shape[0]
    generators = torch.zeros(B, n_gen, dim, dim, dtype=DTYPE, device=params.device)

    # Build upper-triangle index arrays once
    row_idx, col_idx = torch.triu_indices(dim, dim, offset=1, device=params.device)
    # params already ordered as (i,j) for i<j in row-major order
    vals = params.to(DTYPE)  # (B, n_gen, n_params)
    generators[:, :, row_idx, col_idx] = vals
    generators[:, :, col_idx, row_idx] = -vals

    return generators


def antisymmetric_to_params(generators):
    """Extract real parameters from antisymmetric generators (vectorized).

    generators: (B, n_gen, dim, dim)
    Returns: (B, n_gen, dim*(dim-1)//2) real tensor
    """
    B, n_gen, dim, _ = generators.shape

    row_idx, col_idx = torch.triu_indices(dim, dim, offset=1, device=generators.device)
    return generators[:, :, row_idx, col_idx].real


def random_params_so(batch_size, n_gen, dim, device=DEVICE, seeds=None):
    """Random parameters for SO(n) antisymmetric generators."""
    n_params = dim * (dim - 1) // 2
    if seeds is not None:
        params_list = []
        for s in seeds:
            np.random.seed(s)
            p = np.random.randn(1, n_gen, n_params)
            params_list.append(p)
        params_np = np.concatenate(params_list, axis=0)
        return torch.tensor(params_np, dtype=FDTYPE, device=device)
    return torch.randn(batch_size, n_gen, n_params, dtype=FDTYPE, device=device)


def n_params_so(dim):
    """Number of parameters per SO(dim) generator."""
    return dim * (dim - 1) // 2


# =========================
# Sp(2n) PARAMETERIZATION
# =========================
# Compact symplectic algebra sp(2n): X ∈ M(2n,C) with X^†Ω + ΩX = 0
# where Ω = [[0, I_n], [-I_n, 0]].
# Equivalently: X = [[A, B], [C, -A^T]] with B = B^T, C = C^T, A arbitrary,
# and additionally X^† = -X (anti-hermitian → generators are i·X, hermitian).
# For compact Sp(n) ≅ USp(2n): dim = n(2n+1) generators.
#
# Parameterization: we directly enforce the block structure.
#   A (n×n complex): n² real params (real + imag, but A is constrained by anti-hermiticity)
#   B (n×n symmetric complex): n(n+1)/2 complex = n(n+1) real params
#   (C determined by anti-hermiticity: C = -B̄)
#
# Actually for compact form (USp(2n)):
#   X = [[A, B], [-B̄, -Ā]] with A^† = -A (anti-hermitian), B^T = B (symmetric)
#   So A has n²-1 independent real params (anti-hermitian traceless) + 1 trace? 
#   No: anti-hermitian n×n has n² real params (n diagonal imaginary + n(n-1) off-diag pairs)
#   B symmetric complex n×n has n(n+1)/2 complex = n(n+1) real params
#   Total: n² + n(n+1) = n(2n+1) ✓
#
# We parameterize hermitian generators T = iX where X ∈ sp(2n).

def _build_omega(n, device=DEVICE):
    """Build the 2n×2n symplectic form Ω = [[0, I], [-I, 0]]."""
    I = torch.eye(n, dtype=DTYPE, device=device)
    Z = torch.zeros(n, n, dtype=DTYPE, device=device)
    return torch.cat([torch.cat([Z, I], dim=1),
                      torch.cat([-I, Z], dim=1)], dim=0)


def params_to_symplectic(params, n_gen, half_dim):
    """Convert real parameters to USp(2n) hermitian generators.

    half_dim = n, so the matrix dimension is 2n.
    Each generator is a 2n×2n hermitian matrix T such that
    T·Ω + Ω·T^T = 0 (i.e., iT ∈ sp(2n)).

    params: (B, n_gen, n_params_sp) real tensor
    Returns: (B, n_gen, 2n, 2n) complex tensor
    """
    n = half_dim
    dim = 2 * n
    B = params.shape[0]
    generators = torch.zeros(B, n_gen, dim, dim, dtype=DTYPE, device=params.device)

    # Each generator has n² + n(n+1) = n(2n+1) real params
    for g in range(n_gen):
        p = params[:, g, :]
        idx = 0

        # Block A: anti-hermitian n×n (n² real params)
        # → diagonal: n purely imaginary values
        # → upper triangle: n(n-1)/2 complex values (re, im)
        A = torch.zeros(B, n, n, dtype=DTYPE, device=params.device)
        for i in range(n):
            A[:, i, i] = (1j * p[:, idx]).to(DTYPE)
            idx += 1
        for i in range(n):
            for j in range(i + 1, n):
                re = p[:, idx]
                im = p[:, idx + 1]
                A[:, i, j] = (re + 1j * im).to(DTYPE)
                A[:, j, i] = (-re + 1j * im).to(DTYPE)
                idx += 2

        # Block B: symmetric complex n×n (n(n+1) real params)
        B_block = torch.zeros(B, n, n, dtype=DTYPE, device=params.device)
        for i in range(n):
            re = p[:, idx]
            im = p[:, idx + 1]
            B_block[:, i, i] = (re + 1j * im).to(DTYPE)
            idx += 2
        for i in range(n):
            for j in range(i + 1, n):
                re = p[:, idx]
                im = p[:, idx + 1]
                B_block[:, i, j] = (re + 1j * im).to(DTYPE)
                B_block[:, j, i] = (re + 1j * im).to(DTYPE)  # symmetric
                idx += 2

        # Construct X = [[A, B], [-B̄, Ā]] (anti-hermitian in usp(2n))
        # This automatically satisfies X†=-X and X^TΩ+ΩX=0.
        X = torch.zeros(B, dim, dim, dtype=DTYPE, device=params.device)
        X[:, :n, :n] = A
        X[:, :n, n:] = B_block
        X[:, n:, :n] = -B_block.conj()
        X[:, n:, n:] = A.conj()

        # Hermitian generator T = iX
        generators[:, g] = 1j * X

    return generators


def symplectic_to_params(generators, half_dim):
    """Extract parameters from USp(2n) hermitian generators.

    generators: (B, n_gen, 2n, 2n) complex tensor
    Returns: (B, n_gen, n_params_sp) real tensor
    """
    n = half_dim
    B_batch, n_gen, dim, _ = generators.shape
    n_params = n * (2 * n + 1)
    params = torch.zeros(B_batch, n_gen, n_params, dtype=FDTYPE, device=generators.device)

    for g in range(n_gen):
        # T = iX → X = -iT = T/(-i) = -i·T
        X = -1j * generators[:, g]  # (B, 2n, 2n)
        A = X[:, :n, :n]
        B_block = X[:, :n, n:]

        idx = 0
        # A diagonal: imaginary parts
        for i in range(n):
            params[:, g, idx] = A[:, i, i].imag
            idx += 1
        # A upper triangle
        for i in range(n):
            for j in range(i + 1, n):
                params[:, g, idx] = A[:, i, j].real
                params[:, g, idx + 1] = A[:, i, j].imag
                idx += 2
        # B diagonal
        for i in range(n):
            params[:, g, idx] = B_block[:, i, i].real
            params[:, g, idx + 1] = B_block[:, i, i].imag
            idx += 2
        # B upper triangle
        for i in range(n):
            for j in range(i + 1, n):
                params[:, g, idx] = B_block[:, i, j].real
                params[:, g, idx + 1] = B_block[:, i, j].imag
                idx += 2

    return params


def random_params_sp(batch_size, n_gen, half_dim, device=DEVICE, seeds=None):
    """Random parameters for USp(2n) generators."""
    n = half_dim
    n_params = n * (2 * n + 1)
    if seeds is not None:
        params_list = []
        for s in seeds:
            np.random.seed(s)
            p = np.random.randn(1, n_gen, n_params)
            params_list.append(p)
        params_np = np.concatenate(params_list, axis=0)
        return torch.tensor(params_np, dtype=FDTYPE, device=device)
    return torch.randn(batch_size, n_gen, n_params, dtype=FDTYPE, device=device)


def n_params_sp(half_dim):
    """Number of parameters per USp(2·half_dim) generator."""
    return half_dim * (2 * half_dim + 1)


# =========================
# GENERIC GPU OPTIMIZER
# =========================

def evolve_generators_generic(
    n_gen,
    dim,
    steps=10000,
    lr=0.001,
    alpha=5.0,
    beta=1.0,
    killing_target=-1.0,
    batch_size=1,
    seeds=None,
    verbose=True,
    convergence_threshold=1e-6,
    param_to_gen_fn=None,
    gen_to_param_fn=None,
    random_param_fn=None,
    param_to_gen_kwargs=None,
):
    """Generic GPU optimizer that works with any parameterization.

    param_to_gen_fn(params, n_gen, ...) → generators (B, n_gen, dim, dim)
    gen_to_param_fn(generators, ...) → params
    random_param_fn(batch_size, n_gen, ..., device, seeds) → params

    Defaults to SU(N) hermitian traceless parameterization.
    Returns list of (generators_np, history) per seed.
    """
    device = DEVICE
    if param_to_gen_kwargs is None:
        param_to_gen_kwargs = {}

    if param_to_gen_fn is None:
        param_to_gen_fn = params_to_generators
        gen_to_param_fn = generators_to_params
        random_param_fn = lambda bs, ng, device, seeds: random_params(bs, ng, dim, device, seeds)
        if "dim" not in param_to_gen_kwargs:
            param_to_gen_kwargs["dim"] = dim

    # Initialize
    params = random_param_fn(batch_size, n_gen, device=device, seeds=seeds)

    # Normalize
    with torch.no_grad():
        gens = param_to_gen_fn(params, n_gen, **param_to_gen_kwargs)
        norms = (gens.abs() ** 2).sum(dim=(-2, -1)).sqrt()
        norms = norms.clamp(min=1e-10)
        gens = gens / norms.unsqueeze(-1).unsqueeze(-1)
        params = gen_to_param_fn(gens)
    params = params.clone().detach().requires_grad_(True)

    optimizer = torch.optim.Adam([params], lr=lr)
    K_target_mat = killing_target * torch.eye(n_gen, dtype=FDTYPE, device=device)

    per_seed = [{"T_total": [], "T_nc": [], "T_k": [], "T_norm": [],
                 "K_diag": [], "K_offdiag": []} for _ in range(batch_size)]

    prev_totals = torch.zeros(batch_size, device=device)
    active = torch.ones(batch_size, dtype=torch.bool, device=device)

    for step in range(steps):
        optimizer.zero_grad()
        generators = param_to_gen_fn(params, n_gen, **param_to_gen_kwargs)

        # Compute tension components once (with grad) — reuse for logging
        T_nc = batched_non_commutativity_tension(generators)
        K = batched_killing_metric(generators)
        T_k = ((K - K_target_mat.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
        T_norm = batched_norm_tension(generators)
        T_total = T_nc + alpha * T_k + beta * T_norm

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
                per_seed[b]["K_diag"].append(
                    torch.diagonal(K[b]).mean().item())
                per_seed[b]["K_offdiag"].append(
                    (K[b] - torch.diag(torch.diagonal(K[b]))).abs().mean().item())

            converged = T_total < convergence_threshold
            if step > 0:
                plateau = (T_total - prev_totals).abs() < 1e-14
                newly_done = (converged | plateau) & active
                active = active & ~newly_done

            prev_totals = T_total.clone()

            if verbose and step % 2000 == 0:
                n_active = active.sum().item()
                mean_t = T_total[active].mean().item() if n_active > 0 else 0
                print(f"  Step {step}: active={n_active}/{batch_size} "
                      f"mean_T={mean_t:.6f}")

            if not active.any():
                if verbose:
                    print(f"  All seeds converged/plateau at step {step}")
                break

    with torch.no_grad():
        final_gens = param_to_gen_fn(params, n_gen, **param_to_gen_kwargs)
        final_np = final_gens.cpu().numpy()

    results = []
    for b in range(batch_size):
        results.append((final_np[b], per_seed[b]))
    return results
