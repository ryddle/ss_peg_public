"""
Killing-Regularized Neural Networks -- Shared Utilities
=======================================================
Core idea:
  If a neural network layer has n_gen = d^2 - 1 neurons and d^2-dimensional
  inputs, each weight row can be reshaped into a dxd matrix. These d^2-1
  matrices are 'candidate Lie algebra generators'. Adding T_K (Killing tension)
  as a regularizer drives them toward genuine compact Lie algebra structure.

This directly transposes the SS-PEG variational result:
  'T_K + T_closure + gradient descent -> exact Lie algebra'
into a neural network training loop.

Adapted from verification/shared_utils_gpu.py.
"""

import functools
import torch
import numpy as np
from typing import Optional

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# ============================================================
#  CORE ALGEBRAIC FUNCTIONS (all differentiable via autograd)
# ============================================================

def killing_metric(generators: torch.Tensor) -> torch.Tensor:
    """
    K_ab = sum_c Tr([G_a, G_c] @ [G_b, G_c])

    This is the Killing form of the Lie algebra spanned by {G_a}.
    For any simple compact algebra: K = -c * I_n (negative definite, proportional to identity).

    generators : (n_gen, d, d)  real or complex
    Returns    : (n_gen, n_gen) real
    """
    G_i = generators.unsqueeze(1)   # (n, 1, d, d) -- index a
    G_c = generators.unsqueeze(0)   # (1, n, d, d) -- index c
    ad  = G_i @ G_c - G_c @ G_i    # (n, n, d, d) -- ad[a,c] = [G_a, G_c]

    # K[a,b] = Sigma_c Tr(ad[a,c] @ ad[b,c]) = Sigma_c Sigma_{k,l} ad[a,c,k,l] * ad[b,c,l,k]
    K = torch.einsum('ackl,bclk->ab', ad, ad)
    return K.real if generators.is_complex() else K


def killing_tension(generators: torch.Tensor) -> torch.Tensor:
    """
    T_K = || K_normalized + I_n ||^2

    Normalizes K so its diagonal mean equals -1, then penalizes deviation from -I.
    Zero iff K is proportional to -I (signature of a compact simple Lie algebra).

    generators : (n_gen, d, d)
    Returns    : scalar (differentiable)
    """
    K = killing_metric(generators)
    n = K.shape[0]
    # Normalize so eigenvalues target -1
    scale = K.diagonal().abs().mean().clamp(min=1e-8)
    K_norm = K / scale
    target = -torch.eye(n, dtype=K.dtype, device=K.device)
    return ((K_norm - target) ** 2).sum()


def closure_tension(generators: torch.Tensor) -> torch.Tensor:
    """
    T_closure = mean_{i<j} || [G_i,G_j] - proj_V([G_i,G_j]) ||^2 / || [G_i,G_j] ||^2

    Measures the fraction of each commutator that 'escapes' the span of the
    generators. Zero means the generators form a closed algebra.
    Uses SVD for a numerically stable orthonormal basis of the generator span.

    Works for both real and complex generators.

    generators : (n_gen, d, d)  real or complex
    Returns    : scalar in [0, 1] (differentiable, real-valued)

    Projection formula for complex vectors (standard Hermitian inner product):
      Q = Vh[:rank].mH   -- orthonormal columns  (Q.mH @ Q = I)
      coords  = C_flat @ Q.conj()   -- inner products <q_l, c>
      C_proj  = coords  @ Q.T       -- reconstruction in original space
      (equivalent to: C_proj = C_flat @ Q.conj() @ Q.T = C_flat @ P
       where P = Q.conj() @ Q.T is the orthogonal projection matrix)
    For real generators, Q.conj() = Q and Q.T = Q.T, recovering the standard formula.
    """
    n, d, _ = generators.shape

    # Flatten generators: (n, d^2)
    G_flat = generators.reshape(n, d * d)

    # Orthonormal column basis via SVD
    # Vh: (min(n,d^2), d^2) with orthonormal rows -- works for real and complex
    _, S, Vh = torch.linalg.svd(G_flat, full_matrices=False)
    threshold = S.max() * 1e-6
    rank      = int((S > threshold).sum().item())
    Q         = Vh[:rank].mH        # (d^2, rank) -- orthonormal cols (works real+complex)

    # Commutators [G_i, G_j] for i < j
    G_i = generators.unsqueeze(1)   # (n, 1, d, d)
    G_j = generators.unsqueeze(0)   # (1, n, d, d)
    C   = G_i @ G_j - G_j @ G_i    # (n, n, d, d)

    rows, cols = torch.triu_indices(n, n, offset=1, device=generators.device)
    C_flat = C[rows, cols].reshape(-1, d * d)   # (n_pairs, d^2)

    if C_flat.shape[0] == 0:
        return torch.zeros(1, dtype=torch.float64, device=generators.device).squeeze()

    # Hermitian projection onto span(G)
    coords   = C_flat @ Q.conj()    # (n_pairs, rank) -- inner products with columns of Q
    C_proj   = coords @ Q.T         # (n_pairs, d^2) -- reconstruction
    residual = C_flat - C_proj

    # Real squared norms (|*|^2 handles both real and complex uniformly)
    C_sq   = (C_flat.abs()   ** 2).sum(-1).clamp(min=1e-12)  # (n_pairs,) real
    res_sq = (residual.abs() ** 2).sum(-1)                    # (n_pairs,) real

    return (res_sq / C_sq).mean()


def norm_tension(generators: torch.Tensor, target: float = 1.0) -> torch.Tensor:
    """
    T_norm = Sigma_a (||G_a|| - target)^2

    Stabilizes generator norms during optimization (same role as beta in verif. suite).
    """
    norms = generators.reshape(generators.shape[0], -1).norm(dim=-1)
    return ((norms - target) ** 2).sum()


# ============================================================
#  LAYER BRIDGE: neural network weights <-> generators
# ============================================================

def extract_generators(weight: torch.Tensor, d: int) -> torch.Tensor:
    """
    Reshape weight matrix rows into dxd candidate generator matrices.

    Given W in R^(n_gen x n_in) where n_in >= d^2:
      G_a = W[a, :d^2].reshape(d, d)   for a = 0 .. n_gen-1

    This is the key bridge: weight rows ARE the generator parameters.

    weight : (n_gen, n_in)
    d      : target matrix dimension
    Returns: (n_gen, d, d)
    """
    n_gen, n_in = weight.shape
    if n_in < d * d:
        raise ValueError(
            f"Row dimension {n_in} < d^2 = {d*d}: "
            f"cannot embed rows into {d}x{d} matrices."
        )
    return weight[:, : d * d].reshape(n_gen, d, d).contiguous()


# ============================================================
#  COMBINED REGULARIZATION LOSS
# ============================================================

def killing_regularization(
    weight:          torch.Tensor,
    d:               int,
    lambda_killing:  float = 1.0,
    lambda_closure:  float = 0.0,
    lambda_norm:     float = 0.1,
) -> torch.Tensor:
    """
    Full regularization term to add to a task loss for any layer with n_in >= d^2.

      L_reg = lambda_K * T_K + lambda_c * T_closure + lambda_n * T_norm

    weight           : (n_gen, n_in) weight matrix of the layer
    d                : algebra representation dimension (rows -> dxd matrices)
    lambda_killing   : weight for Killing tension (structural driver)
    lambda_closure   : weight for closure tension (algebraic self-consistency)
    lambda_norm      : weight for norm stabilizer
    Returns          : scalar loss (differentiable)
    """
    G   = extract_generators(weight, d)
    reg = torch.zeros(1, dtype=weight.dtype, device=weight.device).squeeze()

    if lambda_killing > 0:
        reg = reg + lambda_killing * killing_tension(G)
    if lambda_closure > 0:
        reg = reg + lambda_closure * closure_tension(G)
    if lambda_norm > 0:
        reg = reg + lambda_norm * norm_tension(G)

    return reg


# ============================================================
#  POST-TRAINING ALGEBRAIC ANALYSIS
# ============================================================

def analyze_generators(generators: torch.Tensor) -> dict:
    """
    Compute the full algebraic fingerprint of a set of matrices:
      - Killing tension (how close to K propto -I)
      - Closure fraction (how closed the algebra is)
      - Killing eigenspectrum
      - Estimated dual Coxeter number h^v

    Used to verify that trained weight rows have crystallized into a Lie algebra.

    generators : (n_gen, d, d)
    Returns    : dict of scalar metrics + numpy arrays
    """
    with torch.no_grad():
        gens_f64 = generators.detach().double()
        K  = killing_metric(gens_f64)
        n  = K.shape[0]

        # Eigenspectrum
        eigvals = torch.linalg.eigvalsh(K.float()).cpu().numpy()

        # Normalized Killing tension
        scale  = K.diagonal().abs().mean().clamp(1e-8)
        K_norm = K / scale
        t_k = ((K_norm + torch.eye(n, dtype=K.dtype, device=K.device)) ** 2).sum().item()

        # Closure fraction
        t_c = closure_tension(gens_f64).item()

        # Eigenvalue uniformity: std/mean(abs)
        eig_abs = np.abs(eigvals)
        uniformity = float(eig_abs.std() / (eig_abs.mean() + 1e-8))

        # Rough h^v estimate via universal formula:
        # |K| / (n^{1/2} * mean||G||^2) ~ h^v  (in fundamental rep, normalized generators)
        frob_mean = float((gens_f64 ** 2).sum(dim=(-2, -1)).mean().item())
        k_frob = float(K.norm().item())
        h_estimate = k_frob / max(n ** 0.5 * frob_mean, 1e-8)

        return {
            "killing_tension":       t_k,
            "closure_fraction":      t_c,
            "killing_eigenvalues":   eigvals,
            "eigenvalue_uniformity": uniformity,         # 0 = perfectly uniform (-cI)
            "h_dual_estimate":       h_estimate,
            "is_compact_algebra":    bool(t_k < 1.0 and uniformity < 0.1),
            "is_closed":             bool(t_c < 0.05),
        }


def print_algebra_report(metrics: dict, title: str = "Algebraic Structure") -> None:
    """Pretty-print the output of analyze_generators."""
    w = 58
    print("\n" + "=" * w)
    print(f"  {title}")
    print("=" * w)
    print(f"  Killing tension T_K          : {metrics['killing_tension']:.4f}")
    print(f"  Closure fraction T_closure   : {metrics['closure_fraction']:.4f}")
    print(f"  Eigenvalue uniformity        : {metrics['eigenvalue_uniformity']:.4f}  (0 = perfect)")
    ev = np.round(metrics["killing_eigenvalues"], 3)
    print(f"  Killing eigenvalues          : {ev}")
    print(f"  K propto -I  (compact algebra)    : {metrics['is_compact_algebra']}")
    print(f"  Algebraically closed         : {metrics['is_closed']}")
    print(f"  Estimated h^v                : {metrics['h_dual_estimate']:.3f}")
    print("=" * w)


# ============================================================
#  SU(d) MANIFOLD PARAMETRIZATION  (for AlgebraicLinear, Exp 3)
# ============================================================

@functools.lru_cache(maxsize=32)
def _su_basis(d: int, device_str: str) -> torch.Tensor:
    """
    Precompute the canonical orthogonal hermitian traceless basis for su(d).

    Uses the generalized Gell-Mann convention:
      - Symmetric off-diagonal : e_{ij} + e_{ji}           for i < j
      - Antisymmetric off-diag : -i(e_{ij} - e_{ji})       for i < j
      - Diagonal Cartan        : H_l = sqrt(2/(l(l+1))) diag(1,...,1,-l,0,...,0)

    All d^2-1 matrices satisfy  Tr(T_a T_b) = 2 delta_{ab}  (GML normalization).
    This orthogonality guarantees that K = -2d*I for any orthonormal rotation
    of the basis -> T_K = 0 is exactly achievable.

    Cached per (d, device) since these are pure constants.
    Returns : (d^2-1, d, d) complex128 on the requested device.
    """
    device = torch.device(device_str)
    dtype  = torch.complex128
    basis  = []

    # Symmetric off-diagonal: e_{ij} + e_{ji}   [Tr = 2 ok]
    for i in range(d):
        for j in range(i + 1, d):
            H = torch.zeros(d, d, dtype=dtype, device=device)
            H[i, j] = 1.0
            H[j, i] = 1.0
            basis.append(H)

    # Antisymmetric off-diagonal: -i*e_{ij} + i*e_{ji}   [Tr = 2 ok]
    for i in range(d):
        for j in range(i + 1, d):
            H = torch.zeros(d, d, dtype=dtype, device=device)
            H[i, j] = -1j
            H[j, i] =  1j
            basis.append(H)

    # Diagonal Cartan generators H_l, l = 1 ... d-1
    # H_l = sqrt(2/(l(l+1))) * diag(1,...,1, -l, 0,...,0)   [l ones, then -l]
    # Tr(H_l^2) = (2/(l(l+1))) * (l*1 + l^2) = (2/(l(l+1))) * l(l+1) = 2 ok
    for l in range(1, d):
        scale      = float(np.sqrt(2.0 / (l * (l + 1))))
        diag_vals  = [scale] * l + [-l * scale] + [0.0] * (d - l - 1)
        H = torch.diag(torch.tensor(diag_vals, dtype=dtype, device=device))
        basis.append(H)

    result = torch.stack(basis)   # (d^2-1, d, d)
    assert result.shape[0] == d * d - 1, f"su({d}) basis: expected {d*d-1} elements, got {result.shape[0]}"
    return result


def params_to_su_generators(params: torch.Tensor, d: int) -> torch.Tensor:
    """
    Convert real parameters to hermitian traceless su(d) generators via a fixed basis.

    params  : (n_gen, d^2-1) real  [learnable, gradient-tracked]
    d       : representation dimension
    Returns : (n_gen, d, d) complex hermitian traceless

    Each generator G_a = sum_k params[a, k] * basis_k
    where {basis_k} is the orthogonal GML basis from _su_basis().

    CRITICAL PROPERTY vs Exp 1 & 2:
      The output ALWAYS lives inside su(d) (hermitian traceless).
      T_K CAN therefore reach exactly 0, since the global minimum
      K = -2d*I is achievable within this parameterization.
      In Exp 1/2, weight rows were arbitrary real matrices -> wrong manifold
      -> T_K stuck at ~26 and could never reach 0.
    """
    basis    = _su_basis(d, str(params.device))          # (d^2-1, d, d) complex, const
    params_c = params.to(torch.complex128)               # (n_gen, d^2-1) complex
    return torch.einsum('ak,kij->aij', params_c, basis)  # (n_gen, d, d)


# ============================================================
#  DATASET UTILITIES
# ============================================================

def generate_cubic_trace_dataset(
    n:      int   = 20_000,
    d:      int   = 3,
    seed:   int   = 42,
    device: torch.device = DEVICE,
) -> tuple:
    """
    SL(d)-variant classification task.

    Input : x in R^{d^2}  (flattened dxd real matrix M)
    Label : sign(Tr(M^3)) > 0  -- a degree-3 polynomial invariant under conjugacy

    This task requires the network to compute a non-linear algebraic function
    of the input. A layer with su(d)-structured weights can represent this
    more efficiently via the adjoint action.

    Returns: (X_tensor, y_tensor) on the given device
    """
    rng = np.random.default_rng(seed)
    X   = rng.standard_normal((n, d * d)).astype(np.float64)
    M   = X.reshape(-1, d, d)
    M3  = np.einsum('nij,njk,nkl->nil', M, M, M)  # M^3 batched
    cubic_trace = np.trace(M3, axis1=1, axis2=2)    # Tr(M^3)
    y   = (cubic_trace > 0).astype(np.int64)

    return (
        torch.tensor(X, dtype=torch.float64, device=device),
        torch.tensor(y, dtype=torch.long,   device=device),
    )


def generate_group_action_dataset(
    n:      int   = 10_000,
    d:      int   = 2,
    seed:   int   = 0,
    device: torch.device = DEVICE,
) -> tuple:
    """
    Group-action dataset for symmetry discovery (Exp 2).

    Samples pairs (x, y) where y = U x and U in SU(d) is a random group element.
    The task: discover the d^2-1 generators of su(d) from these transformation pairs.

    Returns: (X, Y, true_generators) -- X: (n, 2*d^2), Y: (n, 2*d), gens: (d^2-1, d, d) complex
             X encodes (Re(x), Im(x), Re(U-I), Im(U-I)) for the learner.
    """
    import scipy.stats as stats

    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)

    # True su(2) generators (Pauli matrices / 2i)
    if d == 2:
        i = 1j
        G1 = np.array([[0, 1], [1, 0]],  dtype=complex) * 0.5j
        G2 = np.array([[0, -i], [i, 0]], dtype=complex) * 0.5j
        G3 = np.array([[1, 0], [0, -1]], dtype=complex) * 0.5j
        true_generators = np.stack([G1, G2, G3])  # (3, 2, 2) complex
    else:
        raise NotImplementedError(f"group_action_dataset: d={d} not implemented yet (only d=2)")

    # Random unit vectors in C^d
    X_re = rng.standard_normal((n, d)).astype(np.float64)
    X_im = rng.standard_normal((n, d)).astype(np.float64)
    X_c  = X_re + 1j * X_im
    X_c  = X_c / np.linalg.norm(X_c, axis=-1, keepdims=True)

    # Random SU(2) elements via random axis-angle
    angles  = rng.uniform(0, 2 * np.pi, n)
    axes_re = rng.standard_normal((n, 3))
    axes    = axes_re / np.linalg.norm(axes_re, axis=-1, keepdims=True)

    # U = exp(theta * (a1 G1 + a2 G2 + a3 G3))
    Y_c = np.zeros_like(X_c)
    for k in range(n):
        theta = angles[k]
        a     = axes[k]
        H     = a[0] * G1 + a[1] * G2 + a[2] * G3  # element of su(2)
        U     = np.eye(d, dtype=complex)
        Hn    = np.eye(d, dtype=complex)
        for m in range(1, 20):             # exp(theta*H) via Taylor series
            Hn = Hn @ (theta * H) / m
            U  = U + Hn
        Y_c[k] = U @ X_c[k]

    # Encode as real vectors: [Re(x), Im(x), Re(y), Im(y)]
    X_enc = np.concatenate([X_re, X_im, Y_c.real, Y_c.imag], axis=-1).astype(np.float64)

    # Target: (a1, a2, a3, theta) -- the generator coefficients
    Y_enc = np.concatenate([axes * angles[:, None]], axis=-1).astype(np.float64)

    true_gen_t = torch.tensor(true_generators, dtype=torch.complex128, device=device)

    return (
        torch.tensor(X_enc, dtype=torch.float64, device=device),
        torch.tensor(Y_enc, dtype=torch.float64, device=device),
        true_gen_t,
    )
