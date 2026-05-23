"""
Octonion algebra and Jordan algebra J₃(𝕆) for exceptional Lie algebras.

Provides:
  - Full octonion multiplication table (8D: e₀=1, e₁..e₇ imaginary)
  - Conjugation, norm, non-associative product
  - J₃(𝕆): 27-dimensional basis of 3×3 hermitian octonionic matrices
  - Jordan product X∘Y = ½(XY + YX) structure constants
  - Cubic norm tensor N(X) = det(X) for J₃(𝕆)
  - F₄ construction: Der(J₃(𝕆)) via derivation condition
  - E₆ construction: cubic norm preservation δN = 0

References:
  - Baez, "The Octonions" (2002), Bull. AMS 39, 145-205
  - Yokota, "Exceptional Lie Groups" (2009)
"""
import numpy as np


# ===================================================================
# OCTONION MULTIPLICATION TABLE
# ===================================================================
# Convention: e₀ = 1 (identity), e₁..e₇ imaginary units
# Fano plane triples: e_a · e_b = e_c for (a,b,c) in triples
# Same Fano plane as test_g2_gpu.py but with e₀ = 1 added.

FANO_TRIPLES = [
    (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7),
    (5, 6, 1), (6, 7, 2), (7, 1, 3)
]

# Build full multiplication table: e_a * e_b = MULT[a,b,c] * e_c
# MULT[a,b] is an 8-vector of coefficients in basis {e₀..e₇}
MULT_TABLE = np.zeros((8, 8, 8))

# e₀ * e_i = e_i and e_i * e₀ = e_i
for i in range(8):
    MULT_TABLE[0, i, i] = 1.0
    MULT_TABLE[i, 0, i] = 1.0

# e_i * e_i = -e₀ for i=1..7
for i in range(1, 8):
    MULT_TABLE[i, i, 0] = -1.0

# Fano triples: e_a * e_b = e_c, e_b * e_a = -e_c
for (a, b, c) in FANO_TRIPLES:
    MULT_TABLE[a, b, c] = 1.0
    MULT_TABLE[b, a, c] = -1.0
    # cyclic: e_b * e_c = e_a, e_c * e_b = -e_a
    MULT_TABLE[b, c, a] = 1.0
    MULT_TABLE[c, b, a] = -1.0
    # cyclic: e_c * e_a = e_b, e_a * e_c = -e_b
    MULT_TABLE[c, a, b] = 1.0
    MULT_TABLE[a, c, b] = -1.0


def oct_mult(x, y):
    """Multiply two octonions x, y ∈ ℝ⁸. Returns z ∈ ℝ⁸."""
    return np.einsum('ijk,i,j->k', MULT_TABLE, x, y)


def oct_conj(x):
    """Octonion conjugate: x̄ = x₀ - x₁e₁ - ... - x₇e₇."""
    xbar = -x.copy()
    xbar[0] = x[0]
    return xbar


def oct_norm_sq(x):
    """Squared norm |x|² = x · x̄ = Σ xᵢ²."""
    return np.dot(x, x)


def oct_re(x):
    """Real part Re(x) = x₀."""
    return x[0]


def oct_basis(i):
    """Return basis element eᵢ as ℝ⁸ vector."""
    e = np.zeros(8)
    e[i] = 1.0
    return e


# ===================================================================
# J₃(𝕆): 3×3 HERMITIAN OCTONIONIC MATRICES
# ===================================================================
# A general element X ∈ J₃(𝕆) is:
#
#   X = | ξ₁    x₃*   x₂  |
#       | x₃    ξ₂    x₁* |
#       | x₂*   x₁    ξ₃  |
#
# where ξ₁, ξ₂, ξ₃ ∈ ℝ and x₁, x₂, x₃ ∈ 𝕆.
# Real dimension = 3 + 3×8 = 27.
#
# We represent X as a tuple (diag, off) where:
#   diag = (ξ₁, ξ₂, ξ₃) — 3 reals
#   off = (x₁, x₂, x₃)  — 3 octonions (ℝ⁸ each)
#
# But for linear algebra we use a flat ℝ²⁷ representation:
#   indices 0,1,2 → ξ₁, ξ₂, ξ₃
#   indices 3..10 → x₁ (8 components)
#   indices 11..18 → x₂ (8 components)
#   indices 19..26 → x₃ (8 components)


def j3_basis():
    """Return 27 basis elements of J₃(𝕆) as ℝ²⁷ vectors.

    Returns list of 27 basis vectors. Ordering:
      [0..2]   → diagonal: E_{11}, E_{22}, E_{33}
      [3..10]  → x₁ components (off-diag 23 slot): e₀..e₇
      [11..18] → x₂ components (off-diag 13 slot): e₀..e₇
      [19..26] → x₃ components (off-diag 12 slot): e₀..e₇
    """
    basis = []
    for i in range(27):
        v = np.zeros(27)
        v[i] = 1.0
        basis.append(v)
    return basis


def j3_to_components(v):
    """Unpack ℝ²⁷ → (ξ₁, ξ₂, ξ₃, x₁, x₂, x₃)."""
    xi1, xi2, xi3 = v[0], v[1], v[2]
    x1 = v[3:11].copy()
    x2 = v[11:19].copy()
    x3 = v[19:27].copy()
    return xi1, xi2, xi3, x1, x2, x3


def components_to_j3(xi1, xi2, xi3, x1, x2, x3):
    """Pack components → ℝ²⁷."""
    v = np.zeros(27)
    v[0], v[1], v[2] = xi1, xi2, xi3
    v[3:11] = x1
    v[11:19] = x2
    v[19:27] = x3
    return v


def jordan_product(u, v):
    """Jordan product X∘Y = ½(XY + YX) in J₃(𝕆).

    For X = (ξ, x) and Y = (η, y):
      (X∘Y)_{ii} = ξᵢηᵢ + Re(xⱼ · ȳⱼ) + Re(xₖ · ȳₖ)
                    where (i,j,k) are cyclic on {1,2,3} in the off-diagonal assignment.

    Off-diagonal: (X∘Y)_off_m = ½(...) — must use octonion product carefully.

    Convention for off-diag indices: x₁ is in (2,3) slot, x₂ in (1,3), x₃ in (1,2).
    """
    xi1, xi2, xi3, x1, x2, x3 = j3_to_components(u)
    eta1, eta2, eta3, y1, y2, y3 = j3_to_components(v)

    # Diagonal components of X∘Y:
    # (X∘Y)₁₁ = ξ₁η₁ + Re(x₃·ȳ₃) + Re(x₂·ȳ₂)
    # (X∘Y)₂₂ = ξ₂η₂ + Re(x₃·ȳ₃) + Re(x₁·ȳ₁)  — x₃ in (1,2), x₁ in (2,3)
    # (X∘Y)₃₃ = ξ₃η₃ + Re(x₂·ȳ₂) + Re(x₁·ȳ₁)
    #
    # Wait: more carefully. For 3×3 hermitian matrix product:
    # (XY)_{ii} = Σ_j X_{ij} Y_{ji} = Σ_j X_{ij} Y_{ij}*
    #
    # Let's use direct matrix multiplication approach.
    # X_{11}=ξ₁, X_{22}=ξ₂, X_{33}=ξ₃
    # X_{23}=x₁*, X_{32}=x₁, X_{13}=x₂, X_{31}=x₂*, X_{12}=x₃*, X_{21}=x₃

    # (XY)_{11} = X₁₁Y₁₁ + X₁₂Y₂₁ + X₁₃Y₃₁
    #           = ξ₁η₁ + x₃*·y₃ + x₂·y₂*
    # Re part for hermitian: Re(x₃*·y₃) + Re(x₂·y₂*)
    # For Jordan product (X∘Y)_{11} = ½((XY)₁₁ + (YX)₁₁)
    # By hermitianity of X∘Y, diagonal is real.

    # It's simpler to compute XY + YX directly using octonionic matrix multiply.
    # Let's define the full matrix multiply and then symmetrize.

    # XY matrix (3×3 octonionic):
    # (XY)_{ij} = Σ_k X_{ik} · Y_{kj}  (octonion product, non-associative!)

    # We build X and Y as 3×3 arrays of octonions
    X = _j3_to_matrix(xi1, xi2, xi3, x1, x2, x3)
    Y = _j3_to_matrix(eta1, eta2, eta3, y1, y2, y3)

    XY = _oct_mat_mult(X, Y)
    YX = _oct_mat_mult(Y, X)

    # Z = ½(XY + YX) — this should be hermitian
    Z = [[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            Z[i][j] = 0.5 * (XY[i][j] + YX[i][j])

    # Extract components
    z_xi1 = Z[0][0][0]  # real part of diagonal
    z_xi2 = Z[1][1][0]
    z_xi3 = Z[2][2][0]
    z_x1 = Z[2][1]  # (3,2) entry = x₁ component  (row 3, col 2 in 1-indexed)
    z_x2 = Z[0][2]  # (1,3) entry = x₂ component
    z_x3 = Z[1][0]  # (2,1) entry = x₃ component

    return components_to_j3(z_xi1, z_xi2, z_xi3, z_x1, z_x2, z_x3)


def _j3_to_matrix(xi1, xi2, xi3, x1, x2, x3):
    """Build 3×3 octonionic matrix from J₃(𝕆) components.

    M = | ξ₁·e₀   x₃*      x₂    |
        | x₃       ξ₂·e₀   x₁*   |
        | x₂*      x₁       ξ₃·e₀ |
    """
    M = [[None]*3 for _ in range(3)]
    M[0][0] = xi1 * oct_basis(0)
    M[1][1] = xi2 * oct_basis(0)
    M[2][2] = xi3 * oct_basis(0)
    M[0][1] = oct_conj(x3)
    M[1][0] = x3.copy()
    M[0][2] = x2.copy()
    M[2][0] = oct_conj(x2)
    M[1][2] = oct_conj(x1)
    M[2][1] = x1.copy()
    return M


def _oct_mat_mult(A, B):
    """Multiply two 3×3 octonionic matrices: (AB)_{ij} = Σ_k A_{ik}·B_{kj}."""
    C = [[np.zeros(8) for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            s = np.zeros(8)
            for k in range(3):
                s = s + oct_mult(A[i][k], B[k][j])
            C[i][j] = s
    return C


def jordan_product_matrix():
    """Compute Jordan product structure constants L_a: (E_a ∘ X)_b = (L_a)_{b,c} X_c.

    Returns: array of shape (27, 27, 27) where J[a,b,c] means
      (E_a ∘ E_b)_c = J[a,b,c], i.e., E_a ∘ E_b = Σ_c J[a,b,c] E_c
    """
    basis = j3_basis()
    n = 27
    J = np.zeros((n, n, n))
    for a in range(n):
        for b in range(n):
            product = jordan_product(basis[a], basis[b])
            J[a, b, :] = product
    return J


# ===================================================================
# CUBIC NORM ON J₃(𝕆)
# ===================================================================
def cubic_norm(v):
    """Cubic norm N(X) = det(X) for X ∈ J₃(𝕆).

    Using the hermitian matrix determinant formula:
      N(X) = ξ₁ξ₂ξ₃ - ξ₁|x₁|² - ξ₂|x₂|² - ξ₃|x₃|² + 2 Re(x₃*·(x₁*·x₂*))

    The cross-term uses the CONJUGATED off-diagonal entries (the actual upper-triangle
    entries of the hermitian matrix), NOT the raw x_i components.  This distinction
    matters because J₃(𝕆) is an exceptional (non-special) Jordan algebra.

    Note: octonion product is non-associative, so the parenthesization matters.
    """
    xi1, xi2, xi3, x1, x2, x3 = j3_to_components(v)
    # Conjugated entries = upper-triangle entries of the hermitian matrix
    x1s = oct_conj(x1)
    x2s = oct_conj(x2)
    x3s = oct_conj(x3)
    return (xi1 * xi2 * xi3
            - xi1 * oct_norm_sq(x1)
            - xi2 * oct_norm_sq(x2)
            - xi3 * oct_norm_sq(x3)
            + 2.0 * oct_re(oct_mult(x3s, oct_mult(x1s, x2s))))


def cubic_norm_tensor():
    """Compute the fully symmetric cubic form N_{abc} such that
    N(X) = N_{abc} X^a X^b X^c.

    This is the linearization of det(X).
    N_{abc} = (1/6)(N(eₐ+eᵦ+e_c) - N(eₐ+eᵦ) - N(eₐ+e_c) - N(eᵦ+e_c)
               + N(eₐ) + N(eᵦ) + N(e_c) - ... )

    More efficiently via polarization identity:
    6·N_{abc} = N(a+b+c) - N(a+b) - N(a+c) - N(b+c) + N(a) + N(b) + N(c)
    (for N homogeneous degree 3 with N(0)=0).
    """
    basis = j3_basis()
    n = 27
    T = np.zeros((n, n, n))

    # Precompute N on all basis elements and pairs
    N_single = np.array([cubic_norm(basis[a]) for a in range(n)])

    for a in range(n):
        for b in range(a, n):
            ab = basis[a] + basis[b]
            N_ab = cubic_norm(ab)
            for c in range(b, n):
                ac = basis[a] + basis[c]
                bc = basis[b] + basis[c]
                abc = basis[a] + basis[b] + basis[c]
                N_ac = cubic_norm(ac)
                N_bc = cubic_norm(bc)
                N_abc = cubic_norm(abc)

                val = (N_abc - N_ab - N_ac - N_bc
                       + N_single[a] + N_single[b] + N_single[c]) / 6.0

                # Symmetrize over all permutations of (a,b,c)
                for p, q, r in [(a,b,c),(a,c,b),(b,a,c),(b,c,a),(c,a,b),(c,b,a)]:
                    T[p, q, r] = val

    return T


# ===================================================================
# F₄ CONSTRUCTION: Der(J₃(𝕆))
# ===================================================================
def build_f4_generators():
    """Construct the 52 generators of f₄ = Der(J₃(𝕆)).

    A derivation D of J₃(𝕆) satisfies: D(X∘Y) = D(X)∘Y + X∘D(Y)
    for all X, Y ∈ J₃(𝕆).

    In the ℝ²⁷ basis, D is a 27×27 real matrix satisfying:
      D^c_a J_{bc}^d + D^c_b J_{ac}^d = D^d_c J_{ab}^c
    equivalently:
      Σ_c D_{ac} J_{cb}^d + Σ_c D_{bc} J_{ac}^d - Σ_c J_{ab}^c D_{cd} = 0 for all a,b,d

    This gives a linear system for the 27×27 = 729 entries of D.

    Returns: list of 52 orthonormalized 27×27 real matrices.
    """
    print("  Computing Jordan product structure constants (27³)...")
    J = jordan_product_matrix()  # (27, 27, 27)

    print("  Building derivation constraint matrix...")
    # The derivation condition (D acts by left multiplication: (DX)_a = Σ_c D[a,c] X_c):
    #   D(X∘Y) = D(X)∘Y + X∘D(Y)
    # Expanding: Σ_c D[d,c] J[a,b,c] = Σ_c D[c,a] J[c,b,d] + Σ_c D[c,b] J[a,c,d]
    # Rearranging:
    #   Σ_c D[c,a] J[c,b,d] + Σ_c D[c,b] J[a,c,d] - Σ_c D[d,c] J[a,b,c] = 0
    #
    # D has 729 unknowns: D_{ij} for i,j = 0..26
    # Variable index for D_{i,j} = 27*i + j
    #
    # For each (a,b,d) triple, we get one constraint row.
    # Some are redundant due to J symmetry.

    dim = 27
    n_vars = dim * dim  # 729

    # Build constraint matrix more efficiently using symmetry of J
    # Only use (a,b,d) with a <= b to reduce row count (due to J symmetry)
    rows = []
    for a in range(dim):
        for b in range(a, dim):
            for d in range(dim):
                row = np.zeros(n_vars)
                for c in range(dim):
                    # Term 1: D[c,a] * J[c,b,d]
                    if abs(J[c, b, d]) > 1e-15:
                        row[c * dim + a] += J[c, b, d]
                    # Term 2: D[c,b] * J[a,c,d]
                    if abs(J[a, c, d]) > 1e-15:
                        row[c * dim + b] += J[a, c, d]
                    # Term 3: -J[a,b,c] * D[d,c]
                    if abs(J[a, b, c]) > 1e-15:
                        row[d * dim + c] -= J[a, b, c]
                if np.any(np.abs(row) > 1e-15):
                    rows.append(row)

    M = np.array(rows)
    print(f"  Constraint matrix: {M.shape[0]} × {M.shape[1]}")

    # SVD to find null space
    print("  Computing SVD...")
    _, S, Vt = np.linalg.svd(M, full_matrices=True)
    rank = np.sum(S > 1e-8)
    null_dim = n_vars - rank
    print(f"  SVD: rank={rank}, null space dim={null_dim} (expected: 52)")

    if null_dim != 52:
        print(f"  WARNING: null space dimension = {null_dim} (expected 52 for f₄)")

    # Null space vectors → 27×27 matrices
    null_vecs = Vt[rank:]  # shape (null_dim, 729)
    generators = []
    for v in null_vecs:
        D = v.reshape(dim, dim)
        generators.append(D)

    # Orthonormalize w.r.t. Frobenius inner product
    generators = _orthonormalize_frobenius(generators)

    print(f"  Constructed {len(generators)} orthonormal f₄ generators (27×27 real)")
    return generators


# ===================================================================
# E₆ CONSTRUCTION: Cubic norm preservation
# ===================================================================
def build_e6_generators():
    """Construct the 78 generators of e₆ from cubic norm preservation.

    e₆ = {L ∈ gl(27,ℝ) : L preserves the cubic norm N AND is traceless}

    The cubic norm preservation condition:
      N(L·X, X, X) = 0 for all X
    linearized: Σ_a L_{da} N_{abc} + Σ_a L_{da} N_{abc} + ... = 0
    i.e., L^d_a N_{dbc} + L^d_b N_{adc} + L^d_c N_{abd} = 0 for all (a,b,c)

    And Tr(L) = 0.

    Returns: list of 78 orthonormalized 27×27 real matrices.
    """
    print("  Computing cubic norm tensor N_{abc} (27³)...")
    N = cubic_norm_tensor()  # (27, 27, 27)

    print("  Building cubic preservation constraint matrix...")
    dim = 27
    n_vars = dim * dim  # 729

    rows = []
    # For each symmetric triple (a,b,c) with a <= b <= c:
    for a in range(dim):
        for b in range(a, dim):
            for c in range(b, dim):
                row = np.zeros(n_vars)
                for d in range(dim):
                    # L_{da} N[d,b,c]
                    if abs(N[d, b, c]) > 1e-15:
                        row[d * dim + a] += N[d, b, c]
                    # L_{db} N[a,d,c]
                    if abs(N[a, d, c]) > 1e-15:
                        row[d * dim + b] += N[a, d, c]
                    # L_{dc} N[a,b,d]
                    if abs(N[a, b, d]) > 1e-15:
                        row[d * dim + c] += N[a, b, d]
                if np.any(np.abs(row) > 1e-15):
                    rows.append(row)

    # Add tracelessness: Σ_i L_{ii} = 0
    trace_row = np.zeros(n_vars)
    for i in range(dim):
        trace_row[i * dim + i] = 1.0
    rows.append(trace_row)

    M = np.array(rows)
    print(f"  Constraint matrix: {M.shape[0]} × {M.shape[1]}")

    # SVD
    print("  Computing SVD...")
    _, S, Vt = np.linalg.svd(M, full_matrices=True)
    rank = np.sum(S > 1e-8)
    null_dim = n_vars - rank
    print(f"  SVD: rank={rank}, null space dim={null_dim} (expected: 78)")

    if null_dim != 78:
        print(f"  WARNING: null space dimension = {null_dim} (expected 78 for e₆)")

    # Null space → 27×27 matrices
    null_vecs = Vt[rank:]
    generators = []
    for v in null_vecs:
        L = v.reshape(dim, dim)
        generators.append(L)

    # Orthonormalize
    generators = _orthonormalize_frobenius(generators)

    print(f"  Constructed {len(generators)} orthonormal e₆ generators (27×27 real)")
    return generators


def _orthonormalize_frobenius(matrices):
    """Gram-Schmidt orthonormalization with Frobenius inner product."""
    ortho = []
    for g in matrices:
        v = g.copy()
        for u in ortho:
            proj = np.sum(v * u) / np.sum(u * u)
            v = v - proj * u
        norm = np.sqrt(np.sum(v * v))
        if norm > 1e-12:
            ortho.append(v / norm)
    return ortho


# ===================================================================
# VERIFICATION HELPERS
# ===================================================================
def verify_derivation(D, J, tol=1e-8):
    """Check that D is a derivation of Jordan algebra with structure constants J.

    Convention: D acts as (D·v)_a = Σ_c D[a,c] v_c.
    Derivation condition: Σ_c D[c,a] J[c,b,d] + Σ_c D[c,b] J[a,c,d] - Σ_c D[d,c] J[a,b,c] = 0

    Returns max |constraint| over all (a,b,d).
    """
    dim = D.shape[0]
    max_err = 0.0
    for a in range(dim):
        for b in range(a, dim):
            for d in range(dim):
                val = 0.0
                for c in range(dim):
                    val += D[c, a] * J[c, b, d] + D[c, b] * J[a, c, d] - D[d, c] * J[a, b, c]
                max_err = max(max_err, abs(val))
    return max_err


def verify_cubic_preservation(L, N, tol=1e-8):
    """Check that L preserves cubic norm: L^d_a N_{dbc} + cyclic = 0.

    Returns max error over all (a,b,c).
    """
    dim = L.shape[0]
    max_err = 0.0
    for a in range(dim):
        for b in range(a, dim):
            for c in range(b, dim):
                val = 0.0
                for d in range(dim):
                    val += L[d, a] * N[d, b, c]
                    val += L[d, b] * N[a, d, c]
                    val += L[d, c] * N[a, b, d]
                max_err = max(max_err, abs(val))
    return max_err
