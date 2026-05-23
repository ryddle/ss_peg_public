"""
Holonomy Zeta with SU(2) Operators on Edges
=============================================

Extends the toy model: instead of scalar weights on edges, we place
SU(2) matrices (non-commutative holonomies).

Key idea from SS-PEG:
  - Edges carry unitary operators U_e ∈ SU(2)
  - Holonomy around a cycle: W = U₃·U₂·U₁
  - Trace formula: Tr(T^n) involves Tr(holonomy) for each closed walk
  - Non-trivial W ≠ I means physical content (gauge flux)

This script:
  1. Builds a 3-node graph with SU(2) operators on edges
  2. Computes the "block" transfer matrix T (6×6 for 3 nodes × 2 internal)
  3. Verifies the trace formula Tr(T^n) = Σ_{closed walks of length n} Tr(holonomy)
  4. Constructs the generalized zeta and finds poles
  5. Explores how poles move as holonomy varies from trivial (W=I) to maximal (W=random)
  6. Connection to Killing eigenvalues: cyclic edges → compact holonomy

Usage:
    python holonomy_zeta_su2.py [--n-samples 50] [--save]
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import argparse


# ============================================================
# SU(2) utilities
# ============================================================

def su2_generator(axis):
    """Pauli matrices (basis of su(2) Lie algebra)."""
    sigma = {
        'x': np.array([[0, 1], [1, 0]], dtype=np.complex128),
        'y': np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
        'z': np.array([[1, 0], [0, -1]], dtype=np.complex128),
    }
    return sigma[axis]


def random_su2(rng=None):
    """Generate a random SU(2) matrix via exponential map."""
    if rng is None:
        rng = np.random.default_rng()
    # Random direction on S²
    theta = rng.uniform(0, np.pi)
    phi = rng.uniform(0, 2 * np.pi)
    angle = rng.uniform(0, 2 * np.pi)

    n = np.array([np.sin(theta)*np.cos(phi),
                  np.sin(theta)*np.sin(phi),
                  np.cos(theta)])

    # U = exp(i·angle/2·(n·σ))
    sigma = [su2_generator(a) for a in ['x', 'y', 'z']]
    H = sum(n_i * s for n_i, s in zip(n, sigma))
    U = np.cos(angle/2) * np.eye(2) + 1j * np.sin(angle/2) * H
    return U


def su2_from_angle(angle, axis='z'):
    """SU(2) rotation by angle around given axis."""
    sigma = su2_generator(axis)
    return np.cos(angle/2) * np.eye(2, dtype=np.complex128) + \
           1j * np.sin(angle/2) * sigma


def holonomy_class(U):
    """
    Extract the conjugacy class of U ∈ SU(2).
    Returns the holonomy angle θ ∈ [0, π] such that
    eigenvalues are e^{±iθ}.
    """
    tr = np.trace(U)
    cos_theta = np.clip(tr.real / 2, -1, 1)
    return np.arccos(cos_theta)


# ============================================================
# Block transfer matrix with SU(2) on edges
# ============================================================

def build_block_transfer(U_edges, alpha=0.0):
    """
    Build the 6×6 block transfer matrix for a 3-node graph
    with SU(2) operators on directed edges.

    U_edges = [U_AB, U_BC, U_CA] — three 2×2 SU(2) matrices

    Graph: A→B (U_AB), B→C (U_BC), C→A (U_CA)

    The block matrix acts on the space C² ⊗ C³ (internal × node):

    T = [[α·I₂,    0,      U_CA],
         [U_AB,    α·I₂,   0   ],
         [0,       U_BC,   α·I₂]]

    where I₂ = 2×2 identity (self-loop with trivial holonomy).
    """
    I2 = np.eye(2, dtype=np.complex128)
    Z2 = np.zeros((2, 2), dtype=np.complex128)

    U_AB, U_BC, U_CA = U_edges

    T = np.block([
        [alpha * I2, Z2,          U_CA],
        [U_AB,       alpha * I2,  Z2],
        [Z2,         U_BC,        alpha * I2]
    ])
    return T


# ============================================================
# Trace formula for block matrix
# ============================================================

def block_trace_formula(U_edges, alpha=0.0, max_n=15):
    """
    Compute Tr(T^n) for the block transfer matrix and decompose
    into holonomy contributions.

    For n=3: Tr(T³) includes Tr(U_CA·U_BC·U_AB) + cyclic permutations
             = Tr(W) + Tr(W') + Tr(W'') where W is the cycle holonomy.

    For a directed cycle, all cyclic permutations of the holonomy product
    have the same trace (cyclic property of trace), so:
      cycle contribution at n=3 = 3·Re[Tr(U_CA·U_BC·U_AB)]
    """
    T = build_block_transfer(U_edges, alpha)
    eigenvalues = np.linalg.eigvals(T)

    # Full holonomy around the cycle
    U_AB, U_BC, U_CA = U_edges
    W = U_CA @ U_BC @ U_AB  # holonomy A→B→C→A

    results = {
        'n': [], 'tr_block': [], 'tr_spectral': [],
        'holonomy_trace': np.trace(W),
        'holonomy_angle': holonomy_class(W),
        'eigenvalues': eigenvalues,
        'W': W
    }

    T_power = np.eye(6, dtype=np.complex128)
    for n in range(1, max_n + 1):
        T_power = T_power @ T
        tr_block = np.trace(T_power)
        tr_spectral = np.sum(eigenvalues ** n)

        results['n'].append(n)
        results['tr_block'].append(tr_block)
        results['tr_spectral'].append(tr_spectral)

    return results


# ============================================================
# Generalized Ihara zeta for the block graph
# ============================================================

def block_ihara_zeta_reciprocal(u, T_block):
    """
    Generalized zeta: Z(u)^{-1} = det(I - u·T_block)

    This is the most general form: the Ihara zeta of a directed graph
    with matrix-valued edge weights.
    """
    n = T_block.shape[0]
    I = np.eye(n, dtype=np.complex128)
    return np.linalg.det(I - u * T_block)


def find_block_poles(T_block):
    """
    Poles of Z(u) = 1/det(I - uT) are at u = 1/λ_k
    where λ_k are eigenvalues of T.
    """
    eigs = np.linalg.eigvals(T_block)
    # Filter out near-zero eigenvalues
    nonzero = eigs[np.abs(eigs) > 1e-12]
    return 1.0 / nonzero


# ============================================================
# Scan: poles vs holonomy angle
# ============================================================

def scan_holonomy_angle(n_points=100, alpha=0.0):
    """
    Fix U_AB = U_BC = I, vary U_CA = exp(iθσ_z/2) from θ=0 to 2π.
    Track how the poles of Z(u) move.

    At θ=0: W = I (trivial holonomy, no flux)
    At θ=π: W = iσ_z (maximal holonomy class)
    """
    angles = np.linspace(0, 2*np.pi, n_points)
    all_poles = []

    I2 = np.eye(2, dtype=np.complex128)

    for theta in angles:
        U_CA = su2_from_angle(theta, axis='z')
        U_edges = [I2, I2, U_CA]
        T = build_block_transfer(U_edges, alpha)
        poles = find_block_poles(T)
        all_poles.append(poles)

    return angles, all_poles


def scan_random_holonomies(n_samples=200, alpha=0.0, seed=42):
    """
    Random SU(2) on all three edges. Collect pole statistics.
    """
    rng = np.random.default_rng(seed)
    all_poles = []
    all_holonomy_angles = []

    for _ in range(n_samples):
        U_edges = [random_su2(rng) for _ in range(3)]
        T = build_block_transfer(U_edges, alpha)
        poles = find_block_poles(T)
        W = U_edges[2] @ U_edges[1] @ U_edges[0]
        all_poles.append(poles)
        all_holonomy_angles.append(holonomy_class(W))

    return all_poles, all_holonomy_angles


# ============================================================
# Visualization
# ============================================================

def plot_all(n_samples=50, save=False):
    fig = plt.figure(figsize=(18, 16))
    fig.suptitle(
        'SS-PEG: Non-Abelian Holonomy Zeta (SU(2) on 3-Cycle)',
        fontsize=16, fontweight='bold'
    )
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.35)

    # --- Panel 1: Trivial holonomy (W=I) trace formula ---
    ax1 = fig.add_subplot(gs[0, 0])
    I2 = np.eye(2, dtype=np.complex128)
    res_trivial = block_trace_formula([I2, I2, I2], alpha=0.0, max_n=15)
    ax1.plot(res_trivial['n'], [t.real for t in res_trivial['tr_block']],
            'bo-', ms=4, label='Tr(T^n) block')
    ax1.plot(res_trivial['n'], [t.real for t in res_trivial['tr_spectral']],
            'rx', ms=8, label='Σ λ_k^n')
    ax1.set_xlabel('n')
    ax1.set_ylabel('Tr(T^n)')
    ax1.set_title(f'Trivial holonomy W=I\nTr(W)={res_trivial["holonomy_trace"]:.3f}')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # --- Panel 2: Non-trivial holonomy (θ=π/3) ---
    ax2 = fig.add_subplot(gs[0, 1])
    U_CA = su2_from_angle(np.pi/3, 'z')
    res_nontrivial = block_trace_formula([I2, I2, U_CA], alpha=0.0, max_n=15)
    ax2.plot(res_nontrivial['n'], [t.real for t in res_nontrivial['tr_block']],
            'bo-', ms=4, label='Tr(T^n) block')
    ax2.plot(res_nontrivial['n'], [t.real for t in res_nontrivial['tr_spectral']],
            'rx', ms=8, label='Σ λ_k^n')
    theta_hol = res_nontrivial['holonomy_angle']
    ax2.set_xlabel('n')
    ax2.set_ylabel('Tr(T^n)')
    ax2.set_title(f'θ_holonomy = π/3\nTr(W)={res_nontrivial["holonomy_trace"]:.3f}, '
                  f'θ={theta_hol:.3f}')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    # --- Panel 3: Random SU(2) holonomy ---
    ax3 = fig.add_subplot(gs[0, 2])
    rng = np.random.default_rng(7)
    U_rand = [random_su2(rng) for _ in range(3)]
    res_random = block_trace_formula(U_rand, alpha=0.0, max_n=15)
    ax3.plot(res_random['n'], [t.real for t in res_random['tr_block']],
            'bo-', ms=4, label='Tr(T^n) block')
    ax3.plot(res_random['n'], [t.real for t in res_random['tr_spectral']],
            'rx', ms=8, label='Σ λ_k^n')
    theta_rand = res_random['holonomy_angle']
    ax3.set_xlabel('n')
    ax3.set_ylabel('Tr(T^n)')
    ax3.set_title(f'Random SU(2)\nTr(W)={res_random["holonomy_trace"]:.3f}, '
                  f'θ={theta_rand:.3f}')
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)

    # --- Panel 4: Eigenvalues comparison (3 cases) ---
    ax4 = fig.add_subplot(gs[1, 0])
    for label, res, color in [('W=I', res_trivial, 'blue'),
                               ('θ=π/3', res_nontrivial, 'green'),
                               ('random', res_random, 'red')]:
        eigs = res['eigenvalues']
        ax4.scatter(eigs.real, eigs.imag, s=50, c=color, label=label,
                   alpha=0.7, edgecolors='black', linewidths=0.5)
    theta_circle = np.linspace(0, 2*np.pi, 100)
    ax4.plot(np.cos(theta_circle), np.sin(theta_circle), 'k--', alpha=0.2)
    ax4.set_xlabel('Re(λ)')
    ax4.set_ylabel('Im(λ)')
    ax4.set_title('Eigenvalues of T_block')
    ax4.set_aspect('equal')
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)

    # --- Panel 5: Poles comparison ---
    ax5 = fig.add_subplot(gs[1, 1])
    for label, edges, color in [('W=I', [I2, I2, I2], 'blue'),
                                 ('θ=π/3', [I2, I2, U_CA], 'green'),
                                 ('random', U_rand, 'red')]:
        T = build_block_transfer(edges, alpha=0.0)
        poles = find_block_poles(T)
        ax5.scatter(poles.real, poles.imag, s=60, c=color, marker='x',
                   linewidths=2, label=label, alpha=0.8)
    ax5.plot(np.cos(theta_circle), np.sin(theta_circle), 'k--', alpha=0.2, label='|u|=1')
    ax5.set_xlabel('Re(u)')
    ax5.set_ylabel('Im(u)')
    ax5.set_title('Poles of generalized Z(u)')
    ax5.set_aspect('equal')
    ax5.legend(fontsize=8)
    ax5.grid(True, alpha=0.3)

    # --- Panel 6: Pole migration as θ varies ---
    ax6 = fig.add_subplot(gs[1, 2])
    angles, all_poles = scan_holonomy_angle(n_points=150, alpha=0.0)
    for i, (ang, poles) in enumerate(zip(angles, all_poles)):
        color = plt.cm.viridis(ang / (2*np.pi))
        ax6.scatter(poles.real, poles.imag, s=3, c=[color]*len(poles), alpha=0.5)
    ax6.plot(np.cos(theta_circle), np.sin(theta_circle), 'k--', alpha=0.2)
    ax6.set_xlabel('Re(u)')
    ax6.set_ylabel('Im(u)')
    ax6.set_title('Pole migration: θ = 0 → 2π\n(U_CA = exp(iθσ_z/2))')
    ax6.set_aspect('equal')
    ax6.grid(True, alpha=0.3)

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(0, 2*np.pi))
    sm.set_array([])
    plt.colorbar(sm, ax=ax6, label='θ')

    # --- Panel 7: Pole |u| distribution for random SU(2) ---
    ax7 = fig.add_subplot(gs[2, 0])
    all_random_poles, all_hol_angles = scan_random_holonomies(
        n_samples=n_samples, alpha=0.0, seed=42
    )
    all_pole_radii = np.concatenate([np.abs(p) for p in all_random_poles])
    # Adaptive bins: if all radii are nearly identical, use fewer bins
    radii_range = all_pole_radii.max() - all_pole_radii.min()
    n_bins = max(5, min(60, int(radii_range / max(1e-12, radii_range) * 60)))
    if radii_range < 1e-8:
        n_bins = 5
    ax7.hist(all_pole_radii, bins=n_bins, color='steelblue', alpha=0.7, edgecolor='black',
            linewidth=0.5, density=True)
    ax7.axvline(1.0, color='red', linestyle='--', label='|u|=1')
    ax7.set_xlabel('|u| (pole radius)')
    ax7.set_ylabel('Density')
    ax7.set_title(f'Pole radius distribution\n({n_samples} random SU(2) configs)')
    ax7.legend(fontsize=8)
    ax7.grid(True, alpha=0.3)

    # --- Panel 8: Holonomy angle vs min pole radius ---
    ax8 = fig.add_subplot(gs[2, 1])
    min_radii = [np.min(np.abs(p)) for p in all_random_poles]
    ax8.scatter(all_hol_angles, min_radii, s=15, alpha=0.6, c='steelblue',
               edgecolors='none')
    ax8.set_xlabel('Holonomy angle θ_W')
    ax8.set_ylabel('min |u_pole|')
    ax8.set_title('Smallest pole radius vs holonomy')
    ax8.grid(True, alpha=0.3)
    ax8.axhline(1.0, color='red', linestyle='--', alpha=0.5, label='|u|=1')
    ax8.legend(fontsize=8)

    # --- Panel 9: All poles in complex plane (random samples) ---
    ax9 = fig.add_subplot(gs[2, 2])
    for poles, hol_angle in zip(all_random_poles, all_hol_angles):
        color = plt.cm.RdYlBu(hol_angle / np.pi)
        ax9.scatter(poles.real, poles.imag, s=2, c=[color]*len(poles), alpha=0.3)
    ax9.plot(np.cos(theta_circle), np.sin(theta_circle), 'k-', alpha=0.3)
    ax9.set_xlabel('Re(u)')
    ax9.set_ylabel('Im(u)')
    ax9.set_title(f'All poles ({n_samples} configs)\nColored by θ_W')
    ax9.set_aspect('equal')
    ax9.grid(True, alpha=0.3)
    sm2 = plt.cm.ScalarMappable(cmap='RdYlBu', norm=plt.Normalize(0, np.pi))
    sm2.set_array([])
    plt.colorbar(sm2, ax=ax9, label='θ_W')

    plt.tight_layout()

    if save:
        out_path = 'holonomy_zeta_su2_results.png'
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        print(f"Saved to {out_path}")

    plt.show()


# ============================================================
# Console output
# ============================================================

def print_verification():
    print("=" * 70)
    print("NON-ABELIAN HOLONOMY ZETA — SU(2) ON 3-CYCLE")
    print("=" * 70)

    I2 = np.eye(2, dtype=np.complex128)

    configs = [
        ("Trivial (W=I)", [I2, I2, I2]),
        ("θ=π/3 on CA", [I2, I2, su2_from_angle(np.pi/3, 'z')]),
        ("θ=π on CA", [I2, I2, su2_from_angle(np.pi, 'z')]),
        ("Random SU(2)", [random_su2(np.random.default_rng(7)) for _ in range(3)]),
    ]

    for name, U_edges in configs:
        print(f"\n--- {name} ---")
        W = U_edges[2] @ U_edges[1] @ U_edges[0]
        theta_W = holonomy_class(W)
        print(f"  Holonomy W:\n{W}")
        print(f"  Tr(W) = {np.trace(W):.6f}")
        print(f"  θ_W = {theta_W:.6f} ({theta_W/np.pi:.4f}π)")
        print(f"  det(W) = {np.linalg.det(W):.6f} (should be 1)")

        T = build_block_transfer(U_edges, alpha=0.0)
        eigs = np.linalg.eigvals(T)
        poles = find_block_poles(T)

        print(f"  T eigenvalues: {np.sort_complex(eigs)}")
        print(f"  Zeta poles:    {np.sort_complex(poles)}")
        print(f"  |poles|:       {np.sort(np.abs(poles))}")

        # Verify Tr(T³) contains Tr(W)
        T3 = T @ T @ T
        tr_T3 = np.trace(T3)
        print(f"  Tr(T³) = {tr_T3:.6f}")
        print(f"  3·Tr(W) = {3*np.trace(W):.6f}")
        print(f"  [Tr(T³) should contain 3·Tr(W) as cycle contribution]")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Holonomy Zeta with SU(2) on 3-Node Cycle'
    )
    parser.add_argument('--n-samples', type=int, default=50,
                        help='Random SU(2) samples for statistics (default: 50)')
    parser.add_argument('--save', action='store_true',
                        help='Save plots to file')
    args = parser.parse_args()

    print_verification()
    plot_all(n_samples=args.n_samples, save=args.save)
