"""
Landscape Relief Map — topographic visualization of the algebraic optimization landscape.

Generates a 2D PCA slice through the 128-seed parameter space, evaluates
T_total on a dense grid using GPU, and renders as a shaded relief map with
converged solutions marked by category and simplicity band.

Strategy:
  1. Load 128 converged generator sets (128, 52, 27, 27)
  2. Extract real antisymmetric parameters → 128 points in R^{18252}
  3. PCA → top-2 principal components define the "landscape plane"
  4. Dense grid on this plane → reconstruct generators → evaluate T_total on GPU
  5. Render as hillshaded topographic map with contour lines

Usage:
  python _landscape_relief_map.py
  python _landscape_relief_map.py --grid-size 150 --batch-size 50 --padding 0.4
"""

import sys, os, argparse, json, time
import numpy as np
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, Normalize, LogNorm
from matplotlib import cm, patheffects

sys.path.insert(0, os.path.dirname(__file__))
from shared_utils_gpu import (
    batched_killing_metric,
    batched_non_commutativity_tension,
    batched_norm_tension,
    params_to_antisymmetric,
    antisymmetric_to_params,
    DTYPE, FDTYPE, DEVICE,
)
from test_closure_driven import batched_closure_tension

# Constants matching the optimizer
N_GEN = 52
DIM = 27
N_PARAMS_PER_GEN = DIM * (DIM - 1) // 2   # 351
N_PARAMS_TOTAL = N_GEN * N_PARAMS_PER_GEN  # 18252
ALPHA = 5.0
BETA = 1.0
GAMMA_FINAL = 200.0  # annealing-stage gamma


# ============================================================
# DATA LOADING + PCA
# ============================================================

def load_and_flatten(npz_path):
    """Load generators, convert to flat real parameter vectors.

    Returns:
        flat_params: (N, 18252) real array
        results: list of per-seed metric dicts
        raw_gens: (N, 52, 27, 27) complex array
    """
    data = np.load(npz_path, allow_pickle=True)
    results = json.loads(str(data['results']))
    raw_gens = data['generators']  # (N, 52, 27, 27) complex128
    N = raw_gens.shape[0]

    print(f"  Loaded {N} seeds, generators shape: {raw_gens.shape}")

    # Convert to flat params using torch (consistent with optimizer parameterization)
    flat_params = np.zeros((N, N_PARAMS_TOTAL), dtype=np.float64)
    batch = 32
    for start in range(0, N, batch):
        end = min(start + batch, N)
        gens_t = torch.tensor(raw_gens[start:end], dtype=DTYPE, device='cpu')
        params_t = antisymmetric_to_params(gens_t)  # (B, 52, 351) real
        flat_params[start:end] = params_t.numpy().reshape(end - start, -1)

    return flat_params, results, raw_gens


def compute_pca(flat_params, n_components=10):
    """PCA via SVD on centered data. Returns mean, PCs, projections, variance ratios."""
    N, P = flat_params.shape
    mean = flat_params.mean(axis=0)
    centered = flat_params - mean

    # Economy SVD (N << P → U is N×N)
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)

    pcs = Vt[:n_components]  # (n_components, P)
    proj = centered @ pcs.T  # (N, n_components)
    explained = (S ** 2) / (S ** 2).sum()

    print(f"  PCA: top-2 explain {explained[0]:.1%} + {explained[1]:.1%} = {explained[:2].sum():.1%}")
    print(f"  top-10: {explained[:10].sum():.1%}")

    return mean, pcs, proj, explained


# ============================================================
# GPU GRID EVALUATION
# ============================================================

def evaluate_grid_gpu(mean, pc1, pc2, grid_x, grid_y,
                      batch_size=25, gamma=GAMMA_FINAL):
    """Evaluate T_total at each grid point on GPU.

    Returns: log10_T (ny, nx) array, components dict
    """
    nx, ny = len(grid_x), len(grid_y)
    xx, yy = np.meshgrid(grid_x, grid_y)  # (ny, nx)
    points = np.column_stack([xx.ravel(), yy.ravel()])  # (ny*nx, 2)
    n_points = len(points)

    device = DEVICE
    K_target = -1.0 * torch.eye(N_GEN, dtype=FDTYPE, device=device)

    T_total_arr = np.zeros(n_points, dtype=np.float64)
    T_cl_arr = np.zeros(n_points, dtype=np.float64)
    simplicity_arr = np.zeros(n_points, dtype=np.float64)

    print(f"  Evaluating {nx}×{ny} = {n_points} grid points (batch_size={batch_size})...")
    t0 = time.time()

    for start in range(0, n_points, batch_size):
        end = min(start + batch_size, n_points)
        batch_pts = points[start:end]  # (B, 2)
        B = len(batch_pts)

        # Reconstruct flat params: mean + x*pc1 + y*pc2
        flat = (mean[None, :] +
                batch_pts[:, 0:1] * pc1[None, :] +
                batch_pts[:, 1:2] * pc2[None, :])

        # Reshape to (B, n_gen, n_params_per_gen) and send to GPU
        params_t = torch.tensor(
            flat.reshape(B, N_GEN, N_PARAMS_PER_GEN),
            dtype=FDTYPE, device=device
        )

        with torch.no_grad():
            gens = params_to_antisymmetric(params_t, N_GEN, dim=DIM)

            T_nc = batched_non_commutativity_tension(gens)
            K = batched_killing_metric(gens)
            T_k = ((K - K_target.unsqueeze(0)) ** 2).sum(dim=(-2, -1))
            T_norm = batched_norm_tension(gens)
            T_cl = batched_closure_tension(gens)

            T_total = T_nc + ALPHA * T_k + BETA * T_norm + gamma * T_cl

            # Simplicity
            K_diag = K.diagonal(dim1=-2, dim2=-1)
            K_mean = K_diag.mean(dim=-1)
            K_std = K_diag.std(dim=-1)
            S = K_std / K_mean.abs().clamp(min=1e-10)

        T_total_arr[start:end] = T_total.cpu().numpy()
        T_cl_arr[start:end] = T_cl.cpu().numpy()
        simplicity_arr[start:end] = S.cpu().numpy()

        # Clean up GPU memory
        del gens, T_nc, K, T_k, T_norm, T_cl, T_total, params_t
        torch.cuda.empty_cache()

        # Progress
        batch_num = start // batch_size
        if batch_num % 20 == 0:
            elapsed = time.time() - t0
            pct = end / n_points * 100
            print(f"    {pct:5.1f}%  ({end}/{n_points})  {elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"  Grid complete: {elapsed:.1f}s ({n_points / elapsed:.0f} pts/s)")

    return (T_total_arr.reshape(ny, nx),
            T_cl_arr.reshape(ny, nx),
            simplicity_arr.reshape(ny, nx))


# ============================================================
# FIGURE: RELIEF MAP
# ============================================================

def plot_relief_map(grid_x, grid_y, log_T, T_cl_grid, simp_grid,
                    proj_2d, results, explained,
                    save_path, title="Algebraic Landscape"):
    """Render the topographic relief map with seed markers."""

    # -- Dark theme --
    bg = '#0d1117'
    fg = '#c9d1d9'
    grid_color = '#21262d'
    plt.rcParams.update({
        'figure.facecolor': bg, 'axes.facecolor': bg,
        'axes.edgecolor': fg, 'axes.labelcolor': fg,
        'xtick.color': fg, 'ytick.color': fg,
        'text.color': fg, 'legend.facecolor': '#161b22',
        'legend.edgecolor': '#30363d',
        'font.size': 11, 'font.family': 'sans-serif',
    })

    fig = plt.figure(figsize=(16, 12))

    # Main relief map
    ax_main = fig.add_axes([0.08, 0.10, 0.62, 0.82])
    # PCA variance inset
    ax_pca = fig.add_axes([0.74, 0.68, 0.22, 0.22])
    # Band legend inset
    ax_leg = fig.add_axes([0.74, 0.10, 0.22, 0.50])
    ax_leg.axis('off')

    # ------- MAIN PANEL: Hillshaded relief -------
    ls = LightSource(azdeg=315, altdeg=35)
    extent = [grid_x[0], grid_x[-1], grid_y[0], grid_y[-1]]

    # Clip extreme values for better contrast
    Z = log_T.copy()
    vmin, vmax = np.percentile(Z[np.isfinite(Z)], [2, 98])
    Z = np.clip(Z, vmin, vmax)

    # Hillshaded terrain
    rgb = ls.shade(Z, cmap=plt.cm.terrain, vert_exag=0.3,
                   blend_mode='soft', vmin=vmin, vmax=vmax)
    ax_main.imshow(rgb, extent=extent, origin='lower', aspect='auto',
                   interpolation='bilinear')

    # Contour lines (sparse)
    levels = np.linspace(vmin, vmax, 12)
    cs = ax_main.contour(grid_x, grid_y, Z, levels=levels,
                         colors='white', linewidths=0.3, alpha=0.35)

    # ------- SEED MARKERS -------
    # Band colors
    band_colors = {
        'B1': '#e63946',   # red
        'B2': '#f4a261',   # orange
        'B3': '#2a9d8f',   # teal
        'B4': '#457b9d',   # steel blue
        'B5': '#a855f7',   # purple
        'other': '#6b7280' # gray
    }
    band_edges = [0.0, 0.11, 0.17, 0.22, 0.27, 0.35]
    band_names = ['B1', 'B2', 'B3', 'B4', 'B5']

    def get_band(s):
        for i, name in enumerate(band_names):
            if band_edges[i] <= s < band_edges[i + 1]:
                return name
        return 'other'

    # Category markers
    cat_markers = {
        'CLOSED_NON_SIMPLE': 'o',
        'PARTIAL': '^',
        'NOT_CLOSED': 'x',
    }

    # Plot each seed
    for i, r in enumerate(results):
        x, y = proj_2d[i, 0], proj_2d[i, 1]
        cat = r['category']
        s = r['simplicity']
        marker = cat_markers.get(cat, 's')
        band = get_band(s) if cat != 'NOT_CLOSED' else 'other'
        color = band_colors[band]
        edge = 'white' if cat == 'CLOSED_NON_SIMPLE' else '#555555'
        size = 70 if cat == 'CLOSED_NON_SIMPLE' else (50 if cat == 'PARTIAL' else 30)

        ax_main.scatter(x, y, c=color, marker=marker, s=size,
                        edgecolors=edge, linewidths=0.6, zorder=5)

    # Mark rebel seed 100
    rebel_idx = next((i for i, r in enumerate(results) if r['seed'] == 100), None)
    if rebel_idx is not None:
        rx, ry = proj_2d[rebel_idx]
        ax_main.scatter(rx, ry, c='#ff0000', marker='*', s=250,
                        edgecolors='white', linewidths=1.0, zorder=10)
        ax_main.annotate('seed 100\n(rebel)',
                         (rx, ry), (rx + (extent[1] - extent[0]) * 0.04, ry),
                         color='#ff4444', fontsize=8, fontweight='bold',
                         arrowprops=dict(arrowstyle='->', color='#ff4444', lw=0.8),
                         ha='left', va='center', zorder=10)

    ax_main.set_xlabel(f'PC1 ({explained[0]:.1%} variance)', fontsize=12)
    ax_main.set_ylabel(f'PC2 ({explained[1]:.1%} variance)', fontsize=12)
    ax_main.set_title(title, fontsize=15, fontweight='bold', pad=12)

    # Colorbar for elevation
    sm = plt.cm.ScalarMappable(cmap=plt.cm.terrain,
                               norm=Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cax = fig.add_axes([0.08, 0.03, 0.62, 0.02])
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label(r'$\log_{10}\,T_{\mathrm{total}}$  (lower = deeper valley)',
                   fontsize=11, color=fg)
    cbar.ax.tick_params(colors=fg)

    # ------- PCA VARIANCE INSET -------
    ax_pca.bar(range(1, 11), explained[:10] * 100,
               color='#58a6ff', edgecolor='#388bfd', alpha=0.8)
    ax_pca.set_xlabel('PC', fontsize=9)
    ax_pca.set_ylabel('Var. explained (%)', fontsize=9)
    ax_pca.set_title('PCA Spectrum', fontsize=10, fontweight='bold')
    ax_pca.tick_params(labelsize=8)
    # Cumulative line
    ax_pca2 = ax_pca.twinx()
    cumvar = np.cumsum(explained[:10]) * 100
    ax_pca2.plot(range(1, 11), cumvar, 'o-', color='#f0883e',
                 markersize=3, linewidth=1.2)
    ax_pca2.set_ylabel('Cumul. (%)', fontsize=9, color='#f0883e')
    ax_pca2.tick_params(labelsize=8, colors='#f0883e')
    ax_pca2.set_ylim(0, 105)

    # ------- LEGEND PANEL -------
    y0 = 0.95
    dy = 0.055

    ax_leg.text(0.5, y0, 'Simplicity Bands', fontsize=11, fontweight='bold',
                ha='center', va='top', transform=ax_leg.transAxes)
    y0 -= dy

    band_info = [
        ('B1', r'$S\approx 0.08$', '7 seeds (6.8%)'),
        ('B2', r'$S\approx 0.14$', '55 seeds (53.4%)'),
        ('B3', r'$S\approx 0.20$', '18 seeds (17.5%)'),
        ('B4', r'$S\approx 0.25$', '8 seeds (7.8%)'),
        ('B5', r'$S\approx 0.29$', '15 seeds (14.6%)'),
    ]
    for name, s_label, count in band_info:
        ax_leg.scatter(0.08, y0, c=band_colors[name], s=80, marker='o',
                       edgecolors='white', linewidths=0.5,
                       transform=ax_leg.transAxes, zorder=5)
        ax_leg.text(0.18, y0, f'{name}: {s_label}', fontsize=9, va='center',
                    transform=ax_leg.transAxes)
        ax_leg.text(0.18, y0 - 0.025, count, fontsize=7.5, va='center',
                    color='#8b949e', transform=ax_leg.transAxes)
        y0 -= dy * 1.3

    y0 -= dy * 0.5
    ax_leg.text(0.5, y0, 'Category Markers', fontsize=11, fontweight='bold',
                ha='center', va='top', transform=ax_leg.transAxes)
    y0 -= dy

    cat_info = [
        ('o', 'white', 'Closed Non-Simple'),
        ('^', '#555555', 'Partial Closure'),
        ('x', '#555555', 'Not Closed'),
        ('*', 'white', 'Rebel (seed 100)'),
    ]
    cat_colors = ['#f4a261', '#f4a261', '#6b7280', '#ff0000']
    cat_sizes = [80, 60, 40, 150]
    for (mk, ec, label), col, sz in zip(cat_info, cat_colors, cat_sizes):
        ax_leg.scatter(0.08, y0, c=col, s=sz, marker=mk,
                       edgecolors=ec, linewidths=0.6,
                       transform=ax_leg.transAxes, zorder=5)
        ax_leg.text(0.18, y0, label, fontsize=9, va='center',
                    transform=ax_leg.transAxes)
        y0 -= dy

    # Stats footer
    y0 -= dy * 0.5
    n_closed = sum(1 for r in results if r['category'] == 'CLOSED_NON_SIMPLE')
    n_partial = sum(1 for r in results if r['category'] == 'PARTIAL')
    n_notclosed = sum(1 for r in results if r['category'] == 'NOT_CLOSED')
    stats = (
        f"128 seeds in so(27), 52 gens\n"
        f"Closed: {n_closed} | Partial: {n_partial} | Open: {n_notclosed}\n"
        f"F4 found: 0 (P < 2.3% at 95% CI)"
    )
    ax_leg.text(0.5, y0, stats, fontsize=8, va='top', ha='center',
                transform=ax_leg.transAxes, color='#8b949e',
                linespacing=1.5)

    plt.savefig(save_path, dpi=200, facecolor=bg, bbox_inches='tight')
    print(f"\n  Saved: {save_path}")
    plt.close()


# ============================================================
# FIGURE 2: SIMPLICITY MAP (same layout, colored by simplicity)
# ============================================================

def plot_simplicity_map(grid_x, grid_y, simp_grid, proj_2d, results,
                        explained, save_path):
    """Simplicity index as a scalar field over the PCA plane."""

    bg = '#0d1117'
    fg = '#c9d1d9'
    plt.rcParams.update({
        'figure.facecolor': bg, 'axes.facecolor': bg,
        'axes.edgecolor': fg, 'axes.labelcolor': fg,
        'xtick.color': fg, 'ytick.color': fg,
        'text.color': fg, 'font.size': 11, 'font.family': 'sans-serif',
    })

    fig, ax = plt.subplots(figsize=(14, 11))

    extent = [grid_x[0], grid_x[-1], grid_y[0], grid_y[-1]]

    # Clip simplicity to [0, 0.5] for better contrast
    Z = np.clip(simp_grid, 0, 0.5)

    # Hillshaded
    ls = LightSource(azdeg=315, altdeg=35)
    rgb = ls.shade(Z, cmap=plt.cm.RdYlGn_r, vert_exag=0.5, blend_mode='soft')
    ax.imshow(rgb, extent=extent, origin='lower', aspect='auto',
              interpolation='bilinear')

    # Contour lines at band values
    band_levels = [0.08, 0.14, 0.20, 0.247, 0.289]
    cs = ax.contour(grid_x, grid_y, simp_grid, levels=band_levels,
                    colors=['#e63946', '#f4a261', '#2a9d8f', '#457b9d', '#a855f7'],
                    linewidths=1.5, alpha=0.7)
    ax.clabel(cs, fmt='S=%.2f', fontsize=8,
              colors='white')

    # Seeds
    for i, r in enumerate(results):
        x, y = proj_2d[i, 0], proj_2d[i, 1]
        cat = r['category']
        if cat == 'NOT_CLOSED':
            ax.scatter(x, y, c='#6b7280', marker='x', s=25,
                       edgecolors='none', alpha=0.5, zorder=5)
        else:
            ax.scatter(x, y, c='white', marker='o', s=40,
                       edgecolors='#333333', linewidths=0.5, zorder=5)

    ax.set_xlabel(f'PC1 ({explained[0]:.1%} variance)', fontsize=12)
    ax.set_ylabel(f'PC2 ({explained[1]:.1%} variance)', fontsize=12)
    ax.set_title('Simplicity Index Landscape — Band Structure in Parameter Space',
                 fontsize=14, fontweight='bold', pad=12)

    sm = plt.cm.ScalarMappable(cmap=plt.cm.RdYlGn_r,
                               norm=Normalize(vmin=0, vmax=0.5))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label('Simplicity  S = K_std / |K_mean|', fontsize=11, color=fg)
    cbar.ax.tick_params(colors=fg)

    # Mark band levels on colorbar
    for lv, name in zip(band_levels, ['B1', 'B2', 'B3', 'B4', 'B5']):
        cbar.ax.axhline(lv, color='white', linewidth=0.5, alpha=0.5)
        cbar.ax.text(1.15, lv, name, fontsize=8, color='white',
                     va='center', transform=cbar.ax.get_yaxis_transform())

    plt.savefig(save_path, dpi=200, facecolor=bg, bbox_inches='tight')
    print(f"  Saved: {save_path}")
    plt.close()


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Landscape Relief Map — topographic PCA visualization')
    parser.add_argument('--input', type=str, default='landscape_128seeds_results.npz')
    parser.add_argument('--grid-size', type=int, default=120,
                        help='Grid resolution per axis (default: 120)')
    parser.add_argument('--batch-size', type=int, default=25,
                        help='GPU batch size for grid eval (default: 25)')
    parser.add_argument('--padding', type=float, default=0.4,
                        help='Fractional padding beyond seed cloud (default: 0.4)')
    parser.add_argument('--gamma', type=float, default=GAMMA_FINAL,
                        help='Closure weight gamma (default: 200)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory (default: ../figures)')
    args = parser.parse_args()

    print("=" * 70)
    print("  LANDSCAPE RELIEF MAP — Topographic PCA Visualization")
    print("=" * 70)
    print(f"  Grid:    {args.grid_size}×{args.grid_size}")
    print(f"  Batch:   {args.batch_size}")
    print(f"  Gamma:   {args.gamma}")
    print(f"  Device:  {DEVICE}")
    if torch.cuda.is_available():
        print(f"  GPU:     {torch.cuda.get_device_name(0)}")
        mem = torch.cuda.get_device_properties(0).total_memory / 1e9
        print(f"  VRAM:    {mem:.1f} GB")
    print()

    # --- Load and flatten ---
    input_path = os.path.join(os.path.dirname(__file__), args.input)
    flat_params, results, raw_gens = load_and_flatten(input_path)

    # --- PCA ---
    mean, pcs, proj, explained = compute_pca(flat_params)
    pc1, pc2 = pcs[0], pcs[1]
    proj_2d = proj[:, :2]

    # --- Grid bounds ---
    x_min, x_max = proj_2d[:, 0].min(), proj_2d[:, 0].max()
    y_min, y_max = proj_2d[:, 1].min(), proj_2d[:, 1].max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    pad = args.padding
    grid_x = np.linspace(x_min - pad * x_range, x_max + pad * x_range, args.grid_size)
    grid_y = np.linspace(y_min - pad * y_range, y_max + pad * y_range, args.grid_size)

    print(f"\n  Grid range: x=[{grid_x[0]:.2f}, {grid_x[-1]:.2f}], "
          f"y=[{grid_y[0]:.2f}, {grid_y[-1]:.2f}]")

    # --- GPU grid evaluation ---
    T_total_grid, T_cl_grid, simp_grid = evaluate_grid_gpu(
        mean, pc1, pc2, grid_x, grid_y,
        batch_size=args.batch_size, gamma=args.gamma)

    # Log transform
    log_T = np.log10(np.clip(T_total_grid, 1e-10, None))
    log_Tcl = np.log10(np.clip(T_cl_grid, 1e-10, None))

    print(f"\n  T_total range: {T_total_grid.min():.2e} — {T_total_grid.max():.2e}")
    print(f"  log10(T) range: {log_T.min():.2f} — {log_T.max():.2f}")
    print(f"  Simplicity range: {simp_grid.min():.4f} — {simp_grid.max():.4f}")

    # --- Output paths ---
    out_dir = args.output_dir or os.path.join(os.path.dirname(__file__), '..', 'figures')
    os.makedirs(out_dir, exist_ok=True)

    # --- Figure 1: Relief map ---
    print("\n  Rendering relief map...")
    plot_relief_map(
        grid_x, grid_y, log_T, T_cl_grid, simp_grid,
        proj_2d, results, explained,
        os.path.join(out_dir, 'landscape_relief_map.png'),
        title='Algebraic Optimization Landscape — 52 Generators in so(27)')

    # --- Figure 2: Simplicity map ---
    print("  Rendering simplicity map...")
    plot_simplicity_map(
        grid_x, grid_y, simp_grid, proj_2d, results,
        explained,
        os.path.join(out_dir, 'landscape_simplicity_map.png'))

    # --- Save grid data for future use ---
    grid_path = os.path.join(os.path.dirname(__file__),
                             'landscape_relief_grid.npz')
    np.savez_compressed(grid_path,
                        grid_x=grid_x, grid_y=grid_y,
                        T_total=T_total_grid, T_cl=T_cl_grid,
                        simplicity=simp_grid, log_T=log_T,
                        proj_2d=proj_2d, explained=explained[:10],
                        pca_mean=mean, pc1=pc1, pc2=pc2)
    print(f"  Grid data saved: {grid_path}")

    print("\n" + "=" * 70)
    print("  DONE")
    print("=" * 70)


if __name__ == '__main__':
    main()
