"""
Visual dashboard for MC landscape results.

Reads partial or complete .npz results and generates:
  1. Seed grid colored by category (wafer map style)
  2. Scatter: closure_frac vs simplicity
  3. Multi-metric heatmap per seed
  4. Category frequency bar chart

Usage:
  python _landscape_dashboard.py                                    # default file
  python _landscape_dashboard.py landscape_16seeds_results.npz      # specific file
  python _landscape_dashboard.py --watch                            # auto-refresh every 60s
"""

import sys, os, json, argparse, time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec

# Category palette — ordered from "worst" to "best"
CATEGORY_ORDER = ['NOT_CLOSED', 'PARTIAL', 'CLOSED_NON_SIMPLE', 'SIMPLE_OTHER', 'F4']
CATEGORY_COLORS = {
    'NOT_CLOSED':       '#d62728',   # red
    'PARTIAL':          '#ff7f0e',   # orange
    'CLOSED_NON_SIMPLE':'#2ca02c',   # green
    'SIMPLE_OTHER':     '#1f77b4',   # blue
    'F4':               '#9467bd',   # purple
}
CATEGORY_LABELS = {
    'NOT_CLOSED':       'Not Closed',
    'PARTIAL':          'Partial',
    'CLOSED_NON_SIMPLE':'Closed Non-Simple',
    'SIMPLE_OTHER':     'Simple (other)',
    'F4':               'F₄',
}


def load_results(path):
    """Load results from .npz, returns list of dicts."""
    if not os.path.exists(path):
        print(f"File not found: {path}")
        return []
    data = np.load(path, allow_pickle=True)
    results = json.loads(str(data['results']))
    return results


def category_to_int(cat):
    """Map category string to integer for colormap."""
    if cat in CATEGORY_ORDER:
        return CATEGORY_ORDER.index(cat)
    return 0


def make_grid_data(results, n_seeds_total):
    """Build a 2D grid of category indices for the seed wafer map."""
    # Choose grid dimensions close to square
    cols = int(np.ceil(np.sqrt(n_seeds_total)))
    rows = int(np.ceil(n_seeds_total / cols))

    grid = np.full((rows, cols), np.nan)
    seed_map = {r['seed']: r for r in results}

    for seed in range(n_seeds_total):
        row = seed // cols
        col = seed % cols
        if seed in seed_map:
            grid[row, col] = category_to_int(seed_map[seed]['category'])
        # nan = not yet computed (will be grey)

    return grid, rows, cols


def plot_dashboard(results, n_seeds_total, title="", save_path=None):
    """Generate the full dashboard figure."""
    if not results:
        print("No results to plot.")
        return

    fig = plt.figure(figsize=(18, 12), facecolor='#1a1a2e')
    fig.suptitle(
        f"MC Landscape Dashboard — {len(results)}/{n_seeds_total} seeds\n{title}",
        fontsize=16, color='white', fontweight='bold', y=0.98
    )

    gs = gridspec.GridSpec(2, 3, hspace=0.35, wspace=0.3,
                           left=0.06, right=0.96, top=0.92, bottom=0.06)

    dark_bg = '#16213e'
    text_color = '#e0e0e0'

    # ================================================================
    # 1. SEED GRID (wafer map)
    # ================================================================
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax1.set_facecolor(dark_bg)
    ax1.set_title('Seed Grid — Category Map', color=text_color, fontsize=13, pad=10)

    grid, rows, cols = make_grid_data(results, n_seeds_total)

    # Custom colormap: categories + NaN as grey
    cmap_colors = [CATEGORY_COLORS[c] for c in CATEGORY_ORDER]
    cmap = ListedColormap(cmap_colors)
    bounds = np.arange(-0.5, len(CATEGORY_ORDER) + 0.5, 1)
    norm = BoundaryNorm(bounds, cmap.N)

    # Mask NaN for "not computed" seeds
    masked_grid = np.ma.masked_invalid(grid)
    cmap.set_bad(color='#2a2a3a')

    im = ax1.imshow(masked_grid, cmap=cmap, norm=norm, aspect='equal',
                     interpolation='nearest')

    # Seed numbers inside cells
    for seed in range(n_seeds_total):
        row = seed // cols
        col = seed % cols
        color = '#ffffff' if not np.isnan(grid[row, col]) else '#555555'
        fontsize = 7 if n_seeds_total > 64 else 9
        ax1.text(col, row, str(seed), ha='center', va='center',
                 fontsize=fontsize, color=color, fontweight='bold')

    ax1.set_xticks([])
    ax1.set_yticks([])

    # Legend
    legend_elements = [
        Patch(facecolor=CATEGORY_COLORS[c], label=CATEGORY_LABELS[c])
        for c in CATEGORY_ORDER
    ]
    legend_elements.append(Patch(facecolor='#2a2a3a', label='Pending'))
    ax1.legend(handles=legend_elements, loc='upper right',
               fontsize=8, facecolor=dark_bg, edgecolor='#444',
               labelcolor=text_color, framealpha=0.9)

    # ================================================================
    # 2. CATEGORY BAR CHART
    # ================================================================
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.set_facecolor(dark_bg)
    ax2.set_title('Category Frequencies', color=text_color, fontsize=13, pad=10)

    counts = {c: 0 for c in CATEGORY_ORDER}
    for r in results:
        cat = r['category']
        if cat in counts:
            counts[cat] += 1

    cats = CATEGORY_ORDER
    vals = [counts[c] for c in cats]
    colors = [CATEGORY_COLORS[c] for c in cats]
    labels = [CATEGORY_LABELS[c] for c in cats]

    bars = ax2.barh(range(len(cats)), vals, color=colors, edgecolor='#333', height=0.7)
    ax2.set_yticks(range(len(cats)))
    ax2.set_yticklabels(labels, color=text_color, fontsize=10)
    ax2.set_xlabel('Count', color=text_color, fontsize=10)
    ax2.tick_params(axis='x', colors=text_color)
    ax2.invert_yaxis()

    # Percentage labels
    n = len(results)
    for i, (bar, v) in enumerate(zip(bars, vals)):
        pct = v / n * 100 if n > 0 else 0
        ax2.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                 f'{v} ({pct:.0f}%)', va='center', color=text_color, fontsize=9)

    ax2.set_xlim(0, max(vals) * 1.3 + 1 if vals else 1)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_color('#444')
    ax2.spines['left'].set_color('#444')

    # ================================================================
    # 3. SCATTER: closure_frac vs simplicity
    # ================================================================
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.set_facecolor(dark_bg)
    ax3.set_title('Closure vs Simplicity', color=text_color, fontsize=13, pad=10)

    for cat in CATEGORY_ORDER:
        cat_results = [r for r in results if r['category'] == cat]
        if not cat_results:
            continue
        x = [r['closure_frac'] * 100 for r in cat_results]
        y = [min(r['simplicity'], 5.0) for r in cat_results]  # clip for visibility
        ax3.scatter(x, y, c=CATEGORY_COLORS[cat], label=CATEGORY_LABELS[cat],
                    s=40, alpha=0.8, edgecolors='white', linewidth=0.3, zorder=3)

    ax3.set_xlabel('Closure (%)', color=text_color, fontsize=10)
    ax3.set_ylabel('Simplicity (lower=better)', color=text_color, fontsize=10)
    ax3.tick_params(colors=text_color)
    ax3.set_xlim(-5, 105)
    ax3.axhline(y=0.05, color='#666', linestyle='--', linewidth=0.8, label='Simple threshold')
    ax3.axvline(x=99, color='#666', linestyle=':', linewidth=0.8, label='Closed threshold')
    ax3.legend(fontsize=7, facecolor=dark_bg, edgecolor='#444',
               labelcolor=text_color, framealpha=0.9, loc='upper left')
    for spine in ax3.spines.values():
        spine.set_color('#444')

    # ================================================================
    # 4. SCATTER: h∨(I_R=3) vs T_cl
    # ================================================================
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.set_facecolor(dark_bg)
    ax4.set_title('h∨(I_R=3) vs T_closure', color=text_color, fontsize=13, pad=10)

    for cat in CATEGORY_ORDER:
        cat_results = [r for r in results if r['category'] == cat]
        if not cat_results:
            continue
        x = [r['h_dual_ir3'] for r in cat_results]
        y = [r['T_cl'] for r in cat_results]
        ax4.scatter(x, y, c=CATEGORY_COLORS[cat], label=CATEGORY_LABELS[cat],
                    s=40, alpha=0.8, edgecolors='white', linewidth=0.3, zorder=3)

    ax4.set_xlabel('h∨ (I_R=3)', color=text_color, fontsize=10)
    ax4.set_ylabel('T_closure', color=text_color, fontsize=10)
    ax4.set_yscale('log')
    ax4.tick_params(colors=text_color)
    # Mark F4 expected h∨
    ax4.axvline(x=9.0, color='#9467bd', linestyle='--', linewidth=0.8,
                alpha=0.6, label='F₄ h∨=9')
    ax4.legend(fontsize=7, facecolor=dark_bg, edgecolor='#444',
               labelcolor=text_color, framealpha=0.9)
    for spine in ax4.spines.values():
        spine.set_color('#444')

    # ================================================================
    # 5. MULTI-METRIC HEATMAP
    # ================================================================
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.set_facecolor(dark_bg)
    ax5.set_title('Per-Seed Metrics', color=text_color, fontsize=13, pad=10)

    if len(results) > 0:
        # Sort by category for visual coherence
        sorted_results = sorted(results,
                                key=lambda r: (category_to_int(r['category']),
                                               -r['closure_frac']))
        metrics = ['closure_frac', 'simplicity', 'T_cl', 'T_k']
        metric_labels = ['Closure', 'Simplicity', 'T_cl', 'T_k']

        data_matrix = np.zeros((len(sorted_results), len(metrics)))
        for i, r in enumerate(sorted_results):
            for j, m in enumerate(metrics):
                val = r.get(m, 0)
                data_matrix[i, j] = val

        # Normalize each column to [0, 1] for comparable heatmap
        for j in range(len(metrics)):
            col = data_matrix[:, j]
            mn, mx = col.min(), col.max()
            if mx > mn:
                data_matrix[:, j] = (col - mn) / (mx - mn)
            else:
                data_matrix[:, j] = 0.5

        im5 = ax5.imshow(data_matrix.T, aspect='auto', cmap='viridis',
                          interpolation='nearest')
        ax5.set_yticks(range(len(metrics)))
        ax5.set_yticklabels(metric_labels, color=text_color, fontsize=10)
        ax5.set_xlabel('Seeds (sorted by category)', color=text_color, fontsize=10)
        ax5.tick_params(axis='x', colors=text_color)
        ax5.tick_params(axis='y', colors=text_color)

        # Category color bar on top
        cat_colors_bar = [CATEGORY_COLORS[r['category']] for r in sorted_results]
        for i, c in enumerate(cat_colors_bar):
            ax5.plot(i, -0.7, 's', color=c, markersize=3, clip_on=False)

    plt.savefig(save_path or 'landscape_dashboard.png', dpi=150,
                facecolor=fig.get_facecolor())
    print(f"  Dashboard saved: {save_path or 'landscape_dashboard.png'}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='MC Landscape visual dashboard')
    parser.add_argument('file', nargs='?',
                        default='landscape_128seeds_results.npz',
                        help='Results .npz file')
    parser.add_argument('--n-seeds', type=int, default=None,
                        help='Total seeds expected (auto-detect if omitted)')
    parser.add_argument('--out', type=str, default=None,
                        help='Output PNG path')
    parser.add_argument('--watch', action='store_true',
                        help='Auto-refresh every 60s')
    parser.add_argument('--interval', type=int, default=60,
                        help='Refresh interval in seconds (with --watch)')
    args = parser.parse_args()

    save_dir = os.path.dirname(os.path.abspath(__file__))
    results_path = os.path.join(save_dir, args.file)
    out_path = args.out or os.path.join(save_dir,
                                         args.file.replace('.npz', '_dashboard.png'))

    while True:
        results = load_results(results_path)
        if results:
            n_seeds = args.n_seeds
            if n_seeds is None:
                # Try to infer from config or max seed
                try:
                    data = np.load(results_path, allow_pickle=True)
                    if 'config' in data:
                        cfg = json.loads(str(data['config']))
                        n_seeds = cfg.get('n_seeds', None)
                    if n_seeds is None and 'n_seeds_total' in data:
                        n_seeds = int(data['n_seeds_total'])
                except Exception:
                    pass
                if n_seeds is None:
                    n_seeds = max(r['seed'] for r in results) + 1

            ts = time.strftime('%H:%M:%S')
            plot_dashboard(results, n_seeds,
                           title=f"Updated {ts} | {os.path.basename(args.file)}",
                           save_path=out_path)
        else:
            print(f"  No results yet in {results_path}")

        if not args.watch:
            break

        print(f"  Waiting {args.interval}s for next update... (Ctrl+C to stop)")
        try:
            time.sleep(args.interval)
        except KeyboardInterrupt:
            print("\n  Stopped.")
            break


if __name__ == "__main__":
    main()
