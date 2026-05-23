"""Detailed figure: simplicity bands + rebel seed analysis."""
import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ── Load data ──
data = np.load('landscape_128seeds_results.npz', allow_pickle=True)
results = json.loads(str(data['results']))

BG = "#0d1117"
BG_P = "#161b22"
TEXT = "#c9d1d9"
SUBTLE = "#8b949e"
GREEN = "#58d68d"
RED = "#e74c3c"
ORANGE = "#f39c12"
CYAN = "#56c8f5"
PURPLE = "#b388ff"
YELLOW = "#f5e642"

fig = plt.figure(figsize=(18, 12), facecolor=BG)
fig.suptitle("128-Seed Landscape — Band Structure & Rebel Analysis",
             color=TEXT, fontsize=16, fontweight="bold", y=0.97)

# ═══════════════════════════════════════════════════
# Panel 1: Closure vs Simplicity with band highlights
# ═══════════════════════════════════════════════════
ax1 = fig.add_subplot(2, 2, 1, facecolor=BG_P)

cat_colors = {
    "CLOSED_NON_SIMPLE": GREEN,
    "PARTIAL": ORANGE,
    "NOT_CLOSED": RED,
}

for r in results:
    c = cat_colors.get(r['category'], SUBTLE)
    ax1.scatter(r['closure_frac'] * 100, r['simplicity'],
                color=c, s=30, alpha=0.7, edgecolors='none', zorder=3)

# Highlight rebel
rebel = [r for r in results if r['seed'] == 100][0]
ax1.scatter(rebel['closure_frac'] * 100, rebel['simplicity'],
            color=YELLOW, s=120, marker='*', zorder=5, edgecolors='white', linewidths=0.5)
ax1.annotate(f"Seed 100\nhv={rebel['h_dual_ir3']:.1f}",
             xy=(rebel['closure_frac'] * 100, rebel['simplicity']),
             xytext=(75, 0.22), color=YELLOW, fontsize=9, fontweight='bold',
             arrowprops=dict(arrowstyle='->', color=YELLOW, lw=1.5))

# Band lines
bands = [0.0800, 0.1400, 0.2000, 0.2474, 0.2887]
band_labels = ["~1/√156", "~1/√51", "1/5", "~1/√16", "1/√12"]
for b, lbl in zip(bands, band_labels):
    ax1.axhline(b, color=CYAN, alpha=0.25, lw=1, ls='--', zorder=1)
    ax1.text(102, b, f"S={b:.3f} ({lbl})", color=CYAN, fontsize=7,
             va='center', ha='left', alpha=0.7)

ax1.set_xlabel("Closure (%)", color=TEXT, fontsize=10)
ax1.set_ylabel("Simplicity", color=TEXT, fontsize=10)
ax1.set_title("Discrete Simplicity Bands", color=TEXT, fontsize=12)
ax1.set_xlim(-5, 115)
ax1.tick_params(colors=SUBTLE)
for spine in ax1.spines.values():
    spine.set_color(SUBTLE)

# ═══════════════════════════════════════════════════
# Panel 2: Simplicity histogram
# ═══════════════════════════════════════════════════
ax2 = fig.add_subplot(2, 2, 2, facecolor=BG_P)

closed_simps = [r['simplicity'] for r in results if r['closure_frac'] > 0.5]
all_simps = [r['simplicity'] for r in results]

ax2.hist(closed_simps, bins=50, range=(0, 0.65), color=GREEN, alpha=0.7, label="Closed/Partial (>50%)")
ax2.hist([r['simplicity'] for r in results if r['closure_frac'] <= 0.5],
         bins=50, range=(0, 0.65), color=RED, alpha=0.5, label="Not Closed")

for b in bands:
    ax2.axvline(b, color=CYAN, alpha=0.4, lw=1.5, ls='--')

ax2.set_xlabel("Simplicity", color=TEXT, fontsize=10)
ax2.set_ylabel("Count", color=TEXT, fontsize=10)
ax2.set_title("Simplicity Distribution — Sharp Peaks", color=TEXT, fontsize=12)
ax2.legend(fontsize=8, facecolor=BG_P, edgecolor=SUBTLE, labelcolor=TEXT)
ax2.tick_params(colors=SUBTLE)
for spine in ax2.spines.values():
    spine.set_color(SUBTLE)

# ═══════════════════════════════════════════════════
# Panel 3: hv vs closure with rebel highlight
# ═══════════════════════════════════════════════════
ax3 = fig.add_subplot(2, 2, 3, facecolor=BG_P)

for r in results:
    c = cat_colors.get(r['category'], SUBTLE)
    ax3.scatter(r['h_dual_ir3'], r['closure_frac'] * 100,
                color=c, s=30, alpha=0.7, edgecolors='none', zorder=3)

ax3.scatter(rebel['h_dual_ir3'], rebel['closure_frac'] * 100,
            color=YELLOW, s=150, marker='*', zorder=5, edgecolors='white', linewidths=0.5)
ax3.annotate(f"Seed 100\nhv = {rebel['h_dual_ir3']:.2f}\nh_raw = {rebel['h_raw']:.4f}\n"
             f"T_norm = {rebel['T_norm']:.2f}\nK_std = {rebel['K_std']:.4f}",
             xy=(rebel['h_dual_ir3'], rebel['closure_frac'] * 100),
             xytext=(5.5, 60), color=YELLOW, fontsize=8,
             arrowprops=dict(arrowstyle='->', color=YELLOW, lw=1.5),
             bbox=dict(facecolor=BG, edgecolor=YELLOW, alpha=0.8, pad=5))

# Typical cluster
ax3.axvspan(2.5, 3.1, alpha=0.08, color=GREEN, zorder=1)
ax3.text(2.8, 10, "Typical\ncluster", color=GREEN, fontsize=8, ha='center', alpha=0.7)

ax3.set_xlabel("hv (I_R=3)", color=TEXT, fontsize=10)
ax3.set_ylabel("Closure (%)", color=TEXT, fontsize=10)
ax3.set_title("hv Landscape — Seed 100 Rebel (3x typical)", color=TEXT, fontsize=12)
ax3.tick_params(colors=SUBTLE)
for spine in ax3.spines.values():
    spine.set_color(SUBTLE)

# ═══════════════════════════════════════════════════
# Panel 4: Rebel vs Typical — radar-style comparison
# ═══════════════════════════════════════════════════
ax4 = fig.add_subplot(2, 2, 4, facecolor=BG_P)

typicals = [r for r in results if r['category'] == 'CLOSED_NON_SIMPLE' and r['h_dual_ir3'] < 4]

metrics = ['h_raw', 'K_mean', 'K_std', 'max_residual', 'T_norm', 'T_k', 'simplicity']
labels = ['h_raw', '|K_mean|', 'K_std', 'max_resid', 'T_norm', 'T_kill', 'simplicity']

rebel_vals = [abs(rebel[m]) for m in metrics]
typ_means = [abs(np.mean([r[m] for r in typicals])) for m in metrics]

# Normalize to typical
ratios = [rv / tv if tv > 0 else 0 for rv, tv in zip(rebel_vals, typ_means)]

y_pos = np.arange(len(metrics))
bars = ax4.barh(y_pos, ratios, color=[YELLOW if r > 2 or r < 0.5 else GREEN for r in ratios],
                alpha=0.8, height=0.6)

ax4.axvline(1.0, color=TEXT, lw=1, ls='--', alpha=0.5)
ax4.set_yticks(y_pos)
ax4.set_yticklabels(labels, color=TEXT, fontsize=9)
ax4.set_xlabel("Rebel / Typical ratio", color=TEXT, fontsize=10)
ax4.set_title("Seed 100 vs Typical CLOSED_NON_SIMPLE", color=TEXT, fontsize=12)

for i, (ratio, val) in enumerate(zip(ratios, rebel_vals)):
    x_txt = max(ratio + 0.1, 0.3)
    txt = f"{ratio:.2f}x"
    if ratio > 100:
        txt = f"{ratio:.0f}x"
    ax4.text(min(x_txt, ax4.get_xlim()[1] * 0.9), i, txt,
             color=TEXT, fontsize=8, va='center')

ax4.set_xscale('log')
ax4.tick_params(colors=SUBTLE)
for spine in ax4.spines.values():
    spine.set_color(SUBTLE)

plt.tight_layout(rect=[0, 0, 1, 0.94])
out = "landscape_128seeds_bands_rebel.png"
plt.savefig(out, dpi=180, facecolor=BG, bbox_inches='tight')
plt.close()
print(f"Saved: {out}")
