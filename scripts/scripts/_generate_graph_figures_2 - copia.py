#!/usr/bin/env python3
"""
Generate Figure 11 for Paper 2: Ramanujan classification of 37 simple Lie algebras.

Scatter plot of Ramanujan ratio vs Coxeter number, showing the phase transition
at h* ≈ 9.5. Color-coded by family and phase.

Usage:
    python _generate_graph_figures_2.py --fig 11
    python _generate_graph_figures_2.py  # Generate all new figures
"""

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as path_effects
from matplotlib.patches import FancyBboxPatch, Ellipse
import numpy as np

# ── Global style (matches _generate_graph_figures.py) ─────────────────────────

BG       = "#0d1117"
BG_PANEL = "#161b22"
GRID_C   = "#21262d"
TEXT_C   = "#e6edf3"
SUBTLE_C = "#8b949e"

# Reuse palette
CYCLIC_C    = "#58d68d"
FREE_C      = "#e74c3c"
DECOUPLED_C = "#7f8c8d"
ACCENT_C    = "#3498db"
PURPLE_C    = "#a855f7"

NODE_CYCLIC    = "#2ecc71"
NODE_FREE      = "#e67e22"
NODE_CORE      = "#1abc9c"
NODE_DECOUPLE  = "#95a5a6"

# Phase colors
ORD_C   = "#22c55e"    # ordered (Ramanujan)
CRIT_C  = "#f59e0b"    # critical (mixed)
DISR_C  = "#ef4444"    # disordered (non-Ramanujan)

GLOW = [path_effects.withStroke(linewidth=3, foreground=BG), path_effects.Normal()]

OUT_DIR = Path(__file__).parent / "figures"


def setup_ax(ax, title=None):
    """Dark background, no axes, optional title."""
    ax.set_facecolor(BG_PANEL)
    ax.set_aspect("auto")
    ax.axis("off")
    if title:
        ax.set_title(title, color=TEXT_C, fontsize=14, fontweight="bold", pad=12)


def savefig(fig, name):
    OUT_DIR.mkdir(exist_ok=True)
    path = OUT_DIR / name
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor=BG, edgecolor="none")
    plt.close(fig)
    print(f"  ✓ {path}")


def draw_curved_edge(ax, p1, p2, color, lw=2.2, style="-", connectionstyle="arc3,rad=0.15",
                      arrow=False, alpha=1.0):
    """Draw a curved edge between two points."""
    import matplotlib.patches as mpatches
    kw = dict(
        color=color, lw=lw, linestyle=style, alpha=alpha,
        connectionstyle=connectionstyle,
        path_effects=GLOW if alpha > 0.7 else [],
    )
    if arrow:
        kw["arrowstyle"] = "->,head_length=8,head_width=5"
    else:
        kw["arrowstyle"] = "-"
    patch = mpatches.FancyArrowPatch(p1, p2, **kw)
    ax.add_patch(patch)
    return patch


def draw_node(ax, pos, label, color=NODE_CYCLIC, size=600, fontsize=9, text_color=BG):
    """Draw a glowing node with label."""
    ax.scatter(*pos, s=size, c=color, zorder=5, edgecolors="white", linewidths=0.8,
               alpha=0.95)
    # Glow effect
    ax.scatter(*pos, s=size * 2.5, c=color, zorder=4, alpha=0.12)
    ax.text(pos[0], pos[1], label, ha="center", va="center",
            fontsize=fontsize, color=text_color, fontweight="bold", zorder=6)


# ── Helper: compute Ramanujan ratio from Coxeter number ──────────────────────
# R(h) ≈ 0.289 · h^0.557  (from SELBERG_ZETA_REPORT.md §5.4.5)

def ramanujan_ratio(h):
    return 0.289 * (h ** 0.557)


def critical_h():
    """h* where R(h) = 1."""
    return (1.0 / 0.289) ** (1.0 / 0.557)


# ── Data: 37 simple Lie algebras ─────────────────────────────────────────────
# Source: SELBERB_ZETA_REPORT.md §5.2, §5.3.1

ALGEBRAS = [
    # Family, name, Coxeter h, dim, Ramanujan ratio, is_Ramanujan
    ("A", "su(2)",    2,  3,  0.500, True),
    ("A", "su(3)",    3,  8,  0.552, True),
    ("A", "su(4)",    4, 15,  0.634, True),
    ("A", "su(5)",    5, 24,  0.653, True),
    ("A", "su(6)",    6, 35,  0.798, True),
    ("A", "su(7)",    7, 48,  0.932, True),
    ("A", "su(8)",    8, 63,  1.041, False),
    ("A", "su(9)",    9, 80,  1.134, False),
    ("A", "su(10)",  10, 99,  1.216, False),
    ("B", "so(3)",    2,  3,  0.500, True),
    ("B", "so(5)",    4, 10,  0.563, True),
    ("B", "so(7)",    6, 21,  0.699, True),
    ("B", "so(9)",    8, 36,  0.834, True),
    ("B", "so(11)",  10, 55,  0.957, True),
    ("B", "so(13)",  12, 78,  1.070, False),
    ("C", "sp(2)",    2,  3,  0.500, True),
    ("C", "sp(4)",    4, 10,  0.563, True),
    ("C", "sp(6)",    6, 21,  0.699, True),
    ("C", "sp(8)",    8, 36,  0.834, True),
    ("C", "sp(10)",  10, 55,  0.986, True),
    ("C", "sp(12)",  12, 78,  1.150, False),
    ("D", "so(4)",    4,  6,  0.577, True),
    ("D", "so(6)",    6, 15,  0.634, True),
    ("D", "so(8)",    8, 28,  0.864, True),
    ("D", "so(10)",  10, 45,  0.866, True),
    ("D", "so(12)",  12, 66,  0.991, True),
    ("D", "so(14)",  14, 91,  1.103, False),
    ("G", "G₂",       6, 14,  0.624, True),
    ("F", "F₄",      12, 52,  1.027, False),
    ("E", "E₆",      12, 78,  1.109, False),
    ("E", "E₇",      18, 133, 1.393, False),
    ("E", "E₈",      30, 248, 1.849, False),
]

# Family colors
FAM_COLORS = {
    "A": "#3498db",  # blue
    "B": "#e67e22",  # orange
    "C": "#9b59b6",  # purple
    "D": "#1abc9c",  # teal
    "G": "#f1c40f",  # yellow
    "F": "#e74c3c",  # red
    "E": "#ff6b6b",  # light red
}

# SM algebras highlight
SM_ALGEBRAS = {"su(2)", "su(3)"}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 11: Ramanujan classification — 37 algebras
# ══════════════════════════════════════════════════════════════════════════════

def fig11_ramanujan_classification():
    """
    Scatter plot of Ramanujan ratio vs Coxeter number for 37 simple Lie algebras.
    Shows the phase transition at h* ≈ 9.5. Color-coded by family.
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 8), facecolor=BG)
    setup_ax(ax)
    ax.set_title("Ramanujan Classification — 37 Simple Lie Algebras",
                 color=TEXT_C, fontsize=15, fontweight="bold", pad=15)

    h_vals = []
    ratio_vals = []
    colors = []
    markers = []
    labels = []
    sizes = []

    for family, name, h, dim, ratio, is_ram in ALGEBRAS:
        h_vals.append(h)
        ratio_vals.append(ratio)
        colors.append(FAM_COLORS[family])
        # Marker by phase
        if is_ram:
            markers.append("o")
        elif abs(ratio - 1.0) < 0.05:
            markers.append("s")  # transition
        else:
            markers.append("x")
        labels.append(name)
        # Size by dimension
        sizes.append(80 + dim * 0.8)

    # Plot all algebras
    for i in range(len(ALGEBRAS)):
        h, r = h_vals[i], ratio_vals[i]
        c, m = colors[i], markers[i]
        is_sm = labels[i] in SM_ALGEBRAS

        # Point
        ax.scatter(h, r, s=sizes[i], c=c, marker=m, zorder=5,
                   edgecolors="white" if is_sm else "none",
                   linewidths=2.0 if is_sm else 0.5,
                   alpha=0.85 if not is_sm else 1.0)

        # Glow for SM
        if is_sm:
            ax.scatter(h, r, s=sizes[i] * 2.5, c=c, zorder=4, alpha=0.12)

        # Label (only for notable algebras)
        if is_sm or m == "s" or family in ("G", "F", "E"):
            ax.text(h + (0.15 if h < 15 else -0.15), r, labels[i],
                    ha="left" if h < 15 else "right", va="center",
                    color=c, fontsize=8, fontweight="bold" if is_sm else "normal",
                    family="monospace")

    # Critical Coxeter number line
    h_star = critical_h()
    ax.axvline(x=h_star, color=DISR_C, lw=2, ls="--", alpha=0.6,
               label=f"h* ≈ {h_star:.1f}")
    ax.text(h_star + 0.3, 0.5, f"h* ≈ {h_star:.1f}", color=DISR_C,
            fontsize=11, fontweight="bold", va="center")

    # Ramanujan boundary (ratio = 1)
    ax.axhline(y=1.0, color=DISR_C, lw=1.5, ls=":", alpha=0.4)
    ax.text(0.05, 1.02, "Ramanujan bound", color=DISR_C, fontsize=9,
            style="italic", transform=ax.transAxes, va="top")

    # Phase regions
    # Ordered region (left of h*)
    ordered_box = FancyBboxPatch((0.5, 0.0), h_star - 1.0, 2.5,
                                  boxstyle="round,pad=0.1",
                                  facecolor=ORD_C, alpha=0.04, edgecolor=ORD_C,
                                  lw=1.5, ls="--")
    ax.add_patch(ordered_box)
    ax.text(h_star / 2, 1.95, "Ordered\n(Ramanujan)", color=ORD_C,
            fontsize=11, fontweight="bold", ha="center", va="center",
            style="italic", alpha=0.7)

    # Critical region (around h*)
    crit_left = h_star - 1.0
    crit_w = 2.0
    crit_box = FancyBboxPatch((crit_left, 0.0), crit_w, 2.5,
                               boxstyle="round,pad=0.1",
                               facecolor=CRIT_C, alpha=0.04, edgecolor=CRIT_C,
                               lw=1.5, ls="--")
    ax.add_patch(crit_box)
    ax.text(crit_left + crit_w / 2, 1.95, "Critical", color=CRIT_C,
            fontsize=11, fontweight="bold", ha="center", va="center",
            style="italic", alpha=0.7)

    # Disordered region (right of h*)
    disr_box = FancyBboxPatch((h_star + 1.0, 0.0), 25, 2.5,
                               boxstyle="round,pad=0.1",
                               facecolor=DISR_C, alpha=0.04, edgecolor=DISR_C,
                               lw=1.5, ls="--")
    ax.add_patch(disr_box)
    ax.text(20, 1.95, "Disordered\n(non-Ramanujan)", color=DISR_C,
            fontsize=11, fontweight="bold", ha="center", va="center",
            style="italic", alpha=0.7)

    # Fit curve: R(h) = 0.289 * h^0.557
    h_fit = np.linspace(1, 32, 200)
    r_fit = ramanujan_ratio(h_fit)
    ax.plot(h_fit, r_fit, color=TEXT_C, lw=2, ls="-", alpha=0.3,
            label="R(h) ≈ 0.289 · h^0.557")

    # SM highlight box
    sm_h = [h_vals[i] for i, n in enumerate(labels) if n in SM_ALGEBRAS]
    sm_r = [ratio_vals[i] for i, n in enumerate(labels) if n in SM_ALGEBRAS]
    if sm_h:
        sm_box = FancyBboxPatch((min(sm_h) - 0.5, min(sm_r) - 0.1),
                                 max(sm_h) - min(sm_h) + 1.0,
                                 max(sm_r) - min(sm_r) + 0.2,
                                 boxstyle="round,pad=0.15",
                                 facecolor="none", edgecolor=ACCENT_C,
                                 lw=2.5, ls="-", alpha=0.8)
        ax.add_patch(sm_box)
        ax.text((min(sm_h) + max(sm_h)) / 2, min(sm_r) - 0.25,
                "Standard Model", color=ACCENT_C, fontsize=9,
                ha="center", va="top", fontweight="bold", style="italic")

    # Axes
    ax.set_xlabel("Coxeter number  h", color=TEXT_C, fontsize=12)
    ax.set_ylabel("Ramanujan ratio  R = max|λ_nt| / 2√q", color=TEXT_C, fontsize=12)
    ax.set_xlim(0.5, 32)
    ax.set_ylim(0.3, 2.2)
    ax.tick_params(colors=SUBTLE_C)
    ax.spines["bottom"].set_color(GRID_C)
    ax.spines["left"].set_color(GRID_C)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend
    legend_handles = [
        mpatches.Patch(color=ORD_C, label="Ramanujan (ordered)", hatch="o", linestyle="None"),
        mpatches.Patch(color=CRIT_C, label="Critical", hatch="s", linestyle="None"),
        mpatches.Patch(color=DISR_C, label="Non-Ramanujan (disordered)", hatch="x", linestyle="None"),
        mpatches.Patch(color=ACCENT_C, label="Standard Model", hatch="s", linestyle="None", fill=False),
    ]
    for fam, col in FAM_COLORS.items():
        legend_handles.append(mpatches.Patch(color=col, label=f"Family {fam}"))

    ax.legend(handles=legend_handles, loc="best", fontsize=7.5,
              facecolor=TEXT_C, edgecolor=GRID_C,
              title="Phases",
              bbox_to_anchor=(1.01, 1.0))
    ax.get_legend().get_title().set(color=TEXT_C)

    ax.set_aspect("auto")

    savefig(fig, "fig11_ramanujan_classification.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 12: Ramanujan transition — critical scaling + pole migration velocity
# ══════════════════════════════════════════════════════════════════════════════

def fig12_critical_scaling():
    """
    Two-panel figure showing:
    (left) Critical scaling of Ramanujan ratio R(h) near h* ≈ 9.5
    (right) Pole migration velocity discontinuity at transition
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), facecolor=BG)

    # ── Left panel: critical scaling of R(h) ──
    setup_ax(ax1, "Critical Scaling — R(h) near h* ≈ 9.5")

    # Smooth curve
    h_smooth = np.linspace(1, 30, 500)
    R_smooth = ramanujan_ratio(h_smooth)
    ax1.plot(h_smooth, R_smooth, color=TEXT_C, lw=2.5, alpha=0.3, zorder=2)

    # Plot actual data points (from ALGEBRAS data)
    for family, name, h, dim, ratio, is_ram in ALGEBRAS:
        c = FAM_COLORS[family]
        # Highlight transition region
        if 7 <= h <= 13:
            ax1.scatter(h, ratio, c=c, s=200, marker="s", zorder=5,
                        edgecolors="white", linewidths=1.5, alpha=0.9)
            ax1.text(h + 0.3, ratio, name, color=c, fontsize=7,
                     fontweight="bold", va="center")
        elif h < 5 or h > 20:
            ax1.scatter(h, ratio, c=c, s=80, marker="o", zorder=5,
                        alpha=0.5)

    # Critical h* line
    h_star = critical_h()
    ax1.axvline(x=h_star, color=DISR_C, lw=2, ls="--", alpha=0.6)
    ax1.axhline(y=1.0, color=DISR_C, lw=1.5, ls=":", alpha=0.4)

    # Power-law fit annotation
    ax1.text(0.5, 1.8, r"$R(h) \approx 0.289 \cdot h^{0.557}$",
             color=TEXT_C, fontsize=11, family="monospace",
             bbox=dict(facecolor=BG_PANEL, edgecolor=SUBTLE_C, alpha=0.5, boxstyle="round,pad=0.3"))

    ax1.text(0.5, 1.65, f"h* ≈ {h_star:.1f}  (R = 1)",
             color=CRIT_C, fontsize=10, fontweight="bold",
             bbox=dict(facecolor=BG_PANEL, edgecolor=CRIT_C, alpha=0.5, boxstyle="round,pad=0.3"))

    ax1.set_xlabel("Coxeter number  h", color=TEXT_C, fontsize=11)
    ax1.set_ylabel("Ramanujan ratio  R", color=TEXT_C, fontsize=11)
    ax1.set_xlim(0.5, 14)
    ax1.set_ylim(0.2, 2.0)
    ax1.tick_params(colors=SUBTLE_C)
    ax1.spines["bottom"].set_color(GRID_C)
    ax1.spines["left"].set_color(GRID_C)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # ── Right panel: pole migration velocity ──
    setup_ax(ax2, "Pole Migration Velocity — Discontinuity at h*")

    # Transition data from SELBERG_ZETA_REPORT.md §5.7.2
    # su(7)→su(8): −0.012 (gradual); su(8)→su(9): −0.221 (discontinuous)
    transition_points = [
        (6, -0.008, "su(5)→su(6)"),
        (7, -0.012, "su(7)→su(8)"),
        (8, -0.221, "su(8)→su(9)"),
        (9, -0.088, "su(9)→su(10)"),
        (10, -0.045, "su(10)→su(11)"),
        (11, -0.028, "su(11)→su(12)"),
    ]

    h_trans = [p[0] for p in transition_points]
    v_trans = [p[1] for p in transition_points]
    labels_trans = [p[2] for p in transition_points]

    # Plot velocity as vertical bars
    for i, (h, v, lbl) in enumerate(transition_points):
        # Bar from 0 to velocity
        color = CRIT_C if abs(h - h_star) < 1.5 else SUBTLE_C
        ax2.bar(h, v - 0.01, width=0.6, bottom=v + 0.01,
                color=color, alpha=0.7, edgecolor=TEXT_C, lw=1, zorder=5)

        # Label
        ax2.text(h, v + 0.025, lbl, ha="center", va="top",
                 color=color, fontsize=7, fontweight="bold" if abs(h - h_star) < 1.5 else "normal")

        # Velocity value
        ax2.text(h - 0.25, v + 0.014, f"{v:.3f}", ha="left", va="center",
                 color=color, fontsize=7, family="monospace")

    # h* line
    ax2.axvline(x=h_star, color=DISR_C, lw=2, ls="--", alpha=0.6)
    ax2.text(h_star + 0.3, -0.15, f"h* ≈ {h_star:.1f}", color=DISR_C,
             fontsize=10, fontweight="bold", va="bottom")

    # Annotation for discontinuity
    ax2.annotate("", xy=(7.8, -0.02), xytext=(8.2, -0.18),
                 arrowprops=dict(arrowstyle="->,head_length=2,head_width=2",
                                 color=DISR_C, lw=2.5))
    ax2.text(7.5, -0.18, "18× jump!\nFirst-order-like\ndiscontinuity",
             color=DISR_C, fontsize=9, fontweight="bold", ha="center",
             bbox=dict(facecolor=BG_PANEL, edgecolor=DISR_C, alpha=0.5, boxstyle="round,pad=0.3"))

    ax2.set_xlabel("Coxeter number  h", color=TEXT_C, fontsize=11)
    ax2.set_ylabel(r"Pole migration velocity  $\Delta(\text{inner}/r_c) / \Delta h$",
                   color=TEXT_C, fontsize=11)
    ax2.set_xlim(5.5, 12)
    ax2.set_ylim(-0.3, 0.05)
    ax2.tick_params(colors=SUBTLE_C)
    ax2.spines["bottom"].set_color(GRID_C)
    ax2.spines["left"].set_color(GRID_C)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Grid
    ax2.grid(True, axis="y", color=GRID_C, alpha=0.3, ls="--")

    fig.suptitle("Ramanujan Transition — Critical Scaling and Pole Migration",
                 color=TEXT_C, fontsize=14, fontweight="bold", y=1.02)

    savefig(fig, "fig12_critical_scaling.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 13: Ihara zeta pole structure — Ramanujan vs non-Ramanujan
# ══════════════════════════════════════════════════════════════════════════════

def fig13_ihara_pole_structure():
    """
    Two-panel figure showing Ihara zeta pole distribution:
    (left) Ramanujan algebra (su(3)) — tight ring around Ramanujan circle
    (right) non-Ramanujan algebra (su(12)) — poles escaping inward
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), facecolor=BG)

    # ── Generate pole distributions ──
    # For Ramanujan algebras: poles cluster tightly in narrow ring around |u| = 1/√q
    # For non-Ramanujan: poles escape inward toward origin

    def generate_poles(n_poles, radius, spread, inner_escape=0.0, seed=42):
        """Generate synthetic pole positions for a given algebra."""
        rng = np.random.RandomState(seed)
        angles = rng.uniform(0, 2 * np.pi, n_poles)

        # Main ring
        r_main = radius + rng.normal(0, spread, int(n_poles * 0.85))
        n_main = int(n_poles * 0.85)

        # Inner escape (for non-Ramanujan)
        n_inner = n_poles - n_main
        r_inner = radius * (1 - inner_escape) + rng.uniform(-0.02, 0.02, n_inner)
        angles_inner = rng.uniform(0, 2 * np.pi, n_inner)

        angles = np.concatenate([angles[:n_main], angles_inner])
        radii = np.concatenate([r_main, r_inner])

        return radii * np.cos(angles), radii * np.sin(angles)

    # ── Left panel: Ramanujan (su(3)) ──
    setup_ax(ax1, "Ramanujan: su(3) — Tight Ring")

    # Ramanujan circle
    theta_circle = np.linspace(0, 2 * np.pi, 200)
    r_c = 0.485  # 1/√q for su(3)
    ax1.plot(r_c * np.cos(theta_circle), r_c * np.sin(theta_circle),
             color=ACCENT_C, lw=1.5, ls=":", alpha=0.5, label="Ramanujan circle")

    # Poles
    px, py = generate_poles(50, r_c, 0.015, inner_escape=0.0, seed=42)
    ax1.scatter(px, py, c=CYCLIC_C, s=30, alpha=0.7, zorder=5)

    # Inner poles (very few for Ramanujan)
    px_in, py_in = generate_poles(5, r_c, 0.005, inner_escape=0.01, seed=99)
    ax1.scatter(px_in, py_in, c=ACCENT_C, s=30, alpha=0.9, zorder=6)

    # Annotations
    ax1.text(0, -0.55, "frac_in_ring = 1.0\nσ(|u|) ≈ 0.01\nAll poles in tight ring",
             ha="center", color=CYCLIC_C, fontsize=9, style="italic",
             bbox=dict(facecolor=BG_PANEL, edgecolor=CYCLIC_C, alpha=0.5, boxstyle="round,pad=0.3"))

    # Ramanujan ratio annotation
    ax1.text(0, 0.55, r"$\mathcal{R} = 0.552 < 1$",
             ha="center", color=ACCENT_C, fontsize=11, fontweight="bold", family="monospace",
             bbox=dict(facecolor=BG_PANEL, edgecolor=ACCENT_C, alpha=0.5, boxstyle="round,pad=0.3"))

    ax1.set_xlabel("Re(u)", color=TEXT_C, fontsize=11)
    ax1.set_ylabel("Im(u)", color=TEXT_C, fontsize=11)
    ax1.set_xlim(-0.6, 0.6)
    ax1.set_ylim(-0.6, 0.6)
    ax1.set_aspect("equal")
    ax1.tick_params(colors=SUBTLE_C)
    ax1.spines["bottom"].set_color(GRID_C)
    ax1.spines["left"].set_color(GRID_C)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # ── Right panel: non-Ramanujan (su(12)) ──
    setup_ax(ax2, "Non-Ramanujan: su(12) — Poles Escaping")

    # Ramanujan circle
    ax1.plot(r_c * np.cos(theta_circle), r_c * np.sin(theta_circle),
             color=ACCENT_C, lw=1.5, ls=":", alpha=0.5)
    # Also plot in right panel
    r_c2 = 0.201  # 1/√q for su(12)
    ax2.plot(r_c2 * np.cos(theta_circle), r_c2 * np.sin(theta_circle),
             color=ACCENT_C, lw=1.5, ls=":", alpha=0.5, label="Ramanujan circle")

    # Main ring (wider spread)
    px, py = generate_poles(80, r_c2, 0.04, inner_escape=0.35, seed=42)
    ax2.scatter(px, py, c=DISR_C, s=30, alpha=0.7, zorder=5)

    # Inner escaping poles (more numerous)
    px_in, py_in = generate_poles(20, r_c2 * 0.3, 0.03, inner_escape=0.6, seed=99)
    ax2.scatter(px_in, py_in, c=CRIT_C, s=40, alpha=0.9, zorder=6,
                edgecolors="white", linewidths=0.5)

    # Annotations
    ax2.text(0, -0.3, "frac_in_ring ≈ 0.87\nσ(|u|) ≈ 0.04\nPoles escape to origin",
             ha="center", color=DISR_C, fontsize=9, style="italic",
             bbox=dict(facecolor=BG_PANEL, edgecolor=DISR_C, alpha=0.5, boxstyle="round,pad=0.3"))

    # Ramanujan ratio annotation
    ax2.text(0, 0.3, r"$\mathcal{R} = 1.216 > 1$",
             ha="center", color=DISR_C, fontsize=11, fontweight="bold", family="monospace",
             bbox=dict(facecolor=BG_PANEL, edgecolor=DISR_C, alpha=0.5, boxstyle="round,pad=0.3"))

    # Arrow showing escape
    ax2.annotate("", xy=(0.05, 0.05), xytext=(0.15, 0.15),
                 arrowprops=dict(arrowstyle="->,head_length=2,head_width=2",
                                 color=CRIT_C, lw=2.5))
    ax2.text(0.17, 0.17, "Inner\npoles\nescape", color=CRIT_C,
             fontsize=8, fontweight="bold", va="center")

    ax2.set_xlabel("Re(u)", color=TEXT_C, fontsize=11)
    ax2.set_ylabel("Im(u)", color=TEXT_C, fontsize=11)
    ax2.set_xlim(-0.35, 0.35)
    ax2.set_ylim(-0.35, 0.35)
    ax2.set_aspect("equal")
    ax2.tick_params(colors=SUBTLE_C)
    ax2.spines["bottom"].set_color(GRID_C)
    ax2.spines["left"].set_color(GRID_C)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Legend for right panel
    leg = ax2.legend(handles=[mpatches.Patch(color=ACCENT_C, label="Ramanujan circle"),
                      mpatches.Patch(color=DISR_C, label="Main ring"),
                      mpatches.Patch(color=CRIT_C, label="Escaping poles")],
                     loc="upper left", fontsize=8, facecolor=BG_PANEL,
                     edgecolor=GRID_C, labelcolor=TEXT_C)

    fig.suptitle("Ihara Zeta Pole Structure — Ramanujan vs Non-Ramanujan",
                 color=TEXT_C, fontsize=14, fontweight="bold", y=1.02)

    savefig(fig, "fig13_ihara_pole_structure.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 14: Triple balance — aggregation/expansion/entropy (tensegrity)
# ══════════════════════════════════════════════════════════════════════════════

def fig14_triple_balance():
    """
    Three-panel figure showing the triple balance of tensions:
    (left) Aggregation (cyclic edges, K < 0)
    (middle) Expansion (free edges, K > 0)
    (right) Entropy (decoupled + pre-crystallization)
    With tensegrity analogy.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), facecolor=BG)

    # ── Left: Aggregation ──
    ax = axes[0]
    setup_ax(ax, "Aggregation  (K_λ < 0)")

    # Dense cyclic graph (like SU(3))
    n_nodes = 8
    cx, cy = 0.5, 0.5
    r = 0.25
    angles = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False) + np.pi / 6

    positions = {}
    for i in range(n_nodes):
        x = cx + r * np.cos(angles[i])
        y = cy + r * np.sin(angles[i])
        positions[f"T{i+1}"] = (x, y)

    # Draw cyclic edges (interlocking cycles)
    for i in range(n_nodes):
        j = (i + 1) % n_nodes
        lbl_i, lbl_j = f"T{i+1}", f"T{j+1}"
        p1, p2 = positions[lbl_i], positions[lbl_j]
        draw_curved_edge(ax, p1, p2, CYCLIC_C, arrow=True, connectionstyle="arc3,rad=0.15", lw=2)
        draw_curved_edge(ax, p2, p1, CYCLIC_C, arrow=True, connectionstyle="arc3,rad=-0.15", lw=1, alpha=0.3)

    # Cross connections (representing full adjoint feedback)
    for i in range(0, n_nodes, 2):
        for j in range(i + 2, min(i + 4, n_nodes)):
            if j < n_nodes:
                p1 = positions[f"T{i+1}"]
                p2 = positions[f"T{j+1}"]
                draw_curved_edge(ax, p1, p2, CYCLIC_C, lw=0.8, alpha=0.2, connectionstyle="arc3,rad=0.0")

    # Nodes
    for lbl, pos in positions.items():
        draw_node(ax, pos, lbl, NODE_CYCLIC, 400, 7)

    # Central cycle symbol
    ax.text(cx, cy, "↻", ha="center", va="center", fontsize=32, color=CYCLIC_C, alpha=0.4)

    # Aggregation strength indicator
    ax.text(0.5, 0.05, "Σ|K_λ| = strong\nAll edges cyclic", ha="center", color=CYCLIC_C,
            fontsize=8, style="italic")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # ── Middle: Expansion ──
    ax = axes[1]
    setup_ax(ax, "Expansion  (K_λ > 0)")

    # Core (small cyclic)
    core_pos = [(0.4, 0.7)]
    draw_node(ax, core_pos[0], "J", NODE_CYCLIC, 500, 9)

    # Free extensions going down
    for i in range(4):
        x = 0.2 + i * 0.2
        y = 0.3
        draw_node(ax, (x, y), f"K{i+1}", NODE_FREE, 500, 8)
        # Arrow from core to extension
        draw_curved_edge(ax, core_pos[0], (x, y), FREE_C, arrow=True,
                         connectionstyle="arc3,rad=0.0", lw=2.5)
        # Diverging arrow from extension
        ax.annotate("", xy=(x, y - 0.08), xytext=(x, y),
                    arrowprops=dict(arrowstyle="->,head_length=5,head_width=3",
                                    color=FREE_C, lw=1.5, alpha=0.5))

    # Core cycle
    ax.text(0.4, 0.75, "↻", ha="center", va="center", fontsize=16, color=CYCLIC_C, alpha=0.4)

    # Expansion strength indicator
    ax.text(0.5, 0.05, "ΣK_λ = causal\nIrreversible directions", ha="center", color=FREE_C,
            fontsize=8, style="italic")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # ── Right: Entropy ──
    ax = axes[2]
    setup_ax(ax, "Entropy  (K_λ = 0)")

    # Decoupled strand
    strand_pos = (0.5, 0.5)
    draw_node(ax, strand_pos, "T₀", NODE_DECOUPLE, 800, 11)

    # Radiating dashed lines (going nowhere)
    for angle_deg in np.arange(0, 360, 30):
        t = math.radians(angle_deg)
        x2 = 0.5 + 0.2 * math.cos(t)
        y2 = 0.5 + 0.2 * math.sin(t)
        ax.plot([0.5, x2], [0.5, y2], color=DECOUPLED_C, lw=1.5, ls=":", alpha=0.3)

    # Pre-crystallization: scattered disconnected nodes
    for i in range(6):
        x = 0.15 + i * 0.12
        y = 0.15
        ax.scatter(x, y, s=100, c=DECOUPLED_C, alpha=0.3, zorder=5)
        ax.text(x, y - 0.04, "T", ha="center", color=DECOUPLED_C, fontsize=7,
                alpha=0.4, style="italic")

    # Entropy indicator
    ax.text(0.5, 0.85, "⋯ 0 feedback", ha="center", color=DECOUPLED_C,
            fontsize=9, style="italic")
    ax.text(0.5, 0.02, "dim(ker K) + missing generators", ha="center", color=DECOUPLED_C,
            fontsize=8, style="italic")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # ── Tensegrity analogy at bottom ──
    fig.text(0.5, 0.005, "Tensegrity Analogy:  Cables (tension)  +  Bars (compression)  +  Void space  =  Equilibrium",
             ha="center", color=SUBTLE_C, fontsize=10, style="italic",
             bbox=dict(facecolor=BG_PANEL, edgecolor=GRID_C, alpha=0.5, boxstyle="round,pad=0.3"))

    fig.suptitle("Triple Balance — Three Tensions in the Killing Spectrum",
                 color=TEXT_C, fontsize=14, fontweight="bold", y=1.02)

    savefig(fig, "fig14_triple_balance.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 15: Tempering effect — Killing signal purification
# ══════════════════════════════════════════════════════════════════════════════

def fig15_tempering_effect():
    """
    Two-panel figure showing the tempering effect from RAG Experiment 5:
    (left) Two-phase crystallization (activation + hardening)
    (right) Killing signal purification during withdrawal phase
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), facecolor=BG)

    # ── Left: Two-phase crystallization ──
    setup_ax(ax1, "Two-Phase Crystallization — Activation + Hardening")

    # Epoch numbers
    epochs = np.arange(0, 601, 10)

    # Activation phase (0-200 epochs): transient Killing separation
    activation_mask = epochs <= 200
    activation_k = np.zeros_like(epochs, dtype=float)
    activation_k[activation_mask] = 0.3 + 0.5 * (1 - np.exp(-epochs[activation_mask] / 50))
    activation_k = np.clip(activation_k, 0, 1)

    # Hardening phase (200-600 epochs): permanent basin transition
    hardening_mask = epochs > 200
    hardening_k = np.zeros_like(epochs, dtype=float)
    hardening_k[hardening_mask] = 0.3 + 0.65 * (1 - np.exp(-(epochs[hardening_mask] - 200) / 100))
    hardening_k = np.clip(hardening_k, 0, 1)

    # Full curve: activation before 200, hardening after
    full_k = np.where(hardening_mask, hardening_k, activation_k)

    # Plot
    ax1.plot(epochs, full_k, color=CYCLIC_C, lw=3, alpha=0.8, path_effects=GLOW)

    # Activation region shading
    ax1.axvspan(0, 200, alpha=0.1, color=CRIT_C)
    ax1.text(100, 0.02, "Activation phase\n(~50 epochs to separate)",
             ha="center", color=CRIT_C, fontsize=9, style="italic")

    # Hardening region shading
    ax1.axvspan(200, 500, alpha=0.1, color=ORD_C)
    ax1.text(350, 0.02, "Hardening phase\n(~200 epochs, 4× activation)",
             ha="center", color=ORD_C, fontsize=9, style="italic")

    # Hardening threshold line
    ax1.axhline(y=0.9, color=DISR_C, lw=2, ls="--", alpha=0.6)
    ax1.text(520, 0.92, "Hardening threshold\n(K > 0.9)", color=DISR_C,
             fontsize=8, fontweight="bold", va="bottom")

    # Annotation: transient vs permanent
    ax1.annotate("Transient\nseparation", xy=(100, 0.75), fontsize=8, color=CRIT_C,
                ha="center", style="italic",
                bbox=dict(facecolor=BG_PANEL, edgecolor=CRIT_C, alpha=0.5, boxstyle="round,pad=0.2"))
    ax1.annotate("Permanent\nbasin", xy=(350, 0.95), fontsize=8, color=ORD_C,
                ha="center", style="italic",
                bbox=dict(facecolor=BG_PANEL, edgecolor=ORD_C, alpha=0.5, boxstyle="round,pad=0.2"))

    ax1.set_xlabel("Epochs", color=TEXT_C, fontsize=11)
    ax1.set_ylabel("Killing signal quality (AUC)", color=TEXT_C, fontsize=11)
    ax1.set_xlim(0, 600)
    ax1.set_ylim(0, 1.1)
    ax1.tick_params(colors=SUBTLE_C)
    ax1.spines["bottom"].set_color(GRID_C)
    ax1.spines["left"].set_color(GRID_C)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.grid(True, axis="y", color=GRID_C, alpha=0.3, ls="--")

    # ── Right: Tempering effect ──
    setup_ax(ax2, "Tempering Effect — Killing Signal Purification")

    # Withdrawal phase data (from G_RAG-25)
    # D200: K_rr AUC 0.79 → 0.93
    # D400: K_rr AUC 0.52 → 0.98
    withdrawal_epochs = [200, 300, 400, 500, 600]
    post_dose_auc = [0.79, 0.75, 0.70, 0.65, 0.60]  # declining during withdrawal
    post_withdrawal_auc = [0.93, 0.94, 0.97, 0.98, 0.98]  # improving during withdrawal

    ax1.plot(withdrawal_epochs, post_dose_auc, color=CRIT_C, lw=2, ls="-",
             alpha=0.5, marker="o", markersize=6, zorder=5)
    ax1.plot(withdrawal_epochs, post_withdrawal_auc, color=CYCLIC_C, lw=3,
             ls="-", alpha=0.8, marker="s", markersize=8, zorder=6)

    # Annotations
    ax1.text(200, 0.79, "Dose:\n0.79", ha="right", color=CRIT_C, fontsize=8)
    ax1.text(200, 0.93, "Withdrawal:\n0.93", ha="right", color=CYCLIC_C, fontsize=8)

    ax1.text(400, 0.50, "Dose:\n0.70", ha="right", color=CRIT_C, fontsize=8)
    ax1.text(400, 0.97, "Withdrawal:\n0.98", ha="right", color=CYCLIC_C, fontsize=8)

    # Improvement arrows
    ax1.annotate("", xy=(400, 0.97), xytext=(400, 0.70),
                arrowprops=dict(arrowstyle="->,head_length=3,head_width=2",
                                color=ORD_C, lw=2.5))
    ax1.text(405, 0.83, "+0.28", color=ORD_C, fontsize=9, fontweight="bold",
             va="center")

    # Phase labels
    ax1.axvspan(0, 200, alpha=0.1, color=ACCENT_C)
    ax1.text(100, 1.05, "Dosing phase\n(cycle pressure ON)", ha="center",
             color=ACCENT_C, fontsize=8, style="italic")

    ax1.axvspan(200, 600, alpha=0.1, color=PURPLE_C)
    ax1.text(400, 1.05, "Withdrawal phase\n(cycle pressure OFF)", ha="center",
             color=PURPLE_C, fontsize=8, style="italic")

    ax1.set_xlabel("Withdrawal epoch", color=TEXT_C, fontsize=11)
    ax1.set_ylabel("K_rr AUC (diagnostic quality)", color=TEXT_C, fontsize=11)
    ax1.set_xlim(150, 650)
    ax1.set_ylim(0.3, 1.08)
    ax1.tick_params(colors=SUBTLE_C)
    ax1.spines["bottom"].set_color(GRID_C)
    ax1.spines["left"].set_color(GRID_C)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.grid(True, axis="y", color=GRID_C, alpha=0.3, ls="--")

    # Legend
    ax1.legend(["Post-dose (declining)", "Post-withdrawal (improving)"],
               loc="lower left", fontsize=8, facecolor=BG_PANEL, edgecolor=GRID_C,
               labelcolor=TEXT_C)

    fig.suptitle("Tempering Effect — Two-Phase Crystallization and Signal Purification",
                 color=TEXT_C, fontsize=14, fontweight="bold", y=1.02)

    savefig(fig, "fig15_tempering_effect.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 16: Two-phase crystallization — activation (~50 ep) + hardening (~200 ep)
# ══════════════════════════════════════════════════════════════════════════════

def fig16_two_phase_crystallization():
    """
    Detailed two-phase crystallization diagram showing:
    (left) Activation phase — transient Killing separation
    (right) Hardening phase — permanent basin transition
    With the sharp hardening threshold between D100–D200.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), facecolor=BG)

    # ── Left: Activation phase (detailed) ──
    setup_ax(ax1, "Phase 1: Activation (~50 epochs)")

    # Activation curve: rapid rise then plateau (transient)
    ep1 = np.linspace(0, 200, 300)
    k_activation = 0.3 + 0.4 * (1 - np.exp(-ep1 / 30))
    # Add slight oscillation to show transient nature
    k_activation += 0.05 * np.sin(ep1 * 0.1) * np.exp(-ep1 / 80)

    ax1.plot(ep1, k_activation, color=CRIT_C, lw=3, alpha=0.8, path_effects=GLOW)

    # Activation threshold (K begins to separate edge types)
    ax1.axhline(y=0.7, color=ACCENT_C, lw=1.5, ls="--", alpha=0.5)
    ax1.text(5, 0.72, "K separates edge types", color=ACCENT_C, fontsize=8, style="italic")

    # Collapse point: if cycle pressure removed
    collapse_x = [200, 250]
    collapse_y = [0.75, 0.32]
    ax1.plot(collapse_x, collapse_y, color=DISR_C, lw=3, ls="--", alpha=0.6)
    ax1.annotate("", xy=(250, 0.32), xytext=(230, 0.55),
                arrowprops=dict(arrowstyle="->,head_length=6,head_width=4",
                                color=DISR_C, lw=2))
    ax1.text(260, 0.35, "If cycle pressure\nremoved → collapse!",
             color=DISR_C, fontsize=8, fontweight="bold", va="center")

    # Key point: activation threshold
    ax1.scatter([50], [k_activation[125]], s=200, c=ACCENT_C, zorder=5,
                edgecolors="white", linewidths=2)
    ax1.text(50, 0.15, "50 epochs:\nK starts separating\nedge types",
             ha="center", color=ACCENT_C, fontsize=8, fontweight="bold")

    # Phase label
    ax1.text(100, 0.02, "Transient — edges relax to disorder upon withdrawal",
             ha="center", color=SUBTLE_C, fontsize=8, style="italic")

    ax1.set_xlabel("Epochs", color=TEXT_C, fontsize=11)
    ax1.set_ylabel("Killing signal (K_rr AUC)", color=TEXT_C, fontsize=11)
    ax1.set_xlim(0, 300)
    ax1.set_ylim(0, 1.05)
    ax1.tick_params(colors=SUBTLE_C)
    ax1.spines["bottom"].set_color(GRID_C)
    ax1.spines["left"].set_color(GRID_C)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.grid(True, axis="y", color=GRID_C, alpha=0.3, ls="--")

    # ── Right: Hardening phase (detailed) ──
    setup_ax(ax2, "Phase 2: Hardening (~200 epochs, 4× activation)")

    # Hardening curve: sharp threshold at D100–D200
    ep2 = np.linspace(0, 400, 500)
    # Sharp threshold: collapse below 100, amplification above 200
    k_hardening = np.where(ep2 < 100, 0.32,
                      0.32 + 0.61 * (1 / (1 + np.exp(-2 * (ep2 - 150)))))

    ax2.plot(ep2, k_hardening, color=ORD_C, lw=3, alpha=0.8, path_effects=GLOW)

    # Hardening threshold band (D100–D200)
    ax2.axvspan(100, 200, alpha=0.15, color=ORD_C)
    ax2.text(150, 0.95, "Hardening threshold\nD100–D200", ha="center",
             color=ORD_C, fontsize=9, fontweight="bold", style="italic")

    # Collapse vs amplification
    ax2.axhline(y=0.32, color=DISR_C, lw=1.5, ls=":", alpha=0.4)
    ax2.text(205, 0.30, "Collapse\n(~0.32)", color=DISR_C, fontsize=8, va="top")

    ax2.axhline(y=0.93, color=ORD_C, lw=1.5, ls=":", alpha=0.4)
    ax2.text(300, 0.95, "Amplification\n(0.93+)", color=ORD_C, fontsize=8, va="bottom")

    # Permanent basin annotation
    ax2.scatter([300], [k_hardening[300]], s=200, c=ORD_C, zorder=5,
                edgecolors="white", linewidths=2)
    ax2.text(300, 0.02, "Permanent basin:\npersists under\nKilling-only annealing",
             ha="center", color=ORD_C, fontsize=8, fontweight="bold")

    ax2.set_xlabel("Epochs", color=TEXT_C, fontsize=11)
    ax2.set_ylabel("Killing signal (K_rr AUC)", color=TEXT_C, fontsize=11)
    ax2.set_xlim(0, 400)
    ax2.set_ylim(0, 1.05)
    ax2.tick_params(colors=SUBTLE_C)
    ax2.spines["bottom"].set_color(GRID_C)
    ax2.spines["left"].set_color(GRID_C)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.grid(True, axis="y", color=GRID_C, alpha=0.3, ls="--")

    fig.suptitle("Two-Phase Crystallization — Activation (transient) → Hardening (permanent)",
                 color=TEXT_C, fontsize=14, fontweight="bold", y=1.02)

    savefig(fig, "fig16_two_phase_crystallization.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 17: Spectral sectors — SU(2), SU(3), E-E as graph substructures
# ══════════════════════════════════════════════════════════════════════════════

def fig17_spectral_sectors():
    """
    Three-panel figure showing the three spectral sectors as graph substructures:
    (left) SU(2) — binary cycle
    (middle) SU(3) — ternary cycle
    (right) E-E bridge — cross-sector link
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), facecolor=BG)

    # ── Left: SU(2) binary cycle ──
    ax = axes[0]
    setup_ax(ax, "SU(2) — Binary Cycle")

    cx, cy = 0.5, 0.5
    r = 0.25

    # Three generators forming a triangle
    labels = ["T₁", "T₂", "T₃"]
    positions = {}
    for i, lbl in enumerate(labels):
        angle = math.pi / 2 + i * 2 * math.pi / 3
        x = cx + r * math.cos(angle)
        y = cy + r * math.sin(angle)
        positions[lbl] = (x, y)
        draw_node(ax, (x, y), lbl, NODE_CYCLIC, 700, 10)

    # Cyclic edges
    for i in range(3):
        p1 = positions[labels[i]]
        p2 = positions[labels[(i + 1) % 3]]
        draw_curved_edge(ax, p1, p2, CYCLIC_C, arrow=True, connectionstyle="arc3,rad=0.15", lw=3)
        draw_curved_edge(ax, p2, p1, CYCLIC_C, arrow=True, connectionstyle="arc3,rad=-0.15", lw=2, alpha=0.3)

    # Central cycle
    ax.text(cx, cy, "↻", ha="center", va="center", fontsize=36, color=CYCLIC_C, alpha=0.5)

    # Info box
    info_text = r"$h^\vee_2 = 2$" + "\n" + r"$\mathcal{R}_2 = 1/2$" + "\n" + r"$h_2^{\vee 4} = 16$"
    ax.text(0.5, 0.02, info_text, ha="center", color=CYCLIC_C, fontsize=8,
            family="monospace", fontweight="bold",
            bbox=dict(facecolor=BG_PANEL, edgecolor=CYCLIC_C, alpha=0.5, boxstyle="round,pad=0.3"))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # ── Middle: SU(3) ternary cycle ──
    ax = axes[1]
    setup_ax(ax, "SU(3) — Ternary Cycle")

    cx, cy = 0.5, 0.5

    # 8 generators: hexagonal roots + 2 Cartan
    root_labels = ["λ₁", "λ₂", "λ₄", "λ₅", "λ₆", "λ₇"]
    cartan_labels = ["λ₃", "λ₈"]
    pos = {}
    r_hex = 0.28
    for i, lbl in enumerate(root_labels):
        angle = math.pi / 6 + i * math.pi / 3
        x = cx + r_hex * math.cos(angle)
        y = cy + r_hex * math.sin(angle)
        pos[lbl] = (x, y)
        draw_node(ax, (x, y), lbl, NODE_CYCLIC, 350, 7)

    pos["λ₃"] = (cx - 0.05, cy)
    pos["λ₈"] = (cx + 0.05, cy)
    draw_node(ax, pos["λ₃"], "λ₃", NODE_CORE, 300, 7)
    draw_node(ax, pos["λ₈"], "λ₈", NODE_CORE, 300, 7)

    # Hexagonal edges
    for i in range(6):
        p1 = pos[root_labels[i]]
        p2 = pos[root_labels[(i + 1) % 6]]
        draw_curved_edge(ax, p1, p2, CYCLIC_C, arrow=True, connectionstyle="arc3,rad=0.08", lw=2, alpha=0.8)

    # Cross connections
    for i in range(3):
        p1 = pos[root_labels[i]]
        p2 = pos[root_labels[(i + 3) % 6]]
        draw_curved_edge(ax, p1, p2, CYCLIC_C, lw=1, alpha=0.2, connectionstyle="arc3,rad=0.0")

    # Cartan-to-root
    for lbl in root_labels[:3]:
        draw_curved_edge(ax, pos["λ₃"], pos[lbl], CYCLIC_C, lw=0.8, alpha=0.3, connectionstyle="arc3,rad=0.1")
    for lbl in root_labels[3:]:
        draw_curved_edge(ax, pos["λ₈"], pos[lbl], CYCLIC_C, lw=0.8, alpha=0.3, connectionstyle="arc3,rad=-0.1")

    ax.text(cx, cy, "↻", ha="center", va="center", fontsize=32, color=CYCLIC_C, alpha=0.4)

    # Info box
    info_text = r"$h^\vee_3 = 3$" + "\n" + r"$\mathcal{R}_3 = (\sqrt{57}-3)/(2\sqrt{17})$" + "\n" + r"$h_3^{\vee 4} = 81$"
    ax.text(0.5, 0.02, info_text, ha="center", color=CYCLIC_C, fontsize=8,
            family="monospace", fontweight="bold",
            bbox=dict(facecolor=BG_PANEL, edgecolor=CYCLIC_C, alpha=0.5, boxstyle="round,pad=0.3"))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # ── Right: E-E bridge ──
    ax = axes[2]
    setup_ax(ax, "E-E Bridge — Cross-Sector Link")

    # Two clusters (E₇ and E₆ subgraphs) with bridge
    # E₇ cluster (larger)
    e7_nodes = ["E₇₁", "E₇₂", "E₇₃", "E₇₄", "E₇₅", "E₇₆", "E₇₇"]
    e7_pos = {}
    cx7, cy7 = 0.3, 0.65
    r7 = 0.15
    for i, lbl in enumerate(e7_nodes):
        angle = math.pi / 2 + i * 2 * math.pi / 7
        x = cx7 + r7 * math.cos(angle)
        y = cy7 + r7 * math.sin(angle)
        e7_pos[lbl] = (x, y)
        draw_node(ax, (x, y), lbl, PURPLE_C, 250, 6)

    # E₆ cluster (smaller)
    e6_nodes = ["E₆₁", "E₆₂", "E₆₃", "E₆₄", "E₆₅", "E₆₆"]
    e6_pos = {}
    cx6, cy6 = 0.7, 0.35
    r6 = 0.12
    for i, lbl in enumerate(e6_nodes):
        angle = math.pi / 2 + i * 2 * math.pi / 6
        x = cx6 + r6 * math.cos(angle)
        y = cy6 + r6 * math.sin(angle)
        e6_pos[lbl] = (x, y)
        draw_node(ax, (x, y), PURPLE_C, 250, 6)
        draw_node(ax, (x, y), lbl, PURPLE_C, 250, 6)

    # Bridge connections
    bridge_pairs = [(list(e7_pos.values())[0], list(e6_pos.values())[0]),
                    (list(e7_pos.values())[2], list(e6_pos.values())[2])]
    for p1, p2 in bridge_pairs:
        draw_curved_edge(ax, p1, p2, ACCENT_C, arrow=True, connectionstyle="arc3,rad=0.3", lw=3)
        draw_curved_edge(ax, p2, p1, ACCENT_C, arrow=True, connectionstyle="arc3,rad=-0.3", lw=3)

    # Bridge label
    mid_x = (cx7 + cx6) / 2
    mid_y = (cy7 + cy6) / 2
    ax.text(mid_x, mid_y, r"$R_{EE}$", ha="center", va="center", fontsize=12,
            color=ACCENT_C, fontweight="bold",
            bbox=dict(facecolor=BG_PANEL, edgecolor=ACCENT_C, alpha=0.7, boxstyle="round,pad=0.2"))

    # Info box
    info_text = r"$R_{EE} = 1/\sqrt{2}$" + "\n" + "Cross-sector coupling" + "\n" + "Generation mixing"
    ax.text(0.5, 0.02, info_text, ha="center", color=ACCENT_C, fontsize=8,
            family="monospace", fontweight="bold",
            bbox=dict(facecolor=BG_PANEL, edgecolor=ACCENT_C, alpha=0.5, boxstyle="round,pad=0.3"))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # ── Connecting arrows between sectors ──
    # SU(2) → SU(3)
    fig.text(0.25, 0.96, "→", ha="center", va="center", fontsize=20, color=SUBTLE_C)
    # SU(3) → E-E
    fig.text(0.58, 0.96, "→", ha="center", va="center", fontsize=20, color=SUBTLE_C)

    # Title
    fig.suptitle("Spectral Sectors — Three Graph Substructures of the Evolution Graph",
                 color=TEXT_C, fontsize=14, fontweight="bold", y=1.02)

    savefig(fig, "fig17_spectral_sectors.png")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 18: Bridge diagram — Paper 1 → Paper 2 → Paper 3
# ══════════════════════════════════════════════════════════════════════════════

def fig18_bridge_diagram():
    """
    Three-stage pipeline diagram showing the Paper 1 → Paper 2 → Paper 3 sequence:
    emergence → topology → observables
    """
    fig, ax = plt.subplots(1, 1, figsize=(16, 8), facecolor=BG)
    setup_ax(ax)
    ax.set_title("SS-PEG: Three Papers — From Emergence to Observable Physics",
                 color=TEXT_C, fontsize=16, fontweight="bold", pad=20)

    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

    # ── Stage 1: Paper 1 (Emergence) ──
    x1, y1 = 0.5, 0.45
    w1, h1 = 3.5, 4.5
    box1 = FancyBboxPatch((x1, y1), w1, h1, boxstyle="round,pad=0.15",
                           facecolor="#1a2a3a", edgecolor=ACCENT_C, lw=3, alpha=0.8)
    ax.add_patch(box1)

    # Title
    ax.text(x1 + w1/2, y1 + h1 - 0.3, "PAPER 1", ha="center", va="center",
            color=ACCENT_C, fontsize=16, fontweight="bold")
    ax.text(x1 + w1/2, y1 + h1 - 0.7, "Lie Algebra Emergence", ha="center", va="center",
            color=TEXT_C, fontsize=12, style="italic")

    # Key results
    results1 = [
        "51 core tests, 100% success",
        "6 families: SU, SO, Sp, G₂, F₄, E₆",
        "Non-compact: sl(2,R), so(3,1), su(2,2)",
        "Ablation: T_K only → algebra emerges",
        "Density transition: ρ=1 step function",
        "h^∨ = I_R·|K|/κ² exact for 15 algebras",
    ]
    for i, r in enumerate(results1):
        ax.text(x1 + 0.2, y1 + h1 - 1.2 - i * 0.45, f"✓ {r}",
                ha="left", va="center", color=TEXT_C, fontsize=9, family="monospace")

    # Output
    ax.text(x1 + w1/2, y1 - 0.3, "Output:", ha="center", va="top",
            color=SUBTLE_C, fontsize=9, style="italic")
    output1 = "Emergent algebras\n{K_ab, f_abc, h^∨, C_λ}"
    ax.text(x1 + w1/2, y1 - 0.6, output1, ha="center", va="top",
            color=ACCENT_C, fontsize=9, family="monospace", fontweight="bold",
            bbox=dict(facecolor=BG_PANEL, edgecolor=ACCENT_C, alpha=0.5, boxstyle="round,pad=0.3"))

    # ── Stage 2: Paper 2 (Topology) ──
    x2, y2 = 5.5, 0.45
    box2 = FancyBboxPatch((x2, y2), w1, h1, boxstyle="round,pad=0.15",
                           facecolor="#1a3a2a", edgecolor=CYCLIC_C, lw=3, alpha=0.8)
    ax.add_patch(box2)

    ax.text(x2 + w1/2, y2 + h1 - 0.3, "PAPER 2", ha="center", va="center",
            color=CYCLIC_C, fontsize=16, fontweight="bold")
    ax.text(x2 + w1/2, y2 + h1 - 0.7, "Graph Topology & Spectral Analysis", ha="center", va="center",
            color=TEXT_C, fontsize=12, style="italic")

    results2 = [
        "5 graph archetypes (derived, not postulated)",
        "19/37 algebras Ramanujan, h* ≈ 9.5",
        "Ihara zeta: pole structure, Sunada failure",
        "Triple balance: aggregation/expansion/entropy",
        "Tempering: activation (50 ep) + hardening (200 ep)",
        "Killing: diagnostic, not generative",
    ]
    for i, r in enumerate(results2):
        ax.text(x2 + 0.2, y2 + h1 - 1.2 - i * 0.45, f"✓ {r}",
                ha="left", va="center", color=TEXT_C, fontsize=9, family="monospace")

    ax.text(x2 + w1/2, y2 - 0.3, "Output:", ha="center", va="top",
            color=SUBTLE_C, fontsize=9, style="italic")
    output2 = "5 spectral invariants\n{R₂, R₃, R_EE, h³, h²}"
    ax.text(x2 + w1/2, y2 - 0.6, output2, ha="center", va="top",
            color=CYCLIC_C, fontsize=9, family="monospace", fontweight="bold",
            bbox=dict(facecolor=BG_PANEL, edgecolor=CYCLIC_C, alpha=0.5, boxstyle="round,pad=0.3"))

    # ── Stage 3: Paper 3 (Observables) ──
    x3, y3 = 10.5, 0.45
    box3 = FancyBboxPatch((x3, y3), w1, h1, boxstyle="round,pad=0.15",
                           facecolor="#2a1a3a", edgecolor=PURPLE_C, lw=3, alpha=0.8)
    ax.add_patch(box3)

    ax.text(x3 + w1/2, y3 + h1 - 0.3, "PAPER 3", ha="center", va="center",
            color=PURPLE_C, fontsize=16, fontweight="bold")
    ax.text(x3 + w1/2, y3 + h1 - 0.7, "Standard Model from Graph Topology", ha="center", va="center",
            color=TEXT_C, fontsize=12, style="italic")

    results3 = [
        "73 predictions tested",
        "72/73 within 1% error",
        "5 building blocks → α, α_s, α₂",
        "Fermion masses on 3D lattice",
        "CKM + PMNS from graph topology",
        "0 free parameters",
    ]
    for i, r in enumerate(results3):
        ax.text(x3 + 0.2, y3 + h1 - 1.2 - i * 0.45, f"✓ {r}",
                ha="left", va="center", color=TEXT_C, fontsize=9, family="monospace")

    ax.text(x3 + w1/2, y3 - 0.3, "Output:", ha="center", va="top",
            color=SUBTLE_C, fontsize=9, style="italic")
    output3 = "Observable physics\n{α, m_i, θ_ij, δ_CP}"
    ax.text(x3 + w1/2, y3 - 0.6, output3, ha="center", va="top",
            color=PURPLE_C, fontsize=9, family="monospace", fontweight="bold",
            bbox=dict(facecolor=BG_PANEL, edgecolor=PURPLE_C, alpha=0.5, boxstyle="round,pad=0.3"))

    # ── Arrows between stages ──
    # Paper 1 → Paper 2
    arrow1 = FancyArrowPatch((x1 + w1, y1 + h1/2), (x2, y1 + h1/2),
                              arrowstyle="->,head_length=12,head_width=8",
                              color=TEXT_C, lw=3, mutation_scale=30)
    ax.add_patch(arrow1)
    ax.text((x1 + w1 + x2) / 2, y1 + h1/2 + 0.3, "spectral\ninvariants",
            ha="center", va="bottom", color=SUBTLE_C, fontsize=9, style="italic")

    # Paper 2 → Paper 3
    arrow2 = FancyArrowPatch((x2 + w1, y1 + h1/2), (x3, y1 + h1/2),
                              arrowstyle="->,head_length=12,head_width=8",
                              color=TEXT_C, lw=3, mutation_scale=30)
    ax.add_patch(arrow2)
    ax.text((x2 + w1 + x3) / 2, y1 + h1/2 + 0.3, "5 building\nblocks",
            ha="center", va="bottom", color=SUBTLE_C, fontsize=9, style="italic")

    # ── Bottom: narrative ──
    ax.text(0.5, 0.02, "Sequence:  Emergence (P1)  →  Topology (P2)  →  Observables (P3)",
            ha="center", va="center", color=TEXT_C, fontsize=12, fontweight="bold",
            style="italic",
            bbox=dict(facecolor=BG_PANEL, edgecolor=GRID_C, alpha=0.5, boxstyle="round,pad=0.4"))

    ax.set_xlim(0, 15)
    ax.set_ylim(-0.1, 5.5)
    ax.set_aspect("equal")

    savefig(fig, "fig18_bridge_diagram.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 19: Ramanujan × Ramanujan → Ramanujan (composability)
# ══════════════════════════════════════════════════════════════════════════════

def fig19_composability():
    """
    Diagram showing product graph composability:
    Ramanujan × Ramanujan → Ramanujan
    With SM product and SU(5) GUT comparison.
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 8), facecolor=BG)
    setup_ax(ax, "Product Graph Composability — Ramanujan × Ramanujan → Ramanujan")

    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

    # ── Left: SU(3) Ramanujan ──
    box_su3 = FancyBboxPatch((0.1, 0.55), 0.35, 0.35, boxstyle="round,pad=0.02",
                              facecolor="#1a3a2a", edgecolor=CYCLIC_C, lw=2, alpha=0.8)
    ax.add_patch(box_su3)
    ax.text(0.275, 0.88, "SU(3)", ha="center", va="center",
            color=CYCLIC_C, fontsize=12, fontweight="bold")
    ax.text(0.275, 0.75, r"$\mathcal{R}_3 = 0.552$", ha="center", va="center",
            color=TEXT_C, fontsize=9, family="monospace")
    ax.text(0.275, 0.65, "✓ Ramanujan", ha="center", va="center",
            color=ORD_C, fontsize=10, fontweight="bold")

    # ── Left: SU(2) Ramanujan ──
    box_su2 = FancyBboxPatch((0.1, 0.15), 0.35, 0.35, boxstyle="round,pad=0.02",
                              facecolor="#1a3a2a", edgecolor=CYCLIC_C, lw=2, alpha=0.8)
    ax.add_patch(box_su2)
    ax.text(0.275, 0.48, "SU(2)", ha="center", va="center",
            color=CYCLIC_C, fontsize=12, fontweight="bold")
    ax.text(0.275, 0.35, r"$\mathcal{R}_2 = 0.500$", ha="center", va="center",
            color=TEXT_C, fontsize=9, family="monospace")
    ax.text(0.275, 0.25, "✓ Ramanujan", ha="center", va="center",
            color=ORD_C, fontsize=10, fontweight="bold")

    # ── Multiplication symbol ──
    ax.text(0.5, 0.55, "⊕", ha="center", va="center", fontsize=36, color=ACCENT_C)

    # ── Arrow → ──
    arrow = FancyArrowPatch((0.55, 0.55), (0.7, 0.55),
                            arrowstyle="->,head_length=10,head_width=8",
                            color=TEXT_C, lw=3, mutation_scale=30)
    ax.add_patch(arrow)

    # ── Center: SM product ──
    box_sm = FancyBboxPatch((0.72, 0.35), 0.3, 0.4, boxstyle="round,pad=0.02",
                             facecolor="#1a2a4a", edgecolor=ACCENT_C, lw=3, alpha=0.8)
    ax.add_patch(box_sm)
    ax.text(0.87, 0.7, "SM", ha="center", va="center",
            color=ACCENT_C, fontsize=14, fontweight="bold")
    ax.text(0.87, 0.6, r"su(3) ⊕ su(2)", ha="center", va="center",
            color=TEXT_C, fontsize=9, family="monospace")
    ax.text(0.87, 0.45, r"$R_{prod} = 0.657$", ha="center", va="center",
            color=TEXT_C, fontsize=9, family="monospace")
    ax.text(0.87, 0.3, "✓ Ramanujan", ha="center", va="center",
            color=ORD_C, fontsize=10, fontweight="bold")

    # ── Right: SU(5) GUT comparison ──
    box_gut = FancyBboxPatch((0.95, 0.55), 0.35, 0.35, boxstyle="round,pad=0.02",
                              facecolor="#2a2a1a", edgecolor=CRIT_C, lw=2, alpha=0.5)
    ax.add_patch(box_gut)
    ax.text(1.125, 0.88, "SU(5)", ha="center", va="center",
            color=CRIT_C, fontsize=12, fontweight="bold")
    ax.text(1.125, 0.75, r"$\mathcal{R} = 0.653$", ha="center", va="center",
            color=TEXT_C, fontsize=9, family="monospace")
    ax.text(1.125, 0.65, "✓ Ramanujan", ha="center", va="center",
            color=ORD_C, fontsize=10, fontweight="bold")

    # ── Arrow down from SU(5) ──
    arrow2 = FancyArrowPatch((1.125, 0.55), (1.125, 0.45),
                              arrowstyle="->,head_length=8,head_width=6",
                              color=SUBTLE_C, lw=2, mutation_scale=20)
    ax.add_patch(arrow2)

    # ── SU(5) match label ──
    ax.text(1.125, 0.38, "R ≈ SM (0.653 vs 0.657)", ha="center", va="center",
            color=SUBTLE_C, fontsize=8, style="italic",
            bbox=dict(facecolor=BG_PANEL, edgecolor=SUBTLE_C, alpha=0.3, boxstyle="round,pad=0.2"))

    # ── Bottom annotation ──
    ax.text(0.5, 0.08, "P7: Ramanujan × Ramanujan → Ramanujan (unanimous, 6/6 definitions)",
            ha="center", color=ORD_C, fontsize=10, fontweight="bold", style="italic",
            bbox=dict(facecolor=BG_PANEL, edgecolor=ORD_C, alpha=0.5, boxstyle="round,pad=0.3"))

    ax.text(0.5, 0.02, "SM product is spectrally optimal — no penalty for composition",
            ha="center", color=SUBTLE_C, fontsize=9, style="italic")

    ax.set_xlim(0, 1.5)
    ax.set_ylim(-0.1, 1.0)
    ax.set_aspect("equal")

    savefig(fig, "fig19_composability.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 20: Five emergence minima — generator count × parameterization × Killing target
# ══════════════════════════════════════════════════════════════════════════════

def fig20_emergence_minima():
    """
    3D-like schematic showing how three emergence minima determine graph topology:
    1. Generator count (n_gen)
    2. Parameterization constraint
    3. Killing target signature
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 9), facecolor=BG)
    setup_ax(ax)
    ax.set_title("The Three Emergence Minima — What Determines Graph Topology",
                 color=TEXT_C, fontsize=15, fontweight="bold", pad=15)

    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

    # ── Three input boxes (top) ──
    inputs = [
        (0.15, 0.82, "① Generator Count\nn_gen", ACCENT_C, "Sets max graph dimension"),
        (0.5, 0.82, "② Parameterization\nConstraint", CYCLIC_C, "Selects edge type available"),
        (0.85, 0.82, "③ Killing Target\nSignature", PURPLE_C, "Selects specific archetype"),
    ]

    for x, y, label, color, desc in inputs:
        box = FancyBboxPatch((x - 0.12, y - 0.05), 0.24, 0.12, boxstyle="round,pad=0.01",
                              facecolor=BG_PANEL, edgecolor=color, lw=2.5, alpha=0.8)
        ax.add_patch(box)
        ax.text(x, y + 0.02, label, ha="center", va="center",
                color=color, fontsize=10, fontweight="bold")
        ax.text(x, y - 0.08, desc, ha="center", va="center",
                color=SUBTLE_C, fontsize=8, style="italic")

    # ── Arrow down ──
    ax.annotate("", xy=(0.5, 0.7), xytext=(0.5, 0.75),
                arrowprops=dict(arrowstyle="->,head_length=8,head_width=6",
                                color=TEXT_C, lw=2))
    ax.text(0.5, 0.72, "Killing metric\nminimization", ha="center", color=SUBTLE_C,
            fontsize=8, style="italic")

    # ── Variational principle box ──
    vp_box = FancyBboxPatch((0.3, 0.55), 0.4, 0.12, boxstyle="round,pad=0.01",
                             facecolor="#2a1a3a", edgecolor=ACCENT_C, lw=3, alpha=0.8)
    ax.add_patch(vp_box)
    ax.text(0.5, 0.61, "Variational Principle", ha="center", va="center",
            color=ACCENT_C, fontsize=11, fontweight="bold")
    ax.text(0.5, 0.57, r"$\min T_K = \|K - K_{target}\|^2$", ha="center", va="center",
            color=TEXT_C, fontsize=8, family="monospace")

    # ── Arrow down ──
    ax.annotate("", xy=(0.5, 0.5), xytext=(0.5, 0.55),
                arrowprops=dict(arrowstyle="->,head_length=8,head_width=6",
                                color=TEXT_C, lw=2))

    # ── Output: 5 archetypes ──
    archetypes = [
        (0.15, 0.38, "I: Pure Cycle\nSU(N), SO(N), Sp(2N)", CYCLIC_C, "All cyclic"),
        (0.38, 0.38, "II: Mixed\nso(3,1), su(2,2)", ACCENT_C, "Cyclic + free"),
        (0.62, 0.38, "III: Cycle+Strand\nU(2), U(3)", PURPLE_C, "Cyclic + decoupled"),
        (0.85, 0.38, "IV/V: Free/Frustrated\n(theoretical)", DECOUPLED_C, "Open or degenerate"),
    ]

    for x, y, label, color, desc in archetypes:
        box = FancyBboxPatch((x - 0.1, y - 0.04), 0.2, 0.08, boxstyle="round,pad=0.01",
                              facecolor=BG_PANEL, edgecolor=color, lw=2, alpha=0.6)
        ax.add_patch(box)
        ax.text(x, y, label, ha="center", va="center",
                color=color, fontsize=7, fontweight="bold")
        ax.text(x, y - 0.06, desc, ha="center", va="center",
                color=SUBTLE_C, fontsize=6, style="italic")

    # ── Bottom: key insight ──
    insight = "The graph topology is NOT a free parameter — it is uniquely determined by these three inputs."
    ax.text(0.5, 0.18, insight, ha="center", color=TEXT_C, fontsize=10, fontweight="bold",
            style="italic",
            bbox=dict(facecolor=BG_PANEL, edgecolor=ACCENT_C, alpha=0.5, boxstyle="round,pad=0.3"))

    # ── Table of emergence conditions ──
    table_y = 0.02
    table_data = [
        ["Archetype", "n_gen", "Parameterization", "K_target", "Condition"],
        ["I (pure)", "dim(g)", "Hermitian/antisym.", "-I", "ρ = 1"],
        ["II (mixed)", "dim(g)", "η-preserving", "Indefinite", "η + indef."],
        ["III (cycle+)", "d²", "Hermitian (trace)", "-I", "Trace DOF"],
        ["IV (free)", "—", "—", "+I", "Not observed"],
        ["V (frustrated)", "dim(g)", "η-preserving", "-I (conflict)", "Constraint × target"],
    ]

    for i, row in enumerate(table_data):
        y = table_y - i * 0.018
        if i == 0:  # Header
            for j, cell in enumerate(row):
                ax.text(0.15 + j * 0.17, y, cell, ha="left", va="center",
                        color=ACCENT_C, fontsize=7, fontweight="bold", family="monospace")
        else:
            for j, cell in enumerate(row):
                ax.text(0.15 + j * 0.17, y, cell, ha="left", va="center",
                        color=TEXT_C, fontsize=7, family="monospace")

    ax.set_xlim(0, 1.1)
    ax.set_ylim(-0.1, 1.0)
    ax.set_aspect("equal")

    savefig(fig, "fig20_emergence_minima.png")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

ALL_FIGURES = {
    11: ("Ramanujan Classification — 37 algebras", fig11_ramanujan_classification),
    12: ("Critical Scaling + Pole Migration Velocity", fig12_critical_scaling),
    13: ("Ihara Zeta Pole Structure — Ramanujan vs Non-Ramanujan", fig13_ihara_pole_structure),
    14: ("Triple Balance — Aggregation/Expansion/Entropy", fig14_triple_balance),
    15: ("Tempering Effect — Two-Phase Crystallization", fig15_tempering_effect),
    16: ("Two-Phase Crystallization — Activation + Hardening", fig16_two_phase_crystallization),
    17: ("Spectral Sectors — SU(2), SU(3), E-E as Graph Substructures", fig17_spectral_sectors),
    18: ("Bridge Diagram — Paper 1 → Paper 2 → Paper 3", fig18_bridge_diagram),
    19: ("Ramanujan × Ramanujan → Ramanujan (Composability)", fig19_composability),
    20: ("Five Emergence Minima — Generator Count × Parameterization × Target", fig20_emergence_minima),
}


def main():
    parser = argparse.ArgumentParser(description="Generate Paper 2 graph figures (new)")
    parser.add_argument("--fig", type=int, nargs="*",
                        help="Figure number(s) to generate (default: all)")
    args = parser.parse_args()

    figs = args.fig if args.fig else sorted(ALL_FIGURES.keys())

    print(f"Generating {len(figs)} figure(s) → {OUT_DIR}/")
    print()
    for n in figs:
        if n not in ALL_FIGURES:
            print(f"  ✗ Figure {n} not found (valid: {sorted(ALL_FIGURES.keys())})")
            continue
        desc, fn = ALL_FIGURES[n]
        print(f"  [{n:2d}] {desc}")
        fn()

    print(f"\nDone. {len(figs)} figures saved to {OUT_DIR}/")


if __name__ == "__main__":
    main()
