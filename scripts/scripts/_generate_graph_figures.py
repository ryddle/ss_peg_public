#!/usr/bin/env python3
"""
Generate publication-quality graph figures for GRAPH_TOPOLOGY_FROM_KILLING.md

Each figure illustrates a graph archetype or concept from the SS-PEG framework,
showing how Killing eigenvalue spectra determine graph topology.

Usage:
    python _generate_graph_figures.py          # Generate all figures
    python _generate_graph_figures.py --fig 3  # Generate only figure 3
"""

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Arc, Circle
from matplotlib.collections import LineCollection
import matplotlib.patheffects as pe
import numpy as np

# ── Global style ──────────────────────────────────────────────────────────────

BG       = "#0d1117"
BG_PANEL = "#161b22"
GRID_C   = "#21262d"
TEXT_C   = "#e6edf3"
SUBTLE_C = "#8b949e"

# Edge class colours
CYCLIC_C    = "#58d68d"   # green-teal  (compact / cyclic)
FREE_C      = "#e74c3c"   # red-orange  (non-compact / free)
DECOUPLED_C = "#7f8c8d"   # grey        (abelian / decoupled)
ACCENT_C    = "#3498db"   # blue accent
PURPLE_C    = "#a855f7"   # purple for special

# Node colours
NODE_CYCLIC  = "#2ecc71"
NODE_FREE    = "#e67e22"
NODE_DECOUPLE= "#95a5a6"
NODE_CORE    = "#1abc9c"

GLOW = [pe.withStroke(linewidth=5, foreground=BG)]

OUT_DIR = Path(__file__).parent / "figures"


def setup_ax(ax, title=None):
    """Dark background, no axes, optional title."""
    ax.set_facecolor(BG_PANEL)
    ax.set_aspect("equal")
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
    kw = dict(
        color=color, lw=lw, linestyle=style, alpha=alpha,
        connectionstyle=connectionstyle,
        path_effects=GLOW if alpha > 0.7 else [],
    )
    if arrow:
        kw["arrowstyle"] = "->,head_length=8,head_width=5"
    else:
        kw["arrowstyle"] = "-"
    patch = FancyArrowPatch(p1, p2, **kw)
    ax.add_patch(patch)
    return patch


def draw_self_loop(ax, centre, radius=0.18, angle=90, color=CYCLIC_C, lw=2):
    """Draw a small self-loop arc at a node."""
    theta1, theta2 = angle - 120, angle + 120
    arc = Arc(centre, radius * 2, radius * 2, angle=0,
              theta1=theta1, theta2=theta2,
              color=color, lw=lw, path_effects=GLOW)
    ax.add_patch(arc)
    # arrowhead at the end
    t_end = math.radians(theta2)
    x_end = centre[0] + radius * math.cos(t_end)
    y_end = centre[1] + radius * math.sin(t_end)
    t_prev = math.radians(theta2 - 5)
    x_prev = centre[0] + radius * math.cos(t_prev)
    y_prev = centre[1] + radius * math.sin(t_prev)
    ax.annotate("", xy=(x_end, y_end), xytext=(x_prev, y_prev),
                arrowprops=dict(arrowstyle="->,head_length=5,head_width=3",
                                color=color, lw=lw))


def draw_node(ax, pos, label, color=NODE_CYCLIC, size=600, fontsize=9, text_color=BG):
    """Draw a glowing node with label."""
    ax.scatter(*pos, s=size, c=color, zorder=5, edgecolors="white", linewidths=0.8,
               alpha=0.95)
    # Glow effect
    ax.scatter(*pos, s=size * 2.5, c=color, zorder=4, alpha=0.12)
    ax.text(pos[0], pos[1], label, ha="center", va="center",
            fontsize=fontsize, color=text_color, fontweight="bold", zorder=6)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Three Edge Classes
# ══════════════════════════════════════════════════════════════════════════════

def fig1_edge_classes():
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), facecolor=BG)

    # ── Cyclic ──
    ax = axes[0]
    setup_ax(ax, "Cyclic Edge  (K_λ < 0)")
    # Two nodes with a curved bidirectional loop
    p1, p2 = (0.3, 0.5), (0.7, 0.5)
    draw_node(ax, p1, "Tₐ", NODE_CYCLIC, 800, 11)
    draw_node(ax, p2, "T_b", NODE_CYCLIC, 800, 11)
    draw_curved_edge(ax, p1, p2, CYCLIC_C, arrow=True, connectionstyle="arc3,rad=0.25")
    draw_curved_edge(ax, p2, p1, CYCLIC_C, arrow=True, connectionstyle="arc3,rad=0.25")
    draw_self_loop(ax, p1, 0.12, 150, CYCLIC_C)
    ax.text(0.5, 0.15, "Adjoint feedback loop\nCompact direction", ha="center",
            color=CYCLIC_C, fontsize=9, style="italic")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.0)

    # ── Free ──
    ax = axes[1]
    setup_ax(ax, "Free Edge  (K_λ > 0)")
    p1, p2 = (0.2, 0.5), (0.8, 0.5)
    draw_node(ax, p1, "Tₐ", NODE_FREE, 800, 11)
    draw_node(ax, p2, "T_b", NODE_FREE, 800, 11)
    draw_curved_edge(ax, p1, p2, FREE_C, arrow=True, connectionstyle="arc3,rad=0.0", lw=3)
    # small diverging arrows
    for dy in [0.12, -0.12]:
        ax.annotate("", xy=(0.95, 0.5 + dy), xytext=(0.8, 0.5),
                    arrowprops=dict(arrowstyle="->,head_length=6,head_width=4",
                                    color=FREE_C, lw=1.5, alpha=0.5))
    ax.text(0.5, 0.15, "No return path\nNon-compact direction", ha="center",
            color=FREE_C, fontsize=9, style="italic")
    ax.set_xlim(-0.05, 1.1)
    ax.set_ylim(-0.05, 1.0)

    # ── Decoupled ──
    ax = axes[2]
    setup_ax(ax, "Decoupled Edge  (K_λ = 0)")
    p1 = (0.5, 0.5)
    draw_node(ax, p1, "T₀", NODE_DECOUPLE, 800, 11)
    # dashed lines going nowhere
    for angle_deg in [0, 60, 120, 180, 240, 300]:
        t = math.radians(angle_deg)
        x2 = 0.5 + 0.3 * math.cos(t)
        y2 = 0.5 + 0.3 * math.sin(t)
        ax.plot([0.5, x2], [0.5, y2], color=DECOUPLED_C, lw=1.5, ls="--", alpha=0.4)
    ax.text(0.5, 0.08, "Zero adjoint feedback\nAbelian direction", ha="center",
            color=DECOUPLED_C, fontsize=9, style="italic")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.1, 1.0)

    fig.suptitle("The Three Edge Classes — Killing Eigenvalue Classification",
                 color=TEXT_C, fontsize=16, fontweight="bold", y=1.02)
    savefig(fig, "fig1_edge_classes.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Archetype I — Pure Cycle (SU(2) and SU(3))
# ══════════════════════════════════════════════════════════════════════════════

def fig2_pure_cycle():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), facecolor=BG)

    # ── SU(2): 3 generators, triangle ──
    setup_ax(ax1, "su(2) — Minimal Cycle  (3 generators)")
    cx, cy = 0.5, 0.45
    r = 0.28
    positions = {}
    labels = ["J₊", "J₃", "J₋"]
    for i, lbl in enumerate(labels):
        angle = math.pi / 2 + i * 2 * math.pi / 3
        x = cx + r * math.cos(angle)
        y = cy + r * math.sin(angle)
        positions[lbl] = (x, y)
        draw_node(ax1, (x, y), lbl, NODE_CYCLIC, 900, 11)

    # Edges forming a cycle
    keys = list(positions.keys())
    for i in range(3):
        p1 = positions[keys[i]]
        p2 = positions[keys[(i + 1) % 3]]
        draw_curved_edge(ax1, p1, p2, CYCLIC_C, arrow=True, connectionstyle="arc3,rad=0.15")
        draw_curved_edge(ax1, p2, p1, CYCLIC_C, arrow=True, connectionstyle="arc3,rad=0.15",
                         alpha=0.4, lw=1.5)

    # Central cycle symbol
    ax1.text(cx, cy, "↻", ha="center", va="center", fontsize=28, color=CYCLIC_C, alpha=0.6)
    ax1.text(cx, 0.02, "[Jᵢ, Jⱼ] = iεᵢⱼₖ Jₖ\nK = diag(−8, −8, −8)",
             ha="center", color=SUBTLE_C, fontsize=8, family="monospace")
    ax1.set_xlim(0.0, 1.0)
    ax1.set_ylim(-0.08, 0.9)

    # ── SU(3): 8 generators ──
    setup_ax(ax2, "su(3) — Octet Cycle  (8 generators)")
    cx, cy = 0.5, 0.48
    # Root vectors form a hexagon, Cartan generators at center
    root_labels = ["λ₁", "λ₂", "λ₄", "λ₅", "λ₆", "λ₇"]
    cartan_labels = ["λ₃", "λ₈"]
    pos = {}
    r_hex = 0.32
    for i, lbl in enumerate(root_labels):
        angle = math.pi / 6 + i * math.pi / 3
        x = cx + r_hex * math.cos(angle)
        y = cy + r_hex * math.sin(angle)
        pos[lbl] = (x, y)
        draw_node(ax2, (x, y), lbl, NODE_CYCLIC, 500, 8)

    # Cartan at centre
    pos["λ₃"] = (cx - 0.06, cy)
    pos["λ₈"] = (cx + 0.06, cy)
    draw_node(ax2, pos["λ₃"], "λ₃", "#1abc9c", 450, 8)
    draw_node(ax2, pos["λ₈"], "λ₈", "#1abc9c", 450, 8)

    # Hexagonal edges (root connections)
    for i in range(6):
        p1 = pos[root_labels[i]]
        p2 = pos[root_labels[(i + 1) % 6]]
        draw_curved_edge(ax2, p1, p2, CYCLIC_C, arrow=True, connectionstyle="arc3,rad=0.08",
                         lw=1.8, alpha=0.8)

    # Cartan-to-root connections (some)
    for lbl in ["λ₁", "λ₂", "λ₄", "λ₅"]:
        draw_curved_edge(ax2, pos["λ₃"], pos[lbl], CYCLIC_C, lw=1.0, alpha=0.35,
                         connectionstyle="arc3,rad=0.12")
    for lbl in ["λ₆", "λ₇", "λ₁", "λ₂"]:
        draw_curved_edge(ax2, pos["λ₈"], pos[lbl], CYCLIC_C, lw=1.0, alpha=0.35,
                         connectionstyle="arc3,rad=-0.12")

    # Diagonal cross connections
    for i in range(3):
        p1 = pos[root_labels[i]]
        p2 = pos[root_labels[(i + 3) % 6]]
        draw_curved_edge(ax2, p1, p2, CYCLIC_C, lw=1.0, alpha=0.25,
                         connectionstyle="arc3,rad=0.0")

    ax2.text(cx, 0.02, "K = diag(−24, …, −24)\nAll 8 edges cyclic — pure gauge",
             ha="center", color=SUBTLE_C, fontsize=8, family="monospace")
    ax2.set_xlim(0.0, 1.0)
    ax2.set_ylim(-0.08, 0.92)

    fig.suptitle("Archetype I: Pure Cycle Graph — K ∝ −I  (Compact Simple Algebras)",
                 color=TEXT_C, fontsize=15, fontweight="bold", y=0.98)
    savefig(fig, "fig2_archetype_I_pure_cycle.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: Archetype II — Mixed Graph (so(3,1) Lorentz)
# ══════════════════════════════════════════════════════════════════════════════

def fig3_mixed_lorentz():
    fig, ax = plt.subplots(1, 1, figsize=(9, 8), facecolor=BG)
    setup_ax(ax, "Archetype II: Mixed Graph — so(3,1) Lorentz Algebra")

    cx, cy_rot = 0.5, 0.65
    r_rot = 0.22
    # Rotation core (cyclic triangle)
    rot_labels = ["J₁", "J₂", "J₃"]
    rot_pos = {}
    for i, lbl in enumerate(rot_labels):
        angle = math.pi / 2 + i * 2 * math.pi / 3
        x = cx + r_rot * math.cos(angle)
        y = cy_rot + r_rot * math.sin(angle)
        rot_pos[lbl] = (x, y)
        draw_node(ax, (x, y), lbl, NODE_CYCLIC, 900, 11)

    # Cyclic edges
    for i in range(3):
        p1 = rot_pos[rot_labels[i]]
        p2 = rot_pos[rot_labels[(i + 1) % 3]]
        draw_curved_edge(ax, p1, p2, CYCLIC_C, arrow=True, connectionstyle="arc3,rad=0.18")
        draw_curved_edge(ax, p2, p1, CYCLIC_C, arrow=True, connectionstyle="arc3,rad=0.18",
                         alpha=0.35, lw=1.5)

    ax.text(cx, cy_rot, "↻", ha="center", va="center", fontsize=26, color=CYCLIC_C, alpha=0.5)

    # Boost extensions (free arrows going down)
    boost_labels = ["K₁", "K₂", "K₃"]
    boost_pos = {}
    cy_boost = 0.18
    for i, lbl in enumerate(boost_labels):
        x = cx + (i - 1) * 0.25
        boost_pos[lbl] = (x, cy_boost)
        draw_node(ax, (x, cy_boost), lbl, NODE_FREE, 900, 11)

    # Rotation → Boost connections (free arrows)
    connections = [("J₁", "K₁"), ("J₂", "K₂"), ("J₃", "K₃")]
    for jl, kl in connections:
        draw_curved_edge(ax, rot_pos[jl], boost_pos[kl], FREE_C, arrow=True,
                         connectionstyle="arc3,rad=0.0", lw=2.5)

    # Cross-connections: [Jᵢ,Kⱼ] = εᵢⱼₖ Kₖ (rotations drive boosts)
    cross = [("J₁", "K₂"), ("J₁", "K₃"), ("J₂", "K₁"), ("J₂", "K₃"),
             ("J₃", "K₁"), ("J₃", "K₂")]
    for jl, kl in cross:
        draw_curved_edge(ax, rot_pos[jl], boost_pos[kl], FREE_C, arrow=True,
                         connectionstyle="arc3,rad=0.12", lw=1.0, alpha=0.3)

    # [Kᵢ,Kⱼ] = -εᵢⱼₖ Jₖ  (boosts return to rotations — dashed curve back)
    for i in range(3):
        p1 = boost_pos[boost_labels[i]]
        p2 = boost_pos[boost_labels[(i + 1) % 3]]
        # Boosts don't close among themselves — they feed back to rotations
        draw_curved_edge(ax, p1, p2, "#e67e22", lw=1.2, alpha=0.25,
                         connectionstyle="arc3,rad=0.2", style="--")

    # Annotations
    ax.text(0.92, 0.72, "Compact core\n3 cyclic edges\nK_λ < 0", color=CYCLIC_C,
            fontsize=9, ha="right", style="italic")
    ax.text(0.92, 0.22, "Non-compact\nextension\n3 free edges\nK_λ > 0", color=FREE_C,
            fontsize=9, ha="right", style="italic")

    # Bracket around boost layer
    ax.annotate("", xy=(0.13, 0.08), xytext=(0.87, 0.08),
                arrowprops=dict(arrowstyle="-", color=FREE_C, lw=1, ls="--"))
    ax.text(0.5, 0.03, "[Kᵢ, Kⱼ] = −εᵢⱼₖ Jₖ   (boosts → rotations, no self-closure)",
            ha="center", color=SUBTLE_C, fontsize=8, family="monospace")

    # Legend label
    ax.text(0.5, 0.96, "Signature: (3⁺, 0₀, 3⁻) — Balanced mixed",
            ha="center", color=SUBTLE_C, fontsize=9, transform=ax.transAxes)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.05, 1.0)
    savefig(fig, "fig3_archetype_II_mixed_lorentz.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: Archetype II extended — su(2,2) conformal
# ══════════════════════════════════════════════════════════════════════════════

def fig4_conformal():
    fig, ax = plt.subplots(1, 1, figsize=(10, 8), facecolor=BG)
    setup_ax(ax, "Archetype II Extended: su(2,2) = so(4,2) — Conformal Algebra")

    # Compact core: s(u(2)⊕u(2)) — 7 generators, arranged as two triangles + center
    cx, cy_core = 0.5, 0.68
    # Two su(2) triangles side by side
    r = 0.13
    core_pos = {}
    # Left su(2)
    for i, lbl in enumerate(["R₁", "R₂", "R₃"]):
        angle = math.pi / 2 + i * 2 * math.pi / 3
        x = 0.3 + r * math.cos(angle)
        y = cy_core + r * math.sin(angle)
        core_pos[lbl] = (x, y)
        draw_node(ax, (x, y), lbl, NODE_CYCLIC, 500, 8)

    # Right su(2)
    for i, lbl in enumerate(["R₄", "R₅", "R₆"]):
        angle = math.pi / 2 + i * 2 * math.pi / 3
        x = 0.7 + r * math.cos(angle)
        y = cy_core + r * math.sin(angle)
        core_pos[lbl] = (x, y)
        draw_node(ax, (x, y), lbl, NODE_CYCLIC, 500, 8)

    # Central shared
    core_pos["R₇"] = (0.5, cy_core)
    draw_node(ax, (0.5, cy_core), "R₇", "#1abc9c", 500, 8)

    # Cyclic edges within each triangle
    for triple in [["R₁", "R₂", "R₃"], ["R₄", "R₅", "R₆"]]:
        for i in range(3):
            p1 = core_pos[triple[i]]
            p2 = core_pos[triple[(i + 1) % 3]]
            draw_curved_edge(ax, p1, p2, CYCLIC_C, arrow=True,
                             connectionstyle="arc3,rad=0.15", lw=2)

    # Connect central to both
    for lbl in ["R₁", "R₃", "R₄", "R₆"]:
        draw_curved_edge(ax, core_pos["R₇"], core_pos[lbl], CYCLIC_C, lw=1.2, alpha=0.4,
                         connectionstyle="arc3,rad=0.1")

    ax.text(0.3, cy_core, "↻", ha="center", va="center", fontsize=20, color=CYCLIC_C, alpha=0.4)
    ax.text(0.7, cy_core, "↻", ha="center", va="center", fontsize=20, color=CYCLIC_C, alpha=0.4)

    # Non-compact extensions: 8 free generators (Pμ, Kμ, boosts, dilatation)
    free_labels = ["P₀", "P₁", "P₂", "P₃", "K₀", "K₁", "K₂", "D"]
    free_pos = {}
    cy_free = 0.2
    for i, lbl in enumerate(free_labels):
        x = 0.12 + i * (0.76 / 7)
        free_pos[lbl] = (x, cy_free)
        draw_node(ax, (x, cy_free), lbl, NODE_FREE, 400, 7)

    # Arrows from core to free
    core_keys = list(core_pos.keys())
    for fl in free_labels:
        # Connect to nearest core node
        fp = free_pos[fl]
        nearest = min(core_keys, key=lambda k: (core_pos[k][0] - fp[0])**2)
        draw_curved_edge(ax, core_pos[nearest], fp, FREE_C, arrow=True,
                         connectionstyle="arc3,rad=0.0", lw=1.5, alpha=0.5)

    # Diverging arrows from free nodes
    for fl in free_labels:
        fp = free_pos[fl]
        ax.annotate("", xy=(fp[0], fp[1] - 0.08), xytext=fp,
                    arrowprops=dict(arrowstyle="->,head_length=5,head_width=3",
                                    color=FREE_C, lw=1.2, alpha=0.5))

    # Annotations
    ax.text(0.5, 0.90, "7 compact edges (K_λ < 0) — s(u(2) ⊕ u(2))",
            ha="center", color=CYCLIC_C, fontsize=10, style="italic")
    ax.text(0.5, 0.07, "8 non-compact edges (K_λ > 0) — translations, special conformal, dilatation",
            ha="center", color=FREE_C, fontsize=9, style="italic")

    # Bounding boxes
    from matplotlib.patches import FancyBboxPatch
    box_core = FancyBboxPatch((0.1, 0.48), 0.8, 0.38, boxstyle="round,pad=0.02",
                               facecolor="none", edgecolor=CYCLIC_C, lw=1.5, ls="--", alpha=0.4)
    ax.add_patch(box_core)
    box_free = FancyBboxPatch((0.05, 0.1), 0.9, 0.2, boxstyle="round,pad=0.02",
                               facecolor="none", edgecolor=FREE_C, lw=1.5, ls="--", alpha=0.4)
    ax.add_patch(box_free)

    ax.text(0.5, 0.0, "Signature: (8⁺, 0₀, 7⁻) — Free-dominated  |  8 > 7: more spacetime than gauge",
            ha="center", color=SUBTLE_C, fontsize=8, family="monospace")

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.05, 1.0)
    savefig(fig, "fig4_archetype_II_conformal.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5: Archetype III — Cycle + Strand (U(3))
# ══════════════════════════════════════════════════════════════════════════════

def fig5_cycle_strand():
    fig, ax = plt.subplots(1, 1, figsize=(10, 7), facecolor=BG)
    setup_ax(ax, "Archetype III: Cycle + Strand — u(3) = su(3) + u(1)")

    # SU(3) core — hexagonal arrangement
    cx, cy = 0.4, 0.5
    r = 0.25
    root_labels = ["λ₁", "λ₂", "λ₄", "λ₅", "λ₆", "λ₇"]
    cartan_labels = ["λ₃", "λ₈"]
    pos = {}

    for i, lbl in enumerate(root_labels):
        angle = math.pi / 6 + i * math.pi / 3
        x = cx + r * math.cos(angle)
        y = cy + r * math.sin(angle)
        pos[lbl] = (x, y)
        draw_node(ax, (x, y), lbl, NODE_CYCLIC, 450, 8)

    pos["λ₃"] = (cx - 0.04, cy)
    pos["λ₈"] = (cx + 0.04, cy)
    draw_node(ax, pos["λ₃"], "λ₃", NODE_CORE, 400, 8)
    draw_node(ax, pos["λ₈"], "λ₈", NODE_CORE, 400, 8)

    # Hexagonal edges
    for i in range(6):
        p1 = pos[root_labels[i]]
        p2 = pos[root_labels[(i + 1) % 6]]
        draw_curved_edge(ax, p1, p2, CYCLIC_C, arrow=True,
                         connectionstyle="arc3,rad=0.08", lw=1.8, alpha=0.8)

    # Cross connections
    for i in range(3):
        p1 = pos[root_labels[i]]
        p2 = pos[root_labels[(i + 3) % 6]]
        draw_curved_edge(ax, p1, p2, CYCLIC_C, lw=0.8, alpha=0.2, connectionstyle="arc3,rad=0.0")

    # Cartan-to-root
    for lbl in root_labels[:3]:
        draw_curved_edge(ax, pos["λ₃"], pos[lbl], CYCLIC_C, lw=0.8, alpha=0.3,
                         connectionstyle="arc3,rad=0.1")
    for lbl in root_labels[3:]:
        draw_curved_edge(ax, pos["λ₈"], pos[lbl], CYCLIC_C, lw=0.8, alpha=0.3,
                         connectionstyle="arc3,rad=-0.1")

    ax.text(cx, cy, "↻", ha="center", va="center", fontsize=24, color=CYCLIC_C, alpha=0.4)

    # ── U(1) strand — isolated ──
    u1_x = 0.85
    draw_node(ax, (u1_x, 0.5), "T₀", NODE_DECOUPLE, 900, 11)
    # Dashed lines radiating (going nowhere)
    for angle_deg in [0, 45, 90, 135, 180, 225, 270, 315]:
        t = math.radians(angle_deg)
        x2 = u1_x + 0.12 * math.cos(t)
        y2 = 0.5 + 0.12 * math.sin(t)
        ax.plot([u1_x, x2], [0.5, y2], color=DECOUPLED_C, lw=1.5, ls=":", alpha=0.4)

    # Gap indicator
    ax.text(0.67, 0.5, "⋮\n✕\n⋮", ha="center", va="center", fontsize=14,
            color=DECOUPLED_C, alpha=0.6)

    # Labels
    from matplotlib.patches import FancyBboxPatch
    box_su3 = FancyBboxPatch((0.08, 0.15), 0.58, 0.7, boxstyle="round,pad=0.02",
                              facecolor="none", edgecolor=CYCLIC_C, lw=1.5, ls="--", alpha=0.5)
    ax.add_patch(box_su3)
    ax.text(0.37, 0.89, "su(3) — 8 cyclic edges", color=CYCLIC_C, fontsize=10,
            ha="center", style="italic")

    box_u1 = FancyBboxPatch((0.73, 0.32), 0.24, 0.36, boxstyle="round,pad=0.02",
                             facecolor="none", edgecolor=DECOUPLED_C, lw=1.5, ls="--", alpha=0.5)
    ax.add_patch(box_u1)
    ax.text(0.85, 0.27, "u(1) — 1 strand\nK = 0", color=DECOUPLED_C, fontsize=9,
            ha="center", style="italic")

    ax.text(0.5, 0.03, "[T_0, lambda_a] = 0  for all a  =>  Disconnected graph components",
            ha="center", color=SUBTLE_C, fontsize=8, family="monospace")

    ax.set_xlim(0.0, 1.05)
    ax.set_ylim(-0.05, 1.0)
    savefig(fig, "fig5_archetype_III_cycle_strand.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 6: Archetype IV (Pure Free) + Archetype V (Frustrated)
# ══════════════════════════════════════════════════════════════════════════════

def fig6_free_and_frustrated():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), facecolor=BG)

    # ── Archetype IV: Pure Free ──
    setup_ax(ax1, "Archetype IV: Pure Free Graph  (K_λ > 0 ∀λ)")

    # Nodes spread out with only outgoing arrows
    labels = ["T₁", "T₂", "T₃", "T₄"]
    positions = [(0.25, 0.7), (0.55, 0.65), (0.35, 0.4), (0.65, 0.35)]
    for (x, y), lbl in zip(positions, labels):
        draw_node(ax1, (x, y), lbl, NODE_FREE, 700, 10)
        # Outgoing arrows
        for angle_deg in np.linspace(-30, 30, 3):
            t = math.radians(angle_deg)
            x2 = x + 0.18 * math.cos(t)
            y2 = y - 0.15 + 0.18 * math.sin(t)
            ax1.annotate("", xy=(x2, y2 - 0.05), xytext=(x, y),
                        arrowprops=dict(arrowstyle="->,head_length=5,head_width=3",
                                        color=FREE_C, lw=1.5, alpha=0.4))

    ax1.text(0.5, 0.1, "Theoretical — not yet observed\nPure non-compact, no rotations\nρ < 1 regime: pre-crystallization",
             ha="center", color=SUBTLE_C, fontsize=9, style="italic")
    ax1.text(0.5, 0.95, "⚠ No cycles = No symmetry",
             ha="center", color=FREE_C, fontsize=11, fontweight="bold", transform=ax1.transAxes)
    ax1.set_xlim(0.0, 1.0)
    ax1.set_ylim(-0.05, 1.0)

    # ── Archetype V: Frustrated ──
    setup_ax(ax2, "Archetype V: Frustrated Graph — so(3,1) -> so(3)")

    # Surviving compact core (triangle)
    cx, cy = 0.5, 0.7
    r = 0.18
    rot_pos = {}
    for i, lbl in enumerate(["J₁", "J₂", "J₃"]):
        angle = math.pi / 2 + i * 2 * math.pi / 3
        x = cx + r * math.cos(angle)
        y = cy + r * math.sin(angle)
        rot_pos[lbl] = (x, y)
        draw_node(ax2, (x, y), lbl, NODE_CYCLIC, 700, 10)

    for i in range(3):
        keys = list(rot_pos.keys())
        p1 = rot_pos[keys[i]]
        p2 = rot_pos[keys[(i + 1) % 3]]
        draw_curved_edge(ax2, p1, p2, CYCLIC_C, arrow=True, connectionstyle="arc3,rad=0.15")

    ax2.text(cx, cy, "↻", ha="center", va="center", fontsize=20, color=CYCLIC_C, alpha=0.5)

    # Collapsed strands (greyed out, with X through them)
    for i, lbl in enumerate(["K₁", "K₂", "K₃"]):
        x = cx + (i - 1) * 0.2
        y = 0.25
        # Faded node
        ax2.scatter(x, y, s=600, c=DECOUPLED_C, zorder=5, alpha=0.25, edgecolors="white",
                    linewidths=0.5)
        ax2.text(x, y, lbl, ha="center", va="center", fontsize=9, color=DECOUPLED_C,
                 alpha=0.5, zorder=6)
        # X through it
        d = 0.04
        ax2.plot([x - d, x + d], [y - d, y + d], color="#e74c3c", lw=2, alpha=0.6)
        ax2.plot([x - d, x + d], [y + d, y - d], color="#e74c3c", lw=2, alpha=0.6)
        # Dashed connection to core (severed)
        nearest_j = list(rot_pos.values())[i]
        ax2.plot([nearest_j[0], x], [nearest_j[1], y], color=DECOUPLED_C,
                 lw=1.5, ls=":", alpha=0.2)

    ax2.text(0.5, 0.07,
             "Compact target + η constraint → boosts collapse\n"
             "Maximal compact subalgebra survives: so(3)",
             ha="center", color=SUBTLE_C, fontsize=9, style="italic")

    # Arrow showing frustration
    ax2.annotate("K_target = −I\n(incompatible)", xy=(0.85, 0.5), fontsize=8,
                 color="#e74c3c", ha="center", style="italic",
                 bbox=dict(boxstyle="round,pad=0.3", fc="#2d1111", ec="#e74c3c", alpha=0.6))

    ax2.set_xlim(0.0, 1.0)
    ax2.set_ylim(-0.05, 1.0)

    fig.suptitle("Archetypes IV & V — Limiting Cases",
                 color=TEXT_C, fontsize=15, fontweight="bold", y=1.0)
    savefig(fig, "fig6_archetypes_IV_V.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 7: Graph Hierarchy — Levels 0 through 4
# ══════════════════════════════════════════════════════════════════════════════

def fig7_hierarchy():
    fig, ax = plt.subplots(1, 1, figsize=(12, 10), facecolor=BG)
    setup_ax(ax)
    ax.set_title("The Spacetime Hierarchy as Graph Evolution",
                 color=TEXT_C, fontsize=15, fontweight="bold", pad=15)

    levels = [
        (0.12, "Level 0", "Pre-crystallization", "ρ < 1", None),
        (0.32, "Level 1", "Internal symmetry", "ρ = 1, compact", "su(2), su(3), …"),
        (0.52, "Level 2", "Spacetime symmetry", "ρ = 1, η + indef.", "so(3,1)"),
        (0.72, "Level 3", "Conformal symmetry", "ρ = 1, η(2,2)", "su(2,2) ≅ so(4,2)"),
        (0.92, "Level 4", "Standard Model", "Product", "su(3) ⊕ su(2) ⊕ u(1)"),
    ]

    for y_base, title, desc, cond, algebra in levels:
        # Level background
        from matplotlib.patches import FancyBboxPatch
        box = FancyBboxPatch((0.05, y_base - 0.06), 0.9, 0.16,
                              boxstyle="round,pad=0.01",
                              facecolor=BG_PANEL, edgecolor=GRID_C, lw=1, alpha=0.8)
        ax.add_patch(box)

        # Title
        ax.text(0.08, y_base + 0.07, title, color=TEXT_C, fontsize=12, fontweight="bold",
                va="center")
        ax.text(0.08, y_base + 0.01, desc, color=SUBTLE_C, fontsize=9, va="center",
                style="italic")
        ax.text(0.08, y_base - 0.04, cond, color=SUBTLE_C, fontsize=8, va="center",
                family="monospace")

    # Level 0: scattered dots
    y0 = 0.12
    for x in np.linspace(0.45, 0.85, 6):
        ax.scatter(x, y0, s=80, c=DECOUPLED_C, alpha=0.4, zorder=5)
    ax.text(0.65, y0 - 0.04, "No edges · No cycles · No physics", ha="center",
            color=DECOUPLED_C, fontsize=7)

    # Level 1: SU(2) triangle
    y1, cx1 = 0.32, 0.6
    r1 = 0.04
    for i, col in enumerate([NODE_CYCLIC] * 3):
        angle = math.pi / 2 + i * 2 * math.pi / 3
        x = cx1 + r1 * math.cos(angle)
        y = y1 + r1 * math.sin(angle)
        ax.scatter(x, y, s=120, c=col, zorder=5, edgecolors="white", linewidths=0.5)
    # Edges
    for i in range(3):
        a1 = math.pi / 2 + i * 2 * math.pi / 3
        a2 = math.pi / 2 + ((i + 1) % 3) * 2 * math.pi / 3
        ax.plot([cx1 + r1 * math.cos(a1), cx1 + r1 * math.cos(a2)],
                [y1 + r1 * math.sin(a1), y1 + r1 * math.sin(a2)],
                color=CYCLIC_C, lw=2, zorder=4)
    ax.text(cx1, y1, "↻", ha="center", va="center", fontsize=10, color=CYCLIC_C, alpha=0.6)
    ax.text(cx1 + 0.12, y1, "su(2)", color=CYCLIC_C, fontsize=9, va="center")

    # Level 2: SO(3,1)
    y2, cx2 = 0.52, 0.55
    # Triangle top
    for i in range(3):
        angle = math.pi / 2 + i * 2 * math.pi / 3
        x = cx2 + 0.04 * math.cos(angle)
        y = y2 + 0.025 + 0.04 * math.sin(angle)
        ax.scatter(x, y, s=100, c=NODE_CYCLIC, zorder=5, edgecolors="white", linewidths=0.5)
    for i in range(3):
        a1 = math.pi / 2 + i * 2 * math.pi / 3
        a2 = math.pi / 2 + ((i + 1) % 3) * 2 * math.pi / 3
        ax.plot([cx2 + 0.04 * math.cos(a1), cx2 + 0.04 * math.cos(a2)],
                [y2 + 0.025 + 0.04 * math.sin(a1), y2 + 0.025 + 0.04 * math.sin(a2)],
                color=CYCLIC_C, lw=2, zorder=4)
    # Boost arrows
    for dx in [-0.04, 0, 0.04]:
        ax.annotate("", xy=(cx2 + dx, y2 - 0.045), xytext=(cx2 + dx, y2 - 0.01),
                    arrowprops=dict(arrowstyle="->,head_length=3,head_width=2",
                                    color=FREE_C, lw=1.5))
    ax.text(cx2 + 0.15, y2, "so(3,1)", color=TEXT_C, fontsize=9, va="center")

    # Level 3: SU(2,2)
    y3, cx3 = 0.72, 0.55
    # Compact blob
    from matplotlib.patches import Ellipse
    ell = Ellipse((cx3 - 0.03, y3 + 0.02), 0.1, 0.06, facecolor=CYCLIC_C, alpha=0.2, zorder=3)
    ax.add_patch(ell)
    ax.text(cx3 - 0.03, y3 + 0.02, "7↻", ha="center", va="center",
            color=CYCLIC_C, fontsize=9, fontweight="bold")
    # Free arrows
    for i in range(4):
        ax.annotate("", xy=(cx3 + 0.08 + i * 0.025, y3 - 0.035),
                    xytext=(cx3 + 0.08 + i * 0.025, y3 + 0.01),
                    arrowprops=dict(arrowstyle="->,head_length=3,head_width=2",
                                    color=FREE_C, lw=1.2))
    ax.text(cx3 + 0.12, y3 - 0.04, "8→", color=FREE_C, fontsize=9)
    ax.text(cx3 + 0.22, y3, "su(2,2)", color=TEXT_C, fontsize=9, va="center")

    # Level 4: SM product
    y4 = 0.92
    # SU(3) blob
    cx_su3 = 0.48
    ell_su3 = Ellipse((cx_su3, y4), 0.08, 0.06, facecolor=CYCLIC_C, alpha=0.2, zorder=3)
    ax.add_patch(ell_su3)
    ax.text(cx_su3, y4, "8↻", ha="center", va="center", color=CYCLIC_C, fontsize=9, fontweight="bold")
    ax.text(cx_su3, y4 + 0.045, "su(3)", ha="center", color=CYCLIC_C, fontsize=8)

    # SU(2) blob
    cx_su2 = 0.63
    ell_su2 = Ellipse((cx_su2, y4), 0.055, 0.05, facecolor=CYCLIC_C, alpha=0.15, zorder=3)
    ax.add_patch(ell_su2)
    ax.text(cx_su2, y4, "3↻", ha="center", va="center", color=CYCLIC_C, fontsize=9, fontweight="bold")
    ax.text(cx_su2, y4 + 0.04, "su(2)", ha="center", color=CYCLIC_C, fontsize=8)

    # U(1) strand
    cx_u1 = 0.76
    ax.scatter(cx_u1, y4, s=100, c=DECOUPLED_C, zorder=5, edgecolors="white", linewidths=0.5)
    ax.text(cx_u1 + 0.03, y4, "⋯", color=DECOUPLED_C, fontsize=12)
    ax.text(cx_u1, y4 + 0.04, "u(1)", ha="center", color=DECOUPLED_C, fontsize=8)

    # Connecting arrows between levels
    for y_from, y_to in [(0.21, 0.27), (0.41, 0.47), (0.61, 0.67), (0.81, 0.87)]:
        ax.annotate("", xy=(0.35, y_to), xytext=(0.35, y_from),
                    arrowprops=dict(arrowstyle="->,head_length=1,head_width=1",
                                    color=SUBTLE_C, lw=1.5, alpha=0.5))

    # Axis labels
    ax.text(0.02, 0.5, "Complexity  →", rotation=90, ha="center", va="center",
            color=SUBTLE_C, fontsize=11, transform=ax.transAxes)
    ax.text(0.35, -0.02, "ρ (relational density)  →", ha="center",
            color=SUBTLE_C, fontsize=10)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.05, 1.05)
    savefig(fig, "fig7_graph_hierarchy.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 8: The Universe's Graph — Layered composition
# ══════════════════════════════════════════════════════════════════════════════

def fig8_universe_graph():
    fig, ax = plt.subplots(1, 1, figsize=(12, 10), facecolor=BG)
    setup_ax(ax, "The Universe's Graph — Layered Composition of All Archetypes")

    from matplotlib.patches import FancyBboxPatch, Ellipse

    # ── Conformal layer (bottom, largest) ──
    box_conf = FancyBboxPatch((0.05, 0.05), 0.9, 0.25,
                               boxstyle="round,pad=0.02",
                               facecolor="#1a1a2e", edgecolor=PURPLE_C,
                               lw=2, alpha=0.6)
    ax.add_patch(box_conf)
    ax.text(0.5, 0.27, "su(2,2) = so(4,2)  —  Conformal Layer",
            ha="center", color=PURPLE_C, fontsize=11, fontweight="bold")

    # Compact core inside conformal
    ell_conf_compact = Ellipse((0.3, 0.15), 0.2, 0.1,
                                facecolor=CYCLIC_C, alpha=0.1, zorder=3)
    ax.add_patch(ell_conf_compact)
    ax.text(0.3, 0.15, "7 ↻", ha="center", va="center",
            color=CYCLIC_C, fontsize=12, fontweight="bold")

    # Free arrows
    for i in range(8):
        x = 0.5 + i * 0.04
        ax.annotate("", xy=(x, 0.08), xytext=(x, 0.18),
                    arrowprops=dict(arrowstyle="->,head_length=4,head_width=2.5",
                                    color=FREE_C, lw=1.5, alpha=0.6))
    ax.text(0.7, 0.08, "8 →", color=FREE_C, fontsize=10)

    # ── Spacetime layer (middle) ──
    box_st = FancyBboxPatch((0.12, 0.35), 0.76, 0.25,
                             boxstyle="round,pad=0.02",
                             facecolor="#1a2a1a", edgecolor="#e67e22",
                             lw=2, alpha=0.6)
    ax.add_patch(box_st)
    ax.text(0.5, 0.57, "so(3,1)  —  Spacetime Backbone",
            ha="center", color="#e67e22", fontsize=11, fontweight="bold")

    # Rotation core
    cx_r, cy_r = 0.35, 0.47
    r_r = 0.055
    for i in range(3):
        angle = math.pi / 2 + i * 2 * math.pi / 3
        x = cx_r + r_r * math.cos(angle)
        y = cy_r + r_r * math.sin(angle)
        ax.scatter(x, y, s=150, c=NODE_CYCLIC, zorder=5, edgecolors="white", linewidths=0.5)
    for i in range(3):
        a1 = math.pi / 2 + i * 2 * math.pi / 3
        a2 = math.pi / 2 + ((i + 1) % 3) * 2 * math.pi / 3
        ax.plot([cx_r + r_r * math.cos(a1), cx_r + r_r * math.cos(a2)],
                [cy_r + r_r * math.sin(a1), cy_r + r_r * math.sin(a2)],
                color=CYCLIC_C, lw=2.5, zorder=4)
    ax.text(cx_r, cy_r, "↻", ha="center", va="center", fontsize=14, color=CYCLIC_C, alpha=0.5)
    ax.text(cx_r, cy_r - 0.08, "3 rotations", ha="center", color=CYCLIC_C, fontsize=8)

    # Boost arrows
    for j in range(3):
        x_b = 0.55 + j * 0.06
        ax.annotate("", xy=(x_b, 0.39), xytext=(x_b, 0.50),
                    arrowprops=dict(arrowstyle="->,head_length=5,head_width=3",
                                    color=FREE_C, lw=2, alpha=0.7))
    ax.text(0.61, 0.38, "3 boosts", ha="center", color=FREE_C, fontsize=8)

    # ── Gauge layer (top) ──
    box_gauge = FancyBboxPatch((0.08, 0.67), 0.84, 0.28,
                                boxstyle="round,pad=0.02",
                                facecolor="#1a1a30", edgecolor=CYCLIC_C,
                                lw=2, alpha=0.6)
    ax.add_patch(box_gauge)
    ax.text(0.5, 0.93, "Internal Gauge Symmetries  —  Standard Model",
            ha="center", color=CYCLIC_C, fontsize=11, fontweight="bold")

    # SU(3)_color
    cx_su3, cy_su3 = 0.25, 0.80
    ell_su3 = Ellipse((cx_su3, cy_su3), 0.22, 0.14,
                       facecolor=CYCLIC_C, alpha=0.1, zorder=3)
    ax.add_patch(ell_su3)
    # Hexagonal structure
    r_hex = 0.06
    for i in range(6):
        angle = i * math.pi / 3
        x = cx_su3 + r_hex * math.cos(angle)
        y = cy_su3 + r_hex * math.sin(angle)
        ax.scatter(x, y, s=60, c=NODE_CYCLIC, zorder=5, edgecolors="white", linewidths=0.3)
    for i in range(6):
        a1 = i * math.pi / 3
        a2 = ((i + 1) % 6) * math.pi / 3
        ax.plot([cx_su3 + r_hex * math.cos(a1), cx_su3 + r_hex * math.cos(a2)],
                [cy_su3 + r_hex * math.sin(a1), cy_su3 + r_hex * math.sin(a2)],
                color=CYCLIC_C, lw=1.5, zorder=4, alpha=0.8)
    # Cartan
    ax.scatter(cx_su3 - 0.015, cy_su3, s=40, c=NODE_CORE, zorder=5)
    ax.scatter(cx_su3 + 0.015, cy_su3, s=40, c=NODE_CORE, zorder=5)
    ax.text(cx_su3, cy_su3 - 0.09, "su(3)_color", ha="center", color=CYCLIC_C, fontsize=9,
            fontweight="bold")
    ax.text(cx_su3, cy_su3 + 0.09, "8 ↻", ha="center", color=CYCLIC_C, fontsize=10)

    # SU(2)_weak
    cx_su2, cy_su2 = 0.55, 0.80
    r2 = 0.045
    for i in range(3):
        angle = math.pi / 2 + i * 2 * math.pi / 3
        x = cx_su2 + r2 * math.cos(angle)
        y = cy_su2 + r2 * math.sin(angle)
        ax.scatter(x, y, s=100, c=NODE_CYCLIC, zorder=5, edgecolors="white", linewidths=0.4)
    for i in range(3):
        a1 = math.pi / 2 + i * 2 * math.pi / 3
        a2 = math.pi / 2 + ((i + 1) % 3) * 2 * math.pi / 3
        ax.plot([cx_su2 + r2 * math.cos(a1), cx_su2 + r2 * math.cos(a2)],
                [cy_su2 + r2 * math.sin(a1), cy_su2 + r2 * math.sin(a2)],
                color=CYCLIC_C, lw=2, zorder=4)
    ax.text(cx_su2, cy_su2, "↻", ha="center", va="center", fontsize=12, color=CYCLIC_C, alpha=0.5)
    ax.text(cx_su2, cy_su2 - 0.08, "su(2)_weak", ha="center", color=CYCLIC_C, fontsize=9,
            fontweight="bold")
    ax.text(cx_su2, cy_su2 + 0.075, "3 ↻", ha="center", color=CYCLIC_C, fontsize=10)

    # U(1)_Y
    cx_u1, cy_u1 = 0.80, 0.80
    ax.scatter(cx_u1, cy_u1, s=200, c=DECOUPLED_C, zorder=5, edgecolors="white", linewidths=0.6)
    ax.text(cx_u1, cy_u1, "T₀", ha="center", va="center", fontsize=8, color=BG,
            fontweight="bold", zorder=6)
    for angle_deg in np.arange(0, 360, 45):
        t = math.radians(angle_deg)
        x2 = cx_u1 + 0.04 * math.cos(t)
        y2 = cy_u1 + 0.04 * math.sin(t)
        ax.plot([cx_u1, x2], [cy_u1, y2], color=DECOUPLED_C, lw=1, ls=":", alpha=0.5)
    ax.text(cx_u1, cy_u1 - 0.08, "u(1)_Y", ha="center", color=DECOUPLED_C, fontsize=9,
            fontweight="bold")
    ax.text(cx_u1, cy_u1 + 0.065, "1 ⋯", ha="center", color=DECOUPLED_C, fontsize=10)

    # ── Connecting arrows between layers ──
    ax.annotate("", xy=(0.4, 0.60), xytext=(0.4, 0.67),
                arrowprops=dict(arrowstyle="<->", color=SUBTLE_C, lw=1.5, ls="--"))
    ax.text(0.43, 0.635, "embedding", color=SUBTLE_C, fontsize=7, style="italic")

    ax.annotate("", xy=(0.4, 0.30), xytext=(0.4, 0.35),
                arrowprops=dict(arrowstyle="<->", color=SUBTLE_C, lw=1.5, ls="--"))
    ax.text(0.43, 0.325, "embedding", color=SUBTLE_C, fontsize=7, style="italic")

    # Legend
    items = [
        (CYCLIC_C, "↻ Cyclic edges — compact / gauge"),
        (FREE_C, "→ Free edges — non-compact / causal"),
        (DECOUPLED_C, "⋯ Decoupled strand — abelian / phase"),
    ]
    for i, (c, txt) in enumerate(items):
        y_leg = 0.01 - i * 0.025
        ax.plot([0.05, 0.09], [y_leg, y_leg], color=c, lw=3)
        ax.text(0.10, y_leg, txt, color=c, fontsize=8, va="center")

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.08, 1.0)
    savefig(fig, "fig8_universe_graph.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 9: Density transition (percolation)
# ══════════════════════════════════════════════════════════════════════════════

def fig9_density_transition():
    fig, ax = plt.subplots(1, 1, figsize=(10, 6), facecolor=BG)
    ax.set_facecolor(BG_PANEL)

    # Step function: symmetry = 0 for rho < 1, 1 for rho >= 1
    rho = np.linspace(0, 2, 500)
    symmetry = np.where(rho < 1.0, 0.0, 1.0)

    # Smooth version for visual context
    smooth = 1.0 / (1.0 + np.exp(-50 * (rho - 1.0)))

    ax.fill_between(rho, 0, smooth, alpha=0.15, color=CYCLIC_C)
    ax.plot(rho, smooth, color=CYCLIC_C, lw=3, path_effects=GLOW)

    # ρ = 1 vertical line
    ax.axvline(x=1.0, color=ACCENT_C, lw=2, ls="--", alpha=0.7)
    ax.text(1.02, 0.5, "ρ = 1", color=ACCENT_C, fontsize=12, fontweight="bold", rotation=90,
            va="center")

    # Annotations
    ax.text(0.4, 0.15, "ρ < 1\nPre-crystallization\nNo cycles · No symmetry",
            ha="center", color=FREE_C, fontsize=10, style="italic",
            bbox=dict(facecolor=BG, edgecolor=FREE_C, alpha=0.5, boxstyle="round,pad=0.3"))

    ax.text(1.5, 0.85, "ρ ≥ 1\nFull crystallization\nAll cycles close\nGauge symmetry exact",
            ha="center", color=CYCLIC_C, fontsize=10, style="italic",
            bbox=dict(facecolor=BG, edgecolor=CYCLIC_C, alpha=0.5, boxstyle="round,pad=0.3"))

    # Small graph illustrations
    # ρ < 1: scattered dots
    for x in np.linspace(0.2, 0.6, 5):
        ax.scatter(x, 0.92, s=40, c=DECOUPLED_C, alpha=0.5, zorder=5,
                   transform=ax.transData, clip_on=False)

    # ρ = 1: triangle
    cx_t, cy_t = 1.3, 0.92
    r_t = 0.05
    for i in range(3):
        a1 = math.pi / 2 + i * 2 * math.pi / 3
        a2 = math.pi / 2 + ((i + 1) % 3) * 2 * math.pi / 3
        ax.plot([cx_t + r_t * math.cos(a1), cx_t + r_t * math.cos(a2)],
                [cy_t + r_t * 3 * math.sin(a1), cy_t + r_t * 3 * math.sin(a2)],
                color=CYCLIC_C, lw=2, zorder=4)
        ax.scatter(cx_t + r_t * math.cos(a1), cy_t + r_t * 3 * math.sin(a1),
                   s=50, c=NODE_CYCLIC, zorder=5, edgecolors="white", linewidths=0.3)

    ax.set_xlabel("ρ  =  n_gen / (d² − 1)", color=TEXT_C, fontsize=12)
    ax.set_ylabel("Symmetry crystallization", color=TEXT_C, fontsize=12)
    ax.set_title("The Density Threshold — Percolation Transition on the Graph",
                 color=TEXT_C, fontsize=14, fontweight="bold", pad=12)

    ax.set_xlim(0, 2)
    ax.set_ylim(-0.05, 1.1)
    ax.tick_params(colors=SUBTLE_C)
    ax.spines["bottom"].set_color(GRID_C)
    ax.spines["left"].set_color(GRID_C)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    savefig(fig, "fig9_density_transition.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 10: Potential → Cycle → Invariant chain
# ══════════════════════════════════════════════════════════════════════════════

def fig10_chain():
    fig, ax = plt.subplots(1, 1, figsize=(14, 5), facecolor=BG)
    setup_ax(ax, "The Potential  →  Cycle  →  Invariant Chain")

    from matplotlib.patches import FancyBboxPatch

    # Use pixel-friendly absolute coords; axes limits set at end
    box_w, box_h = 3.0, 3.0
    gap = 1.8
    x0 = 0.5

    boxes_data = [
        (x0,                       "POTENTIAL", "Killing matrix K_ab\nAdjoint self-interaction\nGraph metric",        ACCENT_C),
        (x0 + box_w + gap,        "CYCLE",     "Structure constants f_abc\nCommutator relations\nGraph topology",     CYCLIC_C),
        (x0 + 2*(box_w + gap),    "INVARIANT", "h_v, Casimirs, I_R\nTopological invariants\nPhysics",                PURPLE_C),
    ]

    y_base = 0.6
    for bx, title, desc, color in boxes_data:
        box = FancyBboxPatch((bx, y_base), box_w, box_h,
                              boxstyle="round,pad=0.15",
                              facecolor=BG_PANEL, edgecolor=color, lw=2.5,
                              transform=ax.transData, zorder=2)
        ax.add_patch(box)
        # Glow
        box_glow = FancyBboxPatch((bx - 0.05, y_base - 0.05), box_w + 0.1, box_h + 0.1,
                                   boxstyle="round,pad=0.2",
                                   facecolor="none", edgecolor=color, lw=5, alpha=0.15,
                                   transform=ax.transData, zorder=1)
        ax.add_patch(box_glow)

        cx = bx + box_w / 2
        ax.text(cx, y_base + box_h - 0.55, title, ha="center", va="center",
                color=color, fontsize=15, fontweight="bold", zorder=3)
        ax.text(cx, y_base + 0.95, desc, ha="center", va="center",
                color=TEXT_C, fontsize=9, family="monospace", zorder=3)

    y_mid = y_base + box_h / 2

    # Simple triangle-head arrows drawn as lines + filled polygon
    for i, (label_top, label_bot) in enumerate([
        ("adjoint\naction", "Variational\nminimization"),
        ("holonomy", "Mathematical\nconsequence"),
    ]):
        bx = boxes_data[i][0]
        x_start = bx + box_w + 0.12
        x_end   = boxes_data[i + 1][0] - 0.12
        x_tip   = x_end
        head    = 0.25

        # Shaft
        ax.plot([x_start, x_tip - head], [y_mid, y_mid],
                color=TEXT_C, lw=2.5, solid_capstyle="butt", zorder=4)
        # Arrowhead (filled triangle)
        ax.fill([x_tip - head, x_tip - head, x_tip],
                [y_mid + 0.22, y_mid - 0.22, y_mid],
                color=TEXT_C, zorder=4)

        x_label = (x_start + x_end) / 2
        ax.text(x_label, y_mid + 0.50, label_top,
                ha="center", va="bottom", color=SUBTLE_C, fontsize=9, style="italic", zorder=4)
        ax.text(x_label, y_base - 0.35, label_bot,
                ha="center", va="top", color=SUBTLE_C, fontsize=9, style="italic", zorder=4)

    total_w = x0 + 3 * box_w + 2 * gap + 0.5
    ax.set_xlim(-0.1, total_w + 0.1)
    ax.set_ylim(-0.3, box_h + y_base + 0.6)
    ax.set_aspect("equal")
    savefig(fig, "fig10_potential_cycle_invariant.png")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

ALL_FIGURES = {
    1: ("Three Edge Classes", fig1_edge_classes),
    2: ("Archetype I: Pure Cycle (SU(2), SU(3))", fig2_pure_cycle),
    3: ("Archetype II: Mixed — SO(3,1) Lorentz", fig3_mixed_lorentz),
    4: ("Archetype II ext: SU(2,2) Conformal", fig4_conformal),
    5: ("Archetype III: Cycle + Strand — U(3)", fig5_cycle_strand),
    6: ("Archetypes IV & V: Pure Free / Frustrated", fig6_free_and_frustrated),
    7: ("Graph Hierarchy (Levels 0–4)", fig7_hierarchy),
    8: ("Universe's Graph — Layered Composition", fig8_universe_graph),
    9: ("Density Transition (Percolation)", fig9_density_transition),
    10: ("Potential → Cycle → Invariant Chain", fig10_chain),
}


def main():
    parser = argparse.ArgumentParser(description="Generate graph topology figures")
    parser.add_argument("--fig", type=int, nargs="*",
                        help="Figure number(s) to generate (default: all)")
    args = parser.parse_args()

    figs = args.fig if args.fig else sorted(ALL_FIGURES.keys())

    print(f"Generating {len(figs)} figure(s) → {OUT_DIR}/")
    print()
    for n in figs:
        if n not in ALL_FIGURES:
            print(f"  ✗ Figure {n} not found (valid: 1–{max(ALL_FIGURES)})")
            continue
        desc, fn = ALL_FIGURES[n]
        print(f"  [{n:2d}] {desc}")
        fn()

    print(f"\nDone. {len(figs)} figures saved to {OUT_DIR}/")


if __name__ == "__main__":
    main()
