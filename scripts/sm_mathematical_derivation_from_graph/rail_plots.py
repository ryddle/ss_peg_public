#!/usr/bin/env python3
"""
rail_plots.py — Visual exploration of the marble-on-a-rail picture
===================================================================
Plots:
  1. 3D rails in lattice space (one per sector)
  2. V(g) = d²(x(g), lattice) — the potential wells
  3. Speed |v(g)| along each rail
  4. Phase portrait: position vs velocity per coordinate
  5. Conserved quantities / invariants search
  6. Inter-sector comparison: what's universal?
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from fractions import Fraction as F
from collections import OrderedDict
from mpl_toolkits.mplot3d import Axes3D

plt.rcParams.update({
    'figure.facecolor': '#0e1117',
    'axes.facecolor': '#0e1117',
    'axes.edgecolor': '#444',
    'axes.labelcolor': '#ccc',
    'text.color': '#ccc',
    'xtick.color': '#999',
    'ytick.color': '#999',
    'grid.color': '#333',
    'legend.facecolor': '#1a1a2e',
    'legend.edgecolor': '#444',
    'font.size': 10,
})

COLORS = {
    'lepton': '#ff6b6b',
    'up': '#4ecdc4',
    'down': '#ffe66d',
}
GEN_COLORS = ['#ff4444', '#44ff44', '#4488ff']
COORD_NAMES = ['p (ln2)', 'q (lnR₃)', 'r (ln3)']
COORD_SHORT = ['p', 'q', 'r']
PARTICLE_NAMES = {
    'lepton': ['e', 'μ', 'τ'],
    'up': ['u', 'c', 't'],
    'down': ['d', 's', 'b'],
}

# ── Data ──
PARTICLES = OrderedDict([
    ("e",   (F(0),    F(0),    F(0))),
    ("mu",  (F(-1,2), F(-4),   F(3))),
    ("tau", (F(1),    F(-7),   F(3))),
    ("u",   (F(-3,2), F(-6),   F(-1))),
    ("c",   (F(1,2),  F(-7),   F(3))),
    ("t",   (F(6),    F(-7),   F(4))),
    ("d",   (F(-1,2), F(-8),   F(-2))),
    ("s",   (F(3,2),  F(-7),   F(0))),
    ("b",   (F(3,2),  F(-6),   F(4))),
])

SECTORS = OrderedDict([
    ("lepton", ["e", "mu", "tau"]),
    ("up",     ["u", "c", "t"]),
    ("down",   ["d", "s", "b"]),
])

FLAT_COORD = {"lepton": 2, "up": 1, "down": 0}

def fit_parabola(x1, x2, x3):
    x1, x2, x3 = F(x1), F(x2), F(x3)
    a = (x1 - 2*x2 + x3) / 2
    b = (-5*x1 + 8*x2 - 3*x3) / 2
    c = 3*x1 - 3*x2 + x3
    return float(a), float(b), float(c)

coeffs = {}
for sname, names in SECTORS.items():
    coords = [PARTICLES[n] for n in names]
    sector_coeffs = []
    for k in range(3):
        a, b, c = fit_parabola(coords[0][k], coords[1][k], coords[2][k])
        sector_coeffs.append((a, b, c))
    coeffs[sname] = sector_coeffs

def x_vec(sector, g):
    return np.array([a*g**2 + b*g + c for a, b, c in coeffs[sector]])

def v_vec(sector, g):
    return np.array([2*a*g + b for a, b, c in coeffs[sector]])

def dist_to_lattice(x):
    nearest = np.round(2 * x) / 2
    return np.linalg.norm(x - nearest)

g_fine = np.linspace(0.2, 3.8, 1000)
g_gens = [1, 2, 3]

OUT = "output/rail_plots"
import os; os.makedirs(OUT, exist_ok=True)

# ════════════════════════════════════════════════════════════════
# FIGURE 1: 3D Rails in lattice space
# ════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(18, 6))
for idx, sname in enumerate(SECTORS):
    ax = fig.add_subplot(1, 3, idx + 1, projection='3d')
    ax.set_facecolor('#0e1117')

    # Rail
    trail = np.array([x_vec(sname, g) for g in g_fine])
    ax.plot(trail[:, 0], trail[:, 1], trail[:, 2],
            color=COLORS[sname], alpha=0.7, linewidth=1.5, label='rail')

    # Generation points
    for gi, gname in zip(g_gens, PARTICLE_NAMES[sname]):
        x = x_vec(sname, gi)
        ax.scatter(*x, s=120, c=GEN_COLORS[gi-1], zorder=5, edgecolors='white', linewidth=0.5)
        ax.text(x[0]+0.3, x[1]+0.3, x[2]+0.3, f'{gname} (g={gi})',
                color=GEN_COLORS[gi-1], fontsize=9, fontweight='bold')

    # Lattice points near the rail (half-integer grid)
    p_range = np.arange(int(trail[:, 0].min()) - 1, int(trail[:, 0].max()) + 2, 0.5)
    q_range = np.arange(int(trail[:, 1].min()) - 1, int(trail[:, 1].max()) + 2, 0.5)
    r_range = np.arange(int(trail[:, 2].min()) - 1, int(trail[:, 2].max()) + 2, 0.5)

    # Only show lattice points within 1.5 of the rail
    for pp in p_range:
        for qq in q_range:
            for rr in r_range:
                pt = np.array([pp, qq, rr])
                dists = np.linalg.norm(trail - pt, axis=1)
                if dists.min() < 1.0:
                    ax.scatter(pp, qq, rr, s=8, c='#555', alpha=0.3, marker='.')

    ax.set_xlabel('p (ln2)', fontsize=9)
    ax.set_ylabel('q (lnR₃)', fontsize=9)
    ax.set_zlabel('r (ln3)', fontsize=9)
    ax.set_title(f'{sname.upper()} rail', fontsize=12, color=COLORS[sname], pad=10)

plt.tight_layout()
plt.savefig(f'{OUT}/01_3D_rails.png', dpi=150, bbox_inches='tight')
plt.close()
print("  [1/8] 3D rails saved")

# ════════════════════════════════════════════════════════════════
# FIGURE 2: V(g) = d²(x(g), lattice) — the potential wells
# ════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
for idx, sname in enumerate(SECTORS):
    ax = axes[idx]
    V = np.array([dist_to_lattice(x_vec(sname, g))**2 for g in g_fine])

    ax.fill_between(g_fine, V, alpha=0.3, color=COLORS[sname])
    ax.plot(g_fine, V, color=COLORS[sname], linewidth=2)

    for gi, gname in zip(g_gens, PARTICLE_NAMES[sname]):
        ax.axvline(gi, color=GEN_COLORS[gi-1], alpha=0.4, linestyle='--', linewidth=0.8)
        ax.plot(gi, 0, 'v', color=GEN_COLORS[gi-1], markersize=12, zorder=5)
        ax.text(gi, -0.008, gname, ha='center', color=GEN_COLORS[gi-1],
                fontsize=11, fontweight='bold')

    ax.set_xlabel('g (generation parameter)', fontsize=11)
    ax.set_title(f'{sname.upper()} — V(g) = d²(x, lattice)', color=COLORS[sname], fontsize=12)
    ax.set_ylim(-0.015, 0.20)
    ax.grid(True, alpha=0.2)

axes[0].set_ylabel('V(g)', fontsize=12)
plt.tight_layout()
plt.savefig(f'{OUT}/02_potential_wells.png', dpi=150, bbox_inches='tight')
plt.close()
print("  [2/8] Potential wells saved")

# ════════════════════════════════════════════════════════════════
# FIGURE 3: Speed |v(g)| — marble velocity
# ════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(12, 6))
for sname in SECTORS:
    speeds = np.array([np.linalg.norm(v_vec(sname, g)) for g in g_fine])
    ax.plot(g_fine, speeds, color=COLORS[sname], linewidth=2, label=sname)

    for gi, gname in zip(g_gens, PARTICLE_NAMES[sname]):
        s = np.linalg.norm(v_vec(sname, gi))
        ax.plot(gi, s, 'o', color=COLORS[sname], markersize=10,
                markeredgecolor='white', markeredgewidth=1, zorder=5)
        offset = 0.15 if sname != 'down' else -0.3
        ax.annotate(gname, (gi, s), textcoords="offset points",
                    xytext=(8, offset * 30), color=COLORS[sname],
                    fontsize=10, fontweight='bold')

for gi in g_gens:
    ax.axvline(gi, color='#555', alpha=0.3, linestyle='--')

ax.set_xlabel('g (generation parameter)', fontsize=12)
ax.set_ylabel('|v(g)| = |dx/dg|', fontsize=12)
ax.set_title('Marble speed along the rail', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.2)
plt.tight_layout()
plt.savefig(f'{OUT}/03_speed_profile.png', dpi=150, bbox_inches='tight')
plt.close()
print("  [3/8] Speed profile saved")

# ════════════════════════════════════════════════════════════════
# FIGURE 4: Coordinate-by-coordinate x_k(g) — the 9 parabolas
# ════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(3, 3, figsize=(16, 12), sharex=True)

for row, sname in enumerate(SECTORS):
    flat_k = FLAT_COORD[sname]
    for col in range(3):
        ax = axes[row][col]
        a, b, c = coeffs[sname][col]
        x_vals = np.array([a*g**2 + b*g + c for g in g_fine])
        v_vals = np.array([2*a*g + b for g in g_fine])

        # Position curve
        color = COLORS[sname]
        style = '-' if col != flat_k else '--'
        ax.plot(g_fine, x_vals, color=color, linewidth=2, linestyle=style)

        # Velocity as shaded area (normalized)
        v_norm = v_vals / max(abs(v_vals.max()), abs(v_vals.min()), 1) * (x_vals.max() - x_vals.min()) * 0.2
        ax.fill_between(g_fine, x_vals, x_vals + v_norm, alpha=0.15, color=color)

        # Generation points
        for gi, gname in zip(g_gens, PARTICLE_NAMES[sname]):
            xg = a*gi**2 + b*gi + c
            ax.plot(gi, xg, 'o', color=GEN_COLORS[gi-1], markersize=10,
                    markeredgecolor='white', markeredgewidth=1, zorder=5)
            ax.annotate(gname, (gi, xg), textcoords="offset points",
                        xytext=(8, 5), color=GEN_COLORS[gi-1], fontsize=9)

        # Turning point
        if abs(a) > 1e-10:
            g_star = -b / (2 * a)
            if 0 < g_star < 4:
                x_star = a * g_star**2 + b * g_star + c
                ax.axvline(g_star, color='#888', alpha=0.3, linestyle=':')
                ax.plot(g_star, x_star, '*', color='white', markersize=8, alpha=0.6)

        flat_tag = " [FLAT]" if col == flat_k else ""
        db = b + 4 * a
        ax.set_title(f'{sname}/{COORD_SHORT[col]}  a={a:.2f} δb={db:+.2f}{flat_tag}',
                     fontsize=10, color=color)
        ax.grid(True, alpha=0.15)

        if row == 2:
            ax.set_xlabel('g', fontsize=10)
        if col == 0:
            ax.set_ylabel(f'{sname}', fontsize=11, color=color)

plt.suptitle('Nine parabolas x_k(g) — dashed = flat coordinate, star = turning point',
             fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig(f'{OUT}/04_nine_parabolas.png', dpi=150, bbox_inches='tight')
plt.close()
print("  [4/8] Nine parabolas saved")

# ════════════════════════════════════════════════════════════════
# FIGURE 5: Phase portraits v_k vs x_k for each coordinate
# ════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(3, 3, figsize=(16, 12))

for row, sname in enumerate(SECTORS):
    flat_k = FLAT_COORD[sname]
    for col in range(3):
        ax = axes[row][col]
        a, b, c = coeffs[sname][col]

        x_vals = np.array([a*g**2 + b*g + c for g in g_fine])
        v_vals = np.array([2*a*g + b for g in g_fine])

        color = COLORS[sname]
        ax.plot(x_vals, v_vals, color=color, linewidth=1.5, alpha=0.7)

        # Direction arrows
        for i in range(0, len(g_fine), 100):
            dx = x_vals[min(i+1, len(x_vals)-1)] - x_vals[max(i-1, 0)]
            dv = v_vals[min(i+1, len(v_vals)-1)] - v_vals[max(i-1, 0)]
            if abs(dx) + abs(dv) > 0:
                ax.annotate('', xy=(x_vals[i] + dx*0.3, v_vals[i] + dv*0.3),
                           xytext=(x_vals[i], v_vals[i]),
                           arrowprops=dict(arrowstyle='->', color=color, alpha=0.4))

        # Generation points
        for gi, gname in zip(g_gens, PARTICLE_NAMES[sname]):
            xg = a*gi**2 + b*gi + c
            vg = 2*a*gi + b
            ax.plot(xg, vg, 'o', color=GEN_COLORS[gi-1], markersize=12,
                    markeredgecolor='white', markeredgewidth=1.5, zorder=5)
            ax.annotate(gname, (xg, vg), textcoords="offset points",
                        xytext=(8, 5), color=GEN_COLORS[gi-1], fontsize=10,
                        fontweight='bold')

        # v=0 line
        ax.axhline(0, color='#666', alpha=0.3, linewidth=0.5)

        flat_tag = " [FLAT]" if col == flat_k else ""
        ax.set_title(f'{sname}/{COORD_SHORT[col]}{flat_tag}',
                     fontsize=10, color=color)
        ax.grid(True, alpha=0.15)

        if row == 2:
            ax.set_xlabel(f'x_{COORD_SHORT[col]}', fontsize=10)
        if col == 0:
            ax.set_ylabel('v = dx/dg', fontsize=10)

plt.suptitle('Phase portraits (x, v) — each parabola traces a straight line in phase space',
             fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig(f'{OUT}/05_phase_portraits.png', dpi=150, bbox_inches='tight')
plt.close()
print("  [5/8] Phase portraits saved")

# ════════════════════════════════════════════════════════════════
# FIGURE 6: Invariant search — what's conserved along the rail?
# ════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

for idx, sname in enumerate(SECTORS):
    color = COLORS[sname]

    # Row 0: |x(g)|² and |v(g)|² and 2E = |v|² (kinetic energy)
    ax = axes[0][idx]
    x2 = np.array([np.dot(x_vec(sname, g), x_vec(sname, g)) for g in g_fine])
    v2 = np.array([np.dot(v_vec(sname, g), v_vec(sname, g)) for g in g_fine])
    xv = np.array([np.dot(x_vec(sname, g), v_vec(sname, g)) for g in g_fine])

    ax.plot(g_fine, x2, color='#ff6b6b', linewidth=1.5, label='|x|²', alpha=0.8)
    ax.plot(g_fine, v2, color='#4ecdc4', linewidth=1.5, label='|v|²', alpha=0.8)
    ax.plot(g_fine, xv, color='#ffe66d', linewidth=1.5, label='x·v', alpha=0.8)

    # Mark generations
    for gi in g_gens:
        ax.axvline(gi, color='#555', alpha=0.3, linestyle='--')

    ax.set_title(f'{sname.upper()} — Norms', color=color, fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.15)
    ax.set_xlabel('g')

    # Row 1: Candidate invariants
    ax = axes[1][idx]

    # E_k = v²/2 for each coordinate separately
    for k in range(3):
        a, b, c = coeffs[sname][k]
        vk = np.array([2*a*g + b for g in g_fine])
        xk = np.array([a*g**2 + b*g + c for g in g_fine])

        # "Harmonic energy" E = v² + ω²x² — search for ω that makes it flattest
        # For a parabola: v = 2ag+b, x = ag²+bg+c
        # v² = 4a²g² + 4abg + b²  (quartic in g? no, quadratic)
        # x² = a²g⁴ + ... (quartic in g)
        # so v² is quadratic, x² is quartic → no simple conservation with const ω

        label = f'{COORD_SHORT[k]}{"*" if k == FLAT_COORD[sname] else ""}'
        ax.plot(g_fine, vk**2, linewidth=1.5, label=f'v_{label}²', alpha=0.7)

    # Total |v|²
    v2_total = np.array([np.dot(v_vec(sname, g), v_vec(sname, g)) for g in g_fine])
    ax.plot(g_fine, v2_total, 'w--', linewidth=2, label='|v|² total', alpha=0.8)

    for gi in g_gens:
        ax.axvline(gi, color='#555', alpha=0.3, linestyle='--')

    ax.set_title(f'{sname.upper()} — v_k² per coordinate', color=color, fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.15)
    ax.set_xlabel('g')

plt.tight_layout()
plt.savefig(f'{OUT}/06_invariant_search.png', dpi=150, bbox_inches='tight')
plt.close()
print("  [6/8] Invariant search saved")

# ════════════════════════════════════════════════════════════════
# FIGURE 7: Cross-sector comparison — universal features
# ════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (a) Normalized speed: |v(g)| / |v(2)| — does gen 2 universally normalize?
ax = axes[0][0]
for sname in SECTORS:
    speeds = np.array([np.linalg.norm(v_vec(sname, g)) for g in g_fine])
    v2_norm = np.linalg.norm(v_vec(sname, 2.0))
    ax.plot(g_fine, speeds / v2_norm, color=COLORS[sname], linewidth=2, label=sname)
for gi in g_gens:
    ax.axvline(gi, color='#555', alpha=0.3, linestyle='--')
ax.axhline(1, color='#888', alpha=0.3, linestyle=':')
ax.set_xlabel('g'); ax.set_ylabel('|v(g)| / |v(2)|')
ax.set_title('Normalized speed (baseline = gen 2)', fontsize=11)
ax.legend(); ax.grid(True, alpha=0.15)

# (b) "Angular momentum" L = x × v  (3D cross product magnitude)
ax = axes[0][1]
for sname in SECTORS:
    L = np.array([np.linalg.norm(np.cross(x_vec(sname, g), v_vec(sname, g)))
                  for g in g_fine])
    ax.plot(g_fine, L, color=COLORS[sname], linewidth=2, label=sname)
    for gi, gname in zip(g_gens, PARTICLE_NAMES[sname]):
        Lg = np.linalg.norm(np.cross(x_vec(sname, gi), v_vec(sname, gi)))
        ax.plot(gi, Lg, 'o', color=COLORS[sname], markersize=8,
                markeredgecolor='white', zorder=5)
for gi in g_gens:
    ax.axvline(gi, color='#555', alpha=0.3, linestyle='--')
ax.set_xlabel('g'); ax.set_ylabel('|x × v|')
ax.set_title('"Angular momentum" |L| = |x × v|', fontsize=11)
ax.legend(); ax.grid(True, alpha=0.15)

# (c) Runge-Lenz like: v × L direction
ax = axes[1][0]
for sname in SECTORS:
    # Compute x·v / |x||v| (angle between position and velocity)
    angles = []
    for g in g_fine:
        x = x_vec(sname, g)
        v = v_vec(sname, g)
        nx, nv = np.linalg.norm(x), np.linalg.norm(v)
        if nx > 1e-10 and nv > 1e-10:
            angles.append(np.arccos(np.clip(np.dot(x, v) / (nx * nv), -1, 1)) * 180 / np.pi)
        else:
            angles.append(0)
    ax.plot(g_fine, angles, color=COLORS[sname], linewidth=2, label=sname)
for gi in g_gens:
    ax.axvline(gi, color='#555', alpha=0.3, linestyle='--')
ax.set_xlabel('g'); ax.set_ylabel('angle(x, v) [deg]')
ax.set_title('Angle between position and velocity', fontsize=11)
ax.legend(); ax.grid(True, alpha=0.15)

# (d) Virial-like ratio: T/V where T = |v|²/2, V = |x|²/2
ax = axes[1][1]
for sname in SECTORS:
    ratio = []
    for g in g_fine:
        x = x_vec(sname, g)
        v = v_vec(sname, g)
        x2 = np.dot(x, x)
        v2 = np.dot(v, v)
        if x2 > 1e-10:
            ratio.append(v2 / x2)
        else:
            ratio.append(np.nan)
    ax.plot(g_fine, ratio, color=COLORS[sname], linewidth=2, label=sname)
for gi in g_gens:
    ax.axvline(gi, color='#555', alpha=0.3, linestyle='--')
ax.set_xlabel('g'); ax.set_ylabel('|v|² / |x|²')
ax.set_title('Virial ratio |v|²/|x|²', fontsize=11)
ax.legend(); ax.grid(True, alpha=0.15)
ax.set_ylim(0, 10)

plt.tight_layout()
plt.savefig(f'{OUT}/07_cross_sector.png', dpi=150, bbox_inches='tight')
plt.close()
print("  [7/8] Cross-sector comparison saved")

# ════════════════════════════════════════════════════════════════
# FIGURE 8: The delta_b deviation — eccentricity landscape
# ════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# (a) Bar chart of delta_b values
ax = axes[0]
labels = []
db_vals = []
colors_bar = []
for sname in SECTORS:
    flat_k = FLAT_COORD[sname]
    for k in range(3):
        a, b, c = coeffs[sname][k]
        db = b + 4 * a
        labels.append(f'{sname[0]}/{COORD_SHORT[k]}')
        db_vals.append(db)
        if k == flat_k:
            colors_bar.append('#888')
        else:
            colors_bar.append(COLORS[sname])

bars = ax.bar(range(len(labels)), db_vals, color=colors_bar, edgecolor='white', linewidth=0.5)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=9, rotation=45)
ax.axhline(0, color='#666', linewidth=1)
ax.set_ylabel('δb = b + 4a = mean velocity')
ax.set_title('Deviation from kinetic optimum\n(grey = flat, forced by constraint)', fontsize=11)
ax.grid(True, alpha=0.15, axis='y')

# (b) delta_b matrix heatmap
ax = axes[1]
db_matrix = np.zeros((3, 3))
for i, sname in enumerate(SECTORS):
    for k in range(3):
        a, b, c = coeffs[sname][k]
        db_matrix[i, k] = b + 4 * a

im = ax.imshow(db_matrix, cmap='RdBu_r', aspect='auto', vmin=-4, vmax=4)
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(COORD_SHORT, fontsize=11)
ax.set_yticks([0, 1, 2])
ax.set_yticklabels(['lepton', 'up', 'down'], fontsize=10)
for i in range(3):
    for j in range(3):
        flat_k = FLAT_COORD[list(SECTORS.keys())[i]]
        txt = f'{db_matrix[i,j]:.2f}'
        if j == flat_k:
            txt += '\n(flat)'
        ax.text(j, i, txt, ha='center', va='center', fontsize=10,
                color='white' if abs(db_matrix[i,j]) > 2 else 'black',
                fontweight='bold')
plt.colorbar(im, ax=ax, shrink=0.8)
ax.set_title('δB matrix (deviation from kinetic optimum)', fontsize=11)

# (c) Eccentricity vs curvature
ax = axes[2]
for sname in SECTORS:
    flat_k = FLAT_COORD[sname]
    for k in range(3):
        a, b, c = coeffs[sname][k]
        db = b + 4 * a
        if abs(a) > 1e-10:
            e = abs(db) / abs(4 * a)
            marker = 's' if k == flat_k else 'o'
            ax.plot(abs(a), e, marker, color=COLORS[sname], markersize=12,
                    markeredgecolor='white', markeredgewidth=1.5, zorder=5)
            ax.annotate(f'{sname[0]}/{COORD_SHORT[k]}', (abs(a), e),
                        textcoords="offset points", xytext=(8, 5),
                        color=COLORS[sname], fontsize=9)

ax.axhline(0.25, color='#888', alpha=0.5, linestyle='--', label='e=1/4 (flat universal)')
ax.set_xlabel('|a| (curvature)', fontsize=11)
ax.set_ylabel('e = |δb|/|4a| (eccentricity)', fontsize=11)
ax.set_title('Orbital eccentricity vs curvature', fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.15)

plt.tight_layout()
plt.savefig(f'{OUT}/08_eccentricity.png', dpi=150, bbox_inches='tight')
plt.close()
print("  [8/8] Eccentricity landscape saved")

print(f"\n  All plots saved to {OUT}/")

# ════════════════════════════════════════════════════════════════
# NUMERICAL INVARIANT SEARCH
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("  INVARIANT SEARCH — quantities at g=1,2,3")
print("=" * 70)

for sname in SECTORS:
    print(f"\n  {sname.upper()}:")
    for gi, gname in zip(g_gens, PARTICLE_NAMES[sname]):
        x = x_vec(sname, gi)
        v = v_vec(sname, gi)
        L = np.cross(x, v)
        print(f"    g={gi} ({gname:>3}): "
              f"|x|²={np.dot(x,x):8.3f}  "
              f"|v|²={np.dot(v,v):8.3f}  "
              f"x·v={np.dot(x,v):8.3f}  "
              f"|L|={np.linalg.norm(L):8.3f}  "
              f"|v|²/|x|²={'inf' if np.dot(x,x) < 1e-10 else f'{np.dot(v,v)/np.dot(x,x):8.3f}':>8}")

# Check: is |L| constant?
print("\n  |x × v| at generations:")
for sname in SECTORS:
    Ls = [np.linalg.norm(np.cross(x_vec(sname, g), v_vec(sname, g))) for g in g_gens]
    print(f"    {sname}: {Ls[0]:.4f}, {Ls[1]:.4f}, {Ls[2]:.4f}  "
          f"(var = {np.var(Ls):.4f}, mean = {np.mean(Ls):.4f})")

# Check: sum of |L| across sectors?
print("\n  Sum |L| across sectors at each generation:")
for gi, gname in zip(g_gens, ['gen1', 'gen2', 'gen3']):
    total_L = sum(np.linalg.norm(np.cross(x_vec(s, gi), v_vec(s, gi)))
                  for s in SECTORS)
    print(f"    g={gi}: Σ|L| = {total_L:.6f}")

# Check: |x(1)|² + |x(3)|² vs 2|x(2)|² (midpoint rule)
print("\n  Midpoint test: |x(1)|² + |x(3)|² vs 2|x(2)|²")
for sname in SECTORS:
    x1sq = np.dot(x_vec(sname, 1), x_vec(sname, 1))
    x2sq = np.dot(x_vec(sname, 2), x_vec(sname, 2))
    x3sq = np.dot(x_vec(sname, 3), x_vec(sname, 3))
    print(f"    {sname}: |x1|²+|x3|² = {x1sq + x3sq:.4f},  2|x2|² = {2*x2sq:.4f},  "
          f"diff = {x1sq + x3sq - 2*x2sq:.4f}")

# Check: x(1) + x(3) vs 2*x(2) (vectorial)
print("\n  Vector midpoint: x(1) + x(3) - 2x(2)")
for sname in SECTORS:
    diff = x_vec(sname, 1) + x_vec(sname, 3) - 2 * x_vec(sname, 2)
    print(f"    {sname}: ({diff[0]:+.4f}, {diff[1]:+.4f}, {diff[2]:+.4f})  "
          f"|diff| = {np.linalg.norm(diff):.4f}")
    # This should be 2a for each coordinate
    a_check = np.array([coeffs[sname][k][0] for k in range(3)])
    print(f"           2a = ({2*a_check[0]:+.4f}, {2*a_check[1]:+.4f}, {2*a_check[2]:+.4f})  "
          f"match: {np.allclose(diff, 2*a_check)}")

# The KEY: what combination is actually preserved?
print("\n  Discrete derivative test: Δx(g) = x(g+1) - x(g)")
for sname in SECTORS:
    dx12 = x_vec(sname, 2) - x_vec(sname, 1)  # = 3a + b
    dx23 = x_vec(sname, 3) - x_vec(sname, 2)  # = 5a + b
    ddx = dx23 - dx12  # = 2a (constant curvature)
    print(f"    {sname}:")
    print(f"      Δx(1→2) = ({dx12[0]:+.3f}, {dx12[1]:+.3f}, {dx12[2]:+.3f})")
    print(f"      Δx(2→3) = ({dx23[0]:+.3f}, {dx23[1]:+.3f}, {dx23[2]:+.3f})")
    print(f"      Δ²x     = ({ddx[0]:+.3f}, {ddx[1]:+.3f}, {ddx[2]:+.3f}) = 2a (INVARIANT)")
    # Cross product of successive deltas
    Ldisc = np.cross(dx12, dx23)
    print(f"      Δx12 × Δx23 = ({Ldisc[0]:+.3f}, {Ldisc[1]:+.3f}, {Ldisc[2]:+.3f})  "
          f"|L_disc| = {np.linalg.norm(Ldisc):.4f}")
