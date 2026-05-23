"""Deep dive into the rebel seed 100 and the simplicity band structure."""
import numpy as np
import json

data = np.load('landscape_128seeds_results.npz', allow_pickle=True)
results = json.loads(str(data['results']))

# ── Simplicity band analysis ──
# For dim=63 generators, if the K-spectrum splits into two degenerate blocks
# of sizes k and (63-k), simplicity = sqrt(k*(63-k)) / 63 * |delta|/|mean|
# The "pure" geometric factor sqrt(k*(63-k))/63 yields:
print("=" * 60)
print("TWO-BLOCK SPLITTING: sqrt(k*(63-k))/63 for various k")
print("=" * 60)
n = 63
for k in range(1, 32):
    geo = np.sqrt(k * (n - k)) / n
    print(f"  k={k:2d}  (63-k)={n-k:2d}  ->  geo = {geo:.4f}")

print()
print("Observed simplicity bands vs geometric factor:")
observed_bands = [0.0800, 0.1400, 0.2000, 0.2474, 0.2887]
for obs in observed_bands:
    # Find closest k
    best_k = min(range(1, 32), key=lambda k: abs(np.sqrt(k*(n-k))/n - obs))
    geo = np.sqrt(best_k * (n - best_k)) / n
    print(f"  S_obs={obs:.4f}  ->  k_best={best_k}  (geo={geo:.4f}, diff={abs(geo-obs):.4f})")

print()

# ── Alternative: ratio-based interpretation ──
print("=" * 60)
print("RATIO INTERPRETATION: S = 1/sqrt(N)")
print("=" * 60)
for obs in observed_bands:
    N_equiv = 1.0 / (obs * obs)
    print(f"  S={obs:.4f}  ->  N = 1/S^2 = {N_equiv:.1f}")

print()

# ── Rebel seed 100 deep analysis ──
print("=" * 60)
print("REBEL SEED 100 vs TYPICAL CLOSED_NON_SIMPLE")
print("=" * 60)

rebel = [r for r in results if r['seed'] == 100][0]
typicals = [r for r in results if r['category'] == 'CLOSED_NON_SIMPLE' and r['h_dual_ir3'] < 4]

print(f"\nRebel (seed 100):")
for key in sorted(rebel.keys()):
    print(f"  {key:20s} = {rebel[key]}")

print(f"\nTypical (mean of {len(typicals)} seeds):")
numeric_keys = [k for k in rebel.keys() if isinstance(rebel[k], (int, float)) and k != 'seed']
for key in sorted(numeric_keys):
    vals = [r[key] for r in typicals]
    rebel_val = rebel[key]
    ratio = rebel_val / np.mean(vals) if np.mean(vals) != 0 else float('inf')
    print(f"  {key:20s} = {np.mean(vals):12.6f} +/- {np.std(vals):.6f}   rebel/typ = {ratio:.2f}x")

print()

# ── Check: is seed 100 h_raw anomalous or is hv(IR=3) the anomaly? ──
print("=" * 60)
print("hv = h_raw * IR where IR=3 for the standard index")
print("=" * 60)
print(f"Rebel:   h_raw = {rebel['h_raw']:.6f}, hv = {rebel['h_dual_ir3']:.6f}, ratio = {rebel['h_dual_ir3']/rebel['h_raw']:.4f}")
typ_hraw = np.mean([r['h_raw'] for r in typicals])
typ_hv = np.mean([r['h_dual_ir3'] for r in typicals])
print(f"Typical: h_raw = {typ_hraw:.6f}, hv = {typ_hv:.6f}, ratio = {typ_hv/typ_hraw:.4f}")
print(f"\nRebel h_raw / typical h_raw = {rebel['h_raw']/typ_hraw:.4f}")
print(f"Rebel hv    / typical hv    = {rebel['h_dual_ir3']/typ_hv:.4f}")

print()

# ── All seeds sorted by hv ──
print("=" * 60)
print("TOP 10 SEEDS BY hv(IR=3) (closed or partial)")
print("=" * 60)
closed_partial = [r for r in results if r['closure_frac'] > 0.5]
for r in sorted(closed_partial, key=lambda x: -x['h_dual_ir3'])[:10]:
    print(f"  seed={r['seed']:3d}  cat={r['category']:20s}  "
          f"closure={r['closure_frac']:.4f}  S={r['simplicity']:.4f}  "
          f"hv={r['h_dual_ir3']:.4f}  h_raw={r['h_raw']:.4f}  "
          f"K_mean={r['K_mean']:.4f}")
