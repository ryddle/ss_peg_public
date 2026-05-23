"""Analyze simplicity bands and rebel seeds in the 128-seed landscape."""
import numpy as np
import json

data = np.load('landscape_128seeds_results.npz', allow_pickle=True)
results = json.loads(str(data['results']))

# --- Seeds with closure > 50% ---
closed = [r for r in results if r['closure_frac'] > 0.5]
print(f"Seeds con closure > 50%: {len(closed)}")
print()

# --- Discrete simplicity values (horizontal bands) ---
simps_raw = [round(r['simplicity'], 4) for r in closed]
# Group into bands of width 0.005
bands = {}
for r in closed:
    s = r['simplicity']
    band_key = round(s * 200) / 200  # bin to 0.005
    bands.setdefault(band_key, []).append(r)

print(f"Simplicity bands (binned to 0.005):")
for bk in sorted(bands.keys()):
    seeds_in = bands[bk]
    if len(seeds_in) >= 2:
        seed_ids = [r['seed'] for r in seeds_in]
        cats = [r['category'] for r in seeds_in]
        print(f"  S ~ {bk:.4f}  ({len(seeds_in):2d} seeds)  cats: {dict((c, cats.count(c)) for c in set(cats))}")

print()

# --- The unique simplicity values ---
print("Exact simplicity values (all closed/partial):")
simps_unique = sorted(set(round(r['simplicity'], 4) for r in closed))
for s in simps_unique:
    count = sum(1 for r in closed if abs(r['simplicity'] - s) < 0.0005)
    marker = " <== MAJOR BAND" if count >= 5 else ""
    print(f"  {s:.4f}  x{count}{marker}")

print()

# --- Rebel seed: hv > 6 ---
print("=" * 60)
print("REBEL SEEDS (hv_IR3 > 5, closure > 0.90)")
print("=" * 60)
rebels = [r for r in results if r['h_dual_ir3'] > 5 and r['closure_frac'] > 0.90]
for r in sorted(rebels, key=lambda x: -x['h_dual_ir3']):
    print(f"\n  Seed {r['seed']}:")
    print(f"    category     = {r['category']}")
    print(f"    closure_frac = {r['closure_frac']:.6f}")
    print(f"    simplicity   = {r['simplicity']:.6f}")
    print(f"    hv(IR=3)     = {r['h_dual_ir3']:.4f}")
    print(f"    h_raw        = {r['h_raw']:.4f}")
    print(f"    K_mean       = {r['K_mean']:.6f}")
    print(f"    K_std        = {r['K_std']:.6f}")
    print(f"    max_residual = {r['max_residual']:.2e}")
    print(f"    T_total      = {r['T_total']:.4f}")

print()

# --- Compare: typical closed-non-simple vs rebel ---
print("=" * 60)
print("COMPARISON: Typical CLOSED_NON_SIMPLE vs Rebel")
print("=" * 60)
typical = [r for r in results if r['category'] == 'CLOSED_NON_SIMPLE' and r['h_dual_ir3'] < 5]
if typical:
    print(f"\nTypical ({len(typical)} seeds):")
    print(f"  hv(IR=3): mean={np.mean([r['h_dual_ir3'] for r in typical]):.4f}, "
          f"std={np.std([r['h_dual_ir3'] for r in typical]):.4f}")
    print(f"  K_mean:   mean={np.mean([r['K_mean'] for r in typical]):.6f}, "
          f"std={np.std([r['K_mean'] for r in typical]):.6f}")
    print(f"  K_std:    mean={np.mean([r['K_std'] for r in typical]):.6f}")
    print(f"  h_raw:    mean={np.mean([r['h_raw'] for r in typical]):.4f}")

# --- What do the simplicity bands correspond to? ---
print()
print("=" * 60)
print("SIMPLICITY BAND INTERPRETATION")
print("=" * 60)
print("For d=8, dim=63: simplicity = K_std / |K_mean|")
print("Discrete values suggest the optimizer converges to")
print("a finite set of K-spectrum structures.")
print()

# Check if simplicity ~ 1/sqrt(N) for small integers
for s in simps_unique:
    if s > 0:
        inv_sq = 1.0 / (s * s)
        print(f"  S={s:.4f}  ->  1/S^2 = {inv_sq:.1f}  ->  1/S = {1/s:.2f}")
