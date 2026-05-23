"""Quick inspection of best F4 seeds from the 128-campaign."""
import numpy as np, json

d = np.load('landscape_128seeds_results.npz', allow_pickle=True)
results = json.loads(str(d['results']))
gens = d['generators']  # (128, 52, 27, 27) complex128

print(f"Generators shape: {gens.shape}")
print(f"Total seeds: {len(results)}")
print()

# Sort by closure descending, then T_cl ascending
seeds = sorted(results, key=lambda r: (-r['closure_frac'], r['T_cl']))

# Show top 20
print("Top 20 seeds by closure:")
print(f"  {'seed':>4s}  {'closure':>7s}  {'T_cl':>10s}  {'K_std':>8s}  {'S':>7s}  {'h_v(3)':>7s}  category")
print("  " + "-"*70)
for s in seeds[:20]:
    cl_pct = s['closure_frac'] * 100
    print(f"  {s['seed']:4d}  {cl_pct:6.1f}%  {s['T_cl']:10.2e}  {s['K_std']:8.4f}  {s['simplicity']:7.4f}  {s['h_dual_ir3']:7.2f}  {s['category']}")

print()

# Count categories
from collections import Counter
cats = Counter(r['category'] for r in results)
print("Category distribution:")
for cat, n in cats.most_common():
    print(f"  {cat}: {n}")

print()

# Best candidates: CLOSED or CLOSED_NON_SIMPLE with lowest T_cl
best = [r for r in results if r['category'] in ('CLOSED', 'CLOSED_NON_SIMPLE')]
best.sort(key=lambda r: r['T_cl'])
print(f"Best {min(10, len(best))} by T_cl (among CLOSED + CLOSED_NON_SIMPLE):")
for s in best[:10]:
    cl_pct = s['closure_frac'] * 100
    print(f"  seed {s['seed']:3d}: closure={cl_pct:5.1f}%  T_cl={s['T_cl']:.2e}  K_std={s['K_std']:.4f}  S={s['simplicity']:.4f}  h_v(3)={s['h_dual_ir3']:.2f}")
