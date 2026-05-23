import json, numpy as np, sys

path = sys.argv[1] if len(sys.argv) > 1 else r'd:\Desarrollo\projects\pyshics\ss_peg\verification\landscape_control_short.npz'
data = np.load(path, allow_pickle=True)

print(f"Keys: {list(data.keys())}")
results = json.loads(str(data['results']))
config = json.loads(str(data['config']))

print(f"\nConfig: {json.dumps(config, indent=2)}")
print(f"\nSeeds done: {data['n_seeds_done']}, total: {data['n_seeds_total']}")

if 'generators' in data:
    print(f"Generators shape: {data['generators'].shape}")

print(f"\n{'='*80}")
print(f"  RESULTS")
print(f"{'='*80}")
header = f"{'Seed':>4s}  {'Category':20s}  {'Close%':>7s}  {'MaxRes':>9s}  {'Simple':>8s}  {'h_v(3)':>7s}  {'T_total':>8s}  {'T_phys':>8s}  {'T_cl':>10s}"
print(header)
print("-" * len(header))
for r in results:
    print(f"{r['seed']:4d}  {r['category']:20s}  {r['closure_frac']:6.1%}  "
          f"{r['max_residual']:9.2e}  {r['simplicity']:8.4f}  "
          f"{r['h_dual_ir3']:7.2f}  {r['T_total']:8.1f}  "
          f"{r['T_phys']:8.1f}  {r['T_cl']:10.2e}")

# Summary
cats = {}
for r in results:
    cats.setdefault(r['category'], []).append(r)
print(f"\n  CATEGORIES:")
for cat, rs in sorted(cats.items()):
    print(f"    {cat}: {len(rs)}/{len(results)}")

if 'stability' in data:
    stab = json.loads(str(data['stability']))
    print(f"\n  STABILITY TEST:")
    for sr in stab:
        rec = "YES" if sr['category'] == 'F4' else "no"
        print(f"    eps={sr['epsilon']:.0e}: {sr['category']:15s} "
              f"closure={sr['closure_frac']:.1%} simple={sr['simplicity']:.4f} "
              f"h_v(3)={sr['h_dual_ir3']:.2f} recovered={rec}")
