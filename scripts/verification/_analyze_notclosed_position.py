"""Quick analysis: NOT_CLOSED seeds position in PCA plane + band boundaries."""
import numpy as np, json

data = np.load('landscape_128seeds_results.npz', allow_pickle=True)
results = json.loads(str(data['results']))
grid = np.load('landscape_relief_grid.npz', allow_pickle=True)
proj = grid['proj_2d']  # (128, 2)

print('=== NOT_CLOSED seeds: position & simplicity ===')
print(f'{"Seed":>5} {"PC1":>8} {"PC2":>8} {"S":>8} {"T_cl":>10} {"closure":>8}')
print('-' * 55)
for i, r in enumerate(results):
    if r['category'] == 'NOT_CLOSED':
        print(f"{r['seed']:5d} {proj[i,0]:8.3f} {proj[i,1]:8.3f} "
              f"{r['simplicity']:8.4f} {r['T_cl']:10.4e} {r['closure_frac']:8.3f}")

print()

# Band centers
band_centers = [0.08, 0.14, 0.20, 0.247, 0.289]
band_names = ['B1', 'B2', 'B3', 'B4', 'B5']

print('=== Simplicity distribution by category ===')
for cat in ['CLOSED_NON_SIMPLE', 'PARTIAL', 'NOT_CLOSED']:
    seeds = [r for r in results if r['category'] == cat]
    if seeds:
        simps = [r['simplicity'] for r in seeds]
        print(f"  {cat:20s}: n={len(seeds):3d}  S_mean={np.mean(simps):.4f}  "
              f"S_std={np.std(simps):.4f}  range=[{min(simps):.4f}, {max(simps):.4f}]")

print()

# Check inter-band gaps for NOT_CLOSED
print('=== NOT_CLOSED: distance to nearest band center ===')
for r in results:
    if r['category'] == 'NOT_CLOSED':
        s = r['simplicity']
        dists = [(abs(s - bc), bn) for bc, bn in zip(band_centers, band_names)]
        dists.sort()
        print(f"  seed {r['seed']:3d}: S={s:.4f}  "
              f"nearest={dists[0][1]}(d={dists[0][0]:.4f})  "
              f"2nd={dists[1][1]}(d={dists[1][0]:.4f})")

print()

# Angular analysis from centroid
centroid = proj.mean(axis=0)
print(f'Centroid: ({centroid[0]:.3f}, {centroid[1]:.3f})')
print()

print('=== Angular position from centroid (all seeds) ===')
for cat in ['NOT_CLOSED', 'PARTIAL', 'CLOSED_NON_SIMPLE']:
    idxs = [i for i, r in enumerate(results) if r['category'] == cat]
    if not idxs:
        continue
    angles = []
    dists = []
    for i in idxs:
        dx = proj[i, 0] - centroid[0]
        dy = proj[i, 1] - centroid[1]
        angles.append(np.degrees(np.arctan2(dy, dx)))
        dists.append(np.sqrt(dx**2 + dy**2))

    print(f"\n  {cat} ({len(idxs)} seeds):")
    order = np.argsort(angles)
    for j in order:
        i = idxs[j]
        r = results[i]
        print(f"    seed {r['seed']:3d}: angle={angles[j]:7.1f}°  "
              f"dist={dists[j]:.3f}  S={r['simplicity']:.4f}  "
              f"band={'B'+str(1+min(range(5), key=lambda k: abs(r['simplicity']-band_centers[k])))}")

# Summary: angular concentration of NOT_CLOSED
print()
print('=== Summary: angular sectors ===')
nc_angles = []
for i, r in enumerate(results):
    if r['category'] == 'NOT_CLOSED':
        dx = proj[i, 0] - centroid[0]
        dy = proj[i, 1] - centroid[1]
        nc_angles.append(np.degrees(np.arctan2(dy, dx)))

nc_angles = np.array(nc_angles)
print(f"NOT_CLOSED angles: {np.sort(nc_angles)}")
print(f"  Mean: {nc_angles.mean():.1f}°  Std: {nc_angles.std():.1f}°")

# Check if they cluster in ~2 diagonal directions
from collections import Counter
quadrants = Counter()
for a in nc_angles:
    if -45 <= a < 45:
        quadrants['East (0°)'] += 1
    elif 45 <= a < 135:
        quadrants['North (90°)'] += 1
    elif -135 <= a < -45:
        quadrants['South (-90°)'] += 1
    else:
        quadrants['West (180°)'] += 1
print(f"  Quadrants: {dict(quadrants)}")
