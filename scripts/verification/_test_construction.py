"""Quick test: F4 and E6 construction."""
import sys, time
sys.path.insert(0, 'verification')
from octonion_jordan import *

print('=== F4 CONSTRUCTION TEST ===')
t0 = time.time()
gens = build_f4_generators()
t1 = time.time()
print(f'Time: {t1-t0:.1f}s')
print(f'Number of generators: {len(gens)}')
if len(gens) > 0:
    print(f'Generator shape: {gens[0].shape}')
    overlaps = []
    for i in range(min(5, len(gens))):
        for j in range(i+1, min(5, len(gens))):
            overlaps.append(abs(np.sum(gens[i] * gens[j])))
    print(f'Max overlap (first 5): {max(overlaps):.2e}')
    norms = [np.sqrt(np.sum(g*g)) for g in gens[:5]]
    print(f'Norms (first 5): {[f"{n:.6f}" for n in norms]}')

print()
print('=== E6 CONSTRUCTION TEST ===')
t0 = time.time()
gens_e6 = build_e6_generators()
t1 = time.time()
print(f'Time: {t1-t0:.1f}s')
print(f'Number of generators: {len(gens_e6)}')
if len(gens_e6) > 0:
    print(f'Generator shape: {gens_e6[0].shape}')
    overlaps = []
    for i in range(min(5, len(gens_e6))):
        for j in range(i+1, min(5, len(gens_e6))):
            overlaps.append(abs(np.sum(gens_e6[i] * gens_e6[j])))
    print(f'Max overlap (first 5): {max(overlaps):.2e}')
    norms = [np.sqrt(np.sum(g*g)) for g in gens_e6[:5]]
    print(f'Norms (first 5): {[f"{n:.6f}" for n in norms]}')
