"""Quick F4 blind search with gamma=50."""
import sys, time
sys.path.insert(0, r'd:\Desarrollo\projects\pyshics\ss_peg\verification')
from test_closure_driven import evolve_with_closure, analyze_result

print('F4 blind: 52 gens, 27x27, gamma=50, 30k steps, 3 seeds')
t0 = time.time()
results, _ = evolve_with_closure(
    n_gen=52, dim=27, steps=30000, lr=0.001,
    alpha=5.0, beta=1.0, gamma=50.0,
    killing_target=-1.0, batch_size=3,
    seeds=[42, 137, 256], verbose=True,
    param_type='antisymmetric',
)
elapsed = time.time() - t0
print(f'\nElapsed: {elapsed:.1f}s')

for b, (gens_np, hist) in enumerate(results):
    a = analyze_result(list(gens_np), name=f'F4 gamma=50 seed={b}')
    tcl = hist['T_closure'][-1]
    print(f'  T_closure final: {tcl:.6f}')
