"""10-seed diagnostic at tau_Q=3000 to estimate P(frustration)."""
import sys
sys.path.insert(0, 'verification')

import numpy as np
import kibble_zurek_algebraic as kz

# Override config
cfg = kz.QUICK_CONFIG.copy()
cfg['TAU_Q_VALUES'] = [3000]
cfg['N_SEEDS'] = 10

results = kz.run_experiment(cfg)

print('\n' + '='*60)
print('  10-SEED FRUSTRATION DIAGNOSTIC  at tau_Q=3000')
print('='*60)
for r in results:
    s = r['mean_domain_size']
    nd = r['n_domains']
    rho = r['defect_density']
    d = r['mean_fp_distance']
    cl = r['mean_closure']
    tag = 'FRUSTRATED' if nd > 100 else ('SATURATED' if nd <= 2 else 'COARSENING')
    seed = r['seed']
    print(f'  seed={seed}  dom={nd:4d}  <s>={s:7.1f}  rho_d={rho:.3f}  <d>={d:.4f}  cl={cl:.3f}  [{tag}]')

domains = [r['n_domains'] for r in results]
frustrated = sum(1 for d in domains if d > 100)
coarsening = sum(1 for d in domains if 2 < d <= 100)
saturated = sum(1 for d in domains if d <= 2)
print(f'\n  P(frustrated) = {frustrated}/10 = {frustrated/10:.0%}')
print(f'  P(coarsening) = {coarsening}/10 = {coarsening/10:.0%}')
print(f'  P(saturated)  = {saturated}/10 = {saturated/10:.0%}')
print(f'  domains: mean={np.mean(domains):.1f}, std={np.std(domains):.1f}, min={min(domains)}, max={max(domains)}')
