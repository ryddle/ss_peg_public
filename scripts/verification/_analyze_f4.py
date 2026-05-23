"""Re-analyze F4 gamma=50 seed=2 with various tolerances."""
import sys, numpy as np
sys.path.insert(0, r'd:\Desarrollo\projects\pyshics\ss_peg\verification')
from test_closure_driven import evolve_with_closure, analyze_result
from shared_utils_gpu import killing_metric_np, check_in_span_np, commutator_np

# Re-run seed 2 (best result)
print('Re-running F4 gamma=50, seed=256 (best seed)...')
results, _ = evolve_with_closure(
    n_gen=52, dim=27, steps=30000, lr=0.001,
    alpha=5.0, beta=1.0, gamma=50.0,
    killing_target=-1.0, batch_size=1,
    seeds=[256], verbose=False,
    param_type='antisymmetric',
)
gens_np = list(results[0][0])
n = len(gens_np)

# Check closure with different tolerances
residuals = []
for i in range(n):
    for j in range(i+1, n):
        C = commutator_np(gens_np[i], gens_np[j])
        norm_c = np.linalg.norm(C)
        if norm_c < 1e-10:
            residuals.append(0.0)
        else:
            B = np.array([g.flatten() for g in gens_np]).T
            coeffs, _, _, _ = np.linalg.lstsq(B, C.flatten(), rcond=None)
            res = np.linalg.norm(C.flatten() - B @ coeffs)
            residuals.append(res)

residuals = np.array(residuals)
print(f'\nResidual stats:')
print(f'  Max: {residuals.max():.6e}')
print(f'  Mean: {residuals.mean():.6e}')
print(f'  Median: {np.median(residuals):.6e}')
print(f'  99th percentile: {np.percentile(residuals, 99):.6e}')

for tol in [1e-4, 5e-4, 1e-3, 5e-3, 1e-2]:
    closes = np.sum(residuals < tol)
    print(f'  tol={tol:.0e}: {closes}/{len(residuals)} = {closes/len(residuals):.1%}')

# Renormalize to ||T||=1 and check K_diag
norms = np.array([np.linalg.norm(g) for g in gens_np])
print(f'\nGenerator norms: mean={norms.mean():.6f}, std={norms.std():.6e}')
gens_renorm = [g / np.linalg.norm(g) for g in gens_np]
K = killing_metric_np(gens_renorm)
K_diag = np.diag(K)
print(f'K_diag (renorm): mean={K_diag.mean():.6f}, std={K_diag.std():.6e}')
Tr_sq = np.mean([np.trace(g @ g).real for g in gens_renorm])
h_dual = abs(K_diag.mean() / Tr_sq)
print(f'h_dual (renorm, I_R=1): {h_dual:.4f}')
print(f'h_dual (renorm, I_R=3): {3 * h_dual:.4f}  (expected F4: 9.0)')
