"""Quick test: degree-corrected RH using d_max instead of d_mean."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from ihara_zeta_adjoint import *

catalog = build_ihara_catalog()
print("Degree-corrected RH: use d_max instead of d_mean")
print("  " + "=" * 95)
hdr = ("  {:>10s}  {:>4s}  {:>5s}  {:>6s}  {:>15s}  {:>10s}  "
       "{:>5s}  {:>4s}  {:>5s}")
print(hdr.format("Algebra", "dim", "d_max", "d_mean",
                  "1/sqrt(d_max-1)", "min|u_nt|", "dcRH?", "Ram?", "equiv"))
print("  " + "-" * 95)

n_equiv = 0
n_total = 0
for name, rank, all_roots, simple in catalog:
    adj = adjacency_from_roots(rank, all_roots, simple)
    degrees = adj.sum(axis=1)
    d_max = int(np.max(degrees))
    d_mean = np.mean(degrees)
    poles = ihara_poles_bass(adj)
    spec = spectral_analysis(adj)

    # All finite poles excluding near-zero
    radii = np.abs(poles[np.isfinite(poles)])
    radii = radii[radii > 0.001]
    # Exclude trivial poles near 1.0
    radii_nt = radii[np.abs(radii - 1.0) > 0.01]
    if len(radii_nt) == 0:
        continue

    min_nt = np.min(radii_nt)
    rc_dmax = 1.0 / np.sqrt(d_max - 1) if d_max > 1 else float('inf')
    dc_rh = min_nt >= rc_dmax - 1e-6
    ram = spec['is_ramanujan']
    equiv = "Y" if dc_rh == ram else "N"
    if dc_rh == ram:
        n_equiv += 1
    n_total += 1

    print("  {:>10s}  {:4d}  {:5d}  {:6.1f}  {:15.4f}  {:10.4f}  "
          "{:>5s}  {:>4s}  {:>5s}".format(
              name, adj.shape[0], d_max, d_mean,
              rc_dmax, min_nt,
              "YES" if dc_rh else "NO",
              "YES" if ram else "NO", equiv))

print()
print("  Degree-corrected RH equivalence: {}/{}".format(n_equiv, n_total))

# Also test: use d_max for inner boundary, d_min for outer
print()
print("  " + "=" * 95)
print("  Alternative: vertex-weighted Ramanujan circle r = 1/sqrt(q_eff)")
print("  Testing q_eff = sqrt(d_max * d_mean) - 1 (geometric mean correction)")
print("  " + "-" * 95)

n_equiv2 = 0
for name, rank, all_roots, simple in catalog:
    adj = adjacency_from_roots(rank, all_roots, simple)
    degrees = adj.sum(axis=1)
    d_max = int(np.max(degrees))
    d_mean = np.mean(degrees)
    q_eff = np.sqrt(d_max * d_mean) - 1
    poles = ihara_poles_bass(adj)
    spec = spectral_analysis(adj)

    radii = np.abs(poles[np.isfinite(poles)])
    radii = radii[radii > 0.001]
    radii_nt = radii[np.abs(radii - 1.0) > 0.01]
    if len(radii_nt) == 0:
        continue

    min_nt = np.min(radii_nt)
    rc_eff = 1.0 / np.sqrt(q_eff) if q_eff > 0 else float('inf')
    rh_eff = min_nt >= rc_eff - 1e-6
    ram = spec['is_ramanujan']
    equiv = "Y" if rh_eff == ram else "N"
    if rh_eff == ram:
        n_equiv2 += 1

    print("  {:>10s}  q_eff={:6.2f}  1/sqrt(q_eff)={:7.4f}  "
          "min|u_nt|={:7.4f}  effRH={:>3s}  Ram={:>3s}  {}".format(
              name, q_eff, rc_eff, min_nt,
              "YES" if rh_eff else "NO",
              "YES" if ram else "NO", equiv))

print()
print("  Geometric-mean RH equivalence: {}/{}".format(n_equiv2, n_total))
