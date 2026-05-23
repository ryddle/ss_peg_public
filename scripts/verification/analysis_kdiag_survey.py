"""
K_diag Survey: How does the converged Killing diagonal vary across all algebras?

Runs all algebras with uniform optimizer config (where applicable) and collects:
  - K eigenvalues (sorted)
  - K_diag mean, std
  - K_offdiag max
  - T_total final
  - Cross with (n_gen, dim, ρ = n_gen/(dim²-1), h∨)

Categories:
  A. OPTIMIZED compact: SU(2–5), SO(3–6), Sp(2,4,6) — target K = -I
  B. OPTIMIZED non-compact: sl(2,ℝ), so(3,1) — indefinite targets
  C. CONSTRUCTIVE: G₂, F₄, E₆ — prebuilt generators, just measure K

Usage:
    python analysis_kdiag_survey.py
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import torch

from shared_utils_gpu import (
    evolve_generators_gpu,
    evolve_generators_generic,
    killing_metric_np,
    check_in_span_np,
    commutator_np,
    get_device_info,
    params_to_generators,
    generators_to_params,
    random_params,
    params_to_antisymmetric,
    antisymmetric_to_params,
    random_params_so,
    params_to_symplectic,
    symplectic_to_params,
    random_params_sp,
    DTYPE, FDTYPE, DEVICE,
    batched_killing_metric,
    batched_non_commutativity_tension,
    batched_norm_tension,
)

# Non-compact parameterizations
from test_sl2r_gpu import (
    params_to_sl2r, sl2r_to_params, random_params_sl2r,
    evolve_sl2r_indefinite,
)
from test_so31_gpu import (
    params_to_so31, so31_to_params, random_params_so31,
    evolve_so31_indefinite,
)

# Constructive generators
from test_g2_gpu import build_g2_generators
from octonion_jordan import build_f4_generators, build_e6_generators


# ===================================================================
# Uniform optimizer config
# ===================================================================
STEPS = 12000
LR = 0.001
ALPHA = 5.0
BETA = 1.0
N_SEEDS = 3
SEEDS = [42, 137, 256]

# Max n_gen for full closure check (skip for large algebras — too slow)
CLOSURE_MAX_NGEN = 24


def measure_killing(generators_np, check_closure=True):
    """Compute Killing metric properties from numpy generators."""
    n_gen = len(generators_np)
    K = killing_metric_np(generators_np)
    K_diag = np.diag(K).real
    K_offdiag = np.abs(K - np.diag(K_diag)).max()
    K_eigs = np.sort(np.linalg.eigvalsh(K.real))

    # Signature
    n_pos = int(np.sum(K_eigs > 1e-6))
    n_neg = int(np.sum(K_eigs < -1e-6))
    n_zero = n_gen - n_pos - n_neg

    # Closure (sample for large algebras)
    closure = -1.0  # sentinel for "not checked"
    if check_closure and n_gen <= CLOSURE_MAX_NGEN:
        closes = 0
        total_pairs = n_gen * (n_gen - 1) // 2
        for i in range(n_gen):
            for j in range(i + 1, n_gen):
                C = commutator_np(generators_np[i], generators_np[j])
                if np.linalg.norm(C) < 1e-10:
                    closes += 1
                elif check_in_span_np(C, generators_np, tol=1e-4):
                    closes += 1
        closure = closes / total_pairs if total_pairs > 0 else 1.0

    # Frobenius norms
    frob_norms = [np.linalg.norm(g) for g in generators_np]

    return {
        "K_eigs": K_eigs,
        "K_diag_mean": np.mean(K_diag),
        "K_diag_std": np.std(K_diag),
        "K_diag_min": np.min(K_diag),
        "K_diag_max": np.max(K_diag),
        "K_offdiag": K_offdiag,
        "signature": (n_pos, n_zero, n_neg),
        "closure": closure,
        "frob_mean": np.mean(frob_norms),
    }


# ===================================================================
# A. Compact algebras via optimizer
# ===================================================================

def run_su(n):
    """SU(n): n²-1 hermitian traceless generators in dim=n."""
    n_gen = n * n - 1
    dim = n
    results = evolve_generators_generic(
        n_gen=n_gen, dim=dim, steps=STEPS, lr=LR,
        alpha=ALPHA, beta=BETA, killing_target=-1.0,
        batch_size=N_SEEDS, seeds=SEEDS,
        verbose=False, convergence_threshold=1e-6,
    )
    # Pick seed with lowest final T_total
    best = min(results, key=lambda r: r[1]["T_total"][-1])
    return list(best[0])


def run_so(n):
    """SO(n): n(n-1)/2 antisymmetric generators in dim=n."""
    n_gen = n * (n - 1) // 2
    dim = n
    results = evolve_generators_generic(
        n_gen=n_gen, dim=dim, steps=STEPS, lr=LR,
        alpha=ALPHA, beta=BETA, killing_target=-1.0,
        batch_size=N_SEEDS, seeds=SEEDS,
        verbose=False, convergence_threshold=1e-6,
        param_to_gen_fn=lambda p, ng: params_to_antisymmetric(p, ng, dim),
        gen_to_param_fn=antisymmetric_to_params,
        random_param_fn=lambda bs, ng, device, seeds: random_params_so(bs, ng, dim, device, seeds),
        param_to_gen_kwargs={},
    )
    best = min(results, key=lambda r: r[1]["T_total"][-1])
    return list(best[0])


def run_sp(half_n):
    """Sp(2n): n(2n+1) generators in dim=2n."""
    n_gen = half_n * (2 * half_n + 1)
    dim = 2 * half_n
    results = evolve_generators_generic(
        n_gen=n_gen, dim=dim, steps=STEPS, lr=LR,
        alpha=ALPHA, beta=BETA, killing_target=-1.0,
        batch_size=N_SEEDS, seeds=SEEDS,
        verbose=False, convergence_threshold=1e-6,
        param_to_gen_fn=params_to_symplectic,
        gen_to_param_fn=lambda g: symplectic_to_params(g, half_n),
        random_param_fn=lambda bs, ng, device, seeds: random_params_sp(bs, ng, half_n, device, seeds),
        param_to_gen_kwargs={"half_dim": half_n},
    )
    best = min(results, key=lambda r: r[1]["T_total"][-1])
    return list(best[0])


# ===================================================================
# B. Non-compact algebras via custom optimizers
# ===================================================================

def run_sl2r():
    """sl(2,ℝ): 3 real traceless 2×2 generators, indefinite K target."""
    results = evolve_sl2r_indefinite(
        n_gen=3, steps=20000, lr=LR,
        alpha=ALPHA, beta=BETA,
        batch_size=N_SEEDS, seeds=SEEDS,
        verbose=False, convergence_threshold=1e-4,
    )
    best = min(results, key=lambda r: r[1]["loss"][-1])
    return list(best[0])


def run_so31():
    """so(3,1): 6 η-preserving 4×4 generators, indefinite K target."""
    results = evolve_so31_indefinite(
        n_gen=6, steps=30000, lr=LR,
        alpha=ALPHA, beta=BETA,
        batch_size=N_SEEDS, seeds=SEEDS,
        verbose=False, convergence_threshold=1e-4,
    )
    best = min(results, key=lambda r: r[1]["loss"][-1])
    return list(best[0])


# ===================================================================
# C. Constructive (prebuilt) algebras
# ===================================================================

def run_g2_constructive():
    """G₂ constructive: 14 generators in dim=7 from octonion 3-form."""
    return build_g2_generators()


def run_f4_constructive():
    """F₄ constructive: 52 generators in dim=27 from Jordan algebra."""
    return build_f4_generators()


def run_e6_constructive():
    """E₆ constructive: 78 generators in dim=27 from cubic norm."""
    return build_e6_generators()


# ===================================================================
# Known theoretical values
# ===================================================================
# h∨ (dual Coxeter number), I_R (Dynkin index of representation used)
THEORY = {
    "SU(2)":  {"n_gen": 3,  "dim": 2,  "h_dual": 2,  "I_R": 1, "compact": True},
    "SU(3)":  {"n_gen": 8,  "dim": 3,  "h_dual": 3,  "I_R": 1, "compact": True},
    "SU(4)":  {"n_gen": 15, "dim": 4,  "h_dual": 4,  "I_R": 1, "compact": True},
    "SU(5)":  {"n_gen": 24, "dim": 5,  "h_dual": 5,  "I_R": 1, "compact": True},
    "SO(3)":  {"n_gen": 3,  "dim": 3,  "h_dual": 2,  "I_R": 1, "compact": True},
    "SO(4)":  {"n_gen": 6,  "dim": 4,  "h_dual": 2,  "I_R": 1, "compact": True},
    "SO(5)":  {"n_gen": 10, "dim": 5,  "h_dual": 3,  "I_R": 1, "compact": True},
    "SO(6)":  {"n_gen": 15, "dim": 6,  "h_dual": 4,  "I_R": 1, "compact": True},
    "Sp(2)":  {"n_gen": 3,  "dim": 2,  "h_dual": 2,  "I_R": 1, "compact": True},
    "Sp(4)":  {"n_gen": 10, "dim": 4,  "h_dual": 3,  "I_R": 1, "compact": True},
    "Sp(6)":  {"n_gen": 21, "dim": 6,  "h_dual": 4,  "I_R": 1, "compact": True},
    "sl(2,R)":  {"n_gen": 3,  "dim": 2,  "h_dual": 2,  "I_R": 1, "compact": False},
    "so(3,1)":  {"n_gen": 6,  "dim": 4,  "h_dual": 2,  "I_R": 1, "compact": False},
    "G2":     {"n_gen": 14, "dim": 7,  "h_dual": 4,  "I_R": 1, "compact": True},
    "F4":     {"n_gen": 52, "dim": 27, "h_dual": 9,  "I_R": 3, "compact": True},
    "E6":     {"n_gen": 78, "dim": 27, "h_dual": 12, "I_R": 3, "compact": True},
}


# ===================================================================
# Main
# ===================================================================
def main():
    print("=" * 100)
    print("  K_diag SURVEY: Converged Killing diagonal across all tested algebras")
    print("=" * 100)
    print(f"  Device: {get_device_info()}")
    print(f"  Optimizer config: steps={STEPS}, lr={LR}, alpha={ALPHA}, beta={BETA}")
    print(f"  Seeds: {SEEDS} (best of {N_SEEDS})")
    print()

    all_results = {}

    # --- A. Compact (optimized) ---
    print("=" * 100)
    print("  CATEGORY A: Compact algebras (optimized, K_target = -I)")
    print("=" * 100)

    compact_runs = [
        ("SU(2)", lambda: run_su(2)),
        ("SU(3)", lambda: run_su(3)),
        ("SU(4)", lambda: run_su(4)),
        ("SU(5)", lambda: run_su(5)),
        ("SO(3)", lambda: run_so(3)),
        ("SO(4)", lambda: run_so(4)),
        ("SO(5)", lambda: run_so(5)),
        ("SO(6)", lambda: run_so(6)),
        ("Sp(2)", lambda: run_sp(1)),
        ("Sp(4)", lambda: run_sp(2)),
        ("Sp(6)", lambda: run_sp(3)),
    ]

    for name, runner in compact_runs:
        t0 = time.time()
        gens = runner()
        dt = time.time() - t0
        m = measure_killing(gens)
        m["time"] = dt
        all_results[name] = m
        cl_str = f"{m['closure']:.1%}" if m['closure'] >= 0 else "skip"
        print(f"  {name:8s}: K_diag_mean={m['K_diag_mean']:+.6f} "
              f"std={m['K_diag_std']:.2e} "
              f"offdiag={m['K_offdiag']:.2e} "
              f"closure={cl_str} "
              f"frob={m['frob_mean']:.4f} "
              f"({dt:.1f}s)")

    # --- B. Non-compact (optimized) ---
    print()
    print("=" * 100)
    print("  CATEGORY B: Non-compact algebras (optimized, indefinite K_target)")
    print("=" * 100)

    t0 = time.time()
    gens = run_sl2r()
    dt = time.time() - t0
    m = measure_killing(gens)
    m["time"] = dt
    all_results["sl(2,R)"] = m
    print(f"  sl(2,R):  K_eigs={m['K_eigs']} "
          f"offdiag={m['K_offdiag']:.2e} "
          f"closure={m['closure']:.1%} "
          f"({dt:.1f}s)")

    t0 = time.time()
    gens = run_so31()
    dt = time.time() - t0
    m = measure_killing(gens)
    m["time"] = dt
    all_results["so(3,1)"] = m
    print(f"  so(3,1):  K_eigs={m['K_eigs']} "
          f"offdiag={m['K_offdiag']:.2e} "
          f"closure={m['closure']:.1%} "
          f"({dt:.1f}s)")

    # --- C. Constructive ---
    print()
    print("=" * 100)
    print("  CATEGORY C: Constructive (prebuilt generators, no optimization)")
    print("=" * 100)

    constructive_runs = [
        ("G2", run_g2_constructive),
        ("F4", run_f4_constructive),
        ("E6", run_e6_constructive),
    ]

    for name, runner in constructive_runs:
        t0 = time.time()
        gens = runner()
        dt = time.time() - t0
        # Skip closure for large constructive algebras (F4=52, E6=78 gens)
        m = measure_killing(gens, check_closure=(len(gens) <= CLOSURE_MAX_NGEN))
        m["time"] = dt
        all_results[name] = m
        print(f"  {name:4s}: K_diag_mean={m['K_diag_mean']:+.6f} "
              f"std={m['K_diag_std']:.2e} "
              f"K_eigs=[{m['K_eigs'][0]:.4f}..{m['K_eigs'][-1]:.4f}] "
              f"closure={m['closure']:.1%} "
              f"({dt:.1f}s)")

    # ===================================================================
    # MASTER TABLE
    # ===================================================================
    print("\n\n" + "=" * 100)
    print("  MASTER TABLE: K_diag across all algebras")
    print("=" * 100)

    header = (f"  {'Algebra':10s} {'n_gen':>5s} {'dim':>4s} {'rho':>6s} "
              f"{'h_v':>4s} {'|K_eig|_min':>11s} {'|K_eig|_max':>11s} "
              f"{'K_diag_mean':>12s} {'K_diag_std':>11s} "
              f"{'K_off':>9s} {'sig':>10s} {'closure':>8s} {'frob':>6s}")
    print(header)
    print("  " + "-" * (len(header) - 2))

    for name in THEORY:
        if name not in all_results:
            continue
        r = all_results[name]
        th = THEORY[name]
        rho = th["n_gen"] / (th["dim"]**2 - 1) if th["dim"]**2 - 1 > 0 else 0
        eigs = r["K_eigs"]
        abs_eigs = np.abs(eigs)
        abs_eigs_nonzero = abs_eigs[abs_eigs > 1e-8]
        eig_min = np.min(abs_eigs_nonzero) if len(abs_eigs_nonzero) > 0 else 0
        eig_max = np.max(abs_eigs_nonzero) if len(abs_eigs_nonzero) > 0 else 0
        sig = r["signature"]

        cl_str = f"{r['closure']:>7.1%}" if r['closure'] >= 0 else "   skip"
        print(f"  {name:10s} {th['n_gen']:>5d} {th['dim']:>4d} {rho:>6.3f} "
              f"{th['h_dual']:>4d} {eig_min:>11.6f} {eig_max:>11.6f} "
              f"{r['K_diag_mean']:>+12.6f} {r['K_diag_std']:>11.2e} "
              f"{r['K_offdiag']:>9.2e} "
              f"({sig[0]}+,{sig[1]}0,{sig[2]}-) "
              f"{cl_str} {r['frob_mean']:>6.3f}")

    # ===================================================================
    # ANALYSIS: Search for patterns
    # ===================================================================
    print("\n\n" + "=" * 100)
    print("  ANALYSIS: K_diag patterns")
    print("=" * 100)

    # Only compact optimized algebras (where target is -1)
    compact_names = [n for n in THEORY if THEORY[n]["compact"] and n in all_results
                     and n not in ("G2", "F4", "E6")]

    print("\n  --- Compact optimized algebras (target K_diag = -1.0) ---")
    print(f"  {'Algebra':10s} {'n_gen':>5s} {'dim':>4s} {'rho':>6s} "
          f"{'K_diag_mean':>12s} {'ratio K/target':>15s} "
          f"{'h_v':>4s} {'frob':>6s}")
    print("  " + "-" * 65)

    for name in compact_names:
        r = all_results[name]
        th = THEORY[name]
        rho = th["n_gen"] / (th["dim"]**2 - 1)
        ratio = r["K_diag_mean"] / (-1.0)
        print(f"  {name:10s} {th['n_gen']:>5d} {th['dim']:>4d} {rho:>6.3f} "
              f"{r['K_diag_mean']:>+12.6f} {ratio:>15.6f} "
              f"{th['h_dual']:>4d} {r['frob_mean']:>6.3f}")

    # Check correlations
    print("\n  --- Correlation analysis ---")

    rhos = []
    k_means = []
    h_duals = []
    n_gens = []
    dims = []
    frobs = []

    for name in compact_names:
        r = all_results[name]
        th = THEORY[name]
        rho = th["n_gen"] / (th["dim"]**2 - 1)
        rhos.append(rho)
        k_means.append(abs(r["K_diag_mean"]))
        h_duals.append(th["h_dual"])
        n_gens.append(th["n_gen"])
        dims.append(th["dim"])
        frobs.append(r["frob_mean"])

    rhos = np.array(rhos)
    k_means = np.array(k_means)
    h_duals = np.array(h_duals)
    n_gens = np.array(n_gens)
    dims = np.array(dims)
    frobs = np.array(frobs)

    # Pearson correlations
    from numpy import corrcoef
    vars_to_check = [
        ("rho", rhos),
        ("h_dual", h_duals),
        ("n_gen", n_gens),
        ("dim", dims),
        ("frob_mean", frobs),
    ]

    for vname, v in vars_to_check:
        if np.std(v) > 0 and np.std(k_means) > 0:
            corr = corrcoef(v, k_means)[0, 1]
            print(f"  corr(|K_diag_mean|, {vname:12s}) = {corr:+.4f}")

    # Check if K_diag ~ -2*h_dual*I_R / dim  or similar scaling
    print("\n  --- Scaling hypotheses ---")
    for name in compact_names:
        r = all_results[name]
        th = THEORY[name]
        rho = th["n_gen"] / (th["dim"]**2 - 1)
        k_abs = abs(r["K_diag_mean"])
        # Hypothesis 1: K_diag = f(rho)
        # Hypothesis 2: K_diag = f(h_dual)
        # Hypothesis 3: K_diag = f(n_gen, dim)
        # Just print the ratios
        print(f"  {name:10s} |K|/h_v={k_abs/th['h_dual']:.6f}  "
              f"|K|*dim={k_abs*th['dim']:.6f}  "
              f"|K|*n_gen={k_abs*th['n_gen']:.6f}  "
              f"|K|/(h_v/dim)={k_abs/(th['h_dual']/th['dim']):.6f}")

    # Non-compact summary
    print("\n  --- Non-compact algebras ---")
    for name in ["sl(2,R)", "so(3,1)"]:
        if name in all_results:
            r = all_results[name]
            th = THEORY[name]
            eigs = r["K_eigs"]
            neg_eigs = eigs[eigs < -1e-6]
            pos_eigs = eigs[eigs > 1e-6]
            print(f"  {name:10s}: neg_mean={np.mean(neg_eigs):.6f} "
                  f"pos_mean={np.mean(pos_eigs):.6f} "
                  f"|neg|/|pos|={abs(np.mean(neg_eigs)/np.mean(pos_eigs)):.6f}")

    # Constructive summary
    print("\n  --- Constructive algebras ---")
    for name in ["G2", "F4", "E6"]:
        if name in all_results:
            r = all_results[name]
            th = THEORY[name]
            eigs = r["K_eigs"]
            print(f"  {name:4s}: K_eigs range [{eigs[0]:.4f}, {eigs[-1]:.4f}] "
                  f"mean={np.mean(eigs):.4f} "
                  f"h_v={th['h_dual']} "
                  f"|K_mean|/h_v={abs(np.mean(eigs))/th['h_dual']:.6f}")


if __name__ == "__main__":
    main()
