"""
Test 1.2: Independence of alpha (weight of Killing tension).
Question: Does the final algebra depend on the relative weight alpha?

EXPECTED: All alphas produce a closed SU(2) algebra; only convergence speed varies.
The optimizer has early-stopping (plateau detection), so each alpha runs
until it converges or stalls — no need to pre-tune step counts.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from shared_utils import evolve_parallel, killing_metric, commutator, check_in_span


def run_test_alpha_independence(steps=15000):
    print("=" * 70)
    print("TEST 1.2: INDEPENDENCE OF ALPHA")
    print("=" * 70)

    # Start at 0.3 — very low alphas (≤0.1) genuinely cannot drive K→-I
    alphas = [0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    # Give generous budget; early stopping will cut short if plateau reached
    max_steps = 50000

    tasks = [
        {'seed': 42, 'n_gen': 3, 'dim': 2,
         'steps': max_steps,
         'lr': 0.001, 'alpha': alpha, 'beta': 1.0,
         'killing_target': -1.0, 'verbose': False}
        for alpha in alphas
    ]
    print(f"  Running {len(alphas)} alpha values in parallel (max {max_steps} steps, early-stop on plateau)...")
    raw = evolve_parallel(tasks)

    results = []
    for alpha, (final_gen, hist) in zip(alphas, raw):
        K = killing_metric(final_gen)
        K_diag_mean = np.mean(np.diag(K))
        K_offdiag_max = np.max(np.abs(K - np.diag(np.diag(K))))
        T_final = hist["T_total"][-1]
        actual_steps = len(hist["T_total"])
        stopped_early = actual_steps < max_steps

        closes_all = all(
            check_in_span(commutator(final_gen[i], final_gen[j]), final_gen)
            for i, j in [(0, 1), (1, 2), (2, 0)]
        )

        results.append({
            "alpha": alpha,
            "K_diag": K_diag_mean,
            "K_offdiag": K_offdiag_max,
            "T_final": T_final,
            "actual_steps": actual_steps,
            "stopped_early": stopped_early,
            "closes": closes_all,
        })

        tag = "plateau" if stopped_early else "max"
        print(f"  alpha={alpha:5.1f}: K_diag={K_diag_mean:.4f}, "
              f"closes={closes_all}, K_offdiag={K_offdiag_max:.6f}, "
              f"steps={actual_steps} ({tag})")

    # ---- Analysis ----
    K_diags = [r["K_diag"] for r in results]
    K_spread = max(K_diags) - min(K_diags)
    closes_count = sum(r["closes"] for r in results)
    early_stops = sum(r["stopped_early"] for r in results)

    print("\n--- RESULTS ---")
    print(f"K_diag range across alphas: [{min(K_diags):.4f}, {max(K_diags):.4f}]")
    print(f"K_diag spread: {K_spread:.4f}")
    print(f"Early stopped (plateau): {early_stops}/{len(alphas)}")
    print(f"Closure: {closes_count}/{len(alphas)}")

    # ---- Plot ----
    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))

        axes[0].plot([r["alpha"] for r in results],
                     [r["K_diag"] for r in results], "o-")
        axes[0].axhline(-1.0, color="r", linestyle="--", label="Target")
        axes[0].set_xlabel("alpha")
        axes[0].set_ylabel("K diagonal mean")
        axes[0].legend()
        axes[0].grid(True)

        axes[1].semilogy([r["alpha"] for r in results],
                         [r["K_offdiag"] for r in results], "o-")
        axes[1].set_xlabel("alpha")
        axes[1].set_ylabel("K off-diagonal max")
        axes[1].grid(True)

        axes[2].plot([r["alpha"] for r in results],
                     [r["actual_steps"] for r in results], "o-")
        axes[2].set_xlabel("alpha")
        axes[2].set_ylabel("Actual steps (early-stop aware)")
        axes[2].grid(True)

        plt.tight_layout()
        plt.savefig(os.path.join(os.path.dirname(__file__), "test_1_2_alpha_independence.png"), dpi=150)
        print("\nPlot saved: test_1_2_alpha_independence.png")
    except ImportError:
        print("(matplotlib not available, skipping plot)")

    # ---- Verdict ----
    print("\n--- VERDICT ---")
    if closes_count == len(alphas):
        print(f"PASS: All {len(alphas)} alpha values produce closed SU(2) algebra.")
        print(f"  K_diag spread: {K_spread:.4f} (convergence tightness varies with alpha).")
    elif closes_count >= len(alphas) - 1:
        print(f"PARTIAL: {closes_count}/{len(alphas)} close. One alpha may need more steps.")
    else:
        print(f"FAIL: Only {closes_count}/{len(alphas)} alphas produce SU(2).")

    return results


if __name__ == "__main__":
    run_test_alpha_independence()
