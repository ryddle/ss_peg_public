"""
run_all_tests.py — Runs all verification tests from the systematic plan.

Usage:
    python run_all_tests.py          # run all tests
    python run_all_tests.py 1_1 2_3  # run only selected tests
"""
import sys, os, time, traceback, importlib

sys.path.insert(0, os.path.dirname(__file__))

# Registry: (module_name, run_function, description)
TESTS = [
    ("test_1_1_reproducibility",       "run_test_reproducibility",       "Phase 1 — Reproducibility (100 seeds)"),
    ("test_1_2_alpha_independence",     "run_test_alpha_independence",    "Phase 1 — Alpha weight independence"),
    ("test_1_3_noise_robustness",       "run_test_noise_robustness",     "Phase 1 — Noise robustness"),
    ("test_2_1_dim3",                   "run_test_dim3",                 "Phase 2 — dim=3 (spin-1 representation)"),
    ("test_2_2_two_generators",         "run_test_2_generators",         "Phase 2 — 2 generators (U(1) expected)"),
    ("test_2_3_su3_eight_generators",   "run_test_su3",                  "Phase 2 — 8 generators SU(3)"),
    ("test_3_1_only_killing",           "run_test_only_killing",         "Phase 3 — Only T_Killing (no T_nc)"),
    ("test_3_2_only_nc",                "run_test_only_nc",              "Phase 3 — Only T_nc (no Killing)"),
    ("test_3_3_energy_landscape",       "run_test_landscape",            "Phase 3 — Energy landscape scan"),
    ("test_4_1_representation_theory",  "run_test_representation_theory","Phase 4 — Casimir & representation theory"),
    # GPU-accelerated variants (require PyTorch + CUDA)
    ("test_1_1_gpu_reproducibility",    "run_test_reproducibility_gpu",  "GPU — Reproducibility (100 seeds batched)"),
    ("test_2_3_gpu_su3",                "run_test_su3_gpu",              "GPU — 8 generators SU(3) (autograd)"),
    ("test_3_3_gpu_energy_landscape",   "run_test_landscape_gpu",        "GPU — Energy landscape scan (batched)"),
    ("test_su_n_gpu",                   "run_test_su_n_sweep",           "GPU — SU(N) universality sweep"),
]


def run_all(selected=None):
    results = []
    total_t0 = time.time()

    print("=" * 70)
    print("   SYSTEMATIC VERIFICATION SUITE — SU(2) EMERGENCE OPTIMIZER")
    print("=" * 70)
    print()

    for module_name, func_name, desc in TESTS:
        tag = module_name.replace("test_", "")
        if selected and not any(tag.startswith(s) for s in selected):
            continue

        print(f"\n{'#' * 70}")
        print(f"# Running: {desc}")
        print(f"# Module:  {module_name}")
        print(f"{'#' * 70}\n")

        t0 = time.time()
        try:
            mod = importlib.import_module(module_name)
            fn = getattr(mod, func_name)
            fn()
            elapsed = time.time() - t0
            results.append((desc, "OK", elapsed))
        except Exception as e:
            elapsed = time.time() - t0
            results.append((desc, f"ERROR: {e}", elapsed))
            traceback.print_exc()

        print(f"\n[elapsed: {elapsed:.1f}s]")

    # ------ Final summary ------
    total_elapsed = time.time() - total_t0
    print("\n\n" + "=" * 70)
    print("   FINAL SUMMARY")
    print("=" * 70)
    print(f"\n{'Test':50s} {'Status':>10s}  {'Time':>7s}")
    print("-" * 70)
    for desc, status, t in results:
        short_status = "OK" if status == "OK" else "ERROR"
        print(f"{desc:50s} {short_status:>10s}  {t:6.1f}s")
    print("-" * 70)
    print(f"{'TOTAL':50s} {'':>10s}  {total_elapsed:6.1f}s")
    print(f"\nPassed: {sum(1 for _,s,_ in results if s == 'OK')}/{len(results)}")
    errors = [(d, s) for d, s, _ in results if s != "OK"]
    if errors:
        print("\nFailed tests:")
        for d, s in errors:
            print(f"  - {d}: {s}")
    print()


if __name__ == "__main__":
    sel = sys.argv[1:] if len(sys.argv) > 1 else None
    run_all(sel)
