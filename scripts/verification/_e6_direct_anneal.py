"""
E₆ blind discovery using direct+anneal (same strategy as F₄).

E₆: 78 generators in 27-dim (IRREDUCIBLE fundamental rep).
Density: 78/351 = 22.2% (higher than F₄'s 14.8%).

Strategy (proven for F₄):
  Stage 1: gamma=50, lr=0.001, 30k steps → K∝I + partial closure
  Stage 2: gamma=200, lr=0.0001, 30k steps → push to exact closure
"""
import sys, os, time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from _f4_direct_anneal import run_direct_anneal, analyze_detailed


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, nargs="+", default=[42, 137, 256])
    parser.add_argument("--s1-gamma", type=float, default=50.0)
    parser.add_argument("--s2-gamma", type=float, default=200.0)
    parser.add_argument("--s1-steps", type=int, default=30000)
    parser.add_argument("--s2-steps", type=int, default=30000)
    parser.add_argument("--s2-lr", type=float, default=0.0001)
    args = parser.parse_args()

    n_gen = 78  # E₆ dimension
    dim = 27    # fundamental rep
    batch_size = len(args.seeds)

    print(f"E6 direct+anneal: {n_gen} gens in so({dim}), {batch_size} seeds")
    print(f"  Density: {n_gen}/{dim*(dim-1)//2} = {n_gen/(dim*(dim-1)//2):.1%}")
    print(f"  Stage 1: gamma={args.s1_gamma}, lr=0.001, {args.s1_steps} steps")
    print(f"  Stage 2: gamma={args.s2_gamma}, lr={args.s2_lr}, {args.s2_steps} steps")
    print(f"  Seeds: {args.seeds}")

    t0 = time.time()
    all_gens, history = run_direct_anneal(
        n_gen=n_gen, dim=dim,
        stage1_steps=args.s1_steps,
        stage1_lr=0.001,
        stage1_gamma=args.s1_gamma,
        stage2_steps=args.s2_steps,
        stage2_lr=args.s2_lr,
        stage2_gamma=args.s2_gamma,
        batch_size=batch_size,
        seeds=args.seeds,
    )
    elapsed = time.time() - t0
    print(f"\nTotal: {elapsed:.1f}s ({elapsed/60:.1f} min)")

    for b in range(batch_size):
        gens_b = [all_gens[b, k] for k in range(n_gen)]
        analyze_detailed(gens_b, name=f"E6 Seed {args.seeds[b]}")
