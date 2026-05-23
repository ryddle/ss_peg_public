"""
Experiment 5 — Minimum Cycle Dose for Killing Classification
=============================================================
What is the minimum number of cycle-loss epochs needed for K_rr to
separate compositional from atomic relations (AUC > 0.90)?

Additionally: does the Killing separation survive "withdrawal" — continued
training WITHOUT cycle loss — or does T_K drive it back toward abelian
collapse (G_RAG-19)?

Design:
  Same mixed KG as Exp4 (16 Type A compositional + 8 Type B atomic).
  Phase 1: 800 ep triple-only with T_K + T_ortho (identical to Exp4 Phase 1).
  Phase 2: For each dose D in [1, 5, 10, 25, 50, 100, 200, 400]:
    (a) Restore Phase 1 checkpoint
    (b) Train D epochs with cycle loss → measure K_rr AUC + holonomy
    (c) "Withdrawal": 200 ep WITHOUT cycle loss → re-measure

  Key questions:
    Q1: Minimum dose for K_rr AUC > 0.90?
    Q2: Does Killing separation survive withdrawal (T_K re-collapse)?
    Q3: Is there a "permanent inoculation" threshold?

  Motivated by Exp4 findings:
    G_RAG-19: T_K drives abelian collapse without cycle loss
    G_RAG-21: Cycle loss unlocks Killing separation (AUC 0.49 → 0.98)
"""

import sys, os, time, copy
import torch
import numpy as np
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'training_ai_test'))
sys.stdout.reconfigure(encoding='utf-8')

from exp4_killing_classifier import (
    AlgKGModel, compute_cycle_loss, train_phase,
    analyze_killing_spectrum, evaluate_holonomy,
    generate_random_orthogonal_operators,
)
from kg_utils import (
    generate_relation_operators_with_closures,
    generate_cycles_v2,
    generate_entity_grounded_cycles,
    DEVICE,
)
from killing_utils import killing_metric


# ============================================================
#  CONFIG (same KG as Exp4 + Exp5-specific sweep parameters)
# ============================================================

SEED = 42
N_BASE_A = 8
N_CLOSURES_A = 8
N_RANDOM_B = 8
N_TYPE_A = N_BASE_A + N_CLOSURES_A   # 16
N_REL = N_TYPE_A + N_RANDOM_B        # 24
ANGLE_SCALE = 0.25
N_ENTITIES = 200
N_TRAIN, N_VAL = 8000, 2000
PHASE1_EPOCHS = 800
LR = 1e-3
BATCH_SIZE = 256
HEAD_HIDDEN = 32
LAMBDA_K, LAMBDA_ORTHO = 0.5, 1.0
LAMBDA_CYCLE = 1.0
MARGIN = 1.0
D = 3
GRAD_CLIP = 1.0

# Exp5-specific
DOSES = [1, 5, 10, 25, 50, 100, 200, 400]
WITHDRAWAL_EPOCHS = 200


# ============================================================
#  MAIN EXPERIMENT
# ============================================================

def run_experiment(seed=SEED):
    print("=" * 70)
    print("  Exp5 — Minimum Cycle Dose for Killing Classification")
    print("=" * 70)
    print(f"  Device: {DEVICE}")
    print(f"  Seed: {seed}")
    print(f"  Mixed KG: {N_TYPE_A} Type A + {N_RANDOM_B} Type B = {N_REL} relations")
    print(f"  Phase 1: {PHASE1_EPOCHS} ep triple-only (T_K + T_ortho)")
    print(f"  Doses: {DOSES}")
    print(f"  Withdrawal: {WITHDRAWAL_EPOCHS} ep (no cycle loss, T_K + T_ortho)")

    t0 = time.time()

    # ========================================================
    #  PHASE 0: DATA GENERATION (identical to Exp4)
    # ========================================================

    print("\n--- Phase 0: Data generation ---")

    Q_A, alphas, closure_cycles = generate_relation_operators_with_closures(
        n_base=N_BASE_A, n_closures=N_CLOSURES_A, d=D,
        angle_scale=ANGLE_SCALE, seed=seed,
    )
    n_gen = Q_A.shape[1]

    I_np = np.eye(n_gen)
    closure_hols = [np.linalg.norm(Q_A[rc] @ Q_A[rj] @ Q_A[ri] - I_np, 'fro')
                    for ri, rj, rc in closure_cycles]
    print(f"  Type A: {Q_A.shape[0]} ops, closure hol max={max(closure_hols):.2e}")

    Q_B = generate_random_orthogonal_operators(
        N_RANDOM_B, n_gen, angle_scale=ANGLE_SCALE, seed=seed + 5000,
    )
    print(f"  Type B: {Q_B.shape[0]} ops")

    Q_all = np.concatenate([Q_A, Q_B], axis=0)
    type_labels = np.array([0] * N_TYPE_A + [1] * N_RANDOM_B)

    # GT Killing
    Q_all_t = torch.tensor(Q_all, dtype=torch.float64, device=DEVICE)
    with torch.no_grad():
        K_gt = killing_metric(Q_all_t).cpu().numpy()
    K_gt_diag = np.diag(K_gt)
    try:
        gt_k_auc = roc_auc_score(type_labels, K_gt_diag)
    except ValueError:
        gt_k_auc = float('nan')
    print(f"  GT K_rr AUC: {gt_k_auc:.4f}")

    # Synthetic KG
    rng = np.random.default_rng(seed)
    entities = rng.standard_normal((N_ENTITIES, n_gen))
    entities = entities / np.linalg.norm(entities, axis=1, keepdims=True)

    transformed = np.einsum('rij,hj->rhi', Q_all, entities)
    nearest_tail = np.zeros((N_REL, N_ENTITIES), dtype=np.int64)
    nearest_dist = np.zeros((N_REL, N_ENTITIES), dtype=np.float64)
    for r in range(N_REL):
        for h in range(N_ENTITIES):
            diffs = entities - transformed[r, h]
            dists = np.linalg.norm(diffs, axis=1)
            dists[h] = 1e10
            nearest_tail[r, h] = np.argmin(dists)
            nearest_dist[r, h] = dists[nearest_tail[r, h]]

    def _gen_split(n_triples, rng_local):
        triples, labels = [], []
        n_valid = n_triples // 2
        for _ in range(n_valid):
            r_idx = rng_local.integers(N_REL)
            h = rng_local.integers(N_ENTITIES)
            triples.append([h, r_idx, nearest_tail[r_idx, h]])
            labels.append(1)
        for _ in range(n_triples - n_valid):
            r_idx = rng_local.integers(N_REL)
            h = rng_local.integers(N_ENTITIES)
            t = rng_local.integers(N_ENTITIES)
            while t == nearest_tail[r_idx, h] or t == h:
                t = rng_local.integers(N_ENTITIES)
            triples.append([h, r_idx, t])
            labels.append(0)
        perm = rng_local.permutation(len(triples))
        return np.array(triples)[perm], np.array(labels)[perm]

    train_triples, train_labels = _gen_split(N_TRAIN, np.random.default_rng(seed + 100))
    val_triples, val_labels = _gen_split(N_VAL, np.random.default_rng(seed + 200))
    print(f"  Train: {len(train_triples)}, Val: {len(val_triples)}")

    # Cycles
    cycle_data = generate_cycles_v2(
        Q_all, closure_cycles,
        n_consistent=200, n_inconsistent=200, seed=seed,
    )
    con_cycles_t = torch.tensor(cycle_data['consistent_cycles'], dtype=torch.long, device=DEVICE)
    incon_cycles_t = torch.tensor(cycle_data['inconsistent_cycles'], dtype=torch.long, device=DEVICE)

    # Entity-grounded paths
    con_paths, incon_paths = generate_entity_grounded_cycles(
        nearest_tail, closure_cycles,
        N_ENTITIES, N_REL,
        n_consistent=300, n_inconsistent=300, seed=seed + 500,
    )
    gt_con_hol = np.zeros(len(con_paths))
    gt_incon_hol = np.zeros(len(incon_paths))
    for i, (h1, r1, h2, r2, h3, r3) in enumerate(con_paths):
        v = Q_all[r3] @ Q_all[r2] @ Q_all[r1] @ entities[h1]
        gt_con_hol[i] = np.linalg.norm(v - entities[h1])
    for i, (h1, r1, h2, r2, h3, r3) in enumerate(incon_paths):
        v = Q_all[r3] @ Q_all[r2] @ Q_all[r1] @ entities[h1]
        gt_incon_hol[i] = np.linalg.norm(v - entities[h1])

    print(f"  Cycles: {len(cycle_data['consistent_cycles'])}c + "
          f"{len(cycle_data['inconsistent_cycles'])}i, "
          f"Paths: {len(con_paths)}c + {len(incon_paths)}i")

    # ========================================================
    #  PHASE 1: TRIPLE-ONLY TRAINING (identical to Exp4)
    # ========================================================

    print(f"\n{'='*70}")
    print(f"  PHASE 1: Triple-only Training ({PHASE1_EPOCHS} ep)")
    print(f"{'='*70}")

    torch.manual_seed(seed)
    model = AlgKGModel(N_ENTITIES, N_REL, n_gen, HEAD_HIDDEN).to(DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)

    hist_p1 = train_phase(
        model, optimizer,
        train_triples, train_labels, val_triples, val_labels,
        epochs=PHASE1_EPOCHS, batch_size=BATCH_SIZE,
        lambda_k=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
        lambda_cycle=0.0,
        con_cycles_t=con_cycles_t, incon_cycles_t=incon_cycles_t,
        grad_clip=GRAD_CLIP, label="Phase1", log_every=400,
    )

    # Phase 1 baseline metrics
    analysis_p1 = analyze_killing_spectrum(model, type_labels, N_TYPE_A, N_RANDOM_B)
    hol_p1 = evaluate_holonomy(model, cycle_data, con_paths, incon_paths,
                                gt_con_hol, gt_incon_hol, label="Phase1")
    p1_hol_sep = hist_p1['hol_incon'][-1] / max(hist_p1['hol_con'][-1], 1e-10)
    print(f"  Phase 1 K_rr AUC: {analysis_p1['k_diag_auc']:.4f}, "
          f"hol_sep: {p1_hol_sep:.2f}x")

    # Save checkpoint
    state_p1 = copy.deepcopy(model.state_dict())
    adam_p1 = copy.deepcopy(optimizer.state_dict())

    t_phase1 = time.time() - t0

    # ========================================================
    #  PHASE 2: DOSE SWEEP
    # ========================================================

    print(f"\n{'='*70}")
    print(f"  PHASE 2: Cycle Loss Dose Sweep")
    print(f"{'='*70}")

    results = []

    for dose in DOSES:
        t_dose_start = time.time()

        print(f"\n{'─'*60}")
        print(f"  DOSE: {dose} ep cycle loss + {WITHDRAWAL_EPOCHS} ep withdrawal")
        print(f"{'─'*60}")

        # Restore Phase 1 checkpoint
        model.load_state_dict(copy.deepcopy(state_p1))
        optimizer = torch.optim.Adam(model.parameters(), lr=LR)
        optimizer.load_state_dict(copy.deepcopy(adam_p1))

        # --- Dose training (cycle loss ON) ---
        log_every_dose = dose if dose <= 50 else 100
        hist_dose = train_phase(
            model, optimizer,
            train_triples, train_labels, val_triples, val_labels,
            epochs=dose, batch_size=BATCH_SIZE,
            lambda_k=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
            lambda_cycle=LAMBDA_CYCLE,
            con_cycles_t=con_cycles_t, incon_cycles_t=incon_cycles_t,
            margin=MARGIN, grad_clip=GRAD_CLIP,
            label=f"D{dose:03d}", log_every=log_every_dose,
        )

        # Post-dose analysis
        analysis_dose = analyze_killing_spectrum(model, type_labels, N_TYPE_A, N_RANDOM_B)
        hol_dose = evaluate_holonomy(model, cycle_data, con_paths, incon_paths,
                                      gt_con_hol, gt_incon_hol, label=f"D{dose:03d}")
        hol_sep_dose = hist_dose['hol_incon'][-1] / max(hist_dose['hol_con'][-1], 1e-10)

        print(f"  >> Post-dose: K_rr AUC={analysis_dose['k_diag_auc']:.4f}  "
              f"Comm A|A={analysis_dose['comm_AA']:.6f} B|B={analysis_dose['comm_BB']:.6f}")

        # --- Withdrawal (cycle loss OFF, T_K + T_ortho still on) ---
        print(f"\n  --- Withdrawal: {WITHDRAWAL_EPOCHS} ep (no cycle loss) ---")
        hist_withdraw = train_phase(
            model, optimizer,
            train_triples, train_labels, val_triples, val_labels,
            epochs=WITHDRAWAL_EPOCHS, batch_size=BATCH_SIZE,
            lambda_k=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
            lambda_cycle=0.0,
            con_cycles_t=con_cycles_t, incon_cycles_t=incon_cycles_t,
            grad_clip=GRAD_CLIP,
            label=f"W{dose:03d}", log_every=50,
        )

        # Post-withdrawal analysis
        analysis_w = analyze_killing_spectrum(model, type_labels, N_TYPE_A, N_RANDOM_B)
        hol_w = evaluate_holonomy(model, cycle_data, con_paths, incon_paths,
                                   gt_con_hol, gt_incon_hol, label=f"W{dose:03d}")
        hol_sep_w = hist_withdraw['hol_incon'][-1] / max(hist_withdraw['hol_con'][-1], 1e-10)

        print(f"  >> Post-withdrawal: K_rr AUC={analysis_w['k_diag_auc']:.4f}  "
              f"Comm A|A={analysis_w['comm_AA']:.6f} B|B={analysis_w['comm_BB']:.6f}")

        t_dose_total = time.time() - t_dose_start
        print(f"  Dose {dose} total time: {t_dose_total:.1f}s")

        # Withdrawal hol_sep trajectory (at each epoch)
        w_hol_sep = [
            hist_withdraw['hol_incon'][i] / max(hist_withdraw['hol_con'][i], 1e-10)
            for i in range(len(hist_withdraw['hol_con']))
        ]

        results.append({
            'dose': dose,
            # Post-dose
            'k_auc_d': analysis_dose['k_diag_auc'],
            'k_abs_auc_d': analysis_dose['k_diag_abs_auc'],
            'comm_auc_d': analysis_dose['comm_total_auc'],
            'hol_raw_d': hol_dose['raw_auc'],
            'hol_path_d': hol_dose['path_auc'],
            'hol_sep_d': hol_sep_dose,
            'hol_corr_d': hol_dose['correlation'],
            'comm_AA_d': analysis_dose['comm_AA'],
            'comm_AB_d': analysis_dose['comm_AB'],
            'comm_BB_d': analysis_dose['comm_BB'],
            'val_acc_d': hist_dose['val_acc'][-1],
            'tk_d': hist_dose['tk'][-1],
            # Post-withdrawal
            'k_auc_w': analysis_w['k_diag_auc'],
            'k_abs_auc_w': analysis_w['k_diag_abs_auc'],
            'hol_raw_w': hol_w['raw_auc'],
            'hol_path_w': hol_w['path_auc'],
            'hol_sep_w': hol_sep_w,
            'hol_corr_w': hol_w['correlation'],
            'comm_AA_w': analysis_w['comm_AA'],
            'comm_AB_w': analysis_w['comm_AB'],
            'comm_BB_w': analysis_w['comm_BB'],
            'val_acc_w': hist_withdraw['val_acc'][-1],
            'tk_w': hist_withdraw['tk'][-1],
            # Withdrawal trajectory
            'w_hol_sep_traj': w_hol_sep,
            'time_s': t_dose_total,
        })

    # ========================================================
    #  SUMMARY
    # ========================================================

    t_total = time.time() - t0

    print(f"\n{'='*70}")
    print(f"  SUMMARY — Exp5: Minimum Cycle Dose for Killing Classification")
    print(f"{'='*70}")

    print(f"\n  Baseline (Phase 1, dose=0):")
    print(f"    K_rr AUC: {analysis_p1['k_diag_auc']:.4f}  "
          f"path_AUC: {hol_p1['path_auc']:.4f}  hol_sep: {p1_hol_sep:.2f}x")

    # Main dose-response table
    print(f"\n  Dose-Response Table:")
    print(f"  {'Dose':>6s}  {'K(D)':>7s} {'pAUC(D)':>8s} {'sep(D)':>7s}  "
          f"{'K(W)':>7s} {'pAUC(W)':>8s} {'sep(W)':>7s}  {'K_ret':>6s} {'H_ret':>6s}")
    print(f"  {'------':>6s}  {'-------':>7s} {'--------':>8s} {'-------':>7s}  "
          f"{'-------':>7s} {'--------':>8s} {'-------':>7s}  {'------':>6s} {'------':>6s}")

    for r in results:
        # Killing retention: ratio of post-withdrawal to post-dose K_rr AUC
        # Only meaningful if post-dose K_rr AUC > 0.5 (better than random)
        if r['k_auc_d'] > 0.55:
            k_ret = r['k_auc_w'] / r['k_auc_d']
            k_ret_s = f"{k_ret:.3f}"
        else:
            k_ret_s = "  n/a "

        # Holonomy separation retention
        h_ret = r['hol_sep_w'] / max(r['hol_sep_d'], 0.01)
        print(f"  {r['dose']:6d}  {r['k_auc_d']:7.4f} {r['hol_path_d']:8.4f} "
              f"{r['hol_sep_d']:7.2f}  "
              f"{r['k_auc_w']:7.4f} {r['hol_path_w']:8.4f} {r['hol_sep_w']:7.2f}  "
              f"{k_ret_s:>6s} {h_ret:6.3f}")

    # Find thresholds
    threshold = 0.90
    min_dose_k = None
    min_dose_k_persist = None
    min_dose_hol = None
    min_dose_hol_persist = None
    for r in results:
        if r['k_auc_d'] >= threshold and min_dose_k is None:
            min_dose_k = r['dose']
        if r['k_auc_w'] >= threshold and min_dose_k_persist is None:
            min_dose_k_persist = r['dose']
        if r['hol_path_d'] >= threshold and min_dose_hol is None:
            min_dose_hol = r['dose']
        if r['hol_path_w'] >= threshold and min_dose_hol_persist is None:
            min_dose_hol_persist = r['dose']

    print(f"\n  Minimum dose for AUC > {threshold}:")
    print(f"    K_rr  post-dose:       {min_dose_k if min_dose_k else '> 400'} ep")
    print(f"    K_rr  post-withdrawal:  {min_dose_k_persist if min_dose_k_persist else '> 400'} ep")
    print(f"    Hol   post-dose:       {min_dose_hol if min_dose_hol else '> 400'} ep")
    print(f"    Hol   post-withdrawal:  {min_dose_hol_persist if min_dose_hol_persist else '> 400'} ep")

    # Commutator block evolution
    print(f"\n  Commutator block ratios (A|A / B|B):")
    print(f"  {'Dose':>6s}  {'ratio(D)':>9s}  {'ratio(W)':>9s}  {'change':>7s}")
    print(f"  {'------':>6s}  {'---------':>9s}  {'---------':>9s}  {'-------':>7s}")
    for r in results:
        ratio_d = r['comm_AA_d'] / max(r['comm_BB_d'], 1e-10)
        ratio_w = r['comm_AA_w'] / max(r['comm_BB_w'], 1e-10)
        change = ratio_w / max(ratio_d, 1e-10)
        print(f"  {r['dose']:6d}  {ratio_d:9.2f}  {ratio_w:9.2f}  {change:7.3f}")

    print(f"\n  Wall time: {t_total:.1f}s (Phase 1: {t_phase1:.1f}s)")

    # --- Save results ---
    save_path = os.path.join(os.path.dirname(__file__), 'exp5_results.npz')
    save_dict = {
        'doses': np.array(DOSES),
        'type_labels': type_labels,
        'p1_k_auc': analysis_p1['k_diag_auc'],
        'p1_path_auc': hol_p1['path_auc'],
        'p1_hol_sep': p1_hol_sep,
        'gt_k_auc': gt_k_auc,
    }
    for r in results:
        d = r['dose']
        for key, val in r.items():
            if key == 'w_hol_sep_traj':
                save_dict[f'd{d}_w_hol_sep_traj'] = np.array(val)
            elif isinstance(val, (int, float)):
                save_dict[f'd{d}_{key}'] = val
    np.savez(save_path, **save_dict)
    print(f"  Results saved to {save_path}")


if __name__ == "__main__":
    run_experiment()
