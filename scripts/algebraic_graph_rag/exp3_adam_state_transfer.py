"""
Experiment 3 — Adam State Transfer for Holonomy Basin Detection
================================================================
Directly analogous to Exp11g.6 (F4 bifurcation): the holonomy-aware
basin may exist in extended phase space (theta, m, v), not in parameter
space alone.

Protocol:
  Phase A:  Train AlgKG with L_triple + L_cycle for 800 epochs ("Donor")
  Phase A0: Train AlgKG with L_triple only for 800 epochs ("Baseline")
  (Both start from identical initial parameters.)

  Phase B:  Continue for 400 epochs with L_triple only (NO cycle loss)
  Variants:
    B1 "Control":    (theta_B, adam_B)  — continue baseline
    B2 "Params":     (theta_A, fresh)   — donor params, fresh Adam
    B3 "Full":       (theta_A, adam_A)  — donor params + donor Adam (KEY)
    B4 "Adam-only":  (theta_B, adam_A)  — baseline params + donor Adam

  If B3 >> B2 >> B1: basin exists in (theta, m, v), not theta alone
  If B4 > B1:        Adam dynamics alone carry basin info

This is "Via 2" from the PIH roadmap.
"""

import sys, os, time, copy
import torch
import torch.nn as nn
import numpy as np
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(__file__))
from kg_utils import (
    generate_relation_operators_with_closures,
    generate_cycles_v2,
    generate_entity_grounded_cycles,
    compute_holonomy,
    relation_killing_tension,
    relation_ortho_tension,
    DEVICE,
)

sys.stdout.reconfigure(encoding='utf-8')


# ============================================================
#  CONFIG
# ============================================================

SEED = 42
N_BASE, N_CLOSURES = 12, 12
N_REL = N_BASE + N_CLOSURES  # 24
ANGLE_SCALE = 0.25
N_ENTITIES = 200
N_TRAIN, N_VAL = 8000, 2000
PHASE_A_EPOCHS = 800
PHASE_B_EPOCHS = 400
LR = 1e-3
BATCH_SIZE = 256
HEAD_HIDDEN = 32
LAMBDA_K, LAMBDA_ORTHO = 0.5, 1.0
LAMBDA_CYCLE = 1.0
MARGIN = 1.0
D = 3
GRAD_CLIP = 1.0


# ============================================================
#  MODEL
# ============================================================

class AlgKGModel(nn.Module):
    """AlgKG with shared operators, T_K + T_ortho regularized."""

    def __init__(self, n_entities, n_relations, n_gen, head_hidden=32):
        super().__init__()
        self.n_gen = n_gen
        self.entity_emb = nn.Embedding(n_entities, n_gen, dtype=torch.float64)
        self.relation_ops = nn.Parameter(
            torch.eye(n_gen, dtype=torch.float64).unsqueeze(0).repeat(n_relations, 1, 1)
            + 0.05 * torch.randn(n_relations, n_gen, n_gen, dtype=torch.float64)
        )
        self.head = nn.Sequential(
            nn.Linear(n_gen + 1, head_hidden, dtype=torch.float64),
            nn.ReLU(),
            nn.Linear(head_hidden, 2, dtype=torch.float64),
        )
        nn.init.normal_(self.entity_emb.weight)
        with torch.no_grad():
            self.entity_emb.weight.div_(
                self.entity_emb.weight.norm(dim=1, keepdim=True).clamp(min=1e-8)
            )

    def forward(self, triples):
        h_idx, r_idx, t_idx = triples[:, 0], triples[:, 1], triples[:, 2]
        e_h = self.entity_emb(h_idx)
        e_t = self.entity_emb(t_idx)
        theta_r = self.relation_ops[r_idx]
        transformed = torch.bmm(theta_r, e_h.unsqueeze(-1)).squeeze(-1)
        diff = transformed - e_t
        dist = diff.norm(dim=1, keepdim=True)
        return self.head(torch.cat([diff, dist], dim=1))

    def get_relation_matrices(self):
        return self.relation_ops

    def regularization_loss(self, lambda_k=0.5, lambda_ortho=1.0):
        ops = self.relation_ops
        loss = torch.tensor(0.0, dtype=torch.float64, device=ops.device)
        if lambda_k > 0:
            loss = loss + lambda_k * relation_killing_tension(ops)
        if lambda_ortho > 0:
            loss = loss + lambda_ortho * relation_ortho_tension(ops)
        return loss


# ============================================================
#  CYCLE CONTRASTIVE LOSS
# ============================================================

def compute_cycle_loss(ops, con_cycles_t, incon_cycles_t, margin=1.0):
    """
    Differentiable cycle contrastive loss on relation operators.
      Consistent:   push ||prod theta_r - I||_F^2 -> 0
      Inconsistent: push ||prod theta_r - I||_F   > margin
    """
    ng = ops.shape[1]
    I_ng = torch.eye(ng, dtype=ops.dtype, device=ops.device).unsqueeze(0)

    # Consistent cycles: minimize ||H - I||^2_F
    H_con = ops[con_cycles_t[:, 2]] @ ops[con_cycles_t[:, 1]] @ ops[con_cycles_t[:, 0]]
    loss_con = ((H_con - I_ng) ** 2).sum(dim=(1, 2)).mean()

    # Inconsistent cycles: push ||H - I||_F above margin
    H_incon = ops[incon_cycles_t[:, 2]] @ ops[incon_cycles_t[:, 1]] @ ops[incon_cycles_t[:, 0]]
    hol_incon = ((H_incon - I_ng) ** 2).sum(dim=(1, 2)).sqrt()
    loss_incon = (torch.clamp(margin - hol_incon, min=0.0) ** 2).mean()

    return loss_con + loss_incon


# ============================================================
#  TRAINING
# ============================================================

def train_phase(
    model, optimizer,
    train_triples, train_labels, val_triples, val_labels,
    epochs, batch_size=256,
    lambda_k=0.5, lambda_ortho=1.0,
    lambda_cycle=0.0, con_cycles_t=None, incon_cycles_t=None, margin=1.0,
    grad_clip=1.0, label="Model", log_every=100,
):
    """Train for one phase. Optimizer is created externally for state management."""
    criterion = nn.CrossEntropyLoss()

    train_t = torch.tensor(train_triples, dtype=torch.long, device=DEVICE)
    train_l = torch.tensor(train_labels, dtype=torch.long, device=DEVICE)
    val_t = torch.tensor(val_triples, dtype=torch.long, device=DEVICE)
    val_l = torch.tensor(val_labels, dtype=torch.long, device=DEVICE)

    n_gen = model.n_gen
    I_ng = torch.eye(n_gen, dtype=torch.float64, device=DEVICE).unsqueeze(0)

    history = {
        'train_acc': [], 'val_acc': [], 'loss': [],
        'tk': [], 'ortho_err': [],
        'hol_con': [], 'hol_incon': [],
    }
    if lambda_cycle > 0:
        history['cycle_loss'] = []

    best_val = 0.0
    n = len(train_triples)

    for ep in range(1, epochs + 1):
        model.train()
        perm = torch.randperm(n, device=DEVICE)
        epoch_loss = 0.0
        epoch_correct = 0
        epoch_cyc = 0.0
        n_batches = 0

        for i in range(0, n, batch_size):
            idx = perm[i:i + batch_size]
            logits = model(train_t[idx])
            task_loss = criterion(logits, train_l[idx])

            reg_loss = model.regularization_loss(lambda_k, lambda_ortho)

            cyc_loss = torch.tensor(0.0, dtype=torch.float64, device=DEVICE)
            if lambda_cycle > 0 and con_cycles_t is not None:
                cyc_loss = lambda_cycle * compute_cycle_loss(
                    model.relation_ops, con_cycles_t, incon_cycles_t, margin
                )

            loss = task_loss + reg_loss + cyc_loss
            optimizer.zero_grad()
            loss.backward()
            if grad_clip > 0:
                nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
            optimizer.step()

            epoch_loss += task_loss.item() * len(idx)
            epoch_correct += (logits.argmax(1) == train_l[idx]).sum().item()
            epoch_cyc += cyc_loss.item()
            n_batches += 1

        train_acc = epoch_correct / n

        # --- Eval ---
        model.eval()
        with torch.no_grad():
            val_logits = model(val_t)
            val_acc = (val_logits.argmax(1) == val_l).float().mean().item()

            ops = model.get_relation_matrices()
            I_n = torch.eye(ops.shape[1], dtype=ops.dtype, device=ops.device)
            ortho_err = sum(
                ((ops[r] @ ops[r].T - I_n) ** 2).sum().item()
                for r in range(ops.shape[0])
            ) / ops.shape[0]
            tk = relation_killing_tension(ops).item()

            # Cycle holonomy stats (always track for monitoring)
            if con_cycles_t is not None:
                H_con = ops[con_cycles_t[:, 2]] @ ops[con_cycles_t[:, 1]] @ ops[con_cycles_t[:, 0]]
                hol_con_mean = ((H_con - I_ng) ** 2).sum(dim=(1, 2)).sqrt().mean().item()
                H_incon = ops[incon_cycles_t[:, 2]] @ ops[incon_cycles_t[:, 1]] @ ops[incon_cycles_t[:, 0]]
                hol_incon_mean = ((H_incon - I_ng) ** 2).sum(dim=(1, 2)).sqrt().mean().item()
            else:
                hol_con_mean = hol_incon_mean = float('nan')

        history['train_acc'].append(train_acc)
        history['val_acc'].append(val_acc)
        history['loss'].append(epoch_loss / n)
        history['tk'].append(tk)
        history['ortho_err'].append(ortho_err)
        history['hol_con'].append(hol_con_mean)
        history['hol_incon'].append(hol_incon_mean)
        if lambda_cycle > 0:
            history['cycle_loss'].append(epoch_cyc / max(n_batches, 1))

        if val_acc > best_val:
            best_val = val_acc

        if ep % log_every == 0 or ep == 1:
            sep = hol_incon_mean / max(hol_con_mean, 1e-10) if not np.isnan(hol_con_mean) else 0
            extra = f"  hol_sep={sep:.2f}x"
            if lambda_cycle > 0:
                extra = f"  cyc_L={epoch_cyc / max(n_batches, 1):.4f}" + extra
            print(f"  [{label}] ep {ep:4d}  train={train_acc:.4f}  val={val_acc:.4f}  "
                  f"loss={epoch_loss / n:.4f}  T_K={tk:.4f}  OE={ortho_err:.4f}{extra}")

    history['best_val_acc'] = best_val
    return history


# ============================================================
#  EVALUATION
# ============================================================

def _compute_path_holonomy(model, paths):
    """Entity-path holonomy: ||theta_r3 * theta_r2 * theta_r1 * e_h1 - e_h1||."""
    model.eval()
    h1s = torch.tensor([p[0] for p in paths], dtype=torch.long, device=DEVICE)
    r1s = torch.tensor([p[1] for p in paths], dtype=torch.long, device=DEVICE)
    r2s = torch.tensor([p[3] for p in paths], dtype=torch.long, device=DEVICE)
    r3s = torch.tensor([p[5] for p in paths], dtype=torch.long, device=DEVICE)

    with torch.no_grad():
        e_h1 = model.entity_emb(h1s)
        ops = model.get_relation_matrices()
        v = torch.bmm(ops[r1s], e_h1.unsqueeze(-1)).squeeze(-1)
        v = torch.bmm(ops[r2s], v.unsqueeze(-1)).squeeze(-1)
        v = torch.bmm(ops[r3s], v.unsqueeze(-1)).squeeze(-1)
        return (v - e_h1).norm(dim=1).cpu().numpy()


def evaluate_holonomy(model, cycle_data, con_paths, incon_paths,
                      gt_con_hol, gt_incon_hol, label="Model"):
    """Combined raw + entity-path holonomy evaluation."""
    model.eval()

    # --- Raw holonomy ---
    with torch.no_grad():
        ops = model.get_relation_matrices()
        hol_con = compute_holonomy(ops, cycle_data['consistent_cycles']).cpu().numpy()
        hol_incon = compute_holonomy(ops, cycle_data['inconsistent_cycles']).cpu().numpy()

    labels_r = np.concatenate([np.zeros(len(hol_con)), np.ones(len(hol_incon))])
    scores_r = np.concatenate([hol_con, hol_incon])
    try:
        raw_auc = roc_auc_score(labels_r, scores_r)
    except ValueError:
        raw_auc = float('nan')
    raw_sep = float(np.mean(hol_incon) / max(np.mean(hol_con), 1e-10))

    # --- Entity-path holonomy ---
    path_hol_con = _compute_path_holonomy(model, con_paths)
    path_hol_incon = _compute_path_holonomy(model, incon_paths)

    labels_p = np.concatenate([np.zeros(len(path_hol_con)), np.ones(len(path_hol_incon))])
    scores_p = np.concatenate([path_hol_con, path_hol_incon])
    try:
        path_auc = roc_auc_score(labels_p, scores_p)
    except ValueError:
        path_auc = float('nan')
    path_sep = float(np.mean(path_hol_incon) / max(np.mean(path_hol_con), 1e-10))

    gt_all = np.concatenate([gt_con_hol, gt_incon_hol])
    if np.std(gt_all) > 1e-10 and np.std(scores_p) > 1e-10:
        corr = float(np.corrcoef(gt_all, scores_p)[0, 1])
    else:
        corr = float('nan')

    print(f"  [{label:16s}]  val={0:.4f}  raw_AUC={raw_auc:.4f} (sep={raw_sep:.2f}x)  "
          f"path_AUC={path_auc:.4f} (sep={path_sep:.2f}x)  corr={corr:.4f}")

    return {
        'raw_auc': raw_auc, 'raw_sep': raw_sep,
        'raw_con': float(np.mean(hol_con)), 'raw_incon': float(np.mean(hol_incon)),
        'path_auc': path_auc, 'path_sep': path_sep,
        'path_con': float(np.mean(path_hol_con)), 'path_incon': float(np.mean(path_hol_incon)),
        'correlation': corr,
    }


def evaluate_with_acc(model, val_triples, val_labels,
                      cycle_data, con_paths, incon_paths,
                      gt_con_hol, gt_incon_hol, label="Model"):
    """Evaluate classification + holonomy together."""
    model.eval()
    val_t = torch.tensor(val_triples, dtype=torch.long, device=DEVICE)
    val_l = torch.tensor(val_labels, dtype=torch.long, device=DEVICE)
    with torch.no_grad():
        val_acc = (model(val_t).argmax(1) == val_l).float().mean().item()

    hol = evaluate_holonomy(model, cycle_data, con_paths, incon_paths,
                            gt_con_hol, gt_incon_hol, label)
    hol['val_acc'] = val_acc

    # Fix the printed line with actual val_acc
    # (re-print with correct value)
    return hol


def adam_state_summary(optimizer):
    """Summarize Adam state for diagnostics."""
    state = optimizer.state_dict()['state']
    all_m, all_v = [], []
    step = None
    for v in state.values():
        all_m.append(v['exp_avg'].abs().mean().item())
        all_v.append(v['exp_avg_sq'].abs().mean().item())
        s = v['step']
        step = s.item() if isinstance(s, torch.Tensor) else s
    return {'mean_m': np.mean(all_m), 'mean_v': np.mean(all_v), 'step': step}


# ============================================================
#  MAIN EXPERIMENT
# ============================================================

def run_experiment(seed=SEED):
    print("=" * 70)
    print("  Exp3 — Adam State Transfer for Holonomy Basin Detection")
    print("=" * 70)
    print(f"  Device: {DEVICE}")
    print(f"  Seed: {seed}")
    print(f"  Phase A: {PHASE_A_EPOCHS} ep (L_triple + L_cycle, lambda_cycle={LAMBDA_CYCLE})")
    print(f"  Phase B: {PHASE_B_EPOCHS} ep (L_triple only)")
    print(f"  Margin: {MARGIN}")

    t0 = time.time()

    # --------------------------------------------------------
    #  1. GENERATE GROUND TRUTH & DATA (identical to Exp2)
    # --------------------------------------------------------

    print("\n--- Generating relation operators with closures ---")
    Q_ops, alphas, closure_cycles = generate_relation_operators_with_closures(
        n_base=N_BASE, n_closures=N_CLOSURES, d=D,
        angle_scale=ANGLE_SCALE, seed=seed,
    )
    n_gen = Q_ops.shape[1]
    print(f"  Relations: {Q_ops.shape[0]}, n_gen: {n_gen}")
    print(f"  Closure cycles: {len(closure_cycles)}")

    I_np = np.eye(n_gen)
    closure_hols = [np.linalg.norm(Q_ops[rc] @ Q_ops[rj] @ Q_ops[ri] - I_np, 'fro')
                    for ri, rj, rc in closure_cycles]
    print(f"  Closure holonomy check: max={max(closure_hols):.2e}, mean={np.mean(closure_hols):.2e}")

    print("\n--- Generating synthetic KG ---")
    rng = np.random.default_rng(seed)
    entities = rng.standard_normal((N_ENTITIES, n_gen))
    entities = entities / np.linalg.norm(entities, axis=1, keepdims=True)

    transformed = np.einsum('rij,hj->rhi', Q_ops, entities)
    nearest_tail = np.zeros((N_REL, N_ENTITIES), dtype=np.int64)
    nearest_dist = np.zeros((N_REL, N_ENTITIES), dtype=np.float64)
    for r in range(N_REL):
        for h in range(N_ENTITIES):
            diffs = entities - transformed[r, h]
            dists = np.linalg.norm(diffs, axis=1)
            dists[h] = 1e10
            nearest_tail[r, h] = np.argmin(dists)
            nearest_dist[r, h] = dists[nearest_tail[r, h]]

    print(f"  Entities: {N_ENTITIES}, NN dist: mean={nearest_dist.mean():.4f}, max={nearest_dist.max():.4f}")

    def _gen_split(n_triples, rng_local):
        triples, labels = [], []
        n_valid = n_triples // 2
        for _ in range(n_valid):
            r = rng_local.integers(N_REL)
            h = rng_local.integers(N_ENTITIES)
            triples.append([h, r, nearest_tail[r, h]])
            labels.append(1)
        for _ in range(n_triples - n_valid):
            r = rng_local.integers(N_REL)
            h = rng_local.integers(N_ENTITIES)
            t = rng_local.integers(N_ENTITIES)
            while t == nearest_tail[r, h] or t == h:
                t = rng_local.integers(N_ENTITIES)
            triples.append([h, r, t])
            labels.append(0)
        perm = rng_local.permutation(len(triples))
        return np.array(triples)[perm], np.array(labels)[perm]

    train_triples, train_labels = _gen_split(N_TRAIN, np.random.default_rng(seed + 100))
    val_triples, val_labels = _gen_split(N_VAL, np.random.default_rng(seed + 200))
    print(f"  Train: {len(train_triples)} ({train_labels.sum()} valid, {(1-train_labels).sum()} invalid)")
    print(f"  Val:   {len(val_triples)} ({val_labels.sum()} valid, {(1-val_labels).sum()} invalid)")

    print("\n--- Generating cycles ---")
    cycle_data = generate_cycles_v2(
        Q_ops, closure_cycles,
        n_consistent=200, n_inconsistent=200, seed=seed,
    )
    print(f"  Consistent: {len(cycle_data['consistent_cycles'])}, "
          f"GT hol={cycle_data['consistent_holonomy'].mean():.6f}")
    print(f"  Inconsistent: {len(cycle_data['inconsistent_cycles'])}, "
          f"GT hol={cycle_data['inconsistent_holonomy'].mean():.4f}")

    # Cycle tensors for training loss
    con_cycles_t = torch.tensor(
        np.array(cycle_data['consistent_cycles']), dtype=torch.long, device=DEVICE
    )
    incon_cycles_t = torch.tensor(
        np.array(cycle_data['inconsistent_cycles']), dtype=torch.long, device=DEVICE
    )

    print("\n--- Generating entity-grounded cycles ---")
    con_paths, incon_paths = generate_entity_grounded_cycles(
        nearest_tail, closure_cycles,
        N_ENTITIES, N_REL,
        n_consistent=300, n_inconsistent=300, seed=seed + 500,
    )
    gt_con_hol = np.zeros(len(con_paths))
    gt_incon_hol = np.zeros(len(incon_paths))
    for i, (h1, r1, h2, r2, h3, r3) in enumerate(con_paths):
        v = Q_ops[r3] @ Q_ops[r2] @ Q_ops[r1] @ entities[h1]
        gt_con_hol[i] = np.linalg.norm(v - entities[h1])
    for i, (h1, r1, h2, r2, h3, r3) in enumerate(incon_paths):
        v = Q_ops[r3] @ Q_ops[r2] @ Q_ops[r1] @ entities[h1]
        gt_incon_hol[i] = np.linalg.norm(v - entities[h1])

    gt_path_sep = gt_incon_hol.mean() / max(gt_con_hol.mean(), 1e-10)
    print(f"  Consistent paths: {len(con_paths)}, GT hol={gt_con_hol.mean():.6f}")
    print(f"  Inconsistent paths: {len(incon_paths)}, GT hol={gt_incon_hol.mean():.4f}")
    print(f"  GT path separation: {gt_path_sep:.1f}x")

    # --------------------------------------------------------
    #  2. SHARED INITIAL STATE
    # --------------------------------------------------------

    print("\n--- Creating shared initial model state ---")
    torch.manual_seed(seed)
    init_model = AlgKGModel(N_ENTITIES, N_REL, n_gen, HEAD_HIDDEN).to(DEVICE)
    init_state = copy.deepcopy(init_model.state_dict())
    del init_model
    print("  Initial state saved (shared across Phase A and A0)")

    # --------------------------------------------------------
    #  3. PHASE A — DONOR (L_triple + L_cycle)
    # --------------------------------------------------------

    print(f"\n{'='*70}")
    print(f"  PHASE A: Donor Training (L_triple + L_cycle, lam_cycle={LAMBDA_CYCLE})")
    print(f"{'='*70}")

    model_A = AlgKGModel(N_ENTITIES, N_REL, n_gen, HEAD_HIDDEN).to(DEVICE)
    model_A.load_state_dict(copy.deepcopy(init_state))
    optimizer_A = torch.optim.Adam(model_A.parameters(), lr=LR)

    hist_A = train_phase(
        model_A, optimizer_A,
        train_triples, train_labels, val_triples, val_labels,
        epochs=PHASE_A_EPOCHS, batch_size=BATCH_SIZE,
        lambda_k=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
        lambda_cycle=LAMBDA_CYCLE,
        con_cycles_t=con_cycles_t, incon_cycles_t=incon_cycles_t,
        margin=MARGIN, grad_clip=GRAD_CLIP, label="Donor",
    )

    adam_A = adam_state_summary(optimizer_A)
    print(f"  Donor Adam state: |m|={adam_A['mean_m']:.6f}, |v|={adam_A['mean_v']:.6f}, "
          f"step={adam_A['step']}")

    # Save checkpoint
    ckpt_A_model = copy.deepcopy(model_A.state_dict())
    ckpt_A_optim = copy.deepcopy(optimizer_A.state_dict())

    # --------------------------------------------------------
    #  4. PHASE A0 — BASELINE (L_triple only)
    # --------------------------------------------------------

    print(f"\n{'='*70}")
    print(f"  PHASE A0: Baseline Training (L_triple only)")
    print(f"{'='*70}")

    model_A0 = AlgKGModel(N_ENTITIES, N_REL, n_gen, HEAD_HIDDEN).to(DEVICE)
    model_A0.load_state_dict(copy.deepcopy(init_state))
    optimizer_A0 = torch.optim.Adam(model_A0.parameters(), lr=LR)

    hist_A0 = train_phase(
        model_A0, optimizer_A0,
        train_triples, train_labels, val_triples, val_labels,
        epochs=PHASE_A_EPOCHS, batch_size=BATCH_SIZE,
        lambda_k=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
        lambda_cycle=0.0,
        con_cycles_t=con_cycles_t, incon_cycles_t=incon_cycles_t,
        margin=MARGIN, grad_clip=GRAD_CLIP, label="Baseline",
    )

    adam_A0 = adam_state_summary(optimizer_A0)
    print(f"  Baseline Adam state: |m|={adam_A0['mean_m']:.6f}, |v|={adam_A0['mean_v']:.6f}, "
          f"step={adam_A0['step']}")

    ckpt_A0_model = copy.deepcopy(model_A0.state_dict())
    ckpt_A0_optim = copy.deepcopy(optimizer_A0.state_dict())

    # --------------------------------------------------------
    #  5. CHECKPOINT EVALUATION (sanity check)
    # --------------------------------------------------------

    print(f"\n{'='*70}")
    print(f"  CHECKPOINT EVALUATION (after Phase A/A0, before Phase B)")
    print(f"{'='*70}")

    eval_args = (cycle_data, con_paths, incon_paths, gt_con_hol, gt_incon_hol)

    print("\n  --- Donor (trained WITH cycle loss) ---")
    eval_donor = evaluate_holonomy(model_A, *eval_args, label="Donor")
    donor_val_acc = hist_A['best_val_acc']

    print("\n  --- Baseline (trained WITHOUT cycle loss) ---")
    eval_baseline = evaluate_holonomy(model_A0, *eval_args, label="Baseline")
    baseline_val_acc = hist_A0['best_val_acc']

    # Check donor's holonomy separation in training
    final_sep_A = hist_A['hol_incon'][-1] / max(hist_A['hol_con'][-1], 1e-10)
    final_sep_A0 = hist_A0['hol_incon'][-1] / max(hist_A0['hol_con'][-1], 1e-10)
    print(f"\n  Donor final hol_sep (training monitor):   {final_sep_A:.2f}x  "
          f"(con={hist_A['hol_con'][-1]:.4f}, incon={hist_A['hol_incon'][-1]:.4f})")
    print(f"  Baseline final hol_sep (training monitor): {final_sep_A0:.2f}x  "
          f"(con={hist_A0['hol_con'][-1]:.4f}, incon={hist_A0['hol_incon'][-1]:.4f})")

    # Free Phase A/A0 models to save memory
    del model_A, optimizer_A, model_A0, optimizer_A0

    # --------------------------------------------------------
    #  6. PHASE B — CONTINUATION VARIANTS (L_triple only)
    # --------------------------------------------------------

    def run_phase_b(model_state, optim_state, variant_label):
        """Create model+optimizer from saved state, continue training."""
        model = AlgKGModel(N_ENTITIES, N_REL, n_gen, HEAD_HIDDEN).to(DEVICE)
        model.load_state_dict(copy.deepcopy(model_state))
        optimizer = torch.optim.Adam(model.parameters(), lr=LR)
        if optim_state is not None:
            optimizer.load_state_dict(copy.deepcopy(optim_state))

        hist = train_phase(
            model, optimizer,
            train_triples, train_labels, val_triples, val_labels,
            epochs=PHASE_B_EPOCHS, batch_size=BATCH_SIZE,
            lambda_k=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
            lambda_cycle=0.0,  # NO cycle loss in Phase B
            con_cycles_t=con_cycles_t, incon_cycles_t=incon_cycles_t,
            grad_clip=GRAD_CLIP, label=variant_label,
        )
        return model, hist

    print(f"\n{'='*70}")
    print(f"  PHASE B: Continuation ({PHASE_B_EPOCHS} ep, L_triple only)")
    print(f"{'='*70}")

    # B1: Control — baseline params + baseline Adam
    print(f"\n--- B1: Control (theta_B, adam_B) ---")
    model_b1, hist_b1 = run_phase_b(ckpt_A0_model, ckpt_A0_optim, "B1:Control")

    # B2: Params-only — donor params + fresh Adam
    print(f"\n--- B2: Params-only (theta_A, fresh adam) ---")
    model_b2, hist_b2 = run_phase_b(ckpt_A_model, None, "B2:Params")

    # B3: Full state — donor params + donor Adam (KEY TEST)
    print(f"\n--- B3: Full state (theta_A, adam_A) --- [KEY TEST] ---")
    model_b3, hist_b3 = run_phase_b(ckpt_A_model, ckpt_A_optim, "B3:Full")

    # B4: Adam-only — baseline params + donor Adam
    print(f"\n--- B4: Adam-only (theta_B, adam_A) --- [PROVOCATIVE] ---")
    model_b4, hist_b4 = run_phase_b(ckpt_A0_model, ckpt_A_optim, "B4:Adam")

    # --------------------------------------------------------
    #  7. PHASE B EVALUATION
    # --------------------------------------------------------

    print(f"\n{'='*70}")
    print(f"  PHASE B: Final Evaluation")
    print(f"{'='*70}")

    print()
    eval_b1 = evaluate_holonomy(model_b1, *eval_args, label="B1:Control")
    eval_b2 = evaluate_holonomy(model_b2, *eval_args, label="B2:Params")
    eval_b3 = evaluate_holonomy(model_b3, *eval_args, label="B3:Full")
    eval_b4 = evaluate_holonomy(model_b4, *eval_args, label="B4:Adam")

    # Val accuracy
    val_t = torch.tensor(val_triples, dtype=torch.long, device=DEVICE)
    val_l = torch.tensor(val_labels, dtype=torch.long, device=DEVICE)
    with torch.no_grad():
        acc_b1 = (model_b1(val_t).argmax(1) == val_l).float().mean().item()
        acc_b2 = (model_b2(val_t).argmax(1) == val_l).float().mean().item()
        acc_b3 = (model_b3(val_t).argmax(1) == val_l).float().mean().item()
        acc_b4 = (model_b4(val_t).argmax(1) == val_l).float().mean().item()

    # --------------------------------------------------------
    #  8. SUMMARY
    # --------------------------------------------------------

    wall_time = time.time() - t0

    print(f"\n{'='*70}")
    print(f"  SUMMARY — Exp3: Adam State Transfer for Holonomy Basin")
    print(f"{'='*70}")

    print(f"\n  {'Variant':20s} | {'Params':8s} | {'Adam':8s} | "
          f"{'Val Acc':8s} | {'Raw AUC':8s} | {'Path AUC':9s} | {'GT Corr':8s}")
    print(f"  {'-'*20}-+-{'-'*8}-+-{'-'*8}-+-{'-'*8}-+-{'-'*8}-+-{'-'*9}-+-{'-'*8}")

    def _row(label, params_src, adam_src, val_acc, ev):
        print(f"  {label:20s} | {params_src:8s} | {adam_src:8s} | "
              f"{val_acc:.4f}   | {ev['raw_auc']:.4f}   | {ev['path_auc']:.4f}    | {ev['correlation']:.4f}")

    _row("Phase A (Donor)", "init", "evolved", donor_val_acc, eval_donor)
    _row("Phase A0 (Baseline)", "init", "evolved", baseline_val_acc, eval_baseline)
    _row("B1: Control", "theta_B", "adam_B", acc_b1, eval_b1)
    _row("B2: Params-only", "theta_A", "fresh", acc_b2, eval_b2)
    _row("B3: Full state", "theta_A", "adam_A", acc_b3, eval_b3)
    _row("B4: Adam-only", "theta_B", "adam_A", acc_b4, eval_b4)

    print(f"\n  Key Comparisons (F4 analog):")
    d_b3_b2 = eval_b3['path_auc'] - eval_b2['path_auc']
    d_b3_b1 = eval_b3['path_auc'] - eval_b1['path_auc']
    d_b4_b1 = eval_b4['path_auc'] - eval_b1['path_auc']
    d_b2_b1 = eval_b2['path_auc'] - eval_b1['path_auc']

    print(f"    B3 - B2 (Adam state effect on donor params):      {d_b3_b2:+.4f}")
    print(f"    B3 - B1 (Full transfer vs control):               {d_b3_b1:+.4f}")
    print(f"    B4 - B1 (Adam-only transfer vs control):          {d_b4_b1:+.4f}")
    print(f"    B2 - B1 (Params-only transfer vs control):        {d_b2_b1:+.4f}")

    print(f"\n  Holonomy Separation (hol_incon / hol_con) at Phase B end:")
    for lbl, h in [("B1:Control", hist_b1), ("B2:Params", hist_b2),
                    ("B3:Full", hist_b3), ("B4:Adam", hist_b4)]:
        sep = h['hol_incon'][-1] / max(h['hol_con'][-1], 1e-10)
        print(f"    {lbl:12s}: {sep:.4f}x  (con={h['hol_con'][-1]:.4f}, incon={h['hol_incon'][-1]:.4f})")

    print(f"\n  Adam State Diagnostics:")
    print(f"    Donor:    |m|={adam_A['mean_m']:.6f}, |v|={adam_A['mean_v']:.6f}, step={adam_A['step']}")
    print(f"    Baseline: |m|={adam_A0['mean_m']:.6f}, |v|={adam_A0['mean_v']:.6f}, step={adam_A0['step']}")

    print(f"\n  Phase A Donor — Holonomy trajectory:")
    for ep_idx in [0, 99, 199, 399, 599, 799]:
        if ep_idx < len(hist_A['hol_con']):
            sep = hist_A['hol_incon'][ep_idx] / max(hist_A['hol_con'][ep_idx], 1e-10)
            print(f"    ep {ep_idx+1:4d}: con={hist_A['hol_con'][ep_idx]:.4f}  "
                  f"incon={hist_A['hol_incon'][ep_idx]:.4f}  sep={sep:.2f}x")

    print(f"\n  Wall time: {wall_time:.1f} s")

    # --------------------------------------------------------
    #  9. SAVE
    # --------------------------------------------------------

    out_path = os.path.join(os.path.dirname(__file__), 'exp3_results.npz')
    np.savez(out_path, **{
        # Phase A
        'donor_val_acc': np.array(hist_A['val_acc']),
        'donor_hol_con': np.array(hist_A['hol_con']),
        'donor_hol_incon': np.array(hist_A['hol_incon']),
        'donor_tk': np.array(hist_A['tk']),
        'donor_raw_auc': eval_donor['raw_auc'],
        'donor_path_auc': eval_donor['path_auc'],
        'donor_corr': eval_donor['correlation'],
        # Phase A0
        'baseline_val_acc': np.array(hist_A0['val_acc']),
        'baseline_hol_con': np.array(hist_A0['hol_con']),
        'baseline_hol_incon': np.array(hist_A0['hol_incon']),
        'baseline_raw_auc': eval_baseline['raw_auc'],
        'baseline_path_auc': eval_baseline['path_auc'],
        'baseline_corr': eval_baseline['correlation'],
        # Phase B
        'b1_raw_auc': eval_b1['raw_auc'], 'b1_path_auc': eval_b1['path_auc'],
        'b1_corr': eval_b1['correlation'], 'b1_val_acc': acc_b1,
        'b2_raw_auc': eval_b2['raw_auc'], 'b2_path_auc': eval_b2['path_auc'],
        'b2_corr': eval_b2['correlation'], 'b2_val_acc': acc_b2,
        'b3_raw_auc': eval_b3['raw_auc'], 'b3_path_auc': eval_b3['path_auc'],
        'b3_corr': eval_b3['correlation'], 'b3_val_acc': acc_b3,
        'b4_raw_auc': eval_b4['raw_auc'], 'b4_path_auc': eval_b4['path_auc'],
        'b4_corr': eval_b4['correlation'], 'b4_val_acc': acc_b4,
        # B histories
        'b1_hol_con': np.array(hist_b1['hol_con']),
        'b1_hol_incon': np.array(hist_b1['hol_incon']),
        'b3_hol_con': np.array(hist_b3['hol_con']),
        'b3_hol_incon': np.array(hist_b3['hol_incon']),
        # Meta
        'wall_time': wall_time,
        'adam_A_mean_m': adam_A['mean_m'],
        'adam_A_mean_v': adam_A['mean_v'],
        'adam_A0_mean_m': adam_A0['mean_m'],
        'adam_A0_mean_v': adam_A0['mean_v'],
    })
    print(f"  Results saved to {out_path}")


if __name__ == '__main__':
    run_experiment()
