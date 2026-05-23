"""
Experiment 1b — Holonomy Consistency Detection (Revised)
=========================================================
Fixes from Exp1 diagnosis (G_RAG-3):
  1. K = 24 relations (12 base + 12 closure) vs K=6 before
  2. angle_scale = 0.25 (operators close to I) vs 0.8 before
  3. Explicit closure relations with GT holonomy = 0.0 exact

This gives a clean binary signal: consistent cycles have holonomy ≡ 0,
inconsistent cycles have holonomy >> 0. The question is whether the
learned operators can reconstruct this separation.
"""

import sys, os, time
import torch
import torch.nn as nn
import numpy as np
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(__file__))
from kg_utils import (
    generate_synthetic_kg,
    generate_relation_operators_with_closures,
    generate_cycles_v2,
    compute_holonomy,
    relation_killing_tension,
    relation_ortho_tension,
    DEVICE,
)

sys.stdout.reconfigure(encoding='utf-8')


# ============================================================
#  MODELS (identical to Exp1 — import would be cleaner but
#  we keep them self-contained for reproducibility)
# ============================================================

class AlgKGModel(nn.Module):
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


class MLPKGModel(nn.Module):
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


# ============================================================
#  TRAINING
# ============================================================

def train_model(
    model, train_triples, train_labels, val_triples, val_labels,
    epochs=1000, lr=1e-3, batch_size=256,
    lambda_k=0.0, lambda_ortho=0.0,
    grad_clip=1.0, label="Model",
):
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    criterion = nn.CrossEntropyLoss()

    train_t = torch.tensor(train_triples, dtype=torch.long, device=DEVICE)
    train_l = torch.tensor(train_labels, dtype=torch.long, device=DEVICE)
    val_t = torch.tensor(val_triples, dtype=torch.long, device=DEVICE)
    val_l = torch.tensor(val_labels, dtype=torch.long, device=DEVICE)

    history = {
        'train_acc': [], 'val_acc': [], 'loss': [],
        'tk': [], 'ortho_err': [],
    }
    best_val_acc = 0.0
    n = len(train_triples)

    for ep in range(1, epochs + 1):
        model.train()
        perm = torch.randperm(n, device=DEVICE)
        epoch_loss = 0.0
        epoch_correct = 0

        for i in range(0, n, batch_size):
            idx = perm[i:i+batch_size]
            batch_t = train_t[idx]
            batch_l = train_l[idx]

            logits = model(batch_t)
            task_loss = criterion(logits, batch_l)

            reg_loss = torch.tensor(0.0, dtype=torch.float64, device=DEVICE)
            if hasattr(model, 'regularization_loss') and (lambda_k > 0 or lambda_ortho > 0):
                reg_loss = model.regularization_loss(lambda_k, lambda_ortho)

            loss = task_loss + reg_loss
            optimizer.zero_grad()
            loss.backward()
            if grad_clip > 0:
                nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
            optimizer.step()

            epoch_loss += task_loss.item() * len(idx)
            epoch_correct += (logits.argmax(1) == batch_l).sum().item()

        train_acc = epoch_correct / n

        model.eval()
        with torch.no_grad():
            val_logits = model(val_t)
            val_preds = val_logits.argmax(1)
            val_acc = (val_preds == val_l).float().mean().item()

            ops = model.get_relation_matrices()
            I_n = torch.eye(ops.shape[1], dtype=ops.dtype, device=ops.device)
            ortho_err = sum(
                ((ops[r] @ ops[r].T - I_n) ** 2).sum().item()
                for r in range(ops.shape[0])
            ) / ops.shape[0]

            try:
                tk = relation_killing_tension(ops).item()
            except Exception:
                tk = float('nan')

        history['train_acc'].append(train_acc)
        history['val_acc'].append(val_acc)
        history['loss'].append(epoch_loss / n)
        history['tk'].append(tk)
        history['ortho_err'].append(ortho_err)

        if val_acc > best_val_acc:
            best_val_acc = val_acc

        if ep % 100 == 0 or ep == 1:
            print(f"  [{label}] ep {ep:4d}  train={train_acc:.4f}  val={val_acc:.4f}  "
                  f"loss={epoch_loss/n:.4f}  T_K={tk:.4f}  OrthoErr={ortho_err:.4f}")

    history['best_val_acc'] = best_val_acc
    return history


# ============================================================
#  HOLONOMY EVALUATION
# ============================================================

def evaluate_holonomy(model, cycle_data, label="Model"):
    model.eval()
    with torch.no_grad():
        ops = model.get_relation_matrices()
        hol_con = compute_holonomy(ops, cycle_data['consistent_cycles']).cpu().numpy()
        hol_incon = compute_holonomy(ops, cycle_data['inconsistent_cycles']).cpu().numpy()

    n_con = len(hol_con)
    n_incon = len(hol_incon)
    labels = np.concatenate([np.zeros(n_con), np.ones(n_incon)])
    scores = np.concatenate([hol_con, hol_incon])

    try:
        auc = roc_auc_score(labels, scores)
    except ValueError:
        auc = float('nan')

    gt_con = cycle_data['consistent_holonomy']
    gt_incon = cycle_data['inconsistent_holonomy']
    gt_all = np.concatenate([gt_con, gt_incon])
    if len(gt_all) > 2 and np.std(gt_all) > 1e-10 and np.std(scores) > 1e-10:
        corr = np.corrcoef(gt_all, scores)[0, 1]
    else:
        corr = float('nan')

    results = {
        'auc': auc,
        'correlation': corr,
        'hol_consistent_mean': float(np.mean(hol_con)),
        'hol_consistent_std': float(np.std(hol_con)),
        'hol_inconsistent_mean': float(np.mean(hol_incon)),
        'hol_inconsistent_std': float(np.std(hol_incon)),
        'gt_consistent_mean': float(np.mean(gt_con)),
        'gt_inconsistent_mean': float(np.mean(gt_incon)),
        'separation_ratio': float(np.mean(hol_incon) / max(np.mean(hol_con), 1e-10)),
    }

    print(f"\n  [{label}] Holonomy Evaluation:")
    print(f"    Consistent  : Δ_info = {results['hol_consistent_mean']:.4f} ± {results['hol_consistent_std']:.4f}"
          f"  (GT: {results['gt_consistent_mean']:.6f})")
    print(f"    Inconsistent: Δ_info = {results['hol_inconsistent_mean']:.4f} ± {results['hol_inconsistent_std']:.4f}"
          f"  (GT: {results['gt_inconsistent_mean']:.4f})")
    print(f"    Separation ratio: {results['separation_ratio']:.2f}×")
    print(f"    ROC-AUC: {results['auc']:.4f}")
    print(f"    Correlation (learned vs GT holonomy): {results['correlation']:.4f}")

    return results


# ============================================================
#  MAIN EXPERIMENT
# ============================================================

def run_experiment(seed=42):
    # ---- Config (Exp1b changes vs Exp1) ----
    N_BASE = 12          # was: n_relations=6 (all base)
    N_CLOSURES = 12      # new: explicit inverse-product relations
    N_REL = N_BASE + N_CLOSURES  # = 24 total
    ANGLE_SCALE = 0.25   # was: 0.8
    N_ENTITIES = 200
    N_TRAIN = 8000       # slightly more triples for 24 relations
    N_VAL = 2000
    EPOCHS = 800
    LR = 1e-3
    HEAD_HIDDEN = 32
    LAMBDA_K = 0.5
    LAMBDA_ORTHO = 1.0
    D = 3

    print("=" * 70)
    print("  Exp1b — Holonomy Consistency (Revised: closures + small angles)")
    print("=" * 70)
    print(f"  Device: {DEVICE}")
    print(f"  Seed: {seed}")
    print(f"  Config: {N_BASE} base + {N_CLOSURES} closure = {N_REL} relations, "
          f"angle_scale={ANGLE_SCALE}")

    t0 = time.time()

    # --- 1. Generate relation operators with closures ---
    print("\n--- Generating relation operators with closures ---")
    Q_ops, alphas, closure_cycles = generate_relation_operators_with_closures(
        n_base=N_BASE, n_closures=N_CLOSURES, d=D,
        angle_scale=ANGLE_SCALE, seed=seed,
    )
    n_gen = Q_ops.shape[1]
    print(f"  Total relations: {Q_ops.shape[0]}, n_gen: {n_gen}")
    print(f"  Closure cycles (GT holonomy = 0): {len(closure_cycles)}")

    # Verify closure holonomy
    I_n = np.eye(n_gen)
    closure_hols = []
    for ri, rj, rc in closure_cycles:
        H = Q_ops[rc] @ Q_ops[rj] @ Q_ops[ri]
        closure_hols.append(np.linalg.norm(H - I_n, 'fro'))
    print(f"  Closure holonomy check: max={max(closure_hols):.2e}, mean={np.mean(closure_hols):.2e}")

    # --- 2. Generate synthetic KG using these operators ---
    print("\n--- Generating synthetic KG ---")
    kg = generate_synthetic_kg(
        n_entities=N_ENTITIES, n_relations=N_REL,
        n_train_triples=N_TRAIN, n_val_triples=N_VAL,
        d=D, angle_scale=ANGLE_SCALE, seed=seed,
    )
    # Override the relation_ops with our closure-augmented set
    # (generate_synthetic_kg generated its own; we need ours for cycle evaluation)
    # Actually we need to regenerate the KG using OUR operators.
    # Let's rebuild the KG data manually.

    rng = np.random.default_rng(seed)
    entities = rng.standard_normal((N_ENTITIES, n_gen))
    entities = entities / np.linalg.norm(entities, axis=1, keepdims=True)

    # Compute transformed entities for all (relation, head) pairs
    transformed = np.einsum('rij,hj->rhi', Q_ops, entities)  # (n_rel, n_ent, n_gen)

    nearest_tail = np.zeros((N_REL, N_ENTITIES), dtype=np.int64)
    nearest_dist = np.zeros((N_REL, N_ENTITIES), dtype=np.float64)
    for r in range(N_REL):
        for h in range(N_ENTITIES):
            diffs = entities - transformed[r, h]
            dists = np.linalg.norm(diffs, axis=1)
            dists[h] = 1e10
            nearest_tail[r, h] = np.argmin(dists)
            nearest_dist[r, h] = dists[nearest_tail[r, h]]

    print(f"  Entities: {N_ENTITIES}, Relations: {N_REL}, n_gen: {n_gen}")
    print(f"  Nearest-neighbor distance: mean={nearest_dist.mean():.4f}, "
          f"max={nearest_dist.max():.4f}")

    def _generate_split(n_triples, rng_local):
        triples, labels = [], []
        n_valid = n_triples // 2
        n_invalid = n_triples - n_valid
        for _ in range(n_valid):
            r = rng_local.integers(N_REL)
            h = rng_local.integers(N_ENTITIES)
            t = nearest_tail[r, h]
            triples.append([h, r, t])
            labels.append(1)
        for _ in range(n_invalid):
            r = rng_local.integers(N_REL)
            h = rng_local.integers(N_ENTITIES)
            t = rng_local.integers(N_ENTITIES)
            while t == nearest_tail[r, h] or t == h:
                t = rng_local.integers(N_ENTITIES)
            triples.append([h, r, t])
            labels.append(0)
        perm = rng_local.permutation(len(triples))
        return np.array(triples)[perm], np.array(labels)[perm]

    train_rng = np.random.default_rng(seed + 100)
    val_rng = np.random.default_rng(seed + 200)
    train_triples, train_labels = _generate_split(N_TRAIN, train_rng)
    val_triples, val_labels = _generate_split(N_VAL, val_rng)

    print(f"  Train triples: {len(train_triples)} "
          f"(valid: {train_labels.sum()}, invalid: {(1-train_labels).sum()})")
    print(f"  Val triples: {len(val_triples)} "
          f"(valid: {val_labels.sum()}, invalid: {(1-val_labels).sum()})")

    # --- 3. Generate cycles ---
    print("\n--- Generating cycles ---")
    cycle_data = generate_cycles_v2(
        Q_ops, closure_cycles,
        n_consistent=200, n_inconsistent=200, seed=seed,
    )
    print(f"  Consistent cycles: {len(cycle_data['consistent_cycles'])}, "
          f"  GT holonomy mean: {cycle_data['consistent_holonomy'].mean():.6f}")
    print(f"  Inconsistent cycles: {len(cycle_data['inconsistent_cycles'])}, "
          f"  GT holonomy mean: {cycle_data['inconsistent_holonomy'].mean():.4f}")
    gt_sep = cycle_data['inconsistent_holonomy'].mean() / max(cycle_data['consistent_holonomy'].mean(), 1e-10)
    print(f"  GT separation ratio: {gt_sep:.1f}×")

    # --- 4. Train models ---
    print(f"\n--- Training AlgKG (λ_K={LAMBDA_K}, λ_ortho={LAMBDA_ORTHO}) ---")
    algkg = AlgKGModel(N_ENTITIES, N_REL, n_gen, HEAD_HIDDEN).to(DEVICE)
    hist_alg = train_model(
        algkg, train_triples, train_labels,
        val_triples, val_labels,
        epochs=EPOCHS, lr=LR, lambda_k=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
        label="AlgKG",
    )

    print(f"\n--- Training MLP (no algebraic constraints) ---")
    mlp = MLPKGModel(N_ENTITIES, N_REL, n_gen, HEAD_HIDDEN).to(DEVICE)
    hist_mlp = train_model(
        mlp, train_triples, train_labels,
        val_triples, val_labels,
        epochs=EPOCHS, lr=LR,
        label="MLP",
    )

    # --- 5. Holonomy evaluation ---
    print("\n--- Holonomy Evaluation ---")
    hol_alg = evaluate_holonomy(algkg, cycle_data, "AlgKG")
    hol_mlp = evaluate_holonomy(mlp, cycle_data, "MLP")

    # --- 6. Summary ---
    wall_time = time.time() - t0
    print("\n" + "=" * 70)
    print("  SUMMARY — Exp1b")
    print("=" * 70)

    print(f"\n  Config changes vs Exp1:")
    print(f"    Relations: 6 → {N_REL} ({N_BASE} base + {N_CLOSURES} closures)")
    print(f"    angle_scale: 0.8 → {ANGLE_SCALE}")
    print(f"    GT separation: 1.13× → {gt_sep:.1f}×")

    print(f"\n  Triple Classification:")
    print(f"    AlgKG best val acc: {hist_alg['best_val_acc']:.4f}")
    print(f"    MLP   best val acc: {hist_mlp['best_val_acc']:.4f}")
    print(f"    Δ(AlgKG − MLP):    {hist_alg['best_val_acc'] - hist_mlp['best_val_acc']:+.4f}")

    print(f"\n  Holonomy Consistency Detection (ROC-AUC):")
    print(f"    AlgKG: {hol_alg['auc']:.4f}  (separation: {hol_alg['separation_ratio']:.2f}×)")
    print(f"    MLP:   {hol_mlp['auc']:.4f}  (separation: {hol_mlp['separation_ratio']:.2f}×)")

    print(f"\n  Learned-vs-GT Holonomy Correlation:")
    print(f"    AlgKG: {hol_alg['correlation']:.4f}")
    print(f"    MLP:   {hol_mlp['correlation']:.4f}")

    print(f"\n  AlgKG Final Diagnostics:")
    print(f"    T_K (final):       {hist_alg['tk'][-1]:.4f}")
    print(f"    OrthoErr (final):  {hist_alg['ortho_err'][-1]:.4f}")

    print(f"\n  MLP Final Diagnostics:")
    print(f"    OrthoErr (final):  {hist_mlp['ortho_err'][-1]:.4f}")

    print(f"\n  Wall time: {wall_time:.1f} s")

    # --- 7. Save results ---
    out_path = os.path.join(os.path.dirname(__file__), 'exp1b_results.npz')
    np.savez(out_path, **{
        'alg_val_acc': np.array(hist_alg['val_acc']),
        'mlp_val_acc': np.array(hist_mlp['val_acc']),
        'alg_tk': np.array(hist_alg['tk']),
        'alg_ortho': np.array(hist_alg['ortho_err']),
        'mlp_ortho': np.array(hist_mlp['ortho_err']),
        'alg_hol_con': hol_alg['hol_consistent_mean'],
        'alg_hol_incon': hol_alg['hol_inconsistent_mean'],
        'alg_auc': hol_alg['auc'],
        'mlp_auc': hol_mlp['auc'],
        'alg_corr': hol_alg['correlation'],
        'mlp_corr': hol_mlp['correlation'],
        'gt_sep': gt_sep,
    })
    print(f"\n  Results saved to {out_path}")

    return {
        'hist_alg': hist_alg, 'hist_mlp': hist_mlp,
        'hol_alg': hol_alg, 'hol_mlp': hol_mlp,
        'gt_sep': gt_sep, 'wall_time': wall_time,
    }


if __name__ == '__main__':
    run_experiment(seed=42)
