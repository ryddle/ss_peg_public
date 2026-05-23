"""
Experiment 1 — Holonomy Consistency Detection on Synthetic Algebraic KG
========================================================================
Phase 1 of Algebraic Graph RAG: can T_K crystallization on relation
operators produce holonomy that discriminates consistent vs inconsistent
cycles in a knowledge graph?

Design:
  - Synthetic KG with SU(3) adjoint-representation relation operators
  - Task: binary triple classification (valid / invalid)
  - Models: AlgKG (T_K + T_ortho) vs MLP (unconstrained)
  - Post-training: compute holonomy on cycles, measure ROC-AUC

Connects to: Exp4 (SU(3) cubic task), G1 (manifold matters),
             G7 (T_ortho necessary), G18 (decoder matters)
"""

import sys, os, time
import torch
import torch.nn as nn
import numpy as np
from sklearn.metrics import roc_auc_score, accuracy_score

sys.path.insert(0, os.path.dirname(__file__))
from kg_utils import (
    generate_synthetic_kg,
    generate_cycles,
    compute_holonomy,
    relation_killing_tension,
    relation_ortho_tension,
    DEVICE,
)

# Force UTF-8 output
sys.stdout.reconfigure(encoding='utf-8')


# ============================================================
#  MODELS
# ============================================================

class AlgKGModel(nn.Module):
    """
    Algebraic KG model with T_K + T_ortho regularized relation operators.

    Entity embeddings: learnable vectors in R^{n_gen}.
    Relation operators: learnable θ_r ∈ R^{n_gen × n_gen}, regularized toward O(n_gen).
    Scoring: ||θ_r · e_h - e_t||² → 2-layer classifier → valid/invalid.
    """
    def __init__(self, n_entities, n_relations, n_gen, head_hidden=32):
        super().__init__()
        self.n_gen = n_gen
        self.entity_emb = nn.Embedding(n_entities, n_gen, dtype=torch.float64)
        # Relation operators: initialized near identity
        self.relation_ops = nn.Parameter(
            torch.eye(n_gen, dtype=torch.float64).unsqueeze(0).repeat(n_relations, 1, 1)
            + 0.05 * torch.randn(n_relations, n_gen, n_gen, dtype=torch.float64)
        )
        # Head: distance → prediction (following G18: 2-layer + ReLU)
        self.head = nn.Sequential(
            nn.Linear(n_gen + 1, head_hidden, dtype=torch.float64),
            nn.ReLU(),
            nn.Linear(head_hidden, 2, dtype=torch.float64),
        )
        # Initialize entity embeddings as unit vectors
        nn.init.normal_(self.entity_emb.weight)
        with torch.no_grad():
            self.entity_emb.weight.div_(
                self.entity_emb.weight.norm(dim=1, keepdim=True).clamp(min=1e-8)
            )

    def forward(self, triples):
        """
        triples : (batch, 3) int [head_idx, rel_idx, tail_idx]
        Returns : (batch, 2) logits
        """
        h_idx = triples[:, 0]
        r_idx = triples[:, 1]
        t_idx = triples[:, 2]

        e_h = self.entity_emb(h_idx)  # (batch, n_gen)
        e_t = self.entity_emb(t_idx)  # (batch, n_gen)

        # θ_r · e_h for each triple in the batch
        theta_r = self.relation_ops[r_idx]  # (batch, n_gen, n_gen)
        transformed = torch.bmm(theta_r, e_h.unsqueeze(-1)).squeeze(-1)  # (batch, n_gen)

        # Feature: difference vector + distance scalar
        diff = transformed - e_t  # (batch, n_gen)
        dist = diff.norm(dim=1, keepdim=True)  # (batch, 1)
        features = torch.cat([diff, dist], dim=1)  # (batch, n_gen+1)

        return self.head(features)

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
    """
    MLP baseline: unconstrained relation matrices, same architecture otherwise.
    """
    def __init__(self, n_entities, n_relations, n_gen, head_hidden=32):
        super().__init__()
        self.n_gen = n_gen
        self.entity_emb = nn.Embedding(n_entities, n_gen, dtype=torch.float64)
        # Unconstrained relation matrices
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
        h_idx = triples[:, 0]
        r_idx = triples[:, 1]
        t_idx = triples[:, 2]

        e_h = self.entity_emb(h_idx)
        e_t = self.entity_emb(t_idx)
        theta_r = self.relation_ops[r_idx]
        transformed = torch.bmm(theta_r, e_h.unsqueeze(-1)).squeeze(-1)
        diff = transformed - e_t
        dist = diff.norm(dim=1, keepdim=True)
        features = torch.cat([diff, dist], dim=1)
        return self.head(features)

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
    """Train triple classifier and return training history."""
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

        # Validation
        model.eval()
        with torch.no_grad():
            val_logits = model(val_t)
            val_preds = val_logits.argmax(1)
            val_acc = (val_preds == val_l).float().mean().item()

            # Ortho error
            ops = model.get_relation_matrices()
            I_n = torch.eye(ops.shape[1], dtype=ops.dtype, device=ops.device)
            ortho_err = sum(
                ((ops[r] @ ops[r].T - I_n) ** 2).sum().item()
                for r in range(ops.shape[0])
            ) / ops.shape[0]

            # T_K (only if we can compute it)
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

def evaluate_holonomy(
    model, cycle_data, label="Model",
):
    """
    Compute holonomy on consistent and inconsistent cycles.
    Returns ROC-AUC and statistics.
    """
    model.eval()
    with torch.no_grad():
        ops = model.get_relation_matrices()

        con_cycles = cycle_data['consistent_cycles']
        incon_cycles = cycle_data['inconsistent_cycles']

        hol_con = compute_holonomy(ops, con_cycles).cpu().numpy()
        hol_incon = compute_holonomy(ops, incon_cycles).cpu().numpy()

    # Labels: 0 = consistent (low holonomy expected), 1 = inconsistent (high holonomy)
    n_con = len(hol_con)
    n_incon = len(hol_incon)
    labels = np.concatenate([np.zeros(n_con), np.ones(n_incon)])
    scores = np.concatenate([hol_con, hol_incon])

    try:
        auc = roc_auc_score(labels, scores)
    except ValueError:
        auc = float('nan')

    # Ground truth holonomy
    gt_con = cycle_data['consistent_holonomy']
    gt_incon = cycle_data['inconsistent_holonomy']

    # Correlation between learned and ground-truth holonomy
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
          f"  (GT: {results['gt_consistent_mean']:.4f})")
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
    print("=" * 70)
    print("  Exp1 — Holonomy Consistency Detection on Synthetic Algebraic KG")
    print("=" * 70)
    print(f"  Device: {DEVICE}")
    print(f"  Seed: {seed}")

    t0 = time.time()

    # --- 1. Generate synthetic KG ---
    print("\n--- Generating synthetic KG ---")
    kg = generate_synthetic_kg(
        n_entities=200, n_relations=6, n_train_triples=5000,
        n_val_triples=1000, d=3, angle_scale=0.8, seed=seed,
    )
    n_gen = kg['n_gen']
    n_entities = kg['entity_embeddings'].shape[0]
    n_relations = kg['relation_ops'].shape[0]

    print(f"  Entities: {n_entities}, Relations: {n_relations}, n_gen: {n_gen}")
    print(f"  Train triples: {len(kg['train_triples'])} "
          f"(valid: {kg['train_labels'].sum()}, invalid: {(1-kg['train_labels']).sum()})")
    print(f"  Val triples: {len(kg['val_triples'])} "
          f"(valid: {kg['val_labels'].sum()}, invalid: {(1-kg['val_labels']).sum()})")

    # Check nearest distance statistics
    print(f"  Nearest-neighbor distance: mean={kg['nearest_dist'].mean():.4f}, "
          f"max={kg['nearest_dist'].max():.4f}")

    # --- 2. Generate cycles ---
    print("\n--- Generating cycles ---")
    cycle_data = generate_cycles(
        kg['relation_ops'], n_consistent=200, n_inconsistent=200, seed=seed,
    )
    print(f"  Consistent cycles: {len(cycle_data['consistent_cycles'])}, "
          f"  GT holonomy mean: {cycle_data['consistent_holonomy'].mean():.4f}")
    print(f"  Inconsistent cycles: {len(cycle_data['inconsistent_cycles'])}, "
          f"  GT holonomy mean: {cycle_data['inconsistent_holonomy'].mean():.4f}")

    gt_sep = cycle_data['inconsistent_holonomy'].mean() / max(cycle_data['consistent_holonomy'].mean(), 1e-10)
    print(f"  GT separation ratio: {gt_sep:.2f}×")

    # --- 3. Train models ---
    EPOCHS = 800
    LR = 1e-3
    HEAD_HIDDEN = 32
    LAMBDA_K = 0.5
    LAMBDA_ORTHO = 1.0

    print(f"\n--- Training AlgKG (λ_K={LAMBDA_K}, λ_ortho={LAMBDA_ORTHO}) ---")
    algkg = AlgKGModel(n_entities, n_relations, n_gen, HEAD_HIDDEN).to(DEVICE)
    hist_alg = train_model(
        algkg, kg['train_triples'], kg['train_labels'],
        kg['val_triples'], kg['val_labels'],
        epochs=EPOCHS, lr=LR, lambda_k=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
        label="AlgKG",
    )

    print(f"\n--- Training MLP (no algebraic constraints) ---")
    mlp = MLPKGModel(n_entities, n_relations, n_gen, HEAD_HIDDEN).to(DEVICE)
    hist_mlp = train_model(
        mlp, kg['train_triples'], kg['train_labels'],
        kg['val_triples'], kg['val_labels'],
        epochs=EPOCHS, lr=LR,
        label="MLP",
    )

    # --- 4. Holonomy evaluation ---
    print("\n--- Holonomy Evaluation ---")
    hol_alg = evaluate_holonomy(algkg, cycle_data, "AlgKG")
    hol_mlp = evaluate_holonomy(mlp, cycle_data, "MLP")

    # --- 5. Summary ---
    wall_time = time.time() - t0
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
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

    # --- 6. Save results ---
    results = {
        'hist_alg': hist_alg,
        'hist_mlp': hist_mlp,
        'hol_alg': hol_alg,
        'hol_mlp': hol_mlp,
        'cycle_data_gt': {
            'consistent_holonomy': cycle_data['consistent_holonomy'],
            'inconsistent_holonomy': cycle_data['inconsistent_holonomy'],
        },
        'wall_time': wall_time,
        'config': {
            'epochs': EPOCHS, 'lr': LR, 'head_hidden': HEAD_HIDDEN,
            'lambda_k': LAMBDA_K, 'lambda_ortho': LAMBDA_ORTHO,
            'n_entities': n_entities, 'n_relations': n_relations,
            'n_gen': n_gen, 'd': kg['d'], 'seed': seed,
        },
    }
    out_path = os.path.join(os.path.dirname(__file__), 'exp1_results.npz')
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
    })
    print(f"\n  Results saved to {out_path}")

    return results


if __name__ == '__main__':
    results = run_experiment(seed=42)
