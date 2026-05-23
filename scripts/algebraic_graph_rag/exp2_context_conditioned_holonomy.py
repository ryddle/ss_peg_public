"""
Experiment 2 — Context-Conditioned Operators for Holonomy Detection
=====================================================================
G_RAG-8 diagnosis: shared operators can't encode cycle consistency.
Fix: θ_r(h,t) = θ_r + α · MLP(e_h || e_t), where the base θ_r
crystallizes via T_K but the perturbation carries context-specific info.

Hypothesis: entity-path holonomy with context-conditioned operators
should discriminate consistent vs inconsistent cycles, even though
raw operator holonomy collapses (as in Exp1/1b).

Models:
  1. AlgKG_Ctx: base ops (T_K + T_ortho) + context perturbation
  2. AlgKG:     shared ops (T_K + T_ortho) — Exp1b control
  3. MLP:       shared unconstrained ops — baseline

Metrics:
  A. Triple classification (val accuracy)
  B. Raw operator holonomy: ||Π θ_r - I||_F  (Exp1b-compatible)
  C. Entity-path holonomy:  ||θ_{rN}...θ_{r1}·e_{h1} - e_{h1}||_2  (NEW)
"""

import sys, os, time
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
#  MODELS
# ============================================================

class AlgKGCtxModel(nn.Module):
    """AlgKG with context-conditioned operators: θ_r(h,t) = θ_r + α·MLP(e_h‖e_t)."""

    def __init__(self, n_entities, n_relations, n_gen, head_hidden=32, ctx_hidden=16):
        super().__init__()
        self.n_gen = n_gen
        self.entity_emb = nn.Embedding(n_entities, n_gen, dtype=torch.float64)
        self.relation_ops = nn.Parameter(
            torch.eye(n_gen, dtype=torch.float64).unsqueeze(0).repeat(n_relations, 1, 1)
            + 0.05 * torch.randn(n_relations, n_gen, n_gen, dtype=torch.float64)
        )
        self.ctx_mlp = nn.Sequential(
            nn.Linear(2 * n_gen, ctx_hidden, dtype=torch.float64),
            nn.ReLU(),
            nn.Linear(ctx_hidden, n_gen * n_gen, dtype=torch.float64),
        )
        self.alpha = nn.Parameter(torch.tensor(0.01, dtype=torch.float64))
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
            # Start with zero perturbation so early training = pure AlgKG
            nn.init.zeros_(self.ctx_mlp[2].weight)
            nn.init.zeros_(self.ctx_mlp[2].bias)

    def _get_conditioned_op(self, r_idx, e_h, e_t):
        """θ_r(h,t) = θ_r + α · reshape(MLP(e_h ‖ e_t), (ng, ng))."""
        base = self.relation_ops[r_idx]                        # (batch, ng, ng)
        ctx = torch.cat([e_h, e_t], dim=-1)                    # (batch, 2*ng)
        delta = self.ctx_mlp(ctx).reshape(-1, self.n_gen, self.n_gen)
        return base + self.alpha * delta

    def forward(self, triples):
        h_idx, r_idx, t_idx = triples[:, 0], triples[:, 1], triples[:, 2]
        e_h = self.entity_emb(h_idx)
        e_t = self.entity_emb(t_idx)
        theta = self._get_conditioned_op(r_idx, e_h, e_t)
        transformed = torch.bmm(theta, e_h.unsqueeze(-1)).squeeze(-1)
        diff = transformed - e_t
        dist = diff.norm(dim=1, keepdim=True)
        return self.head(torch.cat([diff, dist], dim=1))

    def get_relation_matrices(self):
        """Base operators only (for raw holonomy and T_K diagnostic)."""
        return self.relation_ops

    def regularization_loss(self, lambda_k=0.5, lambda_ortho=1.0):
        ops = self.relation_ops  # Only regularize base operators
        loss = torch.tensor(0.0, dtype=torch.float64, device=ops.device)
        if lambda_k > 0:
            loss = loss + lambda_k * relation_killing_tension(ops)
        if lambda_ortho > 0:
            loss = loss + lambda_ortho * relation_ortho_tension(ops)
        return loss


class AlgKGModel(nn.Module):
    """AlgKG with shared operators (Exp1b control)."""

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
    """MLP with unconstrained shared operators (baseline)."""

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
    epochs=800, lr=1e-3, batch_size=256,
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

        if hasattr(model, 'alpha'):
            history.setdefault('alpha', []).append(model.alpha.item())

        if val_acc > best_val_acc:
            best_val_acc = val_acc

        if ep % 100 == 0 or ep == 1:
            extra = ""
            if hasattr(model, 'alpha'):
                extra = f"  α={model.alpha.item():.6f}"
            print(f"  [{label}] ep {ep:4d}  train={train_acc:.4f}  val={val_acc:.4f}  "
                  f"loss={epoch_loss/n:.4f}  T_K={tk:.4f}  OrthoErr={ortho_err:.4f}{extra}")

    history['best_val_acc'] = best_val_acc
    return history


# ============================================================
#  EVALUATION — Raw Operator Holonomy (Exp1b-compatible)
# ============================================================

def evaluate_raw_holonomy(model, cycle_data, label="Model"):
    """Raw holonomy: ||Π θ_r - I||_F using base operators only."""
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

    sep = float(np.mean(hol_incon) / max(np.mean(hol_con), 1e-10))

    results = {
        'auc': auc,
        'con_mean': float(np.mean(hol_con)),
        'incon_mean': float(np.mean(hol_incon)),
        'separation': sep,
    }

    print(f"  [{label}] Raw Holonomy:  con={results['con_mean']:.4f}  "
          f"incon={results['incon_mean']:.4f}  sep={sep:.2f}×  AUC={auc:.4f}")

    return results


# ============================================================
#  EVALUATION — Entity-Path Holonomy (NEW)
# ============================================================

def _compute_path_holonomy(model, paths):
    """Batched entity-path holonomy computation."""
    model.eval()
    is_ctx = hasattr(model, '_get_conditioned_op')

    h1s = torch.tensor([p[0] for p in paths], dtype=torch.long, device=DEVICE)
    r1s = torch.tensor([p[1] for p in paths], dtype=torch.long, device=DEVICE)
    h2s = torch.tensor([p[2] for p in paths], dtype=torch.long, device=DEVICE)
    r2s = torch.tensor([p[3] for p in paths], dtype=torch.long, device=DEVICE)
    h3s = torch.tensor([p[4] for p in paths], dtype=torch.long, device=DEVICE)
    r3s = torch.tensor([p[5] for p in paths], dtype=torch.long, device=DEVICE)

    with torch.no_grad():
        e_h1 = model.entity_emb(h1s)
        e_h2 = model.entity_emb(h2s)
        e_h3 = model.entity_emb(h3s)

        if is_ctx:
            # Step 1: θ_{r1}(h1, h2) · e_{h1}
            op1 = model._get_conditioned_op(r1s, e_h1, e_h2)
            v = torch.bmm(op1, e_h1.unsqueeze(-1)).squeeze(-1)
            # Step 2: θ_{r2}(h2, h3) · v
            op2 = model._get_conditioned_op(r2s, e_h2, e_h3)
            v = torch.bmm(op2, v.unsqueeze(-1)).squeeze(-1)
            # Step 3: θ_{r3}(h3, h1) · v  (closing the cycle)
            op3 = model._get_conditioned_op(r3s, e_h3, e_h1)
            v = torch.bmm(op3, v.unsqueeze(-1)).squeeze(-1)
        else:
            ops = model.get_relation_matrices()
            v = torch.bmm(ops[r1s], e_h1.unsqueeze(-1)).squeeze(-1)
            v = torch.bmm(ops[r2s], v.unsqueeze(-1)).squeeze(-1)
            v = torch.bmm(ops[r3s], v.unsqueeze(-1)).squeeze(-1)

        holonomies = (v - e_h1).norm(dim=1).cpu().numpy()

    return holonomies


def compute_gt_path_holonomy(Q_ops, entities, paths):
    """GT entity-path holonomy using ground-truth operators."""
    holonomies = np.zeros(len(paths))
    for i, (h1, r1, h2, r2, h3, r3) in enumerate(paths):
        v = Q_ops[r1] @ entities[h1]
        v = Q_ops[r2] @ v
        v = Q_ops[r3] @ v
        holonomies[i] = np.linalg.norm(v - entities[h1])
    return holonomies


def evaluate_path_holonomy(model, con_paths, incon_paths,
                           gt_con_hol, gt_incon_hol, label="Model"):
    """Full entity-path holonomy evaluation with stats."""
    hol_con = _compute_path_holonomy(model, con_paths)
    hol_incon = _compute_path_holonomy(model, incon_paths)

    labels = np.concatenate([np.zeros(len(hol_con)), np.ones(len(hol_incon))])
    scores = np.concatenate([hol_con, hol_incon])

    try:
        auc = roc_auc_score(labels, scores)
    except ValueError:
        auc = float('nan')

    gt_all = np.concatenate([gt_con_hol, gt_incon_hol])
    if np.std(gt_all) > 1e-10 and np.std(scores) > 1e-10:
        corr = np.corrcoef(gt_all, scores)[0, 1]
    else:
        corr = float('nan')

    sep = float(np.mean(hol_incon) / max(np.mean(hol_con), 1e-10))

    results = {
        'auc': auc,
        'correlation': corr,
        'con_mean': float(np.mean(hol_con)),
        'con_std': float(np.std(hol_con)),
        'incon_mean': float(np.mean(hol_incon)),
        'incon_std': float(np.std(hol_incon)),
        'separation': sep,
    }

    print(f"\n  [{label}] Entity-Path Holonomy:")
    print(f"    Consistent:   {results['con_mean']:.4f} ± {results['con_std']:.4f}")
    print(f"    Inconsistent: {results['incon_mean']:.4f} ± {results['incon_std']:.4f}")
    print(f"    Separation:   {results['separation']:.2f}×")
    print(f"    ROC-AUC:      {results['auc']:.4f}")
    print(f"    Corr w/ GT:   {results['correlation']:.4f}")

    return results


# ============================================================
#  MAIN EXPERIMENT
# ============================================================

def run_experiment(seed=42):
    # ---- Config (same GT as Exp1b) ----
    N_BASE = 12
    N_CLOSURES = 12
    N_REL = N_BASE + N_CLOSURES  # = 24
    ANGLE_SCALE = 0.25
    N_ENTITIES = 200
    N_TRAIN = 8000
    N_VAL = 2000
    EPOCHS = 800
    LR = 1e-3
    HEAD_HIDDEN = 32
    CTX_HIDDEN = 16
    LAMBDA_K = 0.5
    LAMBDA_ORTHO = 1.0
    D = 3

    print("=" * 70)
    print("  Exp2 — Context-Conditioned Operators for Holonomy Detection")
    print("=" * 70)
    print(f"  Device: {DEVICE}")
    print(f"  Seed: {seed}")
    print(f"  Config: {N_REL} relations, angle_scale={ANGLE_SCALE}")
    print(f"  Context MLP: hidden={CTX_HIDDEN}, α_init=0.01")

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

    # --- 2. Generate KG ---
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

    # --- 3. Generate relation-level cycles (for raw holonomy) ---
    print("\n--- Generating relation-level cycles ---")
    cycle_data = generate_cycles_v2(
        Q_ops, closure_cycles,
        n_consistent=200, n_inconsistent=200, seed=seed,
    )
    print(f"  Consistent: {len(cycle_data['consistent_cycles'])}, "
          f"GT hol mean: {cycle_data['consistent_holonomy'].mean():.6f}")
    print(f"  Inconsistent: {len(cycle_data['inconsistent_cycles'])}, "
          f"GT hol mean: {cycle_data['inconsistent_holonomy'].mean():.4f}")

    # --- 4. Generate entity-grounded cycles (for path holonomy) ---
    print("\n--- Generating entity-grounded cycles ---")
    con_paths, incon_paths = generate_entity_grounded_cycles(
        nearest_tail, closure_cycles,
        N_ENTITIES, N_REL,
        n_consistent=300, n_inconsistent=300, seed=seed + 500,
    )
    print(f"  Consistent paths: {len(con_paths)}")
    print(f"  Inconsistent paths: {len(incon_paths)}")

    # GT entity-path holonomy
    gt_con_hol = compute_gt_path_holonomy(Q_ops, entities, con_paths)
    gt_incon_hol = compute_gt_path_holonomy(Q_ops, entities, incon_paths)
    print(f"  GT path holonomy — consistent:   {gt_con_hol.mean():.6f} ± {gt_con_hol.std():.6f}")
    print(f"  GT path holonomy — inconsistent: {gt_incon_hol.mean():.4f} ± {gt_incon_hol.std():.4f}")
    gt_path_sep = gt_incon_hol.mean() / max(gt_con_hol.mean(), 1e-10)
    print(f"  GT path separation: {gt_path_sep:.1f}×")

    # --- 5. Train models ---
    print(f"\n--- Training AlgKG_Ctx (λ_K={LAMBDA_K}, λ_ortho={LAMBDA_ORTHO}, ctx_hidden={CTX_HIDDEN}) ---")
    algkg_ctx = AlgKGCtxModel(N_ENTITIES, N_REL, n_gen, HEAD_HIDDEN, CTX_HIDDEN).to(DEVICE)
    hist_ctx = train_model(
        algkg_ctx, train_triples, train_labels,
        val_triples, val_labels,
        epochs=EPOCHS, lr=LR, lambda_k=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
        label="AlgKG_Ctx",
    )

    print(f"\n--- Training AlgKG (λ_K={LAMBDA_K}, λ_ortho={LAMBDA_ORTHO}) ---")
    algkg = AlgKGModel(N_ENTITIES, N_REL, n_gen, HEAD_HIDDEN).to(DEVICE)
    hist_alg = train_model(
        algkg, train_triples, train_labels,
        val_triples, val_labels,
        epochs=EPOCHS, lr=LR, lambda_k=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
        label="AlgKG",
    )

    print(f"\n--- Training MLP (no constraints) ---")
    mlp = MLPKGModel(N_ENTITIES, N_REL, n_gen, HEAD_HIDDEN).to(DEVICE)
    hist_mlp = train_model(
        mlp, train_triples, train_labels,
        val_triples, val_labels,
        epochs=EPOCHS, lr=LR,
        label="MLP",
    )

    # --- 6. Raw operator holonomy (Exp1b-compatible) ---
    print("\n--- Raw Operator Holonomy (Exp1b-compatible) ---")
    raw_ctx = evaluate_raw_holonomy(algkg_ctx, cycle_data, "AlgKG_Ctx")
    raw_alg = evaluate_raw_holonomy(algkg, cycle_data, "AlgKG")
    raw_mlp = evaluate_raw_holonomy(mlp, cycle_data, "MLP")

    # --- 7. Entity-path holonomy (NEW) ---
    print("\n--- Entity-Path Holonomy (KEY METRIC) ---")
    path_ctx = evaluate_path_holonomy(
        algkg_ctx, con_paths, incon_paths, gt_con_hol, gt_incon_hol, "AlgKG_Ctx")
    path_alg = evaluate_path_holonomy(
        algkg, con_paths, incon_paths, gt_con_hol, gt_incon_hol, "AlgKG")
    path_mlp = evaluate_path_holonomy(
        mlp, con_paths, incon_paths, gt_con_hol, gt_incon_hol, "MLP")

    # --- 8. Summary ---
    wall_time = time.time() - t0
    print("\n" + "=" * 70)
    print("  SUMMARY — Exp2")
    print("=" * 70)

    print(f"\n  Triple Classification (val acc):")
    print(f"    AlgKG_Ctx: {hist_ctx['best_val_acc']:.4f}")
    print(f"    AlgKG:     {hist_alg['best_val_acc']:.4f}")
    print(f"    MLP:       {hist_mlp['best_val_acc']:.4f}")

    print(f"\n  Raw Operator Holonomy (ROC-AUC):")
    print(f"    AlgKG_Ctx: {raw_ctx['auc']:.4f}  (base ops only)")
    print(f"    AlgKG:     {raw_alg['auc']:.4f}")
    print(f"    MLP:       {raw_mlp['auc']:.4f}")

    print(f"\n  Entity-Path Holonomy (ROC-AUC):   ← KEY METRIC")
    print(f"    AlgKG_Ctx: {path_ctx['auc']:.4f}  (sep: {path_ctx['separation']:.2f}×)")
    print(f"    AlgKG:     {path_alg['auc']:.4f}  (sep: {path_alg['separation']:.2f}×)")
    print(f"    MLP:       {path_mlp['auc']:.4f}  (sep: {path_mlp['separation']:.2f}×)")

    print(f"\n  Entity-Path ↔ GT Correlation:")
    print(f"    AlgKG_Ctx: {path_ctx['correlation']:.4f}")
    print(f"    AlgKG:     {path_alg['correlation']:.4f}")
    print(f"    MLP:       {path_mlp['correlation']:.4f}")

    print(f"\n  AlgKG_Ctx Final Diagnostics:")
    print(f"    T_K:       {hist_ctx['tk'][-1]:.4f}")
    print(f"    OrthoErr:  {hist_ctx['ortho_err'][-1]:.4f}")
    print(f"    α (final): {hist_ctx['alpha'][-1]:.6f}")

    print(f"\n  AlgKG Final Diagnostics:")
    print(f"    T_K:       {hist_alg['tk'][-1]:.4f}")
    print(f"    OrthoErr:  {hist_alg['ortho_err'][-1]:.4f}")

    print(f"\n  MLP Final Diagnostics:")
    print(f"    OrthoErr:  {hist_mlp['ortho_err'][-1]:.4f}")

    print(f"\n  Wall time: {wall_time:.1f} s")

    # --- 9. Save results ---
    out_path = os.path.join(os.path.dirname(__file__), 'exp2_results.npz')
    np.savez(out_path, **{
        'ctx_val_acc': np.array(hist_ctx['val_acc']),
        'alg_val_acc': np.array(hist_alg['val_acc']),
        'mlp_val_acc': np.array(hist_mlp['val_acc']),
        'ctx_tk': np.array(hist_ctx['tk']),
        'ctx_ortho': np.array(hist_ctx['ortho_err']),
        'ctx_alpha': np.array(hist_ctx['alpha']),
        'alg_tk': np.array(hist_alg['tk']),
        'alg_ortho': np.array(hist_alg['ortho_err']),
        'raw_ctx_auc': raw_ctx['auc'],
        'raw_alg_auc': raw_alg['auc'],
        'raw_mlp_auc': raw_mlp['auc'],
        'path_ctx_auc': path_ctx['auc'],
        'path_alg_auc': path_alg['auc'],
        'path_mlp_auc': path_mlp['auc'],
        'path_ctx_corr': path_ctx['correlation'],
        'path_alg_corr': path_alg['correlation'],
        'path_mlp_corr': path_mlp['correlation'],
        'path_ctx_sep': path_ctx['separation'],
        'path_alg_sep': path_alg['separation'],
        'path_mlp_sep': path_mlp['separation'],
    })
    print(f"\n  Results saved to {out_path}")

    return {
        'hist_ctx': hist_ctx, 'hist_alg': hist_alg, 'hist_mlp': hist_mlp,
        'raw_ctx': raw_ctx, 'raw_alg': raw_alg, 'raw_mlp': raw_mlp,
        'path_ctx': path_ctx, 'path_alg': path_alg, 'path_mlp': path_mlp,
        'wall_time': wall_time,
    }


if __name__ == '__main__':
    run_experiment(seed=42)
