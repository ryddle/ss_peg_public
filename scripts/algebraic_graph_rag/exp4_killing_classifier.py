"""
Experiment 4 — Killing Spectrum as Relation Classifier
=======================================================
Can T_K separate compositional relations (that form closure cycles) from
atomic relations (random orthogonal, no algebraic structure) in a mixed KG,
**without any cycle supervision**?

Setup:
  - Type A (compositional): 8 base + 8 closure ops from SU(3) adjoint
    → participate in algebraic closure cycles
  - Type B (atomic): 8 random SO(n_gen) ops (generic, no closures)
    → valid KG triples but NO algebraic composition

Protocol:
  Phase 1: Train with L_triple + T_K + T_ortho (800 ep, NO cycle loss)
  Phase 2: Killing spectrum analysis → classify relations
  Phase 3: Selective cycle loss on compositional subset (400 ep)

Key question: Does the Killing diagonal K_rr separate Type A from Type B
             after triple-only training?

Motivated by GRAPH_TOPOLOGY_FROM_KILLING §3.1:
  "Diagonalize K → classify edges → draw the graph"
"""

import sys, os, time, copy
import torch
import torch.nn as nn
import numpy as np
from scipy.linalg import expm
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

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'training_ai_test'))
from killing_utils import killing_metric

sys.stdout.reconfigure(encoding='utf-8')


# ============================================================
#  CONFIG
# ============================================================

SEED = 42
N_BASE_A = 8          # base algebraic operators (Type A)
N_CLOSURES_A = 8      # closure operators (Type A)
N_RANDOM_B = 8        # random orthogonal operators (Type B)
N_TYPE_A = N_BASE_A + N_CLOSURES_A  # 16
N_REL = N_TYPE_A + N_RANDOM_B       # 24
ANGLE_SCALE = 0.25
N_ENTITIES = 200
N_TRAIN, N_VAL = 8000, 2000
PHASE1_EPOCHS = 800
PHASE3_EPOCHS = 400
LR = 1e-3
BATCH_SIZE = 256
HEAD_HIDDEN = 32
LAMBDA_K, LAMBDA_ORTHO = 0.5, 1.0
LAMBDA_CYCLE = 1.0
MARGIN = 1.0
D = 3
GRAD_CLIP = 1.0


# ============================================================
#  RANDOM ORTHOGONAL OPERATOR GENERATION (Type B)
# ============================================================

def generate_random_orthogonal_operators(n_ops, n_gen, angle_scale=0.25, seed=9999):
    """
    Generate random SO(n_gen) operators NOT from any specific Lie algebra
    representation.  Uses random antisymmetric n_gen×n_gen matrices,
    exponentiated at the same angular scale as the algebraic operators.

    These are generic elements of SO(n_gen) — they live in the full so(n_gen)
    (dim = n_gen*(n_gen-1)/2 = 28 for n_gen=8), NOT in the 8-dimensional
    adjoint image of su(3).
    """
    rng = np.random.default_rng(seed)
    ops = np.zeros((n_ops, n_gen, n_gen), dtype=np.float64)
    for r in range(n_ops):
        M = rng.standard_normal((n_gen, n_gen))
        A = (M - M.T) / 2  # antisymmetric
        # Normalize so ||A||_F matches typical su(3) adjoint generators
        A = A / np.linalg.norm(A, 'fro') * np.sqrt(n_gen)
        ops[r] = expm(angle_scale * A)
    return ops


# ============================================================
#  MODEL (same as Exp3)
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


# ============================================================
#  CYCLE CONTRASTIVE LOSS (same as Exp3)
# ============================================================

def compute_cycle_loss(ops, con_cycles_t, incon_cycles_t, margin=1.0):
    ng = ops.shape[1]
    I_ng = torch.eye(ng, dtype=ops.dtype, device=ops.device).unsqueeze(0)
    H_con = ops[con_cycles_t[:, 2]] @ ops[con_cycles_t[:, 1]] @ ops[con_cycles_t[:, 0]]
    loss_con = ((H_con - I_ng) ** 2).sum(dim=(1, 2)).mean()
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

        if ep % log_every == 0 or ep == 1:
            sep = hol_incon_mean / max(hol_con_mean, 1e-10) if not np.isnan(hol_con_mean) else 0
            extra = f"  hol_sep={sep:.2f}x"
            if lambda_cycle > 0:
                extra = f"  cyc_L={epoch_cyc / max(n_batches, 1):.4f}" + extra
            print(f"  [{label}] ep {ep:4d}  train={train_acc:.4f}  val={val_acc:.4f}  "
                  f"loss={epoch_loss / n:.4f}  T_K={tk:.4f}  OE={ortho_err:.4f}{extra}")

    return history


# ============================================================
#  KILLING SPECTRUM ANALYSIS
# ============================================================

def analyze_killing_spectrum(model, type_labels, n_type_a, n_random_b):
    """
    Compute and analyze the Killing matrix of learned operators.

    type_labels: array of 0 (Type A) / 1 (Type B) for each relation
    Returns dict with all diagnostic quantities.
    """
    model.eval()
    with torch.no_grad():
        ops = model.get_relation_matrices()  # (n_rel, n_gen, n_gen)

        # --- Full Killing matrix ---
        K = killing_metric(ops)  # (n_rel, n_rel) on DEVICE
        K_np = K.cpu().numpy()

        # --- Eigendecomposition ---
        eigenvalues, eigenvectors = np.linalg.eigh(K_np)  # ascending order
        # Sort descending by absolute value (most significant first)
        order = np.argsort(-np.abs(eigenvalues))
        eigenvalues = eigenvalues[order]
        eigenvectors = eigenvectors[:, order]

        # --- Per-relation Killing diagonal ---
        K_diag = np.diag(K_np)

        # --- Commutator matrix ||[θ_r, θ_r']||_F ---
        ops_np = ops.cpu().numpy()
        n_rel = ops_np.shape[0]
        comm_norms = np.zeros((n_rel, n_rel))
        for a in range(n_rel):
            for b in range(n_rel):
                comm = ops_np[a] @ ops_np[b] - ops_np[b] @ ops_np[a]
                comm_norms[a, b] = np.linalg.norm(comm, 'fro')

        # --- Block analysis: A-A, A-B, B-B mean commutator norms ---
        mask_A = type_labels == 0
        mask_B = type_labels == 1
        AA_mask = np.outer(mask_A, mask_A)
        AB_mask = np.outer(mask_A, mask_B) | np.outer(mask_B, mask_A)
        BB_mask = np.outer(mask_B, mask_B)

        comm_AA = comm_norms[AA_mask].mean() if AA_mask.any() else 0.0
        comm_AB = comm_norms[AB_mask].mean() if AB_mask.any() else 0.0
        comm_BB = comm_norms[BB_mask].mean() if BB_mask.any() else 0.0

        # --- Orthogonality per relation ---
        I_n = np.eye(ops_np.shape[1])
        ortho_err_per_rel = np.array([
            np.linalg.norm(ops_np[r] @ ops_np[r].T - I_n, 'fro')
            for r in range(n_rel)
        ])

        # --- Distance from identity per relation ---
        dist_from_I = np.array([
            np.linalg.norm(ops_np[r] - I_n, 'fro')
            for r in range(n_rel)
        ])

        # --- Classification by K_diag ---
        # Type A should have K_rr < 0 (compact/cyclic); Type B different
        # Try: classify by K_diag threshold at median
        K_diag_A = K_diag[mask_A]
        K_diag_B = K_diag[mask_B]

        # ROC-AUC: can K_diag separate Type A from Type B?
        # Convention: higher K_diag = more likely Type B (less integrated)
        try:
            k_auc = roc_auc_score(type_labels, K_diag)
        except ValueError:
            k_auc = float('nan')

        # Also try with -|K_diag| (lower magnitude = less algebraically integrated)
        try:
            k_auc_abs = roc_auc_score(type_labels, -np.abs(K_diag))
        except ValueError:
            k_auc_abs = float('nan')

        # Try commutator total as classifier: Σ_b ||[θ_a, θ_b]||_F
        comm_total = comm_norms.sum(axis=1)
        try:
            comm_auc = roc_auc_score(type_labels, -comm_total)
        except ValueError:
            comm_auc = float('nan')

    return {
        'K_matrix': K_np,
        'eigenvalues': eigenvalues,
        'eigenvectors': eigenvectors,
        'K_diag': K_diag,
        'K_diag_A': K_diag_A,
        'K_diag_B': K_diag_B,
        'comm_norms': comm_norms,
        'comm_AA': comm_AA, 'comm_AB': comm_AB, 'comm_BB': comm_BB,
        'comm_total': comm_total,
        'ortho_err_per_rel': ortho_err_per_rel,
        'dist_from_I': dist_from_I,
        'k_diag_auc': k_auc,
        'k_diag_abs_auc': k_auc_abs,
        'comm_total_auc': comm_auc,
    }


# ============================================================
#  HOLONOMY EVALUATION (same as Exp3)
# ============================================================

def _compute_path_holonomy(model, paths):
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
    model.eval()
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

    print(f"  [{label:16s}]  raw_AUC={raw_auc:.4f} (sep={raw_sep:.2f}x)  "
          f"path_AUC={path_auc:.4f} (sep={path_sep:.2f}x)  corr={corr:.4f}")

    return {
        'raw_auc': raw_auc, 'raw_sep': raw_sep,
        'raw_con': float(np.mean(hol_con)), 'raw_incon': float(np.mean(hol_incon)),
        'path_auc': path_auc, 'path_sep': path_sep,
        'path_con': float(np.mean(path_hol_con)), 'path_incon': float(np.mean(path_hol_incon)),
        'correlation': corr,
    }


# ============================================================
#  MAIN EXPERIMENT
# ============================================================

def run_experiment(seed=SEED):
    print("=" * 70)
    print("  Exp4 — Killing Spectrum as Relation Classifier")
    print("=" * 70)
    print(f"  Device: {DEVICE}")
    print(f"  Seed: {seed}")
    print(f"  Type A: {N_BASE_A} base + {N_CLOSURES_A} closure = {N_TYPE_A} compositional")
    print(f"  Type B: {N_RANDOM_B} random orthogonal = atomic")
    print(f"  Total relations: {N_REL}")
    print(f"  Phase 1: {PHASE1_EPOCHS} ep (L_triple + T_K + T_ortho, NO cycle loss)")
    print(f"  Phase 3: {PHASE3_EPOCHS} ep (L_triple + L_cycle on Type A only)")

    t0 = time.time()

    # ========================================================
    #  PHASE 0: DATA GENERATION
    # ========================================================

    print("\n--- Generating Type A (compositional) operators with closures ---")
    Q_A, alphas, closure_cycles = generate_relation_operators_with_closures(
        n_base=N_BASE_A, n_closures=N_CLOSURES_A, d=D,
        angle_scale=ANGLE_SCALE, seed=seed,
    )
    n_gen = Q_A.shape[1]
    print(f"  Type A: {Q_A.shape[0]} operators in R^{n_gen}x{n_gen}")
    print(f"  Closure cycles: {len(closure_cycles)}")

    I_np = np.eye(n_gen)
    closure_hols = [np.linalg.norm(Q_A[rc] @ Q_A[rj] @ Q_A[ri] - I_np, 'fro')
                    for ri, rj, rc in closure_cycles]
    print(f"  Closure holonomy check: max={max(closure_hols):.2e}, mean={np.mean(closure_hols):.2e}")

    print("\n--- Generating Type B (atomic) random orthogonal operators ---")
    Q_B = generate_random_orthogonal_operators(
        N_RANDOM_B, n_gen, angle_scale=ANGLE_SCALE, seed=seed + 5000,
    )

    # Verify distance from identity is comparable
    dist_A = np.mean([np.linalg.norm(Q_A[r] - I_np, 'fro') for r in range(Q_A.shape[0])])
    dist_B = np.mean([np.linalg.norm(Q_B[r] - I_np, 'fro') for r in range(Q_B.shape[0])])
    print(f"  Type B: {Q_B.shape[0]} operators in R^{n_gen}x{n_gen}")
    print(f"  Distance from I: Type A mean={dist_A:.4f}, Type B mean={dist_B:.4f}")

    # Check maximum accidental closure among Type B triples
    max_b_closure = 0.0
    for b1 in range(N_RANDOM_B):
        for b2 in range(N_RANDOM_B):
            for b3 in range(N_RANDOM_B):
                # Use absolute indices (Type B starts at N_TYPE_A)
                H = Q_B[b3] @ Q_B[b2] @ Q_B[b1]
                hol = np.linalg.norm(H - I_np, 'fro')
                if hol < max_b_closure or max_b_closure == 0:
                    max_b_closure = hol
    print(f"  Best Type B-only closure: {max_b_closure:.4f} (should be >> 0)")

    # Concatenate all operators
    Q_all = np.concatenate([Q_A, Q_B], axis=0)
    assert Q_all.shape[0] == N_REL

    # Type labels: 0 = Type A (compositional), 1 = Type B (atomic)
    type_labels = np.array([0] * N_TYPE_A + [1] * N_RANDOM_B)

    # --- GT Killing analysis on ground-truth operators ---
    print("\n--- Ground truth Killing analysis ---")
    Q_all_t = torch.tensor(Q_all, dtype=torch.float64, device=DEVICE)
    with torch.no_grad():
        K_gt = killing_metric(Q_all_t).cpu().numpy()
    K_gt_diag = np.diag(K_gt)
    print(f"  GT K_rr  Type A: mean={K_gt_diag[:N_TYPE_A].mean():.4f}, "
          f"std={K_gt_diag[:N_TYPE_A].std():.4f}")
    print(f"  GT K_rr  Type B: mean={K_gt_diag[N_TYPE_A:].mean():.4f}, "
          f"std={K_gt_diag[N_TYPE_A:].std():.4f}")
    try:
        gt_k_auc = roc_auc_score(type_labels, K_gt_diag)
    except ValueError:
        gt_k_auc = float('nan')
    print(f"  GT K_rr  AUC for A/B separation: {gt_k_auc:.4f}")

    # --- Generate synthetic KG ---
    print("\n--- Generating synthetic KG (mixed Type A + B) ---")
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

    dist_A_nn = nearest_dist[:N_TYPE_A].mean()
    dist_B_nn = nearest_dist[N_TYPE_A:].mean()
    print(f"  Entities: {N_ENTITIES}")
    print(f"  NN dist: Type A mean={dist_A_nn:.4f}, Type B mean={dist_B_nn:.4f}")

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
    print(f"  Train: {len(train_triples)} ({train_labels.sum()} valid)")
    print(f"  Val:   {len(val_triples)} ({val_labels.sum()} valid)")

    # --- Cycles (only involving Type A closure relations) ---
    print("\n--- Generating cycles (Type A closures only) ---")
    cycle_data = generate_cycles_v2(
        Q_all, closure_cycles,
        n_consistent=200, n_inconsistent=200, seed=seed,
    )
    print(f"  Consistent: {len(cycle_data['consistent_cycles'])}, "
          f"GT hol={cycle_data['consistent_holonomy'].mean():.6f}")
    print(f"  Inconsistent: {len(cycle_data['inconsistent_cycles'])}, "
          f"GT hol={cycle_data['inconsistent_holonomy'].mean():.4f}")

    con_cycles_t = torch.tensor(cycle_data['consistent_cycles'], dtype=torch.long, device=DEVICE)
    incon_cycles_t = torch.tensor(cycle_data['inconsistent_cycles'], dtype=torch.long, device=DEVICE)

    print("\n--- Generating entity-grounded cycles ---")
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
    gt_path_sep = gt_incon_hol.mean() / max(gt_con_hol.mean(), 1e-10)
    print(f"  Consistent paths: {len(con_paths)}, GT hol={gt_con_hol.mean():.6f}")
    print(f"  Inconsistent paths: {len(incon_paths)}, GT hol={gt_incon_hol.mean():.4f}")
    print(f"  GT path separation: {gt_path_sep:.1f}x")

    # ========================================================
    #  PHASE 1: TRIPLE-ONLY TRAINING (NO CYCLE LOSS)
    # ========================================================

    print(f"\n{'='*70}")
    print(f"  PHASE 1: Triple-only Training (L_triple + T_K + T_ortho)")
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
        grad_clip=GRAD_CLIP, label="Phase1", log_every=100,
    )

    print(f"\n  Phase 1 holonomy evaluation (without cycle loss):")
    hol_p1 = evaluate_holonomy(model, cycle_data, con_paths, incon_paths,
                                gt_con_hol, gt_incon_hol, label="Phase1")

    # ========================================================
    #  PHASE 2: KILLING SPECTRUM ANALYSIS
    # ========================================================

    print(f"\n{'='*70}")
    print(f"  PHASE 2: Killing Spectrum Analysis")
    print(f"{'='*70}")

    analysis = analyze_killing_spectrum(model, type_labels, N_TYPE_A, N_RANDOM_B)

    print(f"\n  Eigenvalue spectrum of learned Killing matrix:")
    eigs = analysis['eigenvalues']
    n_neg = (eigs < -1e-6).sum()
    n_zero = (np.abs(eigs) < 1e-6).sum()
    n_pos = (eigs > 1e-6).sum()
    print(f"    Negative: {n_neg}, Zero: {n_zero}, Positive: {n_pos}")
    print(f"    Top 5: {eigs[:5]}")
    print(f"    Bottom 5: {eigs[-5:]}")

    print(f"\n  K_rr (Killing diagonal) per relation type:")
    print(f"    Type A (compositional):")
    print(f"      mean={analysis['K_diag_A'].mean():.4f}, "
          f"std={analysis['K_diag_A'].std():.4f}, "
          f"range=[{analysis['K_diag_A'].min():.4f}, {analysis['K_diag_A'].max():.4f}]")
    print(f"    Type B (atomic):")
    print(f"      mean={analysis['K_diag_B'].mean():.4f}, "
          f"std={analysis['K_diag_B'].std():.4f}, "
          f"range=[{analysis['K_diag_B'].min():.4f}, {analysis['K_diag_B'].max():.4f}]")

    print(f"\n  Classification AUCs:")
    print(f"    K_rr as classifier:       ROC-AUC = {analysis['k_diag_auc']:.4f}")
    print(f"    |K_rr| as classifier:     ROC-AUC = {analysis['k_diag_abs_auc']:.4f}")
    print(f"    Comm total as classifier:  ROC-AUC = {analysis['comm_total_auc']:.4f}")

    print(f"\n  Commutator block structure (mean ||[θ_a, θ_b]||_F):")
    print(f"    A↔A: {analysis['comm_AA']:.6f}")
    print(f"    A↔B: {analysis['comm_AB']:.6f}")
    print(f"    B↔B: {analysis['comm_BB']:.6f}")

    print(f"\n  Per-relation detail:")
    print(f"    {'Rel':>4s}  {'Type':>6s}  {'K_rr':>10s}  {'CommTotal':>10s}  "
          f"{'OrthoErr':>10s}  {'||θ-I||':>10s}")
    print(f"    {'----':>4s}  {'------':>6s}  {'----------':>10s}  {'----------':>10s}  "
          f"{'----------':>10s}  {'----------':>10s}")
    for r in range(N_REL):
        t = "A_base" if r < N_BASE_A else ("A_clos" if r < N_TYPE_A else "B_rand")
        print(f"    {r:4d}  {t:>6s}  {analysis['K_diag'][r]:10.4f}  "
              f"{analysis['comm_total'][r]:10.4f}  "
              f"{analysis['ortho_err_per_rel'][r]:10.6f}  "
              f"{analysis['dist_from_I'][r]:10.6f}")

    # Save state for Phase 3
    state_after_p1 = copy.deepcopy(model.state_dict())
    adam_after_p1 = copy.deepcopy(optimizer.state_dict())

    # ========================================================
    #  PHASE 3: SELECTIVE CYCLE LOSS (on Type A only)
    # ========================================================

    print(f"\n{'='*70}")
    print(f"  PHASE 3: Selective Cycle Loss (400 ep, on Type A closures)")
    print(f"{'='*70}")

    # Restore to end-of-Phase1 state
    model.load_state_dict(copy.deepcopy(state_after_p1))
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    optimizer.load_state_dict(adam_after_p1)

    hist_p3 = train_phase(
        model, optimizer,
        train_triples, train_labels, val_triples, val_labels,
        epochs=PHASE3_EPOCHS, batch_size=BATCH_SIZE,
        lambda_k=LAMBDA_K, lambda_ortho=LAMBDA_ORTHO,
        lambda_cycle=LAMBDA_CYCLE,
        con_cycles_t=con_cycles_t, incon_cycles_t=incon_cycles_t,
        margin=MARGIN,
        grad_clip=GRAD_CLIP, label="Phase3:CycLoss", log_every=100,
    )

    print(f"\n  Phase 3 holonomy evaluation:")
    hol_p3 = evaluate_holonomy(model, cycle_data, con_paths, incon_paths,
                                gt_con_hol, gt_incon_hol, label="Phase3")

    # --- Phase 3 Killing analysis (how did cycle loss change the spectrum?) ---
    print(f"\n  Phase 3 Killing spectrum (after cycle loss):")
    analysis_p3 = analyze_killing_spectrum(model, type_labels, N_TYPE_A, N_RANDOM_B)
    print(f"    K_rr  Type A: mean={analysis_p3['K_diag_A'].mean():.4f}, "
          f"std={analysis_p3['K_diag_A'].std():.4f}")
    print(f"    K_rr  Type B: mean={analysis_p3['K_diag_B'].mean():.4f}, "
          f"std={analysis_p3['K_diag_B'].std():.4f}")
    print(f"    K_rr AUC: {analysis_p3['k_diag_auc']:.4f}")
    print(f"    Comm A↔A={analysis_p3['comm_AA']:.6f}  "
          f"A↔B={analysis_p3['comm_AB']:.6f}  "
          f"B↔B={analysis_p3['comm_BB']:.6f}")

    # ========================================================
    #  SUMMARY
    # ========================================================

    t_total = time.time() - t0

    print(f"\n{'='*70}")
    print(f"  SUMMARY — Exp4: Killing Spectrum as Relation Classifier")
    print(f"{'='*70}")

    print(f"\n  Ground Truth Operators:")
    print(f"    K_rr separation (GT):  Type A mean={K_gt_diag[:N_TYPE_A].mean():.4f}, "
          f"Type B mean={K_gt_diag[N_TYPE_A:].mean():.4f}")
    print(f"    GT K_rr AUC: {gt_k_auc:.4f}")

    print(f"\n  Phase 1 (triple-only, 800ep):")
    print(f"    Val acc: {hist_p1['val_acc'][-1]:.4f}")
    print(f"    T_K: {hist_p1['tk'][-1]:.4f}")
    print(f"    Holonomy (no cycle loss): raw_AUC={hol_p1['raw_auc']:.4f}, "
          f"path_AUC={hol_p1['path_auc']:.4f}")

    print(f"\n  Killing Classification (Phase 1 — the KEY result):")
    print(f"    K_rr ROC-AUC:        {analysis['k_diag_auc']:.4f}")
    print(f"    |K_rr| ROC-AUC:      {analysis['k_diag_abs_auc']:.4f}")
    print(f"    CommTotal ROC-AUC:    {analysis['comm_total_auc']:.4f}")
    print(f"    K_rr  Type A: {analysis['K_diag_A'].mean():.4f} ± {analysis['K_diag_A'].std():.4f}")
    print(f"    K_rr  Type B: {analysis['K_diag_B'].mean():.4f} ± {analysis['K_diag_B'].std():.4f}")
    print(f"    Comm blocks: A↔A={analysis['comm_AA']:.6f}, "
          f"A↔B={analysis['comm_AB']:.6f}, "
          f"B↔B={analysis['comm_BB']:.6f}")

    print(f"\n  Phase 3 (selective cycle loss, 400ep):")
    print(f"    Val acc: {hist_p3['val_acc'][-1]:.4f}")
    print(f"    Holonomy: raw_AUC={hol_p3['raw_auc']:.4f}, "
          f"path_AUC={hol_p3['path_auc']:.4f}, corr={hol_p3['correlation']:.4f}")
    print(f"    K_rr separation: Type A={analysis_p3['K_diag_A'].mean():.4f}, "
          f"Type B={analysis_p3['K_diag_B'].mean():.4f}")

    print(f"\n  Wall time: {t_total:.1f} s")

    # --- Save results ---
    save_path = os.path.join(os.path.dirname(__file__), 'exp4_results.npz')
    np.savez(save_path,
        type_labels=type_labels,
        # Phase 1
        p1_K_matrix=analysis['K_matrix'],
        p1_eigenvalues=analysis['eigenvalues'],
        p1_K_diag=analysis['K_diag'],
        p1_comm_norms=analysis['comm_norms'],
        p1_hol_raw_auc=hol_p1['raw_auc'],
        p1_hol_path_auc=hol_p1['path_auc'],
        # Phase 3
        p3_K_diag_A_mean=analysis_p3['K_diag_A'].mean(),
        p3_K_diag_B_mean=analysis_p3['K_diag_B'].mean(),
        p3_hol_raw_auc=hol_p3['raw_auc'],
        p3_hol_path_auc=hol_p3['path_auc'],
        p3_hol_corr=hol_p3['correlation'],
        # GT
        gt_K_diag=K_gt_diag,
        gt_k_auc=gt_k_auc,
    )
    print(f"  Results saved to {save_path}")


if __name__ == "__main__":
    run_experiment()
