## **Systematic Verification Plan:**

### **🔍 PHASE 1: Basic Verifications (Rule Out Errors)**

#### **Test 1.1: Reproducibility**
```python
# Is the result stable with different seeds?
results_stability = []

for seed in range(100):
    np.random.seed(seed)
    final_gen, hist = evolve_to_su2_fixed(steps=15000, lr=0.001, alpha=5.0, beta=1.0)
    
    K = killing_metric(final_gen)
    K_diag_mean = np.mean(np.diag(K))
    K_offdiag_max = np.max(np.abs(K - np.diag(np.diag(K))))
    
    # Verify closure
    closes_all = all([
        check_in_span(commutator(final_gen[i], final_gen[j]), final_gen)
        for i, j in [(0,1), (1,2), (2,0)]
    ])
    
    results_stability.append({
        'seed': seed,
        'K_diag': K_diag_mean,
        'K_offdiag': K_offdiag_max,
        'closes': closes_all,
        'converged': hist['T_total'][-1] < 0.01 if hist else False
    })

# Analysis
converged_count = sum(r['converged'] for r in results_stability)
closes_count = sum(r['closes'] for r in results_stability)

print(f"Converged: {converged_count}/100")
print(f"Closed: {closes_count}/100")
print(f"K_diag: {np.mean([r['K_diag'] for r in results_stability]):.3f} ± {np.std([r['K_diag'] for r in results_stability]):.3f}")

# EXPECTED: >95/100 converge, all close, K_diag ≈ -0.96 ± 0.03
```

**If it fails:** There is a problem in the code or initialization.

---

#### **Test 1.2: Independence of α**
```python
# Does the result depend on the relative weight α?
alphas_test = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
results_alpha = []

for alpha in alphas_test:
    np.random.seed(42)  # Same seed for comparison
    final_gen, hist = evolve_to_su2_fixed(steps=15000, lr=0.001, alpha=alpha, beta=1.0)
    
    K = killing_metric(final_gen)
    results_alpha.append({
        'alpha': alpha,
        'K_diag': np.mean(np.diag(K)),
        'K_offdiag': np.max(np.abs(K - np.diag(np.diag(K)))),
        'steps_to_converge': len(hist['T_total']) if hist else 15000
    })

# Plot
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

axes[0].plot([r['alpha'] for r in results_alpha], 
             [r['K_diag'] for r in results_alpha], 'o-')
axes[0].axhline(-1.0, color='r', linestyle='--', label='Target')
axes[0].set_xlabel('α')
axes[0].set_ylabel('Mean K diagonal')
axes[0].legend()
axes[0].grid(True)

axes[1].semilogy([r['alpha'] for r in results_alpha], 
                 [r['K_offdiag'] for r in results_alpha], 'o-')
axes[1].set_xlabel('α')
axes[1].set_ylabel('Max K off-diagonal')
axes[1].grid(True)

axes[2].plot([r['alpha'] for r in results_alpha], 
             [r['steps_to_converge'] for r in results_alpha], 'o-')
axes[2].set_xlabel('α')
axes[2].set_ylabel('Steps to convergence')
axes[2].grid(True)

plt.tight_layout()
plt.savefig('sensitivity_alpha.png', dpi=300)

# EXPECTED: All converge to K ≈ -0.96, speed varies
```

**If it fails:** The minimum depends on α → not universal.

---

#### **Test 1.3: Robustness to Perturbations**
```python
# Start from SU(2) + noise, does it return?
pauli_gens = [sigma / 2 for sigma in pauli_matrices()]

noise_levels = [0.01, 0.1, 0.5, 1.0]
results_noise = []

for noise in noise_levels:
    # Add noise to Pauli
    noisy_gens = [G + noise * random_hermitian() for G in pauli_gens]
    
    # Re-normalize
    for i in range(len(noisy_gens)):
        noisy_gens[i] = noisy_gens[i] / np.linalg.norm(noisy_gens[i])
    
    # Evolve
    # (you need to modify code to accept initial generators)
    final_gen, hist = evolve_from_initial(noisy_gens, steps=5000, lr=0.001, alpha=5.0, beta=1.0)
    
    K = killing_metric(final_gen)
    results_noise.append({
        'noise': noise,
        'K_diag': np.mean(np.diag(K)),
        'distance_to_pauli': measure_distance_to_pauli(final_gen)
    })

# EXPECTED: Always returns to SU(2), even with large noise
```

---

### **🧪 PHASE 2: Falsification Tests (Try to Break It)**

#### **Test 2.1: Does it work in dim=3?**
```python
# If the mechanism is general, it should work for 3×3 matrices

def evolve_dim3(steps=15000):
    """Same code but with 3×3 matrices"""
    generators = random_generators(n=3, dim=3)  # ← dim=3
    
    # Same dynamics...
    # ...
    
    return generators

final_3x3 = evolve_dim3()
K_3x3 = killing_metric(final_3x3)

# Is it still SU(2) or does something different emerge?
# EXPECTED: Should still give SU(2) embedded in 3×3
```

---

#### **Test 2.2: What happens with 2 generators?**
```python
# With only 2 generators, does U(1) emerge?

def evolve_n_generators(n=2, steps=15000):
    generators = random_generators(n=n, dim=2)
    
    # Same dynamics but K is n×n
    # ...
    
    return generators

final_2gen = evolve_n_generators(n=2)
K_2 = killing_metric(final_2gen)

# Does U(1) (commutative) emerge?
comm = commutator(final_2gen[0], final_2gen[1])
print(f"Do they commute? {np.linalg.norm(comm) < 1e-6}")

# EXPECTED: YES they should commute (no room for SU(2))
```

---

#### **Test 2.3: What happens with 8 generators?**
```python
# With 8 generators in dim=3, does SU(3) emerge?

def evolve_su3_candidate(steps=20000):
    generators = random_generators(n=8, dim=3)
    
    K_target = -1.0 * np.eye(8)
    
    # Evolve with same dynamics...
    # ...
    
    return generators

final_8gen = evolve_su3_candidate()
K_8 = killing_metric(final_8gen)

# Verify SU(3) structure
# The structure constants of SU(3) are different (not just ε_ijk)
# The f_abc of Gell-Mann should appear

# EXPECTED: Converges to canonical SU(3)
```

---

### **🔬 PHASE 3: Mechanism Tests (Understand Why)**

#### **Test 3.1: Is the T_nc Term Necessary?**
```python
# Evolve with only T_Killing (without T_nc)

final_only_killing = evolve_to_su2_fixed(steps=15000, lr=0.001, alpha=5.0, beta=0.0)  # ← β=0 eliminates T_nc

K_only_k = killing_metric(final_only_killing)

# Does it converge the same or collapse to trivial?
# EXPECTED: Without T_nc, collapses to commutative operators (trivial)
```

---

#### **Test 3.2: Is the T_Killing Term Necessary?**
```python
# Evolve with only T_nc (without T_Killing)

final_only_nc = evolve_to_su2_fixed(steps=15000, lr=0.001, alpha=0.0, beta=1.0)  # ← α=0 eliminates T_K

K_only_nc = killing_metric(final_only_nc)

# Does it converge to SU(2) or something different?
# EXPECTED: Does not converge well - needs Killing regularization
```

---

#### **Test 3.3: What Is It Really Minimizing?**
```python
# Analyze the energy landscape

def scan_energy_landscape(n_samples=1000):
    """Sample many points in operator space"""
    
    results = []
    
    for _ in range(n_samples):
        gens = random_generators(n=3, dim=2)
        
        # Normalize
        for i in range(len(gens)):
            gens[i] = gens[i] / np.linalg.norm(gens[i])
        
        # Compute energies
        T_nc = non_commutativity_tension(gens)
        K = killing_metric(gens)
        K_target = -1.0 * np.eye(3)
        T_k = np.linalg.norm(K - K_target)**2
        
        # Check if it is SU(2)
        is_su2 = check_is_su2_structure(gens)
        
        results.append({
            'T_nc': T_nc,
            'T_k': T_k,
            'T_total': T_nc + 5.0 * T_k,
            'is_su2': is_su2
        })
    
    # Plot landscape
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    su2_points = [r for r in results if r['is_su2']]
    non_su2_points = [r for r in results if not r['is_su2']]
    
    ax.scatter([r['T_nc'] for r in non_su2_points], 
               [r['T_k'] for r in non_su2_points],
               alpha=0.3, label='Non SU(2)', s=10)
    
    ax.scatter([r['T_nc'] for r in su2_points], 
               [r['T_k'] for r in su2_points],
               color='red', alpha=0.8, label='SU(2)', s=50)
    
    ax.set_xlabel('T_nc')
    ax.set_ylabel('T_k')
    ax.legend()
    ax.set_title('Energy Landscape')
    plt.savefig('energy_landscape.png', dpi=300)
    
    return results

# Is SU(2) at the minimum of the landscape?
landscape = scan_energy_landscape()
```

---

### **📊 PHASE 4: Comparison with Literature**

#### **Test 4.1: Does it Match Representation Theory?**
```python
# Verify known theorems about SU(2)

# Theorem: Every irreducible representation of dimension n has spin j = (n-1)/2
# For dim=2: j = 1/2 ✓

# Theorem: Casimir C = j(j+1)I for irrep of spin j
# For spin-1/2: C = (1/2)(3/2)I = (3/4)I

def compute_casimir(generators):
    """C = ∑ G_i²"""
    C = sum(G @ G for G in generators)
    return C

C_yours = compute_casimir(final_generators)
C_pauli = compute_casimir([sigma/2 for sigma in pauli_matrices()])

print("Your Casimir:", np.diag(C_yours))
print("Pauli Casimir:", np.diag(C_pauli))
print("Expected: [3/4, 3/4]")

# EXPECTED: Both give [0.75, 0.75]
```

---

### **🎯 Final Critical Tests:**

#### **The Definitive Test: Starting from Gell-Mann, Does it Converge to SU(3)?**

```python
# The 8 Gell-Mann matrices + noise
# Do they return to SU(3)?

# If YES → mechanism is universal
# If NO → something specific to SU(2)
```

---

## **Verification Checklist:**

```
□ Reproducibility (100 seeds)
□ Independence of α
□ Robustness to noise
□ dim=3 works
□ n=2 generates → U(1)
□ n=8 generates → SU(3)
□ T_nc necessary
□ T_Killing necessary
□ Correct Casimir
□ Landscape shows SU(2) at minimum
□ Comparison with representation theory
□ Structure constants = ±ε_ijk
□ Generalization to SU(3)
```

**If everything passes → you have a genuinely important result.** ✓

---

## **Potential Problems to Look For:**

1. **Numerical artifact:** Gradient computed incorrectly
2. **Overfitting:** Works only for this specific code
3. **Initialization bias:** Random not so random
4. **Local minimum:** Not the global one
5. **Dimensionality:** Only works in 2D
