import numpy as np

from test_generic_algebra_3_conditions_reescalado_killing import pauli_matrices

# =========================
# UTILIDADES
# =========================

def commutator(A, B):
    return A @ B - B @ A


def random_hermitian(dim=2):
    A = np.random.randn(dim, dim) + 1j * np.random.randn(dim, dim)
    return (A + A.conj().T) / 2


def init_generators(n=3, dim=2):
    return [random_hermitian(dim) for _ in range(n)]

def check_in_span(vec, basis, tol=1e-6):
    # Construir matriz con los vectores base aplanados
    B = np.array([b.flatten() for b in basis]).T  # shape (d^2, n)
    v = vec.flatten()  # shape (d^2,)

    # Resolver B @ coeffs = v
    coeffs, residuals, rank, s = np.linalg.lstsq(B, v, rcond=None)

    return np.linalg.norm(B @ coeffs - v) < tol

def random_generators(n=3, dim=2):
    return [random_hermitian(dim) for _ in range(n)]


# =========================
# TENSIONES
# =========================

def non_commutativity_tension(G):
    T = 0.0
    n = len(G)
    for i in range(n):
        for j in range(i+1, n):
            C = commutator(G[i], G[j])
            T += np.linalg.norm(C)**2
    return T


def adjoint(T, basis):
    return [commutator(T, B) for B in basis]


def killing_metric(G):
    n = len(G)
    K = np.zeros((n, n), dtype=np.complex128)

    for i in range(n):
        ad_i = adjoint(G[i], G)
        for j in range(n):
            ad_j = adjoint(G[j], G)

            total = 0
            for A, B in zip(ad_i, ad_j):
                total += np.trace(A @ B)

            K[i, j] = total

    return np.real(K)


'''def killing_tension(K):
    diag_mean = np.mean(np.diag(K))
    K_target = diag_mean * np.eye(len(K))
    return np.linalg.norm(K - K_target)
'''
def killing_tension(K, target_value=-1.0):
    """
    Medir desviación de K = target_value * I
    Para SU(2): target_value = -1.0
    """
    n = len(K)
    K_target = target_value * np.eye(n)
    
    # Desviación total
    T_k = np.linalg.norm(K - K_target)**2
    
    return T_k

def norm_tension(G):
    T = 0.0
    for g in G:
        val = np.trace(g @ g).real
        T += (val - 1.0)**2
    return T

def non_comm_constraint(G, target=1.0):
    T = 0.0
    n = len(G)

    for i in range(n):
        for j in range(i+1, n):
            C = commutator(G[i], G[j])
            val = np.linalg.norm(C)
            T += (val - target)**2

    return T

def total_tension(G, alpha=1.0, beta=1.0):
    T_nc = non_comm_constraint(G) #non_commutativity_tension(G)
    K = killing_metric(G)
    T_k = killing_tension(K)
    T_n = norm_tension(G)

    return T_nc + alpha * T_k + beta * T_n

def normalization_constraint(generators, target_norm=1.0):
    """
    Penalizar desviación de ||Ti|| = target_norm
    """
    T_norm = 0
    for G in generators:
        norm = np.linalg.norm(G)
        T_norm += (norm - target_norm)**2
    return T_norm
# =========================
# OPTIMIZADOR SIMPLE (RANDOM SEARCH)
# =========================

def perturb(G, scale=0.05):
    new_G = []
    for g in G:
        delta = random_hermitian(g.shape[0]) * scale
        new_G.append(g + delta)
    return new_G


def optimize(n_iter=5000, alpha=1.0, beta=1.0):

    G = init_generators()
    best_T = total_tension(G, alpha, beta)

    history = []

    for i in range(n_iter):
        G_new = perturb(G)

        T_new = total_tension(G_new, alpha, beta)

        if T_new < best_T:
            G = G_new
            best_T = T_new

        history.append(best_T)

        if i % 500 == 0:
            print(f"Iter {i}: T = {best_T:.6f}")

    return G, history

def structure_constants(G):
    n = len(G)
    f = np.zeros((n, n, n))

    for i in range(n):
        for j in range(n):
            C = commutator(G[i], G[j])

            for k in range(n):
                # proyección usando traza
                f[i,j,k] = np.trace(C @ G[k]).real

    return f

# =========================
# ANÁLISIS FINAL
# =========================

def analyze(G):
    print("\n=== RESULTADO FINAL ===\n")

    for i, g in enumerate(G):
        print(f"G[{i}]:\n{np.round(g, 3)}\n")

    print("T_nc:", non_commutativity_tension(G))
    print("T_norm:", norm_tension(G))

    K = killing_metric(G)
    print("Killing:\n", np.round(K, 3))
    print("T_k:", killing_tension(K))

    f = structure_constants(G)
    print("Estructura (f_ijk):")
    print(np.round(f, 3))

# ========================
# EVOLVE TO SU(2)
# ========================

def evolve_to_su2_fixed(steps=10000, lr=0.01, alpha=0.5, beta=1.0):
    """
    Versión corregida con objetivo fijo
    
    minimize: T_nc + α·T_Killing + β·T_norm
    """
    # Inicializar operadores aleatorios NORMALIZADOS
    generators = random_generators(n=3, dim=2)
    
    # Normalización inicial explícita
    for i in range(len(generators)):
        norm = np.linalg.norm(generators[i])
        if norm > 1e-10:
            generators[i] = generators[i] / norm
    
    history = {
        'T_nc': [],
        'T_k': [],
        'T_norm': [],
        'T_total': [],
        'K_diag': [],
        'K_offdiag': []
    }
    
    for step in range(steps):
        # Calcular tensiones
        T_nc = non_commutativity_tension(generators)
        
        K = killing_metric(generators)
        K_target = -1.0 * np.eye(3)  # ← OBJETIVO FIJO
        T_k = np.linalg.norm(K - K_target)**2
        
        # Constraint de normalización
        T_norm = sum((np.linalg.norm(G) - 1.0)**2 for G in generators)
        
        T_total = T_nc + alpha * T_k + beta * T_norm
        
        # Guardar historia
        history['T_nc'].append(T_nc)
        history['T_k'].append(T_k)
        history['T_norm'].append(T_norm)
        history['T_total'].append(T_total)
        history['K_diag'].append(np.mean(np.diag(K)))
        history['K_offdiag'].append(np.mean(np.abs(K - np.diag(np.diag(K)))))

       
        # Calcular gradiente (igual que antes)
        grad = []
        epsilon = 1e-6
        
        for i, G in enumerate(generators):
            grad_G = np.zeros_like(G)
            
            for r in range(G.shape[0]):
                for c in range(G.shape[1]):
                    # Perturbación +
                    G_plus = G.copy()
                    G_plus[r,c] += epsilon
                    gens_plus = generators.copy()
                    gens_plus[i] = G_plus
                    
                    K_plus = killing_metric(gens_plus)
                    T_plus = (non_commutativity_tension(gens_plus) + 
                             alpha * np.linalg.norm(K_plus - K_target)**2 +
                             beta * sum((np.linalg.norm(Gp) - 1.0)**2 for Gp in gens_plus))
                    
                    # Perturbación -
                    G_minus = G.copy()
                    G_minus[r,c] -= epsilon
                    gens_minus = generators.copy()
                    gens_minus[i] = G_minus
                    
                    K_minus = killing_metric(gens_minus)
                    T_minus = (non_commutativity_tension(gens_minus) + 
                              alpha * np.linalg.norm(K_minus - K_target)**2 +
                              beta * sum((np.linalg.norm(Gm) - 1.0)**2 for Gm in gens_minus))
                    
                    grad_G[r,c] = (T_plus - T_minus) / (2 * epsilon)
            
            grad.append(grad_G)
        
        # Update
        for i in range(len(generators)):
            generators[i] = generators[i] - lr * grad[i]
            
            # Proyectar a hermítico
            generators[i] = (generators[i] + generators[i].conj().T) / 2
            
            # Proyectar a traceless
            generators[i] = generators[i] - np.trace(generators[i]) * np.eye(2) / 2
        
        # Print progreso
        if step % 1000 == 0:
            print(f"Step {step}: T_total = {T_total:.6f}, T_nc = {T_nc:.6f}, T_k = {T_k:.6f}, T_norm = {T_norm:.6f}")
            print(f"  K_diag = {np.diag(K)}")
            
            if T_total < 1e-6:
                print(f"\n¡Convergió en step {step}!")
                break
    
    return generators, history

# Ejecutar versión corregida
print("="*70)
print("VERSIÓN CORREGIDA - Objetivo K = -1·I fijo")
print("="*70)

final_generators, history = evolve_to_su2_fixed(steps=10000, lr=0.001, alpha=5, beta=1.0)

# Verificar resultado
K_final = killing_metric(final_generators)
print("\nMétrica de Killing final:")
print(K_final)
print("\nDiagonal:", np.diag(K_final))
print("Off-diagonal max:", np.max(np.abs(K_final - np.diag(np.diag(K_final)))))

# Comparar con Pauli
pauli = [sigma / 2 for sigma in pauli_matrices()]
K_pauli = killing_metric(pauli)

print("\nComparación con Pauli:")
print("K_final diagonal:", np.diag(K_final))
print("K_pauli diagonal:", np.diag(K_pauli))
print("Ratio:", np.diag(K_final) / np.diag(K_pauli))

# Verificar cierre
for i, j in [(0,1), (1,2), (2,0)]:
    comm = commutator(final_generators[i], final_generators[j])
    # Verificar si está en span
    in_span = check_in_span(comm, final_generators)
    print(f"\nCommutator [G{i}, G{j}] in span? {in_span}")
    # ... (código de verificación)


# Rescalar resultados actuales
K_final = killing_metric(final_generators)
current_scale = np.mean(np.abs(np.diag(K_final)))  # ≈ 0.57
target_scale = 1.0
scale_factor = np.sqrt(target_scale / current_scale)  # ≈ √(1/0.57) ≈ 1.32

print(f"\nRescalando por factor: {scale_factor:.3f}")

generators_rescaled = [G * scale_factor for G in final_generators]

# Verificar Killing
K_rescaled = killing_metric(generators_rescaled)
print("\nK rescalado:")
print(K_rescaled)
print("\nDiagonal:", np.diag(K_rescaled))
print("Comparar con Pauli:", np.diag(K_pauli))

# Verificar cierre aún se mantiene
for i, j in [(0,1), (1,2), (2,0)]:
    comm = commutator(generators_rescaled[i], generators_rescaled[j])
    in_span = check_in_span(comm, generators_rescaled)
    print(f"[G{i}, G{j}] in span? {in_span}")

# Verificar normas
for i, G in enumerate(generators_rescaled):
    print(f"||G{i}|| = {np.linalg.norm(G):.3f}")


# Plot evolución
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

axes[0,0].plot(history['T_total'])
axes[0,0].set_title('Tensión Total')
axes[0,0].set_xlabel('Iteración')
axes[0,0].set_yscale('log')
axes[0,0].grid(True)

axes[0,1].plot(history['T_nc'], label='T_nc')
axes[0,1].plot(history['T_k'], label='T_k')
axes[0,1].set_title('Componentes de Tensión')
axes[0,1].set_xlabel('Iteración')
axes[0,1].set_yscale('log')
axes[0,1].legend()
axes[0,1].grid(True)

axes[1,0].plot(history['K_diag'])
axes[1,0].axhline(np.mean(np.diag(K_pauli)), color='r', linestyle='--', label='Pauli')
axes[1,0].set_title('Diagonal Killing')
axes[1,0].set_xlabel('Iteración')
axes[1,0].legend()
axes[1,0].grid(True)

axes[1,1].plot(history['K_offdiag'])
axes[1,1].set_title('Off-diagonal Killing (debería → 0)')
axes[1,1].set_xlabel('Iteración')
axes[1,1].set_yscale('log')
axes[1,1].grid(True)

plt.tight_layout()
plt.savefig('evolution_to_su2.png', dpi=300)
print("\nGráfico guardado: evolution_to_su2.png")

# =========================
# MAIN
# =========================

if __name__ == "__main__":

    np.random.seed(0)

    G, history = optimize(
        n_iter=5000,
        alpha=1.0,
        beta=1.0
    )

    analyze(G)