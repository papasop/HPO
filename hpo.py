import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
from mpmath import zetazero, mp

# Set very high precision for mpmath
mp.dps = 50

# === Final-Precision Extreme Parameters ===
L = 180.0
dz = 0.00125           # Increased spatial resolution
kappa = 1000
beta = 0.1
nmax = 7000            # Higher modal complexity
z = np.arange(0, L + dz, dz)
z_i = z[1:-1]
N_interior = len(z_i)

# Potential function V(z)
def V(z, kappa, beta, nmax):
    Vz = np.zeros_like(z)
    for n in range(1, nmax + 1):
        Vz += (kappa / n) * np.exp(-beta * n * z / L) * np.cos(n * np.pi * z / L)
    return Vz - np.min(Vz) + 1e-6

V_diag = V(z_i, kappa, beta, nmax)
main_diag = 2 / dz**2 + V_diag
off_diag = -1 / dz**2 * np.ones(N_interior - 1)
H = diags([main_diag, off_diag, off_diag], [0, -1, 1], format='csc')

# First 30 Riemann zeta zeros (imaginary parts)
num_tn = 30
true_tn = np.array([float(zetazero(n).imag) for n in range(1, num_tn + 1)])

# σ scan: high resolution + wide range
sigma_list = list(range(200, 11600, 100))
extracted = set()

for sigma in sigma_list:
    try:
        eigvals, _ = eigsh(H, k=40, sigma=sigma, which='LM', tol=1e-10)
        sqrt_vals = np.sqrt(np.sort(eigvals))
        extracted.update(sqrt_vals)
        print(f"σ={sigma}: √λₙ = {sqrt_vals.round(4)}")
    except Exception as e:
        print(f"σ={sigma} failed: {e}")

extracted = np.array(sorted(extracted))

# Greedy matching
matched = []
used = set()
for t in true_tn:
    diffs = np.abs(extracted - t)
    for idx in np.argsort(diffs):
        if idx not in used:
            matched.append(extracted[idx])
            used.add(idx)
            break

matched = np.array(matched)
errors = np.abs(matched - true_tn)

# Output table
print("\n--- HPO Matching: 精度极限构造 (v6) ---")
print(f"{'n':>3} | {'tₙ (ζ)':>10} | {'√λₙ':>10} | {'Error':>10}")
print("-" * 40)
for i in range(num_tn):
    print(f"{i+1:>3} | {true_tn[i]:10.6f} | {matched[i]:10.6f} | {errors[i]:10.6f}")
print("-" * 40)
print(f"L2 error : {np.linalg.norm(errors):.6f}")
print(f"Max error: {np.max(errors):.6f}")

# Plot
plt.figure(figsize=(8, 6))
plt.plot(true_tn, matched, 'o', label='HPO √λₙ')
plt.plot(true_tn, true_tn, '--', color='gray', label='Ideal y = x')
plt.xlabel("Riemann zero tₙ")
plt.ylabel("√λₙ from HPO")
plt.title("HPO Eigenvalue Matching (First 30 Zeta Zeros, v6 Precision)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("hpo_match_first_30_v6.png")
plt.show()
