import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
from scipy.optimize import linear_sum_assignment

# Riemann zeros (imaginary parts)
true_tn = np.array([
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
    67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
    79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
    92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
    103.725538, 105.446623, 107.168611, 111.029536, 111.874659,
    114.320221, 116.226680, 118.790783, 121.370125, 122.946829,
    124.256819, 127.516684, 129.578704, 131.087689, 133.497737,
    134.756510, 138.116042, 139.736209, 141.123707, 143.111846
])

# Fixed operator parameters
L, kappa, beta, nmax, alpha, dz = 400.0, 9000, 0.06, 2000, 6.0, 0.005
z = np.arange(0, L + dz, dz)
z_i = z[1:-1]

def V_regularized(z, L, kappa, beta, nmax, alpha):
    Vz = np.zeros_like(z)
    for n in range(1, nmax + 1):
        weight = np.exp(-alpha * (n / nmax)**2)
        Vz += (kappa / n) * weight * np.exp(-beta * n * z / L) * np.cos(n * np.pi * z / L)
    return Vz - np.min(Vz) + 1e-6

V_diag = V_regularized(z_i, L, kappa, beta, nmax, alpha)
main_diag = 2 / dz**2 + V_diag
off_diag = -1 / dz**2 * np.ones(len(z_i) - 1)
H = diags([main_diag, off_diag, off_diag], [0, -1, 1], format='csc')

# Extract spectrum
all_sqrt_evals = []
for sigma in range(2000, 8001, 100):
    try:
        evals, _ = eigsh(H, k=500, which='LM', sigma=sigma)
        all_sqrt_evals.extend(np.sqrt(np.abs(evals)))
    except:
        continue

# Match spectrum to Riemann zeros
sqrt_evals = np.sort(np.unique(np.array(all_sqrt_evals)))
cost_matrix = np.abs(true_tn[:, None] - sqrt_evals[None, :])
_, col_ind = linear_sum_assignment(cost_matrix)
matched = sqrt_evals[col_ind]
errors = np.abs(matched - true_tn)

# Output results
print("-- HPO Matching (Top 50, Minimal Output) --")
print(f"{'n':>3} | {'t_n':>10} | {'√λ_n':>12} | {'Error':>10}")
print("-" * 44)
print("Done: No files written, no external dependencies, pure terminal test.")
