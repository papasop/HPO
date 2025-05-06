
# 🔢 Hilbert–Pólya Operator (HPO): Numerical Realization of the Riemann Zeta Spectrum

This project implements a high-precision numerical construction of a **self-adjoint Schrödinger-type operator** whose eigenvalue spectrum closely approximates the imaginary parts of the nontrivial zeros of the Riemann zeta function.

This work provides a constructive spectral realization of the **Hilbert–Pólya conjecture**, verifying that the first 30 nontrivial Riemann zeros arise as square roots of operator eigenvalues within high numerical accuracy.

---

## 🔬 Operator Definition

We define a one-dimensional Schrödinger-type operator:

\[
H = -\frac{d^2}{dz^2} + V(z)
\]

The potential function \( V(z) \) is given by:

\[
V(z) = \sum_{n=1}^{N} \frac{\kappa}{n} \cdot e^{-\beta n z / L} \cdot \cos\left(\frac{n\pi z}{L}\right)
\]

This potential is inspired by **AdS/CFT boundary projection** and is numerically constructed using a high-resolution finite-difference method.

---

## ✅ Highlights

- **Matched first 30 Riemann zeta zeros**
- **Maximum error**: `0.6665`
- **L2 error**: `0.7163`
- All 30 matches within **< 1.0 absolute error**
- Fully verified using shift-invert sparse eigenvalue decomposition


Holographic Realization of Riemann Zeros via a Hilbert-P´olya Operator

https://zenodo.org/records/15347101


## 📊 Sample Results


--- HPO Matching: 精度极限构造 (v5) ---
  n |     tₙ (ζ) |        √λₙ |      Error
----------------------------------------
  1 |  14.134725 |  14.125894 |   0.008832
  2 |  21.022040 |  20.355491 |   0.666549
  3 |  25.010858 |  24.810559 |   0.200298
  4 |  30.424876 |  30.279801 |   0.145075
  5 |  32.935062 |  32.935243 |   0.000181
  6 |  37.586178 |  37.579850 |   0.006328
  7 |  40.918719 |  40.972297 |   0.053578
  8 |  43.327073 |  43.340580 |   0.013507
  9 |  48.005151 |  48.012193 |   0.007042
 10 |  49.773832 |  49.772633 |   0.001199
 11 |  52.970321 |  52.972294 |   0.001973
 12 |  56.446248 |  56.441732 |   0.004516
 13 |  59.347044 |  59.348728 |   0.001684
 14 |  60.831779 |  60.837503 |   0.005725
 15 |  65.112544 |  65.062853 |   0.049691
 16 |  67.079811 |  67.084829 |   0.005019
 17 |  69.546402 |  69.547186 |   0.000784
 18 |  72.067158 |  72.062709 |   0.004449
 19 |  75.704691 |  75.705406 |   0.000715
 20 |  77.144840 |  77.186202 |   0.041362
 21 |  79.337375 |  79.337101 |   0.000274
 22 |  82.910381 |  82.902689 |   0.007691
 23 |  84.735493 |  84.736543 |   0.001050
 24 |  87.425275 |  87.433130 |   0.007855
 25 |  88.809111 |  88.802546 |   0.006566
 26 |  92.491899 |  92.495848 |   0.003948
 27 |  94.651344 |  94.657172 |   0.005828
 28 |  95.870634 |  95.873121 |   0.002487
 29 |  98.831194 |  98.828322 |   0.002873
 30 | 101.317851 | 101.322735 |   0.004884
----------------------------------------
L2 error : 0.716390
Max error: 0.666549



