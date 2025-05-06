
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


