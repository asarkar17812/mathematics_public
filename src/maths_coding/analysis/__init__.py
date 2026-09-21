"""Analysis — harmonic analysis and (later) the analytic core.

For now this subpackage centers on **Fourier analysis on finite abelian groups**, the cleanest
case of Pontryagin duality:

- :mod:`fourier` — the discrete Fourier transform as the change of basis into characters, the
  **convolution theorem** (convolution becomes pointwise multiplication), **Parseval/Plancherel**
  (the transform is an isometry), and **Pontryagin duality** ``\\widehat{ℤ/N} ≅ ℤ/N`` with the
  orthogonality of characters.

This is the abelian shadow of :mod:`maths_coding.algebra.representation_theory`: the irreducible
representations of a finite abelian group are its one-dimensional **characters**, and the DFT is
exactly the decomposition of functions into that orthonormal character basis.
"""
