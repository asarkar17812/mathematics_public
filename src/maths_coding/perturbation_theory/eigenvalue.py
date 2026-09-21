"""Rayleigh–Schrödinger perturbation theory for matrix eigenvalues.

Given a matrix that splits as ``H(ε) = H₀ + εV`` with the eigenproblem of ``H₀`` already
solved, how do the eigenvalues and eigenvectors move as ``ε`` turns on? This is the
linear-algebra core of perturbation theory — and, with ``H₀`` a Hamiltonian and ``V`` a
small interaction, exactly the calculation behind quantum-mechanical energy-level shifts.

For a **non-degenerate** eigenpair ``(λₙ⁽⁰⁾, |n⟩)`` of the symmetric matrix ``H₀`` (so the
``|n⟩`` are orthonormal), the corrections are

    λₙ⁽¹⁾ = ⟨n|V|n⟩                                   (first order: just the diagonal of V)
    λₙ⁽²⁾ = Σ_{m≠n} |⟨m|V|n⟩|² / (λₙ⁽⁰⁾ − λₘ⁽⁰⁾)        (second order)
    |n⁽¹⁾⟩ = Σ_{m≠n} ⟨m|V|n⟩ / (λₙ⁽⁰⁾ − λₘ⁽⁰⁾) · |m⟩   (first-order eigenvector)

so the perturbed eigenvalue is ``λₙ(ε) = λₙ⁽⁰⁾ + ε λₙ⁽¹⁾ + ε² λₙ⁽²⁾ + O(ε³)``. The
second-order shift famously pushes neighbouring levels apart ("level repulsion"), since the
``m`` below ``n`` contribute negatively and those above contribute positively.

The denominators ``λₙ⁽⁰⁾ − λₘ⁽⁰⁾`` show why a **degenerate** spectrum needs different
machinery (degenerate perturbation theory); this module requires a simple spectrum and says
so. Computations are numerical (via NumPy); ``H₀`` and ``V`` may be Sage matrices, NumPy
arrays, or nested lists, and should be real symmetric.
"""

import numpy as np

from sage.all import *


def _to_numpy(M):
    """Coerce a Sage matrix, NumPy array, or nested list to a float NumPy array."""
    if hasattr(M, "numpy"):  # Sage matrices expose a .numpy() converter
        return np.array(M.numpy(), dtype=float)
    return np.array(M, dtype=float)


def _eigenbasis(H0, tol):
    """Return ascending eigenvalues and orthonormal eigenvectors (columns) of symmetric H0."""
    lam0, vecs = np.linalg.eigh(H0)
    gaps = [abs(lam0[i] - lam0[j]) for i in range(len(lam0)) for j in range(i)]
    if gaps and min(gaps) < tol:
        raise ValueError(
            "H0 has a (near-)degenerate spectrum; use degenerate perturbation theory")
    return lam0, vecs


def rayleigh_schrodinger(H0, V, order=2, tol=1e-9):
    """
    Eigenvalue corrections of ``H₀ + εV`` to the requested order, per eigenstate.

    Parameters
    ----------
    H0, V : Sage matrix / NumPy array / nested list
        Real symmetric matrices; ``H0`` must have a simple spectrum.
    order : int
        1 or 2 — the highest perturbation order to return.
    tol : float
        Degeneracy threshold for the spectrum of ``H0``.

    Returns
    -------
    list[list[float]]
        For each eigenstate (ordered by ascending unperturbed eigenvalue) the coefficient
        list ``[λ⁽⁰⁾, λ⁽¹⁾, λ⁽²⁾, …]`` up to ``order``.
    """
    H0 = _to_numpy(H0)
    V = _to_numpy(V)
    lam0, vecs = _eigenbasis(H0, tol)
    n = H0.shape[0]

    out = []
    for i in range(n):
        vi = vecs[:, i]
        coeffs = [float(lam0[i])]
        if order >= 1:
            coeffs.append(float(vi @ V @ vi))
        if order >= 2:
            l2 = sum((vecs[:, m] @ V @ vi) ** 2 / (lam0[i] - lam0[m])
                     for m in range(n) if m != i)
            coeffs.append(float(l2))
        out.append(coeffs)
    return out


def perturbed_eigenvalues(H0, V, eps, order=2, tol=1e-9):
    """
    Approximate eigenvalues of ``H₀ + εV`` at the given ``eps``, using corrections to ``order``.

    Returns
    -------
    list[float]
        ``λₙ(ε) ≈ Σ_k λₙ⁽ᵏ⁾ εᵏ``, ascending.
    """
    corr = rayleigh_schrodinger(H0, V, order=order, tol=tol)
    return [sum(c[k] * eps ** k for k in range(len(c))) for c in corr]


def first_order_eigenvectors(H0, V, tol=1e-9):
    """
    First-order eigenvector corrections ``|n⁽¹⁾⟩ = Σ_{m≠n} ⟨m|V|n⟩/(λₙ−λₘ) |m⟩``.

    Returns
    -------
    list[numpy.ndarray]
        One correction vector per eigenstate (ascending order). The perturbed eigenvector is
        ``|n⟩ + ε|n⁽¹⁾⟩ + O(ε²)``.
    """
    H0 = _to_numpy(H0)
    V = _to_numpy(V)
    lam0, vecs = _eigenbasis(H0, tol)
    n = H0.shape[0]

    corrections = []
    for i in range(n):
        vi = vecs[:, i]
        corr = np.zeros(n)
        for m in range(n):
            if m == i:
                continue
            corr = corr + (vecs[:, m] @ V @ vi) / (lam0[i] - lam0[m]) * vecs[:, m]
        corrections.append(corr)
    return corrections


if __name__ == "__main__":
    # H0 has a simple spectrum; V is a symmetric perturbation.
    H0 = matrix(RDF, [[1, 0, 0], [0, 2.5, 0], [0, 0, 4]])
    V = matrix(RDF, [[0, 0.3, 0.2], [0.3, 0, 0.5], [0.2, 0.5, 0]])

    print("Rayleigh-Schrodinger corrections [lambda0, lambda1, lambda2]:")
    for i, c in enumerate(rayleigh_schrodinger(H0, V, order=2)):
        print(f"  state {i}: {[round(v, 5) for v in c]}")

    print("\nApprox vs exact eigenvalues of H0 + eps V:")
    for eps in (0.05, 0.1):
        approx = sorted(perturbed_eigenvalues(H0, V, eps, order=2))
        exact = sorted(np.linalg.eigvalsh(_to_numpy(H0) + eps * _to_numpy(V)))
        err = max(abs(a - b) for a, b in zip(approx, exact))
        print(f"  eps={eps}:  max error = {err:.2e}   (~eps^3 = {eps**3:.1e})")
