"""Three canonical factorizations: the spectral theorem, the SVD, and the Jordan form.

A matrix is understood by the basis in which it becomes simple. Three theorems supply that
basis under three hypotheses:

**Spectral theorem (symmetric ``A``).** A real symmetric matrix has an *orthonormal* basis of
eigenvectors: ``A = Q Λ Qᵀ`` with ``QᵀQ = I`` and ``Λ`` real diagonal. Geometrically ``A``
is a pure scaling along ``n`` perpendicular axes — no rotation, no shear. This is the engine
behind principal axes, the second-derivative test, and the curvature operators elsewhere in
the package.

**Singular value decomposition (any ``A``).** Every ``m × n`` matrix factors as
``A = U Σ Vᵀ`` with ``U, V`` orthogonal and ``Σ ≥ 0`` diagonal. It is the precise meaning of
the "red arrow" of the four-subspaces picture: in the right orthonormal bases, ``A`` is
*rotate → scale → rotate*. Concretely it carries the **unit sphere to an ellipsoid** whose
semi-axes are the singular values ``σᵢ`` (the lengths) along the left singular vectors ``uᵢ``
(the directions), and ``σᵢ² `` are the eigenvalues of ``AᵀA``.

**Jordan form (defective ``A``).** When eigenvectors are too few to form a basis, the matrix
is *not* diagonalizable; the best one can do is ``A = P J P⁻¹`` with ``J`` block-diagonal in
**Jordan blocks** — an eigenvalue on the diagonal and ``1``\\ s on the superdiagonal, the
algebraic-but-not-geometric multiplicity made visible.

- :func:`spectral_decomposition` — ``Q, Λ`` for symmetric ``A``, with orthonormality checked.
- :func:`singular_value_decomposition` — ``U, Σ, Vᵀ``, ``σᵢ``, and the ``σᵢ² = eig(AᵀA)`` link.
- :func:`ellipse_semi_axes` — the semi-axes (``σᵢ``) and axis directions (``uᵢ``) of ``A``·sphere.
- :func:`jordan_decomposition` — ``J, P`` and the algebraic-vs-geometric multiplicity gap.
- :func:`plot_svd_action` — (optional) draw the unit circle mapped to its image ellipse.
"""

import numpy as np
from sage.all import matrix, QQ


# --------------------------------------------------------------------------- #
#  Spectral theorem (symmetric matrices)                                      #
# --------------------------------------------------------------------------- #

def spectral_decomposition(A, tol=1e-9):
    """
    Spectral decomposition ``A = Q Λ Qᵀ`` of a real **symmetric** matrix.

    Uses ``numpy.linalg.eigh`` (exact-arithmetic-free but numerically stable for symmetric
    input), returning an orthonormal eigenbasis. Raises ``ValueError`` if ``A`` is not
    symmetric.

    Returns
    -------
    dict with ``eigenvalues`` (array), ``Q`` (orthonormal eigenvectors as columns),
    ``orthonormal`` (``QᵀQ ≈ I``), and ``reconstructs`` (``QΛQᵀ ≈ A``).
    """
    A = np.array(A, dtype=float)
    if not np.allclose(A, A.T, atol=tol):
        raise ValueError("spectral theorem (this form) needs a symmetric matrix")
    vals, Q = np.linalg.eigh(A)
    recon = Q @ np.diag(vals) @ Q.T
    return {
        "eigenvalues": vals,
        "Q": Q,
        "orthonormal": bool(np.allclose(Q.T @ Q, np.eye(A.shape[0]), atol=tol)),
        "reconstructs": bool(np.allclose(recon, A, atol=tol)),
    }


# --------------------------------------------------------------------------- #
#  Singular value decomposition (any matrix)                                  #
# --------------------------------------------------------------------------- #

def singular_value_decomposition(A, tol=1e-9):
    """
    Singular value decomposition ``A = U Σ Vᵀ`` with the ``σᵢ² = eig(AᵀA)`` cross-check.

    Returns
    -------
    dict with ``U``, ``singular_values`` (``σ``, descending), ``Vt``, ``reconstructs``
    (``UΣVᵀ ≈ A``), and ``sigma_sq_are_AtA_eigs`` (``σᵢ²`` match the eigenvalues of ``AᵀA``).
    """
    A = np.array(A, dtype=float)
    U, s, Vt = np.linalg.svd(A)
    m, n = A.shape
    Sigma = np.zeros((m, n))
    for i in range(min(m, n)):
        Sigma[i, i] = s[i]
    recon = U @ Sigma @ Vt

    ata_eigs = np.sort(np.linalg.eigvalsh(A.T @ A))[::-1]
    sigma_sq = np.zeros(n)
    sigma_sq[: len(s)] = s ** 2
    return {
        "U": U,
        "singular_values": s,
        "Vt": Vt,
        "reconstructs": bool(np.allclose(recon, A, atol=tol)),
        "sigma_sq_are_AtA_eigs": bool(np.allclose(sigma_sq[: n], ata_eigs[: n], atol=1e-6)),
    }


def ellipse_semi_axes(A):
    """
    The semi-axes of the image of the unit sphere under ``A``: lengths ``σᵢ`` along the left
    singular directions ``uᵢ``.

    Returns
    -------
    list of ``(length, direction)`` pairs, longest axis first.
    """
    A = np.array(A, dtype=float)
    U, s, _ = np.linalg.svd(A)
    return [(float(s[i]), U[:, i]) for i in range(len(s))]


# --------------------------------------------------------------------------- #
#  Jordan canonical form (defective matrices)                                 #
# --------------------------------------------------------------------------- #

def jordan_decomposition(A, field=QQ):
    """
    Jordan form ``A = P J P⁻¹`` (exact, over ``field``) and the multiplicity gap that forces it.

    For each eigenvalue, ``algebraic`` multiplicity is its multiplicity as a root of the
    characteristic polynomial and ``geometric`` multiplicity is ``dim ker(A − λI)``; the matrix
    is diagonalizable exactly when these agree for every ``λ``. When they differ, Jordan blocks
    of size ``> 1`` appear.

    Returns
    -------
    dict with ``J``, ``P``, ``reconstructs`` (``P J P⁻¹ = A``), ``diagonalizable`` (bool),
    and ``multiplicities`` (``{λ: (algebraic, geometric)}``).
    """
    A = matrix(field, A)
    J, P = A.jordan_form(transformation=True)
    n = A.nrows()
    I = matrix.identity(field, n)

    mult = {}
    for lam, alg in A.charpoly().roots():
        geom = (A - lam * I).right_kernel().dimension()
        mult[lam] = (int(alg), int(geom))

    diagonalizable = all(alg == geom for alg, geom in mult.values())
    return {
        "J": J,
        "P": P,
        "reconstructs": bool(P * J * P.inverse() == A),
        "diagonalizable": diagonalizable,
        "multiplicities": mult,
    }


def plot_svd_action(A, num_points=200):
    """
    (Optional) Draw the unit circle and its image ellipse ``A·circle`` with the singular axes.

    Requires ``matplotlib``; returns the Matplotlib figure (the caller may ``show()`` it).
    Intended for ``2 × 2`` real matrices.
    """
    import matplotlib.pyplot as plt

    A = np.array(A, dtype=float)
    theta = np.linspace(0, 2 * np.pi, num_points)
    circle = np.vstack([np.cos(theta), np.sin(theta)])
    image = A @ circle

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(circle[0], circle[1], "b--", label="unit circle")
    ax.plot(image[0], image[1], "r-", label="A · circle (ellipse)")
    for length, direction in ellipse_semi_axes(A):
        ax.annotate("", xy=length * direction, xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color="k"))
    ax.set_aspect("equal")
    ax.legend()
    ax.set_title("SVD: the unit circle becomes an ellipse with semi-axes σᵢ")
    return fig


if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)

    print("=" * 70)
    print("Spectral theorem — a symmetric matrix is a scaling along ⟂ axes")
    print("=" * 70)
    S = [[2, 1],
         [1, 2]]
    sp = spectral_decomposition(S)
    print("\nA = [[2,1],[1,2]]  (symmetric)")
    print("  eigenvalues :", sp["eigenvalues"], " (expected 1 and 3)")
    print("  Q orthonormal:", sp["orthonormal"], "   QΛQᵀ = A:", sp["reconstructs"])
    print("  eigen-axes are the diagonals (1,1)/√2 and (1,-1)/√2")

    print("\n" + "=" * 70)
    print("SVD — any matrix sends the unit circle to an ellipse")
    print("=" * 70)
    M = [[3, 0],
         [4, 5]]
    sv = singular_value_decomposition(M)
    print("\nA = [[3,0],[4,5]]")
    print("  singular values σ :", sv["singular_values"])
    print("  UΣVᵀ = A          :", sv["reconstructs"],
          "    σᵢ² = eig(AᵀA) :", sv["sigma_sq_are_AtA_eigs"])
    print("  image-ellipse semi-axes (length along direction):")
    for length, direction in ellipse_semi_axes(M):
        print(f"    σ = {length:.4f}  along  {np.round(direction, 4)}")

    print("\n" + "=" * 70)
    print("Jordan form — a defective matrix cannot be diagonalized")
    print("=" * 70)
    D = [[4, 1],
         [-1, 2]]
    jd = jordan_decomposition(D)
    print("\nA = [[4,1],[-1,2]]   (characteristic poly (λ-3)²)")
    print("  algebraic vs geometric multiplicity:",
          {f"λ={lam}": ag for lam, ag in jd["multiplicities"].items()})
    print("  diagonalizable :", jd["diagonalizable"], " (geometric < algebraic ⇒ defective)")
    print("  Jordan form J =")
    print(jd["J"])
    print("  P J P⁻¹ = A :", jd["reconstructs"], " — one 2×2 Jordan block for λ = 3")
