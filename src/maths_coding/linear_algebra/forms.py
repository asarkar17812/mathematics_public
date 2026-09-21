"""Bilinear and quadratic forms — Sylvester's law of inertia and the signature.

A symmetric matrix ``A`` defines a **bilinear form** ``B(u, v) = uᵀ A v`` and its **quadratic
form** ``q(v) = vᵀ A v``. The decisive invariant is the **signature** ``(p, n, z)`` — the
numbers of positive, negative, and zero eigenvalues. **Sylvester's law of inertia** says this
triple is well-defined: *however* you diagonalize the form by a change of basis (a
**congruence** ``A ↦ Cᵀ A C`` with ``C`` invertible, not a similarity ``C⁻¹ A C``), you always
get the same counts of ``+``, ``−``, and ``0`` on the diagonal. **Positive-definite** —
signature ``(n, 0, 0)`` — is exactly an inner product.

The distinction between congruence and similarity is the whole point: similarity preserves
*eigenvalues* (it is a change of basis for an operator), while congruence preserves only their
*signs* (it is a change of basis for a form). So the signature, not the spectrum, is what a
quadratic form actually sees — which is why the Minkowski form ``−t² + x² + y² + z²`` of
special relativity is characterized by its signature ``(3, 1, 0)``, and why the second-
derivative test reads only the signature of the Hessian.

The bridge to multilinear algebra is the **Gram form**: for vectors ``v₁, …, v_k`` the matrix
``G_{ij} = vᵢ · vⱼ`` is a quadratic form whose determinant is the *squared volume* of the
parallelepiped they span — and for a square array that is ``|det|²``, recovering ``det`` as the
top exterior power ``Λⁿ`` (volume scaling).

- :func:`signature` — ``(p, n, z)`` from the eigenvalue signs.
- :func:`is_positive_definite` — Sylvester's minor criterion, cross-checked with the signature.
- :func:`congruence_diagonalize` — an exact congruence ``Cᵀ A C = D`` to diagonal form.
- :func:`sylvester_law` — same signature from eigenvalues, from a congruence, and from a random
  congruence ``Cᵀ A C`` — the law of inertia in action.
- :func:`gram_volume` — parallelepiped volume ``√det(G)``, equal to ``|det|`` when square.
"""

import numpy as np
from sage.all import matrix, QQ, ZZ, identity_matrix, random_matrix, copy


def signature(A, tol=1e-9):
    """
    The signature ``(num_positive, num_negative, num_zero)`` of a symmetric matrix ``A``.

    Computed from the eigenvalues (real, since ``A`` is symmetric). Raises ``ValueError`` for a
    non-symmetric matrix.
    """
    M = np.array(A, dtype=float)
    if not np.allclose(M, M.T, atol=tol):
        raise ValueError("signature is defined for symmetric matrices")
    evals = np.linalg.eigvalsh(M)
    pos = int(np.sum(evals > tol))
    neg = int(np.sum(evals < -tol))
    zero = int(np.sum(np.abs(evals) <= tol))
    return (pos, neg, zero)


def is_positive_definite(A):
    """
    Test positive-definiteness via **Sylvester's criterion** (all leading principal minors
    ``> 0``), cross-checked against the signature being ``(n, 0, 0)``.

    Returns
    -------
    dict with ``leading_minors`` (exact, over ℚ), ``by_minors`` (bool), ``by_signature``
    (bool), and ``agree``.
    """
    M = matrix(QQ, A)
    n = M.nrows()
    minors = [M.submatrix(0, 0, k, k).det() for k in range(1, n + 1)]
    by_minors = all(m > 0 for m in minors)
    by_sig = signature(A) == (n, 0, 0)
    return {
        "leading_minors": minors,
        "by_minors": by_minors,
        "by_signature": by_sig,
        "agree": by_minors == by_sig,
    }


def congruence_diagonalize(A, field=QQ):
    """
    Diagonalize the symmetric form ``A`` by an exact **congruence**: return ``(D, C)`` with
    ``D`` diagonal, ``C`` invertible, and ``Cᵀ A C = D``.

    Uses simultaneous symmetric row/column operations (the form-theoretic analogue of Gaussian
    elimination): each column operation is mirrored by the transposed row operation, which keeps
    the matrix symmetric and the transform a congruence. Zero pivots are handled by a swap or,
    failing that, by adding in an off-diagonal direction.
    """
    A = matrix(field, A)
    if not A.is_symmetric():
        raise ValueError("congruence diagonalization needs a symmetric matrix")
    n = A.nrows()
    M = copy(A)
    C = copy(identity_matrix(field, n))   # invariant: Cᵀ A C = M throughout

    i = 0
    while i < n:
        if M[i, i] == 0:
            # bring a nonzero diagonal entry to position (i, i) if one exists below
            swap = next((k for k in range(i + 1, n) if M[k, k] != 0), None)
            if swap is not None:
                M.swap_rows(i, swap)
                M.swap_columns(i, swap)
                C.swap_columns(i, swap)
            else:
                # no nonzero pivot on the diagonal: use an off-diagonal entry to make one
                k = next((j for j in range(i + 1, n) if M[i, j] != 0), None)
                if k is None:
                    i += 1            # whole pivot row/col is zero — a radical (zero) direction
                    continue
                M.add_multiple_of_column(i, k, 1)   # col_i += col_k
                M.add_multiple_of_row(i, k, 1)      # row_i += row_k  ⇒ M[i,i] = 2·M[i,k] ≠ 0
                C.add_multiple_of_column(i, k, 1)

        pivot = M[i, i]
        for j in range(i + 1, n):
            if M[i, j] != 0:
                s = -M[i, j] / pivot
                M.add_multiple_of_column(j, i, s)   # col_j += s·col_i
                M.add_multiple_of_row(j, i, s)      # row_j += s·row_i  (keep symmetric)
                C.add_multiple_of_column(j, i, s)
        i += 1

    return M, C


def _diag_signature(D, tol=1e-9):
    """Signature read directly off the diagonal of a (diagonalized) matrix ``D``."""
    diag = [D[i, i] for i in range(D.nrows())]
    pos = sum(1 for d in diag if d > 0)
    neg = sum(1 for d in diag if d < 0)
    zero = sum(1 for d in diag if d == 0)
    return (pos, neg, zero)


def sylvester_law(A):
    """
    Demonstrate Sylvester's law of inertia: the signature is the **same** computed three ways —
    from the eigenvalue signs, from an exact congruence diagonalization, and from a *random*
    congruence ``Cᵀ A C``.

    Returns
    -------
    dict with ``from_eigenvalues``, ``from_congruence``, ``from_random_congruence``,
    ``congruence_valid`` (``Cᵀ A C = D``), and ``all_agree``.
    """
    A_q = matrix(QQ, A)
    n = A_q.nrows()

    sig_eig = signature(A)

    D, C = congruence_diagonalize(A_q)
    sig_cong = _diag_signature(D)
    cong_valid = (C.transpose() * A_q * C == D)

    # A random congruence by a unimodular (det ±1) integer matrix must preserve the signature.
    U = random_matrix(ZZ, n, algorithm="unimodular")
    B = U.transpose() * A_q * U
    sig_random = signature(B)

    return {
        "from_eigenvalues": sig_eig,
        "from_congruence": sig_cong,
        "from_random_congruence": sig_random,
        "congruence_valid": bool(cong_valid),
        "all_agree": sig_eig == sig_cong == sig_random,
    }


def gram_volume(vectors):
    """
    The volume of the parallelepiped spanned by ``vectors`` (given as rows), via the Gram form.

    ``volume = √det(G)`` with ``G_{ij} = vᵢ · vⱼ``. When the vectors form a square array this
    equals ``|det|`` — the determinant as the top exterior-power (``Λⁿ``) volume scaling.

    Returns
    -------
    dict with ``gram_det``, ``volume``, and (when square) ``abs_det`` and ``matches_abs_det``.
    """
    V = np.array(vectors, dtype=float)
    G = V @ V.T
    gram_det = float(np.linalg.det(G))
    volume = float(np.sqrt(max(gram_det, 0.0)))
    out = {"gram_det": gram_det, "volume": volume}
    if V.shape[0] == V.shape[1]:
        abs_det = abs(float(np.linalg.det(V)))
        out["abs_det"] = abs_det
        out["matches_abs_det"] = bool(np.isclose(volume, abs_det))
    return out


if __name__ == "__main__":
    print("=" * 70)
    print("Quadratic forms — Sylvester's law of inertia and the signature")
    print("=" * 70)

    # ---- The Minkowski form of special relativity: signature (3,1,0) ---------
    minkowski = [[-1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 1, 0],
                 [0, 0, 0, 1]]
    print("\nMinkowski form  diag(-1, 1, 1, 1):")
    print("  signature (p, n, z) =", signature(minkowski), " — the (3,1,0) of spacetime")

    # ---- An off-diagonal indefinite form, diagonalized by congruence ---------
    A = [[0, 1],
         [1, 0]]            # q(x,y) = 2xy, indefinite, no nonzero diagonal pivot
    D, C = congruence_diagonalize(matrix(QQ, A))
    print("\nForm q(x,y) = 2xy,  A = [[0,1],[1,0]]  (a zero-diagonal pivot case):")
    print("  congruence diagonal D =", D.diagonal(), "   Cᵀ A C = D:",
          C.transpose() * matrix(QQ, A) * C == D)
    print("  signature =", signature(A), " (indefinite: one +, one −)")

    # ---- Sylvester's law: same signature three ways --------------------------
    print("\nSylvester's law of inertia (signature is congruence-invariant):")
    for name, M in [("[[2,1],[1,2]] (pos-def)", [[2, 1], [1, 2]]),
                    ("[[1,2],[2,1]] (indef)", [[1, 2], [2, 1]]),
                    ("diag(1,1,0) (degenerate)", [[1, 0, 0], [0, 1, 0], [0, 0, 0]])]:
        s = sylvester_law(M)
        print(f"  {name:<26}: eig {s['from_eigenvalues']}  congruence {s['from_congruence']}"
              f"  random CᵀAC {s['from_random_congruence']}   agree: {s['all_agree']}")

    # ---- Positive-definiteness two ways --------------------------------------
    print("\nPositive-definiteness — Sylvester's minor criterion vs the signature:")
    for M in [[[2, 1], [1, 2]], [[1, 2], [2, 1]]]:
        pd = is_positive_definite(M)
        print(f"  {M}:  leading minors {pd['leading_minors']}  → PD by minors {pd['by_minors']},"
              f" by signature {pd['by_signature']}  (agree: {pd['agree']})")

    # ---- Determinant as volume (the Λⁿ bridge) -------------------------------
    print("\nDeterminant as volume — √det(Gram) of the spanning vectors:")
    gv = gram_volume([[2, 0], [1, 3]])
    print(f"  vectors (2,0),(1,3):  volume = {gv['volume']:.4f},  |det| = {gv['abs_det']:.4f},"
          f"  match: {gv['matches_abs_det']}")
