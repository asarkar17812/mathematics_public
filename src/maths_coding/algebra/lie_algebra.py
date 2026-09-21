"""The Lie algebra sl₂ and its representations — the atom of highest-weight theory.

``sl₂`` is the three-dimensional Lie algebra with basis ``H, E, F`` and the **commutation
relations**

    [H, E] = 2E,    [H, F] = −2F,    [E, F] = H.

Here ``H`` is the diagonal "weight" operator, ``E`` raises weights by ``2`` and ``F`` lowers them
by ``2`` — the **ladder operators**. Its finite-dimensional irreducible representations are the
cleanest case of the highest-weight classification: for each integer ``n ≥ 0`` there is a unique
irreducible ``V_n`` of dimension ``n + 1``, with a basis on which ``H`` is diagonal with
**weights** ``n, n−2, …, −n`` (symmetric about ``0``), ``E`` shrinks toward the highest weight,
and ``F`` walks down the ladder. Setting ``n = 2j`` recovers the **spin-``j`` representation of
``su(2)``** familiar from quantum mechanics.

The invariant that detects irreducibility is the **Casimir element**
``Ω = EF + FE + ½H²``: it commutes with the whole algebra, so by Schur's lemma it acts on an
irreducible representation as a single scalar — here ``Ω = ½ n(n+2)·I`` on ``V_n``. This is the
representation-theoretic cousin of :mod:`maths_coding.algebra.representation_theory` (irreducible
``⟺`` the center acts by a scalar) and of the eigenvalue/operator motif that runs through the
package.

Pure-Python and exact (matrices of :class:`fractions.Fraction`), so it runs without SageMath:

- :func:`sl2_irrep` — the matrices ``H, E, F`` of the ``(n+1)``-dimensional irreducible ``V_n``.
- :func:`commutator` — ``[X, Y] = XY − YX``.
- :func:`check_relations` — verify the three defining ``sl₂`` relations.
- :func:`weights` — the diagonal of ``H``: the weights ``n, n−2, …, −n``.
- :func:`casimir` — ``Ω = EF + FE + ½H²`` and its scalar value ``½ n(n+2)``.
"""

from fractions import Fraction


def _zeros(d):
    return [[Fraction(0) for _ in range(d)] for _ in range(d)]


def _matmul(A, B):
    d = len(A)
    C = _zeros(d)
    for i in range(d):
        for k in range(d):
            if A[i][k]:
                aik = A[i][k]
                for j in range(d):
                    C[i][j] += aik * B[k][j]
    return C


def _add(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A))] for i in range(len(A))]


def _scale(c, A):
    c = Fraction(c)
    return [[c * A[i][j] for j in range(len(A))] for i in range(len(A))]


def commutator(X, Y):
    """The Lie bracket ``[X, Y] = XY − YX`` of two square matrices."""
    return _add(_matmul(X, Y), _scale(-1, _matmul(Y, X)))


def sl2_irrep(n):
    """
    The ``(n+1)``-dimensional irreducible representation ``V_n`` of ``sl₂``: matrices ``H, E, F``.

    Basis ``e₀, …, eₙ`` with ``e₀`` the highest-weight vector:
    ``H eₖ = (n−2k) eₖ``, ``F eₖ = e_{k+1}``, ``E eₖ = k(n−k+1) e_{k−1}``. These satisfy the
    ``sl₂`` relations exactly (rational/integer entries).

    Returns
    -------
    (H, E, F) : tuple of ``(n+1)×(n+1)`` matrices of :class:`fractions.Fraction`.
    """
    d = n + 1
    H, E, F = _zeros(d), _zeros(d), _zeros(d)
    for k in range(d):
        H[k][k] = Fraction(n - 2 * k)
        if k + 1 < d:
            F[k + 1][k] = Fraction(1)          # F e_k = e_{k+1}
        if k - 1 >= 0:
            E[k - 1][k] = Fraction(k * (n - k + 1))  # E e_k = k(n-k+1) e_{k-1}
    return H, E, F


def check_relations(n):
    """
    Verify the defining relations ``[H,E]=2E``, ``[H,F]=−2F``, ``[E,F]=H`` for ``V_n``.

    Returns a dict of three booleans plus ``all``.
    """
    H, E, F = sl2_irrep(n)
    r1 = commutator(H, E) == _scale(2, E)
    r2 = commutator(H, F) == _scale(-2, F)
    r3 = commutator(E, F) == H
    return {"[H,E]=2E": r1, "[H,F]=-2F": r2, "[E,F]=H": r3, "all": r1 and r2 and r3}


def weights(n):
    """The weights of ``V_n`` — the diagonal of ``H`` — namely ``n, n−2, …, −n``."""
    H, _, _ = sl2_irrep(n)
    return [int(H[k][k]) for k in range(n + 1)]


def casimir(n):
    """
    The Casimir operator ``Ω = EF + FE + ½H²`` on ``V_n`` and its scalar value.

    By Schur's lemma ``Ω`` acts as a scalar on the irreducible ``V_n``; that scalar is
    ``½ n(n+2)``. Returns a dict with ``matrix`` (``Ω``), ``is_scalar`` (``Ω = c·I``),
    ``scalar`` (the value ``c``), and ``expected`` (``½ n(n+2)``).
    """
    H, E, F = sl2_irrep(n)
    d = n + 1
    omega = _add(_add(_matmul(E, F), _matmul(F, E)), _scale(Fraction(1, 2), _matmul(H, H)))
    c = omega[0][0]
    is_scalar = all(omega[i][j] == (c if i == j else 0) for i in range(d) for j in range(d))
    return {
        "matrix": omega,
        "is_scalar": is_scalar,
        "scalar": c,
        "expected": Fraction(n * (n + 2), 2),
    }


if __name__ == "__main__":
    print("=" * 70)
    print("The Lie algebra sl₂ and its irreducible representations V_n")
    print("=" * 70)

    print("\nDefining relations [H,E]=2E, [H,F]=−2F, [E,F]=H hold in every V_n:")
    for n in range(0, 6):
        chk = check_relations(n)
        print(f"  V_{n} (dim {n+1}): weights {weights(n)!s:<22} relations hold: {chk['all']}")

    print("\nThe Casimir Ω = EF + FE + ½H² acts as the scalar ½·n(n+2) (Schur's lemma):")
    for n in range(0, 6):
        cas = casimir(n)
        print(f"  V_{n}: Ω = {str(cas['scalar']):>4}·I   = ½·{n}·{n+2} = {cas['expected']}   "
              f"scalar: {cas['is_scalar']}  matches: {cas['scalar'] == cas['expected']}")

    print("\n(Setting n = 2j gives the spin-j representation of su(2): "
          "V_0 spin-0, V_1 spin-½, V_2 spin-1, …)")
