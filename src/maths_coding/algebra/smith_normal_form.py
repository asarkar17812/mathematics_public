"""Smith Normal Form: "Gaussian elimination that remembers torsion".

Over a *field*, row and column operations reduce any matrix to a block of ``1``\\ s and
``0``\\ s, and the only invariant is the rank. Over a PID we may no longer divide freely —
only multiply rows/columns by *units* and add integer multiples — so the reduction cannot
kill every entry. What survives on the diagonal is a divisibility chain
``d₁ | d₂ | ⋯ | d_r`` of **invariant factors** that records exactly the torsion of the
module the matrix presents. The two-sided action ``PAQ`` is a simultaneous change of basis
in the *source* and the *target*, chosen to make the linear map as diagonal as the ring
allows:

    Smith Normal Form (PID)
    ───────────────────────
    for A ∈ M_{m×n}(R) there are P ∈ GLₘ(R), Q ∈ GLₙ(R) with
        PAQ = diag(d₁, …, d_r, 0, …, 0),   r = rank(A),   d₁ | d₂ | ⋯ | d_r.

The reason the form matters is the **structure theorem for finitely generated modules over
a PID**. If ``M ≅ Rⁿ / im(A)`` is presented by the relation matrix ``A`` (its columns are
the relations), then

    M ≅ R^{n-r} ⊕ R/(d₁) ⊕ ⋯ ⊕ R/(d_r),

i.e. ``M`` is exactly the **cokernel** of ``A`` (``coker(A) = Rⁿ/im(A)``), computed
explicitly by the Smith form. For ``R = ℤ`` this is the classification of finitely
generated abelian groups; the invariant factors are the torsion, and the number of trailing
zeros is the free rank.

- :func:`smith_normal_form` — the triple ``(D, U, V)`` with ``D = U A V`` in Smith form.
- :func:`invariant_factors` — the divisibility chain ``d₁ | ⋯ | d_r`` (units kept visible).
- :func:`abelian_group_structure` / :func:`cokernel_structure` — read off
  ``ℤ^{free} ⊕ ⊕ ℤ/dᵢ`` from a relation matrix.
- :func:`solve_diophantine_system` — use the Smith form to decouple ``A x = b`` over ℤ into
  independent congruences and return a particular solution plus a kernel basis.

The worked examples in ``__main__`` reproduce the notes: ``ℤ²/⟨(2,2),(4,6)⟩ ≅ ℤ/2 ⊕ ℤ/2``
(the Klein four-group), and ``diag(6,4) ↝ diag(2,12)``, i.e. ``ℤ/6 ⊕ ℤ/4 ≅ ℤ/2 ⊕ ℤ/12`` —
the Chinese Remainder Theorem reading the invariant factors off the diagonal.
"""

from sage.all import matrix, ZZ, vector, gcd


def smith_normal_form(A, ring=ZZ):
    """
    Smith Normal Form of ``A`` over a PID.

    Returns the triple ``(D, U, V)`` with ``D = U * A * V``, where ``U`` and ``V`` are
    invertible over ``ring`` and ``D`` is diagonal with the divisibility chain
    ``D[0,0] | D[1,1] | ⋯``. (Sage's convention names the two-sided transform
    ``U A V`` rather than ``P A Q``; the content is identical.)

    Parameters
    ----------
    A : matrix or list of lists
        The matrix to reduce.
    ring : Sage ring, default ``ZZ``
        The PID to work over.

    Returns
    -------
    (D, U, V) : tuple of Sage matrices
        ``D`` is the Smith form; ``U``, ``V`` are the unimodular transforms.
    """
    A = matrix(ring, A)
    D, U, V = A.smith_form()
    return D, U, V


def invariant_factors(A, ring=ZZ):
    """
    The invariant factors ``d₁ | d₂ | ⋯ | d_r`` (the nonzero diagonal of the Smith form).

    These are the divisor chain that classifies ``coker(A)``. Units (``dᵢ = 1``) are kept
    so the divisibility chain is visible; they contribute a *trivial* summand ``R/(1) = 0``
    to the module and can be dropped when reading off the structure.

    Returns
    -------
    list
        ``[d₁, …, d_r]`` with ``r = rank(A)`` and ``dᵢ | d_{i+1}``.
    """
    D, _, _ = smith_normal_form(A, ring)
    diag = [D[i, i] for i in range(min(D.nrows(), D.ncols()))]
    return [d for d in diag if d != 0]


def cokernel_structure(A, ring=ZZ):
    """
    Structure of ``coker(A) = R^m / im(A)`` for the relation matrix ``A`` (m = #rows).

    Interprets the columns of ``A`` as relations among ``m`` generators of a free module
    ``R^m``, and applies the structure theorem: ``coker(A) ≅ R^{m-r} ⊕ ⊕ᵢ R/(dᵢ)`` where
    the ``dᵢ`` are the non-unit invariant factors.

    Returns
    -------
    dict with keys
        ``free_rank`` : int
            ``m - rank(A)`` — the number of free ``R`` summands.
        ``torsion`` : list
            The invariant factors ``dᵢ > 1`` (the unit factors are dropped).
        ``label`` : str
            A human-readable ``R^k ⊕ R/(d₁) ⊕ ⋯`` description.
    """
    A = matrix(ring, A)
    m = A.nrows()
    r = A.rank()
    factors = invariant_factors(A, ring)
    torsion = [d for d in factors if d != 1 and d != -1]
    free_rank = m - r

    R = repr(ring) if ring is not ZZ else "Z"
    pieces = []
    if free_rank == 1:
        pieces.append(f"{R}")
    elif free_rank > 1:
        pieces.append(f"{R}^{free_rank}")
    pieces.extend(f"{R}/{d}" for d in torsion)
    label = " ⊕ ".join(pieces) if pieces else "0"

    return {"free_rank": free_rank, "torsion": torsion, "label": label}


def abelian_group_structure(relation_matrix):
    """
    Classify the finitely generated abelian group ``ℤ^m / ⟨columns of A⟩``.

    Thin wrapper over :func:`cokernel_structure` for ``R = ℤ``; the columns of
    ``relation_matrix`` are the relations among the generators.

    Examples
    --------
    The Klein four-group from the notes, ``M = ℤ² / ⟨(2,2), (4,6)⟩``::

        abelian_group_structure([[2, 4], [2, 6]])  # -> label 'Z/2 ⊕ Z/2'
    """
    return cokernel_structure(relation_matrix, ZZ)


def solve_diophantine_system(A, b, ring=ZZ):
    """
    Solve the integer linear system ``A x = b`` using the Smith Normal Form.

    With ``D = U A V`` (Smith form), substituting ``y = V⁻¹ x`` turns the coupled system
    into a diagonal one, ``D y = U b`` — a list of *independent congruences* ``dᵢ yᵢ = cᵢ``.
    This is the Smith-form application of the notes: the change of basis decouples the
    system. Solvability is then immediate, and the free coordinates (where ``dᵢ = 0``) give
    a basis for the kernel of ``A`` over ``ℤ``.

    Parameters
    ----------
    A : matrix or list of lists
    b : vector or list
    ring : Sage ring, default ``ZZ``

    Returns
    -------
    dict with keys
        ``solvable`` : bool
        ``particular`` : Sage vector or None
            One integer solution ``x`` with ``A x = b`` (None if unsolvable).
        ``kernel_basis`` : list of Sage vectors
            A ℤ-basis for ``{x : A x = 0}``.
    """
    A = matrix(ring, A)
    b = vector(ring, b)
    m, n = A.nrows(), A.ncols()

    D, U, V = A.smith_form()  # D = U A V
    c = U * b                 # D y = c, with y = V^{-1} x  ->  x = V y

    # Solve the decoupled diagonal system D y = c.
    y = [ring(0)] * n
    free_positions = []
    for i in range(n):
        d = D[i, i] if i < min(m, n) else ring(0)
        ci = c[i] if i < m else ring(0)
        if d == 0:
            if ci != 0:
                return {"solvable": False, "particular": None, "kernel_basis": []}
            free_positions.append(i)  # y_i free
        else:
            if ci % d != 0:
                return {"solvable": False, "particular": None, "kernel_basis": []}
            y[i] = ci // d
    # Any extra rows of c beyond n must also vanish (over-determined consistency).
    for i in range(n, m):
        if c[i] != 0:
            return {"solvable": False, "particular": None, "kernel_basis": []}

    particular = V * vector(ring, y)

    # Kernel basis: free y-directions pushed back through x = V y.
    kernel_basis = []
    for j in free_positions:
        e = [ring(0)] * n
        e[j] = ring(1)
        kernel_basis.append(V * vector(ring, e))

    return {"solvable": True, "particular": particular, "kernel_basis": kernel_basis}


if __name__ == "__main__":
    print("=" * 70)
    print("Smith Normal Form — torsion made visible")
    print("=" * 70)

    # ---- Worked example 1: an abelian group from a presentation matrix --------
    # M = Z^2 / <(2,2), (4,6)>; columns of A are the relations.
    A = [[2, 4],
         [2, 6]]
    D, U, V = smith_normal_form(A)
    print("\nA =")
    print(matrix(ZZ, A))
    print("\nSmith form D = U A V =")
    print(D)
    print("invariant factors d1 | d2 :", invariant_factors(A))
    info = abelian_group_structure(A)
    print("structure of Z^2/im(A)  :", info["label"], " (the Klein four-group)")
    print("|det A| = d1*d2         :", abs(matrix(ZZ, A).det()),
          "=", " * ".join(str(d) for d in invariant_factors(A)))

    # ---- Worked example 2: diag(6,4) is NOT in Smith form --------------------
    # 6 does not divide 4; reducing reveals Z/6 + Z/4 = Z/2 + Z/12 (CRT).
    B = [[6, 0],
         [0, 4]]
    print("\n" + "-" * 70)
    print("diag(6,4) presents Z/6 ⊕ Z/4, but 6 ∤ 4, so it is not yet canonical.")
    print("invariant factors :", invariant_factors(B))
    print("structure         :", cokernel_structure(B)["label"],
          " (so Z/6 ⊕ Z/4 ≅ Z/2 ⊕ Z/12 — CRT)")

    # ---- A non-square example with a free summand ----------------------------
    # 3 generators (rows = ambient Z^3), only 2 relations (columns) -> a free Z survives.
    C = [[2, 0],
         [0, 6],
         [0, 0]]
    print("\n" + "-" * 70)
    print("A 3×2 relation matrix (3 generators, 2 relations) leaves a free summand:")
    print(matrix(ZZ, C))
    print("structure :", cokernel_structure(C)["label"])

    # ---- Smith form solving a linear Diophantine system ----------------------
    print("\n" + "-" * 70)
    print("Solving an integer system A x = b via the decoupled Smith system:")
    Asys = [[2, 4],
            [2, 6]]
    bvec = [10, 14]
    sol = solve_diophantine_system(Asys, bvec)
    print("A =", matrix(ZZ, Asys).rows(), " b =", bvec)
    print("solvable      :", sol["solvable"])
    if sol["solvable"]:
        print("particular x  :", list(sol["particular"]),
              " check A x =", list(matrix(ZZ, Asys) * sol["particular"]))
        print("kernel basis  :", [list(k) for k in sol["kernel_basis"]])
