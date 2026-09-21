"""The four fundamental subspaces — Strang's "big picture" of a matrix.

Every ``m × n`` matrix ``A`` of rank ``r`` organizes two spaces into four pieces:

    domain ℝⁿ   = row space (dim r)      ⊕  null space (dim n − r)
    codomain ℝᵐ = column space (dim r)   ⊕  left null space (dim m − r)

with two structural facts that this module makes explicit and checks:

    Rank–nullity:   dim(null space) + rank = n,    dim(row space) = dim(column space) = r.
    Orthogonality:  row space ⊥ null space  (in ℝⁿ),   column space ⊥ left null space (in ℝᵐ).

The real content is the **red arrow** of the Atlas diagram: ``A`` restricted to the row space
is an **isomorphism onto the column space**. Stripped of its kernels, *every* matrix is a
bijection — and that bijection, measured in orthonormal bases, is exactly the singular value
decomposition (see :mod:`maths_coding.linear_algebra.decompositions`). Solving ``Ax = b`` is
then two questions answered by two subspaces: *solvable?* (is ``b`` in the column space) and
*unique?* (is the null space trivial).

- :func:`four_subspaces` — the four subspaces, their dimensions, and the rank–nullity ledger.
- :func:`check_orthogonality` — row ⊥ null and column ⊥ left-null, verified on bases.
- :func:`row_to_column_isomorphism` — ``A`` carries a row-space basis to a column-space basis.
"""

from sage.all import matrix, QQ


def four_subspaces(A, field=QQ):
    """
    The four fundamental subspaces of ``A`` together with the rank–nullity ledger.

    Returns
    -------
    dict with the four subspaces (``row_space``, ``null_space``, ``column_space``,
    ``left_null_space``), their dimensions, ``rank``, and ``rank_nullity_holds`` (the identity
    ``dim(null) + rank = n`` together with ``dim(row) = dim(col) = rank``).
    """
    A = matrix(field, A)
    m, n = A.nrows(), A.ncols()
    r = A.rank()

    row = A.row_space()
    null = A.right_kernel()
    col = A.column_space()
    left_null = A.left_kernel()

    return {
        "rank": int(r),
        "n": int(n),
        "m": int(m),
        "row_space": row,
        "null_space": null,
        "column_space": col,
        "left_null_space": left_null,
        "dims": {
            "row": int(row.dimension()),
            "null": int(null.dimension()),
            "column": int(col.dimension()),
            "left_null": int(left_null.dimension()),
        },
        "rank_nullity_holds": (
            null.dimension() + r == n
            and row.dimension() == r
            and col.dimension() == r
        ),
    }


def check_orthogonality(A, field=QQ):
    """
    Verify the two orthogonality relations on subspace bases:
    row space ⊥ null space (in the domain) and column space ⊥ left null space (in the codomain).

    The first holds because a null vector is killed by every row of ``A``; the second is the
    same statement for ``Aᵀ``. Returns a dict with ``row_perp_null`` and ``col_perp_leftnull``.
    """
    A = matrix(field, A)

    def all_orthogonal(U, W):
        return all(u.dot_product(w) == 0 for u in U.basis() for w in W.basis())

    return {
        "row_perp_null": all_orthogonal(A.row_space(), A.right_kernel()),
        "col_perp_leftnull": all_orthogonal(A.column_space(), A.left_kernel()),
    }


def row_to_column_isomorphism(A, field=QQ):
    """
    Exhibit the isomorphism ``A : row space → column space`` (the diagram's red arrow).

    Applies ``A`` to a basis of the row space; the images are linearly independent (so the
    restriction is injective) and span the column space (so it is onto) — both of dimension
    ``r = rank(A)``. Returns a dict with the image vectors and ``is_isomorphism`` (bool).
    """
    A = matrix(field, A)
    r = A.rank()
    row_basis = A.row_space().basis()
    images = [A * v for v in row_basis]
    img_matrix = matrix(field, images) if images else matrix(field, 0, A.nrows())

    independent = img_matrix.rank() == r
    spans_column = img_matrix.row_space() == A.column_space()
    return {
        "row_basis": list(row_basis),
        "images": images,
        "is_isomorphism": bool(independent and spans_column),
    }


if __name__ == "__main__":
    print("=" * 70)
    print("The four fundamental subspaces — one matrix, four subspaces")
    print("=" * 70)

    # A 3×4 matrix of rank 2 (rows 1 and 2 are independent; row 3 = row1 + row2).
    A = [[1, 2, 0, 1],
         [0, 1, 1, 0],
         [1, 3, 1, 1]]
    info = four_subspaces(A)
    print("\nA =")
    print(matrix(QQ, A))
    print(f"\nrank r = {info['rank']}   (domain dim n = {info['n']}, codomain dim m = {info['m']})")
    d = info["dims"]
    print("dimensions:")
    print(f"  row space     = {d['row']:>2}      null space      = {d['null']:>2}   "
          f"(row {d['row']} + null {d['null']} = n = {info['n']})")
    print(f"  column space  = {d['column']:>2}      left null space = {d['left_null']:>2}   "
          f"(col {d['column']} + leftnull {d['left_null']} = m = {info['m']})")
    print("rank–nullity holds:", info["rank_nullity_holds"])

    orth = check_orthogonality(A)
    print("\northogonality:")
    print("  row space ⊥ null space        :", orth["row_perp_null"])
    print("  column space ⊥ left null space:", orth["col_perp_leftnull"])

    iso = row_to_column_isomorphism(A)
    print("\nA maps the row space isomorphically onto the column space:")
    for v, Av in zip(iso["row_basis"], iso["images"]):
        print(f"  A·{list(v)} = {list(Av)}")
    print("  is an isomorphism (independent images spanning the column space):",
          iso["is_isomorphism"])
