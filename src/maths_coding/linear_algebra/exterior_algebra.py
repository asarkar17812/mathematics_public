"""Tensor and exterior algebra — multilinear maps, and the determinant as a volume.

Multilinear algebra is built from two constructions on a vector space ``V`` (``dim V = n``):

**Tensor product.** ``V ⊗ W`` is the home of bilinear maps; ``dim(V ⊗ W) = dim V · dim W``.

**Exterior power.** ``Λᵏ V`` is the *antisymmetric* part — alternating multilinear forms — of
dimension ``\\binom{n}{k}`` (Pascal's triangle). Its elements are wedges ``v₁ ∧ ⋯ ∧ vₖ``, which
are **antisymmetric** (swapping two factors flips the sign) and therefore **vanish on
repeats** (``v ∧ v = 0``): a wedge is nonzero exactly when its vectors are linearly
independent. The components of ``v₁ ∧ ⋯ ∧ vₖ`` in the basis ``e_{i₁} ∧ ⋯ ∧ e_{iₖ}`` are the
``k×k`` **minors** (the Plücker coordinates) of the matrix whose rows are the ``vᵢ``.

The punchline is the top power. ``Λⁿ V`` is **one-dimensional**, so ``v₁ ∧ ⋯ ∧ vₙ`` is a single
number times ``e₁ ∧ ⋯ ∧ eₙ`` — and that number is the **determinant**. This is exactly why
``det`` is the signed volume scaling factor of a map: it is the action on the top exterior
power, the oriented ``n``-volume of the parallelepiped on the rows (compare the Gram volume in
:mod:`maths_coding.linear_algebra.forms`). In ``ℝ³`` the three components of ``v ∧ w`` are the
**cross product** (up to the Hodge star), the bridge from this algebra to the differential
forms of geometry.

- :func:`tensor_dimension`, :func:`exterior_power_dimension` — the ``∏ dim`` and ``\\binom{n}{k}`` counts.
- :func:`wedge` — the Plücker coordinates (``k×k`` minors) of ``v₁ ∧ ⋯ ∧ vₖ``.
- :func:`is_decomposable_zero` — whether a wedge vanishes (linear dependence).
- :func:`antisymmetry_check` — ``v ∧ w = −(w ∧ v)`` and ``v ∧ v = 0``.
- :func:`determinant_as_top_wedge` — ``v₁ ∧ ⋯ ∧ vₙ = det · (e₁ ∧ ⋯ ∧ eₙ)``.
- :func:`cross_product_via_wedge` — the ``ℝ³`` wedge ↔ cross product (Hodge star) correspondence.
"""

from itertools import combinations
from sage.all import matrix, QQ, binomial, prod, vector


def tensor_dimension(*dims):
    """``dim(V₁ ⊗ ⋯ ⊗ V_r) = ∏ dim Vᵢ`` — the dimension of the tensor product."""
    return int(prod(dims))


def exterior_power_dimension(n, k):
    """``dim Λᵏ V = \\binom{n}{k}`` for ``0 ≤ k ≤ n`` (and ``0`` otherwise) — note ``Λⁿ`` is 1-dimensional."""
    if k < 0 or k > n:
        return 0
    return int(binomial(n, k))


def wedge(vectors, field=QQ):
    """
    The wedge ``v₁ ∧ ⋯ ∧ vₖ`` as its Plücker coordinates: the ``k×k`` minors of the matrix whose
    rows are the ``vᵢ``.

    Returns
    -------
    list of ``(index_tuple, value)``
        ``index_tuple = (i₁ < ⋯ < iₖ)`` labels the basis element ``e_{i₁} ∧ ⋯ ∧ e_{iₖ}`` and
        ``value`` is the corresponding minor. All-zero ⇔ the ``vᵢ`` are linearly dependent
        (in particular ``k > n`` gives the empty/zero wedge).
    """
    M = matrix(field, [list(v) for v in vectors])
    k, n = M.nrows(), M.ncols()
    if k > n:
        return []
    minors = M.minors(k)
    labels = list(combinations(range(n), k))
    return list(zip(labels, minors))


def is_decomposable_zero(vectors, field=QQ):
    """Whether ``v₁ ∧ ⋯ ∧ vₖ = 0`` — equivalently the vectors are linearly dependent."""
    return all(val == 0 for _, val in wedge(vectors, field))


def antisymmetry_check(v, w, field=QQ):
    """
    Verify the defining antisymmetry: ``v ∧ w = −(w ∧ v)`` and ``v ∧ v = 0``.

    Returns
    -------
    dict with ``swap_negates`` (bool) and ``repeat_is_zero`` (bool).
    """
    vw = [val for _, val in wedge([v, w], field)]
    wv = [val for _, val in wedge([w, v], field)]
    swap_negates = all(a == -b for a, b in zip(vw, wv))
    repeat_is_zero = is_decomposable_zero([v, v], field)
    return {"swap_negates": swap_negates, "repeat_is_zero": repeat_is_zero}


def determinant_as_top_wedge(A, field=QQ):
    """
    Confirm ``v₁ ∧ ⋯ ∧ vₙ = det(A) · (e₁ ∧ ⋯ ∧ eₙ)`` for the rows ``vᵢ`` of a square ``A``.

    The top exterior power ``Λⁿ`` is one-dimensional, so the wedge of the ``n`` rows is a single
    Plücker coordinate — and it equals ``det A``, the oriented volume of the parallelepiped they
    span.

    Returns
    -------
    dict with ``top_wedge`` (the single coordinate), ``det`` (``A.det()``), ``equal``.
    """
    M = matrix(field, A)
    coords = wedge(M.rows(), field)
    top = coords[0][1] if coords else field(0)
    d = M.det()
    return {"top_wedge": top, "det": d, "equal": top == d}


def cross_product_via_wedge(v, w, field=QQ):
    """
    The ``ℝ³`` correspondence ``v ∧ w ↔ v × w`` (the Hodge star).

    The three Plücker coordinates of ``v ∧ w`` (the minors over columns ``{0,1}, {0,2}, {1,2}``)
    rearrange — with one sign — into the cross product ``v × w``. Returns a dict with
    ``plucker`` (the wedge coordinates), ``cross_via_wedge`` (Hodge-star rearrangement), and
    ``cross_builtin`` (Sage's ``cross_product``), and ``equal``.
    """
    coords = wedge([v, w], field)               # [(0,1): m01, (0,2): m02, (1,2): m12]
    m01, m02, m12 = (val for _, val in coords)
    cross_via_wedge = vector(field, [m12, -m02, m01])
    cross_builtin = vector(field, v).cross_product(vector(field, w))
    return {
        "plucker": coords,
        "cross_via_wedge": cross_via_wedge,
        "cross_builtin": cross_builtin,
        "equal": cross_via_wedge == cross_builtin,
    }


if __name__ == "__main__":
    print("=" * 70)
    print("Exterior algebra — antisymmetry, and the determinant as a volume")
    print("=" * 70)

    print("\nDimensions:")
    print("  dim(V⊗W) for dim V = 2, dim W = 3 :", tensor_dimension(2, 3))
    print("  dim Λᵏ(R⁴) for k = 0..4           :",
          [exterior_power_dimension(4, k) for k in range(5)], " (Pascal row; Λ⁴ is 1-dimensional)")

    print("\nA wedge as Plücker coordinates (the k×k minors):")
    for label, val in wedge([[1, 2, 0], [0, 1, 3]]):
        names = "∧".join(f"e{i+1}" for i in label)
        print(f"  coefficient of {names:<7} = {val}")

    print("\nAntisymmetry of the wedge:")
    chk = antisymmetry_check([1, 2, 3], [4, 5, 6])
    print("  v ∧ w = −(w ∧ v):", chk["swap_negates"], "    v ∧ v = 0:", chk["repeat_is_zero"])
    print("  linearly dependent vectors wedge to 0:",
          is_decomposable_zero([[1, 2, 3], [2, 4, 6]]), " (second = 2·first)")

    print("\nThe top wedge is the determinant (Λⁿ is one-dimensional):")
    dt = determinant_as_top_wedge([[1, 2], [3, 4]])
    print(f"  rows (1,2),(3,4):  v₁∧v₂ = {dt['top_wedge']} = det = {dt['det']}   ({dt['equal']})")

    print("\nIn R³, v ∧ w is the cross product (up to the Hodge star):")
    cp = cross_product_via_wedge([1, 2, 3], [4, 5, 6])
    print(f"  Plücker minors → {list(cp['cross_via_wedge'])}   cross_product → {list(cp['cross_builtin'])}"
          f"   match: {cp['equal']}")
