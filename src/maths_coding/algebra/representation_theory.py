"""Representation theory of finite groups — symmetry made into matrices, read off by characters.

A **representation** is a homomorphism ``ρ : G → GL(V)`` — a way to realize the abstract group
as concrete matrices acting on a vector space. Over ``ℂ`` every representation is a direct sum
of **irreducible** ones (Maschke), and the irreducibles are classified by a single, remarkably
small invariant: the **character** ``χ_V(g) = tr ρ(g)``. Characters are **class functions**
(constant on conjugacy classes, since trace is conjugation-invariant), and two representations
are isomorphic iff their characters agree — so the entire theory is encoded in a finite
**character table** whose rows are the irreducible characters and whose columns are the
conjugacy classes.

The structure is governed by **Schur orthogonality**. With the Hermitian inner product

    ⟨φ, ψ⟩ = (1/|G|) Σ_g \\overline{φ(g)} ψ(g) = (1/|G|) Σ_classes |C| \\overline{φ(C)} ψ(C),

the irreducible characters form an **orthonormal basis** of the class functions. Two corollaries
fall out immediately and are checked here:

    #irreducibles = #conjugacy classes      (the table is square),
    Σ dᵢ² = |G|                              (degrees from the identity column).

Decomposing any representation is then just taking inner products: the multiplicity of the
irreducible ``Vᵢ`` in ``W`` is ``⟨χ_W, χᵢ⟩``. The **regular representation** (``G`` acting on
``ℂ[G]``) decomposes with each ``Vᵢ`` appearing exactly ``dᵢ`` times — the cleanest proof of
``Σ dᵢ² = |G|``.

- :func:`character_table`, :func:`conjugacy_class_sizes` — the table and its column weights.
- :func:`irreducible_degrees` — the degrees ``dᵢ`` and the check ``Σ dᵢ² = |G|``.
- :func:`count_check` — ``#irreducibles = #conjugacy classes``.
- :func:`character_inner_product` — the class-function inner product.
- :func:`row_orthogonality`, :func:`column_orthogonality` — the two Schur relations.
- :func:`decompose` — multiplicities of a class function in the irreducible basis.
- :func:`regular_representation_decomposition` — each ``Vᵢ`` appears ``dᵢ`` times.
"""

from sage.all import ZZ


def character_table(G):
    """The character table of ``G`` (rows = irreducible characters, columns = conjugacy classes)."""
    return G.character_table()


def conjugacy_class_sizes(G):
    """The sizes ``|C|`` of the conjugacy classes, in the column order of the character table."""
    return [int(c.cardinality()) for c in G.conjugacy_classes()]


def irreducible_degrees(G):
    """
    The degrees ``dᵢ = χᵢ(e)`` (the identity column of the table), with the sum-of-squares check.

    Returns
    -------
    dict with ``degrees`` (list), ``sum_of_squares``, ``group_order``, and ``holds``
    (``Σ dᵢ² = |G|``).
    """
    T = character_table(G)
    degrees = [int(T[i, 0]) for i in range(T.nrows())]   # column 0 is the identity class
    ssq = sum(d * d for d in degrees)
    return {
        "degrees": degrees,
        "sum_of_squares": ssq,
        "group_order": int(G.order()),
        "holds": ssq == int(G.order()),
    }


def count_check(G):
    """Verify ``#irreducible characters = #conjugacy classes`` (the table is square)."""
    T = character_table(G)
    n_irr = T.nrows()
    n_cls = len(G.conjugacy_classes())
    return {"num_irreducibles": int(n_irr), "num_classes": int(n_cls), "equal": n_irr == n_cls}


def character_inner_product(values1, values2, class_sizes, group_order):
    """
    The class-function inner product ``⟨φ, ψ⟩ = (1/|G|) Σ_classes |C| \\overline{φ(C)} ψ(C)``.

    ``values1``/``values2`` are the values of ``φ``/``ψ`` on the conjugacy classes (table-column
    order). Returns the inner product (an element of the cyclotomic field, equal to a rational
    for genuine characters).
    """
    total = sum(s * v1.conjugate() * v2 for s, v1, v2 in zip(class_sizes, values1, values2))
    return total / group_order


def row_orthogonality(G):
    """
    The first Schur relation: the irreducible characters are **orthonormal**,
    ``⟨χᵢ, χⱼ⟩ = δᵢⱼ``.

    Returns
    -------
    dict with ``gram`` (the matrix of inner products) and ``is_identity`` (bool).
    """
    T = character_table(G)
    sizes = conjugacy_class_sizes(G)
    order = int(G.order())
    n = T.nrows()
    gram = [[character_inner_product(T[i], T[j], sizes, order) for j in range(n)]
            for i in range(n)]
    is_identity = all(gram[i][j] == (1 if i == j else 0) for i in range(n) for j in range(n))
    return {"gram": gram, "is_identity": is_identity}


def column_orthogonality(G):
    """
    The second Schur relation (columns): ``Σᵢ \\overline{χᵢ(g)} χᵢ(h) = |C_G(g)|`` if ``g, h`` are
    conjugate and ``0`` otherwise.

    The diagonal value for a class ``C`` is the centralizer order ``|G| / |C|``. Returns a dict
    with ``holds`` (all column inner products match the expected centralizer pattern).
    """
    T = character_table(G)
    sizes = conjugacy_class_sizes(G)
    order = int(G.order())
    ncols = T.ncols()
    ok = True
    for a in range(ncols):
        for b in range(ncols):
            col_ip = sum(T[i, a].conjugate() * T[i, b] for i in range(T.nrows()))
            expected = (order // sizes[a]) if a == b else 0    # |C_G(g)| = |G|/|class|
            if col_ip != expected:
                ok = False
    return {"holds": ok}


def decompose(G, class_values):
    """
    Decompose a class function (e.g. the character of some representation) into irreducibles.

    The multiplicity of the ``i``-th irreducible is ``mᵢ = ⟨χ, χᵢ⟩``; the representation is a
    genuine one exactly when all ``mᵢ`` are non-negative integers. ``class_values`` gives the
    character's values in table-column order.

    Returns
    -------
    dict with ``multiplicities`` (list) and ``reconstructs`` (``χ = Σ mᵢ χᵢ`` on every class).
    """
    T = character_table(G)
    sizes = conjugacy_class_sizes(G)
    order = int(G.order())
    n = T.nrows()
    mults = [character_inner_product(T[i], class_values, sizes, order) for i in range(n)]
    reconstructs = all(
        sum(mults[i] * T[i, c] for i in range(n)) == class_values[c]
        for c in range(T.ncols())
    )
    return {"multiplicities": mults, "reconstructs": reconstructs}


def regular_representation_decomposition(G):
    """
    Decompose the **regular representation** (``G`` on ``ℂ[G]``): its character is ``|G|`` on the
    identity and ``0`` elsewhere, and it contains each irreducible ``Vᵢ`` with multiplicity
    ``dᵢ`` — the conceptual source of ``Σ dᵢ² = |G|``.

    Returns
    -------
    dict with ``multiplicities`` and ``degrees`` and ``equal`` (multiplicity = degree for each).
    """
    T = character_table(G)
    ncols = T.ncols()
    # χ_reg = (|G|, 0, 0, …): identity class first, value |G|; all other classes 0.
    reg = [int(G.order()) if c == 0 else 0 for c in range(ncols)]
    dec = decompose(G, reg)
    degrees = irreducible_degrees(G)["degrees"]
    equal = [int(m) for m in dec["multiplicities"]] == degrees
    return {"multiplicities": dec["multiplicities"], "degrees": degrees, "equal": equal}


if __name__ == "__main__":
    from sage.all import SymmetricGroup, CyclicPermutationGroup, QuaternionGroup

    print("=" * 70)
    print("Representation theory — characters and Schur orthogonality")
    print("=" * 70)

    for name, G in [("S3", SymmetricGroup(3)),
                    ("C4", CyclicPermutationGroup(4)),
                    ("Q8", QuaternionGroup())]:
        print(f"\n{name}  (order {G.order()}):")
        print("  character table (rows = irreducibles, cols = classes):")
        print(character_table(G))
        deg = irreducible_degrees(G)
        print(f"  degrees {deg['degrees']}   Σdᵢ² = {deg['sum_of_squares']} = |G| = "
              f"{deg['group_order']}  ({deg['holds']})")
        cc = count_check(G)
        print(f"  #irreducibles = {cc['num_irreducibles']} = #conjugacy classes = "
              f"{cc['num_classes']}  ({cc['equal']})")
        print("  characters orthonormal (row orthogonality):", row_orthogonality(G)["is_identity"])
        print("  column orthogonality holds:", column_orthogonality(G)["holds"])
        reg = regular_representation_decomposition(G)
        print(f"  regular rep: each Vᵢ appears dᵢ times → mult {[int(m) for m in reg['multiplicities']]}"
              f" = degrees {reg['degrees']}  ({reg['equal']})")
