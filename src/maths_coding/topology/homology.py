"""Simplicial homology — counting holes, with torsion, by linear algebra over ℤ.

A chain complex is a sequence of free abelian groups joined by **boundary maps**

    ⋯ → C_{k+1} --∂_{k+1}--> C_k --∂_k--> C_{k-1} → ⋯,   with   ∂_k ∂_{k+1} = 0,

and its **homology** ``H_k = ker ∂_k / im ∂_{k+1}`` measures the ``k``-dimensional holes: cycles
(things with no boundary) modulo boundaries (things that bound). The free rank of ``H_k`` is the
**Betti number** ``b_k`` (the count of ``k``-holes), and any finite cyclic summands are
**torsion** (a ``ℤ/d`` records a hole you must go around ``d`` times to bound).

The whole computation is one tool from this package, applied twice: the **Smith normal form**
over ``ℤ`` (cf. :mod:`maths_coding.algebra.smith_normal_form`). For each boundary map its rank
and its invariant factors ``d₁ ∣ d₂ ∣ ⋯`` fall out of the SNF, and then

    b_k = (#k-cells) − rank ∂_k − rank ∂_{k+1},      torsion of H_k = { dᵢ > 1 of ∂_{k+1} }.

The closed surfaces tell the story cleanly: the sphere ``S²`` has homology ``ℤ, 0, ℤ``; the
torus ``ℤ, ℤ², ℤ``; the projective plane ``ℝP²`` has **torsion** ``ℤ, ℤ/2, 0``; and the Klein
bottle ``ℤ, ℤ⊕ℤ/2, 0`` — non-orientability shows up precisely as a ``ℤ/2``. The alternating sum
of Betti numbers is the **Euler characteristic** ``χ = Σ (−1)^k b_k = Σ (−1)^k (#k\\text{-cells})``.

Pure-Python and exact (a from-scratch integer Smith normal form), so it runs without SageMath:

- :func:`smith_invariant_factors` — the nonzero invariant factors of an integer matrix.
- :func:`homology` — Betti numbers and torsion of a chain complex from its boundary maps.
- :func:`euler_characteristic` — ``Σ (−1)^k (#k\\text{-cells})``.
- :data:`SURFACES` — ready-made chain complexes for ``S²``, ``T²``, ``ℝP²``, the Klein bottle.
"""


def smith_invariant_factors(matrix):
    """
    The nonzero **invariant factors** ``d₁ ∣ d₂ ∣ ⋯`` of an integer matrix (its Smith normal
    form diagonal), computed by integer row/column reduction with gcd pivoting.

    ``matrix`` is a list of rows. A row count with empty rows (a ``p×0`` matrix) or no rows
    (``0×q``) has rank ``0`` and returns ``[]``. The length of the result is ``rank(matrix)``.
    """
    A = [row[:] for row in matrix]
    m = len(A)
    n = len(A[0]) if m else 0
    if m == 0 or n == 0:
        return []

    def swap_rows(i, j):
        A[i], A[j] = A[j], A[i]

    def swap_cols(i, j):
        for r in A:
            r[i], r[j] = r[j], r[i]

    factors = []
    t = 0
    while t < min(m, n):
        # Choose the smallest-magnitude nonzero entry in the active submatrix as pivot.
        piv = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i][j] != 0 and (piv is None or abs(A[i][j]) < abs(A[piv[0]][piv[1]])):
                    piv = (i, j)
        if piv is None:
            break
        swap_rows(t, piv[0])
        swap_cols(t, piv[1])

        # Clear the pivot row and column; gcd pivoting guarantees termination.
        cleared = False
        while not cleared:
            cleared = True
            for i in range(t + 1, m):
                if A[i][t] != 0:
                    q = A[i][t] // A[t][t]
                    for j in range(t, n):
                        A[i][j] -= q * A[t][j]
                    if A[i][t] != 0:
                        swap_rows(t, i)
                        cleared = False
            for j in range(t + 1, n):
                if A[t][j] != 0:
                    q = A[t][j] // A[t][t]
                    for i in range(t, m):
                        A[i][j] -= q * A[i][t]
                    if A[t][j] != 0:
                        swap_cols(t, j)
                        cleared = False

        # Enforce the divisibility d_t | (everything below-right); else fold a row in and redo.
        repaired = False
        for i in range(t + 1, m):
            for j in range(t + 1, n):
                if A[i][j] % A[t][t] != 0:
                    for jj in range(n):
                        A[t][jj] += A[i][jj]
                    repaired = True
                    break
            if repaired:
                break
        if repaired:
            continue

        factors.append(abs(A[t][t]))
        t += 1

    return factors


def _rank(matrix):
    """Rank of an integer matrix = number of nonzero invariant factors."""
    return len(smith_invariant_factors(matrix))


def homology(dims, boundaries):
    """
    Homology of a chain complex: Betti numbers and torsion from the boundary maps.

    Parameters
    ----------
    dims : dict ``{k: number of k-cells}``.
    boundaries : dict ``{k: matrix of ∂_k}``, each a ``dims[k-1] × dims[k]`` integer matrix
        (a ``p×0`` map is ``[[] ]*p``; a ``0×q`` map is ``[]``). ``∂_0`` is omitted (zero).

    Returns
    -------
    dict ``{k: {"betti": b_k, "torsion": [d>1 …], "description": "ℤ^b ⊕ ℤ/d ⊕ …"}}``.
    """
    ks = sorted(dims)
    ranks = {k: _rank(boundaries[k]) for k in boundaries}

    def rank_d(k):
        return ranks.get(k, 0)

    out = {}
    for k in ks:
        betti = dims[k] - rank_d(k) - rank_d(k + 1)
        torsion = [d for d in smith_invariant_factors(boundaries.get(k + 1, [])) if d > 1]
        pieces = []
        if betti == 1:
            pieces.append("Z")
        elif betti > 1:
            pieces.append(f"Z^{betti}")
        pieces.extend(f"Z/{d}" for d in torsion)
        out[k] = {
            "betti": betti,
            "torsion": torsion,
            "description": " ⊕ ".join(pieces) if pieces else "0",
        }
    return out


def euler_characteristic(dims):
    """The Euler characteristic ``χ = Σ (−1)^k (#k-cells)`` (equals ``Σ (−1)^k b_k``)."""
    return sum((-1) ** k * dims[k] for k in dims)


# Minimal CW chain complexes of the standard closed surfaces (one 0-cell each).
# Boundary maps are ∂_1 : C_1 → C_0 and ∂_2 : C_2 → C_1, sized dims[k-1] × dims[k].
SURFACES = {
    # Sphere: 1 vertex, 0 edges, 1 face. ∂_1 is 1×0, ∂_2 is 0×1.
    "S^2 (sphere)": (
        {0: 1, 1: 0, 2: 1},
        {1: [[]], 2: []},
    ),
    # Torus: 1 vertex, 2 edges a,b, 1 face a b a⁻¹ b⁻¹ (boundary 0 in each edge).
    "T^2 (torus)": (
        {0: 1, 1: 2, 2: 1},
        {1: [[0, 0]], 2: [[0], [0]]},
    ),
    # Projective plane: 1 vertex, 1 edge a, 1 face a a (boundary 2·a) — torsion appears.
    "RP^2 (proj. plane)": (
        {0: 1, 1: 1, 2: 1},
        {1: [[0]], 2: [[2]]},
    ),
    # Klein bottle: 1 vertex, 2 edges a,b, 1 face a b a b⁻¹ (boundary 2·a + 0·b).
    "Klein bottle": (
        {0: 1, 1: 2, 2: 1},
        {1: [[0, 0]], 2: [[2], [0]]},
    ),
}


if __name__ == "__main__":
    print("=" * 70)
    print("Simplicial homology — Betti numbers and torsion via Smith normal form")
    print("=" * 70)

    print("\nInteger Smith normal form (invariant factors):")
    for M in ([[2, 4], [2, 6]], [[2, 0], [0, 6], [0, 0]], [[2]]):
        print(f"  {M}  →  invariant factors {smith_invariant_factors(M)}")

    print("\nHomology of the standard closed surfaces:")
    for name, (dims, bd) in SURFACES.items():
        H = homology(dims, bd)
        groups = "   ".join(f"H_{k} = {H[k]['description']}" for k in sorted(dims))
        chi = euler_characteristic(dims)
        print(f"  {name:<20}: {groups:<34}  χ = {chi}")

    print("\n  (RP² and the Klein bottle carry ℤ/2 torsion — the algebraic trace of non-orientability;")
    print("   χ: sphere 2, torus 0, ℝP² 1, Klein bottle 0.)")
