"""Orbit–Stabilizer, the class equation, and Burnside counting.

A group action ``G ↻ X`` is a homomorphism ``ρ : G → Sym(X)``; it is the precise way a
symmetry group "is understood through how it permutes something concrete" — the same slogan
that runs through :mod:`maths_coding.algebra.galois_group_and_module` and
:mod:`maths_coding.algebra.rotation_group_polyhedra`. Two attached invariants organise
everything:

    Orbit  Orb(x) = { g·x : g ∈ G } ⊆ X          (orbits partition X)
    Stab(x) = { g ∈ G : g·x = x } ≤ G            (a subgroup, generally not normal)

The hinge is the **Orbit–Stabilizer theorem**: for a finite group,

    |Orb(x)| · |Stab(x)| = |G|,

because ``g·Stab(x) ↦ g·x`` is a bijection ``G/Stab(x) ⟷ Orb(x)`` of ``G``-sets (and the
count is then Lagrange's theorem). Summed over orbit representatives it becomes the
**class equation** ``|X| = Σᵢ [G : Stab(xᵢ)]``; specialised to ``G`` acting on itself by
conjugation it is the familiar ``|G| = Σ [G : C(xᵢ)]`` over conjugacy-class reps.

Averaging fixed points instead of summing orbit sizes gives **Burnside's lemma**
(Cauchy–Frobenius):

    #orbits = (1/|G|) Σ_{g∈G} |X^g|,   X^g = { x : g·x = x },

the foundational counting tool of Pólya enumeration. For ``G`` permuting positions and
``X`` the set of ``k``-colourings, ``|X^g| = k^{c(g)}`` where ``c(g)`` is the number of
cycles of ``g`` — so the average is a sum of ``k``-powers, no enumeration needed.

- :func:`orbit`, :func:`stabilizer`, :func:`orbit_stabilizer_check` — the theorem, verified.
- :func:`class_equation` — partition ``X`` into orbits and check the sum.
- :func:`conjugacy_class_equation` — the conjugation action of ``G`` on itself.
- :func:`burnside_orbit_count` — orbits as an average of fixed-point counts (general action).
- :func:`burnside_colorings` / :func:`necklace_count` — Pólya counting of colourings up to a
  permutation group, via cycle counts and the closed form ``(1/n) Σ_{d|n} φ(d) k^{n/d}``.

The demos reproduce the notes' flavour: ``A₄`` (tetrahedron rotations) on its 4 vertices
gives ``4 · 3 = 12``; its conjugacy class sizes are ``1 + 3 + 4 + 4 = 12``; and Burnside
counts the binary necklaces on ``n`` beads.
"""

from sage.all import (
    AlternatingGroup,
    CyclicPermutationGroup,
    DihedralGroup,
    euler_phi,
    divisors,
)
from itertools import product


def orbit(G, x, action=None):
    """
    The orbit ``Orb(x) = { g·x : g ∈ G }`` as a sorted list.

    Parameters
    ----------
    G : Sage group
        Iterated over its elements (small groups only).
    x : object
        A point of the set being acted on.
    action : callable ``(g, x) -> y``, optional
        How ``g`` acts on ``x``. Defaults to the natural permutation action ``g(x)``
        for ``x`` an integer point and ``g`` a Sage permutation.
    """
    if action is None:
        action = lambda g, p: g(p)
    pts = {action(g, x) for g in G}
    return sorted(pts)


def stabilizer(G, x, action=None):
    """
    The stabilizer ``Stab(x) = { g ∈ G : g·x = x }`` as a list of group elements.
    """
    if action is None:
        action = lambda g, p: g(p)
    return [g for g in G if action(g, x) == x]


def orbit_stabilizer_check(G, x, action=None):
    """
    Verify ``|Orb(x)| · |Stab(x)| = |G|`` and return the data.

    Returns
    -------
    dict with ``orbit_size``, ``stab_size``, ``group_order``, ``holds`` (bool).
    """
    orb = orbit(G, x, action)
    stab = stabilizer(G, x, action)
    o, s, n = len(orb), len(stab), G.order()
    return {
        "orbit_size": o,
        "stab_size": s,
        "group_order": int(n),
        "holds": o * s == int(n),
    }


def orbits(G, X, action=None):
    """Partition ``X`` into orbits; returns a list of orbits (each a sorted list)."""
    if action is None:
        action = lambda g, p: g(p)
    remaining = set(X)
    parts = []
    while remaining:
        x = next(iter(remaining))
        orb = set(action(g, x) for g in G)
        parts.append(sorted(orb & set(X)))
        remaining -= orb
    return parts


def class_equation(G, X, action=None):
    """
    The class equation ``|X| = Σ |Orb(xᵢ)|`` for the action of ``G`` on the finite set ``X``.

    Returns
    -------
    dict with ``orbit_sizes`` (list), ``total`` (int), ``holds`` (bool).
    """
    parts = orbits(G, X, action)
    sizes = [len(p) for p in parts]
    return {"orbit_sizes": sizes, "total": sum(sizes), "holds": sum(sizes) == len(set(X))}


def conjugacy_class_equation(G):
    """
    The class equation for ``G`` acting on itself by conjugation: ``|G| = Σ [G : C(xᵢ)]``.

    Here the orbits are the conjugacy classes and the stabilizers are the centralizers.
    Returns the multiset of class sizes and verifies they sum to ``|G|``.
    """
    elements = list(G)
    action = lambda g, h: g * h * g.inverse()
    parts = orbits(G, elements, action)
    sizes = sorted(len(p) for p in parts)
    return {"class_sizes": sizes, "order": int(G.order()), "holds": sum(sizes) == int(G.order())}


def fixed_points(G_element, X, action=None):
    """The fixed-point set ``X^g = { x ∈ X : g·x = x }``."""
    if action is None:
        action = lambda g, p: g(p)
    return [x for x in X if action(G_element, x) == x]


def burnside_orbit_count(G, X, action=None):
    """
    Number of orbits via Burnside: ``(1/|G|) Σ_{g} |X^g|``.

    Computed by averaging fixed-point counts, then cross-checked against the direct orbit
    partition. Returns a dict with ``count``, ``direct`` (from :func:`orbits`), ``agree``.
    """
    if action is None:
        action = lambda g, p: g(p)
    total = sum(len(fixed_points(g, X, action)) for g in G)
    n = int(G.order())
    assert total % n == 0, "Burnside average must be an integer"
    count = total // n
    direct = len(orbits(G, X, action))
    return {"count": count, "direct": direct, "agree": count == direct}


def burnside_colorings(G, points, num_colors):
    """
    Number of ``num_colors``-colourings of ``points`` up to the action of ``G``, via Burnside.

    Uses ``|X^g| = k^{c(g)}`` where ``c(g)`` is the number of cycles of ``g`` (including
    fixed points) on ``points`` — so no colouring enumeration is needed. ``G`` must act on
    the labels ``points`` by permutation.

    Returns
    -------
    int
        The number of distinct colourings.
    """
    k = num_colors
    total = 0
    for g in G:
        cycles = g.cycle_type()  # partition of len(points), one part per cycle
        total += k ** len(cycles)
    n = int(G.order())
    assert total % n == 0
    return total // n


def burnside_colorings_direct(G, points, num_colors):
    """
    Same count as :func:`burnside_colorings`, but by *enumerating* colourings and applying
    the orbit partition directly — a brute-force cross-check for small cases.
    """
    points = list(points)
    pos_index = {p: i for i, p in enumerate(points)}
    X = list(product(range(num_colors), repeat=len(points)))

    def act(g, coloring):
        # g moves the bead at position p to position g(p); recolour accordingly.
        new = [None] * len(points)
        for p in points:
            new[pos_index[g(p)]] = coloring[pos_index[p]]
        return tuple(new)

    return len(orbits(G, X, act))


def necklace_count(n, k):
    """
    Closed form for ``k``-colour necklaces on ``n`` beads (rotations only, cyclic group):

        N(n, k) = (1/n) Σ_{d | n} φ(d) · k^{n/d}.

    Returns
    -------
    int
    """
    total = sum(euler_phi(d) * k ** (n // d) for d in divisors(n))
    assert total % n == 0
    return total // n


if __name__ == "__main__":
    print("=" * 70)
    print("Orbit–Stabilizer, the class equation, and Burnside counting")
    print("=" * 70)

    # ---- Orbit–Stabilizer: A4 (tetrahedron rotations) on its 4 vertices ------
    # The rotation group of the tetrahedron is A4 (order 12) acting on {1,2,3,4}.
    A4 = AlternatingGroup(4)
    os = orbit_stabilizer_check(A4, 1)
    print("\nA4 ↻ {1,2,3,4}  (tetrahedron rotations on vertices):")
    print(f"  |Orb(1)| = {os['orbit_size']}  (transitive)")
    print(f"  |Stab(1)| = {os['stab_size']}  (rotations fixing a vertex)")
    print(f"  |Orb|·|Stab| = {os['orbit_size']*os['stab_size']} = |G| = {os['group_order']}",
          " ✓" if os["holds"] else " ✗")

    # ---- Class equation via conjugation --------------------------------------
    ce = conjugacy_class_equation(A4)
    print("\nConjugation class equation for A4:")
    print(f"  class sizes {ce['class_sizes']} sum to {sum(ce['class_sizes'])} = |A4| = {ce['order']}",
          " ✓" if ce["holds"] else " ✗")

    # ---- Burnside as an average, cross-checked -------------------------------
    # Count 2-colourings of a square's 4 vertices under the dihedral group D4.
    D4 = DihedralGroup(4)
    square = [1, 2, 3, 4]
    bc = burnside_colorings(D4, square, 2)
    bc_direct = burnside_colorings_direct(D4, square, 2)
    print("\n2-colourings of a square (D4 symmetry):")
    print(f"  Burnside (cycle counts) = {bc}   brute-force = {bc_direct}",
          " ✓" if bc == bc_direct else " ✗")

    # ---- Necklaces: rotation-only colourings, two ways -----------------------
    print("\nBinary necklaces on n beads (rotations C_n), closed form vs Burnside:")
    for n in range(3, 9):
        Cn = CyclicPermutationGroup(n)
        formula = necklace_count(n, 2)
        burnside = burnside_colorings(Cn, list(range(1, n + 1)), 2)
        flag = "✓" if formula == burnside else "✗"
        print(f"  n={n}:  (1/n)Σφ(d)2^(n/d) = {formula:>4}   Burnside = {burnside:>4}  {flag}")
