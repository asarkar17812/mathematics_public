"""Hilbert's Nullstellensatz: the dictionary between geometry and commutative algebra.

Affine algebraic geometry runs on a single contravariant correspondence between *geometry*
(zero loci) and *algebra* (ideals of a polynomial ring). For a field ``k`` it is built from
two maps:

    V(S) = { p ∈ 𝔸ⁿ : f(p) = 0  ∀ f ∈ S }      (geometry from polynomials)
    I(X) = { f : f(p) = 0  ∀ p ∈ X }            (the radical ideal of a set)

Over an **algebraically closed** ``k`` the **Nullstellensatz** makes these mutually inverse,
in three equivalent forms:

    1. Weak:      every maximal ideal is (x₁−a₁, …, xₙ−aₙ)  ⇔  MaxSpec ≅ 𝔸ⁿ.
    2. Solutions: a *proper* ideal I has V(I) ≠ ∅.
    3. Strong:    I(V(I)) = √I.

So the strong form says a polynomial vanishes on V(I) **iff** some power lies in ``I`` —
radical membership, which Sage decides with the Rabinowitsch trick: ``f ∈ √I`` iff
``1 ∈ (I, 1 − t·f)`` in one extra variable ``t``. The resulting bijections are
inclusion-reversing: bigger varieties ↔ smaller ideals, **irreducible** variety ↔ **prime**
ideal, a single **point** ↔ a **maximal** ideal. Decomposing an algebraic set into
irreducible components is exactly the minimal-prime decomposition of ``√I``.

The same Jacobian that governs differential geometry returns here as the **Zariski tangent
space** ``T_pV = ker J_p``: a point is *smooth* when ``dim T_pV = dim V``, and the singular
locus ``V(f, ∂f/∂x, ∂f/∂y)`` is where that fails — the **node** of ``y² = x²(x+1)`` and the
**cusp** of ``y² = x³`` against the **smooth** parabola ``y = x²``.

- :func:`vanishing_locus` — the points of ``V(I)`` (over an algebraically closed field).
- :func:`vanishing_ideal` — ``I(X)`` for a finite set of points.
- :func:`radical_membership` — strong Nullstellensatz via the Rabinowitsch trick.
- :func:`nullstellensatz_identity` — check ``I(V(I)) = √I`` on a rational example.
- :func:`singular_points` / :func:`is_smooth_curve` — the Jacobian smoothness criterion.
- :func:`irreducible_components` — ``V(I) = ⋃ V(Pᵢ)`` via minimal primes.

The Nullstellensatz *fails* over non-closed fields: ``(x²+1) ⊆ ℝ[x]`` is proper yet cuts out
``∅`` in ℝ — the demo shows its two points reappearing over ``\\overline{ℚ}``.
"""

from sage.all import QQ, QQbar, PolynomialRing


def _ring(varnames, base=QQ):
    """Build a polynomial ring and return ``(R, gens)``."""
    R = PolynomialRing(base, list(varnames))
    return R, R.gens()


def vanishing_locus(ideal_gens, field=QQbar):
    """
    The points of ``V(I) = { p : f(p) = 0 ∀ f ∈ I }`` over ``field``.

    Requires ``I`` to be zero-dimensional (finitely many points); over ``QQbar`` this
    realises the *weak Nullstellensatz* — a proper ideal has a nonempty variety. Returns a
    list of coordinate tuples (in the order of the ring's generators).
    """
    I = ideal_gens[0].parent().ideal(ideal_gens)
    gens = I.ring().gens()
    if I.dimension() > 0:
        raise ValueError("V(I) is positive-dimensional; not a finite point set")
    return [tuple(d[g] for g in gens) for d in I.variety(field)]


def vanishing_ideal(points, R):
    """
    The ideal ``I(X) = { f : f(p) = 0 ∀ p ∈ X }`` of a finite set of points.

    Each point ``p = (a₁, …, aₙ)`` contributes the maximal ideal ``(x₁−a₁, …, xₙ−aₙ)``, and
    ``I(X)`` is their intersection — the algebraic incarnation of "a point ↔ a maximal
    ideal", with ``I(∅)`` the whole ring.
    """
    gens = R.gens()
    result = None
    for p in points:
        m = R.ideal([gens[i] - R(p[i]) for i in range(len(gens))])
        result = m if result is None else result.intersection(m)
    return R.ideal([R.one()]) if result is None else result


def radical_membership(ideal_gens, f):
    """
    Decide ``f ∈ √I`` (strong Nullstellensatz) by the **Rabinowitsch trick**.

    ``f ∈ √I`` iff ``1 ∈ (g₁, …, g_r, 1 − t·f)`` in the polynomial ring with one extra
    variable ``t``: a power of ``f`` lies in ``I`` exactly when adjoining an inverse of ``f``
    collapses the ring. The result is cross-checked against Sage's ``f in I.radical()``.

    Returns
    -------
    dict with ``in_radical`` (Rabinowitsch), ``in_ideal`` (plain membership),
    ``cross_check`` (via ``I.radical()``), ``agree`` (bool).
    """
    R = f.parent()
    n = R.ngens()
    base = R.base_ring()
    S = PolynomialRing(base, list(R.variable_names()) + ["t_rabinowitsch"])
    embed = R.hom(list(S.gens()[:n]), S)
    t = S.gens()[-1]

    big = S.ideal([embed(g) for g in ideal_gens] + [1 - t * embed(f)])
    rab = S.one() in big

    I = R.ideal(ideal_gens)
    return {
        "in_radical": rab,
        "in_ideal": f in I,
        "cross_check": f in I.radical(),
        "agree": rab == (f in I.radical()),
    }


def nullstellensatz_identity(ideal_gens, field=QQ):
    """
    Verify the strong Nullstellensatz ``I(V(I)) = √I`` on a zero-dimensional example.

    Computes ``V(I)`` over ``field``, forms ``I(V(I))``, and compares it to ``√I``. For the
    identity to hold the points must be taken over an algebraically closed field; the
    default ``field=QQ`` is correct only when ``V(I)`` is entirely rational (as in the
    examples below), and is used so the recovered ideal stays over ℚ.

    Returns
    -------
    dict with ``points``, ``radical`` (gens), ``recovered`` (gens of I(V(I))), ``equal``.
    """
    R = ideal_gens[0].parent()
    I = R.ideal(ideal_gens)
    rad = I.radical()
    pts = vanishing_locus(ideal_gens, field=field)
    pts = [tuple(field(c) for c in p) for p in pts]
    IV = vanishing_ideal(pts, R)
    return {
        "points": pts,
        "radical": list(rad.gens()),
        "recovered": list(IV.gens()),
        "equal": IV == rad,
    }


def singular_points(f, field=QQbar):
    """
    Singular points of the plane curve ``V(f)``: where ``f`` and all its partials vanish.

    The Zariski tangent space is ``ker J_p``; the curve is singular at ``p`` exactly when the
    Jacobian drops rank, i.e. ``p ∈ V(f, ∂f/∂x, ∂f/∂y)``. Returns the finite list of such
    points (empty for a smooth curve).
    """
    R = f.parent()
    gens = R.gens()
    J = R.ideal([f] + [f.derivative(v) for v in gens])
    if J.dimension() <= 0:
        return [tuple(d[g] for g in gens) for d in J.variety(field)]
    raise ValueError("singular locus is positive-dimensional")


def is_smooth_curve(f):
    """Whether the plane curve ``V(f)`` is smooth (no singular points over ``\\overline{ℚ}``)."""
    return len(singular_points(f)) == 0


def irreducible_components(ideal_gens):
    """
    The irreducible components of ``V(I)`` as the **minimal primes** over ``I``.

    By the dictionary, ``V(I) = V(P₁) ∪ ⋯ ∪ V(P_r)`` where the ``Pᵢ`` are the minimal primes
    containing ``I`` (equivalently the minimal primary decomposition of ``√I``); each prime
    is an irreducible piece. E.g. ``V(xy)`` splits into the two coordinate axes ``V(x)`` and
    ``V(y)``.

    Returns
    -------
    list of prime ideals.
    """
    I = ideal_gens[0].parent().ideal(ideal_gens)
    return I.minimal_associated_primes()


if __name__ == "__main__":
    print("=" * 70)
    print("Hilbert's Nullstellensatz — geometry ⟷ commutative algebra")
    print("=" * 70)

    # ---- Strong Nullstellensatz: radical membership --------------------------
    Rx, (x1,) = _ring(["x"])
    print("\nStrong Nullstellensatz (radical membership), I = (x²) ⊆ ℚ[x]:")
    rm = radical_membership([x1 ** 2], x1)
    print(f"  x ∈ (x²)?  {rm['in_ideal']}    x ∈ √(x²)?  {rm['in_radical']}"
          f"   (Rabinowitsch agrees with I.radical(): {rm['agree']})")
    print("  → x vanishes on V(x²) = {0} but is not in (x²); it lies in √(x²) = (x).")

    # ---- I(V(I)) = √I on a rational example ----------------------------------
    R2, (x, y) = _ring(["x", "y"])
    print("\nI(V(I)) = √I  for I = ((x−1)², y):")
    idn = nullstellensatz_identity([(x - 1) ** 2, y])
    print(f"  V(I) = {idn['points']}")
    print(f"  √I        = {idn['radical']}")
    print(f"  I(V(I))   = {idn['recovered']}")
    print(f"  equal? {idn['equal']}   (a fat point's radical is the reduced point (x−1, y))")

    # ---- Smoothness via the Jacobian: parabola, node, cusp -------------------
    print("\nSmoothness via the Zariski tangent space (singular = J drops rank):")
    parabola = y - x ** 2
    node = y ** 2 - x ** 2 * (x + 1)
    cusp = y ** 2 - x ** 3
    for name, f in [("parabola  y − x²", parabola),
                    ("node      y² − x²(x+1)", node),
                    ("cusp      y² − x³", cusp)]:
        sing = singular_points(f)
        tag = "smooth" if not sing else f"singular at {sing}"
        print(f"  {name:<26} → {tag}")

    # ---- Irreducible decomposition: V(xy) = two axes -------------------------
    print("\nIrreducible components of V(xy)  (reducible f ↔ reducible variety):")
    for P in irreducible_components([x * y]):
        print(f"  prime {P.gens()}  ↔  one coordinate axis")

    # ---- Failure over a non-closed field -------------------------------------
    print("\nFailure of the Nullstellensatz over ℝ:  I = (x² + 1) ⊆ ℚ[x] is proper,")
    print("  yet V(I) has no real/rational point — the two points appear over ℚ̄:")
    Rc, (xc,) = _ring(["x"])
    pts_closed = vanishing_locus([xc ** 2 + 1], field=QQbar)
    print(f"  V(x²+1) over ℚ̄ = {pts_closed}   (i.e. ±i — empty over ℝ)")
