"""Localization: inverting elements to zoom in on a ring.

Localization is the algebraic analogue of restricting attention to a neighborhood — the
ring-theoretic version of "local behavior of a function." Given a multiplicatively closed
set ``D ∋ 1`` of a commutative ring ``R``, the ring of fractions ``D⁻¹R`` is the *smallest*
ring in which every ``d ∈ D`` becomes a unit. Two cases carry almost everything:

**Inverting one element ``f``** (``D = {1, f, f², …}``). Then ``R_f ≅ R[t]/(tf − 1)`` — adjoin
a formal inverse ``t = 1/f``. Geometrically ``Spec(R_f) ≅ D(f)``, the distinguished open
where ``f`` does not vanish: localizing at ``f`` is *deleting the zero set of ``f``*. For
``R = ℤ`` this is ``ℤ[1/n]``, whose primes are exactly the ``(p)`` with ``p ∤ n``.

**Localizing at a prime ``P``** (``D = R ∖ P``). Every element *outside* ``P`` becomes a unit,
leaving ``R_P`` a **local ring**: a unique maximal ideal ``P·R_P``. The prototype is
``ℤ_(p) = { a/b ∈ ℚ : p ∤ b }`` — a **discrete valuation ring** whose only primes are ``(0)``
and ``(p)``, graded by the ``p``-adic valuation ``v_p``. Units are ``v_p = 0``; the maximal
ideal is ``v_p ≥ 1``. Localizing keeps only the arithmetic visible "at ``p``", and the global
ring is recovered as the intersection ``ℤ = ⋂_p ℤ_(p)`` — local data reassembling the global
ring, the same slogan that drives the rest of the package.

This module computes the concrete cases (everything reduces to ``p``-adic valuations and the
quotient construction ``R[t]/(tf−1)``):

- :func:`valuation_p`, :func:`in_local_ring`, :func:`is_unit_in_local_ring` — membership and
  units of ``ℤ_(p)`` by valuation; :func:`maximal_ideal_membership`.
- :func:`primes_of_localization_at_prime` — the two primes of the DVR ``ℤ_(p)``.
- :func:`global_ring_is_intersection` — ``x ∈ ℤ`` iff ``x ∈ ℤ_(p)`` for every prime ``p``.
- :func:`invert_element` — build ``R_f = R[t]/(tf−1)`` and witness that ``f`` is now a unit.
- :func:`spec_away_from` — the primes of ``ℤ[1/n] = D((n))`` (cf. :mod:`spectrum`).
"""

from sage.all import ZZ, QQ, PolynomialRing, prime_range


def valuation_p(x, p):
    """The ``p``-adic valuation ``v_p(a/b) = v_p(a) − v_p(b)`` of a rational ``x``."""
    x = QQ(x)
    if x == 0:
        return ZZ(10) ** 9  # +∞ sentinel: 0 lies in every power of the maximal ideal
    return ZZ(x.numerator().valuation(p) - x.denominator().valuation(p))


def in_local_ring(x, p):
    """Is ``x ∈ ℚ`` in ``ℤ_(p)``? Equivalently ``p`` does not divide the denominator (``v_p ≥ 0``)."""
    return valuation_p(x, p) >= 0


def is_unit_in_local_ring(x, p):
    """Is ``x`` a unit of ``ℤ_(p)``? Equivalently ``v_p(x) = 0`` (numerator and denominator coprime to ``p``)."""
    return valuation_p(x, p) == 0


def maximal_ideal_membership(x, p):
    """Is ``x`` in the unique maximal ideal ``p·ℤ_(p)`` of the local ring? Equivalently ``v_p(x) ≥ 1``."""
    return in_local_ring(x, p) and valuation_p(x, p) >= 1


def primes_of_localization_at_prime(p):
    """
    The prime ideals of the DVR ``ℤ_(p)``: exactly ``(0)`` and ``(p)``.

    The primes of ``R_P`` correspond to the primes of ``R`` *contained in* ``P``; for ``R = ℤ``
    and ``P = (p)`` those are ``(0) ⊂ (p)``. This two-element chain (generic point + one closed
    point) is what makes ``ℤ_(p)`` a discrete valuation ring.
    """
    return [ZZ(0), ZZ(p)]


def global_ring_is_intersection(x, prime_bound=50):
    """
    Verify ``x ∈ ℤ ⟺ x ∈ ℤ_(p)`` for every prime ``p`` (checked up to ``prime_bound``).

    This is the statement ``ℤ = ⋂_p ℤ_(p)``: an integer is precisely a rational that is local at
    *every* prime. Returns ``(is_integer, locally_integral_everywhere, agree)``.
    """
    x = QQ(x)
    is_integer = x.denominator() == 1
    everywhere = all(in_local_ring(x, p) for p in prime_range(prime_bound + 1))
    return is_integer, everywhere, is_integer == everywhere


def invert_element(R, f):
    """
    Build the localization ``R_f = R[t]/(t·f − 1)`` that inverts ``f``, and witness ``f`` a unit.

    Returns
    -------
    dict with ``ring`` (the quotient ``R_f``), ``f_is_unit`` (``f̄ · t̄ = 1`` in ``R_f``), and
    ``note``. For ``R = ℤ`` and ``f = n`` this is ``ℤ[1/n]``.
    """
    S = PolynomialRing(R, "t")
    t = S.gen()
    Q = S.quotient(t * S(f) - 1, names="tbar")
    tbar = Q.gen()
    f_is_unit = (Q(f) * tbar == Q(1))
    return {
        "ring": Q,
        "f_is_unit": bool(f_is_unit),
        "note": f"R_f = R[t]/(t·f - 1); the class of t is 1/f, so f is invertible",
    }


def spec_away_from(n, prime_bound=50):
    """
    The primes of ``ℤ[1/n]`` — the closed points of the distinguished open ``D((n)) ⊆ Spec(ℤ)``.

    Inverting ``n`` deletes exactly the points ``(p)`` with ``p ∣ n``, leaving the primes
    ``(p)`` with ``p ∤ n`` (plus the generic ``(0)``). Mirrors
    :func:`maths_coding.geometry.algebraic_geometry.spectrum.distinguished_open`.
    """
    n = ZZ(n)
    return [p for p in prime_range(prime_bound + 1) if n % p != 0]


if __name__ == "__main__":
    print("=" * 70)
    print("Localization — zooming in on a ring by inverting elements")
    print("=" * 70)

    # ---- Z localized at the prime (5): a discrete valuation ring -------------
    p = 5
    print(f"\nZ_(p) for p = {p}  (the DVR of rationals a/b with p not dividing b):")
    samples = [QQ(3) / 7, QQ(10) / 3, QQ(5) / 1, QQ(1) / 5, QQ(25) / 4]
    for x in samples:
        print(f"  {str(x):>6}:  v_{p} = {valuation_p(x, p):>2}   in Z_(p): {in_local_ring(x, p)!s:>5}"
              f"   unit: {is_unit_in_local_ring(x, p)!s:>5}   in max ideal: {maximal_ideal_membership(x, p)}")
    print("  primes of Z_(5):", [f"({q})" for q in primes_of_localization_at_prime(p)],
          "— only (0) ⊂ (5), so Z_(5) is a DVR")

    # ---- Z is the intersection of all its localizations ----------------------
    print("\nZ = ⋂_p Z_(p):  a rational is an integer iff it is local at every prime")
    for x in [QQ(7), QQ(7) / 10, QQ(22) / 7]:
        isint, everywhere, ok = global_ring_is_intersection(x)
        print(f"  x = {str(x):>5}:  integer? {isint!s:>5}   in every Z_(p)? {everywhere!s:>5}   agree: {ok}")

    # ---- Inverting an element: Z[1/n] and the open it carves out --------------
    print("\n" + "-" * 70)
    print("Inverting f: R_f = R[t]/(t·f − 1) makes f a unit (R = Z):")
    for n in [6, 10]:
        res = invert_element(ZZ, n)
        print(f"  Z[1/{n}]:  f = {n} is a unit -> {res['f_is_unit']}"
              f"    Spec(Z[1/{n}]) = D(({n})) = {[f'({q})' for q in spec_away_from(n, 20)]}")
    print("  (inverting n deletes exactly the points (p) with p | n — 'restrict away from V(f)')")
