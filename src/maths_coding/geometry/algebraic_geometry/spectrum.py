"""Spec of a ring and the Zariski topology: turning a ring into a space.

The prime spectrum of a commutative ring ``R`` is

    Spec(R) = { P ⊴ R : P prime },

topologized by declaring the **closed sets** to be ``V(I) = { P : I ⊆ P }`` for ideals
``I``, with basic **distinguished opens** ``D(f) = { P : f ∉ P } = Spec(R) \\ V(f)``. The
closed points are the maximal ideals; a non-maximal prime is a **generic point** whose
closure ``\\overline{\\{P\\}} = V(P)`` is everything specializing from it. ``Spec`` is a
*contravariant* functor: a ring map ``φ : R → S`` induces a continuous
``φ* : Spec(S) → Spec(R)``, ``P ↦ φ⁻¹(P)`` — the algebraic backbone of scheme theory.

The clarifying image (from the notes) is that **Spec(ℤ) is a curve — the arithmetic line**:
one closed point ``(p)`` per prime, plus the generic point ``(0)`` whose closure is the whole
space. The integers behave like *functions on a line*, and "reduction mod ``p``" is
"evaluation at the point ``(p)``", since the residue field ``ℤ/(p) = 𝔽_p``. The quotient
map ``ℤ → ℤ/n`` realises ``Spec(ℤ/n) ↪ V((n)) ⊂ Spec(ℤ)`` — exactly the closed points
``(p)`` with ``p ∣ n``. The polynomial ring ``k[x]`` carries the *same* shape (closed points
``(x−a)`` plus a generic point), which is the precise sense in which arithmetic and geometry
are analogues — both are Dedekind, so both spectra are smooth curves.

Everything here is computed concretely from **factorization** (the primes of ``ℤ`` and of
``k[x]`` are governed by prime numbers and irreducible polynomials respectively):

- :func:`spec_ZZ` — the points of ``Spec(ℤ)`` up to a bound: closed points + generic point.
- :func:`closed_set`, :func:`distinguished_open` — ``V((n))`` and ``D(n)`` in ``Spec(ℤ)``.
- :func:`closure` — the closure of a point (specialization).
- :func:`spec_Zmod`, :func:`reduction_pullback` — ``Spec(ℤ/n)`` and its image ``V((n))``.
- :func:`spec_polynomial_quotient` — ``Spec(k[x]/(f))`` as the irreducible factors of ``f``,
  the polynomial mirror of ``Spec(ℤ/n)``.
"""

from sage.all import ZZ, prime_divisors, prime_range


# A point is recorded as a tuple: ("closed", p) for a maximal ideal (p), or
# ("generic", 0) for the zero ideal of a domain.


def spec_ZZ(prime_bound=30):
    """
    The points of ``Spec(ℤ)`` with ``(p)`` for ``p ≤ prime_bound``, plus the generic ``(0)``.

    Returns
    -------
    dict with ``closed_points`` (list of primes), ``generic`` (the ideal ``(0)``),
    and ``note`` describing the geometry.
    """
    return {
        "closed_points": list(prime_range(prime_bound + 1)),
        "generic": 0,
        "note": "closed points (p) are maximal (residue field F_p); (0) is the generic point",
    }


def closed_set(n, prime_bound=30):
    """
    The closed set ``V((n)) = { P ∈ Spec(ℤ) : (n) ⊆ P }`` — the primes containing ``n``.

    For ``n ≠ 0`` these are exactly the closed points ``(p)`` with ``p ∣ n`` (the generic
    point ``(0)`` is *not* in ``V((n))``). For ``n = 0``, ``V((0)) = Spec(ℤ)`` (every prime
    contains 0). Returns the list of prime divisors (or the string ``"all of Spec(Z)"``).
    """
    n = ZZ(n)
    if n == 0:
        return "all of Spec(Z) (the generic point (0) and every (p))"
    if n.abs() == 1:
        return []  # (1) is contained in no prime; V((1)) = ∅
    return prime_divisors(n)


def distinguished_open(n, prime_bound=30):
    """
    The distinguished open ``D(n) = Spec(ℤ) \\ V((n))`` up to ``prime_bound``.

    These are the primes ``(p)`` with ``p ∤ n`` (where ``n`` is invertible — "where the
    function ``n`` does not vanish"), together with the generic point. Returns the closed
    points in ``D(n)``; the generic point ``(0)`` lies in every nonempty ``D(n)``.
    """
    divisors_of_n = set(closed_set(n, prime_bound)) if ZZ(n) != 0 else set()
    return [p for p in prime_range(prime_bound + 1) if p not in divisors_of_n]


def closure(point, prime_bound=30):
    """
    The Zariski closure of a single point of ``Spec(ℤ)`` (the specialization order).

    A closed point ``(p)`` is its own closure; the generic point ``(0)`` has closure equal to
    the *entire* space — it is the "point everywhere at once". Pass ``("closed", p)`` or
    ``("generic", 0)``.
    """
    kind, val = point
    if kind == "closed":
        return [val]
    return spec_ZZ(prime_bound)  # closure of the generic point is all of Spec(Z)


def spec_Zmod(n):
    """
    ``Spec(ℤ/n)``: the prime ideals of ``ℤ/n``, which correspond to the prime divisors of ``n``.

    (``ℤ/n`` is a domain — a single generic point — only when ``n`` is prime; in general its
    primes are the ``(p)/(n)`` for ``p ∣ n``, all maximal since ``ℤ/n`` is Artinian.)
    """
    return prime_divisors(ZZ(n))


def reduction_pullback(n):
    """
    The map ``Spec(ℤ/n) → Spec(ℤ)`` induced by the quotient ``ℤ → ℤ/n`` (functoriality of Spec).

    Each prime of ``ℤ/n`` pulls back to the prime ``(p) ⊂ ℤ`` with ``p ∣ n``, so the image is
    exactly ``V((n))`` — the closed subscheme cut out by ``n``. This is "reduction mod n picks
    out the points ``(p)`` dividing ``n``." Returns a dict exhibiting image = ``V((n))``.
    """
    image = spec_Zmod(n)
    return {
        "image": image,
        "equals_V_of_n": image == closed_set(n),
        "note": "Spec(Z/n) ↪ Spec(Z) has image V((n)) = {(p) : p | n}",
    }


def spec_polynomial_quotient(f):
    """
    ``Spec(k[x]/(f))`` as the distinct **irreducible factors** of ``f`` — the polynomial
    mirror of :func:`spec_Zmod` (irreducible polynomials play the role of prime numbers).

    Returns a list of ``(irreducible_factor, multiplicity)`` pairs; a factor with multiplicity
    ``> 1`` is a non-reduced (fat) point, exactly as a repeated prime makes ``ℤ/n`` non-reduced.
    """
    return [(g, ZZ(e)) for g, e in f.factor()]


if __name__ == "__main__":
    print("=" * 70)
    print("Spec(Z) — the arithmetic line, and the Zariski topology")
    print("=" * 70)

    S = spec_ZZ(20)
    print("\nSpec(Z) up to 20:")
    print("  closed points :", [f"({p})" for p in S["closed_points"]])
    print("  generic point : (0)   — its closure is ALL of Spec(Z)")
    print("  (", S["note"], ")")

    print("\nClosed sets V((n)) = {(p) : p | n}  (and distinguished opens D(n)):")
    for n in [12, 30, 7]:
        print(f"  V(({n})) = {[f'({p})' for p in closed_set(n)]:<28}"
              f"   D({n}) up to 20 = {[f'({p})' for p in distinguished_open(n, 20)]}")

    print("\nSpecialization (closure of a point):")
    print("  closure of closed point (5) :", [f"({p})" for p in closure(('closed', 5))])
    print("  closure of generic point (0): all of Spec(Z)  (0) ⤳ (p) for every p")

    print("\n" + "-" * 70)
    print("Reduction mod n is evaluation: Spec(Z/n) ↪ V((n)) ⊂ Spec(Z)")
    for n in [12, 45]:
        pb = reduction_pullback(n)
        print(f"  Spec(Z/{n}) = {[f'({p})' for p in spec_Zmod(n)]}"
              f"  →  image V(({n})) match: {pb['equals_V_of_n']}")

    print("\n" + "-" * 70)
    print("Same shape over k[x]: Spec(F_q[x]/(f)) ↔ irreducible factors of f")
    from sage.all import GF, PolynomialRing
    R = PolynomialRing(GF(2), "x")
    x = R.gen()
    f = (x ** 3 + x + 1) * (x + 1) ** 2          # one cubic prime + a fat point (x+1)^2
    print(f"  over F_2, f = (x^3+x+1)(x+1)^2 factors into 'points':")
    for g, e in spec_polynomial_quotient(f):
        tag = "reduced point" if e == 1 else f"fat point (mult {e})"
        print(f"    ({g})   deg {g.degree()}   — {tag}")
    print("  (irreducible polynomials are the 'primes' of k[x], just as p are the primes of Z)")
