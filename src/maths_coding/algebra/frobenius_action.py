"""Frobenius cycle types: reading a Galois element off the factorization mod p.

There is a beautiful dictionary between number theory and group theory. For an irreducible
``f ∈ Z[x]`` with Galois group ``G`` and an unramified prime ``p`` (one not dividing the
discriminant), the **Frobenius element at ``p``** is a conjugacy class in ``G``, and:

    Dedekind's theorem
    ──────────────────
    the way ``f`` factors mod ``p`` into irreducibles of degrees ``d₁, d₂, …``
    is exactly the *cycle type* of Frobenius acting as a permutation of the roots.

So you can read a group-theoretic invariant (a conjugacy class of ``G ⊂ Sₙ``) straight off
an arithmetic computation (factoring a polynomial mod ``p``) — no Galois theory required at
the keyboard. Chebotarev's density theorem then says every cycle type occurs, with
frequency proportional to its class size.

The example is ``f = x⁵ - 2``, whose Galois group is the order-20 Frobenius group
``F₂₀ ⊂ S₅``; sweeping a few primes exhibits the different cycle types predicted above.

- :func:`frobenius_cycle_type` — factor ``f`` mod ``p`` and return the multiset of factor
  degrees (the Frobenius cycle type).
- :func:`frobenius_action` — illustrate the ``x ↦ xᵖ`` map by raising each root to the
  ``p``-th power in the splitting field.
- :func:`analyze` — print the group data alongside the cycle type at each prime.
"""

from sage.all import *


def frobenius_cycle_type(poly, p):
    """
    Frobenius cycle type at ``p``: the degrees of the irreducible factors of ``f`` mod ``p``.

    By Dedekind's theorem, for an unramified prime ``p`` this multiset of degrees is the
    cycle type of the Frobenius conjugacy class acting on the roots. The reduction is
    done over ``GF(p)`` — factoring over ``Q`` instead would ignore ``p`` entirely.

    Parameters
    ----------
    poly : sage polynomial over Q (or Z)
    p : int
        A prime not dividing the discriminant of ``poly``.

    Returns
    -------
    list[int]
        Sorted degrees of the irreducible factors of ``poly`` mod ``p``.
    """
    R = PolynomialRing(GF(p), 'x')

    f_mod = R(poly)  # reduce the polynomial modulo p
    fac = f_mod.factor()

    cycle = []
    for g, e in fac:
        cycle.extend([g.degree()] * e)

    return sorted(cycle)


def frobenius_action(poly, p):
    """
    Print each root of ``poly`` alongside its ``p``-th power in the splitting field.

    A concrete illustration of the ``x ↦ xᵖ`` map that gives Frobenius its name. Note
    this raises algebraic numbers in a characteristic-0 splitting field to the ``p``-th
    power for display; the genuine Frobenius automorphism lives in characteristic ``p``
    (see :func:`frobenius_cycle_type`).

    Parameters
    ----------
    poly : sage polynomial over Q
    p : int
    """
    R = PolynomialRing(QQ, 'x')
    f_mod = R(poly)
    K = f_mod.splitting_field('a')

    roots = f_mod.roots(ring=K, multiplicities=False)

    print(f"\nFrobenius action mod {p}:")
    for r in roots:
        print(r, "->", r^p)


def analyze(poly, primes):
    """
    Report the Galois group of ``poly`` and its Frobenius cycle type at each prime.

    Prints the group order, the conjugacy classes as cycle types (the *possible*
    Frobenius shapes), then for each prime in ``primes`` either marks it ramified
    (``disc ≡ 0 mod p``) or prints the observed cycle type and ``x ↦ xᵖ`` action.

    Parameters
    ----------
    poly : sage polynomial over Q
    primes : iterable[int]

    Returns
    -------
    None
        Demonstration routine; output is printed.
    """
    G = poly.galois_group()

    print("Polynomial:", poly)
    print("Group order:", G.order())

    print("\nConjugacy classes (cycle types):")
    for cls in G.conjugacy_classes():
        rep = cls.representative()
        print(rep.cycle_type())

    print("\nFrobenius data:")
    for p in primes:
        if poly.discriminant() % p == 0:
            print(f"p={p}: ramified")
            continue

        ct = frobenius_cycle_type(poly, p)
        print(f"p={p}: cycle type {ct}")
        frobenius_action(poly, p)


if __name__ == "__main__":
    R = PolynomialRing(QQ, 'x')
    x = R.gen()
    f = x^5 - 2

    primes = [3,5,7,11,13]

    analyze(f, primes)