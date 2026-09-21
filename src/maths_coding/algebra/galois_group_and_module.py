"""The Galois group of a polynomial, made concrete on its roots.

The Galois group ``G = Gal(f)`` of a polynomial ``f`` is abstract until you watch it *act*.
This module takes a polynomial over ``Q``, builds its splitting field, and exhibits the
group in three increasingly concrete ways — the same group seen from three sides:

1. **As an abstract group** — its order, generators, and conjugacy-class sizes.
2. **As permutations of the roots** — each generator sends the list of roots to a
   permutation of itself; this is the embedding ``G ↪ Sₙ`` that *defines* a Galois group as
   a permutation group on the ``n`` roots.
3. **As a permutation representation** — turning each generator into an ``n × n`` 0/1
   permutation matrix, i.e. the linear action of ``G`` on the ``Q``-vector space spanned by
   the roots. (Hence the *module* in the file name: the roots span a ``Q[G]``-module.)

The worked example is ``f = x⁴ - 2``, whose Galois group is the dihedral group ``D₄`` of
order 8 — large enough to be interesting, small enough to print in full.
"""

from sage.all import *

def analyze_galois(poly):
    """
    Print the Galois group of ``poly`` and its action on the roots.

    Computes the splitting field and Galois group of ``poly`` (a polynomial over Q),
    then reports the group order, generators, and conjugacy-class sizes, followed by
    how each generator permutes the roots and the corresponding permutation matrices.

    Parameters
    ----------
    poly : sage polynomial over Q
        The polynomial to analyse (e.g. ``x**4 - 2``).

    Returns
    -------
    None
        Everything is printed; this is a demonstration routine.
    """
    K = poly.splitting_field('a')
    G = poly.galois_group()

    print("Polynomial:", poly)
    print("Group order:", G.order())
    print("Generators:", G.gens())

    print("\nConjugacy classes:")
    for cls in G.conjugacy_classes():
        print("size", len(cls))

    roots = poly.roots(ring=K, multiplicities=False)

    print("\nAction on roots:")
    for g in G.gens():
        print("\nGenerator:", g)
        for r in roots:
            print(r, "->", g(r))

    n = len(roots)
    V = VectorSpace(QQ, n)

    print("\nPermutation representation matrices:")
    for g in G.gens():
        perm = g
        M = matrix(QQ, n, n)
        for i in range(n):
            j = perm(i+1) - 1
            M[j,i] = 1
        print(M)


if __name__ == "__main__":
    R = PolynomialRing(QQ, 'x')
    x = R.gen()
    f = x^4 - 2
    analyze_galois(f)