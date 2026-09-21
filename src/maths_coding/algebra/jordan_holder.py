"""The Jordan–Hölder theorem: unique "prime factorization" of a finite group.

A **composition series** for a group ``G`` is a chain of subgroups

    {e} = G₀ ◁ G₁ ◁ ⋯ ◁ Gₙ = G

in which every successive quotient ``Gᵢ/Gᵢ₋₁`` is **simple** (no nontrivial proper normal
subgroup). The quotients are the **composition factors**. The **Jordan–Hölder theorem** says
that although the *chain* is far from unique, the **multiset of composition factors is
uniquely determined** by ``G`` — independent of the series chosen. This is the exact group-
theoretic analogue of the Fundamental Theorem of Arithmetic: every finite group "factors"
uniquely into simple groups, just as every integer factors uniquely into primes. (It is what
makes the Classification of Finite Simple Groups a meaningful endpoint.)

For a *solvable* group every composition factor is cyclic of prime order, so the factor data
reduces to a multiset of primes — and for a cyclic group ``ℤ/n`` that multiset is literally
the prime factorization of ``n``. Two genuinely different maximal chains land on the same
multiset; the demo exhibits this for ``ℤ/12`` (both routes give ``{2, 2, 3} = 12``), mirroring
the Z/12 Hasse diagram in the notes.

Where the analogy stops: the factors do **not** determine the group. ``C₆`` and ``S₃`` share
the composition factors ``{C₂, C₃}`` yet are non-isomorphic — the missing information is the
"glue" between factors, the content of the **extension problem**.

- :func:`chain_orders`, :func:`chain_factor_orders` — the orders and factor orders of a chain.
- :func:`is_composition_series` — every step normal with simple (here: prime-order) quotient.
- :func:`jordan_holder_multiset` — the basis-independent multiset of factor orders.
- :func:`verify_jordan_holder` — two chains, same multiset (Jordan–Hölder in action).
"""

from sage.all import (
    ZZ,
    CyclicPermutationGroup,
    SymmetricGroup,
    is_prime,
)


def chain_orders(chain):
    """The orders ``|Gᵢ|`` along a chain ``[G₀, …, Gₙ]`` (ascending, trivial → G)."""
    return [int(H.order()) for H in chain]


def chain_factor_orders(chain):
    """
    The factor orders ``|Gᵢ/Gᵢ₋₁| = |Gᵢ| / |Gᵢ₋₁|`` along the chain.

    For a composition series these are exactly the orders of the composition factors.
    """
    o = chain_orders(chain)
    return [ZZ(o[i + 1] // o[i]) for i in range(len(o) - 1)]


def is_composition_series(chain):
    """
    Test whether ``chain`` is a composition series: each ``Gᵢ₋₁`` normal in ``Gᵢ`` with a
    **simple** quotient.

    Simplicity of ``Gᵢ/Gᵢ₋₁`` is checked here via its order being prime — valid for the
    solvable groups in this module (a prime-order quotient is cyclic of prime order, hence
    simple). General simple non-abelian factors would need a direct simplicity test.
    """
    for i in range(len(chain) - 1):
        lo, hi = chain[i], chain[i + 1]
        if not lo.is_normal(hi):
            return False
        if not is_prime(ZZ(hi.order() // lo.order())):
            return False
    return True


def jordan_holder_multiset(chain):
    """The Jordan–Hölder invariant: the **sorted multiset** of composition-factor orders."""
    return sorted(chain_factor_orders(chain))


def verify_jordan_holder(chain_a, chain_b):
    """
    Check the Jordan–Hölder conclusion for two composition series of the same group:
    equal length and an equal multiset of factor orders.

    Returns
    -------
    dict with ``factors_a``, ``factors_b`` (the ordered factor lists), ``multiset``,
    ``both_valid`` (each is a genuine composition series), and ``agree`` (the theorem).
    """
    fa = chain_factor_orders(chain_a)
    fb = chain_factor_orders(chain_b)
    return {
        "factors_a": fa,
        "factors_b": fb,
        "multiset": jordan_holder_multiset(chain_a),
        "both_valid": is_composition_series(chain_a) and is_composition_series(chain_b),
        "agree": sorted(fa) == sorted(fb),
    }


if __name__ == "__main__":
    print("=" * 70)
    print("Jordan–Hölder — a finite group's unique 'prime factorization'")
    print("=" * 70)

    # ---- Two composition series of Z/12, both giving {2, 2, 3} ---------------
    G = CyclicPermutationGroup(12)
    g = G.gen()                       # the 12-cycle, a generator of Z/12
    triv = G.subgroup([])

    # Chain A:  1 ◁ <g^6>(order 2) ◁ <g^3>(order 4) ◁ Z/12      factors 2, 2, 3
    chain_A = [triv, G.subgroup([g ** 6]), G.subgroup([g ** 3]), G]
    # Chain B:  1 ◁ <g^4>(order 3) ◁ <g^2>(order 6) ◁ Z/12      factors 3, 2, 2
    chain_B = [triv, G.subgroup([g ** 4]), G.subgroup([g ** 2]), G]

    print("\nZ/12, two genuinely different maximal chains:")
    print("  chain A orders :", chain_orders(chain_A), "→ factors", chain_factor_orders(chain_A))
    print("  chain B orders :", chain_orders(chain_B), "→ factors", chain_factor_orders(chain_B))
    res = verify_jordan_holder(chain_A, chain_B)
    print("  both are composition series :", res["both_valid"])
    print("  same factor multiset        :", res["multiset"],
          "✓" if res["agree"] else "✗")
    primes_of_12 = [p for p, e in ZZ(12).factor() for _ in range(e)]
    print("  Fundamental Theorem of Arithmetic agrees: 12 =",
          " · ".join(str(p) for p in primes_of_12))

    # ---- Where the analogy breaks: C6 vs S3, same factors, different group ----
    C6 = CyclicPermutationGroup(6)
    S3 = SymmetricGroup(3)
    a3 = S3.subgroup([S3("(1,2,3)")])           # A3 ◁ S3, index 2
    chain_C6 = [C6.subgroup([]), C6.subgroup([C6.gen() ** 3]), C6]   # 1 ◁ <g^3>(2) ◁ C6
    chain_S3 = [S3.subgroup([]), a3, S3]                              # 1 ◁ A3(3) ◁ S3

    print("\n" + "-" * 70)
    print("Same composition factors need NOT mean isomorphic groups:")
    print("  C6 factor orders :", sorted(chain_factor_orders(chain_C6)))
    print("  S3 factor orders :", sorted(chain_factor_orders(chain_S3)))
    print("  multisets equal  :", jordan_holder_multiset(chain_C6) == jordan_holder_multiset(chain_S3))
    print("  but C6 ≅ S3 ?    :", C6.is_isomorphic(S3),
          "  (the 'glue' — the extension problem — differs)")
