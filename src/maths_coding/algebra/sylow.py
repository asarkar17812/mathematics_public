"""Sylow theory — prime-power subgroups and the counting that constrains group structure.

Lagrange says ``|H|`` divides ``|G|``, but not every divisor is realized by a subgroup. The
**Sylow theorems** repair this for prime powers: writing ``|G| = pᵅ·m`` with ``p ∤ m``, a
**Sylow p-subgroup** is one of order ``pᵅ`` (the largest possible p-power), and

    I   (existence)  a Sylow p-subgroup exists;
    II  (conjugacy)  all Sylow p-subgroups are conjugate, and every p-subgroup sits in one;
    III (counting)   nₚ = |Sylₚ(G)| satisfies  nₚ ≡ 1 (mod p)  and  nₚ ∣ m.

The counting constraint is the workhorse: it often pins ``nₚ`` to a single value, and ``nₚ = 1``
means the Sylow p-subgroup is **normal** — the first foothold for taking a group apart (every
Sylow normal ⟺ ``G`` is nilpotent ⟺ ``G`` is the direct product of its Sylow subgroups). The
engine behind all of it is the **class equation** from the conjugation action of ``G`` on
itself (cf. :mod:`maths_coding.algebra.group_actions`):

    |G| = |Z(G)| + Σ [G : C_G(xᵢ)]   over non-central class representatives,

and **Cauchy's theorem** (``p ∣ |G| ⇒ G`` has an element of order ``p``) is the rank-1 case.

This module verifies the theorems computationally — ``nₚ = [G : N_G(P)]`` by Sylow II — and
reads off normality and nilpotency:

- :func:`class_equation` — the conjugation class equation, split off the center.
- :func:`sylow_data` — ``pᵅ``, ``m``, ``nₚ``, and the two Sylow III congruences, verified.
- :func:`verify_sylow_theorems` — the full ledger over every prime dividing ``|G|``.
- :func:`is_normal_sylow`, :func:`nilpotent_via_sylow` — normality and the nilpotency criterion.
- :func:`cauchy_element` — an explicit element of order ``p``.
"""

from sage.all import ZZ, is_prime, prime_divisors


def class_equation(G):
    """
    The class equation ``|G| = |Z(G)| + Σ [G : C_G(xᵢ)]`` from conjugation.

    Returns
    -------
    dict with ``order``, ``center_order``, ``noncentral_class_sizes`` (each a ``[G:C(x)]``),
    and ``holds`` (the identity checks out).
    """
    order = int(G.order())
    z = int(G.center().order())
    sizes = [int(c.cardinality()) for c in G.conjugacy_classes()]
    noncentral = [s for s in sizes if s > 1]      # central elements form size-1 classes
    return {
        "order": order,
        "center_order": z,
        "noncentral_class_sizes": noncentral,
        "holds": z + sum(noncentral) == order,
    }


def sylow_data(G, p):
    """
    Sylow data at the prime ``p``: the Sylow order ``pᵅ``, the cofactor ``m``, and the count
    ``nₚ = [G : N_G(P)]`` (valid by Sylow II), with the two Sylow III congruences checked.

    Returns
    -------
    dict with ``p``, ``sylow_order`` (``pᵅ``), ``m``, ``n_p``, ``congruence_mod_p`` (``nₚ ≡ 1``),
    ``divides_m`` (``nₚ ∣ m``), and ``normal`` (``nₚ = 1``).
    """
    p = ZZ(p)
    P = G.sylow_subgroup(p)
    sylow_order = int(P.order())
    m = int(G.order()) // sylow_order
    n_p = int(G.order() // G.normalizer(P).order())
    return {
        "p": int(p),
        "sylow_order": sylow_order,
        "m": m,
        "n_p": n_p,
        "congruence_mod_p": n_p % p == 1,
        "divides_m": m % n_p == 0,
        "normal": n_p == 1,
    }


def verify_sylow_theorems(G):
    """Sylow data for every prime dividing ``|G|``; returns a list of :func:`sylow_data` dicts."""
    return [sylow_data(G, p) for p in prime_divisors(ZZ(G.order()))]


def is_normal_sylow(G, p):
    """Whether the Sylow ``p``-subgroup is normal — equivalently ``nₚ = 1``."""
    return sylow_data(G, p)["normal"]


def nilpotent_via_sylow(G):
    """
    Test nilpotency through Sylow theory: ``G`` is nilpotent iff **every** Sylow subgroup is
    normal (equivalently ``G`` is the direct product of its Sylow subgroups).

    Returns
    -------
    dict with ``all_sylow_normal``, ``is_nilpotent`` (Sage's verdict), and ``agree``.
    """
    all_normal = all(is_normal_sylow(G, p) for p in prime_divisors(ZZ(G.order())))
    sage_nilpotent = bool(G.is_nilpotent())
    return {
        "all_sylow_normal": all_normal,
        "is_nilpotent": sage_nilpotent,
        "agree": all_normal == sage_nilpotent,
    }


def cauchy_element(G, p):
    """
    An explicit element of order ``p`` (Cauchy's theorem), found inside a Sylow ``p``-subgroup.

    Returns the element, or ``None`` if ``p ∤ |G|``.
    """
    p = ZZ(p)
    if not is_prime(p) or G.order() % p != 0:
        return None
    P = G.sylow_subgroup(p)
    for x in P:
        if x.order() == p:
            return x
    return None


if __name__ == "__main__":
    from sage.all import (
        SymmetricGroup, AlternatingGroup, CyclicPermutationGroup,
        DihedralGroup, QuaternionGroup,
    )

    print("=" * 70)
    print("Sylow theory — prime-power subgroups and the counting nₚ ≡ 1 (mod p)")
    print("=" * 70)

    groups = [
        ("A4  (order 12 = 2²·3)", AlternatingGroup(4)),
        ("S4  (order 24 = 2³·3)", SymmetricGroup(4)),
        ("Q8  (quaternions, 8)", QuaternionGroup()),
        ("C15 (order 15 = 3·5)", CyclicPermutationGroup(15)),
        ("D6  (order 12)", DihedralGroup(6)),
    ]

    for name, G in groups:
        print(f"\n{name}:")
        ce = class_equation(G)
        print(f"  class equation: |G| = |Z| + Σ[G:C(xᵢ)] = {ce['center_order']} + "
              f"{' + '.join(map(str, ce['noncentral_class_sizes']))} = {ce['order']}  ({ce['holds']})")
        for s in verify_sylow_theorems(G):
            tag = "normal" if s["normal"] else "not normal"
            print(f"  p={s['p']}: Sylow order {s['sylow_order']:>2}, m={s['m']:>2}, "
                  f"n_{s['p']}={s['n_p']:>2}   nₚ≡1 mod p: {s['congruence_mod_p']!s:>5}, "
                  f"nₚ|m: {s['divides_m']!s:>5}  → {tag}")
        nil = nilpotent_via_sylow(G)
        print(f"  all Sylow normal: {nil['all_sylow_normal']!s:>5}  ⇒ nilpotent: "
              f"{nil['is_nilpotent']}  (criterion agrees: {nil['agree']})")

    # Cauchy: an element of each prime order in A4.
    A4 = AlternatingGroup(4)
    print("\nCauchy's theorem in A4 — an explicit element of each prime order:")
    for p in (2, 3):
        g = cauchy_element(A4, p)
        print(f"  order {p}: {g}  (order {g.order()})")
