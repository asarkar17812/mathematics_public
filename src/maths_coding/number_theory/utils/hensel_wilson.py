"""Hensel lifting and Wilson's theorems — local roots and the structure of unit groups.

Two classical prime-by-prime results, and their synthesis (the "Henson" combination of the
notes), all instances of the package's recurring slogan that **local data reconstructs
global structure**.

**Hensel's lemma** is Newton's method done ``p``-adically. If ``f ∈ ℤ[x]`` has a *simple*
root mod ``p`` — meaning ``f(a₁) ≡ 0`` but ``f'(a₁) ≢ 0 (mod p)`` — then that root lifts
**uniquely** all the way up the tower of moduli ``p, p², p³, …`` by the Newton iteration

    a_{k+1} = a_k − f(a_k) · f'(a_k)⁻¹   (mod p^{k+1}).

The non-vanishing derivative is exactly the ``p``-adic inverse function theorem: the linear
Taylor term is invertible, so among the ``p`` candidate lifts of ``a_k`` *exactly one* is
again a root. When ``f'(a₁) ≡ 0`` (the *ramified* case, e.g. ``x² − p`` at ``0``) the
branching is not pruned and the lemma can fail. Taking the inverse limit assembles the lifts
into a genuine root in ``ℤ_p`` — this is how ``√2`` and roots of unity come to live there.

**Wilson's theorem** computes the product of all units. Classically ``(p−1)! ≡ −1 (mod p)``,
proved by pairing each element of ``𝔽_p^×`` with its inverse — only the self-inverses
``±1`` survive — and the converse makes it an (inefficient) primality test. The higher-order
form says the product of *all* units of ``(ℤ/p^kℤ)^×`` is ``≡ −1`` for odd primes ``p``.

The **Henson synthesis** chains them with the CRT: use residuosity to find a root mod ``p``,
Hensel-lift it mod ``p^k``, and glue across primes — recovering, e.g., ``√2`` modulo ``7^k``
and (combined) modulo ``N = ∏ pᵢ^{kᵢ}``.

- :func:`is_simple_root_mod_p` — the Hensel hypothesis ``f(a)≡0, f'(a)≢0 (mod p)``.
- :func:`hensel_lift`, :func:`hensel_lift_tower` — Newton lifting, and the whole tower.
- :func:`wilson_theorem`, :func:`wilson_converse`, :func:`wilson_higher` — the unit products.
- :func:`sqrt_padic` — a square root in ``ℤ_p`` by Hensel lifting (the "Henson" route).
"""

from sage.all import (
    ZZ,
    PolynomialRing,
    is_prime,
    factorial,
    gcd,
    inverse_mod,
    prod,
    CRT_list,
)


def _poly(f):
    """Coerce ``f`` to a univariate integer polynomial and return ``(f, f', x)``."""
    R = PolynomialRing(ZZ, "x")
    f = R(f)
    return f, f.derivative(), R.gen()


def is_simple_root_mod_p(f, p, a):
    """
    Test the Hensel hypothesis at ``a``: ``f(a) ≡ 0`` and ``f'(a) ≢ 0 (mod p)``.

    A *simple* root mod ``p`` is precisely one that Hensel lifts uniquely; a root with
    ``f'(a) ≡ 0`` is a *ramified* (multiple) root where lifting may fail.
    """
    f, df, _ = _poly(f)
    return (int(f(a)) % p == 0) and (int(df(a)) % p != 0)


def hensel_lift(f, p, a1, k):
    """
    Lift a simple root ``a₁`` of ``f`` mod ``p`` to the unique root mod ``p^k``.

    Iterates the Newton formula ``a ← a − f(a)·f'(a)⁻¹ (mod p^{level})``. Raises
    ``ValueError`` if ``a₁`` is not a simple root mod ``p`` (so the hypothesis is honest).

    Parameters
    ----------
    f : polynomial / coercible to ``ZZ[x]``
    p : int (prime)
    a1 : int
        A root of ``f`` mod ``p`` with ``f'(a₁) ≢ 0``.
    k : int
        Target exponent; return value lies in ``ℤ/p^kℤ`` represented in ``[0, p^k)``.

    Returns
    -------
    int
        The unique ``a_k`` with ``f(a_k) ≡ 0 (mod p^k)`` and ``a_k ≡ a₁ (mod p)``.
    """
    if not is_prime(p):
        raise ValueError("p must be prime")
    if not is_simple_root_mod_p(f, p, a1):
        raise ValueError(f"{a1} is not a simple root of f mod {p} (Hensel hypothesis fails)")

    f, df, _ = _poly(f)
    a = int(a1) % p
    for level in range(2, k + 1):
        mod = p ** level
        dfa = int(df(a)) % mod
        inv = inverse_mod(dfa, mod)             # f'(a) is a unit since f'(a₁) ≢ 0 mod p
        a = (a - int(f(a)) * inv) % mod
    return a % (p ** k)


def hensel_lift_tower(f, p, a1, levels):
    """
    Return the whole lifting tower ``[(k, a_k, f(a_k) mod p^k)]`` for ``k = 1..levels``.

    Each row exhibits a root mod ``p^k`` and confirms ``f(a_k) ≡ 0`` there — the "solid
    nodes" of the lifting-tower picture, the unique survivors among the ``p`` candidate
    lifts at each level.
    """
    f_poly, _, _ = _poly(f)
    rows = []
    for k in range(1, levels + 1):
        ak = hensel_lift(f, p, a1, k)
        rows.append((k, ak, int(f_poly(ak)) % (p ** k)))
    return rows


def wilson_theorem(p):
    """
    Verify Wilson's theorem ``(p−1)! ≡ −1 (mod p)`` for a prime ``p``.

    Returns
    -------
    (lhs, holds) : tuple
        ``lhs = (p−1)! mod p`` (which is ``p−1 ≡ −1``) and a boolean.
    """
    lhs = int(factorial(p - 1) % p)
    return lhs, lhs == p - 1


def wilson_converse(n):
    """
    Wilson primality test: ``n > 1`` is prime **iff** ``(n−1)! ≡ −1 (mod n)``.

    Correct but exponentially slow — included to make the *characterisation* concrete, not
    as a practical test. Returns ``True``/``False``.
    """
    n = int(n)
    if n < 2:
        return False
    return int(factorial(n - 1) % n) == n - 1


def wilson_higher(p, k):
    """
    Wilson II: the product of the units of ``(ℤ/p^kℤ)^×`` is ``≡ −1 (mod p^k)`` for an odd
    prime ``p`` (the group is cyclic of order ``φ(p^k) = p^{k-1}(p−1)``).

    Returns
    -------
    (product, expected, holds) : tuple
        ``product`` is the unit product mod ``p^k``; ``expected = p^k − 1 ≡ −1``.

    Notes
    -----
    The case ``p = 2`` is genuinely different: ``(ℤ/2^kℤ)^×`` is non-cyclic for ``k ≥ 3``
    and the unit product is ``+1`` there — so this routine restricts to odd ``p``.
    """
    if p == 2:
        raise ValueError("Wilson II as stated is for odd primes; p = 2 needs special care")
    m = p ** k
    product = int(prod(a for a in range(1, m) if gcd(a, p) == 1) % m)
    expected = m - 1
    return product, expected, product == expected


def sqrt_padic(a, p, levels):
    """
    A square root of ``a`` in ``ℤ_p`` by Hensel-lifting a root of ``f(x) = x² − a``.

    This is the "Henson" route of the notes: residuosity (``a`` a quadratic residue mod
    ``p``, ``p`` odd, ``p ∤ a``) supplies the simple root mod ``p``, then Hensel lifts it up
    the tower. Returns the approximation mod ``p^levels`` together with the tower; ``None``
    if ``a`` is a non-residue (no root to lift).

    Returns
    -------
    dict with ``root`` (int mod ``p^levels`` or None), ``tower`` (list), ``is_residue`` (bool).
    """
    R = PolynomialRing(ZZ, "x")
    x = R.gen()
    f = x ** 2 - a

    # Find a simple root mod p (Euler's criterion decides residuosity for odd p ∤ a).
    a1 = None
    for r in range(p):
        if (r * r - a) % p == 0 and (2 * r) % p != 0:
            a1 = r
            break
    if a1 is None:
        return {"root": None, "tower": [], "is_residue": False}

    tower = hensel_lift_tower(f, p, a1, levels)
    return {"root": tower[-1][1], "tower": tower, "is_residue": True}


if __name__ == "__main__":
    R = PolynomialRing(ZZ, "x")
    x = R.gen()

    print("=" * 70)
    print("Wilson's theorem — the product of units collapses to −1")
    print("=" * 70)
    for p in [5, 7, 11, 13]:
        lhs, ok = wilson_theorem(p)
        print(f"  (p-1)! ≡ {lhs} ≡ -1 (mod {p})", "✓" if ok else "✗")
    print("  Wilson converse as a primality test on 2..15:")
    print("   ", [n for n in range(2, 16) if wilson_converse(n)], "(exactly the primes)")
    print("  Wilson II  (odd prime powers, product of units ≡ -1):")
    for (p, k) in [(3, 2), (5, 2), (7, 2), (3, 3)]:
        prod_, exp_, ok = wilson_higher(p, k)
        print(f"    Π units of (Z/{p**k})^× = {prod_} ≡ -1 (mod {p**k})", "✓" if ok else "✗")

    print("\n" + "=" * 70)
    print("Hensel lifting — a simple root climbs the tower p, p², p³, …")
    print("=" * 70)
    # √2 in Z_7, the worked example from the notes: a1 = 3, a2 = 10, ...
    f = x ** 2 - 2
    print("\nf(x) = x² − 2, prime p = 7  (find √2 in Z_7):")
    print("  3² = 9 ≡ 2 (mod 7), and f'(3) = 6 ≢ 0, so a₁ = 3 is a simple root.")
    res = sqrt_padic(2, 7, 5)
    for k, ak, fk in res["tower"]:
        print(f"   mod 7^{k} = {7**k:>6}:  a_{k} = {ak:>6},  a_{k}² − 2 ≡ {fk} (mod 7^{k})")
    print(f"  → √2 ≈ {res['root']} (mod 7^5);  note 10² = 100 ≡ 2 (mod 49). ✓")

    # Ramified failure: x² − 7 mod 7 has the root 0 but f'(0) = 0.
    print("\nRamified case x² − 7 at p = 7:  0 is a root mod 7 but f'(0) = 0 —")
    print("  is_simple_root_mod_p(x²−7, 7, 0) =", is_simple_root_mod_p(x ** 2 - 7, 7, 0),
          " (Hensel does not apply; the lift is not forced).")

    print("\n" + "=" * 70)
    print("Henson synthesis — Hensel per prime, then CRT across primes")
    print("=" * 70)
    # Solve x² ≡ 2 simultaneously mod 7² and mod 17² and glue with the CRT.
    r7 = sqrt_padic(2, 7, 2)["root"]    # √2 mod 49
    r17 = sqrt_padic(2, 17, 2)["root"]  # √2 mod 289 (2 is a QR mod 17: 6²=36≡2)
    if r7 is not None and r17 is not None:
        N = 49 * 289
        glued = int(CRT_list([r7, r17], [49, 289]))
        print(f"  √2 ≡ {r7} (mod 49),  √2 ≡ {r17} (mod 289)")
        print(f"  CRT ⇒ x = {glued} (mod {N}),  x² − 2 ≡ {(glued*glued - 2) % N} (mod {N}) ✓")
