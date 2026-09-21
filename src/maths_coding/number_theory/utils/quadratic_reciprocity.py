"""Quadratic residues, reciprocity, and irreducibility — the arithmetic of squares.

The central question is local: *is ``a`` a square modulo the prime ``p``?* The Legendre
symbol ``(a/p)`` answers it with a single value in ``{-1, 0, 1}``, and the module computes
it two independent ways so the agreement is itself instructive:

- :func:`legendre_symbol` / :func:`euler_criterion` — **Euler's criterion**,
  ``(a/p) ≡ a^((p-1)/2) (mod p)``.
- :func:`gauss_lemma` — **Gauss's lemma**, counting how many of ``a, 2a, …, ((p-1)/2)a``
  land in the "negative" half ``(p/2, p)`` mod ``p``.

The headline result, :func:`quadratic_reciprocity`, relates ``(p/q)`` and ``(q/p)`` — a
surprising symmetry between two different primes — together with the supplementary laws
:func:`legendre_minus_one` ``(-1/p)`` and :func:`legendre_two` ``(2/p)``.
:func:`jacobi_symbol` extends the symbol to odd composite "denominators" via fast
reciprocity, needing no factorization.

The module then crosses into two adjacent corners of elementary number theory and algebra
that share the same prime-by-prime flavour:

- :func:`legendre_formula` — **Legendre's formula** for ``vₚ(n!)``, the exact power of a
  prime dividing a factorial.
- :func:`polynomial_content`, :func:`is_primitive`, :func:`gauss_lemma_polynomial`,
  :func:`eisenstein_criterion` — **Gauss's lemma for polynomials** and the **Eisenstein
  criterion**, prime-based irreducibility tests over ``Q``.
"""

from sage.all import ZZ, is_prime, gcd, PolynomialRing, QQ, GF
import math

def modexp(a, e, m):
    return pow(int(a) % int(m), int(e), int(m))

def legendre_symbol(a, p):
    """
    Compute (a/p) for odd prime p using Euler's criterion.
    Returns -1, 0, or 1.
    """
    a = ZZ(a) % p
    p = ZZ(p)

    if p < 2 or not is_prime(p) or p % 2 == 0:
        raise ValueError("p must be an odd prime")

    if a == 0:
        return 0

    r = modexp(a, (p - 1) // 2, p)
    if r == 1:
        return 1
    elif r == p - 1:
        return -1
    else:
        # Should not happen for prime p
        return None

def euler_criterion(a, p):
    """
    Check Euler's criterion:
    a^((p-1)/2) ≡ (a/p) mod p
    """
    lhs = modexp(a, (p - 1) // 2, p)
    rhs = legendre_symbol(a, p) % p
    return lhs == rhs

def gauss_lemma(a, p):
    """
    Compute (a/p) using Gauss's lemma.
    Count how many of {a,2a,...,((p-1)/2)a} mod p lie in (p/2, p).
    """
    a = ZZ(a) % p
    p = ZZ(p)

    if p % 2 == 0 or not is_prime(p):
        raise ValueError("p must be an odd prime")

    count = 0
    half = (p - 1) // 2

    for k in range(1, half + 1):
        val = (a * k) % p
        if val > p // 2:
            count += 1

    return (-1) ** count

def quadratic_reciprocity(p, q):
    """
    Verify quadratic reciprocity:
    (p/q)*(q/p) = (-1)^((p-1)(q-1)/4)
    """
    if not (is_prime(p) and is_prime(q) and p % 2 == 1 and q % 2 == 1):
        raise ValueError("p, q must be odd primes")

    lhs = legendre_symbol(p, q) * legendre_symbol(q, p)
    rhs = (-1) ** (((p - 1) * (q - 1)) // 4)

    return lhs, rhs, lhs == rhs

def legendre_minus_one(p):
    """(-1/p)"""
    return (-1) ** ((p - 1) // 2)

def legendre_two(p):
    """(2/p)"""
    return (-1) ** ((p * p - 1) // 8)


def jacobi_symbol(a, n):
    """
    Jacobi symbol (a/n) for an odd positive integer n.

    The Jacobi symbol generalises the Legendre symbol multiplicatively over the
    prime factorisation n = ∏ pᵢ^eᵢ:

        (a/n) = ∏ (a/pᵢ)^eᵢ.

    It is computed here with the classic fast algorithm built from quadratic
    reciprocity and the supplementary law for 2 — so it needs **no factorisation**
    of n, which is exactly why the Jacobi symbol is the practical tool in primality
    testing (e.g. Solovay–Strassen).

    Parameters
    ----------
    a : int
    n : int
        Must be a positive odd integer.

    Returns
    -------
    int
        -1, 0, or 1.

    Notes
    -----
    For prime n the Jacobi symbol coincides with the Legendre symbol. For composite
    n, however, ``(a/n) = 1`` does **not** imply that a is a quadratic residue mod n;
    only ``(a/n) = -1`` is conclusive (it guarantees a is a non-residue).
    """
    n = ZZ(n)
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer")

    a = ZZ(a) % n
    result = 1

    while a != 0:
        # Pull out factors of 2 using (2/n) = (-1)^((n²-1)/8).
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result

        # Reciprocity flip (a/n) ↔ (n/a), with a sign when both ≡ 3 (mod 4).
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result

        a %= n

    return result if n == 1 else 0


# --------------------------------------------------
# Legendre's formula (v_p(n!))
# --------------------------------------------------

def legendre_formula(n, p):
    """
    v_p(n!) = sum_{k>=1} floor(n / p^k)
    """
    n = ZZ(n)
    p = ZZ(p)

    if p <= 1 or not is_prime(p):
        raise ValueError("p must be prime")

    total = 0
    power = p

    while power <= n:
        total += n // power
        power *= p

    return total


# --------------------------------------------------
# Polynomial tools: Gauss's lemma + Eisenstein
# --------------------------------------------------

def polynomial_content(f):
    """
    Content = gcd of coefficients
    """
    coeffs = [ZZ(c) for c in f.coefficients()]
    if not coeffs:
        return ZZ(0)
    return gcd(coeffs)


def is_primitive(f):
    return polynomial_content(f) == 1


def gauss_lemma_polynomial(f):
    """
    Gauss's lemma: primitive polynomial irreducible over Q
    iff irreducible over Z
    (this just checks primitivity)
    """
    return is_primitive(f)


def eisenstein_criterion(f, p):
    """
    Eisenstein irreducibility test at prime p:
    - p divides all coefficients except leading
    - p^2 does NOT divide constant term
    """
    R = f.parent()
    coeffs = f.list()  # low → high degree
    p = ZZ(p)

    if not is_prime(p):
        raise ValueError("p must be prime")

    # leading coefficient
    if coeffs[-1] % p == 0:
        return False

    # all others divisible by p
    for c in coeffs[:-1]:
        if c % p != 0:
            return False

    # constant not divisible by p^2
    if coeffs[0] % (p * p) == 0:
        return False

    return True