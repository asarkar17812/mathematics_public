"""The Lucas–Lehmer test, and why it works — a closed form in Q(√3).

The Lucas–Lehmer test decides primality of a Mersenne number ``M_p = 2ᵖ - 1`` by iterating

    s₀ = 4,   s_{n+1} = s_n² - 2   (mod M_p),

and declaring ``M_p`` prime iff ``s_{p-2} ≡ 0``. That this works is not obvious from the
recurrence alone. The key is that the sequence has a **closed form**,

    s_n = (2 + √3)^{2ⁿ} + (2 - √3)^{2ⁿ},

living in the real quadratic field ``Q(√3)``. (The squaring step is just the identity
``(α^{2ⁿ} + ᾱ^{2ⁿ})² - 2 = α^{2^{n+1}} + ᾱ^{2^{n+1}}`` for ``α = 2 + √3``, ``αᾱ = 1``.) The
primality criterion then comes from analysing this expression in the ring ``Z[√3]/M_p``.

This module deliberately keeps **all three** viewpoints side by side so the connection is
explicit and checkable:

- :func:`lucas_lehmer_value_recursive` — the bare recurrence (exact but exponential cost);
- :func:`lucas_lehmer_value_analytic` — the closed form, evaluated exactly in ``Q(√3)``;
- :func:`lucas_lehmer_test` / :func:`lucas_lehmer_test_with_value` — the practical modular
  test, reducing mod ``M_p`` at every step to keep the numbers small;
- :func:`lucas_lehmer_test_analytic` — the test driven by the closed form, to witness that
  the two agree.

The dataset generator :mod:`maths_coding.number_theory.utils.mersenne_generator` calls
:func:`lucas_lehmer_test_with_value` to sweep across exponents.
"""

from sage.all import *

K = QuadraticField(3, 'a')
a = K.gen()  # generator representing sqrt(3)

def lucas_lehmer_value_recursive(n):
    """
    Compute the n-th term of the Lucas–Lehmer sequence recursively.

    The sequence is defined by:
        s_0 = 4
        s_{n+1} = s_n^2 - 2

    Parameters
    ----------
    n : int
        Index of the sequence term.

    Returns
    -------
    sage.rings.integer.Integer
        The exact value of s_n.

    Notes
    -----
    This implementation is mathematically correct but inefficient
    due to repeated recomputation (exponential recursion).
    """
    if n == 0:
        return ZZ(4)
    else:
        return lucas_lehmer_value_recursive(n - 1)**2 - 2


def lucas_lehmer_value_analytic(n):
    """
    Compute the n-th Lucas–Lehmer value using a closed-form expression.

    The identity used is:
        s_n = (2 + sqrt(3))^(2^n) + (2 - sqrt(3))^(2^n)

    This is evaluated exactly in the number field Q(sqrt(3)),
    avoiding floating-point errors.

    Parameters
    ----------
    n : int
        Index of the sequence term.

    Returns
    -------
    sage.rings.integer.Integer
        Exact integer value of s_n.
    """
    x = 2 + a
    y = 2 - a
    return ZZ(x**(2**n) + y**(2**n))


def mersenne_number_calc(p):
    """
    Compute the Mersenne number M_p = 2^p - 1.

    Parameters
    ----------
    p : int
        Exponent.

    Returns
    -------
    sage.rings.integer.Integer
        The Mersenne number M_p.
    """
    return ZZ(2)**p - 1


def lucas_lehmer_test(p):
    """
    Perform the Lucas–Lehmer primality test for Mersenne numbers.

    The test determines whether M_p = 2^p - 1 is prime
    (for prime p).

    Algorithm:
        s = 4
        repeat p-2 times:
            s = s^2 - 2 mod M_p
        M_p is prime iff s == 0

    Parameters
    ----------
    p : int
        Exponent (should be prime).

    Returns
    -------
    bool
        True if M_p is prime, False otherwise.
    """
    if p == 2:
        return True

    M = mersenne_number_calc(p)
    s = ZZ(4)

    for _ in range(p - 2):
        s = (s*s - 2) % M

    return s == 0


def lucas_lehmer_test_analytic(p):
    """
    Lucas–Lehmer test using the analytic (closed-form) sequence.

    This is not efficient for large p, but demonstrates the
    connection between the analytic formula and the modular test.

    Parameters
    ----------
    p : int

    Returns
    -------
    bool
        True if M_p is prime, False otherwise.
    """
    M = mersenne_number_calc(p)
    s = lucas_lehmer_value_analytic(p - 2)
    return s % M == 0


def lucas_lehmer_sequence_value(p):
    """
    Compute the Lucas–Lehmer sequence value s_{p-2} modulo M_p.

    This is the value actually used in the Lucas–Lehmer test.

    Parameters
    ----------
    p : int

    Returns
    -------
    sage.rings.integer.Integer
        The final sequence value modulo M_p.
        For Mersenne primes, this will be 0.
    """
    M = ZZ(2)**p - 1
    s = ZZ(4)

    for _ in range(p - 2):
        s = (s*s - 2) % M

    return s


def lucas_lehmer_test_with_value(p):
    """
    Perform the Lucas–Lehmer test and return both:
        - primality result
        - final sequence value s_{p-2} mod M_p

    This avoids recomputing the sequence twice.
    """
    if p == 2:
        return True, ZZ(0)

    M = mersenne_number_calc(p)
    s = ZZ(4)

    for _ in range(p - 2):
        s = (s*s - 2) % M

    return s == 0, s

def lucas_lehmer_step(s, M):
    """Apply one Lucas–Lehmer step: return ``(s² - 2) mod M``."""
    return (s*s - 2) % M


if __name__ == "__main__":
    # Demo: the closed-form value s_11 evaluated exactly in Q(√3).
    print(lucas_lehmer_value_analytic(11))