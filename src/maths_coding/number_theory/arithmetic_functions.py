"""Arithmetic functions and Dirichlet convolution — the ring where number theory adds up.

An **arithmetic function** is any ``f : ℤ_{>0} → ℂ``. The cast of characters:

    μ(n)   Möbius:        0 if n is divisible by a square, else (−1)^{#prime factors}
    φ(n)   Euler totient:  the count of units mod n
    τ(n)=σ₀, σ(n)=σ₁:      the number of divisors, and their sum;  σ_k(n) = Σ_{d∣n} d^k
    1(n)=1,  id(n)=n,  ε(n)=[n=1]   (the constant, identity, and unit functions)

What organizes them is the **Dirichlet convolution** ``(f * g)(n) = Σ_{d∣n} f(d) g(n/d)``, which
makes the arithmetic functions a commutative ring with unit ``ε`` — and the **multiplicative**
functions a subgroup under it. In this ring the classical identities are one-line statements:

    1 * id = σ,     1 * 1 = τ,     μ * 1 = ε,     φ * 1 = id,     φ = μ * id.

The inverse of ``1`` is ``μ``, and that single fact is **Möbius inversion**: if
``g = f * 1`` (i.e. ``g(n) = Σ_{d∣n} f(d)``) then ``f = g * μ`` (i.e.
``f(n) = Σ_{d∣n} μ(d) g(n/d)``) — the number-theoretic inclusion–exclusion that recovers a
function from its divisor-sums.

Pure-Python and exact, so it runs without SageMath:

- :func:`mobius`, :func:`euler_phi`, :func:`divisor_sigma`, :func:`num_divisors` — the functions.
- :func:`dirichlet_convolution` — ``(f * g)(n)``.
- :func:`mobius_inversion` — recover ``f`` from its divisor-sum ``g``.
- :func:`is_multiplicative` — the numerical test ``f(ab) = f(a)f(b)`` on coprime ``a, b``.
"""

from math import gcd


def _factorize(n):
    """Prime factorization of ``n`` as a dict ``{p: exponent}`` (trial division)."""
    n = int(n)
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def divisors(n):
    """The sorted list of positive divisors of ``n``."""
    n = int(n)
    small = [d for d in range(1, int(n ** 0.5) + 1) if n % d == 0]
    return sorted(set(small) | {n // d for d in small})


def mobius(n):
    """The Möbius function ``μ(n)``: ``0`` if ``n`` is squareful, else ``(−1)^{#distinct primes}``."""
    if n == 1:
        return 1
    f = _factorize(n)
    if any(e > 1 for e in f.values()):
        return 0
    return (-1) ** len(f)


def euler_phi(n):
    """Euler's totient ``φ(n) = n·∏_{p∣n}(1 − 1/p)`` — the count of units modulo ``n``."""
    if n == 1:
        return 1
    result = n
    for p in _factorize(n):
        result -= result // p
    return result


def divisor_sigma(n, k=1):
    """The divisor power sum ``σ_k(n) = Σ_{d∣n} d^k`` (``k=0`` is the divisor count ``τ``)."""
    return sum(d ** k for d in divisors(n))


def num_divisors(n):
    """The number of divisors ``τ(n) = σ₀(n)``."""
    return divisor_sigma(n, 0)


def dirichlet_convolution(f, g, n):
    """The Dirichlet convolution ``(f * g)(n) = Σ_{d∣n} f(d) g(n/d)``."""
    return sum(f(d) * g(n // d) for d in divisors(n))


def mobius_inversion(g, n):
    """
    Recover ``f`` from ``g = f * 1`` by Möbius inversion: ``f(n) = Σ_{d∣n} μ(d) g(n/d)``.

    ``g`` is the divisor-sum function ``g(n) = Σ_{d∣n} f(d)``; this returns ``f(n)``.
    """
    return sum(mobius(d) * g(n // d) for d in divisors(n))


def is_multiplicative(f, bound=40):
    """
    Numerically test whether ``f`` is multiplicative: ``f(ab) = f(a)f(b)`` for all coprime
    ``a, b`` up to ``bound``. (A true proof needs all coprime pairs; this is a strong check.)
    """
    for a in range(1, bound):
        for b in range(1, bound):
            if gcd(a, b) == 1 and f(a * b) != f(a) * f(b):
                return False
    return True


if __name__ == "__main__":
    print("=" * 70)
    print("Arithmetic functions and Dirichlet convolution")
    print("=" * 70)

    print("\n n :  μ(n)  φ(n)  τ(n)  σ(n)")
    for n in range(1, 13):
        print(f"  {n:>2}: {mobius(n):>4} {euler_phi(n):>5} {num_divisors(n):>5} {divisor_sigma(n):>5}")

    one = lambda n: 1
    idf = lambda n: n
    eps = lambda n: 1 if n == 1 else 0

    print("\nClassical identities as Dirichlet convolutions (checked for n = 1..30):")
    checks = {
        "1 * id = σ": all(dirichlet_convolution(one, idf, n) == divisor_sigma(n) for n in range(1, 31)),
        "1 * 1  = τ": all(dirichlet_convolution(one, one, n) == num_divisors(n) for n in range(1, 31)),
        "μ * 1  = ε": all(dirichlet_convolution(mobius, one, n) == eps(n) for n in range(1, 31)),
        "φ * 1  = id": all(dirichlet_convolution(euler_phi, one, n) == n for n in range(1, 31)),
        "φ = μ * id": all(euler_phi(n) == dirichlet_convolution(mobius, idf, n) for n in range(1, 31)),
    }
    for name, ok in checks.items():
        print(f"  {name:<14}: {ok}")

    print("\nMöbius inversion recovers f from its divisor-sum g = f * 1:")
    f = euler_phi                      # take f = φ
    g = lambda n: sum(f(d) for d in divisors(n))    # g = Σ_{d|n} φ(d)  (which equals n)
    ok = all(mobius_inversion(g, n) == f(n) for n in range(1, 31))
    print(f"  g(n) = Σ_(d|n) φ(d) = n;  inversion gives back φ:  {ok}")

    print("\nMultiplicativity (numerical check on coprime pairs):")
    for name, fn in [("μ", mobius), ("φ", euler_phi), ("τ", num_divisors), ("σ", divisor_sigma)]:
        print(f"  {name} multiplicative: {is_multiplicative(fn)}")
