"""A uniform divisibility test for any prime, from the inverse of 10.

Everyone learns ad-hoc divisibility rules ("a number is divisible by 3 iff its digit sum
is", "...by 7 iff you double the last digit, subtract, and recurse"). These look like a
bag of unrelated tricks, but they are all the *same* trick in disguise.

Write ``n = 10·q + r`` with ``q = ⌊n/10⌋`` and ``r`` the last digit. For a prime ``p`` other
than 2 or 5, ``10`` is invertible mod ``p``; let ``k ≡ 10⁻¹ (mod p)``. Multiplying the
congruence ``n ≡ 0`` by ``k`` gives

    k·n ≡ q + k·r   (mod p),

and since ``k`` is invertible, ``p ∣ n`` **iff** ``p ∣ (q + k·r)``. The map
``n ↦ q + k·r`` therefore preserves divisibility by ``p`` while shrinking the number, so
iterating it reduces any ``n`` to something small enough to check directly. The familiar
"×2 and subtract for 7" rule is exactly this with ``p = 7`` (where ``10⁻¹ ≡ 5 ≡ -2``).

The cases ``p = 2`` and ``p = 5`` are handled directly, since 10 is not invertible there.

- :func:`divisibility_step` — one application of ``n ↦ ⌊n/10⌋ + (10⁻¹ mod p)·r``.
- :func:`divisible_by_p` — iterate to a fixed range, then read off the answer.
- :func:`divisibility_report` — apply the test for every prime up to 31 and print a report.
"""

from sage.all import *

def divisibility_step(n, p):
    """
    Apply one step of the divisibility transformation:
        n -> floor(n/10) + k * (last digit)
    where k ≡ 10^{-1} mod p
    """
    if p in (2, 5):
        return n  
    
    k = inverse_mod(10, p)
    return (n // 10) + k * (n % 10)


def divisible_by_p(n, p):
    """
    General divisibility test using iterative reduction.
    Works for any prime p (except 2, 5 handled directly).
    """
    if p == 2:
        return n % 2 == 0
    if p == 5:
        return n % 5 == 0

    n = ZZ(n)

    # Reduce until small
    while abs(n) >= p:
        n = divisibility_step(n, p)

    return n % p == 0


def divisibility_report(n):
    """
    Check divisibility for primes up to 31 using Sage.
    """
    primes = [p for p in primes_first_n(11) if p <= 31]

    print(f"\nDivisibility report for {n}:\n")

    for p in primes:
        if divisible_by_p(n, p):
            print(f"{p} divides {n}")
        else:
            print(f"{p} does NOT divide {n}")
