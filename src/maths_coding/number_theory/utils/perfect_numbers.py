"""Perfect numbers — the Euclid–Euler bridge to Mersenne primes.

A number is *perfect* if it equals the sum of its proper divisors: ``6 = 1+2+3``,
``28 = 1+2+4+7+14``. The Euclid–Euler theorem makes the even perfect numbers completely
explicit and ties this module to primality testing:

    n is an even perfect number  ⇔  n = 2^(p-1) · (2ᵖ - 1)  with  2ᵖ - 1 prime.

So even perfect numbers are in exact bijection with Mersenne primes ``M_p = 2ᵖ - 1``, and
finding them reduces to deciding primality of Mersenne numbers — which is precisely what
the Lucas–Lehmer test does efficiently. (Whether any *odd* perfect number exists is a
famous open problem.)

The module gives both routes, which is the point: a fast structural one and a slow
honest one to check it against.

- :func:`perfect_numbers_up_to` — generate even perfect numbers structurally, using
  :func:`lucas_lehmer_test` to find the Mersenne primes.
- :func:`is_perfect_bruteforce` — verify perfection directly by summing proper divisors.
- :func:`aliquot_sum` / :func:`classify_number` — the surrounding picture: the sum of
  proper divisors and the deficient / perfect / abundant trichotomy it induces.
"""

from sage.all import ZZ, is_prime

def mersenne_number(p):
    """Return the Mersenne number M_p = 2^p - 1."""
    return ZZ(2)**p - 1

def lucas_lehmer_test(p):
    """
    Decide whether the Mersenne number M_p = 2^p - 1 is prime.

    Iterate the recurrence s ← s² - 2 (mod M_p) starting from s = 4, a total of
    p - 2 times; M_p is prime iff the final residue is 0. A self-contained copy of
    the test lives here so this module is usable without the primality subpackage;
    see :mod:`maths_coding.number_theory.primality.lucas_lehmer_test` for the fuller
    treatment, including the closed form in Q(√3).

    Parameters
    ----------
    p : int
        Exponent (assumed prime; M_p can only be prime when p is).

    Returns
    -------
    bool
        True iff M_p is prime.
    """
    if p == 2:
        return True

    M = mersenne_number(p)
    s = ZZ(4)

    for _ in range(p - 2):
        s = (s*s - 2) % M

    return s == 0

def perfect_numbers_up_to(N):
    """
    Returns all even perfect numbers <= N.
    """
    N = ZZ(N)
    results = []

    p = 2
    while True:
        if is_prime(p):
            if lucas_lehmer_test(p):
                M = mersenne_number(p)
                perfect = ZZ(2)**(p - 1) * M

                if perfect > N:
                    break

                results.append(perfect)

        p += 1

    return results

def is_perfect_bruteforce(n):
    """
    Check if n is perfect by summing proper divisors.
    """
    n = ZZ(n)
    if n <= 1:
        return False

    return sum(d for d in n.divisors() if d < n) == n


def aliquot_sum(n):
    """
    Return the aliquot sum of n: the sum of its proper divisors (σ(n) - n).

    This is the quantity perfection is defined against — n is perfect exactly when
    its aliquot sum equals n.

    Parameters
    ----------
    n : int
        Positive integer.

    Returns
    -------
    sage.rings.integer.Integer
        Sum of the divisors of n that are strictly less than n.
    """
    n = ZZ(n)
    if n < 1:
        raise ValueError("n must be a positive integer")
    return sum(d for d in n.divisors() if d < n)


def classify_number(n):
    """
    Classify n as 'deficient', 'perfect', or 'abundant'.

    Comparing a number to its aliquot sum sorts every integer into one of three
    classes; perfect numbers are the exact knife-edge between deficient (the typical
    case, e.g. all primes and prime powers) and abundant (e.g. 12, whose proper
    divisors 1+2+3+4+6 = 16 overshoot it).

    Parameters
    ----------
    n : int
        Positive integer.

    Returns
    -------
    str
        One of 'deficient', 'perfect', 'abundant'.
    """
    s = aliquot_sum(n)
    n = ZZ(n)
    if s < n:
        return "deficient"
    if s == n:
        return "perfect"
    return "abundant"


if __name__ == "__main__":
    N = 10**8
    nums = perfect_numbers_up_to(N)
    print("Perfect numbers up to", N)
    for n in nums:
        print(n)