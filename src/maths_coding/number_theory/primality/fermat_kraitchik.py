"""Fermat–Kraitchik factorization — turning a factoring problem into a search for squares.

Fermat's idea: any odd ``N = u·v`` can be written as a **difference of two squares**,

    N = a² - b² = (a - b)(a + b),   with  a = (u+v)/2,  b = (v-u)/2.

So to factor ``N`` it suffices to find ``a`` with ``a² - N`` a perfect square. The search
starts at ``a = ⌈√N⌉`` and increments, testing whether ``a² - N`` is a square at each step;
when it is, ``(a - b, a + b)`` is a nontrivial factorization. This is fast when ``N`` has
two factors close to ``√N`` and slow when they are far apart — the complementary regime to
trial division — and it is the conceptual seed of Kraitchik's method and, ultimately, of
the quadratic sieve.

Note
----
The exponent in the ``max_steps`` default is written ``10^6`` in Sage's ``^``-as-power
convention; evaluate this module through Sage (preparsed) rather than as plain Python.
"""

from sage.all import *

def fermat_kraitchik(N, max_steps=10^6):
    """
    Factor an odd integer N using the Fermat–Kraitchik method.

    Parameters
    ----------
    N : integer
        Number to factor (preferably odd)
    max_steps : int
        Maximum number of iterations

    Returns
    -------
    (a-b, a+b) if factorization found, else None
    """

    N = ZZ(N)

    # Handle trivial cases
    if N % 2 == 0:
        return (2, N // 2)

    # Start at ceil(sqrt(N))
    a = ceil(sqrt(N))

    for step in range(max_steps):
        b2 = a*a - N

        if b2.is_square():
            b = sqrt(b2)
            return (a - b, a + b)

        a += 1

    return None