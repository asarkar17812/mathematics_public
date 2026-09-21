"""Continued fractions — the Euclidean algorithm as the best way to approximate a number.

A continued fraction ``[a₀; a₁, a₂, …] = a₀ + 1/(a₁ + 1/(a₂ + ⋯))`` is the most economical
encoding of a real number. For a **rational**, the expansion is finite and its partial
quotients are exactly the quotients of the **Euclidean algorithm** on ``(p, q)`` — so the
continued fraction *is* ``gcd`` computation in disguise. For a **quadratic irrational** like
``√n`` the expansion is eventually **periodic** (Lagrange's theorem): ``√2 = [1; \\overline{2}]``,
``√7 = [2; \\overline{1,1,1,4}]``.

Truncating at ``aₖ`` gives the **convergents** ``hₖ/kₖ``, computed by the recurrence
``hₖ = aₖ hₖ₋₁ + hₖ₋₂`` (same for ``k``). These are the **best rational approximations**: no
fraction with a smaller denominator lies closer to ``x``, and ``|x − hₖ/kₖ| < 1/kₖ²``. They also
solve **Pell's equation** ``x² − n y² = 1``: its fundamental solution is the convergent of
``√n`` at the end of a period — the bridge from approximation theory back to Diophantine number
theory (and the unit group of ``ℤ[√n]``).

Pure-Python and exact (via :class:`fractions.Fraction`), so it runs without SageMath:

- :func:`cf_expansion` — the (finite) continued fraction of a rational, via Euclid.
- :func:`cf_sqrt` — the eventually-periodic continued fraction ``(a₀, period)`` of ``√n``.
- :func:`convergents` — the convergents ``hₖ/kₖ`` (best rational approximations).
- :func:`evaluate` — fold a finite continued fraction back to a rational.
- :func:`pell_fundamental` — the fundamental solution of ``x² − n y² = 1`` from ``√n``.
"""

from fractions import Fraction
from math import isqrt


def cf_expansion(x):
    """
    The continued fraction ``[a₀; a₁, …]`` of a rational ``x`` — the Euclidean-algorithm quotients.
    """
    x = Fraction(x)
    a = []
    while True:
        ai = x.numerator // x.denominator      # floor, correct for negatives too
        a.append(ai)
        x -= ai
        if x == 0:
            break
        x = 1 / x
    return a


def cf_sqrt(n):
    """
    The eventually-periodic continued fraction of ``√n`` as ``(a₀, period)``.

    Uses the standard ``(m, d, a)`` recurrence; the period closes when ``aₖ = 2a₀``. Returns
    ``(a₀, [])`` for a perfect square (where ``√n`` is the integer ``a₀``).
    """
    a0 = isqrt(n)
    if a0 * a0 == n:
        return a0, []
    m, d, a = 0, 1, a0
    period = []
    while a != 2 * a0:
        m = d * a - m
        d = (n - m * m) // d
        a = (a0 + m) // d
        period.append(a)
    return a0, period


def convergents(cf):
    """
    The convergents ``hₖ/kₖ`` of a (finite) continued fraction ``cf = [a₀, a₁, …]``.

    Returns a list of ``Fraction``\\ s via ``hₖ = aₖ hₖ₋₁ + hₖ₋₂``, ``kₖ = aₖ kₖ₋₁ + kₖ₋₂``.
    """
    h2, h1 = 0, 1     # h_{-2}, h_{-1}
    k2, k1 = 1, 0     # k_{-2}, k_{-1}
    out = []
    for a in cf:
        h = a * h1 + h2
        k = a * k1 + k2
        out.append(Fraction(h, k))
        h2, h1 = h1, h
        k2, k1 = k1, k
    return out


def evaluate(cf):
    """Fold a finite continued fraction ``[a₀, a₁, …, aₙ]`` back into a single rational."""
    value = Fraction(cf[-1])
    for a in reversed(cf[:-1]):
        value = a + 1 / value
    return value


def pell_fundamental(n):
    """
    The fundamental solution ``(x, y)`` of Pell's equation ``x² − n y² = 1`` from ``√n``.

    Builds convergents of the periodic continued fraction of ``√n`` and returns the first with
    ``x² − n y² = 1``. Returns ``None`` for a perfect square (no nontrivial solution).
    """
    a0, period = cf_sqrt(n)
    if not period:
        return None
    cf = [a0] + period * 2          # two periods always contain the fundamental solution
    for c in convergents(cf):
        x, y = c.numerator, c.denominator
        if y > 0 and x * x - n * y * y == 1:
            return (x, y)
    return None


if __name__ == "__main__":
    print("=" * 70)
    print("Continued fractions — Euclid, best approximations, and Pell")
    print("=" * 70)

    print("\nRationals: the continued fraction is the Euclidean algorithm:")
    for x in [Fraction(415, 93), Fraction(355, 113)]:
        cf = cf_expansion(x)
        print(f"  {x} = {cf}   folds back to {evaluate(cf)}  ({evaluate(cf) == x})")

    print("\nQuadratic irrationals: eventually periodic (Lagrange):")
    for n in [2, 3, 7, 23]:
        a0, period = cf_sqrt(n)
        print(f"  √{n} = [{a0}; {period} repeating]")

    print("\nConvergents of √2 are the best rational approximations (→ 1.414213…):")
    a0, period = cf_sqrt(2)
    cf = [a0] + period * 5
    for c in convergents(cf)[:7]:
        print(f"  {str(c):>9} = {float(c):.8f}   |√2 − ·| = {abs(2**0.5 - float(c)):.2e}")

    print("\nPell's equation x² − n y² = 1 from the convergents of √n:")
    for n in [2, 3, 7, 13, 61]:
        sol = pell_fundamental(n)
        x, y = sol
        print(f"  n = {n:>2}:  fundamental (x, y) = ({x}, {y})   check x²−{n}y² = {x*x - n*y*y}")
