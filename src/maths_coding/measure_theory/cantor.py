"""The Cantor set and the devil's staircase — measure zero, yet uncountable.

The middle-thirds **Cantor set** ``C`` is built by repeatedly deleting open middle thirds:
start with ``[0,1]``, remove ``(1/3, 2/3)``, then the middle third of each surviving interval,
and so on. After ``n`` stages ``2ⁿ`` intervals of length ``3⁻ⁿ`` remain, total length
``(2/3)ⁿ → 0``. The total deleted length is ``Σ 2ⁿ⁻¹/3ⁿ = 1``, so ``C`` has **Lebesgue measure
zero**. Yet ``C`` is **uncountable**: a point lies in ``C`` exactly when it has a base-3
expansion using only the digits ``0`` and ``2``, and the map ``0↦0, 2↦1`` to binary is a
surjection onto ``[0,1]``. So ``C`` is a measure-zero set with the cardinality of the
continuum — the prototype of "small in measure, large in cardinality."

Its **self-similarity dimension** is ``log 2 / log 3 ≈ 0.6309``: ``C`` is two copies of itself
scaled by ``1/3`` (``N = 2`` pieces, ratio ``r = 1/3``, ``dim = log N / log(1/r)``), a
non-integer Hausdorff dimension — a fractal.

The **Cantor function** ``c(x)`` (the *devil's staircase*) is the CDF of the uniform measure on
``C``: continuous and non-decreasing from ``0`` to ``1``, yet **constant on every deleted
interval**, so ``c'(x) = 0`` almost everywhere. It is the standard example of a **singular
measure** (a Lebesgue–Stieltjes measure with no density) — all its "mass" sits on the
measure-zero set ``C``.

- :func:`cantor_intervals` — the ``2ⁿ`` surviving intervals at stage ``n`` (exact rationals).
- :func:`remaining_length` / :func:`removed_length` — ``(2/3)ⁿ`` and its complement (``→ 1``).
- :func:`in_cantor_set` — the base-3 "no digit 1" membership test.
- :func:`hausdorff_dimension` — ``log 2 / log 3``.
- :func:`cantor_function` — the devil's staircase ``c(x)``.
"""

from fractions import Fraction
from math import log


def cantor_intervals(level):
    """The list of ``2^level`` closed intervals (as ``Fraction`` pairs) remaining at stage ``level``."""
    intervals = [(Fraction(0), Fraction(1))]
    for _ in range(level):
        nxt = []
        for lo, hi in intervals:
            third = (hi - lo) / 3
            nxt.append((lo, lo + third))
            nxt.append((hi - third, hi))
        intervals = nxt
    return intervals


def remaining_length(level):
    """Total length of the surviving intervals at stage ``level``: exactly ``(2/3)^level``."""
    return Fraction(2, 3) ** level


def removed_length(level):
    """Total length removed by stage ``level``: ``1 − (2/3)^level`` (``→ 1``, so ``C`` is null)."""
    return 1 - remaining_length(level)


def in_cantor_set(x, depth=80):
    """
    Membership test for the Cantor set via base-3 digits: ``x ∈ C`` iff it has a ternary
    expansion using only ``0`` and ``2``.

    A terminating digit ``1`` is allowed only when nothing follows it (an endpoint such as
    ``1/3 = 0.0222…₃``), which is handled explicitly. ``x`` is coerced to an exact ``Fraction``.
    """
    x = Fraction(x)
    if x < 0 or x > 1:
        return False
    for _ in range(depth):
        if x == 0 or x == 1:
            return True                 # 0 = 0.000…, 1 = 0.222…
        t = x * 3
        d = int(t)
        frac = t - d
        if d == 1:
            return frac == 0            # endpoint 0.…1 = 0.…0222… ∈ C, else in a deleted gap
        x = frac
    return True


def hausdorff_dimension():
    """The self-similarity (Hausdorff) dimension of the Cantor set: ``log 2 / log 3 ≈ 0.6309``."""
    return log(2) / log(3)


def cantor_function(x, depth=80):
    """
    The Cantor function ``c(x)`` — the devil's staircase — evaluated at ``x ∈ [0, 1]``.

    Algorithm: read the base-3 digits of ``x``; stop at the first ``1`` (keeping it); map the
    preceding ``0/2`` digits via ``0↦0, 2↦1``; read the result as a binary fraction. The result
    is the CDF of the uniform measure on the Cantor set: it climbs only on ``C`` and is flat on
    every deleted middle third.
    """
    x = Fraction(x)
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    bits = []
    for _ in range(depth):
        t = x * 3
        d = int(t)
        if d == 1:
            bits.append(1)
            break
        bits.append(d // 2)             # 0 ↦ 0, 2 ↦ 1
        x = t - d
        if x == 0:
            break
    val = Fraction(0)
    for i, b in enumerate(bits):
        val += Fraction(b, 2 ** (i + 1))
    return float(val)


if __name__ == "__main__":
    print("=" * 70)
    print("The Cantor set — measure zero, uncountable, dimension log2/log3")
    print("=" * 70)

    print("\nThe middle-thirds construction (2ⁿ intervals of length 3⁻ⁿ):")
    for n in range(0, 6):
        print(f"  stage {n}: {2**n:>3} intervals, total length (2/3)^{n} = "
              f"{float(remaining_length(n)):.5f}   removed so far = {float(removed_length(n)):.5f}")
    print("  → remaining length → 0, so the Cantor set has Lebesgue measure 0")

    print("\nMembership (x ∈ C iff base-3 expansion avoids the digit 1):")
    tests = [Fraction(0), Fraction(1, 3), Fraction(2, 3), Fraction(1, 4),
             Fraction(1, 2), Fraction(1, 9), Fraction(4, 9), Fraction(1)]
    for x in tests:
        print(f"  {str(x):>4} ∈ C : {in_cantor_set(x)}")

    print(f"\nHausdorff dimension = log2/log3 = {hausdorff_dimension():.6f}  (a fractal: not an integer)")

    print("\nThe Cantor function c(x) (devil's staircase) — flat on every deleted gap:")
    for x in [Fraction(0), Fraction(1, 9), Fraction(1, 3), Fraction(1, 2),
              Fraction(2, 3), Fraction(7, 9), Fraction(1)]:
        flat = "  (on the deleted middle third [1/3,2/3] ⇒ constant 1/2)" if Fraction(1, 3) <= x <= Fraction(2, 3) else ""
        print(f"  c({str(x):>4}) = {cantor_function(x):.5f}{flat}")
    print("  c is continuous and increases 0→1, yet c'(x)=0 a.e.: a singular measure's CDF")
