"""p-adic numbers — a different notion of "small", and the arithmetic it unlocks.

Fix a prime ``p``. The **p-adic valuation** ``vₚ(x)`` is the exponent of ``p`` in ``x`` (so
``vₚ(p^k · unit) = k``), and the **p-adic absolute value** is ``|x|ₚ = p^{−vₚ(x)}``. The twist:
a number is *small* p-adically when it is **highly divisible by ``p``** — ``|p^k|ₚ = p^{−k} → 0``.
This metric is **ultrametric**: it satisfies the strong triangle inequality

    |x + y|ₚ ≤ max(|x|ₚ, |y|ₚ),   with equality whenever |x|ₚ ≠ |y|ₚ,

which makes p-adic analysis startlingly clean — for instance a series ``Σ aₙ`` converges **iff**
its terms go to zero, ``|aₙ|ₚ → 0``. So ``Σ p^n`` converges, and in fact ``Σ_{n≥0} p^n = 1/(1−p)``
in ``ℚₚ``; for ``p = 2`` that reads ``1 + 2 + 4 + 8 + ⋯ = −1`` in ``ℚ₂``. Every element of ``ℤₚ``
has a unique expansion ``Σ dₙ pⁿ`` with digits ``dₙ ∈ {0, …, p−1}`` — for a negative integer the
expansion is infinite (``−1 = (p−1)(p−1)(p−1)…``), the p-adic shadow of two's complement.

This is the metric completion the Hensel tower lives in (cf.
:mod:`maths_coding.number_theory.utils.hensel_wilson`): a simple root mod ``p`` lifts to a
genuine root in ``ℤₚ``. And **Ostrowski's theorem** says these are not exotic curiosities —
*every* nontrivial absolute value on ``ℚ`` is equivalent to either the usual real one or one of
the ``|·|ₚ``, so ``ℝ`` and the ``ℚₚ`` are the complete siblings of ``ℚ``.

This module is pure-Python (exact via :class:`fractions.Fraction`) and runs without SageMath:

- :func:`valuation`, :func:`p_adic_abs`, :func:`p_adic_distance` — the valuation and metric.
- :func:`p_adic_digits`, :func:`p_adic_expansion` — the base-``p`` digits of a p-adic integer.
- :func:`ultrametric_check` — the strong triangle inequality (with the equality case).
- :func:`geometric_series` — ``Σ p^n → 1/(1−p)`` (e.g. ``Σ 2^n = −1`` in ``ℚ₂``).
- :func:`series_converges` — convergence ⇔ ``|aₙ|ₚ → 0`` (e.g. ``Σ n!`` converges in every ``ℚₚ``).
- :func:`ostrowski_statement` — the classification of absolute values on ``ℚ``.
"""

from fractions import Fraction
from math import inf, factorial


def valuation(x, p):
    """
    The p-adic valuation ``vₚ(x)`` of a rational ``x`` (``+∞`` for ``x = 0``).

    ``vₚ(a/b) = vₚ(a) − vₚ(b)``: the net power of ``p`` dividing ``x``.
    """
    x = Fraction(x)
    if x == 0:
        return inf
    num, den, v = x.numerator, x.denominator, 0
    while num % p == 0:
        num //= p
        v += 1
    while den % p == 0:
        den //= p
        v -= 1
    return v


def p_adic_abs(x, p):
    """The p-adic absolute value ``|x|ₚ = p^{−vₚ(x)}`` (``0`` for ``x = 0``)."""
    v = valuation(x, p)
    return 0.0 if v == inf else float(p) ** (-v)


def p_adic_distance(x, y, p):
    """The p-adic distance ``|x − y|ₚ`` — small means "agree to a high power of ``p``"."""
    return p_adic_abs(Fraction(x) - Fraction(y), p)


def p_adic_digits(x, p, num=12):
    """
    The first ``num`` base-``p`` digits ``d₀, d₁, …`` of a p-adic integer ``x`` (``vₚ(x) ≥ 0``).

    Extracted Hensel-style: ``d₀ = x mod p`` (interpreting ``x`` in ``ℤₚ``), then recurse on
    ``(x − d₀)/p``. For a negative integer or a unit fraction the expansion is genuinely
    infinite. Raises ``ValueError`` if ``|x|ₚ > 1`` (not a p-adic integer).
    """
    x = Fraction(x)
    if valuation(x, p) < 0:
        raise ValueError("not a p-adic integer (|x|_p > 1); factor out a power of p first")
    digits = []
    for _ in range(num):
        den_inv = pow(x.denominator % p, -1, p)        # p ∤ denominator since v_p(x) ≥ 0
        d = (x.numerator * den_inv) % p
        digits.append(d)
        x = (x - d) / p
    return digits


def p_adic_expansion(x, p, num=12):
    """A readable string ``d₀ + d₁·p + d₂·p² + ⋯`` for the p-adic integer ``x``."""
    digits = p_adic_digits(x, p, num)
    terms = [f"{d}·{p}^{i}" if i else f"{d}" for i, d in enumerate(digits)]
    return " + ".join(terms) + " + …"


def ultrametric_check(x, y, p):
    """
    Verify the strong triangle inequality ``|x+y|ₚ ≤ max(|x|ₚ, |y|ₚ)``, and the equality case.

    Returns
    -------
    dict with ``abs_sum``, ``max_abs``, ``holds`` (the inequality), and ``equality``
    (``|x+y|ₚ = max`` — guaranteed by the theorem whenever ``|x|ₚ ≠ |y|ₚ``).
    """
    ax, ay = p_adic_abs(x, p), p_adic_abs(y, p)
    asum = p_adic_abs(Fraction(x) + Fraction(y), p)
    m = max(ax, ay)
    return {
        "abs_sum": asum,
        "max_abs": m,
        "holds": asum <= m + 1e-12,
        "equality": abs(asum - m) < 1e-12,
    }


def geometric_series(p, num_terms=12):
    """
    Partial sums of ``Σ_{n≥0} p^n`` converging p-adically to ``1/(1−p)``.

    The limit is exact: ``S_N − 1/(1−p) = p^{N+1}/(p−1)``, whose p-adic size ``p^{−(N+1)} → 0``.
    For ``p = 2`` the limit is ``1/(1−2) = −1``, so ``1 + 2 + 4 + ⋯ = −1`` in ``ℚ₂``.

    Returns
    -------
    dict with ``limit`` (``Fraction``), and ``rows`` = list of ``(N, partial_sum, |S_N − limit|ₚ)``.
    """
    limit = Fraction(1, 1 - p)
    rows, S = [], Fraction(0)
    for n in range(num_terms):
        S += Fraction(p) ** n
        rows.append((n, S, p_adic_distance(S, limit, p)))
    return {"limit": limit, "rows": rows}


def series_converges(terms, p):
    """
    Test p-adic convergence of ``Σ aₙ`` by the criterion ``|aₙ|ₚ → 0`` (i.e. ``vₚ(aₙ) → ∞``).

    Unlike over ``ℝ``, this term-test is *exact*: a p-adic series converges iff its terms shrink.

    Returns
    -------
    dict with ``valuations`` (``vₚ(aₙ)``) and ``converges`` (whether they are increasing to ∞).
    """
    vals = [valuation(a, p) for a in terms]
    finite = [v for v in vals if v != inf]
    converges = all(finite[i] <= finite[i + 1] for i in range(len(finite) - 1)) and len(finite) > 1
    return {"valuations": vals, "converges": converges}


def ostrowski_statement():
    """Ostrowski's theorem: the classification of the absolute values on ``ℚ`` (returns a string)."""
    return ("Every nontrivial absolute value on Q is equivalent to either the real absolute "
            "value |·|_∞ or a p-adic absolute value |·|_p for some prime p. So R and the Q_p "
            "are exactly the completions of Q.")


if __name__ == "__main__":
    print("=" * 70)
    print("p-adic numbers — 'small' means 'highly divisible by p'")
    print("=" * 70)

    p = 5
    print(f"\nValuation and absolute value at p = {p}:")
    for x in [Fraction(50), Fraction(3), Fraction(1, 25), Fraction(7, 10)]:
        print(f"  x = {str(x):>5}:  v_{p}(x) = {valuation(x, p)!s:>3}   |x|_{p} = {p_adic_abs(x, p):g}")

    print("\nBase-p digit expansions in Z_p:")
    print(f"  -1 in Z_2  : digits {p_adic_digits(-1, 2, 10)}   (= …1111, the 2-adic −1)")
    print(f"  1/3 in Z_2 : digits {p_adic_digits(Fraction(1, 3), 2, 10)}")
    print(f"  100 in Z_3 : digits {p_adic_digits(100, 3, 8)}   (100 = 1·1 + 0·3 + 2·9 + 0·27 + 1·81)")

    print("\nUltrametric (strong triangle) inequality at p = 3:")
    for x, y in [(Fraction(9), Fraction(3)), (Fraction(1), Fraction(2))]:
        u = ultrametric_check(x, y, 3)
        print(f"  |{x}+{y}|_3 = {u['abs_sum']:g} ≤ max = {u['max_abs']:g}   "
              f"holds={u['holds']}  equality={u['equality']}")

    print("\nGeometric series  Σ pⁿ → 1/(1−p)  (p = 2 gives 1 + 2 + 4 + ⋯ = −1):")
    gs = geometric_series(2, 8)
    print(f"  2-adic limit = {gs['limit']}")
    for n, S, dist in gs["rows"]:
        print(f"    Σ_0^{n} 2ⁿ = {int(S):>4}   |S_N − (−1)|_2 = {dist:g}")

    print("\nConvergence ⇔ terms → 0:  Σ n! converges in every ℚₚ (it never does in ℝ):")
    facts = [factorial(n) for n in range(1, 9)]
    sc = series_converges(facts, 3)
    print(f"  v_3(n!) for n=1..8 = {sc['valuations']}  → increasing, so Σ n! converges 3-adically: {sc['converges']}")

    print("\nOstrowski:", ostrowski_statement())
