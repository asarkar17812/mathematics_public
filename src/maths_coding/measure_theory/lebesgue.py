"""Measures, simple functions, and the Lebesgue integral — slicing the range, not the domain.

The motivating contrast: **Riemann** chops the *domain* into thin vertical strips and sums
``f(x)·Δx`` — perfect for continuous functions, but it collapses for wild ones. **Lebesgue**
chops the *range* into horizontal layers and asks, for each value level, "how large is the set
of ``x`` where ``f`` takes roughly this value?" — and "how large is a set" is exactly what a
**measure** answers. Slicing horizontally is what integrates functions too irregular for
strips.

The construction proceeds in stages. A **simple function** ``s = Σ cᵢ·1_{Eᵢ}`` integrates to
``∫ s dμ = Σ cᵢ·μ(Eᵢ)``. Every non-negative measurable ``f`` is the increasing pointwise limit
of simple functions ``sₙ ↗ f``, and the Lebesgue integral is the supremum
``∫ f dμ = sup{ ∫ s dμ : 0 ≤ s ≤ f, s simple }``. The decisive test case is the **Dirichlet
function** ``1_ℚ``: the rationals are countable, hence a Lebesgue-null set, so ``∫₀¹ 1_ℚ dλ = 0``
— while Riemann gives upper sum ``1`` and lower sum ``0`` and simply fails to exist.

Everything here is computed concretely: finite measure spaces exactly, and the continuous
Lebesgue integral by the range-slicing **layer cake** ``∫ f dλ = ∫₀^∞ λ({f > t}) dt``.

- :func:`is_sigma_algebra`, :func:`is_measure` — the axioms, checked on a finite space.
- :func:`counting_measure`, :func:`dirac_measure` — the two simplest measures.
- :func:`simple_function_integral` — ``Σ cᵢ μ(Eᵢ)``.
- :func:`dyadic_simple_integral` — the canonical ``sₙ ↗ f`` approximation, integral from below.
- :func:`lebesgue_integral`, :func:`riemann_integral` — range-slicing vs domain-slicing.
- :func:`dirichlet_demo` — Lebesgue ``= 0`` while Riemann is undefined.
- :func:`pushforward` — the image measure ``(f_*μ)(B) = μ(f⁻¹(B))``.
"""

import math
from bisect import bisect_right
from itertools import combinations


def _grid(a, b, n):
    """A uniform grid of ``n`` points across ``[a, b]`` (pure Python, no numpy)."""
    if n == 1:
        return [a]
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


# --------------------------------------------------------------------------- #
#  σ-algebras and measures on a finite space (exact)                          #
# --------------------------------------------------------------------------- #

def is_sigma_algebra(X, collection):
    """
    Test whether ``collection`` is a σ-algebra on the finite set ``X``: contains ``X``, and is
    closed under complement and (finite ⇒ countable) union.
    """
    X = frozenset(X)
    coll = {frozenset(s) for s in collection}
    if X not in coll:
        return False
    if any((X - A) not in coll for A in coll):
        return False
    return all((A | B) in coll for A in coll for B in coll)


def is_measure(mu):
    """
    Test whether ``mu`` (a dict ``{frozenset: value}`` over a finite collection) satisfies the
    measure axioms: ``μ(∅) = 0`` and finite additivity on disjoint sets present in the domain.
    """
    if mu.get(frozenset(), 0) != 0:
        return False
    sets = list(mu.keys())
    for A in sets:
        for B in sets:
            if (A & B) == frozenset() and (A | B) in mu:
                if abs(mu[A | B] - (mu[A] + mu[B])) > 1e-9:
                    return False
    return True


def counting_measure(S):
    """The counting measure ``μ(S) = |S|`` (``∞`` for infinite ``S``) — the basis of discrete probability."""
    S = list(S)
    return len(S)


def dirac_measure(x0):
    """The Dirac measure ``δ_{x₀}``: returns a function ``A ↦ 1 if x₀ ∈ A else 0``."""
    return lambda A: 1 if x0 in A else 0


# --------------------------------------------------------------------------- #
#  The Lebesgue integral                                                      #
# --------------------------------------------------------------------------- #

def simple_function_integral(coeffs, measures):
    """
    Integrate a simple function ``s = Σ cᵢ·1_{Eᵢ}``: ``∫ s dμ = Σ cᵢ·μ(Eᵢ)`` (convention ``0·∞ = 0``).
    """
    total = 0.0
    for c, m in zip(coeffs, measures):
        total += 0.0 if c == 0 else c * m
    return total


def dyadic_simple_integral(f, a, b, n, cap=None, n_grid=20000):
    """
    Integral of the canonical dyadic simple approximation ``sₙ`` of ``f ≥ 0`` over ``[a, b]``.

    ``sₙ(x) = ⌊2ⁿ f(x)⌋ / 2ⁿ`` (optionally capped at ``cap``) is simple, and ``sₙ ↗ f``, so
    ``∫ sₙ`` increases toward ``∫ f`` from below — the integral *is* that supremum. Computed by
    averaging ``sₙ`` on a fine domain grid.
    """
    xs = _grid(a, b, n_grid)
    scale = 2 ** n
    total = 0.0
    for x in xs:
        s = math.floor(scale * float(f(x))) / scale
        if cap is not None:
            s = min(s, cap)
        total += s
    return total / len(xs) * (b - a)


def lebesgue_integral(f, a, b, n_levels=4000, n_grid=20000, max_value=None):
    """
    The Lebesgue integral of ``f ≥ 0`` over ``[a, b]`` by the **layer cake**
    ``∫ f dλ = ∫₀^∞ λ({f > t}) dt`` — slicing the *range* and measuring superlevel sets.

    The superlevel-set measure ``λ({f > t})`` is estimated by the fraction of a fine domain grid
    where ``f > t``, times ``(b − a)``. This is Lebesgue's construction, not Riemann's.
    """
    fx = sorted(float(f(x)) for x in _grid(a, b, n_grid))
    N = len(fx)
    M = max_value if max_value is not None else fx[-1]
    ts = _grid(0, M, n_levels)
    dt = (ts[1] - ts[0]) if n_levels > 1 else M
    total = 0.0
    for t in ts:
        # λ({f > t}) ≈ (#grid points with f > t)/N · (b−a)
        count = N - bisect_right(fx, t)
        total += count / N * (b - a)
    return total * dt


def riemann_integral(f, a, b, n=20000):
    """The Riemann integral of ``f`` over ``[a, b]`` by a midpoint sum — slicing the *domain*."""
    h = (b - a) / n
    total = sum(float(f(a + (i + 0.5) * h)) for i in range(n))
    return total * h


def dirichlet_demo():
    """
    The Dirichlet function ``1_ℚ`` on ``[0, 1]``: **Lebesgue integrable** (value 0) but **not
    Riemann integrable**.

    The rationals are a countable union of points, each of measure 0, so by countable additivity
    ``λ(ℚ ∩ [0,1]) = 0`` and ``∫ 1_ℚ dλ = 0``. Every subinterval, however, contains both a
    rational and an irrational, so the Riemann upper sum is ``1`` and the lower sum is ``0``.
    """
    return {
        "lebesgue_integral": 0,
        "reason": "Q ∩ [0,1] is countable ⇒ Lebesgue measure 0 ⇒ ∫ 1_Q dλ = 0",
        "riemann_upper": 1,
        "riemann_lower": 0,
        "riemann_integrable": False,
    }


def pushforward(domain_weights, f):
    """
    The pushforward (image) measure ``(f_*μ)(y) = μ(f⁻¹(y)) = Σ_{x : f(x)=y} μ(x)``.

    Parameters
    ----------
    domain_weights : dict ``{x: mass}``  — a (finite) measure on the domain.
    f : callable ``x ↦ y``.

    Returns
    -------
    dict ``{y: mass}`` — the measure transported to the codomain.
    """
    out = {}
    for x, w in domain_weights.items():
        y = f(x)
        out[y] = out.get(y, 0) + w
    return out


if __name__ == "__main__":
    print("=" * 70)
    print("Measures and the Lebesgue integral — slice the range, not the domain")
    print("=" * 70)

    # ---- A finite σ-algebra and a measure on it ------------------------------
    X = {1, 2, 3}
    powerset = [frozenset(s) for r in range(4) for s in combinations(X, r)]
    print("\nFinite measure space on X = {1,2,3}:")
    print("  power set is a σ-algebra:", is_sigma_algebra(X, powerset))
    mu = {s: len(s) for s in powerset}     # counting measure
    print("  counting measure is a measure:", is_measure(mu),
          "  (μ(X) =", mu[frozenset(X)], ")")

    # ---- Simple-function integral --------------------------------------------
    val = simple_function_integral([2, 5], [3, 0.5])   # 2·μ(E1)+5·μ(E2)
    print("\nSimple function s = 2·1_{E1} + 5·1_{E2}, μ(E1)=3, μ(E2)=0.5:")
    print("  ∫ s dμ =", val)

    # ---- sₙ ↗ f, integral from below, → ∫ f -----------------------------------
    f = lambda x: x ** 2
    print("\nf(x) = x² on [0,1]; dyadic simple functions sₙ ↗ f (∫ → 1/3 from below):")
    for n in range(1, 9):
        print(f"  n = {n}:  ∫ sₙ dλ ≈ {dyadic_simple_integral(f, 0, 1, n):.5f}")

    # ---- Lebesgue (range-slice) vs Riemann (domain-slice) agree on nice f ----
    leb = lebesgue_integral(f, 0, 1)
    rie = riemann_integral(f, 0, 1)
    print(f"\n∫₀¹ x² :  Lebesgue (layer cake) ≈ {leb:.5f}   Riemann (strips) ≈ {rie:.5f}"
          f"   (exact 1/3 = {1/3:.5f})")

    # ---- Where they part ways: the Dirichlet function ------------------------
    d = dirichlet_demo()
    print("\nDirichlet function 1_Q on [0,1] — the function that separates the two integrals:")
    print(f"  Lebesgue ∫ = {d['lebesgue_integral']}   ({d['reason']})")
    print(f"  Riemann: upper sum = {d['riemann_upper']}, lower sum = {d['riemann_lower']}"
          f"  ⇒ not Riemann integrable")

    # ---- Pushforward measure -------------------------------------------------
    src = {(-2): 0.1, (-1): 0.2, 0: 0.3, 1: 0.2, 2: 0.2}   # a measure on {-2..2}
    pf = pushforward(src, lambda x: x * x)                  # transport along x ↦ x²
    print("\nPushforward of a measure along x ↦ x²:")
    print("  source masses :", src)
    print("  (x²)_*μ        :", pf, " (masses of ±k collapse onto k²)")
