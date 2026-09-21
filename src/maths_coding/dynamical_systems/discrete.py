"""Discrete dynamical systems — iterated maps, stability, and chaos.

A map ``f : ℝ → ℝ`` generates a dynamical system by iteration, ``x_{n+1} = f(x_n)``. The
questions mirror the continuous case but the answers are richer, because even a smooth
one-dimensional map can behave chaotically.

- **Fixed points** ``x*`` satisfy ``f(x*) = x*``. Their stability is governed by the
  **multiplier** ``f'(x*)``: the fixed point attracts when ``|f'(x*)| < 1`` and repels when
  ``|f'(x*)| > 1`` (compare the continuous criterion ``Re λ < 0``).
- **Periodic orbits** of period ``k`` are fixed points of the ``k``-fold composite ``fᵏ``.
- The **Lyapunov exponent**

      λ = lim (1/N) Σ_{n<N} log |f'(x_n)|

  measures the average exponential rate at which nearby trajectories separate. ``λ < 0``
  means orbits are drawn together (regular behaviour); ``λ > 0`` is the defining signature
  of **chaos** — sensitive dependence on initial conditions.

The recurring example is the **logistic map** ``f(x) = r x (1 - x)``, whose Lyapunov
exponent at ``r = 4`` is exactly ``log 2`` (fully chaotic) while at ``r = 3.2`` it is
negative (a stable period-2 cycle). See :mod:`maths_coding.dynamical_systems.bifurcations`
for how the behaviour changes with ``r``.
"""

import math

from sage.all import *


def iterate(f, x0, n):
    """
    Return the orbit ``[x0, f(x0), f²(x0), …, fⁿ(x0)]`` of length ``n + 1``.

    Parameters
    ----------
    f : callable
        The map ``x ↦ f(x)``.
    x0 : numeric
        Initial condition.
    n : int
        Number of iterations.
    """
    orbit = [x0]
    x = x0
    for _ in range(n):
        x = f(x)
        orbit.append(x)
    return orbit


def fixed_points_map(fexpr, x):
    """
    Solve ``f(x) = x`` for the fixed points of a symbolic map.

    Parameters
    ----------
    fexpr : symbolic expression
        The map as an expression in ``x``.
    x : sage variable

    Returns
    -------
    list
        The real/complex solutions of ``f(x) = x``.
    """
    return solve(fexpr == x, x, solution_dict=True)


def multiplier(fexpr, x, xstar):
    """
    Multiplier ``f'(x*)`` of a fixed (or periodic) point.

    A fixed point is linearly **stable** when ``|f'(x*)| < 1`` and **unstable** when
    ``|f'(x*)| > 1``; ``|f'(x*)| = 1`` is the marginal (bifurcation) case.

    Parameters
    ----------
    fexpr : symbolic expression in ``x``
    x : sage variable
    xstar : value to evaluate the derivative at
    """
    return diff(fexpr, x).subs({x: xstar})


def classify_fixed_point(fexpr, x, xstar, tol=1e-12):
    """Return ``'stable'``, ``'unstable'``, or ``'marginal'`` for the fixed point ``xstar``."""
    m = abs(complex(multiplier(fexpr, x, xstar)))
    if m < 1 - tol:
        return "stable"
    if m > 1 + tol:
        return "unstable"
    return "marginal"


def logistic(r, x):
    """The logistic map ``f(x) = r x (1 - x)``."""
    return r * x * (1 - x)


def logistic_derivative(r, x):
    """Derivative of the logistic map, ``f'(x) = r (1 - 2x)``."""
    return r * (1 - 2 * x)


def lyapunov_exponent(f, fprime, x0, n=100000, burn=1000):
    """
    Numerically estimate the Lyapunov exponent ``λ = lim (1/N) Σ log|f'(x_n)|``.

    A transient of ``burn`` iterations is discarded first so the orbit settles onto its
    attractor. ``λ > 0`` indicates chaos; ``λ < 0`` indicates a stable periodic attractor.

    Parameters
    ----------
    f : callable
        The map ``x ↦ f(x)`` (numeric).
    fprime : callable
        Its derivative ``x ↦ f'(x)`` (numeric).
    x0 : float
        Initial condition.
    n : int
        Number of iterations to average over.
    burn : int
        Transient iterations to discard before averaging.

    Returns
    -------
    float
        The estimated Lyapunov exponent.
    """
    x = float(x0)
    for _ in range(burn):
        x = f(x)
    total = 0.0
    for _ in range(n):
        total += math.log(abs(fprime(x)))
        x = f(x)
    return total / n


if __name__ == "__main__":
    # Fixed points of the logistic map for a chosen r, with stability.
    var('x')
    r0 = QQ(32) / 10  # r = 3.2: the nonzero fixed point has lost stability (period-2 regime)
    f = logistic(r0, x)
    print(f"Logistic map f(x) = {r0} x (1 - x)")
    for sol in fixed_points_map(f, x):
        xs = sol[x]
        print(f"  fixed point x* = {xs}  ->  multiplier f'(x*) = {multiplier(f, x, xs)}"
              f"  ({classify_fixed_point(f, x, xs)})")

    # Lyapunov exponent: chaotic at r=4 (= log 2), regular at r=3.2.
    print("\nLyapunov exponents:")
    for r in (3.2, 3.5, 3.83, 4.0):
        lam = lyapunov_exponent(lambda y, r=r: logistic(r, y),
                                lambda y, r=r: logistic_derivative(r, y),
                                x0=0.1)
        tag = "chaos" if lam > 0 else "regular"
        print(f"  r = {r:4}:  lambda = {lam:+.4f}  ({tag})")
    print(f"\n  (exact value at r=4 is log 2 = {math.log(2):.4f})")
