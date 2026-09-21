"""Bifurcations — how the qualitative picture changes with a parameter.

**This is a demonstration script.** Run it with Sage to see the printout and a plot:

    sage src/maths_coding/dynamical_systems/bifurcations.py

A *bifurcation* is a value of a parameter at which the number or stability of equilibria
suddenly changes — the system reorganises. Three local bifurcations of one-dimensional
flows ``ẋ = f(x; r)`` account for almost everything that can happen generically; they are
captured by **normal forms**:

==================  =====================  =================================================
bifurcation         normal form            what happens as ``r`` increases through 0
==================  =====================  =================================================
saddle-node         ``ẋ = r - x²``         two equilibria ``±√r`` appear out of nothing
transcritical       ``ẋ = r x - x²``       two equilibria cross and *exchange* stability
pitchfork           ``ẋ = r x - x³``       one stable state splits into two (symmetry breaks)
==================  =====================  =================================================

The script first works each normal form out symbolically — equilibria as functions of ``r``
and their stability from the sign of ``∂f/∂x`` — then turns to the **logistic map**
``x ↦ r x(1-x)`` and plots its bifurcation diagram, the iconic period-doubling cascade
(period 1 → 2 → 4 → 8 → …) that accumulates at ``r ≈ 3.5699`` and opens into chaos, with
periodic windows beyond.
"""

import numpy as np
import matplotlib.pyplot as plt

from sage.all import *


def analyze_normal_form(f, x, r, test_values):
    """
    Print the equilibria of ``ẋ = f(x; r)`` and their stability at several ``r`` values.

    For each ``r`` in ``test_values`` the equilibria are found by solving ``f = 0``, and each
    is labelled stable / unstable from the sign of ``∂f/∂x`` there (negative ⇒ stable).

    Parameters
    ----------
    f : symbolic expression in ``x`` and ``r``
    x, r : sage variables
    test_values : iterable
        Parameter values to sample (e.g. ``[-1, 0, 1]``).
    """
    fx = diff(f, x)
    for rv in test_values:
        sols = solve(f.subs({r: rv}) == 0, x, solution_dict=True)
        reals = []
        for s in sols:
            xs = s[x]
            try:
                if abs(complex(xs).imag) < 1e-12:
                    reals.append(xs)
            except TypeError:
                reals.append(xs)
        labels = []
        for xs in reals:
            slope = complex(fx.subs({r: rv, x: xs})).real
            labels.append(f"{xs} ({'stable' if slope < 0 else 'unstable' if slope > 0 else 'marginal'})")
        print(f"  r = {rv:>3}:  equilibria = [{', '.join(labels) if labels else 'none'}]")


def logistic_bifurcation_data(r_min=2.5, r_max=4.0, r_steps=1600,
                              n_iter=600, n_keep=200, x0=0.1):
    """
    Generate ``(r, x)`` points of the logistic-map bifurcation diagram.

    For each ``r`` the map ``x ↦ r x(1-x)`` is iterated ``n_iter`` times from ``x0`` and the
    last ``n_keep`` iterates (the attractor) are recorded.

    Returns
    -------
    (rs, xs) : tuple of numpy arrays
        Flat arrays suitable for a scatter plot.
    """
    r_values = np.linspace(r_min, r_max, r_steps)
    rs, xs = [], []
    for r in r_values:
        x = x0
        for _ in range(n_iter - n_keep):
            x = r * x * (1 - x)
        for _ in range(n_keep):
            x = r * x * (1 - x)
            rs.append(r)
            xs.append(x)
    return np.array(rs), np.array(xs)


if __name__ == "__main__":
    var('x r')

    print("Saddle-node  x' = r - x^2")
    analyze_normal_form(r - x ** 2, x, r, [-1, 0, 1])

    print("\nTranscritical  x' = r x - x^2")
    analyze_normal_form(r * x - x ** 2, x, r, [-1, 0, 1])

    print("\nPitchfork  x' = r x - x^3")
    analyze_normal_form(r * x - x ** 3, x, r, [-1, 0, 1])

    print("\nGenerating logistic bifurcation diagram ...")
    rs, xs = logistic_bifurcation_data()
    plt.figure(figsize=(8, 5))
    plt.plot(rs, xs, ',k', alpha=0.25)
    plt.title("Logistic map bifurcation diagram:  x -> r x (1 - x)")
    plt.xlabel("r")
    plt.ylabel("attractor x")
    plt.show()
