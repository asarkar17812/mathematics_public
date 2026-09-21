"""Singular perturbation — when setting ε = 0 changes everything.

**This is a demonstration script.** Run it with Sage:

    sage src/maths_coding/perturbation_theory/singular.py

A perturbation is *singular* when the ``ε = 0`` problem is fundamentally different from the
``ε ≠ 0`` one — so a naive power series in ``ε`` cannot capture the answer. The cure is to
recognise where the expansion breaks down and to rescale there. This script walks through
the three classic faces of the phenomenon, each self-checked against the exact solution.

1. **A lost root.** ``ε x² + x − 1 = 0`` has *two* roots, but at ``ε = 0`` it degenerates to
   the *linear* ``x − 1 = 0`` with only one. The regular series finds the root near ``x = 1``;
   the second root races off like ``x ∼ −1/ε`` and is recovered only by rescaling
   ``x = X/ε`` (the *dominant balance* that keeps the ``εx²`` term in play).

2. **A boundary layer.** For ``ε y'' + y' + y = 0`` on ``[0, 1]`` with ``y(0)=0, y(1)=1``,
   dropping ``ε`` lowers the order of the ODE from two to one, so one boundary condition must
   be abandoned. The *outer* solution ``y ≈ e^{1−x}`` satisfies the condition at ``x=1``; a
   thin **boundary layer** of width ``ε`` near ``x=0`` (inner variable ``X=x/ε``) restores
   ``y(0)=0``. The uniform **composite** approximation is

       y(x) ≈ e^{1−x} − e · e^{−x/ε}.

3. **Multiple scales.** The weakly damped oscillator ``y'' + 2ε y' + y = 0`` (``y(0)=1,
   y'(0)=0``) decays on the *slow* time ``T = εt`` while oscillating on the *fast* time ``t``.
   A plain expansion produces *secular* terms ``∝ t`` that wrongly blow up; treating ``t`` and
   ``T`` as independent and killing the secular terms gives the correct slow decay

       y(t) ≈ e^{−εt} cos t.
"""

import numpy as np

from sage.all import *


def singular_quadratic_roots(eps_value):
    """
    Both roots of ``ε x² + x − 1 = 0`` at a numeric ``ε``, with their singular asymptotics.

    Returns
    -------
    dict
        ``exact`` (the two exact roots), ``regular_asymptotic`` (``1 − ε``, the root that
        survives ``ε → 0``) and ``singular_asymptotic`` (``−1/ε − 1``, the root that escapes).
    """
    e = float(eps_value)
    roots = sorted(np.roots([e, 1.0, -1.0]))
    return {
        "exact": roots,
        "regular_asymptotic": 1 - e,         # root near x = 1
        "singular_asymptotic": -1 / e - 1,   # root near x = -1/eps
    }


def boundary_layer_composite(x, eps):
    """Uniform composite approximation ``e^{1−x} − e·e^{−x/ε}`` for ``ε y'' + y' + y = 0``."""
    return np.exp(1 - x) - np.e * np.exp(-x / eps)


def _boundary_layer_exact(x, eps):
    """Exact solution of ``ε y'' + y' + y = 0``, ``y(0)=0, y(1)=1`` (for verification)."""
    s = np.sqrt(1 - 4 * eps)
    l1, l2 = (-1 + s) / (2 * eps), (-1 - s) / (2 * eps)
    c1, c2 = np.linalg.solve(np.array([[1, 1], [np.exp(l1), np.exp(l2)]]),
                             np.array([0.0, 1.0]))
    return c1 * np.exp(l1 * x) + c2 * np.exp(l2 * x)


def multiple_scales_approx(t, eps):
    """Leading multiple-scales approximation ``e^{−εt} cos t`` for ``y'' + 2ε y' + y = 0``."""
    return np.exp(-eps * t) * np.cos(t)


def _damped_oscillator_exact(t, eps):
    """Exact solution of ``y'' + 2ε y' + y = 0``, ``y(0)=1, y'(0)=0`` (for verification)."""
    w = np.sqrt(1 - eps ** 2)
    return np.exp(-eps * t) * (np.cos(w * t) + (eps / w) * np.sin(w * t))


if __name__ == "__main__":
    print("1. Lost root of  eps x^2 + x - 1 = 0")
    for eps in (0.1, 0.01):
        d = singular_quadratic_roots(eps)
        print(f"  eps={eps}: exact roots = [{d['exact'][0]:.4f}, {d['exact'][1]:.4f}]")
        print(f"            regular  ~ 1 - eps      = {d['regular_asymptotic']:.4f}")
        print(f"            singular ~ -1/eps - 1   = {d['singular_asymptotic']:.4f}")

    print("\n2. Boundary layer of  eps y'' + y' + y = 0,  y(0)=0, y(1)=1   (eps=0.02)")
    eps = 0.02
    print("     x      composite     exact")
    for x in (0.0, 0.01, 0.05, 0.2, 0.5, 1.0):
        print(f"   {x:5.2f}   {boundary_layer_composite(x, eps):+9.4f}   "
              f"{_boundary_layer_exact(x, eps):+9.4f}")

    print("\n3. Multiple scales for  y'' + 2 eps y' + y = 0,  y(0)=1, y'(0)=0   (eps=0.05)")
    eps = 0.05
    print("     t       e^{-eps t} cos t    exact")
    for t in (0.0, 5.0, 10.0, 20.0, 30.0):
        print(f"   {t:5.1f}      {multiple_scales_approx(t, eps):+9.4f}        "
              f"{_damped_oscillator_exact(t, eps):+9.4f}")
