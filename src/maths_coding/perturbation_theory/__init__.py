"""Perturbation theory — solving hard problems as small deviations from easy ones.

Most equations cannot be solved in closed form. Perturbation theory is the art of solving
the ones that sit *near* a solvable problem: introduce a small parameter ``ε``, expand the
unknown as a power series ``x(ε) = x₀ + x₁ε + x₂ε² + …`` around the ``ε = 0`` answer, and
determine the corrections order by order. It is the workhorse behind quantum mechanics,
celestial mechanics, fluid boundary layers, and much of applied analysis.

The subject splits along a sharp line:

- **Regular** perturbation — the series in ``ε`` converges to the right answer and setting
  ``ε = 0`` recovers a problem of the *same type*. Substitute and match powers.
  (:mod:`regular`, :mod:`eigenvalue`)
- **Singular** perturbation — setting ``ε = 0`` *changes the character* of the problem (a
  root escapes to infinity, the order of a differential equation drops, a boundary
  condition can no longer be met). Naive expansion fails; you need rescaling, boundary
  layers, or multiple time scales. (:mod:`singular`)

Modules
-------
regular
    Regular perturbation series for the roots of an ``ε``-dependent algebraic equation,
    computed by matching powers of ``ε``.
eigenvalue
    **Rayleigh–Schrödinger** perturbation theory for matrix eigenvalues and eigenvectors
    of ``H(ε) = H₀ + εV`` — the linear-algebra heart of quantum perturbation theory.
singular
    The singular world: roots lost at ``ε = 0`` and recovered by rescaling, boundary-layer
    analysis of a stiff ODE, and the method of multiple scales for a weakly damped
    oscillator (a demonstration script).

The natural companion is :mod:`maths_coding.dynamical_systems`: perturbation theory asks how
a dynamical system responds when its defining equations are nudged.
"""
