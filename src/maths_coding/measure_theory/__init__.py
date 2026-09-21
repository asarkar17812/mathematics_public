"""Measure theory — "how large is a set?", and the integral built on the answer.

A measure generalizes length / area / volume: a map ``μ : Σ → [0, ∞]`` on a σ-algebra of
"measurable" sets that vanishes on ``∅`` and is countably additive. The payoff is a theory of
integration that slices the **range** instead of the **domain** — Lebesgue rather than
Riemann — which is exactly what lets us integrate functions far too irregular for vertical
strips.

- :mod:`lebesgue` — the machinery: σ-algebra and measure axioms on a finite space, the
  counting and Dirac measures, the **simple-function integral** and its monotone approximation
  ``sₙ ↗ f``, the Lebesgue integral as a range-slicing layer cake, the Riemann-vs-Lebesgue
  contrast on the Dirichlet function, and the pushforward measure.
- :mod:`cantor` — the canonical counterexample: the middle-thirds **Cantor set**, uncountable
  yet of Lebesgue measure zero, with Hausdorff dimension ``log 2 / log 3``; and the **Cantor
  function** (devil's staircase), the CDF of a singular measure.

This connects to :mod:`maths_coding.geometry.information_geometry` (measures are the points of
a statistical manifold) and to :mod:`maths_coding.geometry.algebraic_geometry.localization`
(both ask, in different categories, "what is true locally / on small sets?").
"""
