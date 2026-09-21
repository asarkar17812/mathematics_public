"""Algebraic geometry — the algebra–geometry dictionary, and Riemann–Roch on curves.

Two complementary views of varieties:

- :mod:`nullstellensatz` — the foundational dictionary itself. The maps ``V`` (zero locus)
  and ``I`` (vanishing ideal) become mutually inverse over an algebraically closed field
  (Hilbert's Nullstellensatz, ``I(V(I)) = √I``), with points ↔ maximal ideals and
  irreducible pieces ↔ prime ideals. It also reads smoothness off the Jacobian (the Zariski
  tangent space) for the parabola, node, and cusp.
- :mod:`spectrum` — the scheme-theoretic upgrade: ``Spec(R)`` with the Zariski topology,
  worked out as the "arithmetic line" ``Spec(ℤ)`` (closed points ``(p)`` + a generic point),
  the functorial ``Spec(ℤ/n) ↪ V((n))``, and the ``ℤ ↔ k[x]`` analogy.
- :mod:`localization` — zooming in by inverting elements: ``ℤ_(p)`` as a discrete valuation
  ring, ``R_f = R[t]/(tf−1)`` realising the distinguished open ``D(f)``, and ``ℤ = ⋂_p ℤ_(p)``.
- :mod:`riemann_roch` — smooth projective curves through their *divisors*. The centrepiece
  is the Riemann–Roch theorem,

      l(D) - l(K - D) = deg(D) + 1 - g,

  relating a curve's genus ``g`` to how many functions have poles bounded by a divisor
  ``D``; it is worked out concretely on an elliptic curve (genus 1).
"""
