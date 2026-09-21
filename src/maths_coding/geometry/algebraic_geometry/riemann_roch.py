"""Riemann–Roch on an elliptic curve — counting functions with prescribed poles.

The Riemann–Roch theorem is the central accounting identity for a smooth projective curve.
For a divisor ``D`` (a formal integer combination of points), let ``ℓ(D)`` be the dimension
of the space ``L(D)`` of rational functions whose poles are bounded by ``D``. The theorem
relates this count to the curve's **genus** ``g`` and a canonical divisor ``K``:

    ℓ(D) - ℓ(K - D) = deg(D) + 1 - g.

The genus is the single topological invariant doing the work — for a curve over ``C`` it is
the number of "holes" of the corresponding surface. Elliptic curves are exactly the genus
``g = 1`` case, where the formula simplifies to ``ℓ(D) - ℓ(K - D) = deg(D)``, and for a
divisor of positive degree ``ℓ(D) = deg(D)``. This module sets up that example concretely:

- :func:`elliptic_curve_data` — build ``E`` and report its genus (always 1).
- :func:`divisor_example` — form a small degree-0 divisor ``D = P - O``.
- :func:`riemann_roch_summary` — print the Riemann–Roch identity specialised to ``E``.
"""

from sage.all import *

def elliptic_curve_data(a, b):
    """
    Build the elliptic curve ``y² = x³ + a x + b`` and return it with its genus.

    Returns
    -------
    (E, g) : tuple
        The Sage ``EllipticCurve`` and its genus (which is 1 for any elliptic curve).
    """
    E = EllipticCurve([a, b])
    g = E.genus()
    return E, g

def divisor_example(E):
    """
    Illustrate a degree-0 divisor ``D = P - O`` on ``E`` and print it.

    Uses the point ``P = (0, 0)`` and the point at infinity ``O`` (the group identity) to
    form ``D = P - O``, a divisor of degree ``1 - 1 = 0``. Divisors of degree 0 modulo
    principal divisors form the group that *is* ``E`` — the bridge between the geometry and
    the group law.

    Note
    ----
    Assumes ``(0, 0)`` lies on ``E`` (true e.g. when ``b = 0``); intended as a minimal
    illustration rather than a general routine.
    """
    P = E(0,0)
    O = E(0)

    # divisor D = P - O
    D = [(P, 1), (O, -1)]

    print("Divisor D = P - O:", D)
    print("Degree of D:", 1 - 1)

def riemann_roch_summary(E):
    """Print the genus of ``E`` and the Riemann–Roch identity it satisfies."""
    g = E.genus()
    print(f"Genus: {g}")
    print("Riemann–Roch predicts:")
    print("l(D) - l(K-D) = deg(D) + 1 - g")