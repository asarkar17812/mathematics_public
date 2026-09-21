"""Curves in space — curvature, torsion, and the Frenet–Serret frame.

A regular space curve ``r(t) = (x(t), y(t), z(t))`` carries, at every point, a natural
orthonormal frame and two scalar invariants that together determine it up to rigid motion
(the *fundamental theorem of space curves*):

- **curvature** ``κ`` — how sharply the curve bends, ``κ = |r' × r''| / |r'|³``;
- **torsion** ``τ`` — how sharply it twists out of its osculating plane,
  ``τ = (r' × r'') · r''' / |r' × r''|²``.

The **Frenet–Serret frame** ``(T, N, B)`` — unit tangent, principal normal, binormal —
evolves along the curve according to

    T' = κ N,   N' = -κ T + τ B,   B' = -τ N        (arc-length derivatives).

A straight line has ``κ = 0``; a plane curve has ``τ = 0``; a circular helix has both
constant, which is the worked example at the bottom of the file.

Everything is symbolic: pass in expressions in a Sage variable ``t`` and you get exact
formulas back. The functions assume the curve is *regular* (``r'(t) ≠ 0``) and, for
torsion, *non-degenerate* (``r' × r'' ≠ 0``).
"""

from sage.all import *


def _diff_vec(v, t):
    """Differentiate a symbolic vector componentwise with respect to ``t``."""
    return vector([diff(c, t) for c in v])


def _norm(v):
    """Symbolic Euclidean norm ``sqrt(v · v)`` (avoids the ``abs`` that ``v.norm()`` emits)."""
    return sqrt(v.dot_product(v))


def _simplify_vec(v):
    """Apply ``simplify_full`` componentwise (it is a scalar method, not a vector one)."""
    return vector([c.simplify_full() for c in v])


def speed(r, t):
    """Return the speed ``|r'(t)|`` of the curve."""
    return _norm(_diff_vec(r, t)).simplify_full()


def curvature(r, t):
    """
    Curvature ``κ(t)`` of the space curve ``r(t)``.

    Uses the parametrisation-independent formula ``κ = |r' × r''| / |r'|³``, valid for any
    regular curve (not just unit-speed ones).

    Parameters
    ----------
    r : sage vector of expressions in ``t``
        The curve, e.g. ``vector([cos(t), sin(t), t])``.
    t : sage symbolic variable
        The curve parameter.

    Returns
    -------
    Symbolic expression for ``κ(t)``.
    """
    r1 = _diff_vec(r, t)
    r2 = _diff_vec(r1, t)
    return (_norm(r1.cross_product(r2)) / _norm(r1) ** 3).simplify_full()


def torsion(r, t):
    """
    Torsion ``τ(t)`` of the space curve ``r(t)``.

    Uses ``τ = (r' × r'') · r''' / |r' × r''|²``. Requires ``r' × r'' ≠ 0`` (i.e. the curve
    is not momentarily straight); a planar curve returns ``0``.

    Parameters
    ----------
    r : sage vector of expressions in ``t``
    t : sage symbolic variable

    Returns
    -------
    Symbolic expression for ``τ(t)``.
    """
    r1 = _diff_vec(r, t)
    r2 = _diff_vec(r1, t)
    r3 = _diff_vec(r2, t)
    cross = r1.cross_product(r2)
    return (cross.dot_product(r3) / cross.dot_product(cross)).simplify_full()


def frenet_frame(r, t):
    """
    Return the Frenet–Serret frame ``(T, N, B)`` of ``r(t)`` as symbolic unit vectors.

    ``T`` is the unit tangent ``r'/|r'|``, ``B`` the unit binormal ``(r' × r'')/|r' × r''|``,
    and ``N = B × T`` the principal normal. The three are mutually orthonormal and form a
    right-handed frame at every regular, non-degenerate point.

    Returns
    -------
    (T, N, B) : tuple of sage vectors
    """
    r1 = _diff_vec(r, t)
    r2 = _diff_vec(r1, t)
    T = (r1 / _norm(r1))
    cross = r1.cross_product(r2)
    B = (cross / _norm(cross))
    N = B.cross_product(T)
    return (_simplify_vec(T), _simplify_vec(N), _simplify_vec(B))


def arc_length(r, t, a, b):
    """
    Arc length of ``r(t)`` for ``t`` in ``[a, b]``: ``∫_a^b |r'(t)| dt``.

    Returns a Sage symbolic integral, evaluated in closed form when possible (and otherwise
    left as an ``integrate(...)`` expression that can be evaluated numerically with ``.n()``).
    """
    return integrate(speed(r, t), t, a, b)


if __name__ == "__main__":
    # Worked example: the circular helix r(t) = (a cos t, a sin t, b t).
    # Classical result: kappa = a/(a^2+b^2) and tau = b/(a^2+b^2), both CONSTANT.
    var('t a b')
    assume(a > 0)
    helix = vector([a * cos(t), a * sin(t), b * t])

    print("Helix r(t) =", helix)
    print("speed     =", speed(helix, t))
    print("curvature =", curvature(helix, t), "   (expected a/(a^2+b^2))")
    print("torsion   =", torsion(helix, t), "   (expected b/(a^2+b^2))")

    T, N, B = frenet_frame(helix, t)
    print("\nFrenet frame:")
    print("  T =", T)
    print("  N =", N)
    print("  B =", B)
