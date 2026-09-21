"""Surfaces in space — fundamental forms and the curvatures they produce.

A regular parametrised surface ``r(u, v)`` is measured by two quadratic forms:

- the **first fundamental form** ``I = E du² + 2F du dv + G dv²`` records *intrinsic*
  geometry — lengths, angles, and areas as felt by an inhabitant of the surface
  (``E = r_u·r_u``, ``F = r_u·r_v``, ``G = r_v·r_v``);
- the **second fundamental form** ``II = L du² + 2M du dv + N dv²`` records *extrinsic*
  bending — how the surface curves away from its tangent plane, measured against the unit
  normal ``n`` (``L = r_uu·n``, etc.).

From these come the curvatures:

    Gaussian   K = (LN - M²) / (EG - F²)        = κ₁ κ₂
    mean       H = (EN - 2FM + GL) / (2(EG - F²)) = (κ₁ + κ₂)/2
    principal  κ₁, κ₂ = H ± √(H² - K)

The deep fact — Gauss's **Theorema Egregium** — is that ``K``, defined here through the
*extrinsic* second form, is in truth *intrinsic*: it depends only on ``E, F, G``. That is
why a flat sheet of paper (``K = 0``) cannot be wrapped onto a sphere (``K = 1/R² > 0``)
without distortion. The worked example confirms ``K = 1/R²`` for a sphere of radius ``R``.

All computations are symbolic in two Sage variables ``u, v``.
"""

from sage.all import *


def _diff_vec(vec_expr, x):
    """Differentiate a symbolic vector componentwise with respect to ``x``."""
    return vector([diff(c, x) for c in vec_expr])


def _norm(v):
    """Symbolic Euclidean norm ``sqrt(v · v)``."""
    return sqrt(v.dot_product(v))


def first_fundamental_form(r, u, v):
    """
    Return ``(E, F, G)`` of the surface ``r(u, v)`` — the first fundamental form.

    These coefficients encode the induced metric on the surface; ``EG - F² > 0`` is the
    regularity condition (the tangent vectors ``r_u, r_v`` are linearly independent).
    """
    ru = _diff_vec(r, u)
    rv = _diff_vec(r, v)
    E = ru.dot_product(ru).simplify_full()
    F = ru.dot_product(rv).simplify_full()
    G = rv.dot_product(rv).simplify_full()
    return E, F, G


def unit_normal(r, u, v):
    """Return the unit normal ``n = (r_u × r_v) / |r_u × r_v|`` of the surface."""
    ru = _diff_vec(r, u)
    rv = _diff_vec(r, v)
    cross = ru.cross_product(rv)
    # simplify componentwise (simplify_full is a scalar method, not a vector one)
    return vector([c.simplify_full() for c in (cross / _norm(cross))])


def second_fundamental_form(r, u, v):
    """
    Return ``(L, M, N)`` of the surface ``r(u, v)`` — the second fundamental form.

    Computed as the normal components of the second derivatives:
    ``L = r_uu·n``, ``M = r_uv·n``, ``N = r_vv·n``, with ``n`` the unit normal.
    """
    ru = _diff_vec(r, u)
    rv = _diff_vec(r, v)
    ruu = _diff_vec(ru, u)
    ruv = _diff_vec(ru, v)
    rvv = _diff_vec(rv, v)
    n = unit_normal(r, u, v)
    L = ruu.dot_product(n).simplify_full()
    M = ruv.dot_product(n).simplify_full()
    N = rvv.dot_product(n).simplify_full()
    return L, M, N


def gaussian_curvature(r, u, v):
    """
    Gaussian curvature ``K = (LN - M²) / (EG - F²)`` of the surface ``r(u, v)``.

    Positive on sphere-like (elliptic) points, negative on saddle-like (hyperbolic) points,
    zero on flat or cylinder-like (parabolic) points. By the Theorema Egregium this is an
    intrinsic invariant.
    """
    E, F, G = first_fundamental_form(r, u, v)
    L, M, N = second_fundamental_form(r, u, v)
    return ((L * N - M ** 2) / (E * G - F ** 2)).simplify_full()


def mean_curvature(r, u, v):
    """
    Mean curvature ``H = (EN - 2FM + GL) / (2(EG - F²))`` of the surface ``r(u, v)``.

    The average of the principal curvatures. Surfaces with ``H ≡ 0`` are *minimal*
    (soap-film) surfaces; the sign depends on the chosen orientation of the normal.
    """
    E, F, G = first_fundamental_form(r, u, v)
    L, M, N = second_fundamental_form(r, u, v)
    return ((E * N - 2 * F * M + G * L) / (2 * (E * G - F ** 2))).simplify_full()


def principal_curvatures(r, u, v):
    """
    Return the principal curvatures ``(κ₁, κ₂)`` of the surface ``r(u, v)``.

    These are the extreme normal curvatures at a point, obtained as the roots
    ``κ = H ± √(H² - K)`` of ``κ² - 2Hκ + K = 0``.
    """
    K = gaussian_curvature(r, u, v)
    H = mean_curvature(r, u, v)
    disc = sqrt((H ** 2 - K).simplify_full())
    return ((H - disc).simplify_full(), (H + disc).simplify_full())


if __name__ == "__main__":
    # Worked example: sphere of radius R, r(u,v) = R(sin u cos v, sin u sin v, cos u).
    # Classical result: K = 1/R^2 everywhere, and H = -1/R (with the outward-ish normal
    # produced by r_u x r_v in these coordinates). The principal curvatures coincide.
    var('u v R')
    assume(R > 0)
    assume(sin(u) > 0)  # 0 < u < pi, the open sphere
    sphere = vector([R * sin(u) * cos(v), R * sin(u) * sin(v), R * cos(u)])

    print("Sphere r(u,v) =", sphere)
    print("First form  (E,F,G) =", first_fundamental_form(sphere, u, v))
    print("Second form (L,M,N) =", second_fundamental_form(sphere, u, v))
    print("Gaussian K  =", gaussian_curvature(sphere, u, v), "   (expected 1/R^2)")
    print("Mean     H  =", mean_curvature(sphere, u, v), "   (expected -1/R)")
    print("Principal   =", principal_curvatures(sphere, u, v))
