"""Intrinsic geometry from a metric — Christoffel symbols, geodesics, and curvature.

This module takes the *intrinsic* view: forget any surrounding space and start from a
**metric** ``g_{ij}(x)`` — a symmetric positive-definite matrix telling you how to measure
infinitesimal length, ``ds² = Σ g_{ij} dxⁱ dxʲ``. Remarkably, the metric alone determines
all of the geometry:

- the **Christoffel symbols**
  ``Γᵏ_{ij} = ½ gᵏˡ(∂ᵢ g_{jl} + ∂ⱼ g_{il} - ∂ₗ g_{ij})``
  define the Levi-Civita connection — the unique sensible way to differentiate vectors on
  a curved space;
- the **geodesic equations** ``ẍᵏ + Γᵏ_{ij} ẋⁱ ẋʲ = 0`` pick out the "straightest possible"
  curves (great circles on a sphere, straight lines in the plane);
- the **Riemann tensor** and, in two dimensions, the **Gaussian curvature**
  ``K = R_{1212} / det(g)`` measure the failure of the space to be flat — computed here
  *without ever leaving the surface*, the computational face of Gauss's Theorema Egregium.

This is exactly the apparatus :mod:`maths_coding.geometry.information_geometry` reuses, with
the Fisher information matrix playing the role of ``g``. The worked example is the round
unit sphere ``g = diag(1, sin²u)``, for which ``K = 1``.

Conventions: ``g`` is a Sage symbolic ``matrix`` and ``coords`` the matching list of Sage
variables; index order in returned arrays is ``Gamma[k][i][j] = Γᵏ_{ij}``.
"""

from sage.all import *


def christoffel_symbols(g, coords):
    """
    Christoffel symbols ``Γᵏ_{ij}`` of the metric ``g`` in the given coordinates.

    Parameters
    ----------
    g : sage symbolic matrix (n×n, symmetric)
        The metric tensor as a function of ``coords``.
    coords : list of sage variables
        The ``n`` coordinate variables.

    Returns
    -------
    list
        Nested list ``Gamma`` with ``Gamma[k][i][j] = Γᵏ_{ij}`` (symmetric in ``i, j``).
    """
    n = len(coords)
    ginv = g.inverse()
    Gamma = [[[SR(0) for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                s = SR(0)
                for l in range(n):
                    s += ginv[k, l] * (
                        diff(g[j, l], coords[i])
                        + diff(g[i, l], coords[j])
                        - diff(g[i, j], coords[l])
                    )
                Gamma[k][i][j] = (s / 2).simplify_full()
    return Gamma


def geodesic_equations(g, coords, param):
    """
    Return the geodesic equations ``ẍᵏ + Γᵏ_{ij} ẋⁱ ẋʲ = 0`` as symbolic expressions.

    Each coordinate becomes a function of ``param`` (e.g. ``x(t)``), and the returned
    expressions are the left-hand sides that a geodesic must set to zero — ready to feed to
    a symbolic or numerical ODE solver.

    Parameters
    ----------
    g : sage symbolic matrix
    coords : list of sage variables
    param : sage variable
        The curve parameter (arc length / affine parameter), e.g. ``var('t')``.

    Returns
    -------
    (eqns, funcs) : tuple
        ``eqns[k]`` is the geodesic expression for coordinate ``k``; ``funcs[k]`` is the
        corresponding symbolic function ``xᵏ(param)``.
    """
    n = len(coords)
    Gamma = christoffel_symbols(g, coords)

    funcs = [function(f"x{k}")(param) for k in range(n)]
    subs = {coords[k]: funcs[k] for k in range(n)}

    eqns = []
    for k in range(n):
        acc = diff(funcs[k], param, 2)
        for i in range(n):
            for j in range(n):
                acc += Gamma[k][i][j].subs(subs) * diff(funcs[i], param) * diff(funcs[j], param)
        eqns.append(acc.simplify_full())
    return eqns, funcs


def riemann_tensor(g, coords):
    """
    Riemann curvature tensor ``R^ρ_{σμν}`` of the metric ``g``.

    Uses ``R^ρ_{σμν} = ∂_μ Γ^ρ_{νσ} - ∂_ν Γ^ρ_{μσ} + Γ^ρ_{μλ} Γ^λ_{νσ} - Γ^ρ_{νλ} Γ^λ_{μσ}``.
    The tensor vanishes identically iff the space is flat.

    Returns
    -------
    list
        Nested list ``R`` with ``R[rho][sig][mu][nu] = R^ρ_{σμν}``.
    """
    n = len(coords)
    G = christoffel_symbols(g, coords)
    R = [[[[SR(0) for _ in range(n)] for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for rho in range(n):
        for sig in range(n):
            for mu in range(n):
                for nu in range(n):
                    term = diff(G[rho][nu][sig], coords[mu]) - diff(G[rho][mu][sig], coords[nu])
                    for lam in range(n):
                        term += G[rho][mu][lam] * G[lam][nu][sig] - G[rho][nu][lam] * G[lam][mu][sig]
                    R[rho][sig][mu][nu] = term.simplify_full()
    return R


def gaussian_curvature_from_metric(g, coords):
    """
    Intrinsic Gaussian curvature ``K`` of a **2-dimensional** metric ``g``.

    Computes ``K = R_{1212} / det(g)`` with ``R_{1212} = g_{1ρ} R^ρ_{212}`` from the Riemann
    tensor — using only the metric, no embedding. Comparing this with
    :func:`maths_coding.geometry.differential_geometry.surfaces.gaussian_curvature` (which
    uses the extrinsic second fundamental form) is Gauss's Theorema Egregium in action: the
    two agree.

    Parameters
    ----------
    g : sage symbolic 2×2 matrix
    coords : list of two sage variables

    Returns
    -------
    Symbolic expression for ``K``.
    """
    if len(coords) != 2:
        raise ValueError("Gaussian curvature here is defined for 2-dimensional metrics")
    R = riemann_tensor(g, coords)
    R_0101 = sum(g[0, rho] * R[rho][1][0][1] for rho in range(2))
    return (R_0101 / g.det()).simplify_full()


if __name__ == "__main__":
    # Worked example: the round unit sphere, metric g = diag(1, sin^2 u) in (u, v).
    # Great circles are geodesics; the intrinsic Gaussian curvature is +1 everywhere.
    var('u v t')
    assume(sin(u) > 0)
    g = matrix(SR, [[1, 0], [0, sin(u) ** 2]])

    Gamma = christoffel_symbols(g, [u, v])
    print("Sphere metric g = diag(1, sin^2 u)")
    print("Gamma^u_vv =", Gamma[0][1][1], "  (expected -sin(u)cos(u))")
    print("Gamma^v_uv =", Gamma[1][0][1], "  (expected cot u)")
    print("Gaussian K =", gaussian_curvature_from_metric(g, [u, v]), "  (expected 1)")

    eqns, funcs = geodesic_equations(g, [u, v], t)
    print("\nGeodesic equations (each = 0):")
    for f, e in zip(funcs, eqns):
        print(f"  {diff(f, t, 2)} + ... :  {e} = 0")
