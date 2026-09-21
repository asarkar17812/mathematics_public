"""Differential geometry — curvature of curves, surfaces, and abstract metrics.

The smooth-geometry counterpart to :mod:`maths_coding.geometry.algebraic_geometry`. The
organising question is *how is a shape curved, intrinsically and extrinsically?* — and the
answer is built up in three stages of increasing abstraction:

Modules
-------
curves
    One-dimensional shapes in space. The **Frenet–Serret** apparatus extracts curvature
    ``κ`` (how fast the curve turns) and torsion ``τ`` (how fast it leaves its osculating
    plane) from a parametrisation — a complete set of bending invariants for a space curve.
surfaces
    Two-dimensional shapes. The **first and second fundamental forms** give the
    **Gaussian**, **mean**, and **principal** curvatures. Gauss's *Theorema Egregium* —
    that ``K`` is intrinsic, computable without leaving the surface — is the headline
    result this module makes tangible.
geodesics
    The intrinsic view. Given only a **metric** ``g``, the **Christoffel symbols** define
    how to differentiate on a curved space, and the **geodesic equations** give its
    "straight lines". This is the abstract Riemannian-geometry language that
    :mod:`maths_coding.geometry.information_geometry` then borrows for statistics.

Everything is done with Sage's symbolic calculus, so the output is exact formulas (e.g.
``K = 1/R²`` for a sphere of radius ``R``) rather than floating-point approximations.
"""
