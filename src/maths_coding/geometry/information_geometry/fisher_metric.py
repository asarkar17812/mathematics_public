"""The Fisher information metric — turning a statistical model into a Riemannian manifold.

For a family ``p(x; θ)`` the **Fisher information matrix** is

    g_{ij}(θ) = E[ ∂ᵢ ℓ · ∂ⱼ ℓ ] = -E[ ∂ᵢ ∂ⱼ ℓ ],     ℓ = log p(x; θ),

where ``E`` is the expectation under ``p(·; θ)`` itself. The two forms agree because
``E[∂ᵢ ℓ] = 0`` (the score has mean zero). This module builds ``g`` three ways:

- :func:`fisher_metric_expected` — straight from the definition, integrating the
  (negative) Hessian of the log-likelihood against the density. Exact but relies on Sage
  evaluating the integral.
- :func:`fisher_metric_discrete` — the same definition for a finite-support model, as a
  sum (always evaluates).
- :func:`fisher_from_log_partition` — for an **exponential family**
  ``p(x; θ) = h(x) exp(θ·T(x) - ψ(θ))``, the metric is simply the Hessian of the
  log-partition function ``ψ``: ``g_{ij} = ∂ᵢ ∂ⱼ ψ``. Pure differentiation, no integration —
  the practical workhorse.

A handful of textbook closed forms (:func:`fisher_gaussian`, :func:`fisher_bernoulli`,
:func:`fisher_exponential`, :func:`fisher_poisson`) are provided to check against.

The punchline (see ``__main__``): the Gaussian family carries the metric
``diag(1/σ², 2/σ²)``, whose Gaussian curvature — computed by the differential-geometry
module — is the constant ``-1/2``. Distribution space is hyperbolic.
"""

from sage.all import *


def fisher_metric_expected(log_pdf, params, x, lo, hi):
    """
    Fisher information metric from the definition, for a continuous model.

    Computes ``g_{ij} = -∫ (∂ᵢ ∂ⱼ log p) · p dx`` over ``x ∈ [lo, hi]``, with
    ``p = exp(log_pdf)``.

    Parameters
    ----------
    log_pdf : symbolic expression
        ``log p(x; θ)`` in the data variable ``x`` and the parameters.
    params : list of sage variables
        The parameters ``θ`` (declare ``assume(sigma > 0)`` etc. as needed so the integral
        converges in closed form).
    x : sage variable
        The data variable to integrate out.
    lo, hi : symbolic
        Integration limits (e.g. ``-oo, oo``).

    Returns
    -------
    sage matrix
        The Fisher information matrix ``g(θ)``.

    Notes
    -----
    Relies on Sage closing the integral symbolically; for models where it cannot, use
    :func:`fisher_from_log_partition` (exponential families) or
    :func:`fisher_metric_discrete` (finite support).
    """
    n = len(params)
    p = exp(log_pdf)
    g = matrix(SR, n, n)
    for i in range(n):
        for j in range(n):
            integrand = -diff(log_pdf, params[i], params[j]) * p
            g[i, j] = integrate(integrand, x, lo, hi).simplify_full()
    return g


def fisher_metric_discrete(log_pmf, params, support):
    """
    Fisher information metric for a finite-support (discrete) model.

    Computes ``g_{ij} = -Σ_x (∂ᵢ ∂ⱼ log p) · p`` over the given support. Always evaluates,
    since the expectation is a finite sum.

    Parameters
    ----------
    log_pmf : callable
        ``log_pmf(x)`` returns the symbolic log-probability of outcome ``x``.
    params : list of sage variables
    support : iterable
        The finite set of outcomes ``x``.

    Returns
    -------
    sage matrix
        The Fisher information matrix.
    """
    n = len(params)
    g = matrix(SR, n, n)
    for xval in support:
        lp = log_pmf(xval)
        p = exp(lp)
        for i in range(n):
            for j in range(n):
                g[i, j] += -diff(lp, params[i], params[j]) * p
    return g.apply_map(lambda e: e.simplify_full())


def fisher_from_log_partition(psi, thetas):
    """
    Fisher metric of an exponential family as the Hessian of the log-partition ``ψ(θ)``.

    For ``p(x; θ) = h(x) exp(θ·T(x) - ψ(θ))`` in natural parameters ``θ``, one has
    ``g_{ij} = ∂ᵢ ∂ⱼ ψ`` (and also ``∇ψ = E[T]``, the mean-value parameters). This needs no
    integration and is the most robust route.

    Parameters
    ----------
    psi : symbolic expression
        The log-partition (cumulant) function ``ψ(θ)``.
    thetas : list of sage variables
        The natural parameters.

    Returns
    -------
    sage matrix
        The Hessian ``∂ᵢ ∂ⱼ ψ``.
    """
    n = len(thetas)
    return matrix(SR, n, n,
                  lambda i, j: diff(psi, thetas[i], thetas[j]).simplify_full())


# --------------------------------------------------------------------------
# Closed-form reference examples (derivations in the docstrings)
# --------------------------------------------------------------------------

def fisher_gaussian(mu, sigma):
    """
    Fisher metric of the normal family ``N(μ, σ²)`` in coordinates ``(μ, σ)``.

    Returns ``diag(1/σ², 2/σ²)``. Famous fact: with this metric the parameter half-plane
    ``{(μ, σ) : σ > 0}`` is a space of constant curvature ``-1/2`` — a scaled hyperbolic
    plane.
    """
    return matrix(SR, [[1 / sigma ** 2, 0], [0, 2 / sigma ** 2]])


def fisher_bernoulli(p):
    """Fisher information of ``Bernoulli(p)`` in the mean parameter ``p``: ``1/(p(1-p))``."""
    return SR(1) / (p * (1 - p))


def fisher_exponential(lam):
    """Fisher information of the exponential family ``Exp(λ)`` in rate ``λ``: ``1/λ²``."""
    return SR(1) / lam ** 2


def fisher_poisson(lam):
    """Fisher information of ``Poisson(λ)`` in the mean ``λ``: ``1/λ``."""
    return SR(1) / lam


if __name__ == "__main__":
    # 1. Exponential-family route (pure differentiation): Bernoulli and Poisson.
    var('theta')
    psi_bernoulli = log(1 + exp(theta))          # natural param theta = log(p/(1-p))
    print("Bernoulli Fisher (natural param) =", fisher_from_log_partition(psi_bernoulli, [theta])[0, 0],
          "  [= p(1-p)]")

    psi_poisson = exp(theta)                      # natural param theta = log(lambda)
    print("Poisson   Fisher (natural param) =", fisher_from_log_partition(psi_poisson, [theta])[0, 0],
          "  [= lambda]")

    # 2. Definition route (integral): the Gaussian family.
    var('x mu sigma')
    assume(sigma > 0)
    log_pdf = -log(sigma * sqrt(2 * pi)) - (x - mu) ** 2 / (2 * sigma ** 2)
    g = fisher_metric_expected(log_pdf, [mu, sigma], x, -oo, oo)
    print("\nGaussian Fisher metric g(mu, sigma) =")
    print(g, "   (expected diag(1/sigma^2, 2/sigma^2))")

    # 3. Cross-link to differential geometry: this metric has constant curvature -1/2.
    #    (Requires the package to be importable, i.e. `pip install -e .` or src on the path.)
    try:
        from src.maths_coding.geometry.differential_geometry.geodesics import (
            gaussian_curvature_from_metric,
        )
        K = gaussian_curvature_from_metric(g, [mu, sigma])
        print("\nGaussian curvature of the statistical manifold:", K, "   (expected -1/2)")
    except ImportError:
        print("\n(Install the package to compute the manifold's curvature: it is -1/2.)")
