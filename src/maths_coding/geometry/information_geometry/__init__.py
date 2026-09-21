"""Information geometry — differential geometry on spaces of probability distributions.

A *parametric statistical model* ``{p(x; θ) : θ ∈ Θ}`` — the Gaussians, the Bernoullis, an
exponential family — is not just a set; it is a smooth manifold whose points are
probability distributions. Information geometry equips that manifold with a metric and
studies it with exactly the tools of
:mod:`maths_coding.geometry.differential_geometry`.

The bridge is the **Fisher information matrix**

    g_{ij}(θ) = E[ ∂ᵢ log p · ∂ⱼ log p ] = -E[ ∂ᵢ ∂ⱼ log p ],

which is a genuine Riemannian metric on parameter space (Rao's insight). Once it is in
hand, "distance between distributions", curvature, and geodesics all make sense.

Two facts give the subject its flavour, and both are demonstrated here:

1. **KL divergence is the metric, infinitesimally.** The Kullback–Leibler divergence is not
   symmetric and is not a distance, yet its second-order Taylor expansion *is* the Fisher
   metric: ``KL(p_θ ‖ p_{θ+dθ}) ≈ ½ dθᵀ g(θ) dθ``. (:mod:`divergence`)
2. **Curved models are genuinely curved.** Feeding the Gaussian family's Fisher metric
   ``g = diag(1/σ², 2/σ²)`` into the differential-geometry curvature machinery gives a
   *constant* Gaussian curvature ``K = -1/2`` — the statistical manifold of normal
   distributions is a hyperbolic space. (:mod:`fisher_metric`)

Modules
-------
fisher_metric
    Build the Fisher information metric — from the definition (an expectation), from the
    log-partition function of an exponential family, or from a discrete model — plus the
    classic closed-form examples.
divergence
    Kullback–Leibler and Bregman divergences, and the demonstration that the Fisher metric
    is the Hessian of the divergence on the diagonal.
"""
