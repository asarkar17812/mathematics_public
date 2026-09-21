"""Regular perturbation — power series for the roots of an ε-dependent equation.

The simplest perturbation problem: a single algebraic equation ``f(x, ε) = 0`` whose
solution at ``ε = 0`` is known, and which we want for small ``ε ≠ 0``. The regular ansatz is

    x(ε) = x₀ + x₁ ε + x₂ ε² + ⋯

Substitute it into the equation, expand in ``ε``, and demand that the coefficient of each
power of ``ε`` vanish. The ``ε⁰`` term reproduces the unperturbed equation (fixing ``x₀``);
each higher power then determines the next correction ``xₖ`` *linearly*, because ``xₖ`` first
appears multiplied by ``∂f/∂x`` evaluated at the leading order. So the whole series unwinds
by solving one linear equation per order.

The textbook example is ``x² + εx − 1 = 0``. Its root near ``+1`` has the expansion

    x(ε) = 1 − ½ ε + ⅛ ε² + ⋯,

which :func:`perturbation_series_root` reproduces. (Contrast ``εx² + x − 1 = 0`` in
:mod:`maths_coding.perturbation_theory.singular`, where the *same-looking* equation is
singular because the ``ε`` multiplies the highest power.)
"""

from sage.all import *


def perturbation_series_root(equation, x, eps, x0, order):
    """
    Compute the regular perturbation series of a root of ``equation`` in powers of ``eps``.

    Substitutes ``x = x0 + Σ_{k≥1} cₖ εᵏ`` into ``equation`` (an expression that should equal
    zero), expands in ``eps``, and solves the resulting triangular system order by order for
    the coefficients ``cₖ``.

    Parameters
    ----------
    equation : symbolic expression
        The left-hand side of ``equation = 0``, in the variables ``x`` and ``eps``.
    x : sage variable
        The unknown.
    eps : sage variable
        The small parameter.
    x0 : value
        The leading-order root; must satisfy ``equation|_{eps=0, x=x0} = 0``.
    order : int
        Highest power of ``eps`` to compute.

    Returns
    -------
    (coeffs, series) : tuple
        ``coeffs`` is the list ``[x0, x1, …, x_order]`` and ``series`` is the assembled
        symbolic expression ``Σ coeffs[k] εᵏ``.

    Examples
    --------
    For ``x**2 + eps*x - 1`` with ``x0 = 1`` the coefficients are ``[1, -1/2, 1/8, …]``.
    """
    if order < 0:
        raise ValueError("order must be non-negative")

    cvars = [SR.var(f"_c{k}") for k in range(1, order + 1)]
    coeffs = [SR(x0)] + cvars
    series = sum(coeffs[k] * eps ** k for k in range(order + 1))

    expanded = taylor(equation.subs({x: series}), eps, 0, order)

    known = {}
    for k in range(1, order + 1):
        # Coefficient of eps^k is linear in c_k once c_1..c_{k-1} are substituted.
        ck_eq = expanded.coefficient(eps, k).subs(known)
        sol = solve(ck_eq == 0, cvars[k - 1], solution_dict=True)
        if not sol:
            raise ValueError(f"could not solve for the order-{k} coefficient")
        known[cvars[k - 1]] = sol[0][cvars[k - 1]].simplify_full()

    result = [SR(x0)] + [known[cv] for cv in cvars]
    return result, sum(result[k] * eps ** k for k in range(order + 1))


if __name__ == "__main__":
    var('x epsilon')

    # x^2 + eps x - 1 = 0, the two roots near +1 and -1.
    eq = x ** 2 + epsilon * x - 1
    for x0 in (1, -1):
        coeffs, series = perturbation_series_root(eq, x, epsilon, x0, order=3)
        print(f"root near {x0:>2}:  coeffs = {coeffs}")
        print(f"            x(eps) = {series}")

    # Numerical check against the exact root for a small epsilon.
    coeffs, series = perturbation_series_root(eq, x, epsilon, 1, order=3)
    e0 = 0.05
    approx = series.subs({epsilon: e0})
    exact = ((-epsilon + sqrt(epsilon ** 2 + 4)) / 2).subs({epsilon: e0})
    print(f"\nepsilon = {e0}:  series = {float(approx):.8f}   exact = {float(exact):.8f}")
