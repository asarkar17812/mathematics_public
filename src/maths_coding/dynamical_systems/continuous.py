"""Continuous dynamical systems — equilibria and linear stability.

An autonomous system of ODEs ``ẋ = F(x)``, ``x ∈ ℝⁿ``, defines a *flow*. Its skeleton is
the set of **equilibria** (points where ``F = 0``, so the system sits still) together with
the stability of each. The central technique is **linearization**: near an equilibrium
``x*``, writing ``x = x* + δ`` gives ``δ̇ ≈ J δ`` where ``J = DF(x*)`` is the Jacobian. The
eigenvalues of ``J`` then classify the equilibrium:

============================  ==================================================
eigenvalues of ``J``          local behaviour
============================  ==================================================
all real parts < 0            asymptotically **stable** (sink)
all real parts > 0            **unstable** (source)
real parts of mixed sign      **saddle**
nonzero imaginary parts       spiralling / oscillatory approach
some real part = 0            **non-hyperbolic** — linearization is inconclusive
============================  ==================================================

The Hartman–Grobman theorem guarantees that at a *hyperbolic* equilibrium (no eigenvalue on
the imaginary axis) the nonlinear flow really does look like its linearization.

Functions take ``F`` as a list/vector of Sage expressions in the variables ``vars``. The
worked example is the **damped pendulum**, whose downward rest state is a stable spiral and
whose inverted state is a saddle.
"""

from sage.all import *


def fixed_points(F, variables):
    """
    Solve ``F(x) = 0`` for the equilibria of the system.

    Parameters
    ----------
    F : list/vector of symbolic expressions
        The right-hand side ``F`` of ``ẋ = F(x)``.
    variables : list of sage variables
        The state variables.

    Returns
    -------
    list of dict
        Each dict maps every state variable to a coordinate of one equilibrium.
    """
    eqns = [Fi == 0 for Fi in F]
    return solve(eqns, list(variables), solution_dict=True)


def jacobian_matrix(F, variables):
    """
    Return the symbolic Jacobian ``DF`` of the vector field, ``J_{ij} = ∂Fᵢ/∂xⱼ``.
    """
    return jacobian(vector(SR, list(F)), list(variables))


def linearize(F, variables, point):
    """
    Evaluate the Jacobian of ``F`` at an equilibrium ``point``.

    Parameters
    ----------
    F : list/vector of symbolic expressions
    variables : list of sage variables
    point : dict
        Substitution ``{var: value}`` for the equilibrium.

    Returns
    -------
    sage matrix
        The constant Jacobian ``J = DF(point)`` governing the linearized flow.
    """
    J = jacobian_matrix(F, variables)
    return J.subs(point)


def classify_equilibrium(J, tol=1e-9):
    """
    Classify an equilibrium from its (constant) Jacobian ``J`` via the eigenvalues.

    Parameters
    ----------
    J : sage matrix
        The Jacobian at the equilibrium (numeric, or symbolic that evaluates to numbers).
    tol : float
        Threshold below which a real or imaginary part is treated as zero.

    Returns
    -------
    (label, eigenvalues) : tuple
        A human-readable classification and the list of (complex, numeric) eigenvalues.
        Classic 2-D names (node / spiral / centre / saddle) are used when ``J`` is 2×2.
    """
    eigs = matrix(CDF, J).eigenvalues()
    reps = [e.real() for e in eigs]
    imps = [e.imag() for e in eigs]

    oscillatory = any(abs(im) > tol for im in imps)
    n = J.nrows()

    if any(abs(re) < tol for re in reps):
        if all(abs(re) < tol for re in reps) and oscillatory:
            label = "centre (non-hyperbolic; linearization inconclusive)"
        else:
            label = "non-hyperbolic (a zero real part; linearization inconclusive)"
        return label, eigs

    all_neg = all(re < 0 for re in reps)
    all_pos = all(re > 0 for re in reps)

    if all_neg:
        stab = "stable"
    elif all_pos:
        stab = "unstable"
    else:
        return ("saddle", eigs)

    if n == 2:
        shape = "spiral" if oscillatory else "node"
        return (f"{stab} {shape}", eigs)
    return (f"{stab}" + (" (oscillatory)" if oscillatory else ""), eigs)


def analyze_system(F, variables, verbose=True):
    """
    Find every equilibrium of ``ẋ = F(x)`` and classify it.

    Returns
    -------
    list of (point, label, eigenvalues)
        One entry per equilibrium. Also prints a report when ``verbose``.
    """
    results = []
    for pt in fixed_points(F, variables):
        J = linearize(F, variables, pt)
        try:
            label, eigs = classify_equilibrium(J)
        except (TypeError, ValueError):
            label, eigs = "indeterminate (non-numeric Jacobian)", None
        results.append((pt, label, eigs))
        if verbose:
            print(f"equilibrium {pt}: {label}")
            if eigs is not None:
                print(f"    eigenvalues: {[complex(e) for e in eigs]}")
    return results


if __name__ == "__main__":
    # Worked example 1: damped pendulum  theta'' + b theta' + sin(theta) = 0,
    # written as a first-order system (x = theta, y = theta').
    #   x' = y
    #   y' = -sin(x) - b y
    # The two equilibria are given explicitly (solve() on sin(x)=0 is unreliable):
    # rest state (0, 0) is a stable spiral; inverted state (pi, 0) is a saddle.
    var('x y')
    b = QQ(1) / 2
    F = [y, -sin(x) - b * y]

    print("Damped pendulum  x'=y,  y'=-sin(x)-(1/2)y\n")
    for pt in ({x: 0, y: 0}, {x: pi, y: 0}):
        J = linearize(F, [x, y], pt)
        label, eigs = classify_equilibrium(J)
        print(f"equilibrium {pt}: {label}")
        print(f"    eigenvalues: {[complex(e) for e in eigs]}")

    # Worked example 2: a polynomial competition model, where fixed_points() enumerates
    # every equilibrium automatically.
    #   x' = x(3 - x - 2y),   y' = y(2 - x - y)
    print("\nCompetition model  x'=x(3-x-2y),  y'=y(2-x-y)\n")
    G = [x * (3 - x - 2 * y), y * (2 - x - y)]
    analyze_system(G, [x, y])
