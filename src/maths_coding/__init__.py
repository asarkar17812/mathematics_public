"""maths_coding — computational pure mathematics, written to illustrate ideas.

This package is a collection of small, self-contained "labs" for concepts in number
theory, algebra, and geometry. The emphasis is pedagogical: each module implements the
relevant objects exactly (via SageMath) and, wherever it sharpens intuition, prints or
plots the structure so the underlying theorem becomes visible.

Subpackages
-----------
algebra
    Groups acting on things, and the structure they reveal: Galois groups, Frobenius cycle
    types, rotation groups of polyhedra, orbit–stabilizer/Burnside counting, Sylow theory,
    Jordan–Hölder composition series, representation theory (characters), the Lie algebra sl₂
    and its root systems, homological algebra (Ext/Tor), and Smith normal form.
linear_algebra
    The geometry of a matrix: the four fundamental subspaces, and the spectral / SVD /
    Jordan decompositions.
geometry
    Algebraic geometry (Riemann–Roch, the Nullstellensatz dictionary, Spec & the Zariski
    topology, localization), differential geometry (curves, surfaces, geodesics), and
    information geometry (Fisher metric, divergences).
topology
    The covering-space / fundamental-group dictionary (the Galois correspondence for covers
    of the circle), and simplicial homology (Betti numbers and torsion of surfaces).
measure_theory
    "How large is a set?" and the integral built on it: the Lebesgue construction (slicing the
    range), and the Cantor set / devil's staircase.
analysis
    Harmonic analysis: the discrete Fourier transform, the convolution theorem, and Pontryagin
    duality on finite abelian groups.
number_theory
    The most developed tower: primality testing, elliptic curves over Q / R / F_p,
    arithmetic dynamics, Hensel/Wilson lifting, p-adic numbers, continued fractions / Pell,
    arithmetic functions (Möbius/Dirichlet), and a kit of utilities and dataset generators.
dynamical_systems, perturbation_theory
    Long-term behaviour of evolving systems (stability, bifurcations, chaos) and solving
    hard problems near easy ones (regular, eigenvalue, and singular perturbation).

A recurring theme ties the modules together: *local data (mod p) constrains global
structure*, with the Frobenius operator linking Galois theory, finite-field dynamics,
and elliptic curves; and *operator + eigenstructure*, linking the Jacobian, the Fisher
metric, curvature operators, and the spectral/SVD decompositions. See the README for the
full picture.

Note
----
Most modules import from :mod:`sage.all`, so a working SageMath installation is required for
them. A growing set of self-contained modules are deliberately dependency-light and run on
plain CPython (standard library only): :mod:`measure_theory.lebesgue`,
:mod:`measure_theory.cantor`, :mod:`number_theory.p_adic`,
:mod:`number_theory.continued_fractions`, :mod:`number_theory.arithmetic_functions`,
:mod:`topology.homology`, :mod:`analysis.fourier`, :mod:`algebra.lie_algebra`,
:mod:`algebra.root_systems`, and :mod:`algebra.homological_algebra`.
"""

__version__ = "0.0.0"
