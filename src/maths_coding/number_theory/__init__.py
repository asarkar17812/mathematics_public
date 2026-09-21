"""Number theory — the package's main tower.

This is the most developed part of ``maths_coding``, built to climb from elementary
arithmetic to the arithmetic of elliptic curves while keeping one idea in view:
*reduction mod p turns hard global questions into finite, computable local ones.*

Subpackages
-----------
primality
    Deciding primality, with a focus on Mersenne numbers (Lucas–Lehmer) and on
    factorization (Fermat–Kraitchik).
elliptic
    Elliptic curves seen simultaneously over Q, R, and F_p — rank, torsion, reduction.
dynamics
    Iterating maps over finite fields and the Frobenius operator: orbits, eigenvalues,
    point counts, the Hasse bound and Sato–Tate angles.
utils
    A toolbox of standalone routines (CRT, quadratic reciprocity, divisibility tests,
    perfect numbers, Hensel/Wilson lifting) plus the dataset generators that populate
    ``data_stores/``.

Modules
-------
p_adic
    The p-adic valuation, absolute value, and digit expansions; the ultrametric inequality;
    p-adic convergence (``Σ pⁿ = 1/(1−p)``); Ostrowski's theorem. The completion the Hensel
    tower lives in. *(Pure Python — runs without SageMath.)*
continued_fractions
    Continued fractions as the Euclidean algorithm, convergents as best rational
    approximations, periodic expansions of ``√n``, and Pell's equation. *(Pure Python.)*
arithmetic_functions
    Multiplicative functions (μ, φ, σ, τ), Dirichlet convolution, and Möbius inversion — the
    classical identities (``φ * 1 = id``, ``μ * 1 = ε``) as ring statements. *(Pure Python.)*

Reading order, roughly: ``utils`` and ``primality`` give the elementary groundwork;
``elliptic`` introduces the central object; ``dynamics`` shows Frobenius as the operator
that controls it; ``p_adic`` is the metric world that ``utils.hensel_wilson`` completes into.
"""
