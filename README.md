# Computational Pure Maths

A computational pure-mathematics package whose goal is **to illustrate ideas, not just
compute answers**. Each module is a small, self-contained "lab" for a concept across number
theory, algebra, geometry, dynamical systems, perturbation theory, and information geometry:
it implements the objects, exposes them through readable functions, and whenever helpful, it prints or plots the structure so you can *see* the theorem at work.

The code utilizes on [SageMath](https://www.sagemath.org/) for exact arithmetic and symbolic
calculus (number fields, finite fields, elliptic curves, Galois groups, polynomial rings,
symbolic differentiation and integration) and on `numpy` / `matplotlib` for numerical
experiments, data generation, and visualization.

---

## The Core Idea Behind the Package:

Most of the modules are different views of one recurring story:

> **Local data constrains global structure, and a single operator — Frobenius — ties the
> views together.**

You will meet the same cast of characters in many disguises:

| Object                        | Where it shows up                                                                      |
| ----------------------------- | -------------------------------------------------------------------------------------- |
| **Groups & operators**  | Galois groups, rotation groups of polyhedra, the group law on elliptic curves          |
| **Eigenvalues**         | Frobenius eigenvalues on `E(F_p)`, the Hasse bound `                                 |
| **Frobenius**           | cycle types in `Gal(f)`, `x ↦ xᵖ` on roots, the trace `a_p` on elliptic curves |
| **Finite fields**       | reductions `mod p`, exact orbit / cycle structure, point counts over `F_{pⁿ}`     |
| **Reduction `mod p`** | divisibility tests, CRT, factorization, local-to-global for torsion                    |

The unifying slogan, stated at the end of [`dynamics.py`](src/maths_coding/number_theory/dynamics/dynamics.py):

```
dynamics = operator + eigenstructure + arithmetic constraints
```

Elliptic curves with their Frobenius action are the deepest instance, and the package
builds up to them through simpler examples.

That same **operator + eigenstructure** motif reaches beyond number theory and is what
links the later, more analytic areas of the package:

| Area                            | The operator, and what its eigenstructure tells you                         |
| ------------------------------- | --------------------------------------------------------------------------- |
| **Dynamical systems**     | the Jacobian at an equilibrium — its eigenvalues decide stability          |
| **Perturbation theory**   | `H₀ + εV` — eigenvalues drift, and the corrections are computable      |
| **Differential geometry** | curvature/shape operators — their eigenvalues are the principal curvatures |
| **Information geometry**  | the Fisher metric — curvature of the manifold of probability distributions |

A second cross-cutting thread connects the geometry: the **metric → curvature** machinery
built in `differential_geometry/geodesics.py` is reused verbatim by
`information_geometry/`, where the Fisher matrix plays the role of the metric and the
Gaussian family turns out to be a hyperbolic plane (`K = −1/2`).

---

## Layout

```
src/maths_coding/
├── algebra/                       # Groups acting on things
│   ├── galois_group_and_module.py     # Galois group: roots, conjugacy classes, perm. representation
│   ├── frobenius_action.py            # Frobenius cycle types vs. factorization mod p
│   ├── rotation_group_polyhedra.py    # Rotation symmetry groups of polyhedra (numerically reconstructed)
│   ├── group_actions.py               # Orbit–Stabilizer, the class equation, Burnside/Pólya counting
│   ├── sylow.py                       # Sylow theorems (nₚ ≡ 1 mod p), normality, nilpotency, Cauchy
│   ├── jordan_holder.py               # Composition series; Jordan–Hölder uniqueness of factors
│   ├── representation_theory.py       # Character tables, Schur orthogonality, Σdᵢ²=|G|, decomposition
│   ├── lie_algebra.py                 # sl₂: ladder operators, irreps Vₙ, Casimir = ½n(n+2)  (pure Python)
│   ├── root_systems.py                # Cartan matrices, root systems, Weyl groups (A–D, G₂)  (pure Python)
│   ├── smith_normal_form.py           # Smith normal form: classify f.g. abelian groups as cokernels
│   └── homological_algebra.py         # Exact sequences; Ext/Tor of cyclic groups = ℤ/gcd  (pure Python)
│
├── linear_algebra/                # The geometry of a matrix
│   ├── four_subspaces.py              # Four fundamental subspaces, rank–nullity, orthogonality
│   ├── decompositions.py             # Spectral theorem, SVD (circle→ellipse), Jordan form
│   ├── forms.py                      # Bilinear/quadratic forms, Sylvester's law, signature
│   └── exterior_algebra.py           # Wedge product, Plücker coords, det as Λⁿ, cross product
│
├── topology/                      # Spaces classified by algebra
│   ├── covering_spaces.py            # Covers of S¹: monodromy, deck groups, Galois correspondence
│   └── homology.py                   # Chain complexes; Betti numbers + torsion of surfaces  (pure Python)
│
├── measure_theory/                # "How large is a set?" + the integral on it  (pure Python)
│   ├── lebesgue.py                   # σ-algebras, simple functions, Lebesgue vs Riemann, pushforward
│   └── cantor.py                     # Cantor set (measure 0, dim log2/log3), the devil's staircase
│
├── analysis/                      # Harmonic analysis  (pure Python)
│   └── fourier.py                    # DFT, convolution theorem, Parseval, Pontryagin duality of ℤ/N
│
├── geometry/
│   ├── algebraic_geometry/
│   │   ├── riemann_roch.py            # Genus, divisors, the Riemann–Roch formula for elliptic curves
│   │   ├── nullstellensatz.py         # V–I dictionary, radical membership, Jacobian smoothness (node/cusp)
│   │   ├── spectrum.py                # Spec(R) + Zariski topology; Spec(Z) as the arithmetic line
│   │   └── localization.py            # Z_(p) as a DVR; R_f = R[t]/(tf−1); Z = ⋂_p Z_(p)
│   ├── differential_geometry/
│   │   ├── curves.py                  # Frenet–Serret: curvature, torsion, the moving frame
│   │   ├── surfaces.py                # Fundamental forms; Gaussian, mean, principal curvatures
│   │   └── geodesics.py               # Metric → Christoffel symbols → geodesics → Riemann curvature
│   └── information_geometry/
│       ├── fisher_metric.py           # Fisher information as a Riemannian metric on a statistical model
│       └── divergence.py              # KL & Bregman divergence; KL's Hessian is the Fisher metric
│
├── dynamical_systems/             # Long-term behaviour of evolving systems
│   ├── continuous.py                  # ODE equilibria, Jacobian, linear-stability classification
│   ├── discrete.py                    # Iterated maps, multipliers, logistic map, Lyapunov exponent
│   └── bifurcations.py                # Normal forms + the logistic period-doubling diagram
│
├── perturbation_theory/           # Solving hard problems near easy ones
│   ├── regular.py                     # Power-series roots of an ε-dependent equation
│   ├── eigenvalue.py                  # Rayleigh–Schrödinger matrix eigenvalue perturbation
│   └── singular.py                    # Lost roots, boundary layers, method of multiple scales
│
└── number_theory/
    ├── p_adic.py                      # p-adic valuation/|·|ₚ, digit expansions, Σpⁿ=1/(1−p)  (pure Python)
    ├── continued_fractions.py        # CF = Euclid, convergents, √n periods, Pell's equation  (pure Python)
    ├── arithmetic_functions.py       # μ, φ, σ, τ; Dirichlet convolution; Möbius inversion  (pure Python)
    ├── primality/
    │   ├── lucas_lehmer_test.py       # Lucas–Lehmer test + closed form in Q(√3)
    │   └── fermat_kraitchik.py        # Fermat–Kraitchik difference-of-squares factoring
    ├── elliptic/
    │   └── elliptic_curves.py         # E over Q vs R vs F_p: rank, torsion, reduction, plot
    ├── dynamics/                      # ARITHMETIC dynamics (cf. top-level dynamical_systems/)
    │   ├── dynamics.py                # Orbits in finite fields, matrix dynamics, Frobenius ↔ eigenvalues
    │   └── frobenius_eigenvalues.py   # Frobenius angles (Sato–Tate), point-count recurrence, Hasse bound
    └── utils/
        ├── primes_generator.py        # Build prime datasets up to N
        ├── prime_factorization.py     # Smallest-prime-factor sieve; ω(n), Ω(n) datasets
        ├── mersenne_generator.py       # Mersenne / Lucas–Lehmer datasets
        ├── perfect_numbers.py         # Even perfect numbers via Mersenne primes
        ├── crt.py                     # Chinese Remainder Theorem (pairwise, coprime, verify)
        ├── divisibility_tests.py      # Iterated divisibility test from 10⁻¹ mod p
        ├── quadratic_reciprocity.py   # Legendre symbol, Euler's criterion, Gauss's lemma, reciprocity, Eisenstein
        └── hensel_wilson.py           # Hensel lifting up the p-power tower; Wilson's theorems; √2 in Z_p via CRT
```

> **Two flavours of "dynamics".** `number_theory/dynamics/` is *arithmetic* dynamics
> (iteration over finite fields, Frobenius); the top-level `dynamical_systems/` is the
> *continuous/real* theory (flows, stability, chaos). They share the operator–eigenvalue
> viewpoint but are otherwise separate.

`data_stores/` holds generated datasets (primes, factorizations, Mersenne tests) as
paired `.csv` (human-readable) and `.npy` (structured NumPy array) files. These are
build artifacts — they are produced by the generators in `number_theory/utils/` and are
git-ignored.

---

## Installing

The main dependency is **SageMath** (most modules start with `from sage.all import ...`).
Sage ships its own Python, so the most reliable setup is to install Sage first and use its
interpreter. A growing set of modules are deliberately dependency-light — `measure_theory/`,
`number_theory/p_adic.py`, `number_theory/continued_fractions.py`, `topology/homology.py`,
`analysis/fourier.py`, and `algebra/lie_algebra.py` use only the Python standard library and
run on plain CPython (no Sage needed).

```bash
# 1. Install SageMath (conda-forge is the easiest cross-platform route)
conda create -n sage -c conda-forge sage python=3.11
conda activate sage

# 2. Install this package (editable) into the Sage environment
pip install -e .
```

`numpy` and `matplotlib` come bundled with Sage; they are also declared in
`pyproject.toml` for completeness.

---

## Quick tour

Run a demo script directly with Sage:

```bash
sage src/maths_coding/number_theory/elliptic/elliptic_curves.py
sage src/maths_coding/number_theory/dynamics/frobenius_eigenvalues.py
```

Or import the libraries in a Sage session:

```python
from maths_coding.number_theory.utils.crt import crt
from maths_coding.number_theory.utils.quadratic_reciprocity import (
    legendre_symbol, quadratic_reciprocity,
)
from maths_coding.number_theory.primality.lucas_lehmer_test import lucas_lehmer_test
from maths_coding.number_theory.utils.perfect_numbers import perfect_numbers_up_to

crt([2, 3, 2], [3, 5, 7])          # (23, 105):  x ≡ 23 mod 105
legendre_symbol(2, 7)              # 1   (2 is a QR mod 7, since 3² = 9 ≡ 2)
quadratic_reciprocity(3, 5)        # ((3/5)(5/3),  (-1)^…,  agree?)
lucas_lehmer_test(13)              # True  (M₁₃ = 8191 is prime)
perfect_numbers_up_to(10**4)       # [6, 28, 496, 8128]
```

Explore the newer areas — geometry, dynamics, and perturbation theory:

```python
from sage.all import var, vector, cos, sin, matrix, SR
from maths_coding.geometry.differential_geometry.curves import curvature, torsion
from maths_coding.geometry.differential_geometry.geodesics import gaussian_curvature_from_metric
from maths_coding.dynamical_systems.discrete import (
    logistic, logistic_derivative, lyapunov_exponent,
)
from maths_coding.perturbation_theory.eigenvalue import perturbed_eigenvalues

var('t a b u v')
helix = vector([a*cos(t), a*sin(t), b*t])
curvature(helix, t)                       # a/(a^2 + b^2)   (constant, as for any helix)
torsion(helix, t)                         # b/(a^2 + b^2)

g = matrix(SR, [[1, 0], [0, sin(u)**2]])  # round-sphere metric
gaussian_curvature_from_metric(g, [u, v]) # 1   (intrinsic curvature, no embedding used)

lyapunov_exponent(lambda x: logistic(4, x),
                  lambda x: logistic_derivative(4, x), x0=0.1)   # ≈ 0.693 = log 2 (chaos)

perturbed_eigenvalues([[1,0,0],[0,2.5,0],[0,0,4]],
                      [[0,0.3,0.2],[0.3,0,0.5],[0.2,0.5,0]], eps=0.1)  # Rayleigh–Schrödinger
```

Run a demonstration script for one of the new topics:

```bash
sage src/maths_coding/dynamical_systems/bifurcations.py      # logistic period-doubling diagram
sage src/maths_coding/perturbation_theory/singular.py        # boundary layers & multiple scales
sage src/maths_coding/geometry/information_geometry/fisher_metric.py   # Gaussian manifold has K = -1/2
```

Regenerate a dataset:

```bash
sage src/maths_coding/number_theory/utils/primes_generator.py   # writes data_stores/primes/
```

---

## Status

The number-theory tower
(primality → elliptic curves → Frobenius dynamics) is the most developed, and the package
now also spans **linear algebra** (four fundamental subspaces, spectral/SVD/Jordan
decompositions, quadratic forms, exterior algebra), **measure theory** (the Lebesgue integral
and the Cantor set), **differential geometry** (curves, surfaces, geodesics), **information
geometry** (Fisher metric, divergences), **dynamical systems** (stability, bifurcations,
chaos), and **perturbation theory** (regular, eigenvalue, and singular). The **algebra** and
**algebraic geometry** corners appear through Smith normal form and the
structure theorem, the Orbit–Stabilizer / Burnside machinery, Sylow theory, Jordan–Hölder
composition series, representation theory (character tables), Hensel + Wilson lifting, p-adic
numbers, continued fractions / Pell, Hilbert's Nullstellensatz, Spec with the Zariski topology,
localization, the Lie algebra `sl₂` with its root systems, and homological algebra (Ext/Tor). The **topology** corner spans the covering-space / fundamental-group dictionary (the Galois
correspondence for covers of the circle) and simplicial homology (Betti numbers and torsion
from Smith normal form), and a new **analysis** corner opens harmonic analysis (the DFT, the
convolution theorem, and Pontryagin duality). One integer Smith-normal-form computation now
recurs across `smith_normal_form`, `homology`, and `homological_algebra` — the same `ℤ/2⊕ℤ/2`
falls out as a cokernel, a torsion homology group, and an `Ext`. These reinforce the
package's two recurring slogans: *local data constrains global structure* (local roots
lifting up a `p`-power tower; invariant factors recording torsion; reduction mod `p` as
evaluation on the arithmetic line `Spec(ℤ)`; `ℤ = ⋂_p ℤ_(p)`) and *operator + eigenstructure*
(the spectral theorem, SVD, and the quadratic-form signature joining the Jacobian, the Fisher
metric, and the curvature operators). The final piece involves *volume and the determinant*, and now runs
from the exterior power `Λⁿ` through the Gram volume to the measure-theoretic length that
underlies the Lebesgue integral.
