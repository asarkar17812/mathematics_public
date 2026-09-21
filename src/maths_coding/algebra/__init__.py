"""Algebra — groups acting on structures.

The thread running through this subpackage is a *group acting on a set*, and what the
action reveals:

- :mod:`galois_group_and_module` — the Galois group of a polynomial acting on its roots,
  realised concretely as permutations and as a permutation representation (matrices).
- :mod:`frobenius_action` — the Frobenius element ``x ↦ xᵖ`` as a *specific* group
  element, read off from how the polynomial factors mod ``p`` (its cycle type).
- :mod:`rotation_group_polyhedra` — the rotation symmetry group of a polyhedron,
  reconstructed numerically from its vertices, then analysed as an abstract group
  (orbits, stabilizers, conjugacy classes).
- :mod:`group_actions` — the abstract machinery behind those examples: orbits, stabilizers,
  the Orbit–Stabilizer theorem, the class equation, and Burnside/Pólya counting.
- :mod:`sylow` — Sylow theory: prime-power subgroups, the counting ``nₚ ≡ 1 (mod p)``, and the
  normality/nilpotency it forces (the first tools for taking a group apart).
- :mod:`jordan_holder` — composition series and the Jordan–Hölder theorem: a finite group's
  unique "prime factorization" into simple composition factors.
- :mod:`representation_theory` — symmetry realized as matrices: characters, the character
  table, and Schur orthogonality (``#irreducibles = #conjugacy classes``, ``Σ dᵢ² = |G|``).
- :mod:`lie_algebra` — the *infinitesimal* version: ``sl₂``, its ladder operators, the
  highest-weight irreducibles ``V_n``, and the Casimir scalar (Schur's lemma in the smooth world).
- :mod:`root_systems` — the combinatorial skeleton behind Lie theory: simple roots, Cartan
  matrices, the root system, and the Weyl group (``A``–``D`` families plus ``G₂``).
- :mod:`smith_normal_form` — the *module*-theoretic counterpart: Smith normal form over a
  PID, which classifies finitely generated abelian groups as the cokernel of a relation
  matrix.
- :mod:`homological_algebra` — exact sequences and the derived functors ``Ext``/``Tor`` of
  cyclic groups (all ``≅ ℤ/gcd``), with exactness tested through the same Smith normal form.

Together they make one point several ways: a symmetry group, or a finitely generated
module, is understood through how it *permutes*, *presents*, or *represents* something concrete.
"""
