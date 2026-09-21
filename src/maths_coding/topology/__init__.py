"""Topology — spaces seen through the algebra that classifies them.

Two complementary tools for telling spaces apart:

- :mod:`covering_spaces` — the **covering-space / fundamental-group dictionary**, the
  topological mirror of Galois theory. Covers of the circle ``S¹``, where ``π₁(S¹) = ℤ`` acts on
  each fiber by *monodromy*. The connected covers ``z ↦ zⁿ`` are classified by the subgroups
  ``nℤ ≤ ℤ`` exactly as field extensions are classified by subgroups of a Galois group: the
  degree of the cover is the index of the subgroup, every cover is regular (``ℤ`` is abelian),
  and the deck-transformation group is the quotient ``ℤ/nℤ``.
- :mod:`homology` — **simplicial homology**: Betti numbers and torsion of a chain complex
  computed by integer Smith normal form, worked out for the sphere, torus, projective plane, and
  Klein bottle (where non-orientability appears as ``ℤ/2`` torsion).

The recurring idea — *a symmetry group realized as transformations of a space* — connects this
to the group actions of :mod:`maths_coding.algebra.group_actions` (deck transformations are a
free group action whose quotient is the base) and to the Galois groups of
:mod:`maths_coding.algebra.galois_group_and_module`.
"""
