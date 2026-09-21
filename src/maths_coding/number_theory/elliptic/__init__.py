"""Elliptic curves — one object seen over many fields.

An elliptic curve ``y² = x³ + ax + b`` is a single equation that lives in several worlds
at once:

- over **Q** it is an arithmetic object — a finitely generated abelian group
  ``E(Q) ≅ Z^r ⊕ torsion`` (Mordell–Weil);
- over **R** it is a smooth one-dimensional shape you can plot;
- over **F_p** it is a finite group, and reduction ``E(Q) → E(F_p)`` is a homomorphism.

The guiding principle — and the reason elliptic curves recur across this package — is that
the finite local data ``E(F_p)`` constrains the global structure: torsion away from ``p``
injects into ``E(F_p)``, so counting points mod small primes already pins down the
rational torsion. See :mod:`maths_coding.number_theory.elliptic.elliptic_curves`, which
makes all three views concrete for ``y² = x³ + 2x + 3``.
"""
