"""Linear algebra — the geometry of a matrix.

A companion in code to the Atlas's Linear Algebra section. The thread is that a single matrix
``A`` carries a complete geometric story, told twice:

- :mod:`four_subspaces` — Strang's "big picture": ``A`` splits its domain into the **row
  space** and **null space** and its codomain into the **column space** and **left null
  space**, mapping the row space *isomorphically* onto the column space and annihilating the
  null space. Rank–nullity and the orthogonality of the pairs are the content.
- :mod:`decompositions` — the three canonical factorizations: the **spectral theorem** for a
  symmetric matrix (``A = QΛQᵀ``, orthonormal eigenbasis), the **singular value
  decomposition** (``A = UΣVᵀ``, which sends the unit sphere to an ellipse with semi-axes the
  singular values), and the **Jordan form** for a defective matrix that no eigenbasis can
  diagonalize.
- :mod:`forms` — bilinear and quadratic forms: **Sylvester's law of inertia** and the
  signature ``(p, n, z)`` (the congruence invariant a form actually sees), positive-
  definiteness, and the determinant as a Gram volume.
- :mod:`exterior_algebra` — tensor and exterior powers: the wedge product (antisymmetry, Plücker
  coordinates), the **determinant as the top exterior power** ``Λⁿ``, and the ``ℝ³`` wedge ↔
  cross-product correspondence.

These reconnect to the rest of the package through the recurring **operator + eigenstructure**
motif: the same eigen/singular machinery is the Jacobian's spectrum in ``dynamical_systems``,
the curvature operators in ``differential_geometry``, and the Fisher metric in
``information_geometry``.
"""
