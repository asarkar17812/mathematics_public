"""Rotation symmetry groups of polyhedra, reconstructed from their vertices.

This module recovers the *rotation group* of a polyhedron — the symmetries you can realise
by physically rotating it — without being told the answer in advance. The strategy is
discovery-by-numbers, then identification-by-algebra:

1. **Find the rotations numerically.** Any rotation that maps the solid to itself must
   permute its vertices. So we hunt for ``3×3`` matrices ``R`` carrying one ordered triple
   of vertices to another; we keep ``R`` only if it is a genuine rotation
   (``RᵀR = I`` and ``det R = +1``) *and* it permutes the whole vertex set.
   (:func:`find_rotations`)
2. **Name the group abstractly.** Each rotation induces a permutation of the vertices;
   feeding those permutations to Sage yields the abstract group, whose order, orbits,
   stabilizers, and conjugacy classes we then print. (:func:`analyze_polyhedron`)

The classical answers fall out: the tetrahedron gives ``A₄`` (order 12), the cube gives
``S₄`` (order 24), and a regular ``n``-gon gives the cyclic group ``C_n`` of in-plane
rotations. The result is also a small live demonstration of the orbit–stabilizer theorem,
``|G| = |orbit| · |stabilizer|``.

Numerical caveat: vertex matching uses a floating-point tolerance (``atol=1e-6``), which is
fine for these exact-coordinate solids.
"""

from sage.all import *
import numpy as np
import itertools
import math

def regular_ngon(n):
    """Return the ``n`` vertices of a unit regular ``n``-gon in the ``z = 0`` plane."""
    return [
        (math.cos(2*math.pi*k/n), math.sin(2*math.pi*k/n), 0.0)
        for k in range(n)
    ]

def cube_vertices():
    """Return the 8 vertices of the cube ``{-1, 1}³``."""
    return [(x, y, z) for x in (-1,1) for y in (-1,1) for z in (-1,1)]

def tetrahedron_vertices():
    """Return 4 vertices of a regular tetrahedron (alternate corners of the cube)."""
    return [(1,1,1), (-1,-1,1), (-1,1,-1), (1,-1,-1)]

def is_rotation_matrix(R, tol=1e-6):
    """
    Test whether ``R`` is a proper rotation: orthogonal with determinant ``+1``.

    Checks ``RᵀR ≈ I`` (orthogonality, so lengths and angles are preserved) and
    ``det R ≈ +1`` (orientation-preserving, ruling out reflections), each within ``tol``.
    """
    return (
        np.allclose(R.T @ R, np.eye(3), atol=tol)
        and np.isclose(np.linalg.det(R), 1, atol=tol)
    )

def rotation_axis_angle(R):
    """
    Recover the axis and angle of a 3-D rotation matrix ``R``.

    The **axis** is the eigenvector of ``R`` with eigenvalue ``1`` (the direction left
    fixed by the rotation); the **angle** comes from the trace via
    ``trace(R) = 1 + 2cos θ``. Returns ``(axis, angle)`` with ``axis`` a unit vector and
    ``angle`` in radians.
    """
    eigvals, eigvecs = np.linalg.eig(R)

    axis = None
    for i, val in enumerate(eigvals):
        if np.isclose(val, 1, atol=1e-6):
            axis = np.real(eigvecs[:, i])
            axis = axis / np.linalg.norm(axis)
            break

    trace = np.trace(R)
    angle = math.acos(max(min((trace - 1)/2, 1), -1))

    return axis, angle


def find_rotations(vertices):
    """
    Find every rotation symmetry of a polyhedron given its vertices.

    A rotational symmetry must send vertices to vertices, so each one is determined by
    where it sends three independent (non-coplanar) vertices. The routine therefore
    tries every ordered triple → ordered triple map ``R = B·A⁻¹``, keeps those that are
    genuine rotations (:func:`is_rotation_matrix`) and that permute the *entire* vertex
    set, then de-duplicates.

    Parameters
    ----------
    vertices : sequence of 3-tuples
        The polyhedron's vertices (e.g. from :func:`cube_vertices`).

    Returns
    -------
    list[tuple[np.ndarray, tuple[int, ...]]]
        Pairs ``(R, perm)`` of a rotation matrix and the vertex permutation it induces
        (``perm[i] = j`` means vertex ``i`` maps to vertex ``j``).
    """
    V = np.array(vertices)
    n = len(V)
    rotations = []

    for i, j, k in itertools.permutations(range(n), 3):
        A = np.stack([V[i], V[j], V[k]], axis=1)
        if np.linalg.matrix_rank(A) < 3:
            continue

        for i2, j2, k2 in itertools.permutations(range(n), 3):
            B = np.stack([V[i2], V[j2], V[k2]], axis=1)

            try:
                R = B @ np.linalg.inv(A)
            except np.linalg.LinAlgError:
                continue

            if not is_rotation_matrix(R):
                continue

            image = (R @ V.T).T

            perm = [-1]*n
            for idx, pt in enumerate(image):
                for jdx, v in enumerate(V):
                    if np.allclose(pt, v, atol=1e-6):
                        perm[idx] = jdx
                        break

            if -1 not in perm:
                rotations.append((R, tuple(perm)))

    # deduplicate
    unique = []
    for R, p in rotations:
        if not any(np.allclose(R, Q, atol=1e-6) for Q,_ in unique):
            unique.append((R, p))

    return unique

def analyze_polyhedron(vertices):
    """
    Identify and describe the rotation group of a polyhedron.

    Calls :func:`find_rotations`, assembles the induced vertex permutations into a Sage
    ``PermutationGroup``, and prints the group's order and generators, an axis/angle
    classification of each rotation, and the orbits, stabilizer sizes, and conjugacy
    classes. The orbit and stabilizer figures together illustrate the orbit–stabilizer
    theorem.

    Parameters
    ----------
    vertices : sequence of 3-tuples

    Returns
    -------
    None
        Demonstration routine; output is printed.
    """
    data = find_rotations(vertices)

    perms = [Permutation(p) for _, p in data]
    G = PermutationGroup(perms)

    print("Group order:", G.order())
    print("Generators:", G.gens())

    print("\nRotation classification (axis + angle):")
    for R, _ in data:
        axis, angle = rotation_axis_angle(R)
        print("axis =", axis, "| angle =", angle)

    print("\nOrbits on vertices:", G.orbits())

    print("\nStabilizer sizes:")
    for v in range(min(5, len(vertices))):
        print(v, "->", G.stabilizer(v).order())

    print("\nConjugacy classes:")
    for cls in G.conjugacy_classes():
        print("size", len(cls))

