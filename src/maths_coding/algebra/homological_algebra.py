"""Homological algebra — exact sequences, and Ext/Tor measuring how exactness fails.

Homological algebra studies **exact sequences** ``⋯ → A → B → C → ⋯`` (image at each step equals
the kernel at the next) and the **derived functors** ``Ext`` and ``Tor`` that quantify how badly
a non-exact situation misses being exact. The starting point (from the Cokernels section of the
notes) is that every map ``φ : A → B`` sits in the four-term exact sequence

    0 → ker φ → A --φ--> B → coker φ → 0,

so ``φ`` is injective iff ``ker φ = 0`` and surjective iff ``coker φ = 0``.

For abelian groups the derived functors are computed from a **free resolution**. A cyclic group
``ℤ/m`` has the length-one resolution ``0 → ℤ --×m--> ℤ → ℤ/m → 0``, and applying
``Hom(−, ℤ/n)`` or ``− ⊗ ℤ/n`` reduces everything to the map ``×m`` on ``ℤ/n``. The result is the
clean, completely computable fact

    Hom(ℤ/m, ℤ/n) ≅ Ext¹(ℤ/m, ℤ/n) ≅ Tor₁(ℤ/m, ℤ/n) ≅ ℤ/gcd(m,n),

where ``Tor₁ = ker(×m)`` is the ``m``-torsion of ``ℤ/n`` and ``Ext¹ = coker(×m) = ℤ/n / m·ℤ/n``.
``Tor`` is also why homology with coefficients and the **universal coefficient theorem** carry a
correction term — the torsion of :mod:`maths_coding.topology.homology` reappears here as a
functor. Exactness of a complex of free abelian groups is just the vanishing of its homology,
which is read off — once again — by the integer **Smith normal form**.

Pure-Python and exact, so it runs without SageMath:

- :func:`hom_cyclic`, :func:`ext1_cyclic`, :func:`tor1_cyclic` — all ``≅ ℤ/gcd(m,n)``, computed
  from the resolution (kernel/cokernel of ``×m`` on ``ℤ/n``).
- :func:`kernel_cokernel_sequence` — the exact sequence ``0→ker→A→B→coker→0`` of a matrix.
- :func:`is_exact` / :func:`homology_ranks` — exactness of a complex of free groups via SNF.
"""

from math import gcd


# --------------------------------------------------------------------------- #
#  Integer Smith normal form (rank + invariant factors), self-contained       #
# --------------------------------------------------------------------------- #

def _invariant_factors(matrix):
    """Nonzero invariant factors of an integer matrix (its Smith normal form diagonal)."""
    A = [row[:] for row in matrix]
    m = len(A)
    n = len(A[0]) if m else 0
    if m == 0 or n == 0:
        return []

    def swap_rows(i, j):
        A[i], A[j] = A[j], A[i]

    def swap_cols(i, j):
        for r in A:
            r[i], r[j] = r[j], r[i]

    factors, t = [], 0
    while t < min(m, n):
        piv = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i][j] != 0 and (piv is None or abs(A[i][j]) < abs(A[piv[0]][piv[1]])):
                    piv = (i, j)
        if piv is None:
            break
        swap_rows(t, piv[0])
        swap_cols(t, piv[1])
        cleared = False
        while not cleared:
            cleared = True
            for i in range(t + 1, m):
                if A[i][t] != 0:
                    q = A[i][t] // A[t][t]
                    for j in range(t, n):
                        A[i][j] -= q * A[t][j]
                    if A[i][t] != 0:
                        swap_rows(t, i)
                        cleared = False
            for j in range(t + 1, n):
                if A[t][j] != 0:
                    q = A[t][j] // A[t][t]
                    for i in range(t, m):
                        A[i][j] -= q * A[i][t]
                    if A[t][j] != 0:
                        swap_cols(t, j)
                        cleared = False
        repaired = False
        for i in range(t + 1, m):
            for j in range(t + 1, n):
                if A[i][j] % A[t][t] != 0:
                    for jj in range(n):
                        A[t][jj] += A[i][jj]
                    repaired = True
                    break
            if repaired:
                break
        if repaired:
            continue
        factors.append(abs(A[t][t]))
        t += 1
    return factors


def _rank(matrix):
    return len(_invariant_factors(matrix))


# --------------------------------------------------------------------------- #
#  Ext, Tor, Hom of cyclic groups via the resolution 0 → Z --×m--> Z → Z/m → 0 #
# --------------------------------------------------------------------------- #

def hom_cyclic(m, n):
    """
    ``Hom(ℤ/m, ℤ/n) ≅ ℤ/gcd(m,n)``: a homomorphism sends the generator to an ``n``-torsion
    element killed by ``m``, i.e. an element of ``ker(×m on ℤ/n)``.

    Returns a dict with the structure (``order``), and the explicit kernel order.
    """
    g = gcd(m, n)
    ker_order = sum(1 for x in range(n) if (m * x) % n == 0)   # |ker(×m on Z/n)|
    return {"group": f"Z/{g}" if g > 1 else "0", "order": g, "kernel_order": ker_order}


def tor1_cyclic(m, n):
    """
    ``Tor₁(ℤ/m, ℤ/n) ≅ ℤ/gcd(m,n)``: the ``m``-torsion ``ker(×m on ℤ/n)`` from tensoring the
    resolution with ``ℤ/n``.
    """
    g = gcd(m, n)
    ker_order = sum(1 for x in range(n) if (m * x) % n == 0)
    return {"group": f"Z/{g}" if g > 1 else "0", "order": g, "kernel_order": ker_order}


def ext1_cyclic(m, n):
    """
    ``Ext¹(ℤ/m, ℤ/n) ≅ ℤ/gcd(m,n)``: the cokernel ``ℤ/n / m·ℤ/n`` of ``×m`` on ``ℤ/n``, from
    applying ``Hom(−, ℤ/n)`` to the resolution.
    """
    g = gcd(m, n)
    image = {(m * x) % n for x in range(n)}                    # image of ×m on Z/n
    coker_order = n // len(image)
    return {"group": f"Z/{g}" if g > 1 else "0", "order": g, "cokernel_order": coker_order}


# --------------------------------------------------------------------------- #
#  Exact sequences of free abelian groups                                     #
# --------------------------------------------------------------------------- #

def kernel_cokernel_sequence(matrix):
    """
    The exact sequence ``0 → ker φ → ℤ^n --φ--> ℤ^m → coker φ → 0`` for an integer matrix ``φ``.

    Returns a dict with ``rank``, ``kernel_rank`` (``n − rank``), ``cokernel_free_rank``
    (``m − rank``), and ``cokernel_torsion`` (invariant factors ``> 1``) — so ``φ`` is injective
    iff ``kernel_rank = 0`` and surjective iff the cokernel vanishes.
    """
    m = len(matrix)
    n = len(matrix[0]) if m else 0
    facts = _invariant_factors(matrix)
    r = len(facts)
    return {
        "rank": r,
        "kernel_rank": n - r,
        "cokernel_free_rank": m - r,
        "cokernel_torsion": [d for d in facts if d > 1],
        "injective": n - r == 0,
        "surjective": (m - r == 0) and all(d == 1 for d in facts),
    }


def homology_ranks(dims, maps):
    """
    Homology of a complex of free abelian groups ``⋯ → ℤ^{dims[k+1]} --maps[k+1]--> ℤ^{dims[k]} → ⋯``.

    ``maps[k]`` is the matrix of ``d_k`` (size ``dims[k-1] × dims[k]``). Returns
    ``{k: {"betti", "torsion"}}`` — a self-contained re-derivation of the homology used to test
    exactness.
    """
    ranks = {k: _rank(maps[k]) for k in maps}
    out = {}
    for k in dims:
        betti = dims[k] - ranks.get(k, 0) - ranks.get(k + 1, 0)
        torsion = [d for d in _invariant_factors(maps.get(k + 1, [])) if d > 1]
        out[k] = {"betti": betti, "torsion": torsion}
    return out


def is_exact(dims, maps):
    """
    Whether a complex of free abelian groups is **exact** — equivalently all its homology
    vanishes (zero Betti numbers and no torsion at every position).
    """
    H = homology_ranks(dims, maps)
    return all(h["betti"] == 0 and not h["torsion"] for h in H.values())


if __name__ == "__main__":
    print("=" * 70)
    print("Homological algebra — exact sequences, Ext, and Tor")
    print("=" * 70)

    print("\nHom(Z/m, Z/n) ≅ Ext¹(Z/m, Z/n) ≅ Tor₁(Z/m, Z/n) ≅ Z/gcd(m,n):")
    print(f"  {'(m,n)':<10}{'gcd':<6}{'Hom':<8}{'Ext¹':<8}{'Tor₁':<8}all agree")
    for m, n in [(2, 2), (2, 3), (4, 6), (6, 9), (12, 8)]:
        h, e, t = hom_cyclic(m, n), ext1_cyclic(m, n), tor1_cyclic(m, n)
        agree = h["order"] == e["order"] == t["order"] == gcd(m, n)
        print(f"  {str((m, n)):<10}{gcd(m, n):<6}{h['group']:<8}{e['group']:<8}{t['group']:<8}{agree}")

    print("\nThe exact sequence 0 → ker φ → Z^n → Z^m → coker φ → 0 of a matrix φ:")
    phi = [[2, 4], [2, 6]]            # the Klein-four presentation again
    ks = kernel_cokernel_sequence(phi)
    print(f"  φ = {phi}:  rank {ks['rank']}, ker rank {ks['kernel_rank']}, "
          f"coker = Z^{ks['cokernel_free_rank']} ⊕ {'⊕'.join(f'Z/{d}' for d in ks['cokernel_torsion']) or '0'}")
    print(f"  injective: {ks['injective']}   surjective: {ks['surjective']}")

    print("\nExactness of a complex of free groups (homology must vanish):")
    # 0 → Z --×2--> Z is NOT exact at the right Z (coker = Z/2); but with Z/2 added it would be.
    seq1 = ({0: 1, 1: 1}, {1: [[2]]})            # Z --2--> Z  (a 2-term complex)
    seq2 = ({0: 1, 1: 1, 2: 1}, {1: [[0]], 2: [[1]]})  # Z --1--> Z --0--> Z  (split, exact in middle)
    print(f"  Z --×2--> Z : homology {homology_ranks(*seq1)}  exact: {is_exact(*seq1)}")
    print(f"  Z --1--> Z --0--> ... at middle: H_1 = {homology_ranks(*seq2)[1]}")
