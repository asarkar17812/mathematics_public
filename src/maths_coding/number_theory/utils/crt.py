"""Chinese Remainder Theorem — recombining an integer from its residues.

The CRT says that knowing ``x`` modulo several moduli is the same as knowing ``x`` modulo
their lcm — and when the moduli are pairwise coprime, every combination of residues
occurs for exactly one ``x`` mod the product. This is the precise sense in which an
integer can be reconstructed from purely *local* data (its residues mod each prime power),
the elementary template for the "local determines global" idea that runs through the rest
of the package.

The module offers three views of the same theorem:

- :func:`crt_pair` — the two-congruence base case, built directly from the extended
  Euclidean algorithm (``xgcd``); it also handles non-coprime moduli, raising when the two
  congruences are inconsistent.
- :func:`crt` — folds :func:`crt_pair` across a whole system, so it works for arbitrary
  (not necessarily coprime) moduli.
- :func:`crt_coprime` — the classical product formula ``x = Σ aᵢ Mᵢ (Mᵢ⁻¹ mod mᵢ)`` with
  ``Mᵢ = M / mᵢ``, faster when the moduli are known to be pairwise coprime.

:func:`verify_crt` independently checks a claimed solution against every congruence.

Examples
--------
>>> crt([2, 3, 2], [3, 5, 7])      # doctest: +SKIP
(23, 105)
>>> crt_coprime([2, 3, 2], [3, 5, 7])  # doctest: +SKIP
(23, 105)
"""

from sage.all import ZZ, xgcd
from functools import reduce

def crt_pair(a1, m1, a2, m2):
    """
    Solve:
        x ≡ a1 (mod m1)
        x ≡ a2 (mod m2)

    Returns (x, lcm(m1,m2)) if solvable, else raises ValueError.
    """
    a1, m1, a2, m2 = ZZ(a1), ZZ(m1), ZZ(a2), ZZ(m2)

    g, s, t = xgcd(m1, m2)

    # consistency condition
    if (a2 - a1) % g != 0:
        raise ValueError("No solution: incompatible congruences")

    lcm = (m1 * m2) // g

    # combine
    x = (a1 + (a2 - a1) // g * s * m1) % lcm

    return x, lcm

def crt(residues, moduli):
    """
    Solve system:
        x ≡ residues[i] (mod moduli[i])

    Returns (x, M) where M = lcm of moduli.
    """
    if len(residues) != len(moduli):
        raise ValueError("residues and moduli must have same length")

    x, m = ZZ(residues[0]), ZZ(moduli[0])

    for a_i, m_i in zip(residues[1:], moduli[1:]):
        x, m = crt_pair(x, m, a_i, m_i)

    return x, m

def crt_coprime(residues, moduli):
    """
    Faster CRT assuming all moduli are pairwise coprime.
    """
    M = reduce(lambda a, b: a*b, moduli, ZZ(1))
    x = ZZ(0)

    for a_i, m_i in zip(residues, moduli):
        M_i = M // m_i
        _, inv, _ = xgcd(M_i, m_i)  # inverse of M_i mod m_i
        x += a_i * inv * M_i

    return x % M, M

def verify_crt(x, residues, moduli):
    return all(ZZ(x) % m == ZZ(a) % m for a, m in zip(residues, moduli))