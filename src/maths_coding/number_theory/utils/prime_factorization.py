"""Factorize every integer up to ``N`` with a smallest-prime-factor sieve.

Trial-dividing each ``n`` separately would be wasteful; instead this module sieves once to
record, for every ``n ≤ N``, its **smallest prime factor** ``spf[n]``. With that table the
full factorization of any ``n`` is read off in ``O(number of prime factors)`` time by
repeatedly dividing out ``spf[n]`` — a classic space-for-time trade that factorizes a whole
range cheaply.

The dataset it builds is aimed at the *statistics* of factorization, recording two
prime-counting functions for each ``n``:

- ``ω(n)`` — the number of **distinct** prime factors;
- ``Ω(n)`` — the number of prime factors **with multiplicity**.

These are the subject of the Hardy–Ramanujan and Erdős–Kac theorems (``ω(n)`` is normally
distributed around ``ln ln n``), so the dataset is a hook for exploring that empirically.

Run as a script to (re)build the dataset:

    sage src/maths_coding/number_theory/utils/prime_factorization.py

Output (written under ``data_stores/factorizations/`` as paired ``.csv`` and ``.npy``):

================  ====================================================
column            meaning
================  ====================================================
``n``             the integer (``2 ≤ n ≤ N``)
``log10_n``       ``log₁₀ n``
``omega``         ``ω(n)``, distinct prime factors
``Omega``         ``Ω(n)``, prime factors with multiplicity
``factorization`` list of ``(prime, exponent)`` pairs, ascending
``time_sec``      per-n wall-clock time for the factorization lookup
================  ====================================================
"""

from sage.all import *
import numpy as np
import csv
import time
import math
import os

def build_spf(N):
    """
    Build the smallest-prime-factor table for ``0 … N``.

    Returns a list ``spf`` where ``spf[n]`` is the smallest prime dividing ``n``
    (with ``spf[n] == n`` left in place for ``n`` prime, ``0``, and ``1``). This is a
    modified Sieve of Eratosthenes: when ``i`` is found prime, it is stamped as the
    smallest factor of its still-unmarked multiples.
    """
    spf = list(range(N + 1))

    for i in range(2, int(N**0.5) + 1):
        if spf[i] == i:  # i is prime
            for j in range(i * i, N + 1, i):
                if spf[j] == j:
                    spf[j] = i

    return spf

def factorize_spf(n, spf):
    """
    Factorize ``n`` using a precomputed smallest-prime-factor table ``spf``.

    Repeatedly divide out ``spf[n]`` until ``n`` reaches 1, accumulating exponents.

    Returns
    -------
    list[tuple[int, int]]
        ``(prime, exponent)`` pairs in ascending order of prime.
    """
    factors = {}

    while n > 1:
        p = spf[n]
        factors[p] = factors.get(p, 0) + 1
        n //= p

    return [(int(p), int(e)) for p, e in sorted(factors.items())]

def generate_dataset(N):
    """Factorize every ``n`` in ``2 … N`` and assemble the records (see module docstring)."""
    records = []
    spf = build_spf(N)

    for n in range(2, N + 1):
        start = time.perf_counter()

        factor_list = factorize_spf(n, spf)

        elapsed = time.perf_counter() - start

        omega = len(factor_list)
        Omega = sum(e for _, e in factor_list)
        log10_n = math.log10(n)

        print(f"n={n:6d} | ω={omega:2d} | Ω={Omega:2d} | t={elapsed:8.6f}s")

        records.append((
            int(n),
            float(log10_n),
            int(omega),
            int(Omega),
            factor_list,  
            float(elapsed)
        ))

    return records

if __name__ == "__main__":
    N = 10000000

    outdir = "data_stores/factorizations"
    os.makedirs(outdir, exist_ok=True)

    records = generate_dataset(N)

    csv_path = os.path.join(outdir, f"factorizations_{N}.csv")
    npy_path = os.path.join(outdir, f"factorizations_{N}.npy")

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["n", "log10_n", "omega", "Omega", "factorization", "time_sec"])

        for row in records:
            n, log10_n, omega, Omega, factor_list, t = row
            writer.writerow([n, log10_n, omega, Omega, str(factor_list), t])

    dtype = np.dtype([
        ("n", np.int64),
        ("log10_n", np.float64),
        ("omega", np.int64),
        ("Omega", np.int64),
        ("factorization", object),  # stores list of tuples
        ("time_sec", np.float64),
    ])

    arr = np.array(records, dtype=dtype)
    np.save(npy_path, arr)

    print(f"\nSaved to: {outdir}")