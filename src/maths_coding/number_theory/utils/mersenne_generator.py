"""Run the Lucas–Lehmer test across prime exponents and record the results.

For each prime ``p ≤ N`` this evaluates the Lucas–Lehmer test on the Mersenne number
``M_p = 2ᵖ - 1`` and stores both the verdict and the *size* of the final residue
``s_{p-2} mod M_p``. The flagged ``is_mp == 1`` rows are exactly the Mersenne-prime
exponents (2, 3, 5, 7, 13, 17, 19, 31, …), which by Euclid–Euler also enumerate the even
perfect numbers (see :mod:`maths_coding.number_theory.utils.perfect_numbers`).

A numerical detail worth noting: the final residue ``s_mod`` can have thousands of digits,
so its base-10 logarithm is computed from the top 53 bits of the Sage integer
(``bits = s_mod.nbits()``; ``log₁₀ s_mod ≈ (bits - 53)·log₁₀2 + log₁₀(mantissa)``) rather
than by converting the whole integer to a float, which would overflow. For Mersenne primes
the residue is 0 and the log is recorded as ``-inf``.

Run as a script to (re)build the dataset:

    sage src/maths_coding/number_theory/utils/mersenne_generator.py

Output (written under ``data_stores/mersenne/`` as paired ``.csv`` and ``.npy``):

==============  ======================================================
column          meaning
==============  ======================================================
``p``           prime exponent
``log10_Mp``    ``log₁₀ M_p = p·log₁₀2``
``log10_s_mod`` ``log₁₀`` of the final Lucas–Lehmer residue (``-inf`` if 0)
``is_mp``       1 iff ``M_p`` is prime
``time_sec``    wall-clock time for the test at this ``p``
==============  ======================================================

Note
----
Run from the ``src/`` directory (or with ``src`` on ``PYTHONPATH``) so that the
``from number_theory.primality...`` import resolves.
"""

from sage.all import *
import numpy as np
import csv
import time
import math
import os
from src.maths_coding.number_theory.primality.lucas_lehmer_test import lucas_lehmer_test_with_value

def generate_dataset(N):
    """Run the Lucas–Lehmer test for every prime ``p ≤ N`` and assemble the records."""
    records = []

    for p in prime_range(2, N + 1):
        start = time.perf_counter()
        is_mp, s_mod = lucas_lehmer_test_with_value(p)
        elapsed = time.perf_counter() - start

        log10_Mp = p * math.log10(2)
        if s_mod == 0:
            log10_s_mod = float("-inf")
        else:
            bits = s_mod.nbits()
            mantissa = int(s_mod >> (bits - 53))
            log10_s_mod = (bits - 53) * math.log10(2) + math.log10(mantissa)

        print(f"p={p:5d} | t={elapsed:8.4f}s | prime={is_mp}")

        records.append((
            int(p),
            float(log10_Mp),
            float(log10_s_mod),
            int(is_mp),
            float(elapsed)
        ))

    return records

if __name__ == "__main__":
    N = 10000
    outdir = "data_stores/mersenne"
    os.makedirs(outdir, exist_ok=True)
    records = generate_dataset(N)
    csv_path = os.path.join(outdir, f"mersenne_dataset_upto_{N}.csv")
    npy_path = os.path.join(outdir, f"mersenne_dataset_upto_{N}.npy")

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["p", "log10_Mp", "log10_s_mod", "is_mp", "time_sec"])
        writer.writerows(records)

    dtype = np.dtype([
        ("p", np.int64),
        ("log10_Mp", np.float64),
        ("log10_s_mod", np.float64),
        ("is_mp", np.int8),
        ("time_sec", np.float64),
    ])

    arr = np.array(records, dtype=dtype)
    np.save(npy_path, arr)

    print(f"\nSaved to: {outdir}")