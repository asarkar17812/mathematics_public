"""Generate a dataset of the primes up to ``N``.

A small data-generation script: it walks the primes ``2 ≤ p ≤ N`` (via Sage's
``prime_range``) and records, for each, its value, its index in the sequence (so ``p`` is
the ``index``-th prime, i.e. ``π(p) = index``), and ``log₁₀ p``. Saved together these let
you plot the prime-counting function ``π(x)`` and eyeball the Prime Number Theorem
(``π(x) ≈ x / ln x``) directly from the data.

Run as a script to (re)build the dataset:

    sage src/maths_coding/number_theory/utils/primes_generator.py

Output (written under ``data_stores/primes/`` as paired ``.csv`` and ``.npy``):

==========  ==========================================================
column      meaning
==========  ==========================================================
``p``       the prime
``index``   1-based index in the prime sequence, equal to ``π(p)``
``log10_p`` ``log₁₀ p``
``time_sec``per-prime wall-clock marker (negligible here; kept for a
            uniform schema with the heavier generators)
==========  ==========================================================
"""

from sage.all import *
import numpy as np
import csv
import time
import math
import os

def generate_dataset(N):
    """Build the prime records for all primes ``p`` with ``2 ≤ p ≤ N`` (see module docstring)."""
    records = []

    for idx, p in enumerate(prime_range(2, N + 1), start=1):
        start = time.perf_counter()
        elapsed = time.perf_counter() - start
        log10_p = math.log10(p)

        print(f"p={p:6d} | idx={idx:6d} | t={elapsed:8.6f}s")

        records.append((
            int(p),
            int(idx),
            float(log10_p),
            float(elapsed)
        ))

    return records

if __name__ == "__main__":
    N = 10000000

    outdir = "data_stores/primes"
    os.makedirs(outdir, exist_ok=True)

    records = generate_dataset(N)

    csv_path = os.path.join(outdir, f"primes_dataset_upto_{N}.csv")
    npy_path = os.path.join(outdir, f"primes_dataset_upto_{N}.npy")

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["p", "index", "log10_p", "time_sec"])
        writer.writerows(records)

    dtype = np.dtype([
        ("p", np.int64),
        ("index", np.int64),
        ("log10_p", np.float64),
        ("time_sec", np.float64),
    ])

    arr = np.array(records, dtype=dtype)
    np.save(npy_path, arr)

    print(f"\nSaved to: {outdir}")