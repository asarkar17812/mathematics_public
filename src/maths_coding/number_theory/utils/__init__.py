"""Number-theory toolbox and dataset generators.

A grab-bag of standalone routines that the rest of the package leans on, plus the scripts
that build the datasets in ``data_stores/``.

Library routines
----------------
crt
    Chinese Remainder Theorem: glue congruences ``x ≡ aᵢ (mod mᵢ)`` into one.
quadratic_reciprocity
    The Legendre symbol three ways (Euler's criterion, Gauss's lemma), the law of
    quadratic reciprocity, Legendre's formula ``vₚ(n!)``, and the Eisenstein
    irreducibility test.
divisibility_tests
    A uniform divisibility test for any prime built from the modular inverse of 10.
perfect_numbers
    Even perfect numbers via the Euclid–Euler correspondence with Mersenne primes.
hensel_wilson
    Hensel lifting (Newton's method ``p``-adically) of a simple root up the tower
    ``p, p², p³, …``; Wilson's theorems on products of units; and their "Henson" synthesis
    with the CRT, e.g. building ``√2`` in ``ℤ_p`` and gluing across primes.

Dataset generators (run as scripts)
----------------------------------
primes_generator
    Enumerate primes up to ``N`` with their index and ``log₁₀``.
prime_factorization
    A smallest-prime-factor sieve, recording ``ω(n)`` and ``Ω(n)`` for every ``n ≤ N``.
mersenne_generator
    Run the Lucas–Lehmer test across prime exponents and record the results.

Each generator writes paired ``.csv`` / ``.npy`` files; those outputs are git-ignored
build artifacts, not source.
"""
