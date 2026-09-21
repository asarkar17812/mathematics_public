"""Primality and factorization.

Two complementary questions about integers, with a Mersenne-number slant:

- :mod:`lucas_lehmer_test` — *is this number prime?* The Lucas–Lehmer test decides
  primality of Mersenne numbers ``M_p = 2ᵖ - 1`` with a short modular recurrence, and the
  module also exhibits the exact closed form of that recurrence inside the number field
  ``Q(√3)``.
- :mod:`fermat_kraitchik` — *if not, what are its factors?* Fermat's difference-of-squares
  method (with Kraitchik's framing) factors an odd ``N`` by searching for
  ``N = a² - b² = (a-b)(a+b)``.

Mersenne numbers tie the two together and lead directly to even perfect numbers (see
:mod:`maths_coding.number_theory.utils.perfect_numbers`).
"""
