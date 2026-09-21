"""Arithmetic dynamics — iterating maps and the Frobenius operator.

What happens when you apply a map over and over on a *finite* set? Because the set is
finite, every trajectory must eventually cycle, so the natural invariants are orbit and
cycle lengths. This subpackage explores that question for several maps and shows how the
same picture — *operator + eigenstructure + arithmetic constraints* — governs all of them.

Modules
-------
dynamics
    A narrated walkthrough: squaring and matrix maps over finite fields, the doubling map
    on ``E(F_p)``, and the Frobenius trace ``a_p`` as the eigenvalue datum tying these
    together.
frobenius_eigenvalues
    The Frobenius eigenvalues of an elliptic curve: the Hasse bound ``|λ| = √p``, the
    Sato–Tate angle ``θ``, and the linear recurrence ``S_{n+1} = a_p S_n - p S_{n-1}``
    that predicts point counts ``#E(F_{pⁿ})``.

Both modules are demonstration scripts — run them to see the printout and plots.
"""
