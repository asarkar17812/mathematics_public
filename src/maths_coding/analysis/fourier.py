"""Fourier analysis on ℤ/N — characters, the convolution theorem, and Pontryagin duality.

A function on the cyclic group ``ℤ/N`` is decomposed into **characters** — the homomorphisms
``χ_k : ℤ/N → ℂˣ``, ``χ_k(n) = e^{2πi kn/N}`` — which are the irreducible representations of the
(abelian) group and form an **orthogonal basis**:

    ⟨χ_j, χ_k⟩ = Σ_{n} χ_j(n) \\overline{χ_k(n)} = N·δ_{jk}.

The **discrete Fourier transform** ``X_k = Σ_n x_n e^{-2πi kn/N}`` is the coordinate map into
that basis, and three facts make it the central tool of harmonic analysis:

    Convolution theorem:   ``\\widehat{x ⊛ y} = \\hat x · \\hat y``   (cyclic convolution ↦ pointwise product),
    Parseval/Plancherel:   ``Σ |x_n|² = (1/N) Σ |X_k|²``            (the DFT is, up to scale, an isometry),
    Pontryagin duality:    ``\\widehat{ℤ/N} ≅ ℤ/N``                  (the characters are themselves a ℤ/N).

The last is the finite, self-dual case of Pontryagin duality; the character ``χ_k`` is indexed by
``k ∈ ℤ/N``, so the dual group is again ``ℤ/N``. This is the abelian shadow of
:mod:`maths_coding.algebra.representation_theory`: the DFT is the decomposition into the
one-dimensional irreducible characters, and character orthogonality is Schur orthogonality.

Pure-Python (``cmath``), so it runs without SageMath:

- :func:`character` — the character ``χ_k`` of ``ℤ/N``.
- :func:`dft`, :func:`idft` — the transform and its inverse (``idft(dft(x)) = x``).
- :func:`cyclic_convolution` — ``(x ⊛ y)_m = Σ_n x_n y_{m−n}``.
- :func:`convolution_theorem_check` — ``\\widehat{x ⊛ y} = \\hat x · \\hat y``.
- :func:`parseval_check` — the Plancherel identity.
- :func:`character_orthogonality` — ``⟨χ_j, χ_k⟩ = N δ_{jk}`` (Pontryagin self-duality of ℤ/N).
"""

import cmath


def character(N, k):
    """The character ``χ_k`` of ``ℤ/N``: returns the function ``n ↦ e^{2πi kn/N}``."""
    return lambda n: cmath.exp(2j * cmath.pi * k * n / N)


def dft(x):
    """The discrete Fourier transform ``X_k = Σ_n x_n e^{-2πi kn/N}``."""
    N = len(x)
    return [sum(x[n] * cmath.exp(-2j * cmath.pi * k * n / N) for n in range(N))
            for k in range(N)]


def idft(X):
    """The inverse DFT ``x_n = (1/N) Σ_k X_k e^{+2πi kn/N}`` (so ``idft(dft(x)) = x``)."""
    N = len(X)
    return [sum(X[k] * cmath.exp(2j * cmath.pi * k * n / N) for k in range(N)) / N
            for n in range(N)]


def cyclic_convolution(x, y):
    """The cyclic convolution ``(x ⊛ y)_m = Σ_n x_n · y_{(m−n) mod N}``."""
    N = len(x)
    return [sum(x[n] * y[(m - n) % N] for n in range(N)) for m in range(N)]


def convolution_theorem_check(x, y, tol=1e-9):
    """
    Verify the **convolution theorem**: ``DFT(x ⊛ y) = DFT(x) · DFT(y)`` (pointwise).

    Returns a dict with ``max_error`` and ``holds``.
    """
    lhs = dft(cyclic_convolution(x, y))
    fx, fy = dft(x), dft(y)
    rhs = [a * b for a, b in zip(fx, fy)]
    err = max(abs(a - b) for a, b in zip(lhs, rhs))
    return {"max_error": err, "holds": err < tol}


def parseval_check(x, tol=1e-9):
    """
    Verify the **Parseval/Plancherel** identity ``Σ |x_n|² = (1/N) Σ |X_k|²``.

    Returns a dict with ``time_energy``, ``freq_energy``, and ``holds``.
    """
    N = len(x)
    X = dft(x)
    time_energy = sum(abs(v) ** 2 for v in x)
    freq_energy = sum(abs(v) ** 2 for v in X) / N
    return {
        "time_energy": time_energy,
        "freq_energy": freq_energy,
        "holds": abs(time_energy - freq_energy) < tol,
    }


def character_orthogonality(N, tol=1e-9):
    """
    Verify character orthogonality ``⟨χ_j, χ_k⟩ = Σ_n χ_j(n) \\overline{χ_k(n)} = N·δ_{jk}``.

    This is Pontryagin duality made concrete: the ``N`` characters are an orthogonal basis (the
    dual group ``\\widehat{ℤ/N} ≅ ℤ/N``). Returns a dict with ``holds`` (all inner products match).
    """
    ok = True
    for j in range(N):
        cj = character(N, j)
        for k in range(N):
            ck = character(N, k)
            ip = sum(cj(n) * ck(n).conjugate() for n in range(N))
            expected = N if j == k else 0
            if abs(ip - expected) > tol:
                ok = False
    return {"holds": ok}


if __name__ == "__main__":
    print("=" * 70)
    print("Fourier analysis on ℤ/N — characters, convolution, duality")
    print("=" * 70)

    x = [1.0, 2.0, 0.0, -1.0]          # a length-4 signal
    X = dft(x)
    print("\nSignal x =", x)
    print("DFT(x)   =", [complex(round(v.real, 4), round(v.imag, 4)) for v in X])
    rec = [round(v.real, 9) for v in idft(X)]
    print("idft(dft(x)) =", rec, " recovers x:", rec == x)

    print("\nConvolution theorem (cyclic convolution ↦ pointwise product):")
    y = [0.0, 1.0, 1.0, 0.0]
    ct = convolution_theorem_check(x, y)
    print(f"  x ⊛ y = {[round(v.real, 4) for v in cyclic_convolution(x, y)]}")
    print(f"  DFT(x⊛y) = DFT(x)·DFT(y)?  max error {ct['max_error']:.2e}  → {ct['holds']}")

    print("\nParseval / Plancherel (the DFT is an isometry up to 1/N):")
    pv = parseval_check(x)
    print(f"  Σ|xₙ|² = {pv['time_energy']:.4f}   (1/N)Σ|X_k|² = {pv['freq_energy']:.4f}  → {pv['holds']}")

    print("\nPontryagin duality of ℤ/N — the characters are an orthogonal basis (⟨χⱼ,χ_k⟩ = N δⱼ_k):")
    for N in (4, 5, 6, 7):
        print(f"  N = {N}: character orthogonality holds: {character_orthogonality(N)['holds']}"
              f"   (so Ẑ/{N} ≅ ℤ/{N})")
