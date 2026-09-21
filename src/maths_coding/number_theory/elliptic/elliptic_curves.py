"""Elliptic curve E: y² = x³ + 2x + 3 seen globally (Q), geometrically (R), and locally (F_p).

**This is a demonstration script** — running it prints a narrated walkthrough and pops up a
matplotlib plot of the real curve with rational points overlaid. Run it with Sage:

    sage src/maths_coding/number_theory/elliptic/elliptic_curves.py

The story it tells, in order:

1. **Global structure over Q.** Compute the rank, torsion subgroup, and generators of
   ``E(Q) ≅ Z^r ⊕ torsion`` (Mordell–Weil).
2. **Local structure mod p.** Reduce ``E`` modulo several primes, count ``#E(F_p)``, and
   read off the group structure; primes dividing the discriminant give *bad reduction*.
3. **Local constrains global.** Track how each rational torsion point reduces into
   ``E(F_p)`` — the injection of prime-to-``p`` torsion is *why* a few small primes already
   bound the global torsion.
4. **Picture.** Plot the real locus implicitly (as the zero contour of
   ``y² - (x³ + ax + b)``, so the whole curve is drawn) with a finite sample of rational
   points — torsion plus small multiples of the generators.

See :mod:`maths_coding.number_theory.dynamics.frobenius_eigenvalues` for the Frobenius /
eigenvalue side of the same curve.
"""

from sage.all import *
import numpy as np
import matplotlib.pyplot as plt

"""
============================================================
ELLIPTIC CURVE: GLOBAL vs LOCAL STRUCTURE (WITH VISUALIZATION)
============================================================

We study an elliptic curve:
    y^2 = x^3 + ax + b

This object lives in multiple "worlds":

1. Over Q (global arithmetic object)
2. Over R (continuous geometric shape)
3. Over F_p (finite, local reductions)

Key philosophy:
    local data (mod p) constrains global structure.

We will:
- compute E(Q)
- reduce mod p
- track torsion behavior
- visualize the real curve + rational points
"""

# --------------------------------------------------
# Define elliptic curve
# --------------------------------------------------
a = 2
b = 3
E = EllipticCurve([a, b])

print("Elliptic curve:", E)

"""
The curve defines a group law geometrically:
- points form an abelian group
- identity = point at infinity
- addition via chord–tangent rule

Algebraically:
    E(Q) ≅ Z^r ⊕ torsion
"""

# --------------------------------------------------
# Global structure over Q
# --------------------------------------------------
rank = E.rank()
torsion = E.torsion_subgroup()
gens = E.gens()

print("\n--- Global structure over Q ---")
print("Rank:", rank)
print("Torsion subgroup:", torsion)
print("Generators:", gens)

"""
Interpretation:

- Rank = number of independent infinite directions
- Torsion = finite "looping" structure
- Generators = basis for infinite part

Think:
    operators (multiplication by n) act on these points
    generating orbits inside the curve
"""

# --------------------------------------------------
# Reduction mod p (local structure)
# --------------------------------------------------
primes_to_check = [5, 7, 11, 13]

print("\n--- Local reductions mod p ---")

for p in primes_to_check:
    if E.discriminant() % p == 0:
        print(f"p={p}: bad reduction (singularity appears)")
        continue

    Ep = E.change_ring(GF(p))

    print(f"\np = {p}")
    print("  #E(F_p) =", Ep.cardinality())
    print("  Group structure:", Ep.abelian_group())

"""
Concept:

Reduction mod p gives:
    E(Q) → E(F_p)

This is a group homomorphism.

Important fact:
    torsion (away from p) injects into E(F_p)

So finite field data constrains global torsion.
"""

# --------------------------------------------------
# Torsion reduction behavior
# --------------------------------------------------
print("\n--- Torsion reduction behavior ---")

for p in primes_to_check:
    if E.discriminant() % p == 0:
        continue

    Ep = E.change_ring(GF(p))

    print(f"\np = {p}")
    for T in E.torsion_points():
        if not T.is_zero():
            red = Ep(T[0] % p, T[1] % p)
            print(f"  {T} → mod {p} → {red}")

"""
This shows how torsion points "descend" into finite fields.

Think of this as:
    a global structure being probed by local projections.
"""

# --------------------------------------------------
# Generate rational points (finite sample)
# --------------------------------------------------
points = []

points.extend(E.torsion_points())

for P in gens:
    for k in range(-5, 6):
        if k != 0:
            points.append(k * P)

points = list(set(points))

# --------------------------------------------------
# Plot using implicit curve (FULL EXTENSION)
# --------------------------------------------------

"""
Instead of solving y = ±sqrt(...), we plot implicitly:

    y^2 - (x^3 + ax + b) = 0

This ensures:
- full curve is drawn
- no artificial truncation
- captures entire geometry inside window
"""

x = np.linspace(-5, 5, 1000)
y = np.linspace(-5, 5, 1000)
X, Y = np.meshgrid(x, y)

Z = Y**2 - (X**3 + a*X + b)

plt.figure(figsize=(6, 6))

# contour at Z = 0 draws the curve
plt.contour(X, Y, Z, levels=[0])

# --------------------------------------------------
# Plot rational points
# --------------------------------------------------
for P in points:
    if not P.is_zero():
        plt.scatter(float(P[0]), float(P[1]), s=40)

# --------------------------------------------------
# Formatting
# --------------------------------------------------
plt.axhline(0)
plt.axvline(0)
plt.title("Elliptic Curve: Real Geometry + Rational Points")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)

plt.show()