"""Frobenius eigenvalues of an elliptic curve: the Hasse bound, Sato–Tate angles, and point counts.

**This is a demonstration script** — run it with Sage to see the printout and a plot:

    sage src/maths_coding/number_theory/dynamics/frobenius_eigenvalues.py

For ``E: y² = x³ + 2x + 3``, reduction mod ``p`` has a Frobenius endomorphism whose
characteristic polynomial is ``x² - a_p x + p``, where ``a_p = p + 1 - #E(F_p)`` is the
*trace of Frobenius*. Its two roots (eigenvalues) ``λ, λ̄`` are complex conjugates with

    |λ| = √p        (the **Hasse bound**, equivalently ``|a_p| ≤ 2√p``),

so they can be written ``λ = √p · e^{iθ}`` with ``cos θ = a_p / (2√p)``. The script:

1. **Frobenius data + angles.** For several primes, print ``a_p``, check
   ``#E(F_p) = p + 1 - a_p``, and extract the Sato–Tate angle ``θ``.
2. **Point-count recurrence.** The power sums ``S_n = λⁿ + λ̄ⁿ`` satisfy the linear
   recurrence ``S_{n+1} = a_p S_n - p S_{n-1}`` (``S₀ = 2``, ``S₁ = a_p``), and
   ``#E(F_{pⁿ}) = pⁿ + 1 - S_n``. The script verifies this against Sage's actual counts
   over ``F_{pⁿ}``.
3. **Plot.** The normalized eigenvalues ``e^{iθ}`` on the unit circle — the empirical data
   whose limiting distribution is the subject of the Sato–Tate conjecture (now a theorem).
"""

from sage.all import *
import matplotlib.pyplot as plt
import math

E = EllipticCurve([2, 3])
print("Elliptic Curve:", E)

primes = [5, 7, 11, 13, 17, 19, 23]

angles = []

print("\n--- Frobenius data ---\n")

for p in primes:
    if E.discriminant() % p == 0:
        continue

    ap = E.ap(p)
    Np = E.change_ring(GF(p)).cardinality()

    print(f"p = {p}")
    print(f"  a_p = {ap}")
    print(f"  #E(F_p) = {Np} (check: {p + 1 - ap})")

    cos_theta = ap / (2 * math.sqrt(p))

    # numerical safety
    cos_theta = max(min(cos_theta, 1), -1)

    theta = math.acos(cos_theta)
    angles.append(theta)

    print(f"  cos(theta) = {cos_theta:.4f}")
    print(f"  theta = {theta:.4f} radians\n")

print("\n--- Recurrence check ---\n")

p = 11
ap = E.ap(p)

print(f"Using p = {p}, a_p = {ap}")

# recurrence:
# S_0 = 2
# S_1 = a_p
# S_{n+1} = a_p S_n - p S_{n-1}

S = [2, ap]

for n in range(2, 6):
    Sn = ap*S[-1] - p*S[-2]
    S.append(Sn)

print("S_n values:", S)

print("\nCompare with actual point counts:")

for n in range(1, 5):
    Epn = E.change_ring(GF(p^n))
    count = Epn.cardinality()
    predicted = p**n + 1 - S[n]
    print(f"n={n}: actual={count}, predicted={predicted}")

print("\n--- Plotting normalized eigenvalues ---\n")

x_vals = [math.cos(t) for t in angles]
y_vals = [math.sin(t) for t in angles]

# unit circle
theta = [i * 2 * math.pi / 200 for i in range(200)]
circle_x = [math.cos(t) for t in theta]
circle_y = [math.sin(t) for t in theta]

plt.figure()
plt.plot(circle_x, circle_y)
plt.scatter(x_vals, y_vals)

plt.title("Normalized Frobenius eigenvalues (angles)")
plt.xlabel("Re")
plt.ylabel("Im")
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True)

plt.show()