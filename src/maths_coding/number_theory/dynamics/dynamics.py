"""One idea, many dynamical systems: operator + eigenstructure + arithmetic constraints.

**This is a demonstration script** — run it with Sage to see the printout and a series of
matplotlib plots:

    sage src/maths_coding/number_theory/dynamics/dynamics.py

It iterates several maps over finite sets and shows that the same picture governs them all.
Because the underlying set is finite, every trajectory eventually cycles, so the natural
data is the distribution of orbit / cycle lengths.

Sections:

1. **Elliptic curves + Frobenius.** For ``E: y² = x³ + 2x + 3`` and small primes, the
   trace of Frobenius ``a_p`` and the Frobenius polynomial ``x² - a_p x + p``.
2. **Orbits in a finite field.** The squaring map ``x ↦ x²`` on ``GF(11)`` decomposed into
   cycles, with a bar chart of cycle lengths.
3. **Linear dynamics.** A fixed matrix acting on ``GF(7)²``, scatter-plotting orbit lengths
   over all start vectors.
4. **Frobenius ↔ eigenvalues.** The characteristic polynomial of Frobenius and the Hasse
   bound ``|λ| = √p`` — the arithmetic mirror of classical eigenvalue dynamics.
5. **Nonlinear group dynamics.** The doubling map ``P ↦ 2P`` on ``E(F_p)``, histogramming
   orbit lengths.

The closing summary states the unifying slogan the whole package circles around:
``dynamics = operator + eigenstructure + arithmetic constraints``, with elliptic curves and
Frobenius as its deepest instance.
"""

from sage.all import *
import matplotlib.pyplot as plt

print("=== EXTENDED STRUCTURE: GROUPS + DYNAMICS + ELLIPTIC CURVES ===\n")

# ==================================================
# SECTION 1: ELLIPTIC CURVES + FROBENIUS
# ==================================================

print("---- ELLIPTIC CURVES + FROBENIUS ----\n")

E = EllipticCurve([2, 3])
print("Elliptic Curve:", E)

primes = [5, 7, 11, 13]

frobenius_data = []

for p in primes:
    if E.discriminant() % p == 0:
        continue
    
    Ep = E.change_ring(GF(p))
    Np = Ep.cardinality()
    ap = E.ap(p)  # Frobenius trace
    
    print(f"p={p}: #E(F_p)={Np}, a_p={ap}")
    
    frob_poly = E.frobenius_polynomial(p)
    print("Frobenius polynomial:", frob_poly)
    
    frobenius_data.append((p, ap))

print("\nInterpretation:")
print("a_p = trace of Frobenius → eigenvalues encode arithmetic structure\n")


# ==================================================
# SECTION 2: ORBIT VISUALIZATION (FINITE FIELD)
# ==================================================

print("---- ORBIT STRUCTURE IN FINITE FIELD ----\n")

p = 11
F = GF(p)

def f(x):
    return x^2

visited = set()
cycles = []

for a in F:
    if a in visited or a == 0:
        continue
    
    orbit = []
    x = a
    
    while x not in orbit:
        orbit.append(x)
        x = f(x)
    
    cycles.append(orbit)
    visited.update(orbit)

print("Cycles in GF(11) under x -> x^2:")
for c in cycles:
    print(c)


# Visualization (cycle lengths)
lengths = [len(c) for c in cycles]

plt.figure()
plt.bar(range(len(lengths)), lengths)
plt.title("Cycle lengths in GF(11) under x^2")
plt.xlabel("Cycle index")
plt.ylabel("Length")
plt.show()


print("\n---- LINEAR DYNAMICS (MATRIX ACTION) ----\n")

F2 = GF(7)
A = Matrix(F2, [[2,1],[0,1]])

points = []
labels = []

for x in F2:
    for y in F2:
        v = vector(F2, [x,y])
        
        orbit = []
        w = v
        
        for _ in range(10):
            if w in orbit:
                break
            orbit.append(w)
            w = A*w
        
        points.append(len(orbit))
        labels.append(v)

# Plot orbit lengths
plt.figure()
plt.scatter(range(len(points)), points)
plt.title("Orbit lengths under matrix action over GF(7)")
plt.xlabel("Vector index")
plt.ylabel("Orbit length")
plt.show()

print("\n---- FROBENIUS ↔ EIGENVALUE LINK ----\n")

for p, ap in frobenius_data:
    print(f"p={p}")
    print(f"Trace a_p = {ap}")
    print(f"Characteristic polynomial: x^2 - {ap}x + {p}")
    print()

print("Eigenvalues λ satisfy |λ| = sqrt(p) (Hasse bound)")
print("This mirrors classical eigenvalue dynamics.\n")

print("---- ELLIPTIC CURVE POINT ORBITS ----\n")

p = 11
Ep = E.change_ring(GF(p))

points = Ep.points()

# Use doubling map
def double(P):
    return 2*P

orbit_lengths = []

for P in points:
    if P.is_zero():
        continue
    
    orbit = []
    Q = P
    
    for _ in range(20):
        if Q in orbit:
            break
        orbit.append(Q)
        Q = double(Q)
    
    orbit_lengths.append(len(orbit))

plt.figure()
plt.hist(orbit_lengths, bins=range(1, max(orbit_lengths)+2))
plt.title("Orbit lengths under doubling map on E(F_p)")
plt.xlabel("Orbit length")
plt.ylabel("Frequency")
plt.show()


print("\n---- SUMMARY ----\n")

print("""
You are now seeing the same structure in multiple forms:

1. Classical groups → operators preserving structure
2. Eigenvalues → control dynamics
3. Cyclotomics → classify periodicity
4. Finite fields → exact cycles
5. Frobenius → arithmetic operator on elliptic curves
6. Elliptic curve maps → nonlinear group dynamics

Unifying principle:

    dynamics = operator + eigenstructure + arithmetic constraints

Elliptic curves + Frobenius are the deepest instance of this.
""")

print("\n=== END ===")