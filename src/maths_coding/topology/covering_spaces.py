"""Covering spaces of the circle: the Galois correspondence, made concrete.

A covering ``p : X̃ → X`` lays disjoint "sheets" over each small open set; for connected ``X``
the number of sheets — the **degree** — is the constant size of a fiber ``p⁻¹(x)``. The two
prototypes both live over the circle ``S¹``:

    universal cover   p : ℝ → S¹,   t ↦ e^{2πit}      (infinite degree, fiber ≅ ℤ)
    finite covers     p : S¹ → S¹,  z ↦ zⁿ            (degree n, fiber ≅ ℤ/n)

The engine is **lifting**: a path downstairs lifts uniquely once a starting sheet is chosen,
and homotopic loops lift to paths with the same endpoint. So the fundamental group
``π₁(S¹) = ℤ`` acts on the fiber by **monodromy** — and for the ``n``-fold cover that action
is the *cyclic shift* of the ``n`` sheets, the generator of the deck-transformation group.

The deepest statement is the **Galois correspondence of covering spaces** — formally identical
to Galois theory, with ``π₁`` in the role of the Galois group:

    connected cover X̃          ⟷   subgroup  p_*π₁(X̃) ≤ π₁(X)
    degree of the cover         ⟷   index of the subgroup
    regular (normal) cover      ⟷   normal subgroup
    deck transformation group   ⟷   quotient  π₁(X) / p_*π₁(X̃)
    universal cover             ⟷   trivial subgroup {e}

For ``X = S¹`` this is the lattice of subgroups ``nℤ ≤ ℤ``: the ``n``-fold cover ``z ↦ zⁿ``
corresponds to ``nℤ`` (index ``n``), is regular because ``ℤ`` is abelian, and has deck group
``ℤ/nℤ``. The universal cover ``ℝ → S¹`` corresponds to ``{0}`` and has deck group ``ℤ``, which
realizes ``π₁(S¹)`` literally as a symmetry group (``ℝ/ℤ ≅ S¹``).

- :func:`fiber` / :func:`degree` — the sheets ``ℤ/n`` of the ``n``-fold cover, and its degree.
- :func:`monodromy` — the permutation of sheets from going once around the base loop.
- :func:`lift_loop_endpoint` — where a loop wound ``w`` times lands (monodromy as a ℤ-action).
- :func:`deck_group` — the deck-transformation group (``ℤ/n``, or ``ℤ`` for the universal cover).
- :func:`galois_correspondence` — the subgroup ``nℤ`` ↔ ``n``-fold cover dictionary, tabulated.
"""

from sage.all import ZZ, CyclicPermutationGroup, Permutation


# The n-fold cover z ↦ z^n of S^1: sheets are labelled 0, 1, …, n-1 (the n preimages of a
# base point, i.e. the n branches of z^{1/n}). The universal cover is encoded by n = 0.


def fiber(n):
    """The fiber ``p⁻¹(x)`` of the ``n``-fold cover: the sheet labels ``{0, …, n−1}`` (``ℤ`` if ``n = 0``)."""
    if n == 0:
        return "Z (the integer fiber of the universal cover R → S^1)"
    return list(range(n))


def degree(n):
    """The degree (number of sheets); ``∞`` for the universal cover ``n = 0``."""
    return "infinity" if n == 0 else int(n)


def monodromy(n):
    """
    The monodromy permutation of the fiber induced by traversing the base loop once.

    Lifting the loop ``t ↦ e^{2πit}`` moves sheet ``k`` to sheet ``k+1 (mod n)``, so the
    monodromy of the generator of ``π₁(S¹) = ℤ`` is the cyclic ``n``-shift. Returned as a Sage
    :class:`Permutation` of ``{1, …, n}`` (1-indexed). For ``n = 0`` the monodromy is the shift
    ``k ↦ k+1`` on ``ℤ`` (described as a string).
    """
    if n == 0:
        return "k ↦ k+1 on Z (free ℤ-action; the universal cover has no finite cycle)"
    if n == 1:
        return Permutation([1])
    return Permutation([(i % n) + 1 for i in range(1, n + 1)])  # (1 2 ... n)


def lift_loop_endpoint(n, start_sheet, winding):
    """
    Lift the loop wound ``winding`` times around ``S¹``, starting on ``start_sheet``; return the
    ending sheet.

    This is monodromy as the ``ℤ``-action on the fiber: ``k ↦ k + winding (mod n)``. That the
    endpoint depends only on the *winding number* (not the particular loop) is homotopy lifting;
    that it can differ from the start is why the cover is nontrivial.
    """
    if n == 0:
        return start_sheet + winding  # on the integer fiber, no wraparound
    return (start_sheet + winding) % n


def deck_group(n):
    """
    The deck-transformation group ``π₁(S¹) / p_*π₁(X̃)`` of the ``n``-fold cover.

    Returns a dict with ``order``, a ``label``, and (for finite ``n ≥ 1``) ``group`` as a Sage
    :class:`CyclicPermutationGroup`. Every cover of ``S¹`` is regular because ``ℤ`` is abelian,
    so the deck group always has order equal to the degree: ``ℤ/nℤ`` (or ``ℤ`` for ``n = 0``).
    """
    if n == 0:
        return {"order": "infinity", "label": "Z", "group": None}
    return {"order": int(n), "label": f"Z/{n}Z", "group": CyclicPermutationGroup(n)}


def galois_correspondence(max_degree=6):
    """
    Tabulate the covering-space ↔ subgroup dictionary for ``S¹``.

    Returns a list of rows, one per connected cover, each a dict with ``cover``, ``subgroup``
    (``nℤ``), ``index`` (= degree), ``regular`` (always ``True`` for ``S¹``), and ``deck_group``.
    Includes the universal cover ``ℝ → S¹`` ↔ ``{0}``.
    """
    rows = [{
        "cover": "R → S^1 (universal)",
        "subgroup": "{0}",
        "index": "infinity",
        "regular": True,
        "deck_group": "Z",
    }]
    for n in range(1, max_degree + 1):
        rows.append({
            "cover": f"z ↦ z^{n}",
            "subgroup": f"{n}Z",
            "index": int(n),
            "regular": True,          # Z is abelian ⇒ every subgroup normal ⇒ every cover regular
            "deck_group": deck_group(n)["label"],
        })
    return rows


if __name__ == "__main__":
    print("=" * 70)
    print("Covering spaces of S¹ — π₁ = ℤ acting on fibers by monodromy")
    print("=" * 70)

    print("\nThe n-fold cover z ↦ z^n:")
    for n in [1, 2, 3, 4]:
        mono = monodromy(n)
        print(f"  n = {n}:  degree {degree(n)},  fiber {fiber(n)},  "
              f"monodromy {mono.cycle_string()},  deck group {deck_group(n)['label']}")

    print("\nMonodromy as the ℤ-action on the fiber (lifting a loop wound w times):")
    n = 4
    print(f"  cover z ↦ z^{n}, start on sheet 0:")
    for w in range(0, 6):
        print(f"    wound {w} time(s)  →  ends on sheet {lift_loop_endpoint(n, 0, w)}"
              + ("   (back to start — loop lifts to a loop)" if lift_loop_endpoint(n, 0, w) == 0 and w > 0 else ""))

    print("\nThe Galois correspondence  {covers of S¹}  ⟷  {subgroups of ℤ}:")
    print(f"  {'cover':<22}{'subgroup':<10}{'index':<10}{'regular':<9}{'deck group'}")
    for r in galois_correspondence(6):
        print(f"  {r['cover']:<22}{r['subgroup']:<10}{str(r['index']):<10}{str(r['regular']):<9}{r['deck_group']}")
    print("\n  (degree = index, deck group = ℤ/nℤ = quotient; universal cover ℝ ↔ {0}, deck group ℤ)")
