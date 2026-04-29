"""
csf_trees.sage
==============

SageMath implementation accompanying the Math 725 Independent Project
report on Martin--Morin--Wagner, "On Distinguishing Trees by Their
Chromatic Symmetric Functions".
"""

from itertools import combinations
from collections import defaultdict


# ---------------------------------------------------------------------------
# Symmetric-function ring (power-sum basis)
# ---------------------------------------------------------------------------

Sym = SymmetricFunctions(QQ)
p   = Sym.power()


# ---------------------------------------------------------------------------
# Type partition of an edge set, computed by union--find
# ---------------------------------------------------------------------------

def _type_partition(vertices, edge_set):
    """
    The type partition of the spanning subgraph (`vertices`, `edge_set`):
    a Sage Partition whose parts are the orders of the connected
    components, sorted weakly decreasingly.

    Implemented via a small dictionary-based union--find rather than by
    constructing a Sage Graph each iteration, since this routine is
    called O(2^|E|) times when computing X_G.
    """
    parent = {v: v for v in vertices}

    def find(x):
        # path compression
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for u, v in edge_set:
        union(u, v)

    sizes = defaultdict(int)
    for v in vertices:
        sizes[find(v)] += 1

    return Partition(sorted(sizes.values(), reverse=True))


# ---------------------------------------------------------------------------
# Chromatic symmetric function (Stanley's power-sum expansion)
# ---------------------------------------------------------------------------

def chromatic_symmetric_function(G):
    """
    Return the chromatic symmetric function X_G of a finite simple
    Sage graph G, as an element of Sym.power() (the power-sum basis).

    Computed via Stanley's expansion (eq. (1) of the report).
    """
    vertices = list(G.vertices(sort=False))
    edges    = [tuple(e) for e in G.edges(labels=False, sort=False)]

    # Accumulate coefficients in a plain dict for speed, then convert
    # to a single symmetric-function element at the end.
    coeffs = defaultdict(int)
    for size in range(len(edges) + 1):
        sign = (-1)**size
        for A in combinations(edges, size):
            coeffs[_type_partition(vertices, A)] += sign

    return p.sum_of_terms(
        ((lam, c) for lam, c in coeffs.items() if c != 0),
        distinct=True,
    )


def csf_coefficient_dict(G):
    """
    The CSF of G as a dict {Partition: c_lambda(G)}, with zero
    coefficients dropped.  Convenience wrapper around
    `chromatic_symmetric_function`.
    """
    return dict(chromatic_symmetric_function(G).monomial_coefficients())


def csf_signature(G):
    """A hashable canonical CSF representation, suitable for equality
    testing across non-isomorphic graphs."""
    return tuple(sorted(
        (tuple(lam), c) for lam, c in csf_coefficient_dict(G).items()
    ))


# ---------------------------------------------------------------------------
# Caterpillars
# ---------------------------------------------------------------------------

def caterpillar(leaf_numbers):
    """
    Build the caterpillar with given leaf-number sequence
    (f_0, f_1, ..., f_s) as a Sage Graph.  Spine has s edges and s+1
    vertices; total vertex count is (s + 1) + sum(f_i).
    """
    if len(leaf_numbers) < 2:
        raise ValueError("a caterpillar's spine must have at least one edge")

    s = len(leaf_numbers) - 1
    G = Graph()
    spine = [('v', i) for i in range(s + 1)]
    for v in spine:
        G.add_vertex(v)
    for i in range(s):
        G.add_edge(spine[i], spine[i + 1])
    for i, fi in enumerate(leaf_numbers):
        for j in range(int(fi)):
            G.add_edge(spine[i], ('leaf', i, j))
    return G


# ---------------------------------------------------------------------------
# Reconstruction (proof of Theorem 5.1 of the report, in code)
# ---------------------------------------------------------------------------

def is_singleton_free(lam):
    """Every part of the partition lam is at least 2."""
    return all(part >= 2 for part in lam)


def reconstruct_caterpillar(coeffs):
    """
    Given the coefficient dictionary `coeffs` of the chromatic symmetric
    function of a caterpillar T whose leaf numbers are strictly positive
    and pairwise distinct, return the leaf-number sequence
    (f_0, ..., f_s) of T (canonicalised so that the smaller endpoint
    comes first).

    Implements the four-step proof of Theorem 5.1 of the report.

    Raises ValueError if the data are inconsistent with the hypotheses.
    """
    sf_nonzero = [(lam, c) for lam, c in coeffs.items()
                  if is_singleton_free(lam) and c != 0]
    if not sf_nonzero:
        raise ValueError("no singleton-free non-zero coefficient found")

    # Step 1: lambda = type(L) is the unique singleton-free partition
    # of maximum length with non-zero coefficient.
    max_len = max(len(lam) for lam, _ in sf_nonzero)
    lam_candidates = [lam for lam, _ in sf_nonzero if len(lam) == max_len]
    if len(lam_candidates) != 1:
        raise ValueError(
            "lambda is not unique; distinctness hypothesis is violated"
        )
    lam = lam_candidates[0]
    s   = len(lam) - 1

    # Step 2: the multiset of leaf numbers is { part - 1 : part in lambda }.
    leaf_numbers_multiset = [int(q) - 1 for q in lam]
    if any(f <= 0 for f in leaf_numbers_multiset):
        raise ValueError("a leaf number is non-positive")
    if len(set(leaf_numbers_multiset)) != len(leaf_numbers_multiset):
        raise ValueError("leaf numbers are not pairwise distinct")

    if s == 0:
        return tuple(leaf_numbers_multiset)

    # Step 3: mu_1, ..., mu_s = type(L u {e_i}) are exactly the
    # singleton-free partitions of length s with non-zero coefficient.
    mus = [lam2 for lam2, _ in sf_nonzero if len(lam2) == s]
    if len(mus) != s:
        raise ValueError(
            f"expected exactly s = {s} singleton-free length-{s} "
            f"partitions, found {len(mus)}"
        )

    # Step 4: each mu_i is obtained from lambda by replacing some pair
    # of parts {a, b} by the single part (a + b).  Build the auxiliary
    # graph H on the leaf-number set, then walk the resulting path.
    lam_count = defaultdict(int)
    for q in lam:
        lam_count[int(q)] += 1

    auxiliary_edges = []
    for mu in mus:
        mu_count = defaultdict(int)
        for q in mu:
            mu_count[int(q)] += 1

        # parts of lambda not in mu (with multiplicity)
        missing = []
        for q, ct in lam_count.items():
            extra = ct - mu_count.get(q, 0)
            missing.extend([q] * max(0, extra))

        # parts of mu not in lambda (with multiplicity)
        new_parts = []
        for q, ct in mu_count.items():
            extra = ct - lam_count.get(q, 0)
            new_parts.extend([q] * max(0, extra))

        if (len(missing) != 2 or len(new_parts) != 1
                or missing[0] + missing[1] != new_parts[0]):
            raise ValueError(
                f"mu = {mu} is not a single 2-part merge of lambda = {lam}"
            )

        a, b = missing
        auxiliary_edges.append((a - 1, b - 1))

    # Build H as a Sage graph.
    H = Graph()
    for f in leaf_numbers_multiset:
        H.add_vertex(f)
    for u, v in auxiliary_edges:
        H.add_edge(u, v)

    deg1 = [v for v in H.vertices(sort=False) if H.degree(v) == 1]
    if len(deg1) != 2 or any(H.degree(v) > 2
                             for v in H.vertices(sort=False)):
        raise ValueError("auxiliary graph H is not a path")

    # Canonical orientation: smaller endpoint first.
    start = min(deg1)
    sequence = [start]
    visited = {start}
    while True:
        nbrs = [w for w in H.neighbors(sequence[-1]) if w not in visited]
        if not nbrs:
            break
        sequence.append(nbrs[0])
        visited.add(nbrs[0])

    return tuple(sequence)


# ---------------------------------------------------------------------------
# Empirical search for collisions (Conjecture 7.1)
# ---------------------------------------------------------------------------

def _compositions(n, k):
    """All k-tuples of non-negative integers summing to n."""
    if k == 1:
        yield (n,)
        return
    for first in range(n + 1):
        for rest in _compositions(n - first, k - 1):
            yield (first,) + rest


def all_caterpillar_sequences(n):
    """
    Yield each non-isomorphic caterpillar on n vertices exactly once,
    represented by its canonicalised leaf-number sequence
    (f_0, ..., f_s) with f_0 >= 1, f_s >= 1, and (f_0, ..., f_s) <=
    (f_s, ..., f_0) lexicographically.

    Total vertex count: (s + 1) + sum(f_i) = n.
    """
    yielded = set()
    for s in range(1, n):
        leaves_total = n - s - 1
        if leaves_total < 2:  # need f_0 >= 1 and f_s >= 1
            continue
        for seq in _compositions(leaves_total, s + 1):
            if seq[0] >= 1 and seq[-1] >= 1:
                key = min(seq, seq[::-1])
                if key not in yielded:
                    yielded.add(key)
                    yield key


def test_caterpillar_conjecture(n_max):
    """
    For each n from 4 to n_max, enumerate all non-isomorphic caterpillars
    on n vertices, compute their CSF signatures, and report whether any
    two distinct caterpillars share a CSF.
    """
    print(f"{'n':>3} | {'caterpillars':>12} | {'distinct CSFs':>14} "
          f"| result")
    print("-" * 60)
    for n in range(4, n_max + 1):
        sig_to_seq = {}
        collisions = []
        total = 0
        for seq in all_caterpillar_sequences(n):
            total += 1
            T = caterpillar(seq)
            sig = csf_signature(T)
            if sig in sig_to_seq:
                collisions.append((sig_to_seq[sig], seq))
            else:
                sig_to_seq[sig] = seq
        result = ("OK (no collisions)" if not collisions
                  else f"COLLISIONS: {collisions}")
        print(f"{n:>3} | {total:>12} | {len(sig_to_seq):>14} | {result}")


# ---------------------------------------------------------------------------
# Demonstrations matching Section 6 of the report
# ---------------------------------------------------------------------------

def demo_section_6_1():
    """Reproduce Example 6.1: full CSF expansion of T_{(2,1,1)}."""
    print("=" * 65)
    print("Example 6.1: full CSF expansion of the caterpillar with")
    print("leaf-number sequence (1, 2)  (== the spider T_{(2,1,1)}).")
    print("=" * 65)

    T  = caterpillar((1, 2))
    XT = chromatic_symmetric_function(T)

    # SageMath prints power-sum elements very nicely.
    print("\nX_T (in the power-sum basis):")
    print(f"   {XT}")

    coeffs = csf_coefficient_dict(T)
    n_v = T.num_verts()
    m   = T.num_edges()

    print("\nProposition 4.1 verifications:")
    print(f"  (i)   |V(T)| = {n_v}")
    c_2_111 = coeffs.get(Partition([2, 1, 1, 1]), 0)
    print(f"  (ii)  -c_(2,1,1,1) = {-c_2_111}  "
          f"(should equal m = {m})")
    c_n1_1 = coeffs.get(Partition([n_v - 1, 1]), 0)
    print(f"  (v)   |c_(n-1,1)| = {abs(c_n1_1)}  (number of leaves)")

    expected = {
        Partition([1, 1, 1, 1, 1]):  1,
        Partition([2, 1, 1, 1]):    -4,
        Partition([3, 1, 1]):        4,
        Partition([2, 2, 1]):        2,
        Partition([4, 1]):          -3,
        Partition([3, 2]):          -1,
        Partition([5]):              1,
    }
    assert coeffs == expected, (
        f"\nMISMATCH:\n  computed: {coeffs}\n  expected: {expected}"
    )
    print("\nAll coefficients match equation (3) of the report.")


def demo_section_6_2():
    """Reproduce Example 6.2: reconstruct caterpillar with f = (1, 2, 3)."""
    print()
    print("=" * 65)
    print("Example 6.2: reconstruct the caterpillar with leaf-number")
    print("sequence (1, 2, 3) from its chromatic symmetric function.")
    print("=" * 65)

    seq = (1, 2, 3)
    T   = caterpillar(seq)
    coeffs = csf_coefficient_dict(T)

    print(f"\nT has n = {T.num_verts()} vertices and "
          f"m = {T.num_edges()} edges.")

    print("Singleton-free coefficients of X_T:")
    sf = [(lam, c) for lam, c in coeffs.items() if is_singleton_free(lam)]
    sf.sort(key=lambda x: (-len(x[0]), tuple(-int(q) for q in x[0])))
    for lam, c in sf:
        parts = ",".join(str(int(q)) for q in lam)
        print(f"   c_({parts}) = {c:+d}")

    recovered = reconstruct_caterpillar(coeffs)
    print(f"\nRecovered leaf-number sequence: {recovered}")
    print(f"Original leaf-number sequence:  {seq}")

    if recovered in (seq, seq[::-1]):
        print("RECONSTRUCTION SUCCESSFUL "
              "(matches original up to reversal).")
    else:
        raise AssertionError("reconstruction failed")


def demo_conjecture_search(n_max=11):
    """Run the empirical test of Conjecture 7.1 up to n = n_max."""
    print()
    print("=" * 65)
    print(f"Empirical test of Conjecture 7.1: every caterpillar on n <= "
          f"{n_max}")
    print("vertices is determined by its chromatic symmetric function.")
    print("=" * 65)
    print()
    test_caterpillar_conjecture(n_max=n_max)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    demo_section_6_1()
    demo_section_6_2()
    demo_conjecture_search(n_max=11)
