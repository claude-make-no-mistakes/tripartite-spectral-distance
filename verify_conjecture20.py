"""Check the small catalogues and polynomial inequalities in the Conjecture 20 proof."""

from itertools import combinations, product
import json

import networkx as nx
from sympy import Matrix, Poly, Rational, Symbol, cancel, expand

from verify_rank_reduction import L1_THRESHOLD, SQUARE_THRESHOLD, root_intervals, verify_catalogue


PAW = "CN"
P4 = "Ck"
BULL = "DNG"
HOUSE = "DNg"


def graph_matrix(graph):
    return Matrix(nx.to_numpy_array(graph, dtype=int).tolist())


def reduced_extensions(seed):
    """Enumerate all twin-free rank-r extensions of an invertible order-r seed, omitting isolates."""
    size = seed.rows
    inverse = seed.inv()
    fixed = [seed[:, j] for j in range(size)]
    extra = []
    for signature in product(range(2), repeat=size):
        if any(signature):
            vector = Matrix(signature)
            if vector not in fixed and (vector.T * inverse * vector)[0] == 0:
                extra.append(vector)
    for count in range(len(extra) + 1):
        for columns in combinations(extra, count):
            signatures = Matrix.hstack(*fixed, *columns)
            result = signatures.T * inverse * signatures
            if all(x in (0, 1) for x in result):
                yield result


def induced_bounds(matrix):
    roots = root_intervals(matrix)
    positive = [max(0, lo) for lo, _ in roots[:-1]]
    negative = [max(0, -hi) for _, hi in roots[2:]]
    return 2 * max(sum(positive), sum(negative)), sum(x*x for x in positive + negative)


def verify_templates():
    atlas = [g for g in nx.graph_atlas_g() if len(g) == 4]
    verify_catalogue(atlas, 4)
    kinds = {g6: nx.from_graph6_bytes(g6.encode("ascii")) for g6 in (PAW, P4, BULL, HOUSE)}
    count = 0
    admitted = set()
    for graph in atlas:
        seed = graph_matrix(graph)
        if seed.det() == 0:
            continue
        for matrix in reduced_extensions(seed):
            count += 1
            assert matrix.rows <= 6  # Ensures the interlacing bounds apply even at n=3.
            l1, square = induced_bounds(matrix)
            if l1 > L1_THRESHOLD and square > SQUARE_THRESHOLD:
                continue
            candidate = nx.Graph()
            candidate.add_nodes_from(range(matrix.rows))
            candidate.add_edges_from((i, j) for i, j in combinations(range(matrix.rows), 2) if matrix[i, j])
            matches = [key for key, base in kinds.items() if nx.is_isomorphic(base, candidate)]
            assert len(matches) == 1, list(candidate.edges())
            admitted.add(matches[0])
    assert admitted == set(kinds)
    triangle = graph_matrix(nx.complete_graph(3))
    assert list(reduced_extensions(triangle)) == [triangle]
    return {"four_vertex_seeds": len(atlas), "extensions_checked": count, "admitted": sorted(admitted)}


def paw_polynomial(x, a, c, d):
    return x**4 - (a*(1+c+d) + c*d)*x*x - 2*a*c*d*x + a*c*d


def distance_intervals(n, sizes, edges):
    order = sum(sizes)
    assert order == 3 * n
    real_adjacency = Matrix(len(sizes), len(sizes),
                            lambda i, j: sizes[j] if (min(i, j), max(i, j)) in edges else 0)
    # This symmetricizable quotient carries all possible nonzero roots.
    roots = sorted(root_intervals(real_adjacency) + [(0, 0)] * (order - len(sizes)))
    targets = [-n, -n] + [0] * (order - 3) + [2 * n]
    lo1 = hi1 = lo2 = hi2 = 0
    for (lo, hi), target in zip(roots, targets, strict=True):
        lo -= target
        hi -= target
        lower = 0 if lo <= 0 <= hi else min(abs(lo), abs(hi))
        upper = max(abs(lo), abs(hi))
        lo1 += lower
        hi1 += upper
        lo2 += lower * lower
        hi2 += upper * upper
    return (lo1, hi1), (lo2, hi2)


def verify_small_n(n):
    paw = {(0, 1), (0, 2), (0, 3), (2, 3)}
    triangle = {(0, 1), (0, 2), (1, 2)}
    incumbent = distance_intervals(n, (n - 1, 1, n, n), paw)
    cases = 0
    for a in range(1, 3 * n - 2):
        for c in range(1, 3 * n - a):
            d = 3 * n - a - c
            if a <= c <= d and (a, c, d) != (n, n, n):
                alternative = distance_intervals(n, (a, c, d), triangle)
                assert all(lo > hi for (lo, _), (_, hi) in zip(alternative, incumbent))
                cases += 1
            d -= 1
            if c <= d and (a, c, d) != (n - 1, n, n):
                alternative = distance_intervals(n, (a, 1, c, d), paw)
                assert all(lo > hi for (lo, _), (_, hi) in zip(alternative, incumbent))
                cases += 1
    return {"n": n, "exact_competitors": cases}


def positive_polynomial(expression, variable):
    coefficients = Poly(expand(expression), variable).all_coeffs()
    assert coefficients[-1] > 0 and all(c >= 0 for c in coefficients), expression


def verify_algebra():
    m = Symbol("m")
    n = m + 3
    r = 2*n - Rational(2, 3)
    area = 3*n - 2
    b = 2*r/(r - 1)
    numerator, denominator = cancel(3*r*r - 2*(b+area)*r + b*area).as_numer_denom()
    positive_polynomial(numerator, m)
    positive_polynomial(denominator, m)
    r = 2*n - Rational(12, 25)
    area = 3*n - 1
    positive_polynomial(3*r*r - 2*(1+area)*r + area, m)
    assert Rational(7, 12) - Rational(12, 25)**2 < Rational(11, 18)**2
    assert Rational(12, 25)**2 + Rational(1, 3)**2 + Rational(1, 2)**2 > Rational(7, 12)

    n = m + 6
    h = n + Rational(11, 18)
    for a, c in ((n + 1, n), (n - 1, n + 2)):
        positive_polynomial(2*a*c*h + a*c - h*h*(a+c), m)
    positive_polynomial(n - Rational(11, 18)*h, m)
    for a, c, d in ((n-2, n, n+1), (n-1, n-1, n+1), (n, n-1, n)):
        positive_polynomial(paw_polynomial(Rational(1, 3), a, c, d), m)
        positive_polynomial(-paw_polynomial(-n + Rational(1, 2), a, c, d), m)
    # Quotient expansion for a triangle with class sizes a,c,d and a leaf at the a-class.
    a, c, d, x = (Symbol(name) for name in ("a", "c", "d", "x"))
    quotient = Matrix([[0, 1, c, d], [a, 0, 0, 0], [a, 0, 0, d], [a, 0, c, 0]])
    assert expand(quotient.charpoly(x).as_expr() - paw_polynomial(x, a, c, d)) == 0


if __name__ == "__main__":
    print(json.dumps(verify_templates()))
    for n in (3, 4, 5):
        print(json.dumps(verify_small_n(n)))
    verify_algebra()
    print(json.dumps({"symbolic_inequalities": "passed"}))
