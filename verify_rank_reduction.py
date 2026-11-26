"""Check rank obstructions for Conjecture 20 and Problem 21 by exact arithmetic."""

import argparse
from itertools import combinations, permutations
import json

import networkx as nx
from sympy import Matrix, Poly, Rational, Symbol, zeros


EPSILON = Rational(1, 10**8)
L1_THRESHOLD = Rational(11, 9)
SQUARE_THRESHOLD = Rational(7, 12)


def root_intervals(adjacency):
    variable = Symbol("x")
    polynomial = Poly(adjacency.charpoly(variable).as_expr(), variable)
    roots = [interval for interval, multiplicity in polynomial.intervals(eps=EPSILON)
             for _ in range(multiplicity)]
    assert len(roots) == adjacency.rows, "A symmetric matrix must have a full real spectrum"
    return roots


def edge_mask(graph, positions, permutation):
    mask = 0
    for u, v in graph.edges():
        edge = tuple(sorted((permutation[u], permutation[v])))
        mask |= 1 << positions[edge]
    return mask


def verify_catalogue(graphs, order):
    """Confirm coverage of all labeled graphs via permutation orbits."""
    positions = {edge: i for i, edge in enumerate(combinations(range(order), 2))}
    covered = set()
    for graph in graphs:
        orbit = {edge_mask(graph, positions, relabel) for relabel in permutations(range(order))}
        assert not covered & orbit, "Graph atlas contains two isomorphic inputs"
        covered.update(orbit)
    assert len(covered) == 1 << len(positions), "Graph atlas misses a labeled graph"


def verify_obstructions(order, square_threshold):
    graphs = [g for g in nx.graph_atlas_g() if len(g) == order]
    verify_catalogue(graphs, order)
    nonsingular = 0
    smallest = [None, None]
    for graph in graphs:
        adjacency = Matrix(nx.to_numpy_array(graph, dtype=int).tolist())
        if adjacency.det() == 0:
            continue
        nonsingular += 1
        roots = root_intervals(adjacency)
        positive = [max(0, lo) for lo, _ in roots[:-1]]
        negative = [max(0, -hi) for _, hi in roots[2:]]
        l1 = 2 * max(sum(positive), sum(negative))
        square = sum(x**2 for x in positive + negative)
        assert l1 > L1_THRESHOLD and square > square_threshold, nx.to_graph6_bytes(graph, header=False)
        for i, bound in enumerate((l1, square)):
            if smallest[i] is None or bound < smallest[i]:
                smallest[i] = bound
    return {"order": order, "isomorphism_classes": len(graphs), "invertible": nonsingular,
            "least_certified_bounds": [str(bound) for bound in smallest]}


def compare_reference(n):
    full = zeros(3 * n)
    sizes = [n - 1, 1, n, n]
    groups = [kind for kind, size in enumerate(sizes) for _ in range(size)]
    links = {(0, 1), (0, 2), (0, 3), (2, 3)}
    for i, j in combinations(range(3 * n), 2):
        if (groups[i], groups[j]) in links:
            full[i, j] = full[j, i] = 1
    spectrum = [-n, -n] + [0] * (3 * n - 3) + [2 * n]
    differences = [max(abs(lo - ref), abs(hi - ref))
                   for (lo, hi), ref in zip(root_intervals(full), spectrum, strict=True)]
    assert sum(x**2 for x in differences) < SQUARE_THRESHOLD
    return {"n": n, "reference_squared_distance_upper": str(sum(x**2 for x in differences))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--problem21", action="store_true", help="Also verify the 2/3 cutoff through order 7")
    args = parser.parse_args()
    for order in (5, 6):
        threshold = Rational(2, 3) if order == 6 else SQUARE_THRESHOLD
        print(json.dumps(verify_obstructions(order, threshold)), flush=True)
    if args.problem21:
        print(json.dumps(verify_obstructions(7, Rational(2, 3))), flush=True)
    for n in (3, 4):
        print(json.dumps(compare_reference(n)), flush=True)
    # Lemma 5 of the source gives Lambda < 17/30 and 2/5 < lambda_2 < 1/2 for n >= 5.
    # The squared Euclidean distance in Theorem 6 is convex in both arguments.
    for x in (0, Rational(17, 30)):
        for y in (Rational(2, 5), Rational(1, 2)):
            assert 2 * (x*x - x*y + y*y) < SQUARE_THRESHOLD


if __name__ == "__main__":
    main()
