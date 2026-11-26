# The Unique Nearest Graph to the Balanced Complete Tripartite Graph

For every $n \geq 3$, the graph $H_n = (n-1)K_1 \vee (K_1 \sqcup K_{n,n})$ is the unique simple graph of order $3n$ other than $K_{n,n,n}$, up to isomorphism, minimizing the Manhattan ($\ell^1$) adjacency-spectral distance to $K_{n,n,n}$, and also the unique minimizer of the Euclidean ($\ell^2$) adjacency-spectral distance. This proves Conjecture 20 of [Jovanović, Sedlar, and Škrekovski (arXiv:2610.18636v2)](https://arxiv.org/html/2610.18636v2#S5). The same induced-subgraph interlacing method resolves Problem 21 of that paper: for every $n \geq 3$, any simple graph of order $3n$ (connected or disconnected) with squared Euclidean adjacency-spectral distance at most $2/3$ from $K_{n,n,n}$ has adjacency rank at most five.

## Files

- [proof.pdf](proof.pdf) — Compiled publication-ready paper ([proof.tex](proof.tex))
- [PROOF.md](PROOF.md) — Markdown version of the proof
- [verify_rank_reduction.py](verify_rank_reduction.py) — Exact rational-interval verification of rank obstructions (orders 5–7) and incumbent bounds ($n = 3, 4$, and corners for $n \geq 5$)
- [verify_conjecture20.py](verify_conjecture20.py) — Exact verification of rank-four templates, small-$n$ competitors ($n = 3, 4, 5$), and symbolic polynomial inequalities
- [rank-certificate.jsonl](rank-certificate.jsonl) and [conjecture20-certificate.jsonl](conjecture20-certificate.jsonl) — Recorded outputs of the verification scripts

## Replaying the finite certificates

Requires Python 3.10+ with `networkx`, `numpy`, and `sympy` (tested on Python 3.10 with NetworkX 3.3 / NumPy 1.25.1 / SymPy 1.13.1, and on Python 3.12 with NetworkX 3.7 / NumPy 2.5.5 / SymPy 1.14.0):

```sh
python3 verify_rank_reduction.py --problem21
python3 verify_conjecture20.py
```

All finite spectral checks use exact integer characteristic polynomials and rational root intervals (no floating-point comparisons). The rank verifier certifies graph-atlas completeness at each order by partitioning all labeled graphs into disjoint permutation orbits. The second script checks every compatible rank-four signature template, compares all remaining competitors for $n \in \{3, 4, 5\}$ against $H_n$, and verifies all symbolic polynomial inequalities used in the proof.
