# The unique nearest graph to the balanced complete tripartite graph

For every $n\geq3$, the graph $(n-1)K_1\vee(K_1\sqcup K_{n,n})$ is the unique graph up to isomorphism nearest to $K_{n,n,n}$ among the other simple graphs of order $3n$, for both Manhattan and Euclidean adjacency-spectral distance. This proves Conjecture 20. The proof also answers Problem 21, since squared Euclidean distance at most $2/3$ forces adjacency rank at most five, including for disconnected graphs. [Original conjecture and problem](https://arxiv.org/html/2610.18636v2#S5)

[Proof note](proof.pdf) · [Editable proof](PROOF.md) · [X draft](POST-DRAFT.md)

## Replay the finite certificates

Use Python 3.10 with NetworkX 3.3, NumPy 1.25.1, and SymPy 1.13.1.

```sh
python verify_rank_reduction.py --problem21
python verify_conjecture20.py
```

The rank verifier checks that the graph atlas covers every labeled input at the required orders. The spectral check uses rational eigenvalue intervals. The second script checks every permitted rank-four template, the small parameter cases, and the polynomial inequalities in the proof.
