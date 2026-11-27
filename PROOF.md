# The Unique Nearest Graph to the Balanced Complete Tripartite Graph

## Abstract

Let $D_1$ and $D_2$ denote the $\ell^1$ (Manhattan) distance and the **squared** $\ell^2$ (Euclidean) distance, respectively, between sorted adjacency spectra of graphs of the same order. We prove that for every $n \geq 3$, the graph
$$
H_n = (n-1)K_1 \vee (K_1 \sqcup K_{n,n})
$$
is the unique simple graph of order $3n$ other than $K_{n,n,n}$, up to isomorphism, minimizing $D_1(K_{n,n,n}, G)$, and likewise the unique minimizer of $D_2(K_{n,n,n}, G)$. This proves Conjecture 20 of [Jovanović, Sedlar, and Škrekovski, arXiv:2610.18636v2](https://arxiv.org/html/2610.18636v2#S5) (whose Theorem 17 already establishes the asserted asymptotic limits). By the same induced-subgraph interlacing method, we also resolve Problem 21 of that paper, showing that for every $n \geq 3$, any simple graph $G$ of order $3n$ (connected or disconnected) satisfying $D_2(K_{n,n,n}, G) \leq 2/3$ has adjacency rank at most five.

## 1. Introduction

For two simple graphs $G$ and $H$ of order $N$ with sorted adjacency spectra $\alpha_1 \geq \cdots \geq \alpha_N$ and $\gamma_1 \geq \cdots \geq \gamma_N$, define
$$
D_1(G, H) = \sum_{i=1}^N |\alpha_i - \gamma_i|, \qquad
D_2(G, H) = \sum_{i=1}^N (\alpha_i - \gamma_i)^2.
$$
Minimizing $D_2$ is equivalent to minimizing the Euclidean spectral distance $\sqrt{D_2}$. All candidate graphs in the proof are simple graphs of order $3n$ ($n \geq 3$) and may be disconnected.

- **Theorem 1 (Conjecture 20).** For every $n \geq 3$, $H_n = (n-1)K_1 \vee (K_1 \sqcup K_{n,n})$ is the unique simple graph of order $3n$ other than $K_{n,n,n}$, up to isomorphism, minimizing $D_1(K_{n,n,n}, G)$, and also the unique minimizer of $D_2(K_{n,n,n}, G)$.
- **Theorem 2 (Problem 21).** For every $n \geq 3$, every simple graph $G$ of order $3n$ (possibly disconnected) with $D_2(K_{n,n,n}, G) \leq 2/3$ has adjacency rank at most five.

## 2. Bounds from induced subgraphs and rank reduction

Write $N = 3n$ with $n \geq 3$. The adjacency spectrum of $K_{n,n,n}$ in nonincreasing order is $(\gamma_1, \dots, \gamma_N) = (2n, 0, \dots, 0, -n, -n)$. Abbreviate $D_1 = D_1(K_{n,n,n}, G)$ and $D_2 = D_2(K_{n,n,n}, G)$, and let $\alpha_1 \geq \cdots \geq \alpha_N$ be the spectrum of $G$.

**Lemma 3 (Interlacing bounds).** Let $F$ be an induced subgraph of $G$ of order $s \leq N - 2$ with spectrum $\beta_1 \geq \cdots \geq \beta_s$, and set $\beta_i^+ = \max(0, \beta_i)$ and $\beta_j^- = \max(0, -\beta_j)$. Then
$$
D_1 \geq 2\max\left(\sum_{i=2}^s \beta_i^+,\; \sum_{j=1}^{s-2} \beta_j^-\right), \qquad
D_2 \geq \sum_{i=2}^s (\beta_i^+)^2 + \sum_{j=1}^{s-2} (\beta_j^-)^2. \tag{1}
$$

*Proof.* Cauchy interlacing gives $\alpha_i \geq \beta_i \geq \alpha_{N-s+i}$ for $1 \leq i \leq s$. Since $s \leq N - 2$, every $i \in \{2, \dots, s\}$ lies in $\{2, \dots, N-2\}$ (where $\gamma_i = 0$), so whenever $\beta_i > 0$ we have $\alpha_i \geq \beta_i > 0$ and $(\alpha_i - \gamma_i)^+ \geq \beta_i^+$. Likewise, for $1 \leq j \leq s - 2$, the index $k = N - s + j$ lies in $\{3, \dots, N-2\}$ (where $\gamma_k = 0$), so whenever $\beta_j < 0$ we have $\alpha_k \leq \beta_j < 0$ and $(\alpha_k - \gamma_k)^- \geq \beta_j^-$. Because $\alpha_i > 0 > \alpha_k$, the positive and negative supports are disjoint, yielding the $D_2$ bound. Since both adjacency matrices have trace zero, $\sum_m (\alpha_m - \gamma_m)^+ = \sum_m (\alpha_m - \gamma_m)^- = D_1/2$, giving the $D_1$ bound. $\square$

**Lemma 4 (Symmetric elimination).** Let $M$ be a real symmetric matrix of rank $r$. Then there exists a nested sequence of invertible principal submatrices of orders $0 = k_0 < k_1 < \cdots < k_t = r$ with $k_\ell - k_{\ell-1} \in \{1, 2\}$. In particular, $M$ has an invertible principal submatrix of order $r$; if $r \geq 5$, it has one of order $5$ or $6$; and if $r \geq 6$, it has one of order $6$ or $7$.

*Proof.* Given an invertible principal submatrix $B$ of order $k < r$ with block form $\bigl(\begin{smallmatrix} B & C \\ C^T & D \end{smallmatrix}\bigr)$, the Schur complement $S = D - C^T B^{-1} C$ is symmetric of rank $r - k \geq 1$. A nonzero diagonal entry $s_{ii}$ yields an invertible principal submatrix of order $k + 1$; otherwise, a nonzero off-diagonal entry $s_{ij} = s_{ji}$ with zero diagonal yields a $2 \times 2$ pivot of determinant $-s_{ij}^2 \neq 0$, enlarging the order to $k + 2$. In both cases the remaining Schur complement has rank $r - \operatorname{order}(B') \geq 0$, so the sequence terminates at $r$. Starting from $k_0 = 0$, the first crossing of $5$ occurs at $5$ or $6$ (for $r \geq 5$), and the first crossing of $6$ occurs at $6$ or $7$ (for $r \geq 6$). $\square$

**Proposition 5 (Rank reduction and Theorem 2).** Theorem 2 holds. Moreover, for every $n \geq 3$,
$$
D_1(K_{n,n,n}, H_n) < 1 + \frac{2}{3n} \leq \frac{11}{9}, \qquad
D_2(K_{n,n,n}, H_n) < \frac{7}{12},
$$
and every graph $G$ of order $3n$ achieving $D_1 \leq D_1(K_{n,n,n}, H_n)$ or $D_2 \leq D_2(K_{n,n,n}, H_n)$ has adjacency rank at most four.

*Proof.* [verify_rank_reduction.py](verify_rank_reduction.py) checks all invertible adjacency matrices on $5$ and $6$ vertices up to isomorphism using exact rational root intervals for their integer characteristic polynomials, and certifies completeness and irredundancy of the graph atlas at each order $s$ by partitioning all $2^{\binom{s}{2}}$ labeled graphs into disjoint permutation orbits. Since $N - 2 = 3n - 2 \geq 7$ for all $n \geq 3$, Lemma 3 applies for $s \in \{5, 6, 7\}$. Every invertible graph of order $5$ or $6$ has lower bounds in (1) strictly exceeding $11/9$ and $7/12$, respectively. With `--problem21`, the script verifies that every invertible graph of order $6$ or $7$ has squared lower bound in (1) strictly exceeding $2/3$. By Lemma 4, any graph of order $3n$ ($n \geq 3$) with adjacency rank at least $6$ satisfies $D_2 > 2/3$, proving Theorem 2.

For the incumbent $H_n$, Theorem 6 and Remark 7 of the source give $D_1(K_{n,n,n}, H_n) < 1 + 2/(3n) \leq 11/9$. For $n \geq 5$, Lemma 5 and Theorem 6 of the source bound $D_2(K_{n,n,n}, H_n)$ by the maximum of the strictly convex function $2(x^2 - xy + y^2)$ on $0 \leq x \leq 17/30$, $2/5 \leq y \leq 1/2$, whose maximum over the four corners is $259/450 < 7/12$ (checked in [verify_rank_reduction.py](verify_rank_reduction.py)). Directly, the quotient polynomial of $H_n$ from Section 6 factors as $(x+n)q(x)$ with $q(x) = x^3 - nx^2 - (2n^2 - n - 1)x + n(n-1)$, whose signs at $-n, 0, 2/5, 1/2, 2n - 1/2 - 1/(3n), 2n - 17/30, 2n$ immediately yield these bounds. For $n \in \{3, 4\}$, [verify_rank_reduction.py](verify_rank_reduction.py) checks $D_2(K_{n,n,n}, H_n) < 7/12$ by rational root intervals. By Lemma 4, any graph of rank at least $5$ has $D_1 > 11/9$ and $D_2 > 7/12$, so any competitor achieving either incumbent distance has rank at most four. $\square$

## 3. Rank-at-most-four structure and small templates

Call a graph $G$ of order $3n$ *competing* if $D_1 < 11/9$ or $D_2 < 7/12$. By Proposition 5, every competing graph has rank at most four.

**Lemma 6.** Every competing graph $G$ contains a triangle, and therefore has adjacency rank $3$ or $4$.

*Proof.* A triangle-free graph of rank at most $4$ is bipartite: otherwise a shortest odd cycle $C_k$ has odd length $k \geq 5$, is induced, and has $k$ nonzero eigenvalues $2\cos(2\pi j/k)$ ($0 \leq j < k$), contradicting the rank bound. A bipartite graph has extreme eigenvalues $\rho, -\rho$ ($\rho \geq 0$, including $\rho = 0$ for the empty graph), giving $D_1 \geq |2n - \rho| + |-\rho - (-n)| \geq n \geq 3 > 11/9$ and $D_2 \geq (\rho - 2n)^2 + (n - \rho)^2 \geq n^2/2 \geq 9/2 > 7/12$. Hence no bipartite graph competes, so $G$ contains a triangle $K_3$ (rank $3$) and has rank $3$ or $4$. $\square$

**Lemma 7.** Let $G$ have adjacency matrix $A$ of rank $r$ with invertible principal submatrix $B$ of order $r$. Then $G$ is obtained from an induced twin-free template $T$ (containing the vertices of $B$) by replacing each vertex of $T$ with a nonempty independent set of false twins and adding zero or more isolated vertices.

*Proof.* Write $A = \bigl(\begin{smallmatrix} B & C \\ C^T & D \end{smallmatrix}\bigr)$. Since $\operatorname{rank}(A) = r$, the Schur complement vanishes, so $D = C^T B^{-1} C$. Each vertex's adjacency signature $v \in \{0,1\}^r$ to the seed satisfies $v^T B^{-1} v = 0$, and two signatures $v, w$ have adjacency $v^T B^{-1} w \in \{0,1\}$. All-zero signatures are isolated vertices; identical nonzero signatures are independent false twins. Since the $r$ columns of $B$ are distinct and nonzero, keeping one vertex per nonzero signature occurring in $G$ (including the seed) yields the induced twin-free template $T$. $\square$

Classifying the templates of competing graphs:

- **Rank 3.** Use a triangle $B = J_3 - I_3$ as the seed. Its inverse is $B^{-1} = \frac{1}{2}J_3 - I_3$, and for $v \in \{0,1\}^3$ of weight $w$, $v^T B^{-1} v = \frac{1}{2}w(w-2)$ vanishes only for $w = 0$ (isolates) and $w = 2$ (the columns of $B$). Thus $G$ is a complete tripartite graph plus zero or more isolated vertices.
- **Rank 4.** By Lemma 4, $G$ has an invertible principal seed $B$ of order $4$. [verify_conjecture20.py](verify_conjecture20.py) verifies the $4$-vertex graph atlas by permutation orbits and enumerates all $30$ compatible twin-free signature extensions (omitting isolates). Every resulting template has order at most $6 \leq N - 2$, so (1) applies for all $n \geq 3$. The only templates surviving at least one metric are the path $P_4$, the paw, the bull, and the house. The path $P_4$ is bipartite (excluded by Lemma 6). The other three consist of a triangle plus one or two extra vertices, each adjacent to a single triangle vertex; when there are two extra vertices, they attach to distinct triangle vertices and are nonadjacent (bull) or adjacent (house).

## 4. Eliminating templates with two or more outside vertices

In a blowup of the paw, bull, or house (with possible isolates), call the three parts corresponding to the triangle vertices the *triangle classes*, and call all other vertices *outside vertices*.

**Proposition 8.** Every competing graph $G$ is either a complete tripartite graph of full order $3n$ or a full-order paw blowup with a single pendant vertex (and no isolated vertices).

*Proof.* We show that every other surviving candidate satisfies $\rho(G) \leq r := 2n - 2/3$.

*Case 1 (paw blowup with at least two outside vertices, or bull/house blowup).* Construct an edge-supergraph $G'$ on the same $3n$ vertices consisting of a complete tripartite graph on $S = 3n - 2$ vertices and two outside vertices, each joined to one full triangle class and optionally to each other:
- For a paw blowup with at least two outside vertices, retain two outside vertices including at least one original pendant vertex (re-attaching the second to the same hub class if it was isolated), move any remaining pendant vertices into a triangle class different from their hub, and assign any remaining isolates to triangle classes.
- For a bull or house blowup, retain one vertex from each outside class. If the first outside class attaches to triangle class $A$ and the second to $B$, move the remaining first-class vertices into $B$, the remaining second-class vertices into $A$, and any isolates into the third triangle class $C$. Every original edge---including edges between the two outside classes in the house---is preserved.

By entrywise domination, $\rho(G) \leq \rho(G')$. Let $a_1, a_2, a_3 \geq 1$ ($\sum a_i = S = 3n - 2$) be the triangle class sizes of $G'$, $v = (\sqrt{a_1}, \sqrt{a_2}, \sqrt{a_3})^T$, $E \in \{0,1\}^{3 \times 2}$ the attachment matrix of the two outside vertices, and $F \in \{0,1\}^{2 \times 2}$ their adjacency matrix. The symmetrized quotient is
$$
Q = \begin{pmatrix}
vv^T - \operatorname{diag}(a_1, a_2, a_3) & \operatorname{diag}(v) E \\
E^T \operatorname{diag}(v) & F
\end{pmatrix}.
$$
Since $\|F\|_{\mathrm{op}} \leq 1$ and $r = 2n - 2/3 \geq 16/3 > 3$, we have $(rI_2 - F)^{-1} \preceq (r - 1)^{-1}I_2$. Writing $b_i \in \{0,1,2\}$ for the number of outside vertices attached to class $i$, we have $\sum_i b_i = 2$ and $EE^T = \operatorname{diag}(b_1, b_2, b_3)$. Taking the Schur complement in $rI - Q$, $rI - Q \succeq 0$ follows from
$$
\operatorname{diag}(r + k_i a_i) - vv^T \succeq 0, \qquad k_i := 1 - \frac{b_i}{r - 1} > 0,
$$
equivalently $\sum_{i=1}^3 a_i / (r + k_i a_i) \leq 1$. Setting $W = \sum_{i=1}^3 1/k_i$ and applying weighted Jensen to the concave function $x/(r+x)$ with weights $(1/k_i)/W$ at $x_i = k_i a_i$ gives the upper bound $WS/(Wr + S)$. Convexity of $t \mapsto (1 - t/(r-1))^{-1}$ on $[0, 2]$ yields
$$
W \leq 2 + \frac{1}{1 - 2/(r - 1)} = 2 + \frac{1}{1 - b/r}, \qquad b := \frac{2r}{r - 1},
$$
with $r - b = r(r - 3)/(r - 1) > 0$. Since $WS/(Wr + S)$ increases with $W$, it is at most $1$ if
$$
3r^2 - 2(b + S)r + bS \geq 0.
$$
At $n = 3 + m$ ($m \geq 0$), the left-hand side equals $4(3m + 2)(3m + 8)/(3(6m + 13)) > 0$ (checked in [verify_conjecture20.py](verify_conjecture20.py)). Thus $\rho(G) \leq \rho(G') \leq 2n - 2/3$.

*Case 2 (complete tripartite graph with at least one isolated vertex).* Here the live order is $S \leq 3n - 1$, and unweighted Jensen ($k_i = 1$, $W = 3$) at $r_0 = 2S/3$ gives $\sum_i a_i / (2S/3 + a_i) \leq 1$, hence $\rho(G) \leq 2S/3 \leq 2n - 2/3$.

*Elimination.* For any competing graph, the errors at the two smallest eigenvalues are strictly less than $\max(11/9, \sqrt{7/12}) < n$, so $\alpha_N \leq \alpha_{N-1} < 0 < \alpha_1$. With rank at most $4$, at most four coordinates $i$ have nonzero error $e_i = \alpha_i - \gamma_i$, including $e_1 = \alpha_1 - 2n \leq -2/3$. Since $\sum_i e_i = 0$, trace zero gives $D_1 \geq 2|e_1| \geq 4/3 > 11/9$, and Cauchy–Schwarz on the remaining at most three nonzero errors gives $D_2 \geq (4/3)e_1^2 \geq 16/27 > 7/12$. $\square$

## 5. Full-order complete tripartite graphs

For $n \in \{3, 4, 5\}$, [verify_conjecture20.py](verify_conjecture20.py) compares every full-order complete tripartite graph (excluding $K_{n,n,n}$) and every full-order single-pendant paw blowup (excluding $H_n$) directly against $H_n$ using exact rational root intervals ($17$, $35$, and $59$ competitors, respectively), confirming strict inequalities in both $D_1$ and $D_2$. Assume henceforth that $n \geq 6$.

**Proposition 9.** For every $n \geq 6$, no complete tripartite graph $K_{a,c,d}$ ($a + c + d = 3n$, $a,c,d \geq 1$) other than $K_{n,n,n}$ achieves $D_1 \leq D_1(K_{n,n,n}, H_n)$ or $D_2 \leq D_2(K_{n,n,n}, H_n)$.

*Proof.* If $K_{a,c,d}$ competes, both negative eigenvalue magnitudes lie strictly between $l := n - 7/9$ and $h := n + 7/9$ (since $11/18 < 7/9$ and $7/12 < (7/9)^2$). For its symmetric quotient $Q = vv^T - \operatorname{diag}(a, c, d)$, $lI + Q$ has two negative eigenvalues and $hI + Q$ is positive definite. Because adding $vv^T \succeq 0$ cannot increase the number of negative eigenvalues, at least two of $a, c, d$ exceed $l$, hence are at least $n$. At most one can exceed $n$, since a vector supported on two classes of size at least $n + 1 > h$ and orthogonal to $v$ would give a negative quadratic form for $hI + Q$.

Thus the sizes are $(n - k, n, n + k)$ in some order ($k \geq 0$). For $k \geq 2$, evaluating $hI + Q \succ 0$ on the vector orthogonal to $v$ supported on the classes of sizes $n$ and $n + k$ gives
$$
h > \frac{2n(n + k)}{2n + k} \geq n + \frac{n}{n + 1} > n + \frac{7}{9} = h,
$$
a contradiction. The case $k = 0$ is $K_{n,n,n}$.

The case $k = 1$ ($K_{n-1,n,n+1}$) is worse than $H_n$ in both metrics by Theorem 15 of the source. Directly, its $3 \times 3$ quotient has characteristic polynomial $f(x) = x^3 - (3n^2 - 1)x - 2n(n^2 - 1)$, with $f(-n + t) = n(1 - 3t^2) + t(t^2 + 1)$. Since $f(-n) = n > 0 > -2n(n^2 - 1) = f(0)$, we have $\alpha_N < -n < \alpha_{N-1} < 0 < \alpha_1$. Evaluating $f(-n + 1/\sqrt{3}) = 4/(3\sqrt{3}) > 0$ and $f(-n - 1/2) = n/4 - 5/8 > 0$ gives $\alpha_{N-1} + n > 1/\sqrt{3}$ and $-n - \alpha_N > 1/2$. For all $n \geq 6$,
$$
D_1 \geq 2(\alpha_{N-1} + n) > \frac{2}{\sqrt{3}} > \frac{10}{9} \geq 1 + \frac{2}{3n} > D_1(K_{n,n,n}, H_n),
$$
and
$$
D_2 \geq (\alpha_{N-1} + n)^2 + (\alpha_N + n)^2 > \frac{1}{3} + \frac{1}{4} = \frac{7}{12} > D_2(K_{n,n,n}, H_n). \quad \square
$$

## 6. Full-order paw blowups with one pendant vertex

Assume $n \geq 6$, and let $a \geq 1$ be the size of the triangle class adjacent to the pendant vertex and $c, d \geq 1$ the other two triangle class sizes, so $a + c + d = 3n - 1$. Here $(a, c, d) = (n - 1, n, n)$ gives $H_n$.

**Proposition 10.** For every $n \geq 6$, every full-order single-pendant paw blowup with $(a, c, d) \neq (n - 1, n, n)$ satisfies $D_1 > 11/9$ and $D_2 > 7/12$.

*Proof.* The induced paw has characteristic polynomial $(x + 1)(x^3 - x^2 - 3x + 1)$ (two positive and two negative roots), so by interlacing and rank at most $4$ the four nonzero eigenvalues satisfy $\alpha_1 \geq \alpha_2 > 0 > \alpha_{N-1} \geq \alpha_N$.

Eliminating the single pendant vertex in $rI - Q$ at $r := 2n - 12/25 > 1$ applies the Jensen bound of Proposition 8 with $S = 3n - 1$, $b = 1$, and $W = 2 + (1 - 1/r)^{-1}$ exactly. Since $3r^2 - 2(1 + S)r + S > 0$ for all $n \geq 3$ (checked by expanding in $n - 3$), we have $\rho(G) \leq 2n - 12/25$.

If $G$ competes, both negative magnitudes lie in $(l, h)$ where $l := n - 11/18$ and $h := n + 11/18$: for $D_1 < 11/9$ this follows from trace zero, and for $D_2 < 7/12$ it follows from the top-error contribution $(12/25)^2$ and $7/12 - (12/25)^2 < (11/18)^2$.

For $u \in \{l, h\}$ ($u > 0$), eliminating the pendant vertex in $uI + Q$ produces the Schur complement
$$
S(u) = \operatorname{diag}\bigl(u - a(1 + 1/u),\; u - c,\; u - d\bigr) + vv^T, \qquad v = (\sqrt{a}, \sqrt{c}, \sqrt{d})^T.
$$
Thus $S(l)$ has two negative eigenvalues and $S(h)$ is positive definite. Hence at least two of $a(1 + 1/l), c, d$ exceed $l$, and at most one of $a(1 + 1/h), c, d$ can exceed $h$. Moreover, evaluating $w^T S(h) w > 0$ on $w = (\sqrt{c}, -\sqrt{a}, 0)^T \perp v$ gives
$$
h > \frac{ac(2 + 1/h)}{a + c}, \tag{2}
$$
and similarly for $(a, d)$, while between $c$ and $d$ we have $h > 2cd/(c + d)$. Splitting on $a$:

- If $a \leq n - 2$, then $a(1 + 1/l) < n - 1 < l$, so both $c, d \geq n$ with at most one above $n$. Up to swapping $c, d$, they are $(n, n + k)$ with $k = n - 1 - a \geq 1$. For $k \geq 2$, $2n(n + k)/(2n + k) \geq n + n/(n + 1) > h$ is impossible, leaving $(a, c, d) = (n - 2, n, n + 1)$.
- If $a \geq n + 1$, then $a(1 + 1/h) > h$, so $c, d \leq n$ and at least one equals $n$ by the condition at $l$. The right-hand side of (2) increases with $a$ and already exceeds $h$ at $(a, c) = (n + 1, n)$ for all $n \geq 6$.
- If $a = n$, then $n(1 + 1/h) > h$ (i.e., $n - \frac{11}{18}h > 0$), so $c, d \leq n$ and $c + d = 2n - 1$, leaving $(n, n - 1, n)$ up to swapping $c, d$.
- If $a = n - 1$, then $c + d = 2n$. The right-hand side of (2) increases with $c$ and already exceeds $h$ at $(a, c) = (n - 1, n + 2)$ for all $n \geq 6$, leaving $(n - 1, n, n)$ (which is $H_n$) and $(n - 1, n - 1, n + 1)$.

The two exclusions from (2) and $n - \frac{11}{18}h > 0$ are verified by clearing positive denominators and expanding in $n - 6$ in [verify_conjecture20.py](verify_conjecture20.py).

The $4 \times 4$ quotient matrix has characteristic polynomial
$$
p(x) = x^4 - \bigl(a(1 + c + d) + cd\bigr)x^2 - 2acdx + acd.
$$
For each of the three surviving nonincumbent triples $(n - 2, n, n + 1)$, $(n - 1, n - 1, n + 1)$, and $(n, n - 1, n)$, expanding in $n - 6$ verifies $p(1/3) > 0$ and $p(-n + 1/2) < 0$ (all six signs are checked in [verify_conjecture20.py](verify_conjecture20.py)). Since $\alpha_1 \geq 1 > 1/3$ (the graph has an edge), $p(1/3) > 0$ forces $\alpha_2 > 1/3$. Since $-n + 1/2 < 0 < \alpha_2$, $p(-n + 1/2) < 0$ places $-n + 1/2 \in (\alpha_N, \alpha_{N-1})$, forcing $\alpha_{N-1} > -n + 1/2$. Both coordinates give positive errors, so trace zero yields
$$
D_1 > 2\left(\frac{1}{3} + \frac{1}{2}\right) = \frac{5}{3} > \frac{11}{9},
$$
and combining coordinates $1$, $2$, and $N - 1$ gives
$$
D_2 > \left(\frac{12}{25}\right)^2 + \left(\frac{1}{3}\right)^2 + \left(\frac{1}{2}\right)^2 = \frac{13309}{22500} > \frac{7}{12}.
$$
This excludes every nonincumbent paw blowup and completes the proof of Theorem 1. $\square$

## Computational verification

All finite checks use exact rational arithmetic and rational root intervals (no floating-point comparisons). Replay with Python 3.10+ (tested on Python 3.10 and 3.12 with NetworkX, NumPy, and SymPy):

```sh
python3 verify_rank_reduction.py --problem21
python3 verify_conjecture20.py
```

The recorded outputs are stored in [rank-certificate.jsonl](rank-certificate.jsonl) and [conjecture20-certificate.jsonl](conjecture20-certificate.jsonl).
