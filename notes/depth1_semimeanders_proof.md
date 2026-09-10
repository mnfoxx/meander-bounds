> **LITERATURE STATUS (added after checking DFGG 1996, hep-th/9607039; Sawada–Li 2012; Albert–Paterson 2005; Jensen 2000).**
> - §10 (escape numbers, Theorem H) is a rediscovery: DFGG define the winding *exactly* as our escape number (their §2.1, Fig. 4), state Theorem H as their eq. (2.11) "R(q) = R̄(q) iff ν(q) < 1", tabulate the winding distribution to n = 24, and estimate ν = 0.52(1). **Cut from any paper; cite.**
> - Lemma F (1 on top + n on bottom = open meander) is folklore (combos.org; Sawada–Li). Its parity remark is also wrong as written: odd n gives open meanders with an odd number of crossings, *not* closed meanders of the same height (closing adds an outer arc and raises the height by one). The counts coincide with A005315 only because A005316(2k−1) = A005315(k).
> - Theorem E is two lines in DFGG's language (bounded height ⇒ bounded winding ⇒ meander growth) plus standard truncation. **Remark, not theorem.**
> - Theorem C (sandwich) appears not to be written down anywhere, but it is elementary.
> - **Believed new:** the exact rational generating functions for height-bounded semi-meanders, closed meanders and open meanders (§3, §6, §11), the closed forms at height 2, and the phase-decomposition structure (§7, Propositions A–B). Jensen's transfer matrix uses width cutoffs only as a device for exact counts; Albert–Paterson's bounds use a word language, not height.

# Semi-meanders of bounded nesting depth: a rationality theorem and an exact formula for depth 1

## Statement

Fold a strip of $n$ stamps, labelled $1,\dots,n$ along the strip, into a single pile, and insist that stamp 1 lands on top. Call such a folding a *semi-meander* (this is the standard identification; the number of them is OEIS A000682, and the total number of foldings A000136 equals $n$ times it, since the top stamp is uniformly distributed).

Every crease between stamps $i$ and $i+1$ becomes a *fold* joining the two stamps. Say a fold is *nested inside* another if the two lie on the same side of the pile and the first sits strictly between the two ends of the second. The *depth* of a folding is the length of the longest chain of successively nested folds, minus one; so depth 0 means no fold sits inside another, and depth $\le 1$ means no three folds are mutually nested.

Write $s_D(n)$ for the number of semi-meanders of $n$ stamps of depth at most $D$.

**Theorem.**
1. For every fixed $D$, the generating function $\sum_n s_D(n)\,x^n$ is a rational function of $x$.
2. For $n\ge 3$,
$$
s_1(n)=\begin{cases}2F_{n+2}-3\cdot 2^{n/2} & n \text{ even},\\[2pt] 2F_{n+2}-2^{(n+3)/2} & n \text{ odd},\end{cases}
$$
where $F_k$ is the Fibonacci number ($F_1=F_2=1$). Equivalently
$$
\sum_{n\ge 3} s_1(n)x^n \;=\; \frac{x^2\,(1+x-x^2+2x^3+2x^4)}{(1-x-x^2)(1-2x^2)} \;-\;1-x^2 .
$$
The first values are $2,4,10,18,36,62,114,192,338,562,964,\dots$

The number of *all* foldings of $n$ labelled stamps with depth $\le 1$ is $n\cdot s_1(n)$.

---

## 1. The arch model

Number the positions in the finished pile $1,\dots,n$ from the top, and let $p(i)$ be the position of stamp $i$; $p(1)=1$. The crease between stamps $i$ and $i+1$ lies on one edge of the pile, and consecutive creases lie on opposite edges. Normalise so that crease $i$ is on the **right** when $i$ is odd and on the **left** when $i$ is even. Draw crease $i$ as an arc on that side joining positions $p(i)$ and $p(i+1)$.

**Fact (Sainte-Laguë; Koehler 1968; Lunnon 1968).** A permutation $p$ arises from a physical folding if and only if, on each side, the arcs are pairwise non-crossing: no two arcs on the same side have exactly one endpoint of one lying strictly inside the other.

So the right side carries the arcs $(1,2),(3,4),\dots$ and the left side carries $(2,3),(4,5),\dots$. Each stamp $2,\dots,n-1$ has exactly one arc on each side. Stamp 1 has only a right arc. Stamp $n$ has only one arc, on the right if $n$ is even and on the left if $n$ is odd. Hence:

* the right arcs form a non-crossing matching of all $n$ positions when $n$ is even, and of all positions except $p(n)$ when $n$ is odd;
* the left arcs form a non-crossing matching of all positions except $1$ when $n$ is odd, and of all positions except $1$ and $p(n)$ when $n$ is even.

**Lemma 1.** Semi-meanders of $n$ stamps are in bijection with pairs $(M_R,M_L)$ of non-crossing matchings on the positions $\{1,\dots,n\}$ such that $M_L$ leaves position 1 unmatched, exactly one further position (call it $q$) is unmatched by exactly one of $M_R,M_L$, every other position is matched by both, and the graph $M_R\cup M_L$ contains no cycle.

*Proof.* Given a folding, take its arcs; the union is the path $1-2-\dots-n$, which is acyclic, and the unmatched positions are as described with $q=p(n)$. Conversely, $M_R\cup M_L$ is a graph on $n$ vertices with maximum degree 2 and exactly $n-1$ edges (each matching contributes half its covered positions; the count above gives $n-1$). An acyclic graph with $n-1$ edges on $n$ vertices is a spanning tree, and with maximum degree 2 it is a Hamiltonian path; its two ends are the two degree-1 vertices, position 1 and $q$. Read the labels $1,2,\dots,n$ along the path starting from position 1. Since every interior vertex has one arc per side, sides alternate along the path, and the first arc (at stamp 1) is on the right; so the arc carrying label $i$ is on the right iff $i$ is odd. Non-crossing on each side then makes this a valid folding by the Fact. The two constructions are inverse. $\square$

**Depth in this picture.** Scan the positions from top to bottom, maintaining on each side a stack of *open* arcs (arcs whose upper end has been seen and whose lower end has not). Non-crossing is exactly the rule that the arc closed at any position is the top of its side's stack. A fold is at depth $d$ iff when it is opened the stack already holds $d$ arcs. So

$$\text{depth}\le D \iff \text{neither stack ever exceeds height } D+1.$$

---

## 2. The scan is a finite automaton (proof of part 1)

Process positions $1,\dots,n$ in order. At each position the stamp there performs one *event* on each side:

* **open**: push a new arc onto that side's stack;
* **close**: pop the top arc of that side's stack;
* **absent**: no arc on that side. This happens on the left at position 1, and at exactly one further position on exactly one side, namely position $q=p(n)$.

A position with both sides absent is impossible (it would be an isolated vertex).

After processing $t$ positions, the arcs seen so far form a graph of maximum degree 2 with no cycles (we reject as soon as a cycle would close), i.e. a disjoint union of paths. Each open arc is a loose end of exactly one such path. The path containing position 1 (the **root**) has one loose end (position 1 itself is a dead end). Once position $q$ has been seen, the path containing it (the **tail**) also has one loose end. Every other path has exactly two loose ends. Two paths are identified for our purposes if they have the same pattern of loose ends, so the **state** after $t$ steps is:

> the two stacks, each entry labelled by the path it belongs to (root, tail, or a two-ended path), up to renaming of the two-ended paths; together with a flag saying whether $q$ has been seen.

The transition rule is local: a close event connects the current stamp to the path owning the popped arc; two closes at one position connect the stamp to two paths, which is forbidden if they are the same path (cycle) and otherwise merges them; an absent event makes the current stamp's path the tail; open events push the current stamp's path label. The scan **accepts** at a position whose events close the root's loose end and the tail's loose end simultaneously (creating the single spanning path) while leaving both stacks empty; by Lemma 1 this happens exactly when we have read a semi-meander, and every semi-meander is read in exactly one way.

If both stacks are bounded by $D+1$, there are at most $2D+2$ loose ends, so finitely many states. Counting accepted scans of length $n$ is counting length-$n$ paths in a finite directed graph, and the generating function of such counts is rational (it is $u^{\!\top}(I-xT)^{-1}v$ for the adjacency matrix $T$). This proves part 1. $\square$

---

## 3. The depth-1 automaton (proof of part 2)

Now both stacks have height $\le 2$. Write $S$ for the side on which the root's loose end is currently open and $\bar S$ for the other side. By the left/right mirror symmetry the count of scans does not depend on which side $S$ is, so we may collapse each mirror pair of states into one.

### 3.1 Before the tail appears

Start (after position 1): the root has one open arc on the right, the left stack is empty. Exhausting the transition rule from this state and its successors while $q$ has not appeared yields exactly three states (stacks listed bottom to top; $c$ is a two-ended path):

| state | $S$ stack | $\bar S$ stack |
|---|---|---|
| $A$ | $[\text{root}]$ | $[\,]$ |
| $B$ | $[\text{root},\,c]$ | $[c]$ |
| $C$ | $[\text{root}]$ | $[c,\,c]$ |

and the transitions (events written as (event on $S$, event on $\bar S$)):

* $A \to A$ via (close, open): the stamp closes the root's arc and opens on the other side; $S$ flips.
* $A \to B$ via (open, open): a new two-ended path $c$ is opened, nested inside the root arc on $S$ and alone on $\bar S$.
* $B \to C$ via (close, open): the inner arc on $S$ closes and $c$ continues on $\bar S$, where it is now nested inside its own earlier arc. The pair (close, close) would join the two ends of $c$ (a cycle) and is rejected; opening on $S$ would give height 3.
* $C \to B$ via (open, close): $c$'s inner arc on $\bar S$ closes and $c$ re-enters $S$ inside the root arc.
* $C \to A$ via (close, close): the stamp joins the root to $c$; the merged root now has its single loose end on $\bar S$, so $S$ flips.

No other non-tail transitions exist (all others either exceed height 2, close a cycle, close an empty stack, or create an isolated vertex).

Let $a_k,b_k,c_k$ be the number of $k$-step scans from the start (i.e. having read positions $2,\dots,k+1$) ending in $A,B,C$, with $a_0=1,\,b_0=c_0=0$. Then
$$
a_{k+1}=a_k+c_k,\qquad b_{k+1}=a_k+c_k,\qquad c_{k+1}=b_k .
$$
With $A(x)=\sum a_kx^k$ etc. this reads $A=1+x(A+C)$, $B=x(A+C)$, $C=xB$, whence
$$
A(x)=\frac{1-x^2}{1-x-x^2},\qquad B(x)=\frac{x\,A(x)}{1-x^2}=\frac{x}{1-x-x^2},\qquad C(x)=\frac{x^2}{1-x-x^2}.
$$
(Indeed $c_k=a_{k-1}$, so $a_{k+1}=a_k+a_{k-1}$: the pre-tail phase is Fibonacci.)

### 3.2 The tail phase

Let $T$ denote the state: root's loose end open on one side, tail's loose end open on the other, nothing else open. From $T$:

* (close, close): root and tail join with empty stacks — **accept**. One position.
* (close on root's side, open on the other): the root re-opens beneath the tail on the tail's side. The only legal continuation is (open, close) on the next position, which returns to $T$. Two positions.
* (open on root's side, close on tail's side): symmetric; returns to $T$ in two positions.

Everything else is dead (a third path could never reconnect, or a stack overflows). So the generating function for completing the scan from $T$ is
$$
\tau(x)=\frac{x}{1-2x^2}.
$$

### 3.3 Entering the tail phase

The tail appears at position $q$, where one side's event is absent. From each pre-tail state, the legal ways to reach $T$ (or accept directly), with the number of positions used:

* From $A$: (close, absent) accepts immediately — one position, contributing $x$. (absent, open) reaches $T$ in one position. (open, absent) puts the tail's arc inside the root's on $S$; the forced next position (close, open) reaches $T$: two positions. Entry series: $x$ (direct accept) $+\;(x+x^2)\,\tau(x)$.
* From $B$: (close, absent) makes $c$ the tail with its remaining end on $\bar S$: reaches $T$ in one position. (absent, open) and (absent, close) each reach $T$ after one forced further position. Entry series: $(x+2x^2)\,\tau(x)$. Note (close, absent) with the root cannot occur here because the closing side's top arc is $c$, not the root.
* From $C$: (absent, close) makes $c$ the tail: one position to $T$. (open, absent) followed by the forced (close, close) reaches $T$ in two positions. (close, absent) would join root and tail while $c$ is still open, so it is rejected. Entry series: $(x+x^2)\,\tau(x)$.

### 3.4 Assembling

Position 1 contributes one factor $x$. Hence
$$
G(x)=\sum_{n} s_1(n)x^n
= x\Big[\,x\,A(x)+\big(A(x)(x+x^2)+B(x)(x+2x^2)+C(x)(x+x^2)\big)\,\tau(x)\Big].
$$
Substituting the closed forms for $A,B,C,\tau$ and simplifying,
$$
G(x)=\frac{x^2\,(1+x-x^2+2x^3+2x^4)}{(1-x-x^2)(1-2x^2)}
= 1+x^2+\frac{2(1+x)}{1-x-x^2}-\frac{3+4x}{1-2x^2}.
$$
Since $\sum_{n\ge0}F_{n+2}x^n=\frac{1+x}{1-x-x^2}$ and $\frac{3+4x}{1-2x^2}=\sum_m\big(3\cdot2^m x^{2m}+4\cdot 2^m x^{2m+1}\big)$, reading off the coefficient of $x^n$ for $n\ge3$ gives
$$
s_1(n)=2F_{n+2}-\begin{cases}3\cdot 2^{n/2}& n\text{ even}\\ 2^{(n+3)/2}& n\text{ odd}\end{cases}
$$
as claimed. (The $1+x^2$ term corrects the degenerate cases $n=0,2$; e.g. $s_1(2)=1$.) $\square$

---

## 4. Checks

* The 18-state unreduced automaton was built mechanically from the transition rule of §2 and, independently, the sequence $s_1(n)$ was computed by a brute-force folding enumerator that knows nothing about the arch model. The two agree for $3\le n\le 26$ (the automaton continues $995690, 1614928, 2627002,\dots$). This validates that the automaton — the one step a reader must check by hand — has been transcribed correctly.
* The reduced three-state derivation of §3 was checked symbolically against the unreduced automaton's generating function; they coincide.
* Asymptotically $s_1(n)\sim \tfrac{2}{\sqrt5}\,\varphi^{\,n+2}$, so depth-1 semi-meanders have growth constant $\varphi\approx1.618$, against the conjectured $\approx 3.5$ for all semi-meanders.

## 5. Remarks

* The same scan with stacks bounded by $D+1$ gives a mechanical route to $\sum_n s_D(n)x^n$ for any $D$; the state count grows quickly, and for $D=2$ nineteen computed terms did not fit any constant-coefficient recurrence of order $\le 8$, so the denominator there has degree at least 9.
* The structure of the proof — a Fibonacci pre-tail phase, a $(1-2x^2)$ tail phase, and a bounded interface between them — suggests that depth-$D$ generating functions may always factor as (a "spiral" part) $\times$ (a "tail" part) in a similar way; this is not proved here.

---

## 6. Addendum: depths 2 and 3, computed rigorously by the same automaton

Running the scan of §2 with stacks bounded by $D+1$ and extracting the minimal linear recurrence of the resulting integer sequence (Berlekamp–Massey on more than twice as many terms as the automaton has states, which certifies minimality) gives:

**Depth $\le 2$** (64 states, order-15 recurrence, all 19 brute-force terms reproduced):
$$
\sum_n s_2(n)x^n=\frac{x^2\,(1+x-9x^2-6x^3+34x^4+23x^5-54x^6-43x^7+42x^8+37x^9-28x^{10}-16x^{11}+14x^{12}+x^{13}-4x^{14})}{(1-x)\,(1-6x^2+8x^4-2x^6)\,(1-5x^2-x^3+5x^4-4x^6+x^8)}
$$
$s_2(n)$: $2,4,10,24,66,156,403,924,2273,5109,12190,\dots$; growth constant $\lambda_2\approx 2.1337$ (root of the degree-8 factor).

**Depth $\le 3$** (236 states, order-37 recurrence, all brute-force terms through $n=18$ reproduced). Denominator:
$$
(1-2x^2)\,(1-14x^2+69x^4-151x^6+164x^8-86x^{10}+16x^{12})\,\cdot\,P_{22}(x),
$$
where $P_{22}(x)=1-2x-14x^2+27x^3+77x^4-142x^5-224x^6+389x^7+397x^8-633x^9-473x^{10}+657x^{11}+414x^{12}-451x^{13}-280x^{14}+209x^{15}+143x^{16}-64x^{17}-51x^{18}+12x^{19}+11x^{20}-x^{21}-x^{22}$ (written with the sign convention $P(0)=1$). $s_3(n)$: $2,4,10,24,66,174,504,1350,3876,10300,28935,\dots$; growth constant $\lambda_3\approx 2.4833$.

**Observed structure at every depth computed.** The denominator factors as (tail factor, a polynomial in $x^2$) $\times$ (pre-tail factor), exactly as in the hand proof for $D=1$, where they are $1-2x^2$ and $1-x-x^2$. The tail factors have degree $1,3,7$ in $x^2$ for $D=1,2,3$. In each case the pre-tail factor carries the dominant root, so the growth constant $\lambda_D$ ($1.618,\ 2.134,\ 2.483$) is governed by the phase before the last stamp is placed. These $\lambda_D$ are rigorous lower bounds on the semi-meander growth constant, though far from the conjectured $\approx 3.5$.

---

## 7. The two phases, identified (Phase 2, step 3 of the project)

Let $q$ be the pile position of stamp $n$. Split the depth-$D$ automaton of §2 into the states before $q$ is read (pre-$q$) and after (post-$q$).

**Proposition A (phase decomposition).** The denominator of $\sum_n s_D(n)x^n$ divides $\chi_{\text{pre}}(x)\cdot\chi_{\text{post}}(x)$, where $\chi$ denotes the reversed characteristic polynomial of the sub-automaton's transfer matrix. Computed exactly:

* $D=1$: $\chi_{\text{pre}}=(1-x-x^2)(1+x-x^2)$, $\chi_{\text{post}}=(1-2x^2)^2$.
* $D=2$: $\chi_{\text{pre}}=(1-x)(1+x)\,P_8(x)P_8(-x)$ with $P_8=1-5x^2-x^3+5x^4-4x^6+x^8$; $\chi_{\text{post}}=(1-2x^2)^4(1-6x^2+8x^4-2x^6)^2$.

So the "pre-tail factor" of §6 is $\chi_{\text{pre}}$ (its mirror image $P(-x)$ cancels against the numerator, because the two mirror-symmetric copies of each pre-$q$ state feed the tail phase identically) and the "tail factor" is $\chi_{\text{post}}$. *Proof of divisibility:* the automaton has no edge from post-$q$ back to pre-$q$, so its transfer matrix is block upper-triangular. $\square$

**Why the tail factor is a polynomial in $x^2$.** After $q$ no event is absent, so every position changes the right-stack height by exactly $\pm1$. The post-$q$ sub-automaton is therefore bipartite (by parity of right-stack height) and its spectrum is symmetric under $\lambda\mapsto-\lambda$. The pre-$q$ automaton is bipartite too; its mirror-symmetric collapse is what produces odd factors like $1-x-x^2$.

**Proposition B (what the tail phase is).** The conjecture that post-$q$ states need no path labels is **false**: at $D=2$, 40 of the 46 post-$q$ states contain a two-ended path. The correct statement comes from reading the pile *upward from the bottom*. In that direction the post-$q$ region contains no single-ended path at all (the root and the tail only acquire their dead ends at positions 1 and $q$, which lie above). Hence the post-$q$ phase is, up to the direction of reading, the automaton

> $\mathcal{M}_D$: pairs of stacks of height $\le D+1$, every open arc belonging to a two-ended path, no absent events, no cycle ever closed,

started from empty stacks. This is precisely the height-bounded form of the transfer matrix used to enumerate **closed meanders** (the pile is the river; the two sides are the two half-planes). Built directly, $\mathcal{M}_D$ has $5,12,34$ states for $D=1,2,3$ and reversed characteristic polynomials
$$
1-2x^2,\qquad 1-6x^2+8x^4-2x^6,\qquad (1-2x^2)(1-4x^2+5x^4-x^6)(1-14x^2+\dots+16x^{12}),
$$
which contain the tail factors of §6 exactly (at $D=3$ the extra factor $1-4x^2+5x^4-x^6$ cancels in the semi-meander generating function; it corresponds to states unreachable from the semi-meander's tail configuration). This is the sense in which "no labels are needed": the tail phase is not label-free, but it is *root-free*, and that is what makes it a meander problem rather than a semi-meander problem.

**Consequence: a reading of the shrinking gap.** Write $\mu_D$ for the tail growth constant and $\lambda_D$ for the pre-tail one:

| $D$ | $\mu_D$ | $\mu_D^2$ | $\lambda_D$ | $\lambda_D-\mu_D$ |
|---|---|---|---|---|
| 1 | 1.4142 | 2 | 1.6180 | 0.204 |
| 2 | 2.0529 | 4.214 | 2.1337 | 0.081 |
| 3 | 2.4488 | 5.997 | 2.4833 | 0.035 |
| 4 | 2.6997 | 7.288 | 2.7219 | 0.022 |

Since $\mathcal{M}_D$ is the height-$(D+1)$ closed-meander matrix, $\mu_D^2$ is the growth constant per pair of positions of height-bounded open-meander systems, and (if the bounded-height constants converge to the unbounded one) $\mu_D^2\to R\approx12.26$, the meander constant, while $\lambda_D\to\bar R\approx3.50$, the semi-meander constant. Di Francesco, Golinelli and Guitter conjectured $\bar R=\sqrt R$. The observed collapse of the gap $\lambda_D-\mu_D$ is exactly what that conjecture predicts inside the bounded-depth family, so the table is (modest) independent numerical support for $\bar R^2=R$, and it isolates *where* the identity would come from: the pre-$q$ phase (one dead end) and the post-$q$ phase (no dead ends) must have the same per-position entropy in the limit.

**Upper bound on the tail phase.** Dropping the no-cycle constraint, $\mathcal{M}_D$ is dominated by two independent Dyck paths of height $\le D+1$, whose joint growth per position is $4\cos^2\!\big(\tfrac{\pi}{D+3}\big)$: $2,\,2.618,\,3,\,3.247$ for $D=1,\dots,4$. The true $\mu_D$ ($1.414,\,2.053,\,2.449,\,2.700$) sit well below, so the cycle constraint is not a perturbation; a closed form for $\mu_D$, if one exists, will not come from Chebyshev-type formulas alone.

**Status of the project bet.** Phase 2 step 3 as posed is refuted; the corrected version (tail phase $\equiv$ bounded-height meander matrix) is proved. Step 4 becomes: any exact result about height-bounded closed-meander transfer matrices transfers immediately to rigorous lower bounds on $\bar R$ via $\mu_D\le\bar R$, and the $\bar R=\sqrt R$ conjecture is equivalent to $\lim\lambda_D=\lim\mu_D$ (given convergence in $D$).

---

## 8. The two phases have the same limiting growth (Theorem C)

Notation. $\mathcal S_D$ is the pre-$q$ automaton of §7 (one single-ended path, the root; stacks of height $\le D+1$) and $\lambda_D$ its spectral radius. $\mathcal M_D$ is the rootless automaton of Proposition B (all paths two-ended; stacks $\le D+1$) and $\mu_D$ its spectral radius. Both automata read the same alphabet: at each position, one event (open/close) per side. A *walk* is a sequence of events valid from a given state.

Computed values:

| $D$ | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| $\mu_D$ | 1.4142 | 2.0529 | 2.4488 | 2.6997 | 2.8739 | 2.9960 |
| $\lambda_D$ | 1.6180 | 2.1337 | 2.4833 | 2.7219 | 2.8867 | — |

**Theorem C.** For every $D\ge 2$,
$$\mu_{D-1}\;\le\;\lambda_D\;\le\;\mu_{D+1}.$$
Consequently $\mu_\infty:=\lim_D\mu_D$ exists and $\lim_D\lambda_D=\mu_\infty$.

*Proof.* Both inequalities are injections of event sequences.

*(i) $\lambda_D\le\mu_{D+1}$.* Let $X$ be any state of $\mathcal S_D$ and $w$ a walk of length $k$ from $X$. Form $X'$ from $X$ by inserting a *frozen* arc at the bottom of the right stack, belonging to the root's path; the root's path now has two loose ends (its live end and the frozen arc), so $X'$ is a state of $\mathcal M_{D+1}$. We claim $w$ is valid from $X'$ in $\mathcal M_{D+1}$. Heights: the right stack is taller by one, so bounded by $D+2$. Closes: $w$ never closes the right stack below its rooted height, so it never touches the frozen arc; every other close pops the same arc as before. Cycles: a cycle in $\mathcal M_{D+1}$ would close both ends of one path; for the root's path one end is frozen, and every other path is unchanged, so no new cycle arises. Merges are identical. Hence $w\mapsto w$ is an injection from walks of $\mathcal S_D$ from $X$ into walks of $\mathcal M_{D+1}$ from $X'$. Since $\mathcal S_D$ is strongly connected (any state can be closed down to the start state and rebuilt), the number of length-$k$ walks from $X$ grows like $\lambda_D^k$, while the number of length-$k$ walks in any finite automaton from any state is $O(k^m\mu^k)$ with $\mu$ its spectral radius. So $\lambda_D\le\mu_{D+1}$.

*(ii) $\mu_{D-1}\le\lambda_D$.* Let $Y$ be a state in the dominant strongly connected class of $\mathcal M_{D-1}$ and $w$ a walk of length $k$ from $Y$. Form $Y'$ by inserting a frozen *single-ended* arc — a root — at the bottom of the right stack. $Y'$ is a state of $\mathcal S_D$ (stacks $\le D+1$). The same three checks apply: heights fine, the frozen root is never popped, and the only new cycle possibility would involve the root, which never closes. So walks of $\mathcal M_{D-1}$ from $Y$ inject into walks of $\mathcal S_D$ from $Y'$, and $\mu_{D-1}\le\lambda_D$.

$\mu_D$ is non-decreasing in $D$ ($\mathcal M_D\subseteq\mathcal M_{D+1}$) and at most $4$ (four event pairs per position), so $\mu_\infty$ exists, and the sandwich gives $\lambda_D\to\mu_\infty$. $\square$

The data respect both inequalities with room: $\lambda_D$ sits inside $(\mu_{D-1},\mu_{D+1})$ at every computed $D$, and in fact $\mu_D<\lambda_D<\mu_{D+1}$, which the proof does not give (the lower inequality $\mu_D\le\lambda_D$ is open; it would follow from a version of (ii) that does not lose a level of height).

### What Theorem C says about the meander constants

Let $\bar R$ be the semi-meander growth constant and $R$ the closed-meander constant per order (so $\sqrt R$ per crossing), with Jensen's estimates $\bar R\approx 3.50$, $R\approx12.26$. Di Francesco–Golinelli–Guitter conjectured $\bar R=\sqrt R$; the easy direction $\sqrt R\le\bar R$ follows from $M_n\le\bar M_{2n}$ (a meander is a semi-meander of zero winding), and the direction $\bar R\le\sqrt R$ is open.

Depth-$\le D$ semi-meanders grow like $\max(\lambda_D,\mu_D)$ and height-$\le D{+}1$ meanders grow like $\mu_D$ per crossing (Proposition B), so
$$\mu_\infty\;\le\;\min\big(\bar R,\sqrt R\big).$$
Theorem C therefore reduces the open direction to one statement:

> **(*)** *Depth-bounded semi-meanders capture the full semi-meander entropy:* $\sup_D\lambda_D=\bar R$.

Indeed (*) gives $\bar R=\mu_\infty\le\sqrt R$, and with the easy direction, $\bar R=\sqrt R$.

(*) is a concrete question about a single family of objects, with no field theory in it. It asks whether a semi-meander of large nesting depth can be replaced by one of bounded depth at sub-exponential cost — a "depth compression" statement of the kind that is a theorem for self-avoiding walks in strips and for bounded-height Dyck paths, and that fails for objects whose entropy lives entirely in a single deep spiral. The data in the table ($\mu_6\approx3.00$ against $\bar R\approx3.50$, with increments shrinking roughly geometrically) are consistent with (*) but far from decisive.

Two smaller points fall out for free. Every $\mu_D$ and $\lambda_D$ is a rigorous algebraic lower bound on $\bar R$ (currently $\bar R\ge 2.996$ from $D=6$, much weaker than $\sqrt{11.38}\approx 3.37$ from the meander bounds of Albert–Paterson via the easy direction, so of no independent use yet). And the mirror-symmetric collapse in §3 shows why the depth-1 pre-tail factor is $1-x-x^2$ rather than a polynomial in $x^2$: the bipartite structure is broken only by the way the root's side alternates.

---

## 9. The depth-bounded constants converge to the meander constant (Theorem E)

Let $\mathcal S_\infty$ be the pre-$q$ automaton with unbounded stacks: a countable, irreducible, non-negative matrix $T$, of which each $\mathcal S_D$ is the principal truncation to states of height $\le D+1$. Let $i$ be the start state ($R=[\text{root}]$, $L=[\,]$).

**Lemma D (truncation).** $\displaystyle\sup_D\lambda_D=\limsup_{n\to\infty}\big((T^n)_{ii}\big)^{1/n}.$

*Proof.* Write $g$ for the right-hand side. Each $\mathcal S_D$ is strongly connected and contains $i$, so $\lambda_D=\limsup_n((\mathcal S_D^n)_{ii})^{1/n}\le g$, since every walk of $\mathcal S_D$ is a walk of $T$. Conversely, fix $\varepsilon>0$ and $n$ with $(T^n)_{ii}\ge(g-\varepsilon)^n$. Finitely many walks are counted, all of height $\le D+1$ for some $D$, so $(\mathcal S_D^{\,n})_{ii}=(T^n)_{ii}$; concatenating, $(\mathcal S_D^{\,kn})_{ii}\ge((\mathcal S_D^{\,n})_{ii})^k$, whence $\lambda_D\ge(g-\varepsilon)$. $\square$

(This is the elementary half of the Vere-Jones truncation theorem for infinite non-negative matrices; no more is needed.)

**What the diagonal counts.** A walk $i\to i$ of length $n-2$ followed by the event (close, absent) is a semi-meander of $n$ stamps in which stamp $n$ sits at the *bottom* of the pile with its arc on the right; the mirror start state gives the left case. Conversely every semi-meander with $p(n)=n$ arises this way. Let $b_n$ be the number of semi-meanders of $n$ stamps with stamp 1 on top and stamp $n$ at the bottom. By mirror symmetry and strong connectivity, $\limsup b_n^{1/n}=g$.

**Lemma F.** $b_n$ equals the number of open meanders with $n-2$ crossings (OEIS A005316):
$$b_n=1,2,3,8,14,42,81,262,538,1828,3926,13820,30694,110954,252939\quad(n=4,\dots,18),$$
verified by brute force; for odd $n$ these are the closed meander numbers A005315.

*Proof.* Delete stamps 1 and $n$. Their arcs were hairpins to positions 2 and $n-1$ (a stamp at the top or bottom of the pile has only one neighbour), so nothing else changes. What remains is the path $2-3-\dots-(n-1)$ drawn as non-crossing arcs alternating sides of the pile line, crossing it $n-2$ times, with both ends free: the end at position 2 can leave upward past the (now empty) top, and the end at $n-1$ downward. That is an open meander with the pile as the river. The construction is reversible: given an open meander with $n-2$ crossings, order the crossings along the river, put a new stamp above the first and below the last, and both are hairpins. Sides alternate along the road, so the labelling is forced. $\square$

**Lemma G (standard).** Open meanders with $m$ crossings grow like $(\sqrt R)^{m}$, where $R\approx12.26$ is the closed meander constant per order. *Proof.* Closed$(k)\le$ Open$(2k)$: cut an outermost arc of a closed meander; both new ends escape to infinity. Open$(m)\le (m{+}1)\cdot$Closed$(\lceil m/2\rceil)$: join the two free ends around the outside (adding one crossing if they lie in opposite half-planes); the open meander is recovered from the closed one together with the marked closing arc. $\square$

**Theorem E.** $\displaystyle\lim_{D\to\infty}\lambda_D=\lim_{D\to\infty}\mu_D=\sqrt R.$

*Proof.* Lemmas D, F, G give $\sup_D\lambda_D=\sqrt R$; $\lambda_D$ is non-decreasing in $D$, so it converges to $\sqrt R$; Theorem C carries this to $\mu_D$. $\square$

So the growth constants of depth-bounded semi-meanders — $1.618,\ 2.134,\ 2.483,\ 2.722,\ 2.887,\dots$ — converge to the *meander* constant $\sqrt R\approx3.50$, not, a priori, to the semi-meander constant $\bar R$. Extended values of $\mu_D$ (rootless automaton, power iteration): $\mu_7=3.0864$ (2528 states), $\mu_8=3.1542$ (8123 states). Convergence is polynomial in $D$, roughly $\sqrt R-\mu_D\sim D^{-p}$ with $p$ drifting upward from 1 to 1.5 over the computed range; unconstrained fits give a limit of $3.6\pm0.1$, consistent with $\sqrt{12.26}=3.501$.

### The conjecture, restated exactly

Since every depth-bounded semi-meander is a semi-meander, $\sqrt R=\sup_D\lambda_D\le\bar R$ (recovering the known easy direction). Statement (*) of §8 — that depth-bounded semi-meanders capture the full entropy — is therefore *literally equivalent* to $\bar R=\sqrt R$; §8's reduction is an identity, not a reduction. What the last two sections buy is a clean and testable form of the conjecture:

> $\bar R=\sqrt R\iff b_n/s(n)$ decays sub-exponentially $\iff$ semi-meanders with stamp $n$ at the bottom of the pile are not exponentially rare.

In matrix language: the pre-$q$ transfer matrix $T$ is *$R$-recurrent-like* at its convergence parameter — the total mass $\sum_X(T^n)_{iX}$ (weighted by the tail phase) grows no faster than the diagonal $(T^n)_{ii}$ — rather than transient. Equivalently, a uniformly random semi-meander does not typically hold linearly many open arcs at the level of its last stamp.

**Evidence.** From exact enumeration for $10\le n\le 18$:
$$\frac{b_n}{s(n)}\sim n^{-\gamma},\qquad \gamma=1.47,1.44,1.45,1.43,1.44,1.42,1.43$$
(local exponents on successive windows), i.e. a clean power law with no sign of an exponential factor. If $s(n)\sim\bar R^{\,n}n^{-\bar\alpha}$ and open meanders scale as $(\sqrt R)^n n^{-\alpha'}$, then under $\bar R=\sqrt R$ one expects $\gamma=\alpha'-\bar\alpha$; with the Di Francesco–Golinelli–Guitter values $\bar\alpha\approx2.05$ and the meander exponent $\approx3.42$ this predicts $\gamma\approx1.37$, which the slowly decreasing local exponents are heading toward. Mean nesting depth of a uniform semi-meander grows from $2.24$ at $n=10$ to $3.90$ at $n=18$, sublinearly but not yet clearly so.

**Where a proof would have to come from.** One needs $b_n\ge s(n)\,e^{-o(n)}$: an injection (or a bounded-to-one map) from semi-meanders of $n$ stamps into open meanders with $n+o(n)$ crossings. Read the semi-meander from both ends: the region above $q=p(n)$ is a rooted walk down from stamp 1, and the region below $q$, read upward, is a rooted walk up from stamp $n$; the two are glued at level $q$ through the arcs spanning it. A proof of the conjecture is exactly a way to "unzip" this gluing at sub-exponential cost — to route the two roots past the spanning arcs so they exit the top and bottom of the pile. The obstruction is that an arc spanning level $q$ blocks the path on its side; the path can only slip past it by alternating sides, and each such step may enter a deeper nest on the other side. Whether the spanning arcs at a typical level can always be escaped in $o(n)$ steps is, after this reduction, the whole content of the conjecture.

---

## 10. Escape numbers: the conjecture as a statement about winding (Theorem H)

**Definitions.** For a semi-meander $F$ of $n$ stamps, an *extension* by $k$ stamps is a semi-meander of $n+k$ stamps whose first $n$ stamps fold as $F$ (deleting the last $k$ stamps of any folding leaves a folding, so restriction is well defined). The *escape number* $e(F)$ is the least $k$ for which some extension has stamp $n+k$ at the bottom of the pile — equivalently (Lemma F) the least number of crossings that must be added to turn $F$ into an open meander. This is the *wind-factor* of Sawada–Li (Electron. J. Combin. 19 (2012) P43), who use it to filter open meanders out of semi-meanders; it is finite for every $F$ (the spiral is the extreme case, with $e=n-1$). Write $s_{\le k}(n)$ for the number of semi-meanders of $n$ stamps with $e\le k$, and recall $b_m=$ open meanders with $m-2$ crossings.

**Theorem H.** For all $n,k$:
$$s_{\le k}(n)\;\le\;\sum_{j=0}^{k} b_{n+j}.$$
Consequently, if there is a sequence $k(n)=o(n)$ with $s_{\le k(n)}(n)\ge s(n)\,e^{-o(n)}$, then $\bar R\le\sqrt R$, and hence $\bar R=\sqrt R$.

*Proof.* Fix $j$ and, for each $F$ with $e(F)=j$, choose one extension $F'$ of $n+j$ stamps with the last stamp at the bottom. Then $F'$ is counted by $b_{n+j}$, and $F$ is recovered from $F'$ by deleting its last $j$ stamps, so $F\mapsto F'$ is injective on $\{e=j\}$. Summing over $j\le k$ gives the inequality. For the consequence: by Lemma G and Fekete's lemma for closed meanders, $b_m\le C_\varepsilon(\sqrt R+\varepsilon)^m$ for every $\varepsilon>0$; so
$$s(n)\,e^{-o(n)}\le s_{\le k(n)}(n)\le (k(n)+1)\,C_\varepsilon(\sqrt R+\varepsilon)^{\,n+k(n)},$$
and taking $n$-th roots with $k(n)=o(n)$ gives $\limsup s(n)^{1/n}\le\sqrt R+\varepsilon$. The reverse inequality $\sqrt R\le\bar R$ is the easy direction of §9. $\square$

**Corollary.** If the escape number of a uniformly random semi-meander of $n$ stamps is $o(n)$ in probability — in particular if its mean is $o(n)$ — then $\bar R=\sqrt R$.

*Proof.* By Markov, $\mathbb P(e\le\varepsilon n)\ge 1-\mathbb E[e]/(\varepsilon n)\to1$, so $s_{\le\varepsilon n}(n)\ge s(n)/2$ for large $n$; apply Theorem H with every $\varepsilon$ and let $\varepsilon\to0$. $\square$

Conversely, $\bar R=\sqrt R$ implies (granting that $s(n)^{1/n}$ converges, as it does for meanders by Fekete) that already $b_n=s_{\le0}(n)\ge s(n)e^{-o(n)}$: the configurations with escape number zero carry the full entropy. So the three statements

* $\bar R=\sqrt R$;
* open meanders are a sub-exponential fraction of semi-meanders;
* some sub-linear escape threshold captures almost all semi-meanders

are equivalent. This makes rigorous the heuristic behind the conjecture: Di Francesco, Golinelli and Guitter wrote that $R=\bar R^2$ "is a consequence of the polymer folding interpretation" (same entropy per monomer with a free end as with a closed loop), and studied the winding of semi-meanders — the number of arches above the source — finding numerically a winding exponent $\nu$ with $\tfrac12\le\nu<1$, i.e. sub-linear typical winding. The escape number is the combinatorial cost of unwinding, and the Corollary is the precise sense in which sub-linear winding *implies* the equality of constants.

**Data (exact, $n\le13$).** Distribution of $e$ over all semi-meanders, $e=0,1,2,\dots$:

| $n$ | $s(n)$ | mean $e$ | distribution |
|---|---|---|---|
| 8 | 174 | 2.517 | 14, 28, 39, 61, 14, 16, 1, 1 |
| 10 | 1406 | 2.883 | 81, 181, 276, 476, 150, 190, 24, 26, 1, 1 |
| 12 | 12198 | 3.209 | 538, 1290, 2098, 3870, 1485, 2031, 374, 436, 36, 38, 1, 1 |
| 13 | 37378 | 3.337 | 1828, 2098, 9894, 6151, 9667, 2835, 3645, 548, 622, 43, 45, 1, 1 |

The mean grows like $n^{\nu}$ with local exponent $0.66,0.63,0.61,0.60,0.59,0.58$ on successive windows, drifting down toward the $\tfrac12$ end of the DFGG range — sub-linear, as the Corollary needs, though small $n$ cannot exclude a late crossover. Several exact patterns are visible and easy to conjecture: $e=n-1$ and $e=n-2$ are attained exactly once each (the spiral and its unique near-relative); the number with $e=n-3$ is $\lfloor (n+1)^2/4\rfloor-4$ for $6\le n\le13$; and the $e=1$ column at $n$ equals the $e=2$ column at $n-1$ ($181,276,1290,2098$), which suggests a simple bijection worth writing down.

**Status.** The unconditional results of this project are Theorems A–E (rational generating functions for every depth bound, the depth-1 closed form, the phase decomposition, and $\lim\lambda_D=\lim\mu_D=\sqrt R$). Theorem H does not prove the conjecture; it identifies the one quantitative statement — sub-linear escape numbers for typical semi-meanders — that would, and shows it is the same statement DFGG measured as $\nu<1$. Proving $\nu<1$ is the problem.


---

## 11. Meanders of bounded height (closed and open): exact results

Height means the maximum number of arcs stacked on either side. Closed meanders of order $k$ (i.e. $2k$ crossings) are counted from the rootless automaton $\mathcal M_{h-1}$ started at empty stacks with a final cycle-closing step; open meanders with $m$ crossings are walks of length $m$ in the rooted pre-$q$ automaton from the start state back to the start state or its mirror. All counts below were verified against independent brute force (closed: all pairs of non-crossing matchings, $k\le7$; open: foldings with stamp 1 on top and $n$ at the bottom, $n\le14$).

**Closed meanders, height $\le2$:** exactly $2^{k-1}$ for $k\ge1$; generating function $y/(1-2y)$ with $y=x^2$.

**Closed meanders, height $\le3$:** $1,2,8,34,144,608,2564,10808,45552,191976,809056,\dots$
$$\sum_k c^{(3)}_k y^k=\frac{y(1-2y)^2}{1-6y+8y^2-2y^3}.$$

**Closed meanders, height $\le4$:** $1,2,8,42,242,1436,8596,51578,309616,1858406,\dots$
$$\frac{y\,(1-6y+8y^2-4y^3)(1-8y+17y^2-11y^3)}{(1-2y)(1-14y+69y^2-151y^3+164y^4-86y^5+16y^6)}.$$
The true meander numbers are $1,2,8,42,262,1828,13820$; height $\le4$ already accounts for all meanders of order $\le4$ and $242/262$ of order 5.

**Open meanders, height $\le2$:** the Fibonacci numbers: $1,1,2,3,5,8,13,\dots$ for $m\ge1$ crossings; generating function $(1-x^2)/(1-x-x^2)$.

**Open meanders, height $\le3$:** $1,1,2,3,8,14,35,65,155,300,694,1382,3126,6349,14126,\dots$
$$\frac{(1-x-2x^2+x^4)(1+x-2x^2+x^4)}{(1-x)(1-5x^2-x^3+5x^4-4x^6+x^8)}.$$

**Open meanders, height $\le4$:** $1,1,2,3,8,14,42,81,243,496,1454,3084,8813,19234,53717,\dots$, with the same degree-22 denominator $P_{22}$ as the depth-3 semi-meander pre-tail factor of §6.

The denominators are the pre-$q$ characteristic polynomials of §7, as they must be: the open-meander generating function is a diagonal entry of $(I-x\mathcal S_D)^{-1}$. The closed-meander denominators are the $\mathcal M_D$ characteristic polynomials. This is the cleanest statement of the phase decomposition: the semi-meander generating function at depth $D$ is assembled from the open-meander and closed-meander generating functions at the same height.

The closed forms at height 2 — $2^{k-1}$ closed meanders, $F_m$ open meanders, and $2F_{n+2}-c\cdot2^{n/2}$ semi-meanders — are the three faces of one five-state automaton and should be presented together, with the semi-meander formula as the glued case.

---

## 12. An improved rigorous lower bound on the meander growth constant (Theorem J)

Let $M_k$ be the number of closed meanders of order $k$ ($2k$ crossings) and $R=\lim_k M_k^{1/k}$ (the limit exists by Fekete's lemma, Lando–Zvonkin). Jensen's estimate is $R\approx12.2629$; the best rigorous bounds in the literature are $11.380\le R\le12.901$ (Albert–Paterson 2005), still cited as best known in Delecroix–Goujard–Zograf–Zorich, Ann. Inst. Fourier (2024).

**Theorem J.** $R\ \ge\ 11.54507$.

### The automaton

Scan the river left to right. Between consecutive crossings, the state is the list of open arcs above the river (bottom to top of the stack) and below it, each labelled by the path it belongs to, up to renaming of paths; every path has exactly two loose ends. At each crossing the road either opens a new arc or closes the innermost open arc on each side (non-crossing forces the innermost). Two closes joining the two ends of the same path make a cycle; that is forbidden except at the last crossing, when both stacks are otherwise empty, which produces a closed meander. The *width* of a state is the total number of open arcs. Let $\mathcal W_W$ be this automaton restricted to states of width $\le W$; walks from the empty state to the accepting step of length $2k$ are in bijection with closed meanders of order $k$ whose width never exceeds $W$. (This is the transfer matrix of Jensen (2000), truncated.) The implementation was checked against a brute-force enumeration over all pairs of non-crossing matchings: widths $\le4,6,8$ give $38,184,890$ / $42,258,1700$ / $42,262,1824$ closed meanders of orders $4,5,6$, matching in every case.

Width truncation converges much faster per state than the per-side height truncation of §§7–9: at $\approx$33 000 states, width $\le16$ gives $3.250$ where height $\le10$ gives $3.207$.

### From a certificate to a bound

Let $T$ be the transition matrix of $\mathcal W_W$ on its non-empty states (the empty state has no incoming edge and is dropped). $T$ is non-negative and every state is reachable from the first state; from the one-path state $j_0=((1),(1))$ the accepting step is available immediately. Suppose $v\ge0$, $v\ne0$, is an integer vector with
$$(T^{\!\top})^2 v\ \ge\ c\,v\quad\text{entrywise}.$$
Then $(T^{\!\top})^{2m}v\ge c^m v$, so for every $m$ some state $i$ has at least $c^m v_{j_0}/\|v\|_1$ walks of length $2m$ from $i$ to $j_0$. Prefixing a fixed walk from the start to $i$ (length $\ell_i\le L$, $L$ depending only on $W$) and appending the accepting step gives, for each $m$, a closed meander count of some order $k\in[m+1,m+L/2+1]$ that is at least $c^m v_{j_0}/\|v\|_1$. Since $M_k^{1/k}\to R$ and width-$\le W$ meanders are meanders, $R\ge c$. (For the odd/even bookkeeping: the automaton is bipartite in the parity of the upper-stack height, so $(T^{\!\top})^2$ is the natural object and $c$ bounds the growth per *order*, i.e. per two crossings.)

The vector $v$ is obtained by power iteration in floating point, then rounded **down** to integers at scale $2^{52}$; the product $(T^{\!\top})^2v$ and the minimum ratio are computed in exact 64/128-bit integer arithmetic. No floating-point quantity enters the certificate.

### Results

| $W$ | states | $\rho(T^2)$ (power iteration) | certified $c$ | $R\ge$ |
|---|---|---|---|---|
| 20 | 478 192 | 11.0416 | 134080992/12143478 | 11.04139 |
| 22 | 1 830 270 | 11.2085 | — | 11.20795 |
| 24 | 7 030 570 | 11.3434 | 32406063/2857120 | 11.34221 |
| 26 | 27 088 870 | 11.4538 | 10994240961/959947840 | 11.45295 |
| 28 | 52 327 345 orbits ($\approx1.05\cdot10^8$ states) | 11.5453 | 4367716392/378318673 | **11.54507** |

The certified increments ($0.167,0.134,0.109,0.092$) shrink by a factor $\approx0.83$ per step, so $W=30$ ($\approx4\cdot10^8$ states) and $W=32$ should certify roughly $11.62$ and $11.68$; both are routine on a machine with 64 GB of memory.

**The mirror quotient.** The $W=28$ computation used the quotient of the automaton by the swap of the two sides. Let $P$ be the swap, which commutes with $T$; orbits have size 1 or 2. For orbits $a,b$ let $E(a,b)$ be the number of transitions from a fixed representative of $a$ into states of $b$; this is what the code stores. By $T[Ps][Pt]=T[s][t]$, the number of transitions from *all* states of $a$ into a fixed $t\in b$ is $E(a,b)\,|a|/|b|$. Hence if $y>0$ on orbits satisfies $(E^{\!\top})^2y\ge c\,y$, the lift $v_t:=y_{b}/|b|$ for $t\in b$ satisfies $(T^{\!\top})^2v\ge c\,v$ on the full automaton: writing $\tilde u_b=\sum_a E(a,b)y_a$ one checks $\tilde u_b=|b|\,(T^{\!\top}v)_t$, and iterating once more gives the claim, the orbit sizes cancelling. So the certificate on the quotient is a certificate on $T$, and the argument of the previous paragraph applies verbatim. As a check, the quotient reproduces the full-automaton certificates at $W=20$ and $W=24$ to the last digit.

The $W=28$ eigenvector was computed in five checkpointed legs of $\approx170$ power iterations; the certified value rose $11.5026\to11.5395\to11.5437\to11.5447\to11.5451$ while the power-iteration estimate settled at $11.54530$.

**Semi-meanders.** Since a closed meander of order $k$ is a semi-meander with $2k$ bridges, $\bar R\ge\sqrt{R}\ge3.39780$; this is the "easy direction" of §9 and inherits the improvement.

### Supplementary facts for §12

* **Width equals twice the number of open paths** in the closed-meander automaton (every path has two loose ends), so truncating by number of paths is the same truncation. The number of states with exactly $p$ open paths is $\binom{2p-1}{p}$ ($1,3,10,35,126,462,1716,\dots$ observed for $p\le12$), i.e. half a central binomial coefficient; this is the count of Jensen's boundary configurations and is presumably known.
* **Where the Perron mass sits.** At $W=24$ the stationary mass (left $\times$ right Perron vectors) by width is $0.3,1.2,3.3,6.9,11.6,15.9,18.3,17.7,13.9,8.2,2.8$ percent for widths $4,6,\dots,24$: it peaks near width $16$ and about $3\%$ sits on the boundary width, which is why each increment of $W$ gains a steady $\approx0.83$ of the previous gain rather than converging fast.
* **Consistency with Jensen.** Fitting $\rho_W=R-aW^{-p}$ to the last four or five values gives $R\approx12.40$–$12.46$, drifting downward as earlier terms are dropped, consistent with Jensen's $12.2629$; the fit is not a claim, only a check that nothing is wrong with the truncation.
* **Exact generating functions at small width** (closed meanders, $y=x^2$ marks the order):
  width $\le4$: $1,2,8,38,184,890,4296,\dots$, $\;y(1-3y)^2/(1-8y+17y^2-8y^3)$ (14 states);
  width $\le6$: $1,2,8,42,258,1700,11500,\dots$, denominator $1-26y+285y^2-1735y^3+6506y^4-15751y^5+25107y^6-26432y^7+18204y^8-7972y^9+2052y^{10}-240y^{11}$ (49 states).
  Compare height $\le2$ (both sides): $2^{k-1}$, and height $\le3$: $y(1-2y)^2/(1-6y+8y^2-2y^3)$ from §11. Width $\le4$ contains height $\le2$ strictly (it allows a $3{+}1$ split), and its growth $\sqrt{4.83}$ already exceeds height $\le3$'s $2.05$ per crossing... (per order: $4.83$ vs $4.21$).
* **Memory.** At $\approx58$ bytes per orbit-state, $W=30$ ($\approx2\cdot10^8$ orbits) needs $\approx12$ GB and $W=32$ about $45$ GB; the code (`src/muq.c`) checkpoints the eigenvector and needs no change.
---

## 13. Adaptive truncation: selecting states by Perron mass instead of by width (Theorem K)

§12 truncated the closed-meander automaton to states of width $\le W$. The mass diagnostic there (Perron mass peaks near width $16$ at $W=24$, only $\approx3\%$ on the boundary) suggested that a width cut keeps many states that contribute almost nothing while discarding, at the boundary, a few that matter. This section replaces the width cut by a *selected* set of states, chosen adaptively, and shows that the resulting certificates are better per state by roughly one full width step: with $30\cdot10^6$ orbit-states it already beats the $W=28$ certificate of §12 ($52\cdot10^6$ orbit-states), and $58\cdot10^6$ states certify $11.5975$, and with $14\cdot10^6$ it beats $W=26$ ($13.5\cdot10^6$) by $0.065$ where the $W=26\to28$ step is $0.092$.

**Theorem K.** $R \ge 11.60805$. *(Certificate: $55\,781\,199$ orbit-states, all of width $\le32$, with the integer vector $v$; $c=803642589184/69231468544$; computed on the author's laptop, see the addendum in §13.3. Files `s58_states.bin`, `s58_vec.bin` (1.6 GB, GitHub release v1.0); verified by the independent checker `src/check.c`. An earlier cloud run gave $R\ge11.59753$ on $58\,398\,791$ states, `results/adaptive/s58_check.log`.)*

### 13.1 Sub-automaton certificates

The argument of §12 never used that the truncated automaton was defined by a width bound. Let $\mathcal A$ be the full closed-meander automaton (mirror quotient, as in §12) and $S$ any finite set of orbit-states containing the one-path state $j_0=((1),(1))$. Write $T_S$ for the transition matrix restricted to $S$ (only transitions with both ends in $S$).

**Lemma K1.** Let $S'\subseteq S$ be the set of states reachable from $j_0$ by walks inside $S$. If $v\ge0$ is an integer vector on $S'$ with $v_{j_0}>0$ and $(T_{S'}^{\!\top})^2v\ge c\,v$ entrywise, then $R\ge c$.

*Proof.* Iterating, $(T_{S'}^{\!\top})^{2m}v\ge c^mv$; the $j_0$ entry says $\sum_i v_i\,\#\{\text{walks } i\to j_0 \text{ of length } 2m \text{ inside } S'\}\ge c^m v_{j_0}$, so some $i\in S'$ has at least $c^mv_{j_0}/\|v\|_1$ such walks. Prefix the walk $\varnothing\to j_0$ (one crossing, event (open,open)) and a fixed walk $j_0\to i$ inside $S'$ of length $\ell_i\le L$ (finite, as $S'$ is finite), and append the accepting step at $j_0$. Each of the counted walks becomes a distinct closed meander of order $k=(2+\ell_i+2m)/2$, and $\ell_i$ is even because the automaton is bipartite (parity of the upper-stack height) and $i$ and $j_0$ lie in the same class whenever a walk of even length joins them. So for every $m$ some order $k\in[m+1,\,m+1+L/2]$ has $M_k\ge c^m v_{j_0}/\|v\|_1$, and $M_k^{1/k}\to R$ gives $R\ge c$. The quotient lemma of §12 lifts a certificate on orbits to one on states unchanged, since $S'$ is a union of orbits. $\square$

The checker therefore needs only the list $S$ and the vector $v$: it re-derives the transitions of each listed state from scratch, computes $S'$ by breadth-first search from $j_0$, evaluates $(T_{S'}^{\!\top})^2v$ in exact integers, and reports the minimum ratio. No property of how $S$ was chosen enters the proof. Any heuristic whatever may be used to choose $S$.

### 13.2 The selection heuristic

For an irreducible non-negative matrix with left and right Perron vectors $a,b$ (normalised $a^{\!\top}b=1$), $\partial\rho/\partial T_{ij}=a_ib_j$, so deleting state $i$ (its row and column) changes $\rho$ to first order by $-(2\rho-T_{ii})\,a_ib_i$. The product $m_i=a_ib_i$ — the Perron *mass* of $i$, already used as a diagnostic in §12 — is thus the natural importance of a state, and a state outside $S$ is scored by the same product with its Perron entries estimated from the neighbours it has inside $S$:
$$\hat a_t=\frac1\rho\sum_{i\in S}a_iT_{it},\qquad \hat b_t=\frac1\rho\sum_{j\in S}T_{tj}b_j,\qquad \text{score}(t)=\hat a_t\hat b_t .$$
This is the selection rule of the CIPSI family of "selected configuration interaction" methods in quantum chemistry (Huron–Malrieu–Rancurel 1973), transplanted to a non-symmetric transfer matrix; the only change is that left and right vectors differ.

The procedure, as run (`src/grow.c`):

1. Start from all states of width $\le W_0$ ($W_0=20$, 239 359 orbit-states).
2. Power-iterate $T_S^2$ from a warm start to get $a,b$ (float32; the certificate is recomputed in exact integers, so floating point only affects how tight $c$ is, never its validity).
3. Enumerate the one-step frontier $F=\{t\notin S:\ T_{it}>0 \text{ for some } i\in S\}$, capped at width $32$, and score it. Keep the top $K$ of $S\cup F$ by score (a log-scale histogram supplies the threshold, so nothing beyond $S$ and one chunk of $F$ is ever held in memory).
4. Rebuild the transitions inside the new $S$ and repeat.

Two things matter in practice. First, since only a one-step frontier is scored, the state set must be grown *and pruned repeatedly*: a schedule of "expand, expand, prune to $K$" cycles with $K$ growing geometrically ($1,2,3.5,3.5,3.5,7,14,14\cdot10^6$) produces a far better core than a single expansion from a width cut — at $13$–$15\cdot10^6$ states, the once-expanded width-24 set certifies $11.450$ while the annealed core certifies $11.485$. Second, the *expansion* of a good core (core plus its full one- or two-step frontier) is itself an excellent state set: it costs no selection and adds exactly the states whose absence caps the pruned core.

### 13.3 Results

All values are rigorous certificates ($c$ from the exact integer computation; the power-iteration estimate $\rho_S$ of the true Perron root of $T_S^2$ is shown for comparison). Width-cut rows are from §12.

| state set | orbit-states | $\rho_S$ | certified $c$ | $R\ge$ |
|---|---|---|---|---|
| width $\le 24$ | 3.52 M | 11.3487 | | 11.342 |
| adaptive, 3 cycles | 3.50 M | 11.4102 | 955416461312/83807846400 | 11.40008 |
| width $\le 26$ | 13.5 M | 11.4538 | | 11.45295 |
| adaptive, pruned | 7.02 M | 11.4616 | | 11.436 (50 it.) |
| adaptive, pruned | 14.05 M | 11.5186 | | 11.489 (50 it.) |
| width $\le 28$ | 52.3 M | 11.5453 | 4367716392/378318673 | **11.54507** (§12) |
| adaptive, core 14 M + 2-step frontier | 55.27 M | 11.5897 | 193002846208/16661409792 | 11.58382 (150 it.) |
| adaptive, pruned | 30.07 M | 11.5529 | 347897499648/30126884864 | 11.54774 |
| adaptive, pruned (final of run) | 40.06 M | 11.5969 | 136407019520/11775848448 | 11.58363 |
| adaptive, core 30 M + frontier, 550 it. (cloud) | 58.40 M | 11.5990 | 3500486230016/301830176768 | 11.59753 |
| adaptive, pruned (laptop reproduction, 10 Sep) | 30.03 M | 11.5677 | 233799364608/20218146816 | 11.56384 |
| adaptive, core 30 M + frontier, 400 it. (laptop) | 55.78 M | 11.6084 | 803642589184/69231468544 | **11.60805** |

(Entries marked "$n$ it." are certificates taken after only $n$ warm-started power iterations at intermediate stages of the run; they are valid bounds but not converged. The cloud row is the stage-10 set of the run re-iterated to convergence, `./grow -1 32 400 s58` from its checkpoint.)

**Addendum, 10 September 2026 (laptop reproduction; this is now Theorem K and Theorem 1.1 of the paper).** The author re-ran the full schedule `grow 20 32 150 s58 -1 -1 1000000 -1 -1 2000000 -1 -1 3500000 -1 -1 3500000 -1 -1 3500000 -1 -1 7000000 -1 -1 14000000 -1 -1 14000000 -1 30000000 …` on an Apple-silicon MacBook Pro in one unbroken pass, interrupted it at the 30 M core (STAGE 26: 30,034,282 states, certified 11.55835 after 150 it.), and resumed with one full-frontier expansion, `grow -1 32 400 s58 -1`. The core expanded to 55,781,199 states (137,423,195 transitions; widths 30 and 32 carry 17.75 M and 16.26 M states against 19.83 M and 13.61 M in the cloud set) and after 400 iterations certified $c = 803642589184/69231468544 = 11.6080534\ldots$, with all 55,781,199 entries positive. The independent checker `check.c` on the same laptop reported all states well formed and canonical, all reachable from $j_0$, the same transition count, and the same rational (`results/adaptive/s58_check_laptop.log`). The improvement over the cloud certificate (fewer states, higher bound) is attributed to the core: the cloud set was assembled over several resumed runs with fewer iterations per stage, so its 30 M core was slightly less well annealed. Two independent runs of the same construction agreeing to within 0.01 is recorded here as a consistency check. Headline: $R \ge 11.60805$, hence $\bar R \ge 3.4070$.

Per state, then, the adaptive sets are worth about one width step: $14\cdot10^6$ adaptive states ($11.519$) sit between width $26$ ($11.454$, $13.5\cdot10^6$) and width $28$ ($11.545$, $52\cdot10^6$), closer to the latter; $58\cdot10^6$ adaptive states give $11.599$, roughly what §12's extrapolation predicted for width $30$ ($\approx11.62$ at $4\cdot10^8$ states) minus a third of a step. The same machine budget therefore reaches the equivalent of $W\approx29.5$.

The $58\cdot10^6$ selected states are distributed over widths as $\dots,\ 22{:}0.67\,\mathrm M,\ 24{:}2.3\,\mathrm M,\ 26{:}6.8\,\mathrm M,\ 28{:}14.9\,\mathrm M,\ 30{:}19.8\,\mathrm M,\ 32{:}13.6\,\mathrm M$. Against the full layer sizes ($0.68,\ 2.6,\ 10,\ 39,\ \approx150,\ \approx600$ M) the set keeps essentially all of width $\le22$, $89\%$ of width 24, $68\%$ of 26, $38\%$ of 28, and — this is where the gain comes from — about $13\%$ of width 30 and $2\%$ of width 32, layers a width cut can only buy wholesale.

### 13.4 What this does and does not change

It does not change the convergence law: the Perron mass migrates outward as the state set grows (a typical meander of large order has many open arcs across a typical cut), so no finite selection converges faster than polynomially, and the selected sets still gain a steady fraction of the previous gain per doubling. What it changes is the constant: for a given memory the bound is one width step better, and — because the method has no $W$ — a 64 GB machine can run it directly at $\approx 4\cdot10^8$ states, which by the table above should certify $R\ge11.66$–$11.68$ rather than the $\approx11.62$ of a plain $W=30$ run.

Reproduce: `gcc -O2 -o grow src/grow.c -lm; ./grow 20 32 150 run -1 -1 1000000 -1 -1 2000000 -1 -1 3500000 -1 -1 3500000 -1 -1 3500000 -1 -1 7000000 -1 -1 14000000 -1 -1 14000000 -1 30000000 40000000 40000000` (about 3 h and 5 GB on two cores; `-1` = add the whole one-step frontier, a positive $K$ = prune to $K$). Verify: `gcc -O2 -o check src/check.c; ./check results/adaptive_55M/s55` (the checker's state model — points on the cut line with a non-crossing partner matching — is written independently of `grow.c`'s stack-of-labels model, and reproduces its edge counts exactly).
---

## 14. Upper bounds by forgetting: $R\le 12.6319$ and the first upper bound on the semi-meander constant (Theorem L)

The transfer matrix of §12 is exact because every open arc carries the label of its path, and the only thing the labels are used for is to refuse the one move that would close a cycle. *Forgetting* labels therefore produces an automaton that accepts every meander and some non-meanders, and its growth rate is a rigorous **upper** bound. With a finite memory the automaton is finite, its Perron root is computable, and — unlike the truncations of §§12–13 — the bound improves as memory grows without any state ever being discarded.

**Theorem L.** Let $R$ be the closed-meander constant and $\bar R=\lim s(n)^{1/n}$ the semi-meander (stamp-folding) constant. Then
$$R\ \le\ 12.63185,\qquad \bar R\ \le\ \sqrt{12.63185}\ =\ 3.55413 .$$
*(Certificates: relaxed automaton with total memory $26$, $72\,808\,296$ orbit-states (5.9 GB), integer vector $w>0$ with $T^2w\le c\,w$, $c=13888861889491/2^{40}$; `src/relax.c`, cross-checked by the independent Python build `src/relax_py.py` at memory $\le16$ and by the simulation check `src/relax_check.py` at width $\le12$.)*

The previous rigorous upper bound on $R$ was $12.901$ (Albert–Paterson 2005, cluster method over forbidden factors of length $\le16$); Jensen's estimate is $12.2629$. For $\bar R$ (Lunnon's conjecture $\bar R = 3.5$, Jensen's estimate $3.5019$) we have found no rigorous upper bound in the literature below the trivial $4$ (Legendre's 2013 survey lists none), and no route from $R$ to $\bar R$ exists without the open conjecture $\bar R=\sqrt R$; the bound above is proved directly.

### 14.1 The relaxed automaton

A *partial state* consists of an upper and a lower stack whose entries are either a label or the symbol $?$; every label occurs exactly twice in total. Read a meander as in §12 with events $(e_U,e_L)\in\{\text{open},\text{close}\}^2$ at each crossing:

* a close pops its stack; it is **refused** only if both events are closes and the two popped entries are the *same label* (a certain cycle);
* if both popped entries are labels $a\neq b$, every remaining $b$ is renamed $a$ (the two paths merge); if one popped entry is $?$, nothing is renamed;
* an opened arc receives the label of the popped arc if exactly one arc was popped and it was labelled, and a fresh label otherwise; two arcs opened at once share a fresh label;
* then the state is **reduced**: every label now occurring once becomes $?$; leading $?$ at the bottom of each stack are deleted; and while the total number of entries exceeds the memory $m$, the deepest entry of the taller stack is deleted (repeating the two previous steps).

The reduced state has at most $m$ entries and at most $m/2$ labels, so the automaton $\mathcal R_m$ is finite (mirror quotient: $58\,043$ orbit-states at $m=16$, $4.23\cdot10^6$ at $m=22$, $17.6\cdot10^6$ at $m=24$, $72.8\cdot10^6$ at $m=26$).

**Lemma L1 (soundness).** Say a partial state $r$ *covers* an exact state $s$ (a state of the full automaton of §12, of any width) if, on each side, the entries of $r$ are the top $|r|$ entries of $s$ with $?$ matching anything and equal labels in $r$ lying on the same path in $s$. Then for every exact transition $s\to t$ under an event $e$ and every $r$ covering $s$, the event $e$ is accepted by $\mathcal R_m$ from $r$ and leads to a state covering $t$.

*Proof.* The event is refused by $\mathcal R_m$ only when both popped entries carry the same label, which by covering means the same path in $s$ — exactly the move the exact automaton refuses. Merging renames $b\mapsto a$ only when $a,b$ are labels of arcs that the road has just joined, so equal labels still mean equal paths; a $?$ pop renames nothing, and an arc whose partner was thereby joined to an unknown path is by construction no longer equal-labelled to anything it is not on a path with. A new arc labelled from a popped labelled arc lies on that arc's path in $t$; a fresh label is equal only to its twin (two arcs opened together, on one path). Reduction only replaces labels by $?$ or deletes entries, which preserves covering. $\square$

The reduction step "a label occurring once becomes $?$" loses nothing: such a label can never trigger a refusal, and its only other role (being renamed on a merge) has no observable effect. It cuts the state count by a factor $4$–$5$ at equal memory.

**Corollary.** For every closed meander of order $k$, the first $2k-1$ events of its event word are accepted by $\mathcal R_m$ from the empty state (start with $r=\varnothing$, which covers everything, and induct with Lemma L1; the last event, which closes the single loop, is the one refusal that is *right* in the exact automaton and is simply not taken as a transition). Distinct meanders have distinct event words, so $M_k\le 4\,\|T_m^{2k-1}\|_\infty$ and $R\le\rho(T_m)^2$. For semi-meanders, Lemma 1 of §1 reads a folding of $n$ stamps as a word of $n$ events with the same alphabet at $n-2$ positions and a single-side event (the *absent* event of §2) at position $1$ and at the position $q$ of the last stamp; the root and the tail are single-ended paths and never close a cycle, so treating them as $?$ is a covering. Hence $s(n)\le n\,\|A\|_\infty^2\,\|T_m^{n-2}\|_\infty$ with $A$ the single-event matrix, and $\bar R\le\rho(T_m)$. $\square$

The certificate is the Collatz–Wielandt inequality: for any strictly positive $w$, $\rho(T^2)\le\max_i (T^2w)_i/w_i$. We take $w=\lceil 2^{40}v\rceil\vee 1$ from the power-iterated right vector and evaluate $T^2w$ in exact integers.

Computationally the simulation of Lemma L1 was also verified mechanically: over all exact states of width $\le12$ and all their transitions, and all covering partial states, for memories $m\le 8$ and per-side memories $\le6$, no violation (`src/relax_check.py`).

### 14.2 Results

| memory rule | orbit-states | $\rho(T_m^2)$ | certified $R\le$ | $\bar R\le$ |
|---|---|---|---|---|
| per side $k=1$ | 3 | 13.9282 | 13.9282 | 3.7320 |
| per side $k=4$ | 103 | 13.1412 | 13.1412 | 3.6251 |
| per side $k=7$ | 4 301 | 12.8899 | 12.8899 | 3.5903 |
| per side $k=10$ | 50 213 | 12.7680 | 12.7680 | 3.5733 |
| per side $k=12$ | 716 861 | 12.7158 | 12.7158 | 3.5660 |
| per side $k=14$ | 10 469 223 | 12.6768 | 12.6768 | 3.5605 |
| total $m=16$ | 58 043 | 12.7425 | 12.7425 | 3.5697 |
| total $m=18$ | 243 356 | 12.7120 | 12.7120 | 3.5654 |
| total $m=20$ | 1 016 670 | 12.6869 | 12.6869 | 3.5619 |
| total $m=22$ | 4 233 616 | 12.6656 | 12.6656 | 3.5589 |
| total $m=24$ | 17 579 062 | 12.6475 | 12.64749 | 3.55633 |
| total $m=26$ | 72 808 296 | 12.6318 | **12.63185** | **3.55413** |

(Per-side rows use the older variant without the singleton reduction for $k\le8$ where both were run; the bounds coincide, as they must.) Already $k=7$ — $4\,301$ states, a few seconds — beats the 2005 bound. The total-memory rule is better per state than the per-side rule (memory $20$ with $10^6$ states beats per-side $12$ with $7\cdot10^5$ by $0.03$; memory $22$ with $4.2\cdot10^6$ beats per-side $14$ with $10.5\cdot10^6$). Two other rules — forgetting the *pair* with the deepest shallow end, or with the largest depth sum — were tried and are worse per state, because $?$ entries then accumulate in the middle of the stacks.

The decrements $0.030,\,0.025,\,0.021,\,0.018,\,0.016$ per two units of memory shrink by a factor $\approx0.85$, so this family converges to something near $12.5$, not to $R$: forgetting deep arcs permits deep cycles, and typical arch systems have a positive density of those. A better relaxation must keep *some* information about deep arcs; that is the natural next question (cf. §14.4).

### 14.3 What the two bounds say together

$$11.6080\ \le\ R\ \le\ 12.6319,\qquad 3.4070\ \le\ \bar R\ \le\ 3.5542,$$
with $R\approx12.2629$, $\bar R\approx3.5019$ (Jensen). The $\bar R$ interval contains $\sqrt R\approx3.5018$, as the conjecture $\bar R=\sqrt R$ requires; note that its upper end is *below* $\sqrt{12.901}=3.592$, i.e. the semi-meander bound could not have been obtained from the old meander bound even if $\bar R\le\sqrt R$ were known.

### 14.4 Remarks and directions

* **Adaptive forgetting.** Any deterministic, information-losing post-processing gives a sound relaxation (Lemma L1 needs only that reduction preserves covering). One can therefore forget more in states of low Perron mass and less in states of high mass, the upper-bound analogue of §13; we did not pursue it, because the geometric slowdown above suggests the gain is a constant factor in memory, not a change of limit.
* **Prime meanders (remark, probably known — check La Croix 2003).** Lando–Zvonkin's sewing $A\cdot B$ (cut the upper arc at the last crossing of $A$ and the upper arc at the first crossing of $B$; reconnect crosswise), together with its lower-arc mirror, has unique factorisation: writing $P(x)=\sum P_kx^k$ for meanders that are not a product ($P_k=1,0,4,18,126,916,7204,59618,\dots$), every meander is uniquely a sequence of primes joined by upper or lower sewings, so
$$M(x)=\frac{P(x)}{1-2P(x)}\qquad(\text{verified for }k\le7\text{ by brute force}).$$
Since $M(1/R)<\infty$ ($M_k\sim R^kk^{-3.42}$), $2P(1/R)<1$: the pole is never reached, $P$ and $M$ have the same radius, and $P_k/M_k\to$ a constant ($0.64$ at $k=28$, rising). Consequently concatenation carries no entropy, and the lower bound $R\ge1/x^*$ with $2\sum_{k\le K}P_kx^{*k}=1$ that this identity yields is weak ($8.45$ at $K=28$) — a small negative result worth one sentence in the paper.

### 14.5 Attempts to keep partial information about deep arcs (all negative or marginal)

Recorded so they are not repeated.

* **Tilting by width.** Since a meander's event word returns to width $0$, $M_k$ is also bounded by the weighted walk count in which every (open,open) carries $\theta^{2}$ and every (close,close) $\theta^{-2}$, for any $\theta>0$; on the relaxed automaton this is *not* a diagonal similarity (the true width is not a function of the reduced state), so $\min_\theta\rho(T_\theta)$ could in principle beat $\theta=1$. It does not: at $m=14$, $\rho^2=12.78$ at $\theta=1$ and $12.89/12.92$ at $\theta=0.95/1.05$ (`src/tilt.py`). The forgetting and the phantom-close boundary effects balance exactly at $\theta=1$.
* **Coarse identities cannot forbid.** Any relaxation of the form "superset automaton" can only refuse a move when it is *certain* the two closed arcs are one path; colours, parities, side-of-partner bits, or nesting information about forgotten arcs never certify equality, so they forbid nothing. Weighting uncertain closures by a loop fugacity $q\in(0,1]$ (an upper bound since $M_k\le Z_k(q)/q$) is monotone in $q$ and optimal at $q=0$, i.e. the plain relaxation.
* **Forgetting rules.** Deleting the pair with the deepest shallow end, or with the largest depth sum, instead of the deepest single entry, is worse per state at every memory tried (holes of $?$ in the middle of a stack multiply the state count faster than the retained pairs pay for). Per-side and total budgets: total is better (§14.2).
* **Adaptive coarsening** (the upper-bound analogue of §13): build $\mathcal R_m$, compute the Perron mass, and additionally truncate every state whose mass is in the lowest fraction $f$. At $m=18$, $f=0.7$ gives $73\,000$ states and $12.7288$, against $12.7425$ for $m=16$ ($58\,000$ states) and $12.7120$ for the full $m=18$ ($243\,000$): a gain of about a third of a memory step at fixed size. Extending memory to $m+2$ only in high-mass regions (mass taken from the $m$ reference) is worse than uniform memory at equal size. Marginal; `src/relax_adapt.c`.
* **Where the Perron mass of $\mathcal R_m$ sits.** Not at the boundary: at $m=16$, $47\%$ of the mass is on fully-labelled states and only $12\%$ on states with $m$ entries; $5.6\%$ sits on the *empty* reduced state, from which all four events are phantom-allowed. The excess of $\rho(\mathcal R_m)^2$ over $R$ is therefore not a boundary effect but the accumulated ignorance of pairings created during wide excursions (`src/relaxmass.py`).
* **Why the convergence is slow, heuristically.** A loop is invisible to memory $m$ as soon as the width exceeded $m$ at some time during its life. A spurious system can therefore hide one extra loop per $O(m)$ crossings at bounded entropy cost, which suggests $\rho(\mathcal R_m)^2-R\asymp 1/m$ rather than geometric decay; the observed decrements $0.021,0.018,0.016$ at $m=22,24,26$ are consistent with $\sim m^{-2}$ ($(24/26)^2=0.85$). Whether $\rho(\mathcal R_m)^2\to R$ at all is not proved here; it is the natural theoretical question left open by this section. (For the lower bounds, $\rho_W\to R$ is immediate from Lemma D.)
* **Block insertion / prime structure** (lower-bound side, for completeness): sewing closed sub-meanders into a truncated walk is not injective unless one restricts to prime blocks, and then it is the identity $M=P/(1-2P)$ of §14.4, whose pole is never reached. Wide excursions of real meanders are not self-contained blocks; the entropy is in structures that connect to the outside, which is exactly what neither a width cut nor a memory cut can summarise cheaply.

---

## 15. Window certificates: lower-bound vectors with infinite support (tried; a constant-factor gain only)

The truncation certificates of §§12–13 put weight zero on every state outside the chosen set $S$. Since the relaxation of §14 is, in the same language, an *upper*-bound vector supported on all states and depending only on the top of the stacks, it is natural to ask whether a *lower*-bound vector of the same kind — positive on all states — can do better than truncation. It can, but only by a constant factor, and the reason is instructive.

### 15.1 What an infinite-support certificate must satisfy

**Lemma M1 (reachability).** Every configuration of the exact automaton of width $w$ is reachable from $\varnothing$ in at most $w$ crossings, and all $\binom{2p+1}{p}$ non-crossing configurations with $p$ paths ($w=2p$) are reachable. *(Verified for $w\le10$ by breadth-first search — `src/prelim.py`; the general construction opens each cross pair with one (open,open) and each same-side pair with an (open,open) followed later by the transfer (open, close), processing nested blocks inside out.)*

**Proposition M2.** Let $v\ge0$ be a function on all exact states with $v_{j_0}>0$ and $(T^{\!\top}v)_t\ge c\,v_t$ for every state $t$ (sum over all exact predecessors). If moreover $v_i\le C\,x^{\,w(i)}$ for some $x<1/c$ and all $i$ (where $w(i)$ is the width), then $R\ge c^2$.

*Proof.* Iterating, $\sum_i v_i N_n(i\to j_0)\ge c^nv_{j_0}$. Each walk $i\to j_0$ of length $n$, prefixed by a fixed shortest walk $\varnothing\to i$ of length $\ell_i\le w(i)$ (Lemma M1) and followed by the accepting step, is a closed meander; for fixed $n$ and $\ell$ the map is injective on $\{i:\ell_i=\ell\}$ (the state at position $\ell$ determines $i$), and a meander of order $k$ arises from at most $2k$ pairs $(\ell,n)$. Hence $\sum_k 2k\,M_kx^{2k}\ \ge\ \sum_n x^{n+1}\sum_i x^{\ell_i}N_n(i\to j_0)\ \ge\ \sum_n x^{n+1}c^n v_{j_0}/C$, which diverges if $xc\ge1$; so $\sum_kM_kx^{2k}$ diverges for every $x\in(1/c,\,1/\sqrt R)$ if that interval is non-empty, a contradiction. $\square$

The decay condition is not a technicality. Walks of the exact automaton that do *not* return to $\varnothing$ may well grow faster than $\sqrt R$ per crossing (walks ending at width $\varepsilon n$ are bounded only by $R^{(1+\varepsilon)/2}$ per crossing), so a slowly decaying $v$ could certify a "growth" that is not the growth of meanders. Consequently any infinite-support ansatz $v=g(\text{window})\,\theta^{\#\text{deep arcs}}$ is forced to take $\theta<1/c$, i.e. $\theta<1/\sqrt R\approx0.29$, whereas the true Perron vector decays at roughly $0.67$ per arc.

### 15.2 The window automaton

Take a window of the top $k$ entries of each stack. The *pattern* of a state is the window with each entry labelled if its partner is inside the window and marked $?$ otherwise, plus one flag per side recording whether arcs exist below the window. Patterns are finite in number ($21,64,212,728,2593$ for $k=2,\dots,6$, mirror-quotiented; the cut-line sequence of a pattern is a word in the atoms "$?$" and "chord enclosing a complete non-crossing matching", so their number grows like a central binomial coefficient rather than a Catalan number squared). The ansatz $v(s)=g(\pi(s))\,\theta^{d(s)}$ is verified pattern by pattern: for each pattern $t$ and each reverse event, the predecessor's pattern is determined up to the labels of entries that *enter* the window from below, and the certificate takes the **minimum of $g$ over all completions** consistent with non-crossing (at most one candidate partner per side plus "$?$"). Predecessors that cannot be certified from the window alone — a (close,close) that split a path with a deep end — are simply not counted. At $\theta=0$ the scheme is exactly the per-side height truncation of §7.

The best $g$ is a max–min problem (a concave, monotone, homogeneous operator $\Phi$; $\sup_g\min_t\Phi g(t)/g(t)$). Plain power iteration converges to a spurious fixed point with zeros; the correct value is found by bisection on $c$ with the decreasing iteration $g\leftarrow\min(g,\Phi g/c)$ from $g\equiv1$, whose limit is the largest feasible $g$ (or $0$). Results (`src/window.py`, iteration not fully converged, so "certified" is a lower bound of the optimum):

| $k$ | patterns | $\theta=0$ (= $\mu_{k-1}$) | best $\theta<1/c$ | gain |
|---|---|---|---|---|
| 2 | 21 | 1.4142 | 1.802 ($\theta=0.29$) | +0.39 |
| 3 | 64 | 2.0529 | 2.235 | +0.18 |
| 4 | 212 | 2.4488 | 2.553 | +0.10 |
| 5 | 728 | 2.6997 | 2.777 | +0.08 |
| 6 | 2593 | 2.8739 | 2.929 ($\theta=0.30$) | +0.06 |

(per crossing). The gain is real but shrinks with $k$, settling at about half of one step of the truncation family. The cause is the adversarial minimum: whenever a forgotten arc re-enters the window, the certificate must assume the least favourable pairing, and the value of information about deep arcs, taken in the worst case, is small. This is the exact mirror image of §14.5's finding for upper bounds (where the *maximum* over completions is what makes forgetting cheap). A total-window variant, or the same ansatz layered on the adaptive set of §13, would add the same constant fraction of a step (an estimated $+0.04$ on $11.5975$); we did not build it because it would not change the convergence law, which is what would have justified the effort. Recorded so that the natural question "why not weight the states you throw away?" has an answer.

---

## 16. Do the forgetting bounds converge to $R$? (partial results and an open problem)

Write $\rho_m$ for the Perron root (per crossing) of the relaxed automaton $\mathcal R_m$ of §14 with total memory $m$, so $R\le\rho_m^2$.

**Proposition N1.** $\rho_m$ is non-increasing in $m$, hence $\rho_\infty=\lim_m\rho_m$ exists and $\rho_\infty\ge\sqrt R$.

*Proof.* The language accepted by $\mathcal R_{m'}$ for $m'>m$ is contained in that of $\mathcal R_m$: the reduced state of $\mathcal R_m$ covers that of $\mathcal R_{m'}$ (it is obtained from it by further forgetting), so by the simulation argument of Lemma L1 every word $\mathcal R_{m'}$ accepts, $\mathcal R_m$ accepts. Growth rates of nested languages are ordered. $\square$

**Proposition N2 (what the relaxed language is).** A word (a pair of arch systems, i.e. a system of closed curves) is accepted by $\mathcal R_m$ iff none of its loops is *visible*: a loop is visible when, at the crossing that closes it, both closed arcs carry the same label, and this can only fail if some arc of the loop was, at some earlier time, among the entries deleted by the memory rule — which requires the total number of open arcs to have exceeded $m$ at that time. In particular every loop of an accepted system contains an arc that was open at a moment when the width exceeded $m$ ("$m$-covered systems"). Meanders are the systems with one loop, and their single loop closes only at the last crossing, so they are accepted for every $m$.

**What a convergence proof would need.** $\rho_\infty=\sqrt R$ is equivalent to: for every $\varepsilon>0$ there is $m$ with $\#\{m\text{-covered systems of length }n\}\le(\sqrt R+\varepsilon)^n$ for large $n$. Two natural strategies both stall:

* *Surgery.* Reconnecting the $L$ loops of a covered system into one loop by $L-1$ local crosswise swaps, and recording the swaps, injects covered systems into meanders times $n^{O(L)}$; this proves convergence if $L=o(n/\log n)$ for all accepted systems, but it is not so: one excursion to width $m+1$ ("a cover") makes *every* older arc forgettable at once, so a single cover of cost $O(m)$ can hide $\Theta(n)$ loops.
* *Dilution.* In the cover construction the hidden loops live among arcs that were opened during an *exact* (loop-free) prefix and are closed during a phase whose entropy per crossing is that of the relaxed dynamics with a finite supply of forgotten arcs. Whether such mixtures beat $\sqrt R$ depends on two quantities we cannot control: the large-deviation cost $I(\varepsilon)$ of an exact prefix ending at width $\varepsilon n$, and the entropy bonus per forgotten arc in the closing phase. If $I(\varepsilon)=o(\varepsilon)$ — as it would be if the width of a random meander prefix fluctuated like a random walk, $I(\varepsilon)\asymp\varepsilon^2$ — and the bonus per arc is bounded below, then $\rho_\infty>\sqrt R$ and the family does **not** converge; if $I$ is linear with a large enough slope, it does. Neither is known, and proving $\rho_\infty>\sqrt R$ would in any case require exhibiting a sub-language with growth provably above $R$, which is impossible while $R$ itself is only known to lie in $[11.5975,\,12.6319]$.

**Numerics.** $\rho_m^2=12.7425,\ 12.7120,\ 12.6869,\ 12.6656,\ 12.6475,\ 12.6319$ for $m=16,\dots,26$; the successive decrements $0.0305,0.0251,0.0213,0.0181,0.0156$ fall like $m^{-1.65}$ (log–log slope over the range), i.e. $\rho_m^2-\rho_\infty^2\sim m^{-0.65}$, and the tail sum from $m=26$ is then $\approx0.28$, giving $\rho_\infty^2\approx12.35$ with an uncertainty of at least $\pm0.1$ that comfortably contains $R\approx12.26$. A geometric fit of the same decrements (ratio $0.86$) would give $12.53$. The data cannot decide between convergence and a limit slightly above $R$; the decay is clearly slower than $1/m$, which rules out the naive "one hidden loop per $m$ crossings" heuristic of §14.5.

**Open problem.** Determine $\rho_\infty$; in particular decide whether $\lim_m\rho(\mathcal R_m)^2=R$. A positive answer would make the forgetting relaxation a *complete* method (bounds converging to $R$ from above), complementing the truncations, which converge from below (Lemma D). We state this in the paper as a question, with Propositions N1–N2 and the numerical evidence, and we do not claim convergence.
