# HANDOFF — stamp folding / meanders project (as of 9 Sept 2026, evening)

Upload this whole folder (the zip) into the new conversation. Everything below is self-contained;
`notes/depth1_semimeanders_proof.md` is the full lab journal (read its LITERATURE STATUS block first),
`results/` has the certificates, `src/` has all code.

## 0. One-paragraph status
We set out to "make progress on stamp folding". The problem itself (closed form for A000136 / A000682)
remains open and nothing here touches that. What we have is (a) one certified improvement of a
20-year-old rigorous bound on the meander growth constant, (b) a family of exact enumerative results
for meanders/semi-meanders of bounded height or width, and (c) some structural theorems. A large
middle chapter of the journal (escape numbers, winding, "R̄ = √R iff ν < 1") turned out to be
Di Francesco–Golinelli–Guitter 1996 verbatim and must be cited, not claimed.

## 0b. NEWEST (notes §14, Theorem L): UPPER BOUNDS. R ≤ 12.63185 (was 12.901, Albert–Paterson 2005) and
R̄ ≤ 3.55413 — to our knowledge the FIRST rigorous upper bound on the semi-meander / stamp-folding constant
(conjectured 3.50; Legendre's 2013 survey lists none). Method: forget path labels beyond a memory budget m
(deepest arcs first; singleton labels ≡ '?'); the relaxed automaton accepts every meander and every folding,
so its Perron root (certified by Collatz–Wielandt with an exact integer vector) bounds both constants from
above. Soundness = Lemma L1 (simulation/covering argument), also checked mechanically to width 12.
Code src/relax.c (memory ≤ 28 fits 7 GB up to m=26), independent Python build src/relax_py.py, simulation check
src/relax_check.py. m=7 per side (4301 states, seconds) already beats 12.901. Increments shrink ×0.85 per two
units of memory → this family stalls near 12.5; keeping partial information about deep arcs is the next idea.
Also §14.4: prime-meander identity M = P/(1−2P) (sewing has unique factorisation) — small remark, likely known.
Current two-sided picture: 11.6080 ≤ R ≤ 12.6319 (est. 12.2629); 3.4070 ≤ R̄ ≤ 3.5542 (est. 3.5019).
(10 Sep 2026: headline lower bound is now the author's laptop certificate, R ≥ 11.60805 on 55,781,199 states; see notes §13.3 addendum. The cloud certificate R ≥ 11.59753 on 58,398,791 states remains as an independent confirmation.)

## 1. The headline theorem (NEW, publishable)
**UPDATE (notes §13, Theorem K): R ≥ 11.60805** (laptop, 55,781,199 states; earlier cloud certificate 11.59753 on 58,398,791
orbit-states (widths ≤ 32), computed on a 2-core / 7 GB machine. Method: select states by Perron mass
a_i·b_i (first-order sensitivity of ρ to deleting the state), score frontier states by the same product
estimated from their in-set neighbours (CIPSI-style "selected configuration interaction"), expand/prune in
cycles. Worth ≈ one width step per state: 14M adaptive states ≈ 11.519 vs 11.453 for width ≤ 26 (13.5M).
Code src/grow.c; independent checker src/check.c (different state model; verified). Certificate files
s58_states.bin / s58_vec.bin (1.6 GB) delivered separately. Lemma K1 (notes §13.1) is the two-line
extension of the §12 argument to arbitrary state sets: only states reachable from j0 inside S count.
Everything below about W=28 is still true, it is now the second-best certificate.

Original headline: **R ≥ 11.54507**, where R = lim M_k^{1/k}, M_k = number of closed meanders of order k (2k crossings).
Previous best rigorous lower bound: 11.380 (Albert–Paterson, JCTA 2005), still cited as best known in
Delecroix–Goujard–Zograf–Zorich, Ann. Inst. Fourier 2024. Jensen's numerical estimate: 12.2629.
Upper bound 12.901 (Albert–Paterson) untouched.

How: Jensen-style transfer matrix for closed meanders, truncated to WIDTH ≤ W (total open arcs ≤ W),
quotiented by the mirror symmetry (swap the two sides). Any positive integer vector y with
(Eᵀ)² y ≥ c·y entrywise certifies R ≥ c (argument in notes §12: reachability + Fekete). y from power
iteration, rounded DOWN at scale 2^52; products and min-ratio in exact int64/int128 arithmetic.
Certified values: W=20: 11.04139; 22: 11.20795; 24: 11.34221; 26: 11.45295; 28: 11.54507
(rational 4367716392/378318673; W=28 has 52,327,345 orbit-states ≈ 1.05e8 states).
Automaton validated against brute force over all pairs of non-crossing matchings (widths 4,6,8; k≤6).
Quotient validated: reproduces full-automaton certificates at W=20, 24 to the last digit.
Corollary: semi-meander constant R̄ ≥ √R ≥ 3.3978 (easy direction M_k ≤ M̄_{2k}).

Reproduce: `gcc -O2 -o muq src/muq.c -lm; ./muq 28 170 v28.bin` repeated ~5 times (checkpoints the
eigenvector in v28.bin; each leg ~4.5 min, needs ~2.5 GB). W=30 needs ~12 GB, W=32 ~45 GB;
expected certified ≈ 11.62 and 11.68 (increments shrink by ≈0.83 per step).

## 2. Exact enumerative results (NEW)
All generating functions rational; all verified by independent brute force.
Semi-meanders = foldings of n stamps with stamp 1 on top (A000682). "Depth ≤ D" = at most D+1 arcs
stacked on either side of the pile (= height ≤ D+1). y = x².
- Semi-meanders, height ≤ 2: s(n) = 2F(n+2) − 3·2^{n/2} (n even), 2F(n+2) − 2^{(n+3)/2} (n odd), n≥3.
  GF x²(1+x−x²+2x³+2x⁴)/((1−x−x²)(1−2x²)). Hand proof: notes §3 (three-state automaton).
- Semi-meanders height ≤3,4,5: order 15/37/109 recurrences, denominators factored in notes §6.
- Closed meanders height ≤2: 2^{k−1}. Height ≤3: 1,2,8,34,144,608,2564,…; y(1−2y)²/(1−6y+8y²−2y³).
  Height ≤4: explicit (notes §11).
- Open meanders height ≤2: Fibonacci F_m. Height ≤3: 1,1,2,3,8,14,35,65,155,300,694,1382,…; GF in §11.
- Closed meanders width ≤4: 1,2,8,38,184,890,4296; y(1−3y)²/(1−8y+17y²−8y³). Width ≤6: §12 supp.
- The three height-2 formulas come from ONE five-state automaton (present them together).

## 3. Structural theorems (NEW but elementary)
- Rationality for every fixed height/width bound (finite automaton; notes §2).
- Phase decomposition (notes §7): the semi-meander automaton splits at the position of the last stamp;
  pre-phase = rooted automaton (its diagonal = open meanders), post-phase read bottom-up = rootless
  closed-meander automaton. Denominators factor accordingly; tail factor is a polynomial in x² because
  the post-phase automaton is bipartite (right-stack height changes by ±1 each step).
- Sandwich: μ_D ≤ λ_D ≤ μ_{D+1} (λ = semi-meander bounded-height constant, μ = closed-meander one).
  Lower half: cut an outermost arc of a closed meander (height doesn't increase). Upper half: freeze a
  partner arc at the bottom of a stack (notes §8, Theorem C; needs the μ_D ≤ λ_D update written in).
- Consequence lim λ_D = lim μ_D = √R — TRUE but essentially known (bounded height ⇒ bounded winding
  ⇒ meander growth, DFGG language). State as a remark with citation.
- States with p open paths in the closed-meander automaton number C(2p−1,p) (probably known).

## 4. KNOWN — cite, do not claim (see LITERATURE STATUS block at top of notes)
- DFGG, hep-th/9607039 (Nucl. Phys. B 1997?) "Meanders: a direct enumeration approach": winding w
  defined exactly as our escape number (their §2.1, Fig. 4); "R(q) = R̄(q) iff ν(q) < 1" (eq. 2.11);
  winding distribution to n=24; ν(0) = 0.52(1); scaling function. => our §10 is a rediscovery.
- DFGG hep-th/9602025 (TL algebra): R = R̄² conjecture as "same entropy per monomer".
- 1-on-top + n-on-bottom folding = open meander: folklore (combos.org; Sawada–Li EJC 2012 wind-factor).
- Our Lemma F parity remark was WRONG: odd n gives open meanders with odd crossings, not closed meanders
  of the same height (closing raises height by one). Counts coincide only because A005316(2k−1)=A005315(k).
- A000136 = n·A000682 (uniform position of stamp 1): known.
- Jensen 2000 (cond-mat/0008178): the transfer matrix itself (used for exact counts, not bounds).
- Albert–Paterson 2005: bounds via a word language (shifts/jumps), not width truncation.

## 5. Unproved observations (data only) — candidates for small lemmas or exclusion
- #{semi-meanders with exactly n−k hairpins}: n−2; (n−3)(3n−10)/2; cubic quasi-polynomial (period 2)
  for k=2,3,4. Conjecture: degree k−1 quasi-polynomial. (Flavour of DFGG's large-q polynomials.)
- Semi-meanders with stamp 2j at bottom = with stamp 2j+1 at bottom, all j (n ≤ 18).
- Winding/escape distribution: #{w=n−1}=#{w=n−2}=1; #{w=n−3}=⌊(n+1)²/4⌋−4 (6≤n≤20);
  E_n(1)=E_{n−1}(2) PROVED (notes §10 discussion: remove last stamp; unique gap). But DFGG territory.
  Mean w ≈ 1.07√n − 0.48 for 8≤n≤20 (ν = 1/2 numerically; DFGG say 0.52(1)).
- Extrapolation of width-truncated ρ_W gives R ≈ 12.4 drifting down (sanity only).

## 6. Paper plan
Title idea: "Bounded-width transfer matrices for meanders: an improved lower bound on the meander
constant, and exact generating functions at small height."
Title now: "Transfer-matrix bounds for meanders and stamp foldings: 11.6080 ≤ R ≤ 12.6319, R̄ ≤ 3.5542,
and exact generating functions at small height" (or similar). Sections: 1 intro & statement;
3c upper bounds by forgetting (§14: relaxed automaton, Lemma L1, table, semi-meander corollary); 2 the automaton and its verification; 3 the certificate
argument (Fekete + reachability + quotient lemma + Lemma K1 for arbitrary state sets) with the rational
and reproducibility; 3b adaptive state selection (§13: heuristic, schedule, table vs width cuts); 4 exact
generating functions (height 2 closed forms; tables for height ≤4, width ≤6); 5 structure (phase
decomposition, sandwich); 6 remarks & citations (DFGG winding; convergence to √R; R̄ ≥ 3.3978).
Model for presenting a certified computation: Barequet–Rote–Shalah, "λ > 4" (CACM), which certifies a
transfer-matrix bound for polyominoes.
Before submission: (i) someone with journal access re-checks 2005–2026 for any transfer-matrix
lower bound on R; (ii) run W=30 (and 32) on a big machine; (iii) consider posting code + certificates.

## 7. Open research directions still on the table
- DONE: adaptive (Perron-mass) truncation, notes §13 — one width step per state for free. Next: run
  src/grow.c on a 64 GB machine with K ≈ 3–4·10^8 (schedule in §13); expect R ≥ 11.66–11.68. Possible
  refinements: second-order (Epstein–Nesbet-like) scores, a two-step frontier, OpenMP on the gathers.
- Bigger W as a plain width cut is now dominated by the adaptive method.
- Upper bounds: a relaxation automaton (forget deep labels) is easy to write but likely gives ≈16;
  AP's 12.901 via cluster method looks hard to beat this way.
- Semi-meander-specific bounds: λ_W ≤ μ_{W+1} so nothing beyond √R this way.
- Prove ν < 1 (DFGG's condition ⇔ R̄ = √R): the real prize; nothing here attacks it.
- Recursion in the height bound for the exact GFs (à la Bousquet-Mélou F_k/F_{k+1} for Dyck paths)?
  Denominator degrees 1,3,7,16 (tail) and 2,8,22,76 (pre-tail) — no pattern found; not 2^D−1.
- Prove the small quasi-polynomial and involution observations of §5 if a home is found for them.

## 8. File map
src/fold*.c, list1.c      brute-force folding enumerators (fold5: n-at-bottom + depth distribution)
src/tm.py, tmD.py         depth-D semi-meander automaton → certified minimal recurrence + factored GF
src/phases.py, tailgf.py, meanderphase.py   phase decomposition experiments
src/mu.py, mu2.py, lam.py  μ_D, λ_D by eigenvalues (height truncation)
src/bounded_meanders.py, openmeanders.py, closedcheck.py, bcheck.c, wcount.py, wbrute.py, wgf2.py
                           bounded-height/width meander GFs and brute-force checks
src/mubig.c, muw.c        height-/width-truncated μ in C with float CW bound
src/mucert.c              width-truncated, EXACT integer certificate (full automaton)
src/muq.c                 mirror-quotient, packed states, checkpointing — produced the W=28 headline
src/mass.c                Perron-mass diagnostics by width/#paths
src/relax.c               label-forgetting relaxation → certified UPPER bounds (R ≤ 12.6319, R̄ ≤ 3.5542)
src/relax_py.py           independent Python build of the same automaton (agrees to all digits, m ≤ 16)
src/relax_check.py        mechanical simulation check of Lemma L1 (exact width ≤ 12)
src/prime.py, irred.py    prime-meander decomposition experiments (M = P/(1−2P))
results/relax/            run logs
src/adapt.c               pilot: top-K-by-mass subsets of a fixed-W automaton vs width baselines (§13)
src/grow.c                ADAPTIVE selection (expand/prune cycles, chunked frontier, checkpoints) — produced R ≥ 11.60805 (laptop) / 11.59753 (cloud)
src/check.c               independent certificate checker (cut-line/partner model, BFS reachability, exact integers)
results/adaptive/         stage logs of the runs, converged s58 log, checker output
src/escape.c, escape2.c   winding (escape number) distributions (DFGG territory)
results/certificate_W26.txt, certificate_W28.txt, depth4_gf.txt, summary.md, escape_distribution_n14-20.txt
notes/depth1_semimeanders_proof.md   the full journal, §§1–12 + literature ledger at top
