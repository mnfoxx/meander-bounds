# Results so far (all exact / certified unless marked)

| depth D | states | recurrence order | growth λ_D | tail-factor degree in x² | pre-tail degree |
|---|---|---|---|---|---|
| 1 | 18  | 5   | φ ≈ 1.6180 | 1  | 2  |
| 2 | 64  | 15  | 2.1337     | 3  | 8  |
| 3 | 236 | 37  | 2.4833     | 7  | 22 |
| 4 | 836 | 109 | 2.7219     | 16 | 76 |
| 5 | 2976 | ? | 2.8867 | ? | ? |

Tail constants μ_D (rootless/meander automaton): 1.4142, 2.0529, 2.4488, 2.6997, 2.8739, 2.9960 (D=1..6).
Theorem C (notes §8): μ_{D−1} ≤ λ_D ≤ μ_{D+1}, hence lim λ_D = lim μ_D.
Theorem E (notes §9): lim λ_D = lim μ_D = √R exactly (meander constant), via
  b_n := #{1 on top, n on bottom} = open meanders(n−2) = A005316 (Lemma F).
μ_7 = 3.0864, μ_8 = 3.1542.
DFGG conjecture R̄ = √R  ⟺  b_n/s(n) decays sub-exponentially. Data: b_n/s(n) ~ n^−1.44 (n≤18).

s_D(n) = number of semi-meanders (stamp 1 on top) of n stamps with nesting depth ≤ D.
Labelled foldings with depth ≤ D = n · s_D(n).

Closed form, D=1 (proved in notes/): s_1(n) = 2F(n+2) − 3·2^(n/2) (n even), 2F(n+2) − 2^((n+3)/2) (n odd).

Tail-phase growth vs pre-tail growth: D=2: 2.053 / 2.134; D=3: 2.449 / 2.483; D=4: 2.700 / 2.722.
Gap shrinking: 0.08, 0.03, 0.02.

Other verified facts (brute force, n ≤ 18):
- Position of stamp 1 in the pile is uniform (⇒ A000136 = n·A000682). Known.
- Semi-meanders with stamp 2j at bottom = with stamp 2j+1 at bottom, all j. Unproved, likely easy involution.
- Semi-meanders with exactly n−k hairpins: k=2: n−2; k=3: (n−3)(3n−10)/2; k=4: cubic quasi-polynomial, period 2.
  Conjecture: degree k−1 quasi-polynomial for all k.

Failed guess: tail degrees are NOT 2^D − 1 (1,3,7,16).

Escape numbers (notes §10): R̄ = √R ⟸ E[e(F_n)] = o(n). Mean e: 2.52 (n=8), 2.88 (10), 3.21 (12), 3.34 (13); local exponent ~0.6 falling.

## Meander constant lower bound (certified)
W=20: 11.04139; W=22: 11.20795; W=24: 11.34221; W=26: 11.45295; W=28: R >= 11.54507 (Albert–Paterson 2005: 11.380).
Adaptive selection (notes §13): 3.5M states 11.40008; 30M 11.54774; 40M 11.58363; 58.4M **R >= 11.59753** (Theorem K).
Upper bounds (notes §14): per-side k=7: 12.8899; total memory m=20: 12.6869; m=22: 12.6656; m=24: **R <= 12.63185**, hence **R̄ <= 3.55413**.
Two-sided: 11.5975 <= R <= 12.6319 ; 3.4055 <= R̄ <= 3.5542.
