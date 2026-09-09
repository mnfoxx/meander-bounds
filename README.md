# Transfer-matrix bounds for meanders and stamp foldings

Code, certificates and notes for the paper

> M. Fox, *Transfer-matrix bounds for meanders and stamp foldings: 11.5975 ≤ R ≤ 12.6319 and R̄ ≤ 3.5542* (2026).

`R` is the growth constant of closed meanders (OEIS A005315) and `R̄` that of semi-meanders / stamp foldings (A000682, A000136). The previous rigorous bounds were 11.380 ≤ R ≤ 12.901 (Albert–Paterson 2005) and 3.065 ≤ R̄ ≤ 4 (Uehara 2011).

## Layout

| path | contents |
|---|---|
| `paper/` | the manuscript (`meanders.tex`, `meanders.pdf`) |
| `notes/` | the full lab journal (`depth1_semimeanders_proof.md`, §§1–16) and the project handoff |
| `src/` | all programs (C and Python; see below) |
| `results/` | run logs, certificate summaries, exact generating functions |

The lower-bound certificate itself (`s58_states.bin`, 58,398,791 state records of 20 bytes, and `s58_vec.bin`, the integer vector) is 1.6 GB and is distributed as a release asset, not in the tree.

## Reproducing the two headline bounds

Everything needs only `gcc`; the Python scripts need NumPy and SciPy. A 2-core machine with 7 GB of RAM suffices.

**Lower bound, R ≥ 11.59753** (adaptive selection, notes §13, paper §3):

    gcc -O2 -o grow src/grow.c -lm
    ./grow 20 32 150 run -1 -1 1000000 -1 -1 2000000 -1 -1 3500000 -1 -1 3500000 -1 -1 3500000 \
                        -1 -1 7000000 -1 -1 14000000 -1 -1 14000000 -1 30000000 40000000 40000000
    # ~3 h; each stage prints a certified bound. Re-iterate the best stage to convergence:
    ./grow -1 32 400 run            # resumes from run_ck_*.bin, writes run_states.bin / run_vec.bin

**Independent check of a lower-bound certificate** (different state model, exact integers):

    gcc -O2 -o check src/check.c
    ./check s58                      # reads s58_states.bin, s58_vec.bin; prints "CERTIFICATE OK: R >= 11.597535632"

**Upper bounds, R ≤ 12.63185 and R̄ ≤ 3.55414** (forgetting relaxation, notes §14, paper §4):

    gcc -O2 -o relax src/relax.c -lm
    ./relax 40 600 26                # total memory 26: 72.8M states, ~6 GB, ~1.5 h
    ./relax 7 300                    # per-side memory 7: 4,301 states, seconds, already beats 12.901
    python3 src/relax_py.py 40 16    # independent Python build of the same automaton
    python3 src/relax_check.py 12 4  # mechanical check of the soundness lemma to width 12

**Classical width truncation, R ≥ 11.54507** (notes §12): `gcc -O2 -o muq src/muq.c -lm; ./muq 28 170 v28.bin` repeated ~5 times (2.5 GB, checkpoints in `v28.bin`).

## Program index

| file | purpose |
|---|---|
| `grow.c` | adaptive (Perron-mass) state selection with exact integer certificates; checkpoints |
| `check.c` | independent certificate checker (cut-line/partner model, BFS reachability) |
| `relax.c`, `relax_py.py`, `relax_check.py`, `relax_adapt.c` | forgetting relaxation for upper bounds; Python twin; soundness check; adaptive-coarsening experiment |
| `muq.c`, `mucert.c`, `muw.c`, `mubig.c`, `mass.c` | width/height truncations, mirror quotient, Perron-mass diagnostics |
| `adapt.c` | pilot: top-K-by-mass subsets of a fixed-width automaton |
| `window.py`, `prelim.py` | window certificates (notes §15) and the reachability lemma |
| `tm.py`, `tmD.py`, `bounded_meanders.py`, `openmeanders.py`, `phases.py`, `tailgf.py`, `meanderphase.py`, `mu.py`, `mu2.py`, `lam.py` | exact generating functions at bounded height; phase decomposition |
| `fold*.c`, `list1.c`, `wbrute.py`, `wcount.py`, `wgf2.py`, `closedcheck.py`, `bcheck.c` | brute-force enumerators used as cross-checks |
| `prime.py`, `irred.py` | prime-meander decomposition M = P/(1−2P) |
| `escape.c` | winding / escape-number distributions (DFGG territory) |
| `tilt.py`, `relaxmass.py`, `guess*.py` | experiments recorded in notes §14.5 |

## Verification chain

- Automata validated against brute-force enumeration of pairs of non-crossing matchings (closed meanders, widths ≤ 8) and of foldings (semi-meanders, n ≤ 26).
- Every certificate is an exact integer inequality (`(Tᵀ)²v ≥ c·v` for lower bounds, `T²w ≤ c·w` for upper bounds); floating point only produces the candidate vector.
- The lower-bound checker (`check.c`) shares no code with the producer (`grow.c`) and uses a different representation of states.
- The upper-bound automaton was built independently in C and Python with identical state counts and eigenvalues.

## Licence

Code: MIT. Text of the notes and paper: CC BY 4.0.
