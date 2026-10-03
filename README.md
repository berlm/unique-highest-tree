# Highest trees of random mappings: exact values and simulations

Code and data for Section 5 of the paper *Highest Trees of Random Mappings*
(arXiv:1504.04532). Repository: <https://github.com/berlm/unique-highest-tree>.

A random mapping of `n` elements is viewed as a functional graph. The *depth* of
a vertex is its distance to the cycles. For `c >= 0`, a *c-branch* is the subtree
rooted at a vertex of depth `c` (0-branches are the trees), and its *height* is
the largest distance from its root to one of its vertices.

* `gap_c` is the difference between the heights of the highest and the second
  highest c-branch. It is 0 if the maximal height is attained at least twice or
  there is no c-branch, and infinite if there is exactly one c-branch.
* If `gap_c >= 1`, the *crown* `H` consists of the vertices of the highest
  c-branch that lie above the height of the second highest one, and `r` is the
  number of its roots. The crown is *sparse* for `alpha` if `|H| <= alpha * r`.

The paper proves that `P(gap_c < j) ~ (2j - 1) * sqrt(pi/8) / sqrt(n)` and that
`|H| > alpha * r > 0` holds with probability `1 - Theta(n^(-1/2))`.

## Contents

| Path | Purpose |
|---|---|
| `exact/exact_gap.py` | Exact values of `P(gap_c < j)` from the generating functions (power series arithmetic with NumPy). |
| `sim/mapping_sim.c` | Monte Carlo simulation of gap and crown statistics; also enumerates all `n^n` mappings for small `n`. |
| `scripts/make_tables.py` | Builds the tables (`results/tables.md`, `results/tables.tex`) from the raw results. |
| `scripts/run_all.sh` | Reproduces everything in `results/`. |
| `tests/test_small.py` | Checks both programs against an independent enumeration for `n = 3, 4, 5`. |
| `results/` | Raw results and tables used in the paper. |

## Requirements

Python 3 with NumPy, and a C99 compiler with `make`. Nothing else.

## Usage

```sh
make -C sim                      # build the simulator
python3 tests/test_small.py      # consistency tests (a few seconds)
sh scripts/run_all.sh            # reproduce results/ (about 11 minutes on one core)
```

Individual programs:

```sh
# exact P(gap_c < j) for c = 0..2, j = 1..3
python3 exact/exact_gap.py 100 1000 10000 -o out.csv

# simulation: n, trials, seed, cmax, alpha values
./sim/mapping_sim sim 1000 100000 12345 2 1.5 2 3

# exhaustive enumeration of all n^n mappings (n <= 9)
./sim/mapping_sim enum 5 2 2
```

The simulator prints one CSV line per `c` and `alpha`:

```
n,c,alpha,trials,seed,gap0,gap1,gap2,single,sparse
```

`gap0`, `gap1`, `gap2` count the mappings with `gap_c` equal to 0, 1, 2;
`single` counts those with exactly one c-branch; `sparse` counts those with
`gap_c >= 1` and a sparse crown. The number of mappings violating
`|H| > alpha * r > 0` is `gap0 + sparse`.

## Notes on the methods

**Exact values.** The exponential generating function of the mappings with
`gap_c >= j` is `sum_H Delta_H * psi_c(max(H - j, -1))`, where `t_h` is the
generating function of trees of height at most `h`, `Delta_H = t_H - t_{H-1}`
and `psi_c(k) = t_{k+1} ... t_{k+c} / (1 - t_{k+c})^2`. The series are truncated
at degree `n` and computed in the variable `e*z`, in double precision. The sum
over `H` stops when the probability that a tree with `n` vertices is higher
than `H` is below `--tol` (default `1e-13`); `--full` sums over all heights.

**Simulation.** Mappings are generated with xoshiro256** (seeded through
splitmix64) and unbiased bounded integers. Each mapping is processed in linear
time by repeatedly removing leaves (Kahn's algorithm). Runs are reproducible
from the seeds in `scripts/run_all.sh`.

## Results

See `results/tables.md`. In short: the exact tie probabilities for
`c = 0, 1, 2` behave like `sqrt(pi/8)/sqrt(n) + (c - 1)/(2n)`; the simulation
agrees with the exact gap probabilities (largest deviation 1.76 standard
deviations in 54 comparisons); and the scaled crown failure probabilities level
off as `n` grows.

## License

MIT; see `LICENSE`.
