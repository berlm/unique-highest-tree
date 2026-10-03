#!/usr/bin/env python3
"""Consistency tests on small mappings.

1. A direct (slow, independent) enumeration of all n^n mappings in Python is
   compared with the exact generating-function computation (exact/exact_gap.py).
2. The same enumeration is compared with the C program in enumeration mode
   (sim/mapping_sim enum), including the crown statistics.

Run from the repository root:  python3 tests/test_small.py   (or: pytest tests)
"""
import itertools
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "exact"))
from exact_gap import gap_probabilities  # noqa: E402

BINARY = os.path.join(os.path.dirname(__file__), "..", "sim", "mapping_sim")
INF = 10**9


def statistics(f, c):
    """Return (gap_c, crown size, crown roots) of the mapping f (a tuple)."""
    n = len(f)
    cyclic = set()
    for start in range(n):
        x = start
        for _ in range(n):
            x = f[x]
        cyclic.add(x)

    def depth(x):
        d = 0
        while x not in cyclic:
            x, d = f[x], d + 1
        return d

    def root_of(x):  # the ancestor of x at depth c
        while depth(x) > c:
            x = f[x]
        return x

    levels = {}  # branch root -> list of distances of its vertices from the root
    for x in range(n):
        if depth(x) >= c:
            levels.setdefault(root_of(x), []).append(depth(x) - c)
    heights = sorted((max(v), r) for r, v in levels.items())
    if not heights or (len(heights) > 1 and heights[-1][0] == heights[-2][0]):
        return 0, 0, 0
    h2 = heights[-2][0] if len(heights) > 1 else -1
    gap = heights[-1][0] - h2 if len(heights) > 1 else INF
    top = levels[heights[-1][1]]
    return gap, sum(1 for d in top if d >= h2 + 1), sum(1 for d in top if d == h2 + 1)


def enumerate_counts(n, cmax, alphas):
    counts = {(c, a): dict(gap0=0, gap1=0, gap2=0, single=0, sparse=0) for c in range(cmax + 1) for a in alphas}
    for f in itertools.product(range(n), repeat=n):
        for c in range(cmax + 1):
            gap, size, roots = statistics(f, c)
            for a in alphas:
                rec = counts[(c, a)]
                if gap == 0:
                    rec["gap0"] += 1
                    continue
                if gap == INF:
                    rec["single"] += 1
                elif gap <= 2:
                    rec[f"gap{gap}"] += 1
                if size <= a * roots:
                    rec["sparse"] += 1
    return counts


def test_exact_series_against_enumeration():
    for n in (3, 4, 5):
        counts = enumerate_counts(n, 2, (2.0,))
        exact = gap_probabilities([n], cs=(0, 1, 2), jmax=3, full=True)
        for c in range(3):
            rec, total = counts[(c, 2.0)], n**n
            below = {1: rec["gap0"], 2: rec["gap0"] + rec["gap1"], 3: rec["gap0"] + rec["gap1"] + rec["gap2"]}
            for j in (1, 2, 3):
                assert math.isclose(1 - exact[(n, c, j)], below[j] / total, abs_tol=1e-12), (n, c, j)


def test_c_program_against_enumeration():
    if not os.path.exists(BINARY):
        print("sim/mapping_sim not built; skipping the comparison with the C program")
        return
    alphas = (1.5, 2.0, 3.0)
    for n in (3, 4, 5):
        counts = enumerate_counts(n, 2, alphas)
        out = subprocess.run([BINARY, "enum", str(n), "2"] + [str(a) for a in alphas],
                             capture_output=True, text=True, check=True).stdout
        for line in out.strip().splitlines():
            _, c, a, trials, _, gap0, gap1, gap2, single, sparse = line.split(",")
            rec = counts[(int(c), float(a))]
            assert int(trials) == n**n
            assert (int(gap0), int(gap1), int(gap2), int(single), int(sparse)) == \
                (rec["gap0"], rec["gap1"], rec["gap2"], rec["single"], rec["sparse"]), line


if __name__ == "__main__":
    test_exact_series_against_enumeration()
    test_c_program_against_enumeration()
    print("all tests passed")
