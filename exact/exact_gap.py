#!/usr/bin/env python3
"""Exact gap probabilities for the highest c-branches of a random mapping.

For a uniformly random mapping g of {1..n}, a c-branch is the subtree rooted at a
vertex at distance c from the cycles (0-branches are the trees).  gap_c(g) is the
difference between the heights of the highest and the second highest c-branch
(0 if the maximum is attained twice or there is no c-branch, infinity if there
is exactly one c-branch).

This script evaluates P(gap_c(g) >= j) from the exact generating function

    G_j^{(c)}(z) = sum_{H >= 0} Delta_H(z) * psi_c(max(H - j, -1)),
    psi_c(k)     = t_{k+1} ... t_{k+c} / (1 - t_{k+c})^2,

where t_h is the EGF of rooted labelled trees of height <= h (t_{-1} = 0,
t_h = z exp(t_{h-1})) and Delta_H = t_H - t_{H-1}.  All series are truncated at
degree N and computed in the scaled variable w = e*z, so that the coefficients
are of moderate size and double precision is sufficient.

The sum over H is stopped once the probability that a tree with N vertices has
height > H is below `tol`; `--full` disables this and sums up to H = N.
"""
import argparse
import csv
import math
import sys

import numpy as np


class Truncated:
    """Arithmetic of power series truncated at degree N (FFT multiplication)."""

    def __init__(self, N):
        self.N = N
        self.size = 1 << (2 * N + 1).bit_length()

    def mul(self, a, b):
        fa = np.fft.rfft(a, self.size)
        fb = np.fft.rfft(b, self.size)
        return np.fft.irfft(fa * fb, self.size)[: self.N + 1]

    def exp_small(self, d):
        """exp(d) for a series d with non-negative coefficients (Taylor series in d)."""
        norm = float(np.abs(d).sum())
        out = np.zeros(self.N + 1)
        out[0] = 1.0
        term, bound, m = out.copy(), 1.0, 0
        while True:
            m += 1
            term = self.mul(term, d) / m
            bound *= norm / m
            out += term
            if bound < 1e-19:
                return out

    def refine_inverse(self, a, q):
        """Newton iteration for 1/a, starting from the approximation q."""
        for _ in range(60):
            resid = -self.mul(a, q)
            resid[0] += 1.0
            if np.abs(resid).max() < 1e-15:
                return q
            q = q + self.mul(q, resid)
        raise RuntimeError("inverse did not converge")


def gap_probabilities(ns, cs=(0, 1, 2), jmax=3, tol=1e-13, full=False, progress=False):
    """Return {(n, c, j): P(gap_c >= j)} for n in ns, c in cs, 1 <= j <= jmax."""
    N, cmax = max(ns), max(cs)
    S = Truncated(N)
    n_arr = np.arange(1, N + 1)
    lg = np.array([math.lgamma(k + 1) for k in n_arr])
    tree = np.zeros(N + 1)                       # scaled coefficients of t(z)
    tree[1:] = np.exp((n_arr - 1) * np.log(n_arr) - n_arr - lg)
    total = np.zeros(N + 1)                      # scaled coefficients of M(z) = 1/(1-t)
    total[1:] = np.exp(n_arr * np.log(n_arr) - n_arr - lg)

    one = np.zeros(N + 1)
    one[0] = 1.0
    w = np.zeros(N + 1)
    w[1] = math.exp(-1.0)                        # z = w/e

    keep = cmax + jmax + 3
    t = {-1: np.zeros(N + 1)}                    # t_h
    delta = {}                                   # Delta_h
    q = {-1: one.copy()}                         # 1/(1 - t_h)
    psi = {c: {} for c in cs}                    # psi_c(k)
    G = {(c, j): np.zeros(N + 1) for c in cs for j in range(1, jmax + 1)}
    next_H = {key: 0 for key in G}
    expt = one.copy()                            # exp(t_{h-1})

    def make_psi(c, k):
        p = S.mul(q[k + c], q[k + c])
        for i in range(1, c + 1):
            p = S.mul(p, t[k + i])
        return p

    for c in cs:
        if c == 0:
            psi[c][-1] = one.copy()

    h, last = 0, N + cmax + 1
    while h <= last:
        t[h] = S.mul(w, expt)
        delta[h] = t[h] - t[h - 1]
        expt = S.mul(expt, S.exp_small(np.maximum(delta[h], 0.0)))
        q[h] = S.refine_inverse(one - t[h], q[h - 1])
        for c in cs:
            k = h - c
            if k >= -1 and k not in psi[c]:
                psi[c][k] = make_psi(c, k)
        for (c, j), acc in G.items():
            while True:
                H = next_H[(c, j)]
                k = max(H - j, -1)
                if H > h or k not in psi[c]:
                    break
                acc += S.mul(delta[H], psi[c][k])
                next_H[(c, j)] = H + 1
        for store in [t, delta, q] + [psi[c] for c in cs]:
            for key in [x for x in store if x < h - keep and x != -1]:
                del store[key]
        tail = (tree[N] - t[h][N]) / tree[N]     # P(height > h | tree size N)
        if progress and h % 200 == 0:
            print(f"  height {h:6d}   P(height > h | size N) = {tail:.3e}", file=sys.stderr)
        if not full and tail < tol and last > h + keep:
            last = h + keep                      # flush the pending terms and stop
        h += 1

    return {(n, c, j): G[(c, j)][n] / total[n] for n in ns for (c, j) in G}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("n", type=int, nargs="+", help="mapping sizes")
    ap.add_argument("--cmax", type=int, default=2)
    ap.add_argument("--jmax", type=int, default=3)
    ap.add_argument("--tol", type=float, default=1e-13)
    ap.add_argument("--full", action="store_true", help="sum over all heights up to n")
    ap.add_argument("-o", "--output", default="-")
    args = ap.parse_args()
    res = gap_probabilities(sorted(args.n), tuple(range(args.cmax + 1)), args.jmax,
                            args.tol, args.full, progress=True)
    out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
    wr = csv.writer(out, lineterminator="\n")
    wr.writerow(["n", "c", "j", "prob_gap_lt_j", "scaled_by_sqrt_n", "limit"])
    for (n, c, j), p in sorted(res.items()):
        fail = 1.0 - p
        wr.writerow([n, c, j, f"{fail:.12e}", f"{fail * math.sqrt(n):.8f}",
                     f"{(2 * j - 1) * math.sqrt(math.pi / 8):.8f}"])


if __name__ == "__main__":
    main()
