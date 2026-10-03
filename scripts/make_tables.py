#!/usr/bin/env python3
"""Build the tables of the paper from results/exact_gap.csv and results/simulation.csv.

Writes results/tables.md (for reading) and results/tables.tex (LaTeX tabulars).
Run from the repository root.
"""
import csv
import math
from collections import defaultdict

LIMIT = math.sqrt(math.pi / 8)
Z95 = 1.959963984540054


def wilson(k, m, z=Z95):
    """Wilson score interval for a binomial proportion k/m."""
    p = k / m
    denom = 1 + z * z / m
    centre = (p + z * z / (2 * m)) / denom
    half = z * math.sqrt(p * (1 - p) / m + z * z / (4 * m * m)) / denom
    return centre - half, centre + half


def read_exact(path="results/exact_gap.csv"):
    exact = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            exact[(int(row["n"]), int(row["c"]), int(row["j"]))] = float(row["prob_gap_lt_j"])
    return exact


def read_sim(path="results/simulation.csv"):
    """Sum the counts over all runs with the same (n, c, alpha)."""
    sim = defaultdict(lambda: defaultdict(int))
    with open(path) as fh:
        for row in csv.DictReader(fh):
            key = (int(row["n"]), int(row["c"]), float(row["alpha"]))
            for field in ("trials", "gap0", "gap1", "gap2", "single", "sparse"):
                sim[key][field] += int(row[field])
    return sim


def fmt_n(n):
    return f"{n:,}"


def main():
    exact, sim = read_exact(), read_sim()
    ns = sorted({n for (n, _, _) in exact})
    cs = sorted({c for (_, c, _) in exact})
    js = sorted({j for (_, _, j) in exact})
    alphas = sorted({a for (_, _, a) in sim})
    md, tex = [], []

    # ---- Table 1: exact tie probabilities --------------------------------
    md += ["## Table 1. Exact tie probabilities: sqrt(n) * P(gap_c = 0)", "",
           "| n | " + " | ".join(f"c = {c}" for c in cs) + " | "
           + " | ".join(f"n*(P - lim/sqrt(n)), c = {c}" for c in cs) + " |",
           "|---" * (1 + 2 * len(cs)) + "|"]
    tex += ["% Table 1: exact tie probabilities", r"\begin{tabular}{r" + "c" * (2 * len(cs)) + "}", r"\hline",
            r"$n$ & " + " & ".join(rf"$\sqrt{{n}}\,p_n^{{({c})}}$" for c in cs) + " & "
            + " & ".join(rf"$n\,\bigl(p_n^{{({c})}} - \sqrt{{\pi/8n}}\bigr)$" for c in cs) + r" \\ \hline"]
    for n in ns:
        scaled = [exact[(n, c, 1)] * math.sqrt(n) for c in cs]
        second = [n * (exact[(n, c, 1)] - LIMIT / math.sqrt(n)) for c in cs]
        md.append(f"| {fmt_n(n)} | " + " | ".join(f"{x:.5f}" for x in scaled) + " | "
                  + " | ".join(f"{x:+.4f}" for x in second) + " |")
        tex.append(f"{fmt_n(n).replace(',', '{,}')} & " + " & ".join(f"{x:.5f}" for x in scaled) + " & "
                   + " & ".join(f"${x:+.4f}$" for x in second) + r" \\")
    md += [f"| limit | " + " | ".join(f"{LIMIT:.5f}" for _ in cs) + " |" + " |" * len(cs), ""]
    tex += [r"\hline", r"limit & " + " & ".join(f"{LIMIT:.5f}" for _ in cs) + " &" * len(cs) + r" \\ \hline",
            r"\end{tabular}", ""]

    # ---- Table 2: exact gap probabilities ---------------------------------
    cols = [(c, j) for j in js if j >= 2 for c in cs]
    md += ["## Table 2. Exact gap probabilities: sqrt(n) * P(gap_c < j)", "",
           "| n | " + " | ".join(f"j = {j}, c = {c}" for c, j in cols) + " |", "|---" * (1 + len(cols)) + "|"]
    tex += ["% Table 2: exact gap probabilities", r"\begin{tabular}{r" + "c" * len(cols) + "}", r"\hline",
            r"$n$ & " + " & ".join(rf"$j={j},\ c={c}$" for c, j in cols) + r" \\ \hline"]
    for n in ns:
        vals = [exact[(n, c, j)] * math.sqrt(n) for c, j in cols]
        md.append(f"| {fmt_n(n)} | " + " | ".join(f"{x:.4f}" for x in vals) + " |")
        tex.append(f"{fmt_n(n).replace(',', '{,}')} & " + " & ".join(f"{x:.4f}" for x in vals) + r" \\")
    lims = [(2 * j - 1) * LIMIT for _, j in cols]
    md += ["| limit | " + " | ".join(f"{x:.4f}" for x in lims) + " |", ""]
    tex += [r"\hline", "limit & " + " & ".join(f"{x:.4f}" for x in lims) + r" \\ \hline", r"\end{tabular}", ""]

    # ---- Consistency of the simulation with the exact values --------------
    md += ["## Check. Simulation against exact values (standard scores)", "",
           "For each n, c and j the simulated number of mappings with gap_c < j is compared with",
           "trials * (exact probability); z = (observed - expected) / standard deviation.", "",
           "| n | trials | " + " | ".join(f"z (c = {c}, j = {j})" for c in cs for j in js) + " |",
           "|---" * (2 + len(cs) * len(js)) + "|"]
    zmax = 0.0
    for n in ns:
        zs = []
        for c in cs:
            s = sim[(n, c, alphas[0])]
            m = s["trials"]
            observed = {1: s["gap0"], 2: s["gap0"] + s["gap1"], 3: s["gap0"] + s["gap1"] + s["gap2"]}
            for j in js:
                p = exact[(n, c, j)]
                z = (observed[j] - m * p) / math.sqrt(m * p * (1 - p))
                zs.append(z)
                zmax = max(zmax, abs(z))
        md.append(f"| {fmt_n(n)} | {fmt_n(sim[(n, cs[0], alphas[0])]['trials'])} | "
                  + " | ".join(f"{z:+.2f}" for z in zs) + " |")
    md += ["", f"Largest |z|: {zmax:.2f} (over {len(ns) * len(cs) * len(js)} comparisons).", ""]

    # ---- Table 3: crown failure probabilities -----------------------------
    for c in cs:
        md += [f"## Table 3 (c = {c}). Simulated sqrt(n) * P(not |H| > alpha r > 0), with 95% Wilson intervals", "",
               "| n | trials | " + " | ".join(f"alpha = {a:g}" for a in alphas) + " |",
               "|---" * (2 + len(alphas)) + "|"]
        tex += [f"% Table 3 (c = {c}): crown failure probabilities, simulation",
                r"\begin{tabular}{rr" + "c" * len(alphas) + "}", r"\hline",
                r"$n$ & trials & " + " & ".join(rf"$\alpha = {a:g}$" for a in alphas) + r" \\ \hline"]
        for n in ns:
            cells_md, cells_tex = [], []
            for a in alphas:
                s = sim[(n, c, a)]
                k, m, rt = s["gap0"] + s["sparse"], s["trials"], math.sqrt(n)
                lo, hi = wilson(k, m)
                cells_md.append(f"{k / m * rt:.3f} [{lo * rt:.3f}, {hi * rt:.3f}]")
                cells_tex.append(f"{k / m * rt:.3f} $\\pm$ {(hi - lo) / 2 * rt:.3f}")
            trials = fmt_n(sim[(n, c, alphas[0])]["trials"])
            md.append(f"| {fmt_n(n)} | {trials} | " + " | ".join(cells_md) + " |")
            tex.append(f"{fmt_n(n).replace(',', '{,}')} & {trials.replace(',', '{,}')} & "
                       + " & ".join(cells_tex) + r" \\")
        md.append("")
        tex += [r"\hline", r"\end{tabular}", ""]

    with open("results/tables.md", "w") as fh:
        fh.write("\n".join(md))
    with open("results/tables.tex", "w") as fh:
        fh.write("\n".join(tex))
    print("\n".join(md))


if __name__ == "__main__":
    main()
