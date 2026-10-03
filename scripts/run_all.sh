#!/bin/sh
# Reproduces the files in results/.  Run from the repository root.
set -e
make -C sim
mkdir -p results

# exact gap probabilities (about one minute)
python3 exact/exact_gap.py 100 250 1000 2500 10000 25000 -o results/exact_gap.csv

# Monte Carlo: gap and crown statistics for c = 0, 1, 2 and alpha = 1.5, 2, 3
# (about ten minutes on one core)
OUT=results/simulation.csv
echo "n,c,alpha,trials,seed,gap0,gap1,gap2,single,sparse" > "$OUT"
./sim/mapping_sim sim   100 4000000 2026100301 2 1.5 2 3 >> "$OUT"
./sim/mapping_sim sim   250 4000000 2026100302 2 1.5 2 3 >> "$OUT"
./sim/mapping_sim sim  1000 2000000 2026100303 2 1.5 2 3 >> "$OUT"
./sim/mapping_sim sim  2500 1000000 2026100304 2 1.5 2 3 >> "$OUT"
./sim/mapping_sim sim 10000  400000 2026100305 2 1.5 2 3 >> "$OUT"
./sim/mapping_sim sim 25000  200000 2026100306 2 1.5 2 3 >> "$OUT"

python3 scripts/make_tables.py
