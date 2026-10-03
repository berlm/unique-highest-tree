/*
 * mapping_sim.c -- gap and crown statistics of c-branches of random mappings.
 *
 * A mapping g of {0..n-1} is viewed as a functional graph.  The depth of a
 * vertex is its distance to the cycles.  For c >= 0, a c-branch is the subtree
 * rooted at a vertex of depth c; its height is the largest distance from its
 * root to one of its vertices.
 *
 *   gap_c = 0          if there is no c-branch or the maximal height is
 *                      attained by at least two c-branches,
 *   gap_c = H1 - H2    if there is a unique highest c-branch (height H1) and
 *                      H2 is the largest height of the other c-branches,
 *   gap_c = infinity   if there is exactly one c-branch (then H2 := -1).
 *
 * If gap_c >= 1, the crown consists of the vertices of the highest c-branch at
 * distance >= H2 + 1 from its root; r is the number of those at distance
 * exactly H2 + 1.  The crown is called sparse (for alpha) if |crown| <= alpha*r.
 *
 * Usage:
 *   mapping_sim sim  <n> <trials> <seed> <cmax> <alpha> [<alpha> ...]
 *   mapping_sim enum <n> <cmax> <alpha> [<alpha> ...]      (all n^n mappings)
 *
 * Output (CSV, one line per c and alpha):
 *   n,c,alpha,trials,seed,gap0,gap1,gap2,single,sparse
 * where gap0/gap1/gap2 count the mappings with gap_c = 0/1/2, single counts the
 * mappings with exactly one c-branch, and sparse counts the mappings with
 * gap_c >= 1 and a sparse crown.  In enum mode, trials = n^n and seed = 0.
 *
 * Each mapping is processed in O(n * (cmax + 1)) time: leaves are removed
 * repeatedly (Kahn's algorithm), which yields the cyclic vertices, the height
 * of the subtree below every vertex and an order in which every vertex comes
 * after its descendants.
 */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ALPHA 16

/* ---- xoshiro256** (Blackman, Vigna), seeded with splitmix64 ---- */
static uint64_t rng_s[4];

static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static void rng_seed(uint64_t seed) {
    for (int i = 0; i < 4; i++) rng_s[i] = splitmix64(&seed);
}

static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static inline uint64_t rng_next(void) {
    uint64_t result = rotl(rng_s[1] * 5, 7) * 9, t = rng_s[1] << 17;
    rng_s[2] ^= rng_s[0];
    rng_s[3] ^= rng_s[1];
    rng_s[1] ^= rng_s[2];
    rng_s[0] ^= rng_s[3];
    rng_s[2] ^= t;
    rng_s[3] = rotl(rng_s[3], 45);
    return result;
}

/* unbiased integer in [0, n) (Lemire's method) */
static inline uint32_t rng_below(uint32_t n) {
    uint64_t m = (uint64_t)(uint32_t)(rng_next() >> 32) * n;
    uint32_t lo = (uint32_t)m;
    if (lo < n) {
        uint32_t thresh = (uint32_t)(-n) % n;
        while (lo < thresh) {
            m = (uint64_t)(uint32_t)(rng_next() >> 32) * n;
            lo = (uint32_t)m;
        }
    }
    return (uint32_t)(m >> 32);
}

/* ---- per-mapping analysis ---- */
static int n, cmax, n_alpha;
static double alpha[MAX_ALPHA];
static int32_t *f, *indeg, *order, *hgt, *depth, *anc;
static uint64_t *cnt_gap0, *cnt_gap1, *cnt_gap2, *cnt_single, *cnt_sparse; /* [c], [c][alpha] */

static void analyse(void) {
    int head = 0, tail = 0;
    memset(indeg, 0, (size_t)n * sizeof *indeg);
    memset(hgt, 0, (size_t)n * sizeof *hgt);
    for (int x = 0; x < n; x++) indeg[f[x]]++;
    for (int x = 0; x < n; x++)
        if (indeg[x] == 0) order[tail++] = x;
    while (head < tail) {                     /* remove leaves; children before parents */
        int x = order[head++], p = f[x];
        if (hgt[x] + 1 > hgt[p]) hgt[p] = hgt[x] + 1;
        if (--indeg[p] == 0) order[tail++] = p;
    }
    /* vertices never removed (indeg > 0) are cyclic */
    for (int x = 0; x < n; x++) depth[x] = (indeg[x] > 0) ? 0 : -1;
    for (int i = tail - 1; i >= 0; i--) {     /* parents before children */
        int x = order[i];
        depth[x] = depth[f[x]] + 1;
    }
    for (int c = 0; c <= cmax; c++) {
        int h1 = -1, h2 = -1, top = -1, branches = 0, ties = 0;
        for (int x = 0; x < n; x++) {
            if (depth[x] != c) continue;
            branches++;
            if (hgt[x] > h1) { h2 = h1; h1 = hgt[x]; top = x; ties = 1; }
            else if (hgt[x] == h1) { ties++; }
            else if (hgt[x] > h2) { h2 = hgt[x]; }
        }
        if (branches == 0 || ties >= 2) { cnt_gap0[c]++; continue; }
        if (branches == 1) cnt_single[c]++;   /* h2 = -1 */
        else if (h1 - h2 == 1) cnt_gap1[c]++;
        else if (h1 - h2 == 2) cnt_gap2[c]++;
        /* crown of the branch rooted at top: vertices of depth >= c + h2 + 1 */
        for (int x = 0; x < n; x++) anc[x] = (depth[x] == c) ? x : -1;
        for (int i = tail - 1; i >= 0; i--) {
            int x = order[i];
            if (depth[x] > c) anc[x] = anc[f[x]];
        }
        int64_t size = 0, roots = 0;
        int level = c + h2 + 1;
        for (int x = 0; x < n; x++) {
            if (anc[x] != top || depth[x] < level) continue;
            size++;
            if (depth[x] == level) roots++;
        }
        for (int a = 0; a < n_alpha; a++)
            if ((double)size <= alpha[a] * (double)roots) cnt_sparse[c * MAX_ALPHA + a]++;
    }
}

static void *xcalloc(size_t k, size_t s) {
    void *p = calloc(k, s);
    if (!p) { fprintf(stderr, "out of memory\n"); exit(1); }
    return p;
}

int main(int argc, char **argv) {
    int sim = argc >= 2 && strcmp(argv[1], "sim") == 0;
    int enumerate = argc >= 2 && strcmp(argv[1], "enum") == 0;
    if ((sim && argc < 7) || (enumerate && argc < 5) || (!sim && !enumerate)) {
        fprintf(stderr,
                "usage: %s sim <n> <trials> <seed> <cmax> <alpha> [<alpha> ...]\n"
                "       %s enum <n> <cmax> <alpha> [<alpha> ...]\n", argv[0], argv[0]);
        return 2;
    }
    uint64_t trials = 0, seed = 0;
    int arg = 2;
    n = atoi(argv[arg++]);
    if (sim) {
        trials = strtoull(argv[arg++], NULL, 10);
        seed = strtoull(argv[arg++], NULL, 10);
    }
    cmax = atoi(argv[arg++]);
    n_alpha = argc - arg;
    if (n < 1 || cmax < 0 || n_alpha < 1 || n_alpha > MAX_ALPHA) {
        fprintf(stderr, "bad arguments\n");
        return 2;
    }
    for (int a = 0; a < n_alpha; a++) alpha[a] = atof(argv[arg + a]);

    f = xcalloc((size_t)n, sizeof *f);
    indeg = xcalloc((size_t)n, sizeof *indeg);
    order = xcalloc((size_t)n, sizeof *order);
    hgt = xcalloc((size_t)n, sizeof *hgt);
    depth = xcalloc((size_t)n, sizeof *depth);
    anc = xcalloc((size_t)n, sizeof *anc);
    cnt_gap0 = xcalloc((size_t)cmax + 1, sizeof *cnt_gap0);
    cnt_gap1 = xcalloc((size_t)cmax + 1, sizeof *cnt_gap1);
    cnt_gap2 = xcalloc((size_t)cmax + 1, sizeof *cnt_gap2);
    cnt_single = xcalloc((size_t)cmax + 1, sizeof *cnt_single);
    cnt_sparse = xcalloc(((size_t)cmax + 1) * MAX_ALPHA, sizeof *cnt_sparse);

    if (sim) {
        rng_seed(seed);
        for (uint64_t it = 0; it < trials; it++) {
            for (int x = 0; x < n; x++) f[x] = (int32_t)rng_below((uint32_t)n);
            analyse();
        }
    } else {
        if (n > 9) { fprintf(stderr, "enum mode needs n <= 9\n"); return 2; }
        for (;;) {                             /* f runs through all n^n mappings */
            analyse();
            trials++;
            int x = 0;
            while (x < n && ++f[x] == n) f[x++] = 0;
            if (x == n) break;
        }
    }
    for (int c = 0; c <= cmax; c++)
        for (int a = 0; a < n_alpha; a++)
            printf("%d,%d,%g,%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 "\n",
                   n, c, alpha[a], trials, seed, cnt_gap0[c], cnt_gap1[c], cnt_gap2[c],
                   cnt_single[c], cnt_sparse[c * MAX_ALPHA + a]);
    return 0;
}
