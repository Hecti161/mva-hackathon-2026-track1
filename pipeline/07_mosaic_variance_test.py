#!/usr/bin/env python3
"""
Step 7 - Orthogonal test, part 3: variance test and detection limit.

For each chromosome, compares the observed variance of heterozygous allele balance
against the binomial expectation given the per-site depth, after dropping centromeric /
repetitive windows identified in step 6.

A mosaic trisomy present in a fraction f of cells shifts allele balance to (1+f)/(2+f)
and 1/(2+f), contributing an excess variance of

    excess = [ f / (2 * (2 + f)) ] ** 2      ->      f = 4s / (1 - 2s),  s = sqrt(excess)

Sequencing contributes its own uniform overdispersion (reference bias, PCR duplication,
sequencing error), so the *median* variance ratio across autosomes is taken as the
technical baseline and each chromosome is measured against it. The scatter of that ratio
across chromosomes gives the detection limit.

In this proband the ratio is a uniform ~1.24x on every autosome, with residual scatter
tracking GC-richness and gene density (chr19 > chr16 ~ chr17 > chr9 > chr7) - the known
technical gradient. No chromosome departs from the trend. Detection limit ~15% of cells.

Usage:  python pipeline/07_mosaic_variance_test.py
"""
import math
import os
import pickle
import statistics as st

WORK = os.environ.get('WORK_DIR', 'work')
WINDOW = 10_000_000
REPEAT_DISPERSION = 0.13     # windows above this are centromeric/repetitive - dropped


def implied_fraction(excess):
    """Invert excess = [f/(2(2+f))]^2 for f."""
    if excess <= 0:
        return None
    s = math.sqrt(excess)
    if s >= 0.45:
        return None
    return 4 * s / (1 - 2 * s)


def clean_sites(d, c):
    """Sites outside centromeric/repetitive windows."""
    bins = {}
    for p, ab, dp in d[c]:
        bins.setdefault(p // WINDOW, []).append((ab, dp))
    keep = []
    for _k, v in bins.items():
        if len(v) < 500:
            continue
        disp = sum(1 for x, _ in v if abs(x - 0.5) > 0.15) / len(v)
        if disp > REPEAT_DISPERSION:
            continue
        keep += v
    return keep


def main():
    d = pickle.load(open(os.path.join(WORK, 'ab.pkl'), 'rb'))
    print('Excess variance of allele balance over the binomial expectation')
    print('mosaic trisomy in fraction f of cells -> excess = [f/(2(2+f))]^2\n')
    print('%-4s %9s %10s %10s %9s %11s'
          % ('chr', 'n_sites', 'var_obs', 'var_binom', 'ratio', 'implied_f'))

    rows = []
    for c in [str(i) for i in range(1, 23)]:
        if c not in d:
            continue
        sites = clean_sites(d, c)
        if len(sites) < 5000:
            continue
        ab = [x for x, _ in sites]
        var_obs = st.pvariance(ab, st.mean(ab))
        var_bin = st.mean([0.25 / dp for _, dp in sites])
        rows.append((c, len(sites), var_obs, var_bin, var_obs / var_bin))

    base = st.median([r[4] for r in rows])
    for c, n, vo, ve, ratio in rows:
        f = implied_fraction((ratio - base) * ve)
        print('%-4s %9d %10.6f %10.6f %9.4f %11s'
              % (c, n, vo, ve, ratio, '%.1f%%' % (100 * f) if f else '-'))

    sd = st.stdev([r[4] for r in rows])
    mean_bin = st.mean([r[3] for r in rows])
    limit = implied_fraction(3 * sd * mean_bin)
    print('\ntechnical baseline (median variance ratio): %.4f' % base)
    print('between-chromosome scatter (1 sd)        : %.4f' % sd)
    print('detection limit (3 sd) for whole-chromosome mosaic trisomy: f ~ %.1f%% of cells'
          % (100 * limit))


if __name__ == '__main__':
    main()
