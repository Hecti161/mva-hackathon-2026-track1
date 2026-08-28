#!/usr/bin/env python3
"""
Step 6 - Orthogonal test, part 2: per-chromosome and windowed dispersion.

Prints, per chromosome, the median allele balance and the fraction of heterozygous sites
with |AB - 0.5| > 0.15, then breaks each chromosome into 10 Mb windows.

The windowed view is the point of this step, and it is what prevents a false positive.
A chromosome-wide summary makes small, centromere-heavy chromosomes look aneuploid,
because centromeric and pericentromeric repeat inflates dispersion locally. Binning
separates the two cases:

  - uniform elevation across all windows -> genuine whole-chromosome event
  - elevation confined to a few windows  -> repeat/mappability artefact

In this proband every apparent signal localised to centromeric windows (chr1 120-140 Mb,
chr9 40-60 Mb spanning the 9q12 heterochromatin block, chr20 20-40 Mb, chr17 20-30 Mb,
chr22 0-20 Mb) and is reported as the artefact it is.

Usage:  python pipeline/06_mosaic_windows.py
"""
import os
import pickle
import statistics as st

WORK = os.environ.get('WORK_DIR', 'work')
WINDOW = 10_000_000
THRESH = 0.15
CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y']


def dispersion(abs_):
    return 100.0 * sum(1 for x in abs_ if abs(x - 0.5) > THRESH) / len(abs_)


def main():
    d = pickle.load(open(os.path.join(WORK, 'ab.pkl'), 'rb'))
    gdp = st.median([dp for c in CHROMS if c in d for _p, _ab, dp in d[c]])
    n = sum(len(d[c]) for c in CHROMS if c in d)
    print('heterozygous SNVs: %d | global median depth: %.1fx\n' % (n, gdp))

    print('%-4s %9s %8s %8s %9s %9s' % ('chr', 'n_het', 'medAB', 'disp%', 'medDP', 'DPratio'))
    summary = []
    for c in CHROMS:
        if c not in d or len(d[c]) < 200:
            continue
        ab = [x for _p, x, _dp in d[c]]
        dp = [z for _p, _x, z in d[c]]
        mdp = st.median(dp)
        summary.append((c, len(ab), st.median(ab), dispersion(ab), mdp, mdp / gdp))
    base = st.median([r[3] for r in summary])
    for c, cnt, med, disp, mdp, ratio in summary:
        flag = ' <== check' if abs(disp - base) > 3 else ''
        print('%-4s %9d %8.4f %8.2f %9.1f %9.3f%s' % (c, cnt, med, disp, mdp, ratio, flag))
    print('\nmedian autosomal dispersion: %.2f%%' % base)

    print('\n\nPer-window dispersion (%d Mb windows) - artefact discrimination' % (WINDOW // 10**6))
    print('uniform = technical | localised = repeat | whole-chromosome = biological\n')
    for c in CHROMS[:22]:
        if c not in d:
            continue
        bins = {}
        for p, ab, _dp in d[c]:
            bins.setdefault(p // WINDOW, []).append(ab)
        vals = [dispersion(v) for k, v in sorted(bins.items()) if len(v) >= 500]
        if not vals:
            continue
        print('chr%-3s n=%2d  min=%5.1f%%  median=%5.1f%%  max=%5.1f%%  | %s' % (
            c, len(vals), min(vals), st.median(vals), max(vals),
            ' '.join('%.0f' % v for v in vals)))


if __name__ == '__main__':
    main()
