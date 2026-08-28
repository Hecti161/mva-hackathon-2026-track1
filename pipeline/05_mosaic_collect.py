#!/usr/bin/env python3
"""
Step 5 - Orthogonal test of the cellular phenotype, part 1: collection.

MVA is defined by a *cellular* phenotype - mosaic aneuploidy across multiple chromosomes -
not only by a genotype. That phenotype is measurable from the VCF alone, without
alignments: if a chromosome is trisomic in a fraction f of cells, heterozygous SNVs on it
split their allele balance away from 0.5 toward (1+f)/(2+f) and 1/(2+f), inflating the
*variance* of allele balance while leaving its median at 0.5.

This step harvests allele balance and depth for every high-confidence heterozygous SNV
genome-wide. Steps 6 and 7 test it.

Filters: PASS, biallelic SNV, GT het, GQ >= 60, 15 <= DP <= 80.
The depth window excludes collapsed repeats and extreme copy-number regions, which
would otherwise dominate the variance.

Usage:  python pipeline/05_mosaic_collect.py <path/to/proband.vcf.gz>
Output: work/ab.pkl   (NOT for redistribution - derived from proband genotypes)
"""
import collections
import gzip
import os
import pickle
import sys

WORK = os.environ.get('WORK_DIR', 'work')
HET = {'0/1', '0|1', '1|0'}
GQ_MIN, DP_MIN, DP_MAX = 60, 15, 80


def norm(c):
    return c[3:] if c.lower().startswith('chr') else c


def main(vcf_path):
    data = collections.defaultdict(list)
    n = 0
    for line in gzip.open(vcf_path, 'rt'):
        if line[0] == '#':
            continue
        f = line.split('\t', 10)
        if f[6] != 'PASS':
            continue
        if len(f[3]) != 1 or len(f[4]) != 1:      # biallelic SNV only
            continue
        fmt, sample = f[8].split(':'), f[9].split(':')
        if sample[0] not in HET:
            continue
        try:
            ad = sample[fmt.index('AD')].split(',')
            gq = int(sample[fmt.index('GQ')])
            ref, alt = int(ad[0]), int(ad[1])
        except (ValueError, IndexError):
            continue
        dp = ref + alt
        if gq < GQ_MIN or not (DP_MIN <= dp <= DP_MAX):
            continue
        data[norm(f[0])].append((int(f[1]), alt / dp, dp))
        n += 1

    path = os.path.join(WORK, 'ab.pkl')
    pickle.dump(dict(data), open(path, 'wb'))
    print('high-confidence heterozygous SNVs: %d' % n)
    print('wrote %s' % path)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
