#!/usr/bin/env python3
"""
Step 2 - Extraction.

Single streaming pass over the proband VCF, emitting only records that intersect the
candidate gene intervals (padded by +/-25 kb to retain promoter, UTR and near-regulatory
space).

Pure Python over gzip: bcftools/tabix are unavailable on Windows, and a full pass over
~5 M records takes a few minutes, which is acceptable for a one-shot analysis.

Contig naming is normalised, so a VCF using either '1' or 'chr1' works unchanged.

Usage:  python pipeline/02_extract_regions.py <path/to/proband.vcf.gz>
Output: work/candidates.vcf   (NOT for redistribution - contains proband genotypes)
        work/chromcount.json
"""
import collections
import gzip
import json
import os
import sys

PAD = 25_000
WORK = os.environ.get('WORK_DIR', 'work')


def norm(c):
    """'chr1' and '1' both normalise to '1'."""
    return c[3:] if c.lower().startswith('chr') else c


def main(vcf_path):
    coords = json.load(open(os.path.join(WORK, 'gene_coords.json')))
    by_chrom = collections.defaultdict(list)
    for gene, (c, s, e, _desc) in coords.items():
        by_chrom[norm(str(c))].append((int(s) - PAD, int(e) + PAD, gene))

    out_path = os.path.join(WORK, 'candidates.vcf')
    n_total = n_hit = 0
    chromcount = collections.Counter()

    with gzip.open(vcf_path, 'rt') as fh, open(out_path, 'w') as out:
        for line in fh:
            if line[0] == '#':
                out.write(line)
                continue
            n_total += 1
            t1 = line.index('\t')
            chrom = norm(line[:t1])
            chromcount[chrom] += 1
            intervals = by_chrom.get(chrom)
            if not intervals:
                continue
            t2 = line.index('\t', t1 + 1)
            pos = int(line[t1 + 1:t2])
            for start, end, _gene in intervals:
                if start <= pos <= end:
                    out.write(line)
                    n_hit += 1
                    break

    json.dump(dict(chromcount), open(os.path.join(WORK, 'chromcount.json'), 'w'), indent=1)
    print('total variants in VCF : %d' % n_total)
    print('extracted to %s : %d' % (out_path, n_hit))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
