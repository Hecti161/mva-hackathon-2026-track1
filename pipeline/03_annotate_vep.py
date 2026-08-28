#!/usr/bin/env python3
"""
Step 3 - Annotation.

Annotates the extracted candidate variants through the Ensembl VEP REST API, requesting:
consequence terms, MANE Select transcript mapping, HGVS c./p. nomenclature, exon numbering,
protein domains, gnomAD exome + genome allele frequencies, and SIFT / PolyPhen-2 /
AlphaMissense predictions.

Usage:  python pipeline/03_annotate_vep.py GENE[,GENE...]
        python pipeline/03_annotate_vep.py BUB1B,CEP57,TRIP13,BUB1,CENPE
Output: work/vep_<FIRSTGENE>.json
        work/rows.pkl   (parsed candidate records, reused by step 4)
"""
import collections
import json
import os
import pickle
import sys
import time
import urllib.request

WORK = os.environ.get('WORK_DIR', 'work')
BATCH = 150

URL = ('https://rest.ensembl.org/vep/human/region?content-type=application/json'
       '&canonical=1&mane=1&hgvs=1&pick_order=mane_select,canonical,biotype'
       '&per_gene=1&af_gnomade=1&af_gnomadg=1&sift=b&polyphen=b&numbers=1'
       '&domains=1&AlphaMissense=1&vcf_string=1')


def norm(c):
    return c[3:] if c.lower().startswith('chr') else c


def load_rows():
    """Parse work/candidates.vcf into (gene, chrom, pos, ref, alt, qual, filter, GT, sample)."""
    coords = json.load(open(os.path.join(WORK, 'gene_coords.json')))
    pad = 25_000
    by_chrom = collections.defaultdict(list)
    for gene, (c, s, e, _d) in coords.items():
        by_chrom[norm(str(c))].append((int(s) - pad, int(e) + pad, gene))

    rows = []
    for line in open(os.path.join(WORK, 'candidates.vcf')):
        if line[0] == '#':
            continue
        f = line.rstrip('\n').split('\t')
        chrom, pos = norm(f[0]), int(f[1])
        hits = [g for s, e, g in by_chrom[chrom] if s <= pos <= e]
        rows.append((hits[0], chrom, pos, f[3], f[4], f[5], f[6],
                     f[9].split(':')[0], f[9]))
    pickle.dump(rows, open(os.path.join(WORK, 'rows.pkl'), 'wb'))
    return rows


def post(batch):
    body = json.dumps({'variants': batch}).encode()
    req = urllib.request.Request(
        URL, data=body,
        headers={'Content-Type': 'application/json', 'Accept': 'application/json'})
    for attempt in range(4):
        try:
            return json.load(urllib.request.urlopen(req, timeout=180))
        except Exception as e:
            print('  retry %d: %s' % (attempt, e), file=sys.stderr)
            time.sleep(5 * (attempt + 1))
    return []


def main(genes):
    rows = load_rows()
    sel = [r for r in rows if r[0] in genes]
    print('annotating %d variants across %s' % (len(sel), ','.join(genes)), file=sys.stderr)

    res = []
    for i in range(0, len(sel), BATCH):
        chunk = sel[i:i + BATCH]
        res.extend(post(['%s %d . %s %s . . .' % (r[1], r[2], r[3], r[4]) for r in chunk]))
        print('  batch %d done' % (i // BATCH), file=sys.stderr)

    path = os.path.join(WORK, 'vep_%s.json' % genes[0])
    json.dump(res, open(path, 'w'))
    print('wrote %s (%d annotations)' % (path, len(res)), file=sys.stderr)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1].split(','))
