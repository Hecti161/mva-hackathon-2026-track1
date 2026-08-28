#!/usr/bin/env python3
"""
BLIND step 2 - Variant extraction over the phenotype-ranked candidate set.

Takes the top K genes from step 1 (phenotype ranking only - no human gene picking),
resolves their GRCh38 coordinates from Ensembl, and extracts every variant in those
loci from the proband VCF in a single streaming pass.

K is fixed at 200 (top 3.8% of the 5,268 ranked genes) before looking at any variant.

Usage:  python pipeline_blind/03_extract_topgenes.py <path/to/proband.vcf.gz>
"""
import collections
import gzip
import json
import os
import pickle
import sys
import urllib.request

WORK = os.environ.get('WORK_DIR', 'work')
K = 200
PAD = 25_000


def norm(c):
    return c[3:] if c.lower().startswith('chr') else c


def lookup_batch(symbols):
    """Ensembl POST /lookup/symbol - up to 1000 symbols per call."""
    url = 'https://rest.ensembl.org/lookup/symbol/homo_sapiens?content-type=application/json'
    body = json.dumps({'symbols': symbols}).encode()
    req = urllib.request.Request(
        url, data=body,
        headers={'Content-Type': 'application/json', 'Accept': 'application/json'})
    return json.load(urllib.request.urlopen(req, timeout=180))


def main(vcf_path):
    ranking = pickle.load(open(os.path.join(WORK, 'gene_ranking.pkl'), 'rb'))
    top = [(s, g) for s, g in ranking[:K] if g and g != '-']
    pheno = {g: s for s, g in top}
    print('candidate genes from phenotype ranking: %d' % len(top))

    coords, missing = {}, []
    symbols = [g for _s, g in top]
    for i in range(0, len(symbols), 200):
        chunk = symbols[i:i + 200]
        try:
            res = lookup_batch(chunk)
        except Exception as e:
            print('  lookup failed for a chunk: %s' % e, file=sys.stderr)
            res = {}
        for g in chunk:
            d = res.get(g)
            if d and d.get('seq_region_name') and len(str(d['seq_region_name'])) <= 2:
                coords[g] = (norm(str(d['seq_region_name'])), int(d['start']), int(d['end']))
            else:
                missing.append(g)
    print('resolved to GRCh38: %d   unresolved (alias/scaffold): %d' % (len(coords), len(missing)))
    if missing:
        print('  unresolved: %s' % ', '.join(missing[:15]) + (' ...' if len(missing) > 15 else ''))

    by_chrom = collections.defaultdict(list)
    for g, (c, s, e) in coords.items():
        by_chrom[c].append((s - PAD, e + PAD, g))

    out_path = os.path.join(WORK, 'blind_candidates.vcf')
    rows, n_total = [], 0
    with gzip.open(vcf_path, 'rt') as fh, open(out_path, 'w') as out:
        for line in fh:
            if line[0] == '#':
                out.write(line)
                continue
            n_total += 1
            t1 = line.index('\t')
            c = norm(line[:t1])
            iv = by_chrom.get(c)
            if not iv:
                continue
            t2 = line.index('\t', t1 + 1)
            p = int(line[t1 + 1:t2])
            for s, e, g in iv:
                if s <= p <= e:
                    out.write(line)
                    f = line.rstrip('\n').split('\t')
                    rows.append((g, c, p, f[3], f[4], f[6], f[9]))
                    break

    pickle.dump(rows, open(os.path.join(WORK, 'blind_rows.pkl'), 'wb'))
    json.dump(pheno, open(os.path.join(WORK, 'pheno_scores.json'), 'w'))
    print('\ntotal variants in VCF : %d' % n_total)
    print('extracted             : %d  across %d genes' % (len(rows), len(coords)))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
