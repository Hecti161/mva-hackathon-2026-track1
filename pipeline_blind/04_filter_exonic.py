#!/usr/bin/env python3
"""
BLIND step 3 - Exonic/splice-region restriction.

The 200-gene candidate set yields ~32k variants, the vast majority deep intronic.
This step fetches the canonical transcript exon structure for each candidate gene
from Ensembl and keeps only variants falling in an exon or within 10 bp of a
splice site.

This is an exome-style triage, applied blind: no variant is inspected, only
coordinates. Its cost is that deep-intronic and regulatory variation is dropped,
which is declared as a limitation rather than hidden.
"""
import collections
import json
import os
import pickle
import sys
import time
import urllib.request

WORK = os.environ.get('WORK_DIR', 'work')
SPLICE_PAD = 10


def post_expand(symbols):
    url = ('https://rest.ensembl.org/lookup/symbol/homo_sapiens'
           '?content-type=application/json&expand=1')
    body = json.dumps({'symbols': symbols}).encode()
    req = urllib.request.Request(
        url, data=body,
        headers={'Content-Type': 'application/json', 'Accept': 'application/json'})
    for attempt in range(4):
        try:
            return json.load(urllib.request.urlopen(req, timeout=240))
        except Exception as e:
            print('  retry %d: %s' % (attempt, e), file=sys.stderr)
            time.sleep(5 * (attempt + 1))
    return {}


def main():
    rows = pickle.load(open(os.path.join(WORK, 'blind_rows.pkl'), 'rb'))
    genes = sorted({r[0] for r in rows})
    print('fetching exon structure for %d genes' % len(genes))

    exons = collections.defaultdict(list)
    for i in range(0, len(genes), 40):
        chunk = genes[i:i + 40]
        res = post_expand(chunk)
        for g, d in res.items():
            transcripts = d.get('Transcript') or []
            canon = next((t for t in transcripts if t.get('is_canonical')), None)
            if canon is None and transcripts:
                canon = max(transcripts, key=lambda t: len(t.get('Exon') or []))
            for ex in (canon or {}).get('Exon') or []:
                exons[g].append((int(ex['start']) - SPLICE_PAD, int(ex['end']) + SPLICE_PAD))
        print('  %d/%d genes' % (min(i + 40, len(genes)), len(genes)), file=sys.stderr)

    n_exons = sum(len(v) for v in exons.values())
    print('exons collected: %d across %d genes' % (n_exons, len(exons)))

    kept = []
    for r in rows:
        gene, _c, pos = r[0], r[1], r[2]
        for s, e in exons.get(gene, ()):
            if s <= pos <= e:
                kept.append(r)
                break

    pickle.dump(kept, open(os.path.join(WORK, 'blind_exonic.pkl'), 'wb'))
    print('\nvariants before exonic filter : %d' % len(rows))
    print('variants after exonic filter  : %d' % len(kept))
    by_gene = collections.Counter(r[0] for r in kept)
    print('genes retaining variants      : %d' % len(by_gene))


if __name__ == '__main__':
    main()
