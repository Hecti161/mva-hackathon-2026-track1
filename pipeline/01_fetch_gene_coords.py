#!/usr/bin/env python3
"""
Step 1 - Candidate locus definition.

Retrieves GRCh38 coordinates for the candidate gene panel from the Ensembl REST API.
Coordinates are fetched live rather than hard-coded so the panel cannot silently drift
out of sync with the assembly.

Two gene groups are queried:
  - MVA / mitotic spindle assembly checkpoint (SAC) genes
  - childhood-cancer and chromosomal-instability differential diagnosis

Output: work/gene_coords.json
"""
import json
import os
import sys
import urllib.request

# MVA / spindle assembly checkpoint
SAC = ['BUB1B', 'CEP57', 'TRIP13', 'BUB1', 'CENPE', 'BUB3', 'MAD1L1', 'MAD2L1',
       'CDC20', 'PLK4', 'TUBB', 'KNL1', 'TTK', 'ZWILCH', 'CENPF']
# differential diagnosis: childhood rhabdomyosarcoma / chromosomal instability
DIFF = ['TP53', 'DICER1', 'ATM', 'NBN', 'MRE11']

GENES = SAC + DIFF
WORK = os.environ.get('WORK_DIR', 'work')


def lookup(symbol):
    url = ('https://rest.ensembl.org/lookup/symbol/homo_sapiens/%s'
           '?content-type=application/json' % symbol)
    d = json.load(urllib.request.urlopen(url, timeout=30))
    assert d.get('assembly_name') == 'GRCh38', 'unexpected assembly for %s' % symbol
    return d['seq_region_name'], d['start'], d['end'], d.get('description', '')[:60]


def main():
    os.makedirs(WORK, exist_ok=True)
    out = {}
    for g in GENES:
        try:
            out[g] = lookup(g)
            print('%-8s %-3s %10d %10d' % (g, out[g][0], out[g][1], out[g][2]))
        except Exception as e:
            print('%-8s FAILED: %s' % (g, e), file=sys.stderr)
    path = os.path.join(WORK, 'gene_coords.json')
    json.dump(out, open(path, 'w'), indent=1)
    print('\nwrote %s (%d genes)' % (path, len(out)))


if __name__ == '__main__':
    main()
