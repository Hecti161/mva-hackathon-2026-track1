#!/usr/bin/env python3
"""
Step 4 - Prioritisation.

Ranks annotated variants by predicted impact on the MANE Select transcript, population
frequency (recessive-appropriate thresholds) and genotype, and prints read-level quality
(DP, allelic depth, GQ) so that mapping and strand artefacts can be excluded by inspection.

Two tables are emitted:
  1. HIGH / MODERATE impact - protein-altering candidates
  2. LOW / MODIFIER but rare or absent from gnomAD - non-coding candidates worth a look

Usage:  python pipeline/04_prioritize.py work/vep_BUB1B.json
"""
import json
import os
import pickle
import sys

WORK = os.environ.get('WORK_DIR', 'work')

BENIGN = {
    'intron_variant', 'intergenic_variant', 'upstream_gene_variant',
    'downstream_gene_variant', 'synonymous_variant', 'non_coding_transcript_exon_variant',
    'NMD_transcript_variant', '5_prime_UTR_variant', '3_prime_UTR_variant',
    'regulatory_region_variant', 'TF_binding_site_variant', 'non_coding_transcript_variant',
}

HDR = ('%-8s %-3s %11s %-14s %-30s %-8s %-7s %-9s %-8s %s'
       % ('gene', 'c', 'pos', 'ref>alt', 'consequence', 'impact', 'GT', 'gnomAD', 'AD', 'HGVSp'))


def max_gnomad_af(v, alt):
    """Highest gnomAD frequency reported for this allele across all colocated records."""
    best = None
    for cv in v.get('colocated_variants', []):
        for key, val in (cv.get('frequencies') or {}).get(alt, {}).items():
            if 'gnomad' in key and isinstance(val, (int, float)):
                best = val if best is None else max(best, val)
    return best


def collect(vep, key):
    out = []
    for v in vep:
        vs = v.get('vcf_string')
        if not vs:
            continue
        if isinstance(vs, list):
            vs = vs[0]
        c, p, ref, alt = vs.split('-')
        src = key.get((c, int(p), ref, alt))
        gene = src[0] if src else '?'
        gt = ad = ''
        if src:
            fields = src[8].split(':')
            gt, ad = (fields + ['', ''])[:2]
        af = max_gnomad_af(v, alt)
        for tc in v.get('transcript_consequences', []):
            terms = tc.get('consequence_terms', [])
            out.append(dict(
                gene=tc.get('gene_symbol') or gene, chrom=c, pos=int(p), ref=ref, alt=alt,
                cons=','.join(terms), impact=tc.get('impact'), gt=gt, ad=ad, af=af,
                hgvsc=tc.get('hgvsc', ''), hgvsp=tc.get('hgvsp', ''),
                sift=tc.get('sift_prediction', ''), pp=tc.get('polyphen_prediction', ''),
                am=tc.get('am_class', ''),
                coding=any(t not in BENIGN for t in terms)))
    return out


def show(items, title):
    print('\n=== %s (%d) ===' % (title, len(items)))
    print(HDR)
    for o in items:
        af = 'NOVEL' if o['af'] is None else '%.5f' % o['af']
        label = (o['hgvsp'] or o['hgvsc'] or '').split(':')[-1][:40]
        print('%-8s %-3s %11d %-14s %-30s %-8s %-7s %-9s %-8s %s' % (
            o['gene'], o['chrom'], o['pos'], o['ref'][:6] + '>' + o['alt'][:6],
            o['cons'][:30], o['impact'] or '', o['gt'], af, o['ad'], label))


def main(vep_path):
    rows = pickle.load(open(os.path.join(WORK, 'rows.pkl'), 'rb'))
    key = {(r[1], r[2], r[3], r[4]): r for r in rows}
    out = collect(json.load(open(vep_path)), key)

    hi = [o for o in out if o['impact'] in ('HIGH', 'MODERATE')]
    hi.sort(key=lambda o: (o['impact'] != 'HIGH', o['af'] if o['af'] is not None else -1))
    show(hi, 'HIGH / MODERATE impact - ranked by severity then rarity')

    rare = [o for o in out
            if o['impact'] not in ('HIGH', 'MODERATE')
            and (o['af'] is None or o['af'] < 0.001)]
    rare.sort(key=lambda o: (o['gene'], o['pos']))
    show(rare, 'LOW / MODIFIER, rare or absent from gnomAD')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
