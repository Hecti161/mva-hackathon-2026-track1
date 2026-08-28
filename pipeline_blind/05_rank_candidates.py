#!/usr/bin/env python3
"""
BLIND step 4 - Annotate and rank candidate findings.

Annotates the exonic candidate set through Ensembl VEP, then scores every gene
under explicit inheritance models and ranks the resulting findings by

    score = phenotype_score  x  variant_evidence

Nothing here names a gene or a disease. The ranking is produced from the
phenotype similarity of step 1 and the variant evidence in the proband's VCF.

Inheritance models evaluated per gene:
  homozygous     one rare damaging variant, genotype 1/1
  compound_het   two or more rare damaging heterozygous variants
  dominant       one rare damaging heterozygous variant (de novo status unknown)

Recessive models allow gnomAD AF < 1%; the dominant model requires AF < 0.1%,
reflecting that a dominant childhood-onset allele cannot be common.
"""
import collections
import json
import os
import pickle
import sys
import time
import urllib.request

WORK = os.environ.get('WORK_DIR', 'work')

URL = ('https://rest.ensembl.org/vep/human/region?content-type=application/json'
       '&canonical=1&mane=1&hgvs=1&per_gene=1&pick_order=mane_select,canonical,biotype'
       '&af_gnomade=1&af_gnomadg=1&sift=b&polyphen=b&numbers=1&AlphaMissense=1&vcf_string=1')

IMPACT = {'HIGH': 1.0, 'MODERATE': 0.5, 'LOW': 0.1, 'MODIFIER': 0.02}
AF_RECESSIVE, AF_DOMINANT = 0.01, 0.001
HET = {'0/1', '0|1', '1|0'}
HOM = {'1/1', '1|1'}


def post(batch):
    body = json.dumps({'variants': batch}).encode()
    req = urllib.request.Request(
        URL, data=body,
        headers={'Content-Type': 'application/json', 'Accept': 'application/json'})
    for attempt in range(4):
        try:
            return json.load(urllib.request.urlopen(req, timeout=240))
        except Exception as e:
            print('  retry %d: %s' % (attempt, e), file=sys.stderr)
            time.sleep(5 * (attempt + 1))
    return []


def max_af(v, alt):
    best = None
    for cv in v.get('colocated_variants', []):
        for k, val in (cv.get('frequencies') or {}).get(alt, {}).items():
            if 'gnomad' in k and isinstance(val, (int, float)):
                best = val if best is None else max(best, val)
    return best


def deleteriousness(tc):
    """0-1 boost from concordant in-silico predictions, for missense calls."""
    hits = sum([tc.get('sift_prediction', '').startswith('deleterious'),
                tc.get('polyphen_prediction', '') in ('probably_damaging', 'possibly_damaging'),
                tc.get('am_class', '') == 'likely_pathogenic'])
    return hits / 3.0


MIN_DP, MIN_GQ, MIN_MQ = 15, 50, 50
READ_LEN = 150


def load_quality():
    """QUAL / DP / GQ / MQ per variant, parsed from the extracted VCF."""
    q = {}
    path = os.path.join(WORK, 'blind_candidates.vcf')
    for line in open(path):
        if line[0] == '#':
            continue
        f = line.rstrip('\n').split('\t')
        fmt, s = f[8].split(':'), f[9].split(':')
        def g(tag, cast=int):
            try:
                return cast(s[fmt.index(tag)])
            except (ValueError, IndexError):
                return 0
        mq = 0.0
        for kv in f[7].split(';'):
            if kv.startswith('MQ='):
                try:
                    mq = float(kv[3:])
                except ValueError:
                    pass
        def gs(tag):
            try:
                return s[fmt.index(tag)]
            except (ValueError, IndexError):
                return ''
        ad = gs('AD')
        vaf = 0.0
        if ad and ',' in ad:
            try:
                parts = [int(x) for x in ad.split(',')]
                vaf = parts[1] / sum(parts) if sum(parts) else 0.0
            except ValueError:
                pass
        q[(f[0], int(f[1]), f[3], f[4])] = dict(
            qual=f[5], filt=f[6], dp=g('DP'), gq=g('GQ'), mq=mq,
            pgt=gs('PGT'), pid=gs('PID'), vaf=vaf)
    return q


def passes_quality(v, q):
    r = q.get((v['chrom'], v['pos'], v['ref'], v['alt']))
    if not r:
        return False
    if not (r['filt'] == 'PASS' and r['dp'] >= MIN_DP
            and r['gq'] >= MIN_GQ and r['mq'] >= MIN_MQ):
        return False
    # allele balance: a real het sits near 0.5, a real hom-alt near 1.0
    if v['gt'] in HOM:
        return r['vaf'] >= 0.90
    return 0.20 <= r['vaf'] <= 0.80


def phase_relation(a, b, q):
    """
    Physical phase from GATK's PGT/PID fields.

    HaplotypeCaller emits PID (phase-set id) and PGT (phased genotype) whenever two
    variants are close enough to be seen on the same reads. Two variants sharing a
    PID with the SAME PGT are on the same haplotype - in cis, one allele, not two.
    Sharing a PID with OPPOSITE PGT places them in trans and confirms a genuine
    compound heterozygote. Absent a shared PID they are too far apart to phase with
    short reads and the relationship is unknown.

    This is read-level evidence from the caller, not an inference from coordinates.
    """
    ra = q.get((a['chrom'], a['pos'], a['ref'], a['alt']))
    rb = q.get((b['chrom'], b['pos'], b['ref'], b['alt']))
    if not ra or not rb:
        return 'unknown'
    if ra['pid'] and ra['pid'] == rb['pid']:
        return 'cis' if ra['pgt'] == rb['pgt'] else 'trans'
    return 'unknown'


def main():
    rows = pickle.load(open(os.path.join(WORK, 'blind_exonic.pkl'), 'rb'))
    pheno = json.load(open(os.path.join(WORK, 'pheno_scores.json')))
    key = {(r[1], r[2], r[3], r[4]): r for r in rows}
    qual = load_quality()

    cache = os.path.join(WORK, 'blind_vep.json')
    if os.path.exists(cache):
        vep = json.load(open(cache))
    else:
        vep = []
        for i in range(0, len(rows), 150):
            chunk = rows[i:i + 150]
            vep.extend(post(['%s %d . %s %s . . .' % (r[1], r[2], r[3], r[4]) for r in chunk]))
            print('  annotated %d/%d' % (min(i + 150, len(rows)), len(rows)), file=sys.stderr)
        json.dump(vep, open(cache, 'w'))
    print('annotations: %d' % len(vep))

    # one best consequence per variant per gene
    per_gene = collections.defaultdict(list)
    for v in vep:
        vs = v.get('vcf_string')
        if not vs:
            continue
        if isinstance(vs, list):
            vs = vs[0]
        c, p, ref, alt = vs.split('-')
        src = key.get((c, int(p), ref, alt))
        if not src:
            continue
        gt = src[6].split(':')[0]
        af = max_af(v, alt)
        for tc in v.get('transcript_consequences', []):
            g = tc.get('gene_symbol')
            if g not in pheno:
                continue
            imp = IMPACT.get(tc.get('impact'), 0.0)
            if tc.get('impact') == 'MODERATE':
                imp *= (0.4 + 0.6 * deleteriousness(tc))
            per_gene[g].append(dict(
                chrom=c, pos=int(p), ref=ref, alt=alt, gt=gt, af=af, impact=imp,
                cons=','.join(tc.get('consequence_terms', [])),
                hgvsp=(tc.get('hgvsp') or tc.get('hgvsc') or '').split(':')[-1]))

    findings = []
    n_dropped = 0
    for g, vs in per_gene.items():
        clean = [v for v in vs if passes_quality(v, qual)]
        n_dropped += len(vs) - len(clean)
        rec = [v for v in clean if v['impact'] >= 0.3 and (v['af'] is None or v['af'] < AF_RECESSIVE)]
        dom = [v for v in clean if v['impact'] >= 0.3 and (v['af'] is None or v['af'] < AF_DOMINANT)]
        homs = [v for v in rec if v['gt'] in HOM]
        hets = sorted([v for v in rec if v['gt'] in HET], key=lambda v: -v['impact'])
        cand = []
        if homs:
            best = max(homs, key=lambda v: v['impact'])
            cand.append(('homozygous', [best], best['impact']))
        pair, rel = None, 'unknown'
        for i, a in enumerate(hets):
            for b in hets[i + 1:]:
                r = phase_relation(a, b, qual)
                if r == 'cis':          # same haplotype: one allele, not a pair
                    continue
                pair, rel = (a, b), r
                break
            if pair:
                break
        if pair:
            base = (pair[0]['impact'] + pair[1]['impact']) / 2
            # trans is proven by read-level phasing; unknown phase is not penalised
            cand.append(('compound_het' + ('(trans)' if rel == 'trans' else ''),
                         list(pair), base * (1.25 if rel == 'trans' else 1.0)))
        dhets = [v for v in dom if v['gt'] in HET]
        if dhets:
            best = max(dhets, key=lambda v: v['impact'])
            cand.append(('dominant', [best], best['impact'] * 0.5))
        for model, vlist, vscore in cand:
            findings.append((pheno[g] * vscore, pheno[g], vscore, g, model, vlist))
    print('variants failing the read-level quality gate: %d' % n_dropped)

    findings.sort(reverse=True)
    pickle.dump(findings, open(os.path.join(WORK, 'blind_findings.pkl'), 'wb'))

    print('\nBLIND RANKED CANDIDATES (phenotype score x variant evidence)')
    print('%-4s %-9s %-9s %6s %6s  %s' % ('#', 'gene', 'model', 'combo', 'pheno', 'variants'))
    for i, (combo, ph, vs_, g, model, vlist) in enumerate(findings[:15], 1):
        desc = ' + '.join('%s:%d %s>%s %s [%s]' % (
            v['chrom'], v['pos'], v['ref'][:4], v['alt'][:4],
            v['hgvsp'][:14], 'novel' if v['af'] is None else '%.1e' % v['af']) for v in vlist)
        print('%-4d %-9s %-9s %6.3f %6.3f  %s' % (i, g, model, combo, ph, desc))


if __name__ == '__main__':
    main()
