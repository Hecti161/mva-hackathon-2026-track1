#!/usr/bin/env python3
"""
BLIND step 1 - Phenotype-driven gene ranking.

Input:  the proband's 8 HPO terms. Nothing else. No gene list, no disease name,
        no mention of MVA or any candidate gene.
Output: every disease gene in the HPO corpus, ranked by how well it explains
        those 8 terms.

Similarity measure: Resnik with best-match average, as used by Phenomizer and
Exomiser. For each query term the most informative common ancestor (MICA) shared
with each of the gene's annotated terms is found, and the gene scores the mean
over query terms of its best match:

    score(gene) = mean_i [ max_j IC(MICA(q_i, g_j)) ]
    IC(t)       = -log( genes annotated with t or a descendant / all genes )

A gene therefore scores highly only by explaining *several* of the proband's
features, and rare features count for more than common ones.
"""
import collections
import math
import os
import pickle

WORK = os.environ.get('WORK_DIR', 'work')

# The proband's reported phenotype - the ONLY input to this step.
QUERY = {
    'HP:0002859': 'Rhabdomyosarcoma',
    'HP:0000121': 'Nephrocalcinosis',
    'HP:0004322': 'Short stature',
    'HP:0001508': 'Failure to thrive',
    'HP:0003202': 'Skeletal muscle atrophy',
    'HP:0001622': 'Premature birth',
    'HP:0001518': 'Small for gestational age',
    'HP:0200067': 'Recurrent spontaneous abortion',
}


def parse_obo(path):
    """term -> (name, direct parents). Obsolete terms dropped; alt_ids mapped."""
    parents, names, alt = {}, {}, {}
    cur, obsolete = None, False
    for line in open(path, encoding='utf-8'):
        line = line.rstrip('\n')
        if line == '[Term]':
            cur, obsolete = None, False
        elif line.startswith('id: HP:'):
            cur = line[4:]
            parents.setdefault(cur, set())
        elif cur and line.startswith('name: '):
            names[cur] = line[6:]
        elif cur and line.startswith('is_a: '):
            parents[cur].add(line[6:].split(' !')[0].strip())
        elif cur and line.startswith('alt_id: '):
            alt[line[8:].strip()] = cur
        elif cur and line.startswith('is_obsolete: true'):
            obsolete = True
            parents.pop(cur, None)
            names.pop(cur, None)
    return parents, names, alt


def ancestors_of(term, parents, cache):
    """Transitive closure of is_a, including the term itself."""
    if term in cache:
        return cache[term]
    seen, stack = {term}, [term]
    while stack:
        t = stack.pop()
        for p in parents.get(t, ()):
            if p not in seen:
                seen.add(p)
                stack.append(p)
    cache[term] = seen
    return seen


def main():
    parents, names, alt = parse_obo(os.path.join(WORK, 'hp.obo'))
    cache = {}
    print('ontology: %d terms' % len(parents))

    # gene -> annotated terms
    gene_terms = collections.defaultdict(set)
    with open(os.path.join(WORK, 'genes_to_phenotype.txt'), encoding='utf-8') as fh:
        next(fh)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 4:
                continue
            term = alt.get(f[2], f[2])
            if term in parents:
                gene_terms[f[1]].add(term)
    print('genes with annotations: %d' % len(gene_terms))

    # IC of every term: how many genes are annotated with it or any descendant
    gene_closure = {g: set().union(*(ancestors_of(t, parents, cache) for t in ts))
                    for g, ts in gene_terms.items()}
    freq = collections.Counter()
    for anc in gene_closure.values():
        freq.update(anc)
    total = len(gene_terms)
    ic = {t: -math.log(c / total) for t, c in freq.items() if c > 0}

    # For every query term, the IC of its MICA with every ontology term
    mica = {}
    for q in QUERY:
        qa = ancestors_of(q, parents, cache)
        row = {}
        for t in parents:
            common = qa & ancestors_of(t, parents, cache)
            if common:
                row[t] = max((ic.get(c, 0.0) for c in common), default=0.0)
        mica[q] = row

    # score each gene: mean over query terms of its best-matching annotation
    scored = []
    for g, ts in gene_terms.items():
        per_q = []
        for q in QUERY:
            row = mica[q]
            per_q.append(max((row.get(t, 0.0) for t in ts), default=0.0))
        scored.append((sum(per_q) / len(QUERY), g, per_q))
    scored.sort(reverse=True)

    pickle.dump([(s, g) for s, g, _ in scored],
                open(os.path.join(WORK, 'gene_ranking.pkl'), 'wb'))

    qkeys = list(QUERY)
    print('\nTop 30 genes explaining the proband phenotype')
    print('%-4s %-10s %7s   %s' % ('#', 'gene', 'score', 'per-term IC (' + ' '.join(
        QUERY[q][:11] for q in qkeys) + ')'))
    for i, (s, g, per_q) in enumerate(scored[:30], 1):
        print('%-4d %-10s %7.3f   %s' % (i, g, s, ' '.join('%5.2f' % v for v in per_q)))

    print('\nRank of every gene later found to carry a candidate variant:')
    ranks = {g: i for i, (_s, g, _p) in enumerate(scored, 1)}
    for g in ['BUB1B', 'CEP57', 'TRIP13', 'TP53', 'DICER1', 'ATM']:
        print('  %-8s rank %5s of %d' % (g, ranks.get(g, 'n/a'), len(scored)))


if __name__ == '__main__':
    main()
