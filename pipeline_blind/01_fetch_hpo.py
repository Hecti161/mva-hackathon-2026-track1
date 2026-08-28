#!/usr/bin/env python3
"""
BLIND step 0 - Fetch the Human Phenotype Ontology and its gene annotations.

These are the only knowledge inputs the blind pipeline uses. No gene list, no
disease name, no candidate panel.

Output: work/hp.obo, work/genes_to_phenotype.txt
"""
import os
import urllib.request

WORK = os.environ.get('WORK_DIR', 'work')
URLS = {
    'genes_to_phenotype.txt': 'https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt',
    'hp.obo': 'https://purl.obolibrary.org/obo/hp.obo',
}


def main():
    os.makedirs(WORK, exist_ok=True)
    for name, url in URLS.items():
        path = os.path.join(WORK, name)
        if os.path.exists(path) and os.path.getsize(path) > 0:
            print('%-24s already present' % name)
            continue
        req = urllib.request.Request(url, headers={'User-Agent': 'hpo-fetch'})
        data = urllib.request.urlopen(req, timeout=300).read()
        open(path, 'wb').write(data)
        print('%-24s %6.1f MB' % (name, len(data) / 1e6))


if __name__ == '__main__':
    main()
