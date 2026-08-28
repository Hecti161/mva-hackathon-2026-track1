# MVA Hackathon 2026 — Track 1 (Variant Prediction)
## Methods description and report — Model 2: blind phenotype-driven pipeline

**Team name:** Hecti161
**Model number:** 2 of up to 6
**Proband:** `PROBAND01` (challenge identifier; the sequencing sample ID inside the VCF is `WGS_EX2312012`)
**Submission file:** `Hecti161_blind-hpo-pipeline.csv`
**Assembly:** GRCh38
**Code:** https://github.com/Hecti161/mva-hackathon-2026-track1 (`pipeline_blind/`)

---

## 1. What this model is for

Model 1 identified the causal `BUB1B` compound heterozygote by starting from a
20-gene panel that a human selected *because they already recognised the syndrome*.
That is a legitimate way to solve one case and a poor way to build a method, and we
declared it as the principal limitation of model 1.

**Model 2 removes the human from the loop and asks whether the answer survives.**

The only input is the eight HPO terms from the clinical document. No disease name, no
gene panel, no candidate list, no mention of Mosaic Variegated Aneuploidy or of
`BUB1B` anywhere in the code or its inputs.

**Result: the same compound heterozygote is recovered at rank 1.**

| Blind rank | Gene | Model | Score | Variants |
|---|---|---|---|---|
| **1** | **`BUB1B`** | **compound_het** | **1.293** | **chr15:40,209,701 T>G `p.Leu737Ter` + chr15:40,220,612 T>G `p.Asn1002Lys`** |
| 2 | `BUB1B` | dominant | 0.924 | chr15:40,209,701 T>G `p.Leu737Ter` |
| 3 | `RAI1` | dominant | 0.710 | chr17:17,793,784 AGC>A `p.Gln280AlafsTer` |

Only three findings survive the full chain, and the top two are the same gene.

---

## 2. The funnel, with numbers at every stage

| Stage | In | Out |
|---|---|---|
| Phenotype similarity over the HPO corpus | 8 HPO terms | 5,268 genes ranked — `BUB1B` **15th**, `CEP57` 14th, `TRIP13` 62nd |
| Top *K*=200 genes → Ensembl GRCh38 coordinates | 200 | 199 resolved |
| Streaming extraction from the proband VCF | 5,012,204 variants | 32,309 |
| Exonic / splice-region restriction (±10 bp) | 32,309 | 884 |
| Read-level quality gate | 884 | 43 rejected |
| Rarity + impact + inheritance modelling | — | **3 findings** |

*K* was fixed at 200 — the top 3.8% of the ranked corpus — before any variant was
examined.

**The two other MVA genes rank adjacent to `BUB1B`.** `CEP57` (MVA2) is 14th and
`TRIP13` (MVA3) 62nd, with `CEP57` and `BUB1B` scoring identically at 1.848. This is
the expected and correct behaviour: phenotype alone identifies the *disease*, and it
cannot distinguish which gene in that disease family is responsible. The genomic
evidence makes that call. Genes ranking above `BUB1B` are legitimate confounders that
also cause childhood rhabdomyosarcoma — `PTCH1` (Gorlin syndrome) and `HRAS`
(Costello syndrome) among them — not noise.

### Phenotype similarity

Resnik similarity with best-match average, as used by Phenomizer and Exomiser:

    score(gene) = mean_i [ max_j IC(MICA(q_i, g_j)) ]
    IC(t)       = -log( genes annotated with t or a descendant / all genes )

A gene scores highly only by explaining *several* of the proband's features, and rare
features count for more than common ones. Rhabdomyosarcoma (HP:0002859) carries an
information content of 5.52 against `BUB1B`, the single largest contribution to its
rank.

### Inheritance modelling

Three models are evaluated per gene: homozygous, compound heterozygous, and dominant.
Recessive models admit gnomAD AF < 1%; the dominant model requires AF < 0.1%, since a
dominant childhood-onset allele cannot be common. Findings are ranked by
`phenotype_score × variant_evidence`.

---

## 3. Two methodological contributions

### Physical phasing from PGT/PID rejects cis pairs automatically

The naive pipeline placed a `RAI1` compound heterozygote above `BUB1B`: two novel
frameshifts, `p.Gln280AlafsTer` and `p.Gln280HisfsTer`, 3 bp apart, both PASS, both
at GQ 99 and DP > 45. It is exactly the pattern a compound-heterozygote filter is
built to find, and it is wrong.

**The VCF already contained the disproof.** Both calls carry

    PGT:PID = 0|1 : 17793784_AGC_A

GATK's HaplotypeCaller emits `PID` (phase-set identifier) and `PGT` (phased genotype)
whenever two variants are close enough to be observed on the same reads. An identical
`PID` with an identical `PGT` places both variants **on the same haplotype** — one
complex indel emitted as two records, not two alleles in trans.

The pipeline now reads phase directly:

| Relationship | Evidence | Action |
|---|---|---|
| **cis** | shared `PID`, same `PGT` | reject the pair — it is one allele |
| **trans** | shared `PID`, opposite `PGT` | accept and **upweight** — compound het proven at read level |
| **unknown** | no shared `PID` | accept without penalty — too far apart to phase |

This is read-level evidence from the caller, not an inference from coordinates, and it
works in both directions: it kills false pairs *and* confirms real ones. `BUB1B`'s two
alleles fall in the third category — 10.9 kb apart, unphaseable by short reads — which
is precisely why parental testing remains the outstanding experiment.

We are not aware of another submission using `PGT`/`PID` this way. The rule is
generic: it applies to any recessive candidate in any short-read callset.

### What did *not* work, reported because it is useful

Three widely-recommended remedies were tested against these two false positives and
**none of them would have helped**:

1. **Broad Institute hard filters** (`QD < 2.0`, `FS > 60.0`, `MQ < 40.0`,
   `SOR > 3.0`) were already applied by the data provider — the VCF header records the
   exact `VariantFiltration` invocation — and we already require `FILTER == PASS`.
   Re-applying them is a no-op. Measured directly: `RAI1` has QD 12.14 and 12.68,
   FS 4.28 and 2.65, SOR 1.40 and 1.18, MQ 60. The low-depth `PEX5` call has QD 41.67.
   Both pass every Broad threshold comfortably.
2. **Allele-fraction filtering of homozygous calls** would not have caught the `PEX5`
   artefact: at AD 0,10 its VAF is 1.00. What exposes it is sample-level depth, DP 10.
3. **VCF normalisation** (`bcftools norm`, `vt`) left-aligns indels and splits or joins
   multiallelic records *at the same position*. It does not merge adjacent variants into
   a single MNP, so it would not have collapsed the `RAI1` pair. GATK's
   `--max-mnp-distance` does, but it is a HaplotypeCaller argument that requires
   re-calling from alignments.

What did work was the sample-level depth and genotype-quality gate (DP ≥ 15, GQ ≥ 50,
MQ ≥ 50), an allele-balance requirement (0.20–0.80 for heterozygous calls, ≥ 0.90 for
homozygous), and the phasing rule above.

---

## 4. Methods description form

**Automated output or manually curated?** This submission is **the automated output**
of the pipeline, which is what distinguishes it from model 1. Row 1 of the CSV is the
pipeline's rank-1 finding verbatim. Human input is confined to the `notes` column and
to the decision to omit the pipeline's rank-2 row, which restates a variant already
present in row 1 under a different inheritance model.

**Compound heterozygous output?** Yes, natively, and with phase evaluated from
read-level evidence rather than assumed — see section 3.

**Secondary and incidental findings.** The pipeline's rank-3 finding, `RAI1`
`p.Gln280AlafsTer`, is reported at EPCR 0.08 as `secondary` and annotated in the
submission as **not contributory**, with the reason: it is the surviving single
representation of the complex indel described above, it sits in a repetitive tract,
and the proband's phenotype does not match RAI1 haploinsufficiency. It is included so
that the automated output is visible in full rather than silently trimmed.

**Data provenance. Public data only; no proprietary data.** Human Phenotype Ontology
(`hp.obo`) and its gene annotations (`genes_to_phenotype.txt`); Ensembl REST for
coordinates, exon structures, MANE Select transcripts, HGVS and VEP consequences;
gnomAD v4 exome and genome frequencies, SIFT, PolyPhen-2 and AlphaMissense as served
through VEP. The only restricted input is the challenge VCF.

**Runtime and cost.**

| Stage | Wall time |
|---|---|
| HPO download | ~1 min |
| Phenotype ranking of 5,268 genes | ~3 min |
| VCF extraction over 199 loci | ~3 min |
| Exon structures + exonic filter | ~2 min |
| VEP annotation of 884 variants | ~2 min |
| **Total** | **~11 min, $0** |

Python standard library only. No reference genome, no `bcftools`/`tabix`/`samtools`,
no cluster, no GPU, no FASTQ download. Runs unmodified on Windows.

---

## 5. Strengths and limitations

**Strengths.** The result is reached without human disease knowledge, which is the
property that makes a method reusable on a case nobody has recognised yet. Every stage
reports its numbers, so the funnel is auditable rather than asserted. Phase is taken
from read-level evidence instead of assumed. Negative methodological results are
reported alongside positive ones. It runs in eleven minutes at no cost.

**Limitations.**

1. **It can only find genes the HPO already associates with the phenotype.** This is
   the fundamental ceiling of any phenotype-driven approach: a genuinely novel disease
   gene has no annotations, scores zero, and is invisible. This pipeline accelerates
   diagnosis of *recognised* disease; it does not discover new disease genes.
2. **The result depends on the true gene ranking within the top *K*.** `BUB1B` at 15th
   left ample margin at *K* = 200, but a case whose phenotype is less specific — or
   less completely recorded — could push the causal gene below any practical cutoff.
   The eight HPO terms here are unusually informative.
3. **The exonic restriction discards deep-intronic and regulatory variation.** It is
   what makes annotating 32,309 candidates tractable, and it is immaterial here because
   both causal alleles are coding. A case whose second allele were deep-intronic would
   be missed.
4. **The quality thresholds were formalised after observing the `PEX5` and `RAI1`
   failures.** The values themselves are conventional rather than tuned to the answer,
   but a pre-registered analysis would have fixed them in advance. Stated plainly
   because it is the kind of thing that is easy to leave out.
5. **Phase between the two `BUB1B` alleles remains unproven** — 10.9 kb apart, no
   shared `PID`, no parental data. Parental Sanger sequencing of the two positions
   remains the single most informative next experiment.

---

## 6. Method abstract (up to 500 words)

Model 1 of this submission identified `BUB1B` `c.2210T>G` `p.Leu737Ter` and `c.3006T>G`
`p.Asn1002Lys` as a causal compound heterozygote, but did so from a gene panel chosen
by an analyst who had already recognised the syndrome from the clinical description.
Model 2 tests whether the answer survives when that knowledge is removed.

The only input is the proband's eight HPO terms. Gene annotations from the Human
Phenotype Ontology are scored by Resnik similarity with best-match average — the
measure used by Phenomizer and Exomiser — weighting each matched term by its
information content, so that a gene ranks highly only by explaining several of the
proband's features and rare features count for more than common ones. This ranks 5,268
annotated disease genes, placing `BUB1B` 15th, with `CEP57` 14th and `TRIP13` 62nd.
That the three MVA genes cluster at the top while remaining mutually indistinguishable
is the correct behaviour: phenotype identifies the disease, genomics identifies the
gene.

The top 200 genes were resolved to GRCh38 coordinates through the Ensembl REST API and
their loci extracted from the 5,012,204-variant VCF in one streaming pass, yielding
32,309 variants. Canonical-transcript exon structures reduced this to 884 exonic or
splice-region calls, annotated through Ensembl VEP. A read-level quality gate
(FILTER PASS, DP ≥ 15, GQ ≥ 50, MQ ≥ 50, heterozygous allele balance 0.20–0.80,
homozygous VAF ≥ 0.90) rejected 43. Three inheritance models were then evaluated per
gene and findings ranked by phenotype score times variant evidence. `BUB1B` compound
heterozygous emerged at rank 1, recovering both causal alleles.

The decisive component proved to be phase. A naive run ranked a `RAI1` compound
heterozygote first: two novel frameshifts 3 bp apart, both PASS, GQ 99, DP > 45 — a
textbook false positive. The VCF contained its own disproof in GATK's `PGT`/`PID`
fields, which carry physical phase whenever variants are close enough to share reads.
Both `RAI1` calls report `PID 17793784_AGC_A` with identical `PGT 0|1`, placing them in
cis: one complex indel emitted as two records. The pipeline now rejects cis pairs and
upweights proven-trans pairs from these fields directly. The rule is generic to any
short-read callset. Notably, three standard remedies would not have helped, and this is
reported: Broad hard filters were already applied by the data provider and both
artefacts pass them (`RAI1` QD 12.14/12.68); allele-fraction filtering misses the
low-depth `PEX5` call at VAF 1.00; and `bcftools norm` does not merge adjacent variants
into MNPs.

The approach runs in eleven minutes at zero cost on Python's standard library, with no
reference genome, no bioinformatics toolchain and no FASTQ download. Its principal
limitation is structural: a phenotype-driven method can only find genes the ontology
already associates with the presentation, so it accelerates diagnosis of recognised
disease but cannot discover novel disease genes. Phase between the two `BUB1B` alleles
remains unproven at 10.9 kb separation.

---

*Released under CC BY 4.0, per challenge terms.*
