# MVA Hackathon 2026 — Track 1 (Variant Prediction)
## Methods description and report

**Team name:** Hecti161
**Model number:** 1 of up to 6
**Proband:** `PROBAND01` (challenge identifier; the sequencing sample ID inside the VCF is `WGS_EX2312012`)
**Submission file:** `Hecti161_bub1b-compound-het.csv`
**Assembly:** GRCh38
**Code:** https://github.com/Hecti161/mva-hackathon-2026-track1 (`pipeline/` — seven scripts, end-to-end reproducible)

---

## 1. Result summary

The proband is **compound heterozygous for two variants in `BUB1B`**, consistent with
**Mosaic Variegated Aneuploidy syndrome type 1** (MVA1; MIM 257300; autosomal recessive).

| | Allele 1 | Allele 2 |
|---|---|---|
| **GRCh38** | chr15:40,209,701 T>G | chr15:40,220,612 T>G |
| **HGVS** (NM_001211.6 / ENST00000287598, MANE Select) | `c.2210T>G` `p.Leu737Ter` | `c.3006T>G` `p.Asn1002Lys` |
| **Consequence** | stop_gained, exon 17/23 | missense, exon 23/23 |
| **dbSNP** | rs759242053 | — (novel) |
| **gnomAD** | exomes 7.87e-05 / genomes 3.29e-05 | **absent** |
| **Genotype** | 0/1 | 0/1 |
| **Depth / allelic depth** | DP 46, AD 21,25 (AB 0.54) | DP 28, AD 15,13 (AB 0.46) |
| **GQ / MQ / FILTER** | 99 / 60.0 / PASS | 99 / 60.0 / PASS |
| **In silico** | n/a (nonsense) | SIFT deleterious; PolyPhen-2 probably_damaging; AlphaMissense likely_pathogenic |

### Why this is the answer

**Allelic architecture matches the syndrome.** `p.Leu737Ter` introduces a premature
termination codon in exon 17 of 23 — well upstream of the penultimate exon-junction
complex — and is therefore predicted to trigger nonsense-mediated decay, producing a
**null allele**. `p.Asn1002Lys` falls in the **final exon (23/23)**, so the transcript
escapes NMD and is translated, and the substituted residue lies within the **BubR1 kinase
domain** (Pfam/CDD `cd14029`). Complete biallelic loss of `BUB1B` is not compatible with
live birth; surviving MVA1 patients characteristically carry **one null allele plus one
hypomorphic allele retaining partial function**. That is precisely the configuration
observed here.

**Phenotype is gene-specific, not merely compatible.** Rhabdomyosarcoma is the signature
malignancy of `BUB1B`-associated MVA rather than a generic cancer-predisposition feature.
Around it, the reported HPO cluster — intrauterine growth restriction (~1 kg at 32 weeks,
HP:0001518), short stature (HP:0004322), failure to thrive and skeletal muscle atrophy
(HP:0001508, HP:0003202), nephrocalcinosis (HP:0000121) and parental recurrent pregnancy
loss (HP:0200067) — reconstructs the classical MVA1 presentation. Parental recurrent
miscarriage is additionally consistent with both parents being heterozygous carriers.

**The gene region is otherwise clean.** Across the full `BUB1B` locus ±25 kb, 106 variants
were called; exactly **two are coding**, and they are the two above. There are no
splice-donor, splice-acceptor, frameshift or additional missense candidates competing for
the causal role.

**Competing explanations were tested and excluded.** No pathogenic variant was found in
`TP53` or `DICER1`, the two principal alternative causes of childhood rhabdomyosarcoma.
The only `TP53` coding call is the common benign polymorphism `p.Pro72Arg`
(rs1042522, AF 0.75). The rhabdomyosarcoma in this proband is therefore attributable to
`BUB1B`, not to Li-Fraumeni or DICER1 syndrome.

### Variant classification (ACMG/AMP 2015)

| Variant | Criteria applied | Classification |
|---|---|---|
| `p.Leu737Ter` | PVS1 (nonsense, NMD-predicted, LoF is the established mechanism) + PM2_Supporting + PP4 | **Pathogenic** |
| `p.Asn1002Lys` | PM1 (kinase domain) + PM2_Supporting (absent from gnomAD) + PP3 (concordant in silico) + PP4 + PM3_Supporting (in trans, *inferred*) | **VUS, favour Likely Pathogenic** |

`p.Asn1002Lys` is deliberately **not** called Likely Pathogenic outright, because the
evidence code that would move it there (PM3) depends on phase, which this dataset cannot
establish. See section 4.

---

## 2. Approach in detail

The pipeline is deliberately small, dependency-light and auditable. It runs on a
consumer laptop with no cluster, no reference genome download and no licensed software.

**Step 1 — Candidate locus definition.** Gene coordinates were retrieved programmatically
from the **Ensembl REST API** (`/lookup/symbol`) rather than hard-coded, so the panel is
declarative and cannot drift out of sync with the assembly. Two groups were queried: the
MVA / mitotic spindle assembly checkpoint genes (`BUB1B`, `CEP57`, `TRIP13`, `BUB1`,
`CENPE`, `BUB3`, `MAD1L1`, `MAD2L1`, `CDC20`, `PLK4`, `TUBB`, `KNL1`, `TTK`, `ZWILCH`,
`CENPF`) and the chromosomal-instability / childhood-cancer differential (`TP53`,
`DICER1`, `ATM`, `NBN`, `MRE11`). Each interval was padded by +/-25 kb to retain promoter,
UTR and near-regulatory space.

**Step 2 — Extraction.** A single streaming pass over the 5,012,204-variant VCF
(`pipeline/02_extract_regions.py`) emitted the 3,259 records intersecting those intervals. `bcftools`
and `tabix` are unavailable on the target platform (Windows), so the reader is pure Python
over `gzip` — slower than a tabix seek but portable and free of a toolchain install.

**Step 3 — Annotation.** Extracted variants were annotated via the **Ensembl VEP REST API**
in batches (`pipeline/03_annotate_vep.py`), requesting consequence terms, MANE Select transcript
mapping, HGVS c./p. nomenclature, exon numbering, protein domains, gnomAD exome and genome
allele frequencies, and SIFT / PolyPhen-2 / AlphaMissense predictions.

**Step 4 — Prioritisation.** Variants were ranked by (a) predicted impact on the MANE
transcript, (b) population frequency, with a recessive-appropriate threshold, and
(c) genotype consistency with a recessive model. Quality was then inspected at read level
(DP, allelic balance, GQ, MQ, FS, SOR) to exclude mapping and strand artefacts.

**Step 5 — Orthogonal test of the cellular phenotype.** See section 3.

### Was the submission automatically generated or manually curated?

**Both, and the distinction is stated explicitly.** The *identification and ranking* of the
two causal variants is the automated output of steps 1-4: they emerge as the only two
coding variants in the locus and the only HIGH/MODERATE-impact rare calls among the MVA
gene set. No manual variant hunting was involved.

Manual curation was applied to three things: (i) **assembly of the compound-heterozygous
pair** into a single submission row, (ii) **ACMG criterion assignment and the decision not
to over-call `p.Asn1002Lys`**, and (iii) the **written rationale** in the `notes` column.
The gene panel in step 1 also embeds human prior knowledge — see the limitations in
section 5, where this is treated as the principal weakness of the approach rather than
glossed over.

### Compound heterozygous output

**Yes.** The submission format's paired-variant encoding is used natively: row 1 carries
both alleles as a single compound-heterozygous prediction. Rows 2 and 3 re-state each
allele individually, at lower EPCR, purely as a defensive hedge against the ground truth
being encoded per-variant rather than per-pair. They are the same finding, not additional
candidates.

### Handling of secondary and incidental findings

Two variants are reported with `finding_type = secondary` and deliberately **very low
EPCR** (0.05 and 0.03): `ATM` `p.Ser978Pro` (rs139552233) and `MAD1L1` `p.Arg59Cys`
(rs121908982). Both are surfaced for transparency of the cancer-predisposition sweep and
both are annotated in the submission as **not contributory** — `ATM` at gnomAD AF 0.5% is
far too common for a pathogenic allele of that gene and is monoallelic; `MAD1L1` is
SIFT-tolerated, monoallelic, and at AF 0.57%.

The intent is the opposite of padding the list. Reporting them at near-zero confidence,
with the reasons they were rejected, documents that the differential was actually
examined. No variant was included to inflate recall.

### Data provenance

**Public data only. No proprietary data were used.** Sources: Ensembl REST API (gene
coordinates, VEP consequences, MANE Select transcripts, HGVS); gnomAD v4 exome and genome
allele frequencies as served through VEP; dbSNP identifiers; SIFT, PolyPhen-2 and
AlphaMissense predictions as served through VEP; Pfam/CDD domain assignments; HPO terms as
supplied in the challenge phenotype document; and published clinical literature on
`BUB1B`-associated MVA for phenotype-gene reasoning. The only restricted input is the
challenge dataset itself.

### Runtime and cost

| Stage | Wall time | Cost |
|---|---|---|
| VCF streaming extraction (5.0 M variants) | ~3 min | — |
| VEP REST annotation (3,259 variants) | ~6 min | — |
| Genome-wide allele-balance scan (section 3) | ~5 min | — |
| **Total** | **< 20 min on a consumer laptop** | **$0** |

No GPU, no cluster, no compute account, ~300 MB of input. The **FASTQ files were not
required and were not downloaded**: the causal alleles are cleanly called in the supplied
VCF, and the ~10.9 kb separation between them puts phasing out of reach of short reads
regardless. That saves the entire ~85 GB download for this analysis.

---

## 3. Orthogonal analysis: measuring the cellular phenotype from the VCF

MVA is defined by a *cellular* phenotype — mosaic aneuploidy across multiple chromosomes —
not only by a genotype. That phenotype is directly measurable from the VCF alone, without
alignments: if a chromosome is trisomic in a fraction *f* of cells, heterozygous SNVs on
that chromosome split their allele balance away from 0.5 toward (1+f)/(2+f) and 1/(2+f),
inflating the variance of allele balance without shifting its median.

2,226,786 high-quality heterozygous SNVs (PASS, biallelic, GQ >= 60, 15 <= DP <= 80; median
depth 44x) were extracted genome-wide and tested per chromosome against the binomial
expectation, after excluding centromeric and pericentromeric windows.

**Result: no mosaic aneuploidy detected, with a stated detection limit of ~15% of cells.**

Observed variance exceeds the binomial expectation by a uniform factor of ~1.24x on every
autosome. That excess is technical overdispersion — reference bias, PCR duplication,
sequencing error — and it is flat across the genome. Residual between-chromosome scatter
tracks GC-richness and gene density (`chr19` > `chr16` ~ `chr17` > `chr9` > `chr7`), the
known technical gradient; **no chromosome departs from that trend.**

An earlier, cruder pass appeared to flag chromosomes 17, 20 and 22. Binning at 10 Mb
resolution localised the entire signal to centromeric windows (chr1 120-140 Mb, chr9
40-60 Mb spanning the 9q12 heterochromatin block, chr20 20-40 Mb, chr17 20-30 Mb, chr22
0-20 Mb), with flat 6.5-8.5% dispersion everywhere else. The small chromosomes had merely
looked elevated because a larger fraction of their length is centromeric repeat. **This
was an artefact and is reported as one.**

**This negative result does not weaken the diagnosis, and the reason matters.** MVA
aneuploidy is classically scored by karyotype on *cultured* lymphocytes or fibroblasts.
Aneuploid cells are depleted in vivo and further selected against during culture and
library preparation, and a bulk blood WGS at 44x is a comparatively insensitive assay for
them. Absence of a detectable signal at >=15% cell fraction is therefore expected and
uninformative against `BUB1B` causality — but it is worth measuring and reporting, because
a *positive* signal would have been strong independent confirmation.

---

## 4. The single most informative next experiment

**Parental segregation testing by Sanger sequencing of the two `BUB1B` positions.**

It is inexpensive, it is the only outstanding piece of evidence, and it resolves the one
genuine weakness in this submission. The two variants lie 10.9 kb apart, which is beyond
short-read phasing, and no parental data are provided. *Trans* configuration is therefore
**inferred, not proven.** Confirming that each parent carries one allele would activate
PM3 and move `p.Asn1002Lys` from VUS to Likely Pathogenic, completing the biallelic
diagnosis formally. It would additionally establish 25% recurrence risk for genetic
counselling — directly relevant given the documented history of recurrent pregnancy loss.

The inference is nevertheless strong: a rare nonsense allele and a novel, computationally
deleterious kinase-domain missense allele in the same autosomal-recessive gene, in a
proband whose phenotype is specific to that gene, with the allelic architecture the
syndrome requires.

---

## 5. Strengths and limitations

**Strengths.** Sub-20-minute, zero-cost, laptop-scale runtime with no reference genome or
licensed toolchain. Gene coordinates and annotations are pulled live from Ensembl rather
than hard-coded, so the method does not silently decay. Read-level quality was verified
rather than assumed. The differential diagnosis was actively excluded, not ignored. A
second, orthogonal assay of the cellular phenotype was performed, and its result reported
honestly as negative with an explicit sensitivity bound. Uncertainty is stated where it
exists: `p.Asn1002Lys` is left as a VUS rather than over-called.

**Limitations.**

1. **The gene panel encodes human prior knowledge.** The strongest criticism of this
   submission is that step 1 selected ~20 genes *because the analyst already knew that
   rhabdomyosarcoma plus IUGR plus parental miscarriage implies MVA implies `BUB1B`*. The
   pipeline confirmed and characterised that hypothesis rigorously; it did not generate it
   from first principles. For a method intended to be reused on *other* undiagnosed
   individuals, this is the component that must be replaced by automated, HPO-driven,
   genome-wide gene prioritisation. It is the priority for a subsequent model.
2. **Secondary findings are not exhaustive.** Only 20 genes were screened. A complete
   incidental-findings analysis would require the full ACMG SF v3.2 panel (~80 genes).
   This does not affect the primary finding.
3. **Only small variants were considered.** Structural variants, copy-number changes and
   mobile-element insertions were not called. This is immaterial here, since two coding
   alleles fully explain the phenotype, but a CNV caller would be required for a case where
   only one hit were found.
4. **Phase is inferred, not proven** (section 4).
5. **Non-coding variation was annotated but not systematically modelled.** Twelve intronic
   `BUB1B` variants were annotated; none was in a splice-relevant position. No
   deep-intronic splicing predictor (e.g. SpliceAI) was applied, as it was unnecessary once
   two coding alleles were identified.

---

## 6. Method abstract (up to 500 words)

We identify the proband as compound heterozygous for `BUB1B` `c.2210T>G` `p.Leu737Ter`
(chr15:40,209,701 T>G) and `c.3006T>G` `p.Asn1002Lys` (chr15:40,220,612 T>G), establishing
Mosaic Variegated Aneuploidy syndrome type 1 (MIM 257300).

Our approach inverts the usual genome-wide filtering cascade. Rather than progressively
narrowing five million variants by frequency and impact, we treat the clinical phenotype as
a strong prior that collapses the search space before any variant is examined. The
co-occurrence of rhabdomyosarcoma, intrauterine growth restriction, short stature,
nephrocalcinosis and parental recurrent pregnancy loss is not a generic
cancer-predisposition pattern; it is the presentation of a mitotic spindle assembly
checkpoint disorder, and rhabdomyosarcoma specifically points to the `BUB1B` subtype. Gene
coordinates for the MVA/SAC gene set and for the childhood-cancer differential were pulled
live from the Ensembl REST API; a single streaming pass over the VCF extracted 3,259 of
5,012,204 variants; these were annotated through the Ensembl VEP REST API for consequence,
MANE transcript mapping, gnomAD frequency, protein domain and SIFT/PolyPhen-2/AlphaMissense
prediction.

Exactly two coding variants exist across the entire `BUB1B` locus, and they are the two
reported. Their configuration is the one MVA1 requires: a nonsense allele in exon 17 of 23,
predicted to undergo nonsense-mediated decay and therefore null, paired with a novel
last-exon missense allele that escapes NMD, sits within the BubR1 kinase domain and is
deleterious by three independent predictors. Complete biallelic `BUB1B` loss is not
survivable; a null plus a hypomorph is the canonical surviving genotype. Both calls are of
high quality (GQ 99, MQ 60, allelic balance 0.54 and 0.46, PASS). `TP53` and `DICER1`, the
main alternative causes of childhood rhabdomyosarcoma, were screened and are negative.

We additionally tested the *cellular* phenotype directly from the VCF, without alignments,
by measuring per-chromosome excess variance in heterozygous allele balance against the
binomial expectation across 2.2 million high-quality sites. No mosaic aneuploidy is
detectable above ~15% cell fraction; an apparent signal on chromosomes 17, 20 and 22
localised entirely to centromeric repeat and is reported as the artefact it is. We regard
this negative as expected — MVA aneuploidy is scored on cultured cells, and aneuploid cells
are selected against in blood and in library preparation — and report it with its
sensitivity bound rather than omitting it.

The method's principal strength is economy and auditability: under 20 minutes, $0, one
laptop, no reference genome, no licensed software, no FASTQ download. Its principal
limitation is that the gene panel encodes analyst prior knowledge; the pipeline confirmed
and characterised the hypothesis rigorously but did not generate it automatically.
Replacing that step with HPO-driven genome-wide prioritisation is the priority for
generalising this to other undiagnosed cases. Phase is inferred rather than proven, and
parental Sanger sequencing of the two positions is the single most informative next
experiment — it would upgrade `p.Asn1002Lys` from VUS to Likely Pathogenic and establish
recurrence risk.

---

*Released under CC BY 4.0, per challenge terms.*
