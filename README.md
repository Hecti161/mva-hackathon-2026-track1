# MVA Hackathon 2026 — Track 1 (Variant Prediction)

**Team:** Hecti161
**Proband:** `PROBAND01` (challenge identifier; the sequencing sample ID inside the VCF is `WGS_EX2312012`)
**Assembly:** GRCh38

A phenotype-first pipeline that identifies the causal variants in an undiagnosed child with
suspected Mosaic Variegated Aneuploidy, and independently tests the *cellular* phenotype of
the syndrome from the VCF alone.

Runs in **under 20 minutes on a laptop**, at **zero cost**, with no reference genome, no
cluster, no licensed software and no FASTQ download.

---

## Result

The proband is **compound heterozygous in `BUB1B`** — Mosaic Variegated Aneuploidy
syndrome type 1 (MVA1, [MIM 257300](https://omim.org/entry/257300), autosomal recessive).

| | Allele 1 | Allele 2 |
|---|---|---|
| GRCh38 | chr15:40,209,701 T>G | chr15:40,220,612 T>G |
| HGVS (NM_001211.6, MANE Select) | `c.2210T>G` `p.Leu737Ter` | `c.3006T>G` `p.Asn1002Lys` |
| Consequence | stop_gained, exon 17/23 | missense, exon 23/23 |
| gnomAD | 7.87e-05 (rs759242053) | **absent** |
| Genotype | 0/1 (DP 46, AD 21/25, GQ 99) | 0/1 (DP 28, AD 15/13, GQ 99) |
| ACMG | **Pathogenic** | **VUS, favour Likely Pathogenic** |

`p.Leu737Ter` is a premature termination codon in exon 17 of 23, predicted to trigger
nonsense-mediated decay — a **null allele**. `p.Asn1002Lys` sits in the final exon, escapes
NMD, falls inside the BubR1 **kinase domain**, and is deleterious by SIFT, PolyPhen-2 and
AlphaMissense alike. Biallelic null `BUB1B` is not survivable; **one null plus one
hypomorph** is the allelic architecture MVA1 requires, and it is what this proband carries.

Across the entire `BUB1B` locus ±25 kb, **exactly two variants are coding**, and they are
these two. `TP53` and `DICER1` — the main alternative causes of childhood rhabdomyosarcoma —
were screened and are negative.

Full reasoning, ACMG criteria and limitations: [`results/Hecti161_track1_report.md`](results/Hecti161_track1_report.md).

---

## What is interesting here beyond the answer

**MVA is a cellular phenotype, and we measured it.** The syndrome is defined by mosaic
aneuploidy across multiple chromosomes, not only by a genotype. That is measurable from the
VCF alone, with no alignments: if a chromosome is trisomic in a fraction *f* of cells,
heterozygous SNVs on it split their allele balance away from 0.5 toward (1+f)/(2+f) and
1/(2+f), inflating the **variance** of allele balance while leaving its **median** at 0.5.

Steps 5–7 implement that test across 2.2 M high-confidence heterozygous sites.

**The result is a negative, and it is reported as one.** No mosaic aneuploidy is detectable
above **~15% of cells**. Observed variance exceeds the binomial expectation by a uniform
~1.24× on every autosome — technical overdispersion, flat across the genome — and the
residual scatter tracks GC-richness and gene density (`chr19` > `chr16` ≈ `chr17` > `chr9`),
the known technical gradient. No chromosome departs from the trend.

**A chromosome-level summary alone would have produced a false positive.** Chromosomes 17,
20 and 22 initially appeared elevated. Binning at 10 Mb resolution localised the entire
signal to centromeric windows — chr1 120–140 Mb, chr9 40–60 Mb across the 9q12
heterochromatin block, chr20 20–40 Mb, chr17 20–30 Mb, chr22 0–20 Mb — with flat 6.5–8.5%
dispersion everywhere else. The small chromosomes had merely looked elevated because a
larger fraction of their length is centromeric repeat. Step 6 exists to catch exactly this,
and the windowed output is printed so a reader can check it rather than trust it.

The negative does not weaken the diagnosis: MVA aneuploidy is scored by karyotype on
*cultured* cells, and aneuploid cells are depleted in vivo and further selected against
during library preparation. A bulk blood WGS at 44× is an insensitive assay for them. It is
worth measuring anyway, because a *positive* would have been strong independent
confirmation.

---

## Model 2 — the blind pipeline

The targeted pipeline above answers one case. It does not generalise, because a human
chose its gene panel after recognising the syndrome. Model 2 removes that human and asks
whether the answer survives.

**The only input is the proband's eight HPO terms.** No disease name, no gene panel, no
candidate list. The strings "MVA" and "BUB1B" appear nowhere in `pipeline_blind/` or its
inputs.

**It recovers the same compound heterozygote at rank 1.**

| Stage | In | Out |
|---|---|---|
| Phenotype similarity over the HPO corpus | 8 HPO terms | 5,268 genes ranked — `BUB1B` **15th** |
| Top *K*=200 genes → GRCh38 coordinates | 200 | 199 resolved |
| Streaming extraction from the VCF | 5,012,204 variants | 32,309 |
| Exonic / splice-region restriction | 32,309 | 884 |
| Read-level quality gate | 884 | 43 rejected |
| Rarity + impact + inheritance modelling | — | **3 findings, `BUB1B` compound het at #1** |

*K* was fixed at 200 before any variant was examined. Gene ranking uses Resnik similarity
with best-match average, weighting each matched HPO term by its information content.
`CEP57` (MVA2) ranks 14th and `TRIP13` (MVA3) 62nd — phenotype identifies the *disease*
and cannot separate genes within it; the genomic evidence makes that call.

### Physical phasing from PGT/PID

A naive run ranked a `RAI1` compound heterozygote above `BUB1B`: two novel frameshifts
3 bp apart, both PASS, GQ 99, DP > 45. A textbook false positive — and the VCF contained
its own disproof. Both calls carry `PID 17793784_AGC_A` with identical `PGT 0|1`.

GATK emits `PID` (phase-set id) and `PGT` (phased genotype) whenever two variants are
close enough to be seen on the same reads. Identical `PID` **and** identical `PGT` places
them on one haplotype — a single complex indel emitted as two records, not two alleles in
trans.

| Relationship | Evidence | Action |
|---|---|---|
| **cis** | shared `PID`, same `PGT` | reject — one allele, not a pair |
| **trans** | shared `PID`, opposite `PGT` | accept and upweight — proven at read level |
| **unknown** | no shared `PID` | accept without penalty — unphaseable |

Read-level evidence from the caller, not an inference from coordinates, and it works both
ways: it kills false pairs and confirms real ones. Generic to any short-read callset.

Three standard remedies were tested against these artefacts and **none would have helped**,
which is worth recording: Broad hard filters were already applied by the data provider (the
VCF header carries the `VariantFiltration` invocation) and both artefacts pass them
comfortably — `RAI1` scores QD 12.14/12.68, FS 4.28/2.65, SOR 1.40/1.18, MQ 60; allele-
fraction filtering of homozygous calls misses the low-depth `PEX5` artefact, whose VAF is
1.00 at DP 10; and `bcftools norm` left-aligns and splits or joins multiallelics at one
position but does not merge adjacent variants into an MNP.

Full write-up: [`results/Hecti161_track1_report_model2.md`](results/Hecti161_track1_report_model2.md).

---

## Honest limitations

Stated plainly, because the panel cannot assess what it cannot see.

1. ~~**The gene panel encodes human prior knowledge.**~~ **Resolved by model 2.** The
   targeted pipeline selects ~20 genes *because the analyst already knew* that
   rhabdomyosarcoma + IUGR + parental miscarriage implies MVA implies `BUB1B`. The blind
   pipeline in [`pipeline_blind/`](pipeline_blind/) removes that knowledge entirely and
   recovers the same compound heterozygote at rank 1 — see below.
2. **Phase is inferred, not proven.** The two alleles are 10.9 kb apart — beyond short-read
   phasing — and no parental data are provided. Parental Sanger sequencing of the two
   positions is the single most informative next experiment: it would activate ACMG PM3,
   upgrade `p.Asn1002Lys` to Likely Pathogenic, and establish 25% recurrence risk.
3. **Secondary findings are not exhaustive** — 20 genes screened, not the full ACMG SF v3.2
   panel (~80 genes).
4. **Small variants only.** No SV, CNV or mobile-element calling. Immaterial here, since two
   coding alleles fully explain the phenotype, but required for a case with only one hit.
5. **Non-coding variation annotated but not modelled.** No SpliceAI; unnecessary once two
   coding alleles were identified.

---

## Reproducing

### Data access

The input data are **controlled access** and are **not** in this repository. Request access
at [`SageBio/mva-hackathon-2026-data`](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data),
then:

```bash
hf download SageBio/mva-hackathon-2026-data \
  WGS_EX2312012_HGWCNDSX7.vcf.gz \
  --repo-type dataset --local-dir ./data
```

Only the VCF is needed (~300 MB). The FASTQ files are not required by this pipeline.

### Requirements

Python 3.8+ and the standard library. No third-party packages, no bioinformatics toolchain.
`bcftools`/`tabix` are deliberately not used, so the pipeline runs unmodified on Windows.

Network access is required: gene coordinates and variant annotations are fetched live from
the [Ensembl REST API](https://rest.ensembl.org).

### Run

```bash
export WORK_DIR=work           # optional, defaults to ./work
VCF=data/WGS_EX2312012_HGWCNDSX7.vcf.gz

python pipeline/01_fetch_gene_coords.py                       # ~10 s
python pipeline/02_extract_regions.py       $VCF              # ~3 min
python pipeline/03_annotate_vep.py          BUB1B,CEP57,TRIP13,BUB1,CENPE
python pipeline/04_prioritize.py            work/vep_BUB1B.json
python pipeline/03_annotate_vep.py          TP53,DICER1,ATM,NBN,MRE11,MAD1L1,MAD2L1,BUB3,KNL1,TTK,CENPF,ZWILCH,CDC20,TUBB,PLK4
python pipeline/04_prioritize.py            work/vep_TP53.json

python pipeline/05_mosaic_collect.py        $VCF              # ~5 min
python pipeline/06_mosaic_windows.py
python pipeline/07_mosaic_variance_test.py
```

Model 2, the blind pipeline — takes no gene list, only the HPO terms hard-coded in step 02:

```bash
python pipeline_blind/01_fetch_hpo.py                         # ~1 min
python pipeline_blind/02_rank_genes_by_phenotype.py           # ~3 min
python pipeline_blind/03_extract_topgenes.py $VCF             # ~3 min
python pipeline_blind/04_filter_exonic.py                     # ~2 min
python pipeline_blind/05_rank_candidates.py                   # ~2 min
```

Step 4 on `vep_BUB1B.json` prints the two causal variants at the top of the
HIGH/MODERATE table. Contig naming is normalised internally, so a VCF using either `1` or
`chr1` works unchanged — note that the challenge VCF uses `1` while the submission format
requires `chr15`.

### Everything under `work/` is derived patient data

`work/` holds proband genotypes and is excluded by [`.gitignore`](.gitignore). Do not commit
or redistribute it. The same applies to the VCF, the FASTQ files and the clinical phenotype
document, all of which are covered by the challenge's controlled-access terms.

---

## Contents

```
pipeline/
  01_fetch_gene_coords.py      candidate loci from Ensembl REST (not hard-coded)
  02_extract_regions.py        single streaming pass over the VCF
  03_annotate_vep.py           Ensembl VEP REST: consequence, MANE, gnomAD, SIFT/PolyPhen/AlphaMissense
  04_prioritize.py             impact + frequency + genotype ranking, read-level quality
  05_mosaic_collect.py         genome-wide allele balance harvest
  06_mosaic_windows.py         per-chromosome and 10 Mb windowed dispersion
  07_mosaic_variance_test.py   variance test vs binomial, detection limit
pipeline_blind/                model 2 - no gene panel, HPO terms only
  01_fetch_hpo.py              Human Phenotype Ontology + gene annotations
  02_rank_genes_by_phenotype.py  Resnik/best-match-average over 5,268 disease genes
  03_extract_topgenes.py       top-200 loci from Ensembl, one pass over the VCF
  04_filter_exonic.py          canonical-transcript exons, +/-10 bp splice pad
  05_rank_candidates.py        VEP, quality gate, PGT/PID phasing, inheritance models
results/
  Hecti161_bub1b-compound-het.csv     model 1 submission (targeted panel)
  Hecti161_track1_report.md           model 1 report
  Hecti161_blind-hpo-pipeline.csv     model 2 submission (blind, HPO-driven)
  Hecti161_track1_report_model2.md    model 2 report
```

---

## Licence

Code, report and submission are released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), per the challenge terms, so that
the methods here can be reused for other undiagnosed individuals.

No patient data are included in this repository.
