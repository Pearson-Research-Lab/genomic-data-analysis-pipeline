# Quality Control (QC)

This directory contains scripts and documentation for genotype quality control (QC) prior to imputation and downstream GWAS analyses.

---

## Overview

Ancestry assignment and sample stratification were performed upstream and are **not part of this repository**.

This QC workflow starts from **per-chromosome bgzipped VCF files (`.vcf.gz`) with harmonised variant IDs** and processes two datasets separately:

1. **European-only dataset**
2. **Mixed-ancestry dataset** (including European samples)

Both datasets are processed independently through QC and imputation.

---

## Analysis strategy

The mixed-ancestry dataset was used for the **primary GWAS analysis**.

The European-only dataset was analysed separately as a **sensitivity analysis**, to assess the robustness of results to population structure and ancestry-specific effects.

---

## Pre-imputation QC

QC was performed using PLINK2, starting from indexed VCF files.

### Input requirements

- Per-chromosome VCF files (`chr1–22.vcf.gz`)
- Indexed with tabix-compatible `.tbi` files  
  (e.g. generated using `bcftools index -t`)
- Variant IDs harmonised upstream (e.g. `chr:pos:ref:alt` format)
- Biallelic SNPs only
- Autosomes only (chr1–22)

---

## QC filters applied

### European-only dataset

- **Hardy-Weinberg equilibrium**  
  Exclude variants with:  
  `p < 1e-6`

- **Minor allele frequency (MAF)**  
  Exclude variants with:  
  `MAF < 0.01`

- **Sample missingness**  
  Exclude samples with:  
  `missingness > 0.02`

- **Variant missingness**  
  Exclude variants with:  
  `missingness > 0.03`

---

### Mixed-ancestry dataset

- **Minor allele frequency (MAF)**  
  Exclude variants with:  
  `MAF < 0.01`

- **Sample missingness**  
  Exclude samples with:  
  `missingness > 0.02`

- **Variant missingness**  
  Exclude variants with:  
  `missingness > 0.03`

- **No Hardy-Weinberg equilibrium filtering applied**

---

## Rationale

Hardy-Weinberg equilibrium filtering was applied **only to the European-only dataset**.

This avoids inappropriate exclusion of variants in the mixed-ancestry dataset, where population structure can lead to deviations from HWE that are not due to genotyping error.

---

## Imputation

The filtered datasets were uploaded separately to the **Michigan Imputation Server** using:

- **Genome build:** GRCh37  
- **Reference panel:** 1000 Genomes Phase 3 v5  

European-only and mixed-ancestry datasets were imputed independently.

---

## Post-imputation QC

Post-imputation QC was performed separately for the European-only and mixed-ancestry datasets.

Imputed genotype data were filtered prior to GWAS analysis.

### Filters applied

- **Imputation quality (R2)**  
  Exclude variants with:  
  `R2 < 0.5`

- **Minor allele frequency (MAF)**  
  Exclude variants with:  
  `MAF < 0.01`

- **Variant type**  
  Retain only biallelic variants

- **Monomorphic variants**  
  Remove monomorphic variants

- **Variant identifiers**  
  Harmonise variant IDs to:  
  `chr:pos:ref:alt`

---

### Notes

- No Hardy-Weinberg equilibrium filtering was applied post-imputation
- Post-imputation QC thresholds were based on the study analysis plan and collaborator agreement

---

## Output

The QC pipeline produces:

- Filtered per-chromosome VCF files (`.vcf.gz`)
- Corresponding tabix index files (`.tbi`)
- Summary logs of variant counts before and after filtering

---

## Pipeline implementation

- Pre-imputation QC:  
  [`run_preimputation_qc.sh`](run_preimputation_qc.sh)

- Post-imputation QC:  
  [`run_postimputation_qc.sh`](run_postimputation_qc.sh)
