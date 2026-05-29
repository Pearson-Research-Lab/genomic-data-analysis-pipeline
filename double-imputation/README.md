# Double Imputation Pipeline

This directory contains scripts and documentation for a two-phase double imputation workflow designed to maximise overlap between independently imputed genotype datasets.

When multiple cohorts or batches are imputed separately, differences in variant coverage, QC filtering, and imputation quality can substantially reduce the number of variants shared across all datasets. This pipeline identifies high-confidence common variants, harmonises datasets, and generates a scaffold for a second imputation phase to improve cross-dataset compatibility and downstream meta-analysis.

The common variant scaffold is intentionally restricted to very high-confidence variants (typically INFO/R² ≥ 0.99) to minimise propagation of imputation uncertainty into the second imputation phase.

## Overview

The pipeline performs:

1. Raw genotype cleaning
2. Variant ID harmonisation to `chr:pos:ref:alt`
3. Per-chromosome splitting
4. Biallelic SNP restriction
5. Pre-imputation QC
6. Will Rayner strand/reference checks
7. Michigan Imputation Server input preparation
8. Post-imputation QC
9. Cross-batch harmonisation and variant intersection
10. Common variant scaffold generation
11. Second-phase imputation input preparation

## QC thresholds

### Pre-imputation

- Variant call rate filter
- Sample genotyping rate filter
- HWE filter where appropriate
- Biallelic SNPs only
- Variant IDs standardised to `chr:pos:ref:alt`
- Strand/reference alignment checked using Will Rayner tools
- Thresholds may vary by project, ancestry composition, and study design.

### Post-imputation

- INFO/R2 threshold applied depending on stage:
  - Standard GWAS-ready data: `R2 >= 0.5`
  - Double-imputation scaffold/intersection: stricter, typically `R2 >= 0.99`
- MAF filter, typically `MAF >= 0.001` or analysis-specific
- Biallelic SNPs only
- Monomorphic variants removed

## Inputs

Expected input format:

```text
${base}/${batch}/raw_${batch}.bed
${base}/${batch}/raw_${batch}.bim
${base}/${batch}/raw_${batch}.fam
```

## Workflow

```text
Raw Genotypes
    |
    v
Pre-imputation QC
    |
    v
Michigan Imputation
    |
    v
Post-imputation QC
    |
    v
Cross-batch Harmonisation & Variant Intersection
    |
    v
Common Variant Scaffold
    |
    v
Second-phase Imputation
    |
    v
Final Harmonised Dataset
```

## Outputs

The pipeline generates:

```text
harmonised/
second_phase_input/
summary_statistics/
logs/
```

Key outputs include:

- Harmonised VCF files
- Lists of shared variants across batches
- Common variant scaffold files
- Second-phase imputation inputs
- Variant count and QC summaries

## Software Requirements

- PLINK2
- bcftools
- tabix
- Michigan Imputation Server
- Will Rayner HRC/1000G checking tools
- 1000 Genomes Phase 3 reference panel (GRCh37)
  
## Current Status

This workflow is under active development and documentation is being expanded as scripts are generalised for reuse across projects.
