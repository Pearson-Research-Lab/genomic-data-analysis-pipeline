# GWAS Analysis Pipeline

This directory contains scripts for genome-wide association analyses (GWAS) using genotype data derived from the post-imputation QC pipeline.

The pipeline implements mixed-model GWAS using a genetic relationship matrix (GRM), and supports both:

- **Mixed-ancestry analyses** *(primary analysis)*
- **European-only analyses** *(sensitivity analysis)*

Analyses are performed separately by ancestry, trial arm, and sex stratum, consistent with the study design described in the analysis plan. 

---

## Overview

Within each ancestry set, GWAS analyses are run separately for:

- **Trial arm**
  - Intervention (`albi`; albiglutide)
  - Placebo

- **Sex stratum**
  - All participants
  - Female-only
  - Male-only

The scripts in this directory support:

1. PCA and GRM generation from GWAS-ready genotype data
2. All-participant GWAS within each trial arm
3. Sex-stratified GWAS within each trial arm

The implementation reflects the study pipeline used in practice, including the use of logistic mixed models for outcomes 26–28 so that relatedness could be handled within the same GRM-based framework. 

---

## Pipeline Components

| Script | Description |
|--------|-------------|
| `config_example.R` | Example configuration file defining ancestry-specific paths and phenotype mappings |
| `pcair.R` | Generates principal components (PCs) and a GRM from GWAS-ready PLINK input |
| `gwas_all.R` | Runs GWAS for all participants within each trial arm |
| `gwas_sexstratified.R` | Runs GWAS for sex-stratified cohorts |

---

## Input Data

### 1. Genotype data

These scripts expect a **GWAS-ready PLINK binary dataset** (`.bed/.bim/.fam`) for each ancestry.

These PLINK datasets are derived from the post-imputation QC outputs and should already have been:

- filtered on imputation quality (`R2 >= 0.5`)
- filtered on minor allele frequency (`MAF >= 0.01`)
- restricted to biallelic variants
- cleared of monomorphic variants
- harmonised to `chr:pos:ref:alt` variant IDs
- merged across chromosomes
- converted to PLINK binary format

Separate genotype datasets are used for:

- mixed-ancestry analysis
- European-only sensitivity analysis

This matches the QC workflow described upstream in the repository. :contentReference[oaicite:5]{index=5}

---

### 2. Phenotype data

Phenotype inputs are supplied as RDS files containing cohort-specific data.

#### Full cohort directory

This directory should contain the all-participant trial-arm cohorts, for example:

- `albi_all.rds`
- `placebo_all.rds`

These files are used by `gwas_all.R`.

### Cohort naming requirement

Cohort RDS filenames must follow this convention:

- Intervention arm: `albi_*`
- Comparator arm: `placebo_*`

- Sex strata:
  - `_all`
  - `_female`
  - `_male`

Examples:
- `albi_all.rds`
- `placebo_female.rds`

These names are used to construct SAP-compliant output filenames.

#### Sex-stratified directory

This directory should contain arm- and sex-specific cohorts, for example:

- `albi_female.rds`
- `albi_male.rds`
- `placebo_female.rds`
- `placebo_male.rds`

These files are used by `gwas_sexstratified.R`.

The scripts use the RDS filename stem as the cohort label in output filenames.

---

## Configuration

All paths, ancestry-specific inputs, and phenotype column mappings are controlled through a config file.

Copy and edit:

```bash
cp config_example.R config_study.R
