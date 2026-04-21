# GWAS Analysis Pipeline

This directory contains scripts for genome-wide association analyses (GWAS) using post-imputation QC genotype data.

The pipeline implements mixed-model GWAS using a genetic relationship matrix (GRM), and supports both **primary (mixed-ancestry)** and **sensitivity (European-only)** analyses.

---

## Overview

GWAS analyses are performed separately for:

- **Mixed-ancestry dataset** *(primary analysis)*
- **European-only dataset** *(sensitivity analysis)*

Within each ancestry set, analyses are further stratified by:

- **Trial arm**
  - Intervention (albiglutide; `albi`)
  - Placebo

- **Sex**
  - All participants
  - Female-only
  - Male-only

---

## Pipeline Components

| Script                 | Description                                          |
|------------------------|------------------------------------------------------|
| `pcair.R`              | Generates principal components (PCs) and GRM         |
| `gwas_all.R`           | Runs GWAS for all participants within each trial arm |
| `gwas_sexstratified.R` | Runs GWAS for sex-stratified cohorts                 |

---

## Input Data

### Genotype Data

GWAS-ready genotype data derived from post-imputation QC outputs.

Post-imputation QC VCF files were:

- Filtered (R2 ≥ 0.5, MAF ≥ 0.01)
- Restricted to biallelic SNPs
- Harmonised to `chr:pos:ref:alt` variant IDs
- Merged across chromosomes
- Converted to PLINK binary format (`.bed/.bim/.fam`)

The resulting PLINK dataset is used as input to the GWAS scripts.

Separate datasets are used for:
- mixed-ancestry analysis
- European-only sensitivity analysis

---

### Phenotype Data

RDS files containing cohort-specific phenotype data.

#### Full cohort directory

Contains:
- `albi_all.rds`
- `placebo_all.rds`

#### Sex-stratified directory

Contains:
- `albi_female.rds`
- `albi_male.rds`
- `placebo_female.rds`
- `placebo_male.rds`

---

## Configuration

All paths and analysis settings are controlled via a configuration file.

Copy and edit:

```bash
cp config_example.R config_study.R
