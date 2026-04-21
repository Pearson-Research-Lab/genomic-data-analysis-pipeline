# GWAS Analysis Pipeline

This directory contains scripts for genome-wide association analyses (GWAS) using genotype data derived from the post-imputation QC pipeline.

The pipeline implements mixed-model GWAS using a genetic relationship matrix (GRM), and supports both:

- **Mixed-ancestry analyses** *(primary analysis)*
- **European-only analyses** *(sensitivity analysis)*

Analyses are performed separately by ancestry, trial arm, and sex stratum, consistent with the study design.

---

## Overview

Within each ancestry set, GWAS analyses are run separately for:

- **Trial arm**
  - Intervention (`albi`)
  - Placebo

- **Sex stratum**
  - All participants
  - Female-only
  - Male-only

The pipeline consists of:

1. PCA and GRM generation  
2. All-participant GWAS  
3. Sex-stratified GWAS  
4. SAP-compliant result formatting  
5. (Optional) downstream visualisation  

---

## Directory Structure

```
gwas/
├── README.md
├── config_example.R
├── pcair.R
├── gwas_all.R
├── gwas_sexstratified.R
├── format_and_rename_gwas.R
└── run_gwas_pipeline.sh
```

---


---

## Input Data

### Genotype data

GWAS-ready PLINK binary files (`.bed/.bim/.fam`) derived from post-imputation QC:

- INFO filtered (R² ≥ 0.5)
- MAF ≥ 0.01
- Biallelic variants only
- Monomorphic variants removed
- Variant IDs formatted as `chr:pos:ref:alt`

Separate datasets are used for:
- mixed ancestry (primary)
- European-only (sensitivity)

---

### Phenotype data

Cohort data are provided as RDS files.

Each cohort must include:

#### All-participant cohorts
- `ID`
- `Age`
- `Sex`
- `Binary.Insulin`
- `PC1–PC10`
- phenotype columns defined in config

#### Sex-stratified cohorts
- `ID`
- `Age`
- `Binary.Insulin`
- `PC1–PC10`
- phenotype columns defined in config

---

## Results visualisation

Formatted GWAS result files can be visualised using `results_visualisation.R`.

Example:

```bash
Rscript results_visualisation.R \
  --input /path/to/formatted_results/E_A1CP_I_ALL.csv \
  --outdir /path/to/plots
```
