# Genomic Data Analysis Pipeline

End-to-end computational genomics workflow for large-scale human genomic data (QC, imputation, association analysis, visualisation)

## Overview

This repository contains a modular, reproducible workflow for processing and analysing large-scale human genomic data in high-performance computing (HPC) environments.

The workflow spans the full analysis lifecycle from genotype-level data through to association testing and result visualisation, and is designed for scalable, flexible, and reproducible genomic analyses in large cohort settings. It reflects real-world pipelines used in GWAS and cohort-based genomic studies, with compatibility for extension to multi-omics data integration.

This repository contains a **generalised example workflow** and does not include any proprietary or sensitive data.

## Status

This repository contains the structure and documentation for an active computational genomics workflow. Scripts and example implementations are currently being finalised and will be uploaded shortly.

---

## Key Features

* End-to-end genomic data analysis workflow
* Designed for HPC systems (SGE / SLURM)
* Modular, script-based pipeline structure
* Reproducible and parameter-driven execution
* Suitable for large, heterogeneous cohort datasets
  
---

## Workflow Structure

### 1. Genotype Quality Control

* Sample and variant filtering (missingness, MAF, HWE)
* Removal of low-quality samples and variants
* Preparation of high-quality datasets for downstream analysis

### 2. Pre-Imputation Processing

* Variant harmonisation and allele alignment
* Strand checks and reference alignment
* Conversion to imputation-ready formats (VCF)

### 3. Imputation and Data Integration

* Integration with external imputation workflows (e.g. Michigan Imputation Server)
* Handling of phased genotype data
* Harmonisation across cohorts or genotyping platforms, including double-imputation workflows
  
### 4. Post-Imputation Quality Control

* Filtering based on imputation quality metrics (e.g. INFO score)
* Generation of analysis-ready datasets
* Preparation for downstream statistical analysis

### 5. Association Analysis

* Flexible implementation of genome-wide association analyses across multiple outcomes  
* Parameter-driven model specification enabling reproducible and scalable analyses  
* Support for running multiple models within a single workflow  
* Designed for cohort-based analyses with structured phenotype inputs  
* Supports regression and mixed-model approaches depending on study design
  
### 6. Visualisation and Summary Outputs

* Generation of Manhattan and QQ plots
* Summary statistics and QC visualisations

---

## Repository Structure

```
genomic-data-analysis-pipeline/
│
├── README.md
├── scripts/
│   ├── 01_qc.sh
│   ├── 02_pre_imputation.sh
│   ├── 03_post_imputation_qc.sh
│   ├── 04_gwas.sh
│   └── 05_plots.R
│
├── config/
│   └── parameters.txt
│
├── example_data/
│
└── results/
```

---

## Tools & Technologies

* PLINK / PLINK2
* bcftools
* R (data.table, qqman; mixed-model association methods using GMMAT and SAIGE)
* Bash / Linux scripting
* High-performance computing (SGE / SLURM)

---

## External Tools and Resources

This workflow incorporates established community tools for genomic data processing, including strand checking and reference alignment workflows developed by the Wellcome Centre for Human Genetics (Oxford), including scripts by Will Rayner.

All external tools are used in accordance with their respective documentation and guidelines.

---

## Reproducibility

The workflow is designed to support reproducible research practices:

* Modular scripts for each stage of analysis
* Clear separation of input, processing, and output
* Parameter-driven configuration via `config/parameters.txt`
* Version-controlled codebase

Large genomic data files and intermediate outputs are excluded from version control via `.gitignore`.

---

## Example Usage

```bash
bash scripts/01_qc.sh
bash scripts/02_pre_imputation.sh
bash scripts/03_post_imputation_qc.sh
bash scripts/04_gwas.sh
Rscript scripts/05_plots.R
```

---

## Notes

* This repository contains example scripts and workflow structure only
* No real genomic or patient data are included
* File names and paths are generalised for demonstration purposes

---

## Author

Eram Haider-McInnes
Computational Genomics | Genetic Epidemiology

---

## License

This project is licensed under the MIT License.

