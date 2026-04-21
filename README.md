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
bash qc/run_preimputation_qc.sh
bash qc/run_postimputation_qc.sh
bash gwas/run_gwas_pipeline.sh
```

---

## Notes

* This repository contains example scripts and workflow structure only
* No real genomic or patient data are included
* File names and paths are generalised for demonstration purposes

---

## Citations

If you use this pipeline, please cite the following software and resources:

### GWAS and mixed model analysis

- **GMMAT**  
  Chen H, Wang C, Conomos MP, et al.  
  *Control for population structure and relatedness for binary traits in genetic association studies via logistic mixed models.*  
  American Journal of Human Genetics. 2016;98(4):653–666.

- **GENESIS**  
  Conomos MP, Miller MB, Thornton TA.  
  *Robust inference of population structure for ancestry prediction and correction of stratification in the presence of relatedness.*  
  Genetic Epidemiology. 2015;39(4):276–293.

---

### Genotype processing and QC

- **PLINK / PLINK2**  
  Chang CC, Chow CC, Tellier LC, et al.  
  *Second-generation PLINK: rising to the challenge of larger and richer datasets.*  
  GigaScience. 2015;4:7.

- **bcftools**  
  Danecek P, Bonfield JK, Liddle J, et al.  
  *Twelve years of SAMtools and BCFtools.*  
  GigaScience. 2021;10(2):giab008.

---

### Reference panels and imputation

- **1000 Genomes Project (Phase 3)**  
  The 1000 Genomes Project Consortium.  
  *A global reference for human genetic variation.*  
  Nature. 2015;526:68–74.

- **Michigan Imputation Server**  
  Das S, Forer L, Schönherr S, et al.  
  *Next-generation genotype imputation service and methods.*  
  Nature Genetics. 2016;48:1284–1287.

---

### GWAS visualisation

- **qqman (R package)**  
  Turner SD.  
  *qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots.*  
  bioRxiv. 2014. doi:10.1101/005165

---

## Acknowledgements

This pipeline was developed based on a study-specific statistical analysis plan (SAP) and internal methodological guidance.  
Implementation reflects practical decisions made during analysis, including model selection and QC thresholds.

---

## Author

Eram Haider-McInnes
Computational Genomics | Genetic Epidemiology

---

## License

This project is licensed under the MIT License.

