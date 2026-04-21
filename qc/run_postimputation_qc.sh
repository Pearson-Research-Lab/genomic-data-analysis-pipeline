#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# run_postimputation_qc.sh
#
# Purpose:
#   Run post-imputation QC separately on:
#     1) European-only imputed data
#     2) Mixed-ancestry imputed data
#
# Assumptions:
#   - Ancestry assignment was performed upstream
#   - European-only and mixed-ancestry datasets were imputed separately
#   - Input consists of per-chromosome bgzipped VCFs (.vcf.gz), chr1-chr22
#   - Michigan Imputation Server output includes INFO/R2
#
# QC applied:
#   - Exclude variants with R2 < 0.5
#   - Exclude variants with MAF < 0.01
#   - Retain biallelic SNPs only
#   - Remove monomorphic variants (captured by MAF filter)
#   - Harmonise variant IDs to chr:pos:ref:alt
#   - No HWE filtering
#
# Output:
#   - Post-QC VCFs per chromosome
#   - GWAS-ready PLINK2 pfiles per chromosome
###############################################################################

############################
# User-configurable paths
############################

BASE_DIR="/path/to/project"

EUR_INPUT_DIR="${BASE_DIR}/imputed/eur"
MIXED_INPUT_DIR="${BASE_DIR}/imputed/mixed"

EUR_OUTPUT_DIR="${BASE_DIR}/qc/postimputation/eur_qced"
MIXED_OUTPUT_DIR="${BASE_DIR}/qc/postimputation/mixed_qced"

TMP_DIR="${BASE_DIR}/qc/postimputation/tmp"

BCFTOOLS="bcftools"
PLINK2="plink2"

############################
# QC thresholds
############################

R2_MIN="0.5"
MAF_MIN="0.01"

############################
# Create output directories
############################

mkdir -p "${EUR_OUTPUT_DIR}/"{vcf,pfile,logs}
mkdir -p "${MIXED_OUTPUT_DIR}/"{vcf,pfile,logs}
mkdir -p "${TMP_DIR}"

############################
# Check required software
############################

command -v "${BCFTOOLS}" >/dev/null 2>&1 || {
    echo "[ERROR] bcftools not found in PATH" >&2
    exit 1
}

command -v "${PLINK2}" >/dev/null 2>&1 || {
    echo "[ERROR] plink2 not found in PATH" >&2
    exit 1
}

############################
# Function: process dataset
############################
process_dataset() {
    local dataset_label="$1"   # eur or mixed
    local input_dir="$2"
    local output_dir="$3"

    echo "[INFO] Processing dataset: ${dataset_label}"

    for chr in {1..22}; do
        local in_vcf="${input_dir}/chr${chr}.vcf.gz"
        local stage1_vcf="${output_dir}/vcf/${dataset_label}_chr${chr}.r2_ids.vcf.gz"
        local final_vcf="${output_dir}/vcf/${dataset_label}_chr${chr}.postqc.vcf.gz"
        local pfile_prefix="${output_dir}/pfile/${dataset_label}_chr${chr}.postqc"
        local log_file="${output_dir}/logs/${dataset_label}_chr${chr}.log"

        if [[ ! -f "${in_vcf}" ]]; then
            echo "[WARNING] Missing input file: ${in_vcf}. Skipping chr${chr}."
            continue
        fi

        if [[ ! -f "${in_vcf}.tbi" ]]; then
            echo "[INFO] Index not found for ${in_vcf}; creating with bcftools index -t"
            "${BCFTOOLS}" index -f -t "${in_vcf}"
        fi

        {
            echo "[INFO] chr${chr}: filtering on R2, retaining biallelic SNPs, setting IDs"

            "${BCFTOOLS}" annotate \
                --set-id +'%CHROM:%POS:%REF:%ALT' \
                "${in_vcf}" -Ou | \
            "${BCFTOOLS}" view \
                -i "INFO/R2>=${R2_MIN}" \
                -m2 -M2 \
                -v snps \
                -Oz \
                -o "${stage1_vcf}"

            "${BCFTOOLS}" index -f -t "${stage1_vcf}"

            echo "[INFO] chr${chr}: importing to PLINK2 and applying MAF filter"

            "${PLINK2}" \
                --vcf "${stage1_vcf}" dosage=DS \
                --double-id \
                --allow-extra-chr \
                --maf "${MAF_MIN}" \
                --make-pgen \
                --out "${pfile_prefix}"

            echo "[INFO] chr${chr}: exporting filtered VCF"

            "${PLINK2}" \
                --pfile "${pfile_prefix}" \
                --export vcf-4.2 bgz \
                --out "${pfile_prefix}"

            mv "${pfile_prefix}.vcf.gz" "${final_vcf}"
            "${BCFTOOLS}" index -f -t "${final_vcf}"

            echo "[INFO] chr${chr}: complete"
        } > "${log_file}" 2>&1
    done
}

############################
# Run EUR and mixed QC
############################

process_dataset "eur" "${EUR_INPUT_DIR}" "${EUR_OUTPUT_DIR}"
process_dataset "mixed" "${MIXED_INPUT_DIR}" "${MIXED_OUTPUT_DIR}"

echo "[INFO] All post-imputation QC completed successfully."
