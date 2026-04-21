#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# run_preimputation_qc.sh
#
# Purpose:
#   Run pre-imputation QC separately on:
#     1) European-only genotype data
#     2) Mixed-ancestry genotype data
#
# Assumptions:
#   - Ancestry assignment was performed upstream
#   - Variant IDs were harmonised upstream
#   - Input consists of per-chromosome bgzipped VCFs (.vcf.gz), chr1-chr22
#   - Files are indexed with .tbi indices
#
# QC applied:
#   EUR:
#     - HWE p < 1e-6 excluded
#     - MAF < 0.01 excluded
#     - sample missingness > 0.02 excluded
#     - variant missingness > 0.03 excluded
#
#   Mixed:
#     - MAF < 0.01 excluded
#     - sample missingness > 0.02 excluded
#     - variant missingness > 0.03 excluded
#     - no HWE filter
#
# Output:
#   QCed VCFs for each chromosome in separate EUR and mixed output folders
###############################################################################

############################
# User-configurable paths
############################

BASE_DIR="/path/to/project"

EUR_INPUT_DIR="${BASE_DIR}/eur_vcf"
MIXED_INPUT_DIR="${BASE_DIR}/mixed_vcf"

EUR_OUTPUT_DIR="${BASE_DIR}/qc/eur_qced"
MIXED_OUTPUT_DIR="${BASE_DIR}/qc/mixed_qced"

TMP_DIR="${BASE_DIR}/qc/tmp"

PLINK2="plink2"
BCFTOOLS="bcftools"

THREADS=4

############################
# QC thresholds
############################

HWE_EUR="1e-6"
MAF_MIN="0.01"

# Interpreted as:
# sample call rate >= 98%  => missingness <= 0.02
# variant genotyping rate >= 97% => missingness <= 0.03
MIND_MAX="0.02"
GENO_MAX="0.03"

############################
# Create output directories
############################

mkdir -p "${EUR_OUTPUT_DIR}" "${MIXED_OUTPUT_DIR}" "${TMP_DIR}"
mkdir -p "${EUR_OUTPUT_DIR}/logs" "${MIXED_OUTPUT_DIR}/logs"

############################
# Check required software
############################

command -v "${PLINK2}" >/dev/null 2>&1 || {
    echo "[ERROR] plink2 not found in PATH" >&2
    exit 1
}

command -v "${BCFTOOLS}" >/dev/null 2>&1 || {
    echo "[ERROR] bcftools not found in PATH" >&2
    exit 1
}

############################
# Function: process dataset
############################
process_dataset() {
    local dataset_label="$1"     # eur or mixed
    local input_dir="$2"
    local output_dir="$3"
    local apply_hwe="$4"         # yes or no

    echo "[INFO] Processing dataset: ${dataset_label}"

    for chr in {1..22}; do
        local in_vcf="${input_dir}/chr${chr}.vcf.gz"
        local out_prefix="${output_dir}/chr${chr}"
        local out_vcf="${output_dir}/chr${chr}.qc.vcf.gz"
        local log_file="${output_dir}/logs/chr${chr}.log"

        if [[ ! -f "${in_vcf}" ]]; then
            echo "[WARNING] Missing input file: ${in_vcf}. Skipping chr${chr}."
            continue
        fi

        if [[ ! -f "${in_vcf}.tbi" ]]; then
            echo "[INFO] Index not found for ${in_vcf}; creating with bcftools index -t"
            "${BCFTOOLS}" index -f -t "${in_vcf}"
        fi

        echo "[INFO] chr${chr}: running PLINK2 QC"
        {
            if [[ "${apply_hwe}" == "yes" ]]; then
                "${PLINK2}" \
                    --vcf "${in_vcf}" \
                    --allow-extra-chr \
                    --double-id \
                    --snps-only just-acgt \
                    --max-alleles 2 \
                    --maf "${MAF_MIN}" \
                    --mind "${MIND_MAX}" \
                    --geno "${GENO_MAX}" \
                    --hwe "${HWE_EUR}" midp \
                    --make-pgen \
                    --out "${out_prefix}"
            else
                "${PLINK2}" \
                    --vcf "${in_vcf}" \
                    --allow-extra-chr \
                    --double-id \
                    --snps-only just-acgt \
                    --max-alleles 2 \
                    --maf "${MAF_MIN}" \
                    --mind "${MIND_MAX}" \
                    --geno "${GENO_MAX}" \
                    --make-pgen \
                    --out "${out_prefix}"
            fi
        } > "${log_file}" 2>&1

        echo "[INFO] chr${chr}: exporting filtered VCF"
        "${PLINK2}" \
            --pfile "${out_prefix}" \
            --export vcf-4.2 bgz \
            --out "${out_prefix}.qc"

        mv "${out_prefix}.qc.vcf.gz" "${out_vcf}"

        echo "[INFO] chr${chr}: indexing QCed VCF"
        "${BCFTOOLS}" index -f -t "${out_vcf}"

        echo "[INFO] chr${chr}: complete"
    done
}

############################
# Run EUR and mixed QC
############################

process_dataset "eur" "${EUR_INPUT_DIR}" "${EUR_OUTPUT_DIR}" "yes"
process_dataset "mixed" "${MIXED_INPUT_DIR}" "${MIXED_OUTPUT_DIR}" "no"

echo "[INFO] All pre-imputation QC completed successfully."
