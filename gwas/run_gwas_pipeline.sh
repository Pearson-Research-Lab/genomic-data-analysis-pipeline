#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# run_gwas_pipeline.sh
#
# Run GWAS workflow for one ancestry set:
#   1) PCA/GRM generation
#   2) all-participant GWAS
#   3) sex-stratified GWAS
#   4) SAP-formatting / renaming of result files
#
# Usage:
#   bash run_gwas_pipeline.sh mixed /path/to/config_study.R
#   bash run_gwas_pipeline.sh eur   /path/to/config_study.R
###############################################################################

if [[ $# -ne 2 ]]; then
  echo "Usage: bash run_gwas_pipeline.sh <mixed|eur> <config_file>"
  exit 1
fi

ANCESTRY="$1"
CONFIG_FILE="$2"

if [[ "${ANCESTRY}" != "mixed" && "${ANCESTRY}" != "eur" ]]; then
  echo "Error: ancestry must be 'mixed' or 'eur'"
  exit 1
fi

if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "Error: config file not found: ${CONFIG_FILE}"
  exit 1
fi

# -----------------------------------------------------------------------------
# User-editable paths
# -----------------------------------------------------------------------------
R_BIN="Rscript"

# These are expected to exist outside the config for now.
# Edit to match your environment.
INFO_FILE="/path/to/extracted_info.txt"
HWE_FILE="/path/to/ho_allchrplusxhardy.hwe"
CALLRATE_FILE="/path/to/ho_allchrplusxhardy.lmiss"

# Output roots
GWAS_BASE_OUT="/path/to/gwas_outputs"

PCAIR_OUT="${GWAS_BASE_OUT}/${ANCESTRY}/pcair"
GWAS_ALL_OUT="${GWAS_BASE_OUT}/${ANCESTRY}/results_all"
GWAS_SEX_OUT="${GWAS_BASE_OUT}/${ANCESTRY}/results_sex"

FORMATTED_ALL_OUT="${GWAS_BASE_OUT}/${ANCESTRY}/formatted_all"
FORMATTED_SEX_OUT="${GWAS_BASE_OUT}/${ANCESTRY}/formatted_sex"

# -----------------------------------------------------------------------------
# Pull genotype prefix from config using R
# -----------------------------------------------------------------------------
GENOTYPE_PREFIX=$(
  "${R_BIN}" -e "env <- new.env(); sys.source('${CONFIG_FILE}', envir = env); cat(env\$cfg\$ancestry[['${ANCESTRY}']]\$genotype_prefix)" \
  2>/dev/null
)

if [[ -z "${GENOTYPE_PREFIX}" ]]; then
  echo "Error: could not read genotype_prefix for ancestry '${ANCESTRY}' from config"
  exit 1
fi

echo "========================================"
echo "Running GWAS pipeline for ancestry: ${ANCESTRY}"
echo "Config: ${CONFIG_FILE}"
echo "Genotype prefix: ${GENOTYPE_PREFIX}"
echo "========================================"

mkdir -p "${PCAIR_OUT}" "${GWAS_ALL_OUT}" "${GWAS_SEX_OUT}" "${FORMATTED_ALL_OUT}" "${FORMATTED_SEX_OUT}"

# -----------------------------------------------------------------------------
# 1. PCA / GRM
# -----------------------------------------------------------------------------
echo "[1/5] Running PCA / GRM generation"
"${R_BIN}" pcair.R \
  --genotype-prefix "${GENOTYPE_PREFIX}" \
  --outdir "${PCAIR_OUT}"

# -----------------------------------------------------------------------------
# 2. All-participant GWAS
# -----------------------------------------------------------------------------
echo "[2/5] Running all-participant GWAS"
"${R_BIN}" gwas_all.R \
  --config "${CONFIG_FILE}" \
  --ancestry "${ANCESTRY}" \
  --outdir "${GWAS_ALL_OUT}"

# -----------------------------------------------------------------------------
# 3. Sex-stratified GWAS
# -----------------------------------------------------------------------------
echo "[3/5] Running sex-stratified GWAS"
"${R_BIN}" gwas_sexstratified.R \
  --config "${CONFIG_FILE}" \
  --ancestry "${ANCESTRY}" \
  --outdir "${GWAS_SEX_OUT}"

# -----------------------------------------------------------------------------
# 4. Format / rename all-participant GWAS results
# -----------------------------------------------------------------------------
echo "[4/5] Formatting and renaming all-participant GWAS results"
"${R_BIN}" format_and_rename_gwas.R \
  --config "${CONFIG_FILE}" \
  --ancestry "${ANCESTRY}" \
  --results-dir "${GWAS_ALL_OUT}/results" \
  --info-file "${INFO_FILE}" \
  --hwe-file "${HWE_FILE}" \
  --callrate-file "${CALLRATE_FILE}" \
  --outdir "${FORMATTED_ALL_OUT}"

# -----------------------------------------------------------------------------
# 5. Format / rename sex-stratified GWAS results
# -----------------------------------------------------------------------------
echo "[5/5] Formatting and renaming sex-stratified GWAS results"
"${R_BIN}" format_and_rename_gwas.R \
  --config "${CONFIG_FILE}" \
  --ancestry "${ANCESTRY}" \
  --results-dir "${GWAS_SEX_OUT}/results" \
  --info-file "${INFO_FILE}" \
  --hwe-file "${HWE_FILE}" \
  --callrate-file "${CALLRATE_FILE}" \
  --outdir "${FORMATTED_SEX_OUT}"

echo "========================================"
echo "GWAS pipeline completed for ancestry: ${ANCESTRY}"
echo "Outputs:"
echo "  PCA/GRM:            ${PCAIR_OUT}"
echo "  GWAS all:           ${GWAS_ALL_OUT}"
echo "  GWAS sex:           ${GWAS_SEX_OUT}"
echo "  Formatted all:      ${FORMATTED_ALL_OUT}"
echo "  Formatted sex:      ${FORMATTED_SEX_OUT}"
echo "========================================"
