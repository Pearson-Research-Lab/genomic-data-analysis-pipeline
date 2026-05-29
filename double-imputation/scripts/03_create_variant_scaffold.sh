```bash
#!/bin/bash
#$ -cwd
#$ -V
#$ -pe smp 4
#$ -t 1-22

# get_commonsnps.sh
#
# Post-imputation QC and common variant scaffold generation for double imputation.
#
# This script:
#   1. Filters each imputed batch by INFO/R2
#   2. Calculates/fills MAF and filters by MAF
#   3. Identifies variants shared across all batches using bcftools isec
#   4. Verifies that the common-site set is present in every batch
#   5. Extracts common variants from each batch
#   6. Merges batches across samples to create a second-phase imputation scaffold
#
# Optional/testing:
#   Will Rayner checks on the merged scaffold can be enabled by setting RUN_WR=true.
#   This is disabled by default because merged post-imputation VCFs may not always
#   convert cleanly back to PLINK format.

set -euo pipefail

############################################
# User configuration
############################################

base="/cluster/inspired/Go_data/eram/new_doubleimputation_lenientpreimputationqc"

batches=(gdaffy gdctrl gdillumina gsfreeze1 gsfreeze3)

chr="${SGE_TASK_ID:-}"

# Input VCF pattern:
# Expected file:
#   ${base}/${batch}/imputed/chr${chr}.dose.vcf.gz
input_pattern='${base}/${batch}/imputed/chr${chr}.dose.vcf.gz'

# Double-imputation scaffold thresholds
r2_threshold="0.99"
maf_threshold="0.001"

# Merge behaviour
# If sample IDs overlap between batches, bcftools merge normally stops.
# Set this to true only if duplicate sample IDs are expected and safe.
allow_force_samples=false

# Optional Will Rayner section
run_wr=false
wrperl_script="${base}/wr_stuff/HRC-1000G-check-bim-NoReadKey.pl"
legend="${base}/wr_stuff/1000GP_Phase3_combined.legend"
fastafile="/cluster/ezp_lab/ehaider001/1g_fasta/human_g1k_v37.fasta"
population="EUR"

############################################
# Checks
############################################

if [[ -z "$chr" ]]; then
    echo "ERROR: SGE_TASK_ID is not set. Run as an array job or set chr manually." >&2
    exit 1
fi

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found." >&2; exit 1; }
command -v plink >/dev/null 2>&1 || { echo "ERROR: plink not found." >&2; exit 1; }

isec_dir="${base}/isec_chr${chr}"
summary_dir="${base}/summary_statistics"
logs_dir="${base}/logs"

mkdir -p "$isec_dir" "$summary_dir" "$logs_dir"

step1_counts="${summary_dir}/postimpqc_step1_nvar_r2_${r2_threshold}_chr${chr}.txt"
step2_counts="${summary_dir}/postimpqc_step2_nvar_maf_${maf_threshold}_chr${chr}.txt"
common_counts="${summary_dir}/common_site_check_chr${chr}.txt"

: > "$step1_counts"
: > "$step2_counts"
: > "$common_counts"

echo "Starting common SNP scaffold generation for chromosome ${chr}"
echo "R2 threshold: ${r2_threshold}"
echo "MAF threshold: ${maf_threshold}"

############################################
# Step 1: Filter each batch by R2 and MAF
############################################

filtered_files=()

for batch in "${batches[@]}"; do
    batchdir="${base}/${batch}/imputed"

    in_vcf="${batchdir}/chr${chr}.dose.vcf.gz"

    r2_vcf="${batchdir}/chr${chr}_${batch}_r2_${r2_threshold}.vcf.gz"
    maf_tagged_vcf="${batchdir}/chr${chr}_${batch}_r2_${r2_threshold}_maftagged.vcf.gz"
    maf_vcf="${batchdir}/chr${chr}_${batch}_r2_${r2_threshold}_maf_${maf_threshold}.vcf.gz"

    if [[ ! -f "$in_vcf" ]]; then
        echo "ERROR: Missing input VCF: $in_vcf" >&2
        exit 1
    fi

    if [[ ! -f "${in_vcf}.tbi" && ! -f "${in_vcf}.csi" ]]; then
        echo "Indexing input VCF: $in_vcf"
        bcftools index -t "$in_vcf"
    fi

    echo "Filtering ${batch}, chr${chr} by INFO/R2 >= ${r2_threshold}"

    bcftools view \
        -i "INFO/R2>=${r2_threshold}" \
        "$in_vcf" \
        -Oz -o "$r2_vcf"

    bcftools index -t "$r2_vcf"
    bcftools index -n "$r2_vcf" >> "$step1_counts"

    echo "Filling MAF tag for ${batch}, chr${chr}"

    bcftools +fill-tags "$r2_vcf" \
        -Oz -o "$maf_tagged_vcf" \
        -- -t MAF

    bcftools index -t "$maf_tagged_vcf"

    echo "Filtering ${batch}, chr${chr} by MAF >= ${maf_threshold}"

    bcftools view \
        -i "INFO/MAF>=${maf_threshold}" \
        "$maf_tagged_vcf" \
        -Oz -o "$maf_vcf"

    bcftools index -t "$maf_vcf"
    bcftools index -n "$maf_vcf" >> "$step2_counts"

    filtered_files+=("$maf_vcf")
done

awk '{sum += $1} END {print sum}' "$step1_counts" > "${summary_dir}/totalnvar_step1_r2_${r2_threshold}_chr${chr}.txt"
awk '{sum += $1} END {print sum}' "$step2_counts" > "${summary_dir}/totalnvar_step2_maf_${maf_threshold}_chr${chr}.txt"

############################################
# Step 2: Identify variants shared by all batches
############################################

echo "Identifying variants present in all ${#filtered_files[@]} batches"

bcftools isec \
    -n="${#filtered_files[@]}" \
    -c none \
    -p "$isec_dir" \
    "${filtered_files[@]}"

common_vcf="${isec_dir}/0000.vcf.gz"
common_sites="${isec_dir}/common_sites.tsv"

if [[ ! -f "$common_vcf" ]]; then
    echo "ERROR: Expected common variant file not found: $common_vcf" >&2
    exit 1
fi

if [[ ! -f "${common_vcf}.tbi" && ! -f "${common_vcf}.csi" ]]; then
    bcftools index -t "$common_vcf"
fi

common_n=$(bcftools index -n "$common_vcf")

if [[ "$common_n" -eq 0 ]]; then
    echo "ERROR: No common variants found for chromosome ${chr}." >&2
    exit 1
fi

echo "Common variants identified for chr${chr}: ${common_n}"

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    "$common_vcf" > "$common_sites"

############################################
# Step 3: Verify common sites are present in every batch
############################################

echo "Verifying common sites are present in every filtered batch"

common_filtered_files=()

for idx in "${!batches[@]}"; do
    batch="${batches[$idx]}"
    in_vcf="${filtered_files[$idx]}"
    out_vcf="${base}/${batch}/imputed/chr${chr}_${batch}_common_scaffold.vcf.gz"

    bcftools view \
        -R "$common_sites" \
        "$in_vcf" \
        -Oz -o "$out_vcf"

    bcftools index -t "$out_vcf"

    extracted_n=$(bcftools index -n "$out_vcf")

    echo -e "${batch}\t${extracted_n}\t${common_n}" >> "$common_counts"

    if [[ "$extracted_n" -ne "$common_n" ]]; then
        echo "ERROR: Common-site check failed for ${batch}, chr${chr}." >&2
        echo "Expected ${common_n}, extracted ${extracted_n}." >&2
        echo "See: $common_counts" >&2
        exit 1
    fi

    common_filtered_files+=("$out_vcf")
done

echo "Common-site verification passed for all batches."

############################################
# Step 4: Check sample overlap before merge
############################################

sample_list="${isec_dir}/sample_ids_chr${chr}.txt"
duplicate_samples="${isec_dir}/duplicate_sample_ids_chr${chr}.txt"

: > "$sample_list"

for vcf in "${common_filtered_files[@]}"; do
    bcftools query -l "$vcf" >> "$sample_list"
done

sort "$sample_list" | uniq -d > "$duplicate_samples"

if [[ -s "$duplicate_samples" ]]; then
    echo "WARNING: Duplicate sample IDs detected across batches."
    echo "Duplicate sample list: $duplicate_samples"

    if [[ "$allow_force_samples" == true ]]; then
        echo "allow_force_samples=true, so bcftools merge will use --force-samples."
        merge_extra_args=(--force-samples)
    else
        echo "ERROR: Duplicate sample IDs found and allow_force_samples=false." >&2
        echo "Either fix sample IDs before merging, or set allow_force_samples=true if this is expected." >&2
        exit 1
    fi
else
    merge_extra_args=()
fi

############################################
# Step 5: Merge batches across samples
############################################

merged_vcf="${isec_dir}/r2_${r2_threshold}_maf_${maf_threshold}_chr${chr}_second_phase_scaffold.vcf.gz"

echo "Merging common scaffold VCFs for chr${chr}"

bcftools merge \
    "${merge_extra_args[@]}" \
    "${common_filtered_files[@]}" \
    -Oz -o "$merged_vcf"

bcftools index -t "$merged_vcf"

merged_n=$(bcftools index -n "$merged_vcf")

if [[ "$merged_n" -ne "$common_n" ]]; then
    echo "ERROR: Merged VCF variant count does not match common variant count." >&2
    echo "Common variants: ${common_n}" >&2
    echo "Merged variants: ${merged_n}" >&2
    exit 1
fi

echo "Merged scaffold complete: $merged_vcf"
echo "Merged scaffold variant count: ${merged_n}"

############################################
# Optional/testing: Will Rayner check
############################################

if [[ "$run_wr" == true ]]; then
    echo "Running optional Will Rayner check on merged scaffold."

    if [[ ! -f "$wrperl_script" ]]; then
        echo "ERROR: Will Rayner script not found: $wrperl_script" >&2
        exit 1
    fi

    if [[ ! -f "$legend" ]]; then
        echo "ERROR: Legend file not found: $legend" >&2
        exit 1
    fi

    if [[ ! -f "$fastafile" ]]; then
        echo "ERROR: FASTA file not found: $fastafile" >&2
        exit 1
    fi

    wr_prefix="${isec_dir}/r2_${r2_threshold}_maf_${maf_threshold}_chr${chr}_second_phase_scaffold"

    plink \
        --vcf "$merged_vcf" \
        --make-bed \
        --threads 4 \
        --out "$wr_prefix"

    plink \
        --bfile "$wr_prefix" \
        --freq \
        --threads 4 \
        --out "${isec_dir}/chr${chr}_scaffold"

    cd "$isec_dir"

    perl "$wrperl_script" \
        -b "$(basename "${wr_prefix}.bim")" \
        -f "chr${chr}_scaffold.frq" \
        -r "$legend" \
        -g "$fastafile" \
        -p "$population"

    bash Run-plink.sh

    cd - > /dev/null

    echo "Optional Will Rayner check completed for chr${chr}."
else
    echo "Skipping optional Will Rayner check. Set run_wr=true to enable."
fi

echo "Finished chromosome ${chr}."
```
