#!/bin/bash
#$ -cwd
#$ -V
#$ -pe smp 4
#$ -t 1-22

# preimputation_qc.sh
#
# Run chromosome-specific pre-imputation QC for each batch,
# perform Will Rayner HRC/1000G checks, and export bgzipped VCFs
# ready for Michigan Imputation Server.
#
# Expected input:
#   ${base}/${batch}/split_by_chr/chr${chr}_${batch}.bed
#   ${base}/${batch}/split_by_chr/chr${chr}_${batch}.bim
#   ${base}/${batch}/split_by_chr/chr${chr}_${batch}.fam
#
# Main output:
#   ${base}/${batch}/qc_chr${chr}/chr${chr}_${batch}_readyforMIS.vcf.gz

############################################
# User configuration
############################################

base="/cluster/inspired/Go_data/eram/new_doubleimputation_lenientpreimputationqc"

batches=(gdaffy gdctrl gdillumina gsfreeze1 gsfreeze3)

wrperl_script="${base}/wr_stuff/HRC-1000G-check-bim-NoReadKey.pl"
legend="${base}/wr_stuff/1000GP_Phase3_combined.legend"
fastafile="/cluster/ezp_lab/ehaider001/1g_fasta/human_g1k_v37.fasta"

population="EUR"

threads=4

geno_threshold="0.03"
mind_threshold="0.05"
hwe_threshold="1e-8"

chr="${SGE_TASK_ID:-}"

############################################
# Checks
############################################

if [[ -z "$chr" ]]; then
    echo "ERROR: SGE_TASK_ID is not set. Submit as an array job or set chr manually." >&2
    exit 1
fi

command -v plink2 >/dev/null 2>&1 || { echo "ERROR: plink2 not found." >&2; exit 1; }
command -v plink >/dev/null 2>&1 || { echo "ERROR: plink not found." >&2; exit 1; }
command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found." >&2; exit 1; }
command -v perl >/dev/null 2>&1 || { echo "ERROR: perl not found." >&2; exit 1; }

[[ -f "$wrperl_script" ]] || { echo "ERROR: Will Rayner script not found: $wrperl_script" >&2; exit 1; }
[[ -f "$legend" ]] || { echo "ERROR: Legend file not found: $legend" >&2; exit 1; }
[[ -f "$fastafile" ]] || { echo "ERROR: FASTA file not found: $fastafile" >&2; exit 1; }

echo "Starting pre-imputation QC for chromosome ${chr}"
echo "Base directory: ${base}"
echo "Population for Will Rayner check: ${population}"

############################################
# Run QC per chromosome by batch
############################################

for batch in "${batches[@]}"; do
    batchdir="${base}/${batch}"
    split_dir="${batchdir}/split_by_chr"
    qc_dir="${batchdir}/qc_chr${chr}"

    mkdir -p "$qc_dir"

    input_prefix="${split_dir}/chr${chr}_${batch}"
    qc_prefix="${qc_dir}/chr${chr}_${batch}_preimpqc"
    freq_prefix="${qc_dir}/chr${chr}_${batch}"
    ready_prefix="${qc_dir}/chr${chr}_${batch}_readyforMIS"

    echo "----------------------------------------"
    echo "Batch: ${batch}"
    echo "Chromosome: ${chr}"

    for ext in bed bim fam; do
        if [[ ! -f "${input_prefix}.${ext}" ]]; then
            echo "ERROR: Missing input file: ${input_prefix}.${ext}" >&2
            exit 1
        fi
    done

    echo "Running PLINK2 pre-imputation QC"

    plink2 \
        --bfile "$input_prefix" \
        --max-alleles 2 \
        --snps-only just-acgt \
        --geno "$geno_threshold" \
        --mind "$mind_threshold" \
        --hwe "$hwe_threshold" \
        --make-bed \
        --threads "$threads" \
        --out "$qc_prefix"

    echo "Generating allele frequency file for Will Rayner check"

    plink \
        --bfile "$qc_prefix" \
        --freq \
        --threads "$threads" \
        --out "$freq_prefix"

    echo "Running Will Rayner HRC/1000G check"

    cd "$qc_dir"

    perl "$wrperl_script" \
        -b "$(basename "${qc_prefix}.bim")" \
        -f "$(basename "${freq_prefix}.frq")" \
        -r "$legend" \
        -g "$fastafile" \
        -p "$population"

    if [[ ! -f "Run-plink.sh" ]]; then
        echo "ERROR: Will Rayner check did not create Run-plink.sh in ${qc_dir}" >&2
        exit 1
    fi

    bash Run-plink.sh

    updated_prefix="chr${chr}_${batch}_preimpqc-updated-chr${chr}"

    for ext in bed bim fam; do
        if [[ ! -f "${updated_prefix}.${ext}" ]]; then
            echo "ERROR: Expected Will Rayner updated file missing: ${updated_prefix}.${ext}" >&2
            exit 1
        fi
    done

    echo "Exporting Will Rayner-updated data to bgzipped VCF for MIS"

    plink \
        --bfile "$updated_prefix" \
        --recode vcf bgz \
        --threads "$threads" \
        --out "$(basename "$ready_prefix")"

    if [[ ! -f "$(basename "${ready_prefix}.vcf.gz")" ]]; then
        echo "ERROR: VCF export failed for ${batch}, chr${chr}" >&2
        exit 1
    fi

    bcftools index -t "$(basename "${ready_prefix}.vcf.gz")"

    cd - > /dev/null

    echo "Finished batch ${batch}, chromosome ${chr}"
done

echo "Finished pre-imputation QC for chromosome ${chr}."
