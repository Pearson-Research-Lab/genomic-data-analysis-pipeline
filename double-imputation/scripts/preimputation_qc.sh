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

base="/cluster/inspired/Go_data/eram/new_doubleimputation_lenientpreimputationqc"
batches=(gdaffy gdctrl gdillumina gsfreeze1 gsfreeze3)

wrperl_script="${base}/wr_stuff/HRC-1000G-check-bim-NoReadKey.pl"
legend="${base}/wr_stuff/1000GP_Phase3_combined.legend"
fastafile="/cluster/ezp_lab/ehaider001/1g_fasta/human_g1k_v37.fasta"

chr="${SGE_TASK_ID}"

for i in "${batches[@]}"; do
    batchdir="${base}/${i}"
    split_dir="${batchdir}/split_by_chr"
    qc_dir="${batchdir}/qc_chr${chr}"

    mkdir -p "$qc_dir"

    input_prefix="${split_dir}/chr${chr}_${i}"
    qc_prefix="${qc_dir}/chr${chr}_${i}_preimpqc"

    echo "Running pre-imputation QC for batch ${i}, chromosome ${chr}"

    plink2 \
        --bfile "$input_prefix" \
        --max-alleles 2 \
        --snps-only just-acgt \
        --geno 0.03 \
        --mind 0.05 \
        --hwe 1e-8 \
        --make-bed \
        --threads 4 \
        --out "$qc_prefix"

    echo "Generating allele frequency file for Will Rayner check"

    plink \
        --bfile "$qc_prefix" \
        --freq \
        --threads 4 \
        --out "${qc_dir}/chr${chr}_${i}"

    echo "Running Will Rayner check"

    cd "$qc_dir"

    perl "$wrperl_script" \
        -b "chr${chr}_${i}_preimpqc.bim" \
        -f "chr${chr}_${i}.frq" \
        -r "$legend" \
        -g "$fastafile" \
        -p EUR

    bash Run-plink.sh

    echo "Exporting Will Rayner-updated files to bgzipped VCF for MIS"

    plink \
        --bfile "chr${chr}_${i}_preimpqc-updated-chr${chr}" \
        --recode vcf bgz \
        --threads 4 \
        --out "chr${chr}_${i}_readyforMIS"

    bcftools index -t "chr${chr}_${i}_readyforMIS.vcf.gz"

    cd - > /dev/null
done

echo "Finished pre-imputation QC for chromosome ${chr}."
