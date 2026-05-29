#!/bin/bash
#$ -cwd
#$ -V
#$ -pe smp 4

# process_raw.sh
#
# Rename raw PLINK binary variant IDs to chr:pos:ref:alt format
# and split each batch into autosomal chromosome-specific PLINK files.

base="/cluster/inspired/Go_data/eram/firstphase_preimputationqc"
batches=(gdaffy gdctrl gdillumina gsfreeze1 gsfreeze3)

for i in "${batches[@]}"; do
    batchdir="${base}/${i}"
    raw_prefix="${batchdir}/raw_${i}"
    split_dir="${batchdir}/split_by_chr"

    mkdir -p "$split_dir"

    echo "Renaming variants for batch: ${i}"

    plink2 \
        --bfile "$raw_prefix" \
        --set-all-var-ids @:#:\$r:\$a \
        --new-id-max-allele-len 1000 missing \
        --make-bed \
        --threads 4 \
        --out "${batchdir}/raw_${i}_renamed"

    echo "Splitting ${i} into chromosomes 1-22"

    for chr in $(seq 1 22); do
        plink2 \
            --bfile "${batchdir}/raw_${i}_renamed" \
            --chr "$chr" \
            --make-bed \
            --threads 4 \
            --out "${split_dir}/chr${chr}_${i}"
    done
done

echo "Finished processing raw PLINK files."
