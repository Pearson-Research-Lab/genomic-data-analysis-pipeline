#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GENESIS)
  library(SNPRelate)
  library(GWASTools)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)

genotype_prefix <- args[which(args == "--genotype-prefix") + 1]
outdir <- args[which(args == "--outdir") + 1]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bed <- paste0(genotype_prefix, ".bed")
bim <- paste0(genotype_prefix, ".bim")
fam <- paste0(genotype_prefix, ".fam")

gdsfile <- file.path(outdir, "gwas_input.gds")

snpgdsBED2GDS(bed, bim, fam, gdsfile)

geno <- snpgdsOpen(gdsfile)

# KING kinship
king <- snpgdsIBDKING(geno, remove.monosnp = TRUE)
kinship <- kingToMatrix(king)

# LD pruning
pruned <- snpgdsLDpruning(geno)
pruned_snps <- unlist(pruned, use.names = FALSE)

# PCAir
pcs <- pcair(
  geno,
  kinobj = kinship,
  divobj = kinship,
  eigen.cnt = 10,
  snp.include = pruned_snps
)

pcs_df <- as.data.frame(pcs$vectors)
colnames(pcs_df) <- paste0("PC", 1:ncol(pcs_df))

saveRDS(pcs_df, file.path(outdir, "principalcomponents.rds"))

# PC-Relate GRM
geno_reader <- GdsGenotypeReader(gdsfile)
geno_data <- GenotypeData(geno_reader)
geno_iter <- GenotypeBlockIterator(geno_data, snpInclude = pruned_snps)

pcrel <- pcrelate(
  geno_iter,
  pcs = pcs$vectors[, 1:10],
  training.set = pcs$unrels
)

grm <- pcrelateToMatrix(pcrel)
saveRDS(grm, file.path(outdir, "grm.rds"))

snpgdsClose(geno)

cat("PCA and GRM complete\n")
