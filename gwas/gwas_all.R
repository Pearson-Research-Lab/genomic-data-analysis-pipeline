#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GMMAT)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

config_file <- args[which(args == "--config") + 1]
ancestry <- args[which(args == "--ancestry") + 1]
outdir <- args[which(args == "--outdir") + 1]

source(config_file)

cfg_a <- cfg$ancestry[[ancestry]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "fits"), showWarnings = FALSE)
dir.create(file.path(outdir, "results"), showWarnings = FALSE)

grm <- readRDS(cfg_a$grm_rds)
genotype_prefix <- cfg_a$genotype_prefix

cohort_files <- list.files(cfg_a$full_cohort_dir, pattern = "\\.rds$", full.names = TRUE)

run_gwas <- function(pheno, cohort_name, outcome, formula, family) {

  fit <- glmmkin(
    formula,
    data = pheno,
    kins = grm,
    id = "ID",
    family = family
  )

  saveRDS(fit, file.path(outdir, "fits", paste0("fit_", outcome, "_", cohort_name, ".rds")))

  glmm.score(
    fit,
    infile = genotype_prefix,
    outfile = file.path(outdir, "results", paste0("gwas", outcome, "_", cohort_name, ".txt"))
  )
}

for (file in cohort_files) {

  pheno <- readRDS(file)
  name <- tools::file_path_sans_ext(basename(file))

  for (o in cfg$outcomes$continuous_no_baseline) {
    f <- as.formula(paste0("pheno", o, " ~ Age + Sex + Binary.Insulin + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
    run_gwas(pheno, name, o, f, gaussian())
  }

  for (o in cfg$outcomes$continuous_with_baseline) {
    f <- as.formula(paste0("pheno", o, " ~ BL_pheno_", o, " + Age + Sex + Binary.Insulin + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
    run_gwas(pheno, name, o, f, gaussian())
  }

  for (o in cfg$outcomes$binary) {
    f <- as.formula(paste0("pheno_", o, " ~ Age + Sex + Binary.Insulin + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
    run_gwas(pheno, name, o, f, binomial())
  }
}

cat("GWAS (all participants) complete\n")
