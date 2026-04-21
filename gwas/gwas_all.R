#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(GMMAT)
  library(R.utils)
  library(data.table)
  library(SNPRelate)
  library(GWASTools)
  library(fs)
})

parse_args <- function(args) {
  res <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (!startsWith(key, "--")) stop("Arguments must be supplied as --key value")
    key <- sub("^--", "", key)

    if (i == length(args) || startsWith(args[i + 1], "--")) {
      res[[key]] <- TRUE
      i <- i + 1
    } else {
      res[[key]] <- args[i + 1]
      i <- i + 2
    }
  }
  res
}

load_config <- function(config_file) {
  env <- new.env(parent = emptyenv())
  sys.source(config_file, envir = env)
  if (!exists("cfg", envir = env, inherits = FALSE)) {
    stop("Config file must define an object named 'cfg'")
  }
  get("cfg", envir = env, inherits = FALSE)
}

build_pc_terms <- function(pc_count = 10) {
  paste0("PC", seq_len(pc_count), collapse = " + ")
}

build_formula_continuous_no_baseline <- function(pheno_col, include_sex = TRUE, pc_count = 10) {
  base_terms <- c("Age", "Binary.Insulin")
  if (include_sex) base_terms <- c(base_terms, "Sex")
  rhs <- paste(c(base_terms, build_pc_terms(pc_count)), collapse = " + ")
  as.formula(paste(pheno_col, "~", rhs))
}

build_formula_continuous_with_baseline <- function(pheno_col, baseline_col, include_sex = TRUE, pc_count = 10) {
  base_terms <- c("Age", "Binary.Insulin")
  if (include_sex) base_terms <- c(base_terms, "Sex")
  rhs <- paste(c(baseline_col, base_terms, build_pc_terms(pc_count)), collapse = " + ")
  as.formula(paste(pheno_col, "~", rhs))
}

build_formula_binary <- function(pheno_col, include_sex = TRUE, pc_count = 10) {
  base_terms <- c("Age", "Binary.Insulin")
  if (include_sex) base_terms <- c(base_terms, "Sex")
  rhs <- paste(c(base_terms, build_pc_terms(pc_count)), collapse = " + ")
  as.formula(paste(pheno_col, "~", rhs))
}

check_columns_exist <- function(data, cols, cohort_label) {
  missing_cols <- setdiff(cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in cohort '", cohort_label, "': ", paste(missing_cols, collapse = ", "))
  }
}

run_score_test <- function(fit, genotype_prefix, outfile, maf_range, missing_method, nperbatch) {
  glmm.score(
    obj = fit,
    infile = genotype_prefix,
    outfile = outfile,
    MAF.range = maf_range,
    missing.method = missing_method,
    nperbatch = nperbatch,
    verbose = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("config", "ancestry", "outdir")
missing_required <- required[!required %in% names(args)]
if (length(missing_required) > 0) {
  stop("Missing required arguments: ", paste(missing_required, collapse = ", "))
}

cfg <- load_config(args[["config"]])
ancestry_name <- args[["ancestry"]]
if (!ancestry_name %in% names(cfg$ancestry)) stop("Unknown ancestry: ", ancestry_name)

analysis_cfg <- cfg$ancestry[[ancestry_name]]
outdir <- args[["outdir"]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "fits"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "results"), recursive = TRUE, showWarnings = FALSE)

grm <- as.matrix(readRDS(analysis_cfg$grm_rds))
genotype_prefix <- analysis_cfg$genotype_prefix

pc_count <- cfg$gwas$pc_count
maf_range <- cfg$gwas$maf_range
missing_method <- cfg$gwas$missing_method
nperbatch <- cfg$gwas$nperbatch

for (cohort in analysis_cfg$cohorts_all) {
  if (!file.exists(cohort$file)) stop("Missing cohort file: ", cohort$file)

  pheno <- readRDS(cohort$file)
  cohort_label <- cohort$label

  for (outcome_id in names(cfg$outcomes$continuous_no_baseline)) {
    pheno_col <- cfg$outcomes$continuous_no_baseline[[outcome_id]]
    required_cols <- c(pheno_col, "Age", "Sex", "Binary.Insulin", paste0("PC", seq_len(pc_count)), "ID")
    check_columns_exist(pheno, required_cols, cohort_label)

    fit <- glmmkin(
      build_formula_continuous_no_baseline(pheno_col, include_sex = TRUE, pc_count = pc_count),
      data = pheno, kins = grm, id = "ID", verbose = FALSE,
      family = gaussian(link = "identity")
    )

    saveRDS(fit, file.path(outdir, "fits", paste0("fit_", outcome_id, "_", cohort_label, ".rds")))
    run_score_test(
      fit, genotype_prefix,
      file.path(outdir, "results", paste0("gwas", outcome_id, "_", cohort_label, ".txt")),
      maf_range, missing_method, nperbatch
    )
  }

  for (outcome_id in names(cfg$outcomes$continuous_with_baseline)) {
    outcome_spec <- cfg$outcomes$continuous_with_baseline[[outcome_id]]
    pheno_col <- outcome_spec$pheno
    baseline_col <- outcome_spec$baseline
    required_cols <- c(pheno_col, baseline_col, "Age", "Sex", "Binary.Insulin", paste0("PC", seq_len(pc_count)), "ID")
    check_columns_exist(pheno, required_cols, cohort_label)

    fit <- glmmkin(
      build_formula_continuous_with_baseline(pheno_col, baseline_col, include_sex = TRUE, pc_count = pc_count),
      data = pheno, kins = grm, id = "ID", verbose = FALSE,
      family = gaussian(link = "identity")
    )

    saveRDS(fit, file.path(outdir, "fits", paste0("fit_", outcome_id, "_", cohort_label, ".rds")))
    run_score_test(
      fit, genotype_prefix,
      file.path(outdir, "results", paste0("gwas", outcome_id, "_", cohort_label, ".txt")),
      maf_range, missing_method, nperbatch
    )
  }

  for (outcome_id in names(cfg$outcomes$binary)) {
    pheno_col <- cfg$outcomes$binary[[outcome_id]]
    required_cols <- c(pheno_col, "Age", "Sex", "Binary.Insulin", paste0("PC", seq_len(pc_count)), "ID")
    check_columns_exist(pheno, required_cols, cohort_label)

    fit <- glmmkin(
      build_formula_binary(pheno_col, include_sex = TRUE, pc_count = pc_count),
      data = pheno, kins = grm, id = "ID", verbose = FALSE,
      family = binomial(link = "logit")
    )

    saveRDS(fit, file.path(outdir, "fits", paste0("fit_", outcome_id, "_", cohort_label, ".rds")))
    run_score_test(
      fit, genotype_prefix,
      file.path(outdir, "results", paste0("gwas", outcome_id, "_", cohort_label, ".txt")),
      maf_range, missing_method, nperbatch
    )
  }
}

message("All-participant GWAS complete for ancestry: ", ancestry_name)
