#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
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

sap_outcome_code <- function(outcome_id) {
  map <- c(
    "1"  = "E_A1CP",
    "2"  = "E_A1CF",
    "3"  = "E_FGP",
    "4"  = "E_FGF",
    "5"  = "E_BMIP",
    "6"  = "E_BMIF",
    "7"  = "E_WTP",
    "8"  = "E_WTF",
    "9"  = "E_SBPP",
    "10" = "E_SBPF",
    "11" = "E_DBPP",
    "12" = "E_DBPF",
    "13" = "E_HRP",
    "14" = "E_HRF",
    "15" = "E_EGFRP",
    "16" = "E_EGFRF",
    "17" = "E_UACRP",
    "18" = "E_UACRF",
    "19" = "E_HCTP",
    "20" = "E_HCTF",
    "21" = "AE_HG",
    "22" = "AE_GI",
    "23" = "AE_NAU",
    "24" = "AE_THR",
    "25" = "AE_CYS",
    "26" = "O_PRIMARY",
    "27" = "O_HFHOSP",
    "28" = "O_PRIMHF"
  )
  if (!outcome_id %in% names(map)) {
    stop("No SAP phenotype code found for outcome ID: ", outcome_id)
  }
  unname(map[[outcome_id]])
}

make_join_id <- function(x) {
  sub("(:[^:]+){2}$", "", x)
}

extract_outcome_id <- function(filename) {
  m <- str_match(filename, "gwas([0-9]+)_")
  if (is.na(m[1, 2])) stop("Could not extract outcome ID from filename: ", filename)
  m[1, 2]
}

extract_cohort_label <- function(filename) {
  m <- str_match(filename, "gwas[0-9]+_(.*)\\.txt$")
  if (is.na(m[1, 2])) stop("Could not extract cohort label from filename: ", filename)
  m[1, 2]
}

build_cohort_map <- function(analysis_cfg) {
  cohorts <- c(analysis_cfg$cohorts_all, analysis_cfg$cohorts_sex)
  out <- lapply(cohorts, function(x) {
    list(
      label = x$label,
      file = x$file,
      arm = x$arm,
      sex = x$sex
    )
  })
  names(out) <- vapply(cohorts, function(x) x$label, character(1))
  out
}

get_event_count <- function(cohort_file, outcome_id, cfg) {
  pheno <- readRDS(cohort_file)

  pheno_col <- cfg$outcomes$binary[[outcome_id]]

  if (is.null(pheno_col) || !(pheno_col %in% colnames(pheno))) {
    stop("Cannot find binary phenotype column for outcome ", outcome_id, " in ", cohort_file)
  }

  sum(pheno[[pheno_col]] == 1, na.rm = TRUE)
}

format_continuous <- function(x) {
  x$CALLRATE <- 1 - x$F_MISS
  x$STRAND <- "+"
  x$AC <- 2 * x$AF * x$N
  x$BETA <- (x$SCORE / sqrt(x$VAR)) / sqrt(x$AC)
  x$SE <- abs(x$BETA / qnorm(x$PVAL / 2))
  x$QUAL_TYPE <- ifelse(x$TYPE == "IMPUTED", 1, 0)

  x %>%
    mutate(
      MARKER = SNP.x,
      CHR = CHR.x,
      POS = POS.x,
      EFFECT_ALLELE = A2.x,
      OTHER_ALLELE = A1.x,
      N_SAMPLE = N,
      EAF = AF,
      P = PVAL,
      IMPUTED = TYPE,
      P_HWE = P,
      QUAL_SCORE = R2
    ) %>%
    select(
      MARKER, STRAND, CHR, POS, EFFECT_ALLELE, OTHER_ALLELE,
      N_SAMPLE, EAF, BETA, SE, P, IMPUTED, P_HWE,
      CALLRATE, QUAL_SCORE, QUAL_TYPE
    ) %>%
    filter(!duplicated(MARKER))
}

format_binary <- function(x, n_events) {
  x$CALLRATE <- 1 - x$F_MISS
  x$STRAND <- "+"
  x$BETA <- x$SCORE / x$VAR
  x$SE <- abs(x$BETA / qnorm(x$PVAL / 2))
  x$QUAL_TYPE <- ifelse(x$TYPE == "IMPUTED", 1, 0)
  x$N_EVENTS <- n_events

  x %>%
    mutate(
      MARKER = SNP.x,
      CHR = CHR.x,
      POS = POS.x,
      EFFECT_ALLELE = A2.x,
      OTHER_ALLELE = A1.x,
      N_SAMPLE = N,
      EAF = AF,
      P = PVAL,
      IMPUTED = TYPE,
      P_HWE = P,
      QUAL_SCORE = R2
    ) %>%
    select(
      MARKER, STRAND, CHR, POS, EFFECT_ALLELE, OTHER_ALLELE,
      N_SAMPLE, N_EVENTS, EAF, BETA, SE, P, IMPUTED, P_HWE,
      CALLRATE, QUAL_SCORE, QUAL_TYPE
    ) %>%
    filter(!duplicated(MARKER))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("config", "ancestry", "results-dir", "info-file", "hwe-file", "callrate-file", "outdir")
missing_required <- required[!required %in% names(args)]
if (length(missing_required) > 0) {
  stop("Missing required arguments: ", paste(missing_required, collapse = ", "))
}

cfg <- load_config(args[["config"]])
ancestry_name <- args[["ancestry"]]
if (!ancestry_name %in% names(cfg$ancestry)) {
  stop("Unknown ancestry: ", ancestry_name)
}

analysis_cfg <- cfg$ancestry[[ancestry_name]]
cohort_map <- build_cohort_map(analysis_cfg)

results_dir <- args[["results-dir"]]
info_file <- args[["info-file"]]
hwe_file <- args[["hwe-file"]]
callrate_file <- args[["callrate-file"]]
outdir <- args[["outdir"]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

info <- fread(info_file)
hwe <- fread(hwe_file)
cr <- fread(callrate_file)

hwe$ID <- make_join_id(hwe$SNP)
cr$ID <- make_join_id(cr$SNP)

result_files <- list.files(results_dir, pattern = "^gwas[0-9]+_.*\\.txt$", full.names = TRUE)
if (length(result_files) == 0) {
  stop("No GWAS result files found in: ", results_dir)
}

for (f in result_files) {
  a <- fread(f)
  a$ID <- make_join_id(a$SNP)

  outcome_id <- extract_outcome_id(path_file(f))
  cohort_label <- extract_cohort_label(path_file(f))

  if (!cohort_label %in% names(cohort_map)) {
    stop("Cohort label not found in config: ", cohort_label)
  }

  cohort_info <- cohort_map[[cohort_label]]
  sap_code <- sap_outcome_code(outcome_id)
  suffix <- paste(cohort_info$arm, cohort_info$sex, sep = "_")
  out_name <- paste0(sap_code, "_", suffix, ".csv")

  merged <- a %>%
    left_join(info, by = "ID") %>%
    filter(!is.na(TYPE)) %>%
    left_join(hwe, by = "ID") %>%
    left_join(cr, by = "ID")

  if (as.integer(outcome_id) <= 20) {
    out <- format_continuous(merged)
  } else {
    n_events <- get_event_count(cohort_info$file, outcome_id, cfg)
    out <- format_binary(merged, n_events = n_events)
  }

  fwrite(out, file.path(outdir, out_name))
}
