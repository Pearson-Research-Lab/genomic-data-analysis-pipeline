#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(qqman)
  library(fs)
})

# -------------------------------
# Argument parser
# -------------------------------
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

# -------------------------------
# Detect binary/outcome files
# -------------------------------
is_binary_file <- function(filename) {
  startsWith(filename, "AE_") || startsWith(filename, "O_")
}

# -------------------------------
# Compute lambda (inflation)
# -------------------------------
compute_lambda <- function(pvals) {
  chisq <- qchisq(1 - pvals, 1)
  median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
}

# -------------------------------
# Build plot title
# -------------------------------
make_plot_title <- function(file_name, dat) {
  n_sample <- NA
  n_events <- NA

  if ("N_SAMPLE" %in% names(dat)) {
    vals <- unique(na.omit(dat$N_SAMPLE))
    if (length(vals) > 0) n_sample <- vals[1]
  }

  if ("N_EVENTS" %in% names(dat)) {
    vals <- unique(na.omit(dat$N_EVENTS))
    if (length(vals) > 0) n_events <- vals[1]
  }

  title <- file_name

  if (!is.na(n_sample)) {
    title <- paste0(title, " | N=", n_sample)
  }

  if (is_binary_file(file_name) && !is.na(n_events)) {
    title <- paste0(title, " | Events=", n_events)
  }

  title
}

# -------------------------------
# Main
# -------------------------------
args <- parse_args(commandArgs(trailingOnly = TRUE))

required <- c("input", "outdir")
missing_required <- required[!required %in% names(args)]
if (length(missing_required) > 0) {
  stop("Missing required arguments: ", paste(missing_required, collapse = ", "))
}

input_file <- args[["input"]]
outdir <- args[["outdir"]]

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- fread(input_file)
file_name <- path_file(input_file)
file_stem <- path_ext_remove(file_name)

# -------------------------------
# Validate columns
# -------------------------------
required_cols <- c("CHR", "POS", "P")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# -------------------------------
# Clean data
# -------------------------------
plot_dat <- dat[, .(CHR, POS, P)]
plot_dat <- plot_dat[!is.na(CHR) & !is.na(POS) & !is.na(P)]
plot_dat <- plot_dat[P > 0 & P <= 1]

if (nrow(plot_dat) == 0) {
  stop("No valid rows for plotting")
}

# -------------------------------
# Title + lambda
# -------------------------------
title_text <- make_plot_title(file_name, dat)
lambda <- compute_lambda(plot_dat$P)

# -------------------------------
# Manhattan plot
# -------------------------------
manhattan_file <- file.path(outdir, paste0(file_stem, "_manhattan.png"))

png(manhattan_file, width = 1400, height = 900, res = 120)
manhattan(
  plot_dat,
  chr = "CHR",
  bp = "POS",
  p = "P",
  main = title_text,
  genomewideline = -log10(5e-8),
  suggestiveline = -log10(1e-5)
)
dev.off()

# -------------------------------
# QQ plot
# -------------------------------
qq_file <- file.path(outdir, paste0(file_stem, "_qq.png"))

png(qq_file, width = 900, height = 900, res = 120)
qq(plot_dat$P, main = paste0(title_text, " | lambda=", round(lambda, 3)))
dev.off()

# -------------------------------
# Done
# -------------------------------
message("Saved:")
message("  Manhattan: ", manhattan_file)
message("  QQ:        ", qq_file)
