#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(qqman)
  library(fs)
})

parse_args <- function(args) {
  res <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (!startsWith(key, "--")) {
      stop("Arguments must be supplied as --key value")
    }
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

is_binary_or_outcome <- function(filename) {
  startsWith(filename, "AE_") || startsWith(filename, "O_")
}

make_plot_title <- function(file_name, dat) {
  n_sample <- NA
  n_events <- NA

  if ("N_SAMPLE" %in% names(dat)) {
    n_sample <- unique(na.omit(dat$N_SAMPLE))
    if (length(n_sample) > 0) n_sample <- n_sample[1] else n_sample <- NA
  }

  if ("N_EVENTS" %in% names(dat)) {
    n_events <- unique(na.omit(dat$N_EVENTS))
    if (length(n_events) > 0) n_events <- n_events[1] else n_events <- NA
  }

  if (is_binary_or_outcome(file_name) && !is.na(n_events)) {
    return(paste0(file_name, " | N=", n_sample, " | Events=", n_events))
  }

  if (!is.na(n_sample)) {
    return(paste0(file_name, " | N=", n_sample))
  }

  file_name
}

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

required_cols <- c("CHR", "POS", "P")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Input file is missing required columns: ", paste(missing_cols, collapse = ", "))
}

plot_dat <- dat[, .(CHR, POS, P)]
plot_dat <- plot_dat[!is.na(CHR) & !is.na(POS) & !is.na(P)]
plot_dat <- plot_dat[P > 0 & P <= 1]

if (nrow(plot_dat) == 0) {
  stop("No valid rows available for plotting after filtering")
}

title_text <- make_plot_title(file_name, dat)

manhattan_png <- file.path(outdir, paste0(file_stem, "_manhattan.png"))
qq_png <- file.path(outdir, paste0(file_stem, "_qq.png"))

png(manhattan_png, width = 1400, height = 900, res = 120)
manhattan(
  plot_dat,
  chr = "CHR",
  bp = "POS",
  p = "P",
  snp = NULL,
  main = title_text
)
dev.off()

png(qq_png, width = 900, height = 900, res = 120)
qq(
  plot_dat$P,
  main = title_text
)
dev.off()

message("Saved Manhattan plot: ", manhattan_png)
message("Saved QQ plot: ", qq_png)
