#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GENESIS)
  library(SNPRelate)
  library(GWASTools)
  library(Matrix)
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

args <- parse_args(commandArgs(trailingOnly = TRUE))

required <- c("genotype-prefix", "outdir")
missing_required <- required[!required %in% names(args)]
if (length(missing_required) > 0) {
  stop("Missing required arguments: ", paste(missing_required, collapse = ", "))
}

genotype_prefix <- args[["genotype-prefix"]]
outdir <- args[["outdir"]]
pc_count <- if ("pc-count" %in% names(args)) as.integer(args[["pc-count"]]) else 10

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bed <- paste0(genotype_prefix, ".bed")
bim <- paste0(genotype_prefix, ".bim")
fam <- paste0(genotype_prefix, ".fam")

for (f in c(bed, bim, fam)) {
  if (!file.exists(f)) {
    stop("Missing input file: ", f)
  }
}

gdsfile <- file.path(outdir, "gwas_input.gds")

snpgdsBED2GDS(
  bed.fn = bed,
  bim.fn = bim,
  fam.fn = fam,
  out.gdsfn = gdsfile
)

geno_gds <- snpgdsOpen(gdsfile, allow.duplicate = TRUE)

kingmat <- snpgdsIBDKING(
  gdsobj = geno_gds,
  remove.monosnp = TRUE,
  maf = 0.05,
  type = "KING-robust"
)
saveRDS(kingmat, file.path(outdir, "kingmat.rds"))

kincoef <- kingToMatrix(kingmat)

prunedsnps <- snpgdsLDpruning(
  gdsobj = geno_gds,
  sample.id = NULL,
  autosome.only = FALSE,
  remove.monosnp = TRUE,
  maf = 0.05,
  missing.rate = NaN,
  method = "corr",
  slide.max.bp = 500000,
  slide.max.n = NA,
  ld.threshold = 0.316228,
  num.thread = 1
)

pruned <- unlist(prunedsnps, use.names = FALSE)
saveRDS(prunedsnps, file.path(outdir, "prunedsnps.rds"))

pcs <- pcair(
  gdsobj = geno_gds,
  kinobj = kincoef,
  divobj = kincoef,
  eigen.cnt = pc_count,
  kin.thresh = 0.0442,
  snp.include = pruned
)

principalcomponents <- as.data.frame(pcs$vectors)
colnames(principalcomponents) <- paste0("PC", seq_len(ncol(principalcomponents)))
saveRDS(principalcomponents, file.path(outdir, "principalcomponents.rds"))

geno_reader <- GdsGenotypeReader(filename = gdsfile)
geno_data <- GenotypeData(geno_reader)
geno_iter <- GenotypeBlockIterator(geno_data, snpInclude = pruned)

mypcrel <- pcrelate(
  geno_iter,
  pcs = pcs$vectors[, seq_len(pc_count), drop = FALSE],
  training.set = pcs$unrels,
  ibd.probs = FALSE,
  maf.thresh = 0.05,
  maf.bound.method = "filter"
)

grm <- pcrelateToMatrix(mypcrel)
saveRDS(grm, file.path(outdir, "grm.rds"))

snpgdsClose(geno_gds)

message("PCA/GRM generation complete: ", outdir)
