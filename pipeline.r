#!/usr/bin/env Rscript
#full-pipeline.r - orchestrates the full pipeline: DESeq -> WGCNA -> cleanup

options(repos = c(CRAN = "https://cloud.r-project.org"))

ensure_pkg <- function(pkg, bioc = FALSE) {
  if (!base::requireNamespace(pkg, quietly = TRUE)) {
    if (!base::requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager")
    }
    if (bioc) BiocManager::install(pkg, update = FALSE, ask = FALSE)
    else utils::install.packages(pkg)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

msg_step <- function(...) cat("\033[1;36m", ..., "\033[0m\n", sep = "")
msg_ok   <- function(...) cat("\033[1;32m", ..., "\033[0m\n", sep = "")
red      <- function(x) paste0("\033[1;31m", x, "\033[0m")

ensure_pkg("optparse")

option_list <- list(
  make_option(c("-o", "--outdir"), type = "character", default = "results",
              help = "Output directory [default: %default]"),
  make_option(c("-d", "--database"), type = "character", default = "org.Mmu.eg.db",
              help = "Bioconductor annotation package [default: %default]"),
  make_option(c("-k", "--keytype"), type = "character", default = "ENSEMBL",
              help = "Gene ID type [default: %default]"),
  make_option(c("-c", "--condition"), type = "character", default = "condition",
              help = "Metadata column defining the comparison [default: %default]"),
  make_option(c("-r", "--reference"), type = "character", default = "Control",
              help = "Reference (baseline) level [default: %default]"),
  make_option(c("-f", "--mincount"), type = "integer", default = 10,
              help = "Minimum total count to keep a gene [default: %default]"),
  make_option(c("-s", "--softpower"), type = "integer", default = NULL,
              help = "WGCNA soft-thresholding power. If unset, auto-selects [default: auto]"),
  make_option(c("-n", "--ngenes"), type = "integer", default = 5000,
              help = "WGCNA: number of top-variable genes [default: %default]"),
  make_option(c("-t", "--kmethreshold"), type = "double", default = NULL,
              help = "WGCNA: kME cutoff for hub genes. If unset, uses --maxhubs [default: top-N]"),
  make_option(c("-m", "--maxhubs"), type = "integer", default = 30,
              help = "WGCNA: max hub genes to report [default: %default]"),
  make_option(c("-a", "--ont"), type = "character", default = "BP",
              help = "GO ontology: BP, MF, CC, or ALL [default: %default]"),
  make_option(c("-p", "--pvalue"), type = "double", default = 0.05,
              help = "GO adjusted p-value cutoff [default: %default]"),
  make_option(c("-g", "--showgoterms"), type = "integer", default = 20,
              help = "GO: number of top terms to show in plots [default: %default]")
)

parser <- OptionParser(
  usage = "Rscript run_all.r [options] <counts.csv> <metadata.csv>",
  option_list = option_list,
  description = "\nRuns the full pipeline: DESeq2 differential expression, then WGCNA co-expression, then removes intermediate files.",
  epilogue = "\nExample:\n  Rscript run_all.r data/counts.csv data/meta.csv -o results/cohort1 -c condition -r Control\n"
)

opt <- parse_args(parser, positional_arguments = 2)

counts_file <- opt$args[1]
meta_file   <- opt$args[2]
o <- opt$options

# ---- helper: run a script, halt the whole pipeline if it fails ----
run_step <- function(script, arglist, label) {
  msg_step(paste0("\n===== Running ", label, " ====="))
  status <- system2("Rscript", args = c(script, arglist))
  if (status != 0) {
    stop(red(paste0(label, " failed (exit code ", status,
                    "). Pipeline halted; intermediate files kept for debugging.")),
         call. = FALSE)
  }
}

# ---- Step 1: DESeq pipeline ----
run_step("deseq-pipeline.r",
         c(counts_file, meta_file,
           "-o", o$outdir, "-d", o$database, "-k", o$keytype,
           "-c", o$condition, "-r", o$reference, "-f", o$mincount),
         "DESeq2 differential expression")

# ---- Step 2: WGCNA (build args, adding optional flags only when set) ----
wgcna_args <- c("-o", o$outdir, "-c", o$condition,
                "-n", o$ngenes, "-m", o$maxhubs)
if (!is.null(o$softpower))    wgcna_args <- c(wgcna_args, "-s", o$softpower)
if (!is.null(o$kmethreshold)) wgcna_args <- c(wgcna_args, "-t", o$kmethreshold)

run_step("wgcna.r", wgcna_args, "WGCNA co-expression")

run_step("go.r",
         c("-o", o$outdir, "-d", o$database, "-k", o$keytype,
           "-a", o$ont, "-p", o$pvalue, "-g", o$showgoterms),
         "GO enrichment")

# ---- Step 3: cleanup, only reached on full success ----
intermediate <- file.path(o$outdir, "intermediate")
if (dir.exists(intermediate)) {
  unlink(intermediate, recursive = TRUE)
  msg_ok(paste0("\nPipeline complete. Removed intermediate files in ", intermediate))
} else {
  msg_ok("\nPipeline complete.")
}
