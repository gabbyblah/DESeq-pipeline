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

msg_step <- function(...) cat("\033[36m", ..., "\033[0m\n", sep = "")
msg_ok   <- function(...) cat("\033[32m", ..., "\033[0m\n", sep = "")
red      <- function(x) paste0("\033[1;31m", x, "\033[0m")

ensure_pkg("optparse")

option_list <- list(
  make_option(c("-o", "--outdir"), type = "character", default = "results",
              help = "Output directory (must contain WGCNA go_input.rds) [default: %default]"),
  make_option(c("-d", "--database"), type = "character", default = "org.Mmu.eg.db",
              help = "Bioconductor annotation package [default: %default]"),
  make_option(c("-k", "--keytype"), type = "character", default = "ENSEMBL",
              help = "Gene ID type in the counts data [default: %default]"),
  make_option(c("-a", "--ont"), type = "character", default = "BP",
              help = "GO ontology: BP, MF, CC, or ALL [default: %default]"),
  make_option(c("-p", "--pvalue"), type = "double", default = 0.05,
              help = "Adjusted p-value cutoff for enrichment [default: %default]"),
  make_option(c("-g", "--showgoterms"), type = "integer", default = 20,
              help = "Number of top terms to show in the dotplot [default: %default]")
)

parser <- OptionParser(
  usage = "Rscript go.r [options]",
  option_list = option_list,
  description = "\nGO enrichment analysis on the WGCNA condition-associated module. Reads go_input.rds from --outdir/intermediate, runs over-representation analysis, and writes an enrichment table and dotplot.",
  epilogue = "\nRun the WGCNA stage first - this reads its module output.\n\nExample:\n  Rscript go.r -o results/cohort1 -d org.Mmu.eg.db -a BP\n"
)

opt <- parse_args(parser)

outdir    <- opt$outdir
org_db    <- opt$database
keytype   <- opt$keytype
ontology  <- opt$ont
pvalue    <- opt$pvalue
show_terms <- opt$showgoterms

intermediate <- file.path(outdir, "intermediate")
go_dir <- file.path(outdir, "GO")
dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)

ensure_pkg("clusterProfiler", bioc = TRUE)
ensure_pkg("ggplot2")
ensure_pkg(org_db, bioc = TRUE)

if (!keytype %in% keytypes(get(org_db))) {
  stop(red(paste0("Keytype '", keytype, "' not available for ", org_db,
                  ". Available: ", paste(keytypes(get(org_db)), collapse = ", "))),
       call. = FALSE)
}

if (!ontology %in% c("BP", "MF", "CC", "ALL")) {
  stop(red(paste0("Invalid --ont value '", ontology, "'. Must be one of: BP, MF, CC, ALL")),
       call. = FALSE)
}

go_input_file <- file.path(intermediate, "go_input.rds")
if (!file.exists(go_input_file)) {
  stop(red(paste0("Required input not found: ", go_input_file, "\n",
                  "Run the WGCNA stage first with the same --outdir.")), call. = FALSE)
}

res_df <- read.csv(file.path(outdir, "DESeq", "DE_results.csv"), row.names = 1)
go_input   <- readRDS(go_input_file)
gene_list  <- go_input$selected      # module genes to test
background <- go_input$background     # all analyzed genes (the universe)
module     <- go_input$module        # which module (for labeling)

msg_step(paste0("Running GO enrichment on module '", module, "' (",
                length(gene_list), " genes, ontology: ", ontology, ")"))

ego <- tryCatch(
  enrichGO(
    gene          = gene_list,
    universe      = background,
    OrgDb         = get(org_db),
    keyType       = keytype,
    ont           = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff  = pvalue,
    readable      = TRUE
  ),
  error = function(e) {
    stop(red(paste0(
      "GO enrichment failed.\n",
      "This usually means the database (--database ", org_db,
      ") or keytype (--keytype ", keytype, ") doesn't match your gene IDs.\n",
      "Your genes look like: ", paste(head(gene_list, 3), collapse = ", "), "\n",
      "Original error: ", conditionMessage(e)
    )), call. = FALSE)
  }
)

if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
  msg_step(paste0("No GO terms enriched at padj < ", pvalue,
                  ". This is common for sparsely-annotated organisms. No outputs written."))
  quit(save = "no")
}

ego_df <- as.data.frame(ego)
write.csv(ego_df, file.path(go_dir, "GO_results.csv"), row.names = FALSE)
msg_ok("GO results saved to: ", file.path(go_dir, "GO_results.csv"))

fc_vector <- res_df$log2FoldChange
names(fc_vector) <- res_df$symbol
fc_vector <- fc_vector[!is.na(names(fc_vector))]

png(file.path(go_dir, "GO_dotplot.png"), width = 9, height = 8, units = "in", res = 300)
print(dotplot(ego, showCategory = show_terms) +
        ggtitle(paste0("GO enrichment (", ontology, "): module ", module)) +
        theme_bw(base_size = 13) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              panel.grid = element_blank()))
invisible(dev.off())
msg_ok("GO dotplot saved to: ", file.path(go_dir, "GO_dotplot.png"))

# ---- barplot ----
tryCatch({
  png(file.path(go_dir, "GO_barplot.png"), width = 9, height = 8, units = "in", res = 300)
  print(barplot(ego, showCategory = show_terms) +
          ggtitle(paste0("GO enrichment (", ontology, "): module ", module)) +
          theme_bw(base_size = 13) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                panel.grid = element_blank()))
  invisible(dev.off())
  msg_ok("GO barplot saved to: ", file.path(go_dir, "GO_barplot.png"))
}, error = function(e) {
  invisible(dev.off())
  msg_step(paste0("Barplot skipped: ", conditionMessage(e)))
})

# ---- cnetplot (gene-term network) ----
tryCatch({
  png(file.path(go_dir, "GO_cnetplot.png"), width = 10, height = 9, units = "in", res = 300)
  print(cnetplot(ego, showCategory = min(show_terms, 5)) +
          ggtitle(paste0("Gene-term network (", ontology, "): module ", module)))
  invisible(dev.off())
  msg_ok("GO cnetplot saved to: ", file.path(go_dir, "GO_cnetplot.png"))
}, error = function(e) {
  invisible(dev.off())
  msg_step(paste0("Cnetplot skipped: ", conditionMessage(e)))
})

# ---- heatplot (gene-term membership grid) ----
tryCatch({
  png(file.path(go_dir, "GO_heatplot.png"), width = 12, height = 5, units = "in", res = 300)
  print(heatplot(ego, showCategory = show_terms, foldChange = fc_vector) +
          ggtitle(paste0("Gene-term heatmap (", ontology, "): module ", module)) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  invisible(dev.off())
  msg_ok("GO heatplot saved to: ", file.path(go_dir, "GO_heatplot.png"))
}, error = function(e) {
  invisible(dev.off())
  msg_step(paste0("Heatplot skipped: ", conditionMessage(e)))
})
