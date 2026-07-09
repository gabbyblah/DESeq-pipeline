# Version 3
options(repos = c(CRAN = "https://cloud.r-project.org"))
ensure_pkg <- function(pkg, bioc = FALSE) {
  if (!base::requireNamespace(pkg, quietly = TRUE)) {
    if (!base::requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager")
    }
    if (bioc) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      utils::install.packages(pkg)
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

msg_step <- function(...) cat("\033[36m", ..., "\033[0m\n", sep = "")      # cyan: step in progress
msg_ok   <- function(...) cat("\033[32m", ..., "\033[0m\n", sep = "")      # green: success
red      <- function(x) paste0("\033[1;31m", x, "\033[0m") 

ensure_pkg("optparse")

option_list <- list(
  make_option(c("-o", "--outdir"), type = "character", default = "results",
              help = "Output directory for results & plots [default: %default]"),
  make_option(c("-d", "--database"), type = "character", default = "org.Mmu.eg.db",
              help = "Bioconductor annotation package for your data's organism [default: %default]"),
  make_option(c("-k", "--keytype"), type = "character", default = "ENSEMBL",
              help = "Gene ID type in your counts file [default: %default]"),
  make_option(c("-c", "--condition"), type = "character", default = "condition",
              help = "Name of metadata column for grouping [default: %default]"),
  make_option(c("-r", "--reference"), type = "character", default = "Control", # can we make this not case sensitive? if not we need to emphasize everything is case sensitive & address in errors
              help = "Reference (baseline) level within the condition column [default: %default]"),
  make_option(c("-f", "--mincount"), type = "integer", default = 10,
            help = "Minimum total count across samples to keep a gene [default: %default]")
)

parser <- OptionParser(
  usage = "Rscript deseq-pipeline.r <counts.csv> <metadata.csv> [options]",
  option_list = option_list,
  description = "\nRNA-seq differential expression analysis (DESeq2). Produces an annotated results table, PCA, heatmap, and volcano plot.",
  epilogue = "\nExample:\n Rscript deseq-pipeline.r data/counts.csv data/meta.csv -d org.Hs.eg.db -c grouping -r Untreated\n"
)

opt <- parse_args(parser, positional_arguments = 2)

counts_file    <- opt$args[1]
meta_file      <- opt$args[2]
outdir         <- opt$options$outdir
org_db         <- opt$options$database
keytype        <- opt$options$keytype
condition_col  <- opt$options$condition
reference      <- opt$options$reference
min_count      <- opt$options$mincount

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
intermediate <- file.path(outdir, "intermediate")
dir.create(intermediate, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(counts_file)) stop(red(paste("Counts file does not exist:",
                                          counts_file)))
if (!file.exists(meta_file)) stop(red(paste("Metadata file does not exist:",
                                        meta_file)))

ensure_pkg("BiocManager")
ensure_pkg("DESeq2", bioc = TRUE)
ensure_pkg(org_db, bioc = TRUE)
ensure_pkg("ggplot2")
ensure_pkg("pheatmap")
ensure_pkg("RColorBrewer")
ensure_pkg("ggrepel")

counts <- read.csv(counts_file, header = TRUE, row.names = 1)
meta <- read.csv(meta_file, header = TRUE, row.names = 1)

if (!condition_col %in% colnames(meta)) {
  stop(red(paste0("Condition column '", condition_col, "' not found in metadata.\n",
              "Available columns: ", paste(colnames(meta), collapse = ", "))),
       call. = FALSE)
}

# validate the reference level exists in that column
if (!reference %in% meta[[condition_col]]) {
  stop(red(paste0("Reference level '", reference, "' not found in column '", condition_col, "'.\n",
              "Values present: ", paste(unique(meta[[condition_col]]), collapse = ", "),
              "\n(Note: matching is case-sensitive.)")),
       call. = FALSE)
}

counts <- counts[, rownames(meta)]
stopifnot(all(colnames(counts) == rownames(meta)))

meta[[condition_col]] <- relevel(factor(meta[[condition_col]]), ref = reference)

counts <- counts[which(rowSums(counts) >= min_count), ]

msg_step("Starting differential expression analysis...")
design_formula <- as.formula(paste("~", condition_col))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design_formula)
dds <- suppressMessages(DESeq(dds)) 
msg_ok("Differential expression analysis completed.")
msg_step("Extracting results for annotation...")
res <- results(dds) 

res <- res[order(res$padj), ]

# Annotate results with gene symbols
msg_step("Annotating results...")
res_df <- as.data.frame(res) # Convert results to a data frame

orgdb_obj <- base::get(org_db)

if (!keytype %in% keytypes(orgdb_obj)) {
  stop(red(paste0("Keytype '", keytype, "' is not available for ", org_db,
             ". Available: ", paste(keytypes(orgdb_obj), collapse = ", "))))
}

res_df$symbol <- tryCatch(
  suppressMessages(mapIds(orgdb_obj,
                       keys = rownames(res_df),
                       keytype = keytype,
                       column = "SYMBOL",
                       multiVals = "first")),
  error = function(e) {
    stop(red(paste0(
      "\nGene ID annotation failed.\n",
      "This usually means the annotation database (--db ", org_db,
      ") or the keytype (--keytype ", keytype, ")\n",
      "does not match the gene IDs in your counts file.\n\n",
      "Your first few gene IDs look like: ",
      paste(head(rownames(res_df), 3), collapse = ", "), "\n",
      "Check that --db is the right organism and --keytype matches this ID format.\n",
      "Original error: ", conditionMessage(e)
    )), call. = FALSE)
  }
)

res_df$label <- ifelse(is.na(res_df$symbol), rownames(res_df), res_df$symbol)
msg_step("Annotation completed.")

write.csv(as.data.frame(res_df), file = file.path(outdir, "DE_results.csv"))
msg_ok("Annotated results saved to: ", file.path(outdir, "DE_results.csv"))

# Generate PCA plot
msg_step("Generating PCA plot...")
vsdata <- vst(dds, blind = FALSE) # Variance stabilizing transformation
pca_plot <- suppressMessages(plotPCA(vsdata, intgroup = condition_col)) + # Create PCA plot
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = name), size = 3) +
  theme_bw(base_size = 14) + # set specific theme
  ggtitle("PCA: Disease vs Control") + # title
  scale_color_manual(values = c("Control" = "#2c7fb8",
                                "Disease" = "#d95f0e")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        panel.grid = element_blank()) # remove grid lines

ggsave(file.path(outdir, "PCA.png"), pca_plot, width = 7, height = 6, dpi = 300)
msg_ok("PCA plot saved to: ", file.path(outdir, "PCA.png"))

# Generate heatmap of top differentially expressed genes
msg_step("Generating heatmap of top 30 differentially expressed genes...")
sig_genes <- rownames(res_df)[which(res_df$padj < 0.05)]
top_genes <- head(sig_genes, 30) # Select top 30 significant genes for heatmap

mat <- assay(vsdata)[top_genes, ] # Extract normalized counts for top genes

row_labels <- res_df[top_genes, "label"] # Get gene symbols for top genes

annotation_col <- data.frame(condition = meta[[condition_col]], row.names = rownames(meta))

png(file.path(outdir, "heatmap.png"), width = 8, height = 10, units = "in",
    res = 300)
pheatmap(mat,
         labels_row = row_labels,
         annotation_col = annotation_col,
         scale = "row",
         show_colnames = TRUE,
         main = "Top 30 Differentially Expressed Genes")
invisible(dev.off())
msg_ok("Heatmap saved to: ", file.path(outdir, "heatmap.png"))

# Generate volcano plot
msg_step("Generating volcano plot...")

res_df$diffexpressed <- "NS" # not significant
res_df$diffexpressed[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <-
  "Upregulated"
res_df$diffexpressed[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <-
  "Downregulated"

# pick which genes to label on the volcano plot (top 15 up and down regulated)
sig <- res_df[which(res_df$padj < 0.05), ] # significant only - which drops NAs
sig <- sig[order(abs(sig$log2FoldChange), decreasing = TRUE), ]
sig <- sig[!is.na(sig$symbol), ] # remove NAs in symbol column
label_genes <- head(sig, 15) # top 15 genes to label

# build the volcano plot
volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj),
                              color = diffexpressed)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_text_repel(data = label_genes, aes(label = label), size = 3,
                  max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = c("Upregulated" = "#d95f0e",
                                "Downregulated" = "#2c7fb8", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  theme_bw(base_size = 14) +
  ggtitle("Volcano: Disease vs Control") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank())

suppressWarnings(ggsave(file.path(outdir, "volcano.png"), volcano, width = 8,
       height = 7, dpi = 300))
msg_ok("Volcano plot saved to: ", file.path(outdir, "volcano.png"))

saveRDS(vsdata, file.path(intermediate, "vsdata.rds"))
saveRDS(meta, file.path(intermediate, "meta.rds"))
