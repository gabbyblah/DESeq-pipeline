args <- commandArgs(trailingOnly = TRUE) # Get command line arguments

if (length(args) > 0 && args[1] %in% c("-h", "--help")) {
  cat("
RNA-seq Differential Expression Analysis (DESeq2)
  
Usage: 
  Rscript deseq-pipeline.r <counts.csv> <metadata.csv> <output_dir> <org_db> <keytype>

Arguments (all positional):
  counts.csv:     Path to counts matrix: genes as rows, samples as columns        [default: data/counts.csv]
  metadata.csv:   Path to sample metadata with a 'condition' column               [default: data/meta.csv]
  output_dir:     Directory to save results & plots                               [default: /results/]
  org_db:         Bioconductor organism annotation package (eg: org.Mmu.eg.db)    [default: org.Mmu.eg.db]
  keytype:        Gene ID type in counts file (eg: ENSEMBLE, ENTRZID, SYMBOL)     [default: ENSEMBL]
    
Outputs (written to output_dir):
  DE_results.csv:     Annotated differential expression results
  PCA.png:            PCA plot of samples
  heatmap.png:        Heatmap of top 30 differentially expressed genes
  volcano.png:        Volcano plot of differential expression results
    
Examples:
  Rscript deseq-pipeline.r data/counts.csv data/meta.csv results/cohort1 
  Rscript deseq-pipeline.r Desktop/counts.csv downloads/meta.csv results/human org.Hs.eg/db ENSEMBL
")
  quit(save = "no")
}

counts_file <- if (length(args) >= 1) args[1] else "data/counts.csv"
meta_file <- if (length(args) >= 2) args[2] else "data/meta.csv"
outdir <- if (length(args) >= 3) args[3] else "results"
org_db <- if (length(args) >= 4) args[4] else "org.Mmu.eg.db"
keytype <- if (length(args) >= 5) args[5] else "ENSEMBL"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(counts_file)) stop(paste("Counts file does not exist:",
                                          counts_file))
if (!file.exists(meta_file)) stop(paste("Metadata file does not exist:",
                                        meta_file))

ensure_pkg <- function(pkg, bioc = FALSE) {
  if(!base::requireNamespace(pkg, quietly = TRUE)) {
    if (!base::requireNamespace("BiocManager", quietly = TRUE)) 
      utils::install.packages("BiocManager")
    if (bioc) BiocManager::install(pkg, update = FALSE, ask = FALSE)
    else utils::install.packages(pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }

ensure_pkg("BiocManager")
ensure_pkg("DESeq2", bioc = TRUE)
ensure_pkg(org_db, bioc = TRUE)
ensure_pkg("ggplot2")
ensure_pkg("pheatmap")
ensure_pkg("RColorBrewer")
ensure_pkg("ggrepel")

# Read counts & metadata
counts <- read.csv(counts_file, header = TRUE, row.names = 1) # Read counts data
meta <- read.csv(meta_file, header = TRUE, row.names = 1) # Read metadata
counts <- counts[, rownames(meta)] # Ensure counts columns match metadata rows
stopifnot(all(colnames(counts) == rownames(meta)))

meta$condition <- relevel(factor(meta$condition), ref = "Control")

counts <- counts[which(rowSums(counts) > 10), ]  # Filter out low count genes

# Create DESeq2 dataset and run analysis
cat("Starting differential expression analysis...\n")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta,
                              design = ~ condition)
dds <- suppressMessages(DESeq(dds)) 
cat("Differential expression analysis completed.\n")
cat("Extracting results for annotation...\n")
res <- results(dds) # Extract results from the DESeq analysis

res <- res[order(res$padj), ]

# Annotate results with gene symbols
cat("Annotating results...\n")
res_df <- as.data.frame(res) # Convert results to a data frame

orgdb_obj <- base::get(org_db)

if (!keytype %in% keytypes(orgdb_obj)) {
  stop(paste0("Keytype '", keytype, "' is not available for ", org_db,
             ". Available: ", paste(keytypes(orgdb_obj), collapse = ", ")))
}

res_df$symbol <- suppressMessages(mapIds(orgdb_obj,
                       keys = rownames(res_df),
                       keytype = keytype,
                       column = "SYMBOL",
                       multiVals = "first"))

res_df$label <- ifelse(is.na(res_df$symbol), rownames(res_df), res_df$symbol)
cat("Annotation completed.\n")

write.csv(as.data.frame(res_df), file = file.path(outdir, "DE_results.csv"))
cat("Annotated results saved to:", file.path(outdir, "DE_results.csv"), "\n")

# Generate PCA plot
cat("Generating PCA plot...\n  ")
vsdata <- vst(dds, blind = FALSE) # Variance stabilizing transformation
pca_plot <- plotPCA(vsdata, intgroup = "condition") + # Create PCA plot
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
cat("PCA plot saved to:", file.path(outdir, "PCA.png"), "\n")

# Generate heatmap of top differentially expressed genes
cat("Generating heatmap of top 30 differentially expressed genes...\n")
sig_genes <- rownames(res_df)[which(res_df$padj < 0.05)]
top_genes <- head(sig_genes, 30) # Select top 30 significant genes for heatmap

mat <- assay(vsdata)[top_genes, ] # Extract normalized counts for top genes

row_labels <- res_df[top_genes, "label"] # Get gene symbols for top genes

annotation_col <- data.frame(condition = meta$condition,
                             row.names = rownames(meta))

png(file.path(outdir, "heatmap.png"), width = 8, height = 10, units = "in",
    res = 300)
pheatmap(mat,
         labels_row = row_labels,
         annotation_col = annotation_col,
         scale = "row",
         show_colnames = TRUE,
         main = "Top 30 Differentially Expressed Genes")
dev.off()
cat("Heatmap saved to:", file.path(outdir, "heatmap.png"), "\n")

# Generate volcano plot
cat("Generating volcano plot...\n")

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
cat("Volcano plot saved to:", file.path(outdir, "volcano.png"), "\n")
