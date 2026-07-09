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

msg_step <- function(...) cat("\033[36m", ..., "\033[0m\n", sep = "")
msg_ok   <- function(...) cat("\033[32m", ..., "\033[0m\n", sep = "")
red      <- function(x) paste0("\033[1;31m", x, "\033[0m")

ensure_pkg("optparse")

option_list <- list(
        make_option(c("-o", "--outdir"), type = "character", default = "results",
                    help = "Output directory (must contain DESeq pipeline output) [default: %default]"),
        make_option(c("-c", "--condition"), type = "character", default = "condition",
                    help = "Name of metadata column defining the comparison [default: %default]"),
        make_option(c("-p", "--softpower"), type = "integer", default = NULL,
                    help = "Soft-thresholding power. If unset, auto-selects (pickSoftThreshold estimate, else sample-size default)"),
        make_option(c("-n", "--ngenes"), type = "integer", default = 5000,
                    help = "Number of top-variable genes to use [default: %default]"),
        make_option(c("-t", "--kmethreshold"), type = "double", default = NULL,
                    help = "Module membership (kME) cutoff for hub genes. If unset, uses --maxhubs top genes"),
        make_option(c("-m", "--maxhubs"), type = "integer", default = 30,
                    help = "Max number of hub genes to report [default: %default]")
)

parser <- OptionParser(
        usage = "Rscript wgcna.r [options]",
        option_list = option_list,
        description = "\nWGCNA co-expression analysis. Reads DESeq-pipeline output from --outdir, builds a co-expression network, finds the module most correlated with condition, and reports hub genes and figures.",
        epilogue = "\nRun the DESeq pipeline first - this script reads its output from the same --outdir. \n\nExample:\n Rscript wgcna.r -o results/cohort1 -c condition\n"
)

opt <- parse_args(parser)

outdir        <- opt$outdir
condition_col <- opt$condition
soft_power    <- opt$softpower
n_genes       <- opt$ngenes
kme_threshold <- opt$kmethreshold
max_hubs      <- opt$maxhubs

intermediate <- file.path(outdir, "intermediate")

ensure_pkg("SummarizedExperiment", bioc = TRUE)
ensure_pkg("WGCNA", bioc = TRUE)
cor <- WGCNA::cor

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

if (!file.exists(file.path(intermediate, "vsdata.rds"))) {
  stop(red(paste0("Required input not found: ", file.path(intermediate, "vsdata.rds"), "\n",
                  "Run the DESeq pipeline first with the same --outdir.")), call. = FALSE)
}

msg_step("Starting WGCNA analysis...")

meta <- readRDS(file.path(intermediate, "meta.rds"))
vsdata <- readRDS(file.path(intermediate, "vsdata.rds"))
res_df <- read.csv(file.path(outdir, "DE_results.csv"), row.names = 1)
vst_mat <- assay(vsdata) 
vars <- apply(vst_mat, 1, var)
top <- names(sort(vars, decreasing = TRUE))[1:n_genes]
datExpr <- t(vst_mat[top, ])

if (!condition_col %in% colnames(meta)) {
  stop(red(paste0("Condition column '", condition_col, "' not found in metadata.\n",
                  "Available columns: ", paste(colnames(meta), collapse = ", "))),
       call. = FALSE)
}

gsg <- goodSamplesGenes(datExpr, verbose = 0)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)

if (!is.null(soft_power)) {
  softPower <- soft_power
  msg_step(paste("Using user-specified soft power:", softPower))
} else if (!is.na(sft$powerEstimate)) {
  softPower <- sft$powerEstimate
  msg_step(paste("Auto-selected soft power (scale-free fit):", softPower))
} else {
  softPower <- 9
  msg_step("No power reached scale-free fit; using fallback power 9")
}

png(file.path(outdir, "soft_threshold.png"), width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))   # two plots side by side
# left: scale-free fit vs power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)", ylab = "Scale-free topology fit (signed R^2)",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.8, col = "blue", lty = 2)   # the 0.8 fit target
# right: mean connectivity vs power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft threshold (power)", ylab = "Mean connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
invisible(dev.off())
msg_ok("Soft-threshold plot saved to: ", file.path(outdir, "soft_threshold.png"))

net <- blockwiseModules(datExpr,
                        power = softPower,
                        networkType = "signed",
                        TOMType = "signed",
                        minModuleSize = 30,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        saveTOMs = FALSE,
                        verbose = 0)
moduleColors <- labels2colors(net$colors)

# Figure 3: gene dendrogram with module colors
png(file.path(outdir, "module_dendrogram.png"), width = 20, height = 6, units = "in", res = 300)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene clustering dendrogram and module assignment"
)
invisible(dev.off())
msg_ok("Module dendrogram saved to: ", file.path(outdir, "module_dendrogram.png"))

MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs)
condition_num <- as.numeric(meta[[condition_col]]) - 1
moduleTraitCor <- cor(MEs, condition_num, use = "p")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nrow(datExpr))   # ADDED: p-values

msg_step("=== Module-condition correlation ===")
print(cbind(cor = moduleTraitCor, p = moduleTraitP))

# Figure 2: module-trait relationship heatmap
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitP, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)

png(file.path(outdir, "module_trait_heatmap.png"), width = 6, height = 10, units = "in", res = 300)
par(mar = c(6, 10, 3, 3))   # margins: bottom, left, top, right — room for labels
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = condition_col,
  yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "Module-condition relationships"
)
invisible(dev.off())
msg_ok("Module-trait heatmap saved to: ", file.path(outdir, "module_trait_heatmap.png"))

mtc <- moduleTraitCor[rownames(moduleTraitCor) != "MEgrey", , drop = FALSE]
mod <- sub("^ME", "", rownames(mtc)[which.max(abs(mtc[, 1]))])
msg_step("Auto-selected module of interest:", mod)

MM <- cor(datExpr, MEs[, paste0("ME", mod)], use = "p")
rownames(MM) <- colnames(datExpr)

# rank all genes in the module by module membership (kME)
if (!is.null(kme_threshold)) {
  # threshold mode: keep genes above the kME cutoff, capped at max_hubs
  above <- rownames(MM)[abs(MM[, 1]) >= kme_threshold]
  above <- above[order(abs(MM[above, 1]), decreasing = TRUE)]   # strongest first
  hub <- head(above, max_hubs)
  msg_step(paste0("Hub selection: kME >= ", kme_threshold,
                  " (", length(hub), " genes, capped at ", max_hubs, ")"))
} else {
  # default: top max_hubs by kME
  ranked <- rownames(MM)[order(abs(MM[, 1]), decreasing = TRUE)]
  hub <- head(ranked, max_hubs)
  msg_step(paste0("Hub selection: top ", length(hub), " genes by kME"))
}

if (length(hub) == 0) {
  msg_step(paste0("Warning: no genes passed kME threshold ", kme_threshold,
                  " in module '", mod, "'. hub_genes.csv will be empty."))
}

hub_df <- res_df[hub, c("symbol", "log2FoldChange", "padj")]
hub_df$moduleMembership <- MM[hub, 1]
write.csv(hub_df, file.path(outdir, "hub_genes.csv"))
msg_ok("Hub genes saved to:", file.path(outdir, "hub_genes.csv"))

module_df <- data.frame(gene_id = colnames(datExpr), module = moduleColors)
write.csv(module_df, file.path(outdir, "module_assignments.csv"), row.names = FALSE)
msg_ok("Module assignments saved to:", file.path(outdir, "module_assignments.csv"))

selected_genes <- colnames(datExpr)[moduleColors == mod]
background_genes <- colnames(datExpr)
saveRDS(list(selected = selected_genes, background = background_genes, module = mod), 
        file.path(intermediate, "go_input.rds"))
msg_ok("GO input saved. Selected module '", mod, "' has ", 
    length(selected_genes), " genes.")