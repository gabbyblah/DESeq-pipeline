# RNAseq-analysis-pipeline
A command-line RNA-seq analysis pipeline in R. Runs DESeq2 differential expression, then WGCNA co-expression network analysis, producing annotated tables and publication-style figures at each stage. 

Three scripts:
- **`deseq-pipeline.r`** - differential expression (DESeq2): results table + PCA plot, heatmap, & volcano plot
- **`wgcna.r`** - co-expression network (WGCNA): module detection, hub genes + soft-threshold, dendrogram, & module-trait figures
- **`pipeline.r`** - runs both scripts in sequence from raw counts matrix, then cleans up intermediate files

Works with any organism that has a Bioconductor org.*.eg.db annotation package - defaults to rhesus macaque (*Macaca mulatta*), fully configurable via command-line flags.

**Quick Start (full pipeline)**

Run differential expression *and* co-expression in one command:
```bash
Rscript pipeline.r <counts.csv> <metadata.csv> [options]
```
This runs `deseq-pipeline.r`, then `wgcna.r`, then removes the intermediate files. To run a single stage instead, call that script directly (see below). 

```bash
# full pipeline, all defaults (macaque, Ensembl IDs)
Rscript pipeline.r data/counts.csv data/meta.csv -o results/cohort1
```
```bash
# full pipeline, human data with custom column/reference and WGCNA tuning
Rscript pipeline.r data/counts.csv data/meta.csv -o results/human -d org.Hs.eg.db -c Treatment -r Untreated -m 20
```
  Note: `wgcna.r` reads the output of `deseq-pipeline.r`, so it must run *after* DESeq with the same `--outdir`. `pipeline.r` handles this ordering automatically. If a stage fails, the pipeline halts and keeps the intermediate files for debugging.


# DESeq-pipeline.r - Differential Expression
```bash
Rscript deseq-pipeline.r <counts.csv> <metadata.csv> [options]
```

**Positional Arguments**

| Position | Argument | Default | Description |
|---|---|---|---|
|1|counts.csv|data/counts.csv|Counts matrix|
|2|metadata.csv|data/meta.csv|Sample metadata|

**Flag Arguments**

| Flag | Default | Description |
|---|---|---|
|-o / --outdir|results/|Where outputs are written|
|-d / --database| org.Mmu.eg.db|Bioconductor annotation package for specified organism|
|-k / --keytype|ENSEMBL|Gene ID type used in counts file|
|-c / --condition|condition|Title of your grouping column in metadata file|
|-r / --reference|Control|Baseline to which all conditions will be compared|
|-f / --mincount|10|Minimum total count across samples to keep a gene|

**Outputs**

All written to the specified output directory with -o / --outdir (default: results/)

| Output | Description |
|---|---|
|DE_results.csv|Full results table, sorted by adjusted p-value with gene symbols |
|PCA.png|PCA of variance-stabilized counts, colored by condition|
|heatmap.png|Heatmap of the top 30 significant genes, row-scaled with condition annotation|
|volcano.png|Volcano plot with the top 15 significant named genes by fold change magnitude|

# Examples
```bash
# macaque, Ensembl IDs (all defaults)
Rscript deseq-pipeline.r data/counts.csv data/meta.csv 
```
```bash
# human data, Ensembl IDs, condition column titled 'Treatment' with reference level 'Untreated'
Rscript deseq-pipeline.r ~/counts_sample.csv ~/metadata.csv -o human -d org.Hs.eg.db -k ENSEMBL -c Treatment -r Untreated
```
```bash
# mouse data, Entrez IDs, condition column titled 'grouping' with reference level 'control'
Rscript deseq-pipeline.r data/counts.csv Downloads/meta.csv --database org.Mm.eg.db --keytype ENTREZID --condition grouping -r control
```

# WGCNA.r - Co-Expression Network Analysis
Builds a WGCNA co-expression network from the DESeq output, detects gene modules, identifies the module most correlated with condition, and reports its hub genes.

```bash
Rscript wgcna.r [options]
```

# Flag Arguments

|Flag|Default|Description|
|---|---|---|
|-o / --outdir|results/|Output directory (must contain DESeq output)|
|-c / --condition|condition|Metadata column defining the comparison|
|-p / --softpower|auto|Soft-thresholding power; auto-selects if unset|
|-n / --ngenes|5000|Number of top-variable genes to use|
|-t / --kmethreshold|(none)|kME cutoff for hub genes; uses --maxhubs if unset|
|-m / --maxhubs|30|Max hub genes to report|

# Outputs

|Output|Description|
|---|---|
|module_trait_heatmap.png|Correlation of each module with condition|
|module_dendrogram.png|Gene clustering with module colors|
|soft_threshold.png|Scale-free fit diagnostic for power selection|
|hub_genes.csv|Top hop genes of condition-associated module with DE stats + membership|
|module_assignments.csv|Every gene and its assigned module|

**Notes on soft power & hub selection**

- If --softpower is unset, the script uses WGCNA's automatic scale-free-fit estimate. If no power reaches the fit threshold (common with small sample sizes), it falls back to a signed-network default of 9. The soft_threshold.png figure documents the choice.
- Hub genes default to the top --maxhubs genes by module membership (kME). If --kmethreshold is set, hubs are instead genes above that kME cutoff, capped at --maxhubs.
- The condition-associated module is selected automatically as the module with the strongest correlation to condition (excluding the unassigned "grey" module).

# Pipeline.r - Full Pipeline Orchestrator
Runs `deseq-pipeline.r` then `wgcna.r`, and then removes intermediate files. Accepts all DESeq flags + WGCNA parameters, passed through to the appropriate stage. 

```bash
Rscript pipeline.r <counts.csv> <metadata.csv> [options]
```

|Flag|Default|Description|
|-o / --outdir|results/|Where outputs will be written|
|-d / --database|org.Mmu.eg.db|Bioconductor annotation package (DESeq stage)|
|-k / --keytype|ENSEMBL|Gene ID type (DESeq stage)|
|-c / --condition|condition|Title of your grouping column in metadata file (DESeq stage)|
|-r / --reference|Control|Baseline to which all conditions will be compared (DESeq stage)|
|-f / --mincount|10|Minimum total count across samples to keep a gene (DESeq stage)|
|-p / --softpower|auto|Soft-thresholding power; auto-selects if unset (WGCNA stage)|
|-n / --ngenes|5000|Number of top-variable genes to use (WGCNA stage)|
|-t / --kmethreshold|(none)|kME cutoff for hub genes; uses --maxhubs if unset (WGCNA stage)|
|-m / --maxhubs|30|Max hub genes to report (WGCNA stage)|

If any stage fails, the pipeline halts and keeps the intermediate files so you may debug or resume. Intermediate files are only removed after successful run.

# Input Formats
**Counts CSV** - genes as rows, samples as columns. First column carries gene ID (used as row names); first row is sample names.
|gene_id|C1|C2|C3|C4|D1|D2|D3|D4|
|---|---|---|---|---|---|---|---|---|
|ENSMMUG000000001|12|40|8|2|145|55|67|8|

**Metadata CSV** - samples as rows, one row per sample. First column is the sample ID (must match the column names).
|Samples|condition|
|---|---|
|C1|Control|
|C2|Control|
|C3|Control|
|C4|Control|
|D1|Disease|
|D2|Disease|
|D3|Disease|
|D4|Disease|

If your grouping is not named `condition`, set it with `-c`. This is **case-sensitive** - a column named `Condition` (capital C) requires `-c Condition`.

**Notes**
- Sample names in the counts columns exactly match the metadata sample IDs (case sensitive).
- The gene ID type in the counts file matches the keytype argument and is valid for the chosen annotation package.
- Gene symbol coverage varies by organism; for macaque it is incomplete, and many IDs will appear unannotated in the figures. Each figure only shows gene names that are annotated.
- Genes with total count below --mincount are filtered before analysis.
- The annotation database and counts ID type must be consistent - passing an annotation package for the wrong organism will error during symbol mapping.
- Passing a distinct directory per run keeps results from overwriting each other.
- The script installs any missing R packages automatically on first run (*requires an internet connection*)
- Arguments are case sensitive, especially -r & -c, so check your metadata file & ensure your cases match your arguments.
- Be sure to set -c if your metadata column for condition is named something else, i.e. 'group,' 'genotype,' 'treatment'.

**Scope**
- Single-factor designs only (one condition column). Multifactor models (i.e. controlling for batch) will be supported in the future.
- Plot titles and colors in the DESeq figures assume two groups. Other group counts still run but may use default colors.

# Requirements
R(>=4.0). The script auto-installs missing packages on the first run, but you can pre-install them:

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "org.Mmu.eg.db", "WGCNA", "SummarizedExperiment"))
install.packages(c("optparse", "ggplot2", "pheatmap", "RColorBrewer", "ggrepel"))
```
*You must install the matching annotation package for your samples (eg. org.Hs.eg.db for human, org.Mm.eg.db for mouse) and pass it as the org_db argument - annotationDbi packages do NOT auto-install. 

# Credits
This pipeline relies on the following open-source R/Bioconductor packages:
|Package|Description|
|---|---|
|DESeq2|Differential expression analysis|
|WGCNA|Weighted gene co-expression network analysis|
|AnnotationDbi + org.*.eg.db|Gene ID to symbol matching|
|ggplot2, ggrepel|PCA & volcano plots|
|pheatmap, RColorBrewer|Heatmap|
|optparse|Command-line argument parsing|

# Citation
If you use this pipeline, please cite the underlying methods!

  DESeq2 - Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)

  WGCNA - Langfelder, P., Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9:559.

# License
Released under MIT License - see LICENSE for details.
