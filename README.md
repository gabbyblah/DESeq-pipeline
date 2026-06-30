# DESeq-pipeline
A single command-line R script that runs a DESeq2 differential expression analysis from a counts table &amp; sample metadata file and produces an annotated results table plus three figures (PCA, heatmap, volcano plot)

Works with any organism that has a Bioconductor org.*.db.eg annotation package - defaults to rhesus macaque (*Macaca mulatta*), but the annotation database & gene-ID type are configurable from the command line.

# What it Does
1. Reads a count matrix and sample metadata from CSV
2. Aligns samples by name & verifies match
3. Runs DESeq2 with condition as experimental variable using Control as reference level
4. Annotates gene IDs with gene symbols
5. Writes a results table & three figures to an output directory

The script installs any missing R packages automatically on first run (*requires an internet connection*)

# Outputs
All written to the specified output directory (default: results/)

| Output | Description |
|---|---|
|DE_results.csv|Full results table, sorted by adjusted p-value with gene symbols |
|PCA.png|PCA of variance-stabilized counts, colored by condition|
|heatmap.png|Heatmap of the top 30 significant genes, row-scaled with condition annotation|
|volcano.png|Volcano plot with the top 15 significant named genes by fold change magnitude|

# Usage
```bash
Rscript deseq-pipeline.r <counts.csv> <metadata.csv> <output_dir> <org_db> <keytype>
```

| Position | Argument | Default | Description |
|---|---|---|---|
|1|counts.csv|data/counts.csv|Counts matrix|
|2|metadata.csv|data/meta.csv|Sample metadata|
|3|output_dir|results/|Where outputs are written|
|4|org_db| org.Mmu.eg.db|Bioconductor annotation package for specified organism|
|5|keytype|ENSEMBL|Gene ID type used in counts file|

*Arugments are positional, so in order to set a later one you must supply all earlier arguments*

# Examples
```bash
# macaque, Ensembl IDs (all defaults)
Rscript deseq-pipeline.r data/counts.csv data/meta.csv results/cohort1
```
```bash
# human data, Ensembl IDs
Rscript deseq-pipeline.r ~/counts_sample.csv ~/metadata.csv results/human org.HS.eg.db EMSEMBL
```
```bash
# mouse data, Entrez IDs
Rscript deseq-pipeline.r data/counts.csv Downloads/meta.csv results/mouse org.Mm.eg.db ENTREZID
```

# Input Formats
Counts CSV  - genes as rows, samples as columns, first column has gene ID used as row names & first row is sample names:

|gene_id|C1|C2|C3|C4|D1|D2|D3|D4|
|---|---|---|---|---|---|---|---|---|
|ENSEMMUG00000000001|12|40|8|2|145|55|67|8|

Metadata CSV -  samples as rows, one row per sample, first column is sample ID (must match counts column names); must include a column named *condition*:
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

# Additional Notes
- The metadata file must have a column literally named *condition*
- That condition column has a level literally named *Control*, used as the baseline (so log2foldchange is relative to Control)
- Sample names in the counts columns exactly match the metadata sample IDs (case sensitive)
- The gene ID type in the counts file matches the keytype argument and is valid for the chosen annotation package
- Gene symbol coverage varies by organism; for macaque it is incomplete, and many IDs will appear unannotated in the figures. Each figure only shows gene names that are annotated.
- Low-count genes (row sum <= 10) are filtered before analysis
- The annotation database and counts ID type must be consistent - passing an annotation package for the wrong organism will error during symbol mapping.
- Passing a distinct directory per run keeps results from overwriting each other

*To use a different baseline or condition column, edit the relevel() line & design formula in the script*

# Requirements
R(>=4.0). The script auto-installs missing packages on the first run, but you can pre-install them:
```r
install.packages("BiocManager)
BiocManager::install(c("DESeq2", "org.*.eg.db"))
install.packages(c("ggplot2", "pheatmap", "RColorBrewer", "ggrepel"))
```
Install the matching annotation package for your samples (eg. org.Hs.eg.db for human, org.Mm.eg.db for mouse) and pass it as the org_db argument.

# Credits
This pipeline relies on the following open-source R/Bioconductor packages:
|Package|Description|
|---|---|
|DESeq2|differential expression analysis|
|AnnotationDbi + org.*.eg.db|gene ID to symbol matching|
|ggplot2, ggrepel|PCA & volcano plots|
|pheatmap, RColorBrewer|heatmap|

# Citation
If you use this pipeline, please cite the underlying DESeq2 method!

  Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and
  dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)

# License
Released under MIT License - see LICENSE for details.
