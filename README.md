# nash-transcriptomics-analysis
RNA-seq differential expression and pathway enrichment analysis in NASH.

This repository contains a bioinformatics workflow for the analysis of RNA-seq data from NASH samples:
- SRA download and FASTQ format generation;
- Transcript quantification using Salmon;
- PCA analysis;
- Differential expression analysis with DESeq2;
- Volcano plot visualization;
- Gene Set Enrichment Analysis (GSEA);
- Over-Representation Analysis (ORA);
- KEGG pathway overlap analysis.

# Dataset
GEO accession: GSE126848
Organism: Homo sapiens

# Repository structure
scripts/
metadata/
results/
plots/

# Scripts
01_download_and_quantification.sh: download SRA files and perform transcript quantification using Salmon.

01_import_quantification_and_pca.R: import Salmon quantification files and perform PCA analysis.

02_differential_expression_and_volcano.R: perform differential expression analysis using DESeq2 and generate volcano plots.

03_enrichment_and_kegg_analysis.R: perform GSEA, ORA, and KEGG pathway analyses.

# Required annotation files
The following files should be manually downloaded and placed inside the metadata/ directory before running enrichment analyses:
c5.go.bp.v2024.1.Hs.symbols.gmt
Human.GRCh38.p13.annot.tsv

# Main R packages
tximport
DESeq2
biomaRt
clusterProfiler
enrichplot
ggplot2
KEGGREST

# Output
The pipeline generates:
- PCA plots
- Volcano plots
- Differentially expressed gene tables
- GSEA enrichment plots
- ORA pathway analyses
- KEGG overlap analyses

# Author

Lucrezia Banci
