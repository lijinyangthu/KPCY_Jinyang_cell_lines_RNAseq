# KPCY_Jinyang_cell_lines
RNA-Seq analysis of Jinyang's clonal KPCY cell lines
David Balli - Stanger lab
balli.dave @ gmail

## Bioinformatic pipeline
1. Upenn HPC job script workflow
  1. Trim poor quality reads with seqtk trimmer, fastqc report of trimmed.fastq.gz files, and alignment/quantification of gene-level counts against mm10.
  - seqtk_trim.sh --> fastqc.sh --> star_alignment.sh

2. Analysis workflow
  1. Preprocessing of RNAseq counts, Consensus Clustering versus PDA types, Immune gene set enrichment, and RNAseq/FLOW correlations
  - STAR_preprocessing_normalization.R --> ConsesnsusClustering_PDAsubtyping.R --> Immune_GSVA_enrichment.R --> RNAseq_Flow_correlations.R

  2. generalized R functions for common analyses are included in 'functions.R' script
  - PCA_analysis(), pca_theme(), STAR_to_TPM(), gene_exon_length(), gene_set_statistic_test()

3. Results per 2016-11-17
  1. Consensus cluster results frozen and saved in 2016-11-15-annotaiton.xlsx/.csv
  2. Consensus cluster results for each signature stored in subfolders in results/
  3. Enrichment versus Immunome gene signature done (2016-11-15-Immunome-GSVA.pdf)
  4. Moffitt Stromal level analysis versus KEGG, Biocarta, Reactome, CGP-Immune genesets from broad
