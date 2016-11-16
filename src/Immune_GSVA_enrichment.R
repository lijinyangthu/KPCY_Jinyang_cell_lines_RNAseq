#----------------------------------------------------------------------------------------------------------------------------------------------------------
# script 3 for analysis of KPCY cell lines
# gene set variation analysis of Immunome gene signature + other selected populations from Bindea et al
# references for both gene sets taken from publications:
# Bindea manuscript: Cell, Spatiotemporal Dynamics of Intratumoral Immune Cells Reveal the Immune Landscape in Human Cancer, 2013
#----------------------------------------------------------------------------------------------------------------------------------------------------------

library("tidyverse") # loads ggplot2, tibble, tidyr, readr, purr, dplyr
library("genefilter")
library("edgeR")
library("ggrepel")
library("gplots")
library("reshape2")
library("data.table")
library("magrittr")
library("GSVA")
library("GSA")
library("limma")
library("edgeR")
source("src/functions.R")
library("rio")

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
# set graphics parameters
op <- par(mar = c(8, 5, 4, 2)+ 0.1)
options(op)
# set seed for reproducibility
set.seed(123)
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
#
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
## Immunome signature from Bindea et al. versus BULK tumor gene signatures
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------
signature <- read.table("data/bindea_Immunome_Signature.txt",
                        header = TRUE, sep = "\t",
                        col.names = c("cell_type", "hgnc_symbol"))
names(signature)
# adding MDSCs, Th17, Cytolytic index (Rooney Cell, 2015), Immunsuppression index (Balli in submission), and Treg cells
mdsc <- data.frame(cell_type = "MDSC", hgnc_symbol = c("CD11B", "CD33", "CD14", "CD15"))
Th17_cells <- data.frame(cell_type = "Th17_cells", hgnc_symbol = c("BTLA", "CD200", "CD99"))
cyt_index <- data.frame(cell_type = "CYT_index", hgnc_symbol = c("PRF1", "GZMA"))
immune_suppresion_index <- data.frame(cell_type = "Immune_suppression_index", hgnc_symbol =  c("CTLA4", "PDCD1", "CD274", "PDCD1LG2", "IDO1",
                                                                "IDO2", "LAG3", "ADORA2A", "HAVCR2", "TIGIT", "VTCN1", "C10orf54"))
treg <- data.frame(cell_type = "TReg", hgnc_symbol = c("IL2RA", "FOXP3"))

signature <- rbind(signature, mdsc, treg, Th17_cells, cyt_index, immune_suppresion_index)

immunome <- right_join(human_mouse_convesion, signature, by = "hgnc_symbol") %>%
    dplyr::select(Geneid, hgnc_symbol, cell_type) %>%
    unique() %>% na.omit() %>% data.table()

# using Bindea's gene signature
# but adding markers for MDSCs (CD11b (Itgam), CD33, CD14, CD15 (Fut4))
# adding CYT index genes (GZMA, PRF1)
# adding Immune suppression index genes ()

# add each cell type gene list into an R list for GSVA input
Immune_PDA_Signature <- list(aDC = immunome[immunome$cell_type == "aDC", ]$Geneid,
                             B_cells = immunome[immunome$cell_type == "B_cells", ]$Geneid,
                             blood_vessels = immunome[immunome$cell_type == "Blood_vessels", ]$Geneid,
                             CD8_T_cells = immunome[immunome$cell_type == "CD8_t_cells", ]$Geneid,
                             Cytotoxic_cells = immunome[immunome$cell_type == "Cytotoxic_cells", ]$Geneid,
                             DC = immunome[immunome$cell_type == "DC", ]$Geneid,
                             Eosinophils = immunome[immunome$cell_type == "Eosinophils", ]$Geneid,
                             iDC = immunome[immunome$cell_type == "iDC", ]$Geneid,
                             Lymph_vessels = immunome[immunome$cell_type == "Lymph_vessels", ]$Geneid,
                             Macrophages = immunome[immunome$cell_type == "Macrophages", ]$Geneid,
                             Mast_cells = immunome[immunome$cell_type == "Mast_cells", ]$Geneid,
                             Neutrophils = immunome[immunome$cell_type == "Neutrophils", ]$Geneid,
                             NK_CD56bright_cells = immunome[immunome$cell_type == "NK_CD56bright_cells", ]$Geneid,
                             NK_CD56dim_cells = immunome[immunome$cell_type == "NK_CD56dim_cells", ]$Geneid,
                             NK_cells = immunome[immunome$cell_type == "Nk_cells", ]$Geneid,
                             helper_T_cells = immunome[immunome$cell_type == "T_helper_cells", ]$Geneid,
                             Tcm = immunome[immunome$cell_type == "Tcm", ]$Geneid,
                             Tem = immunome[immunome$cell_type == "Tem", ]$Geneid,
                             Tfh = immunome[immunome$cell_type == "TFH", ]$Geneid,
                             Tgd = immunome[immunome$cell_type == "Tgd", ]$Geneid,
                             Th1_cells = immunome[immunome$cell_type == "Th1_cells", ]$Geneid,
                             Th17_cells = immunome[immunome$cell_type == "Th17_cells", ]$Geneid,
                             Th2_cells = immunome[immunome$cell_type == "Th2_cells", ]$Geneid,
                             Treg = immunome[immunome$cell_type == "TReg", ]$Geneid,
                             MDSC = immunome[immunome$cell_type == "MDSC", ]$Geneid,
                             CYT_index = immunome[immunome$cell_type == "CYT_index", ]$Geneid,
                             Immune_suppression_index = immunome[immunome$cell_type == "Immune_suppression_index", ]$Geneid)

# filter TPM values for Geneid and samples labeled as BULK
immunome_tpm <- tpm %>% filter(Geneid %in% immunome$Geneid) %>% dplyr::select(matches("Geneid|BULK")) %>% data.frame()
row.names(immunome_tpm) <- immunome_tpm$Geneid
immunome_tpm$Geneid <- NULL

# isexpr <- rowSums(immunome_tpm) >= 1 ) >= ncol(immunome_tpm) * 0.3
isexpr <- rowSums(immunome_tpm) >= 1

immunome_tpm <- log2(immunome_tpm[isexpr,] + 1)
dim(immunome_tpm)

# gene set variation analysis
gsva_immunome <- gsva(as.matrix(immunome_tpm), Immune_PDA_Signature,
                      rnaseq = TRUE, verbose = TRUE, mx.diff = TRUE)$es.obs

# heatmap of same genes as in PCA
anno_df <- data.frame(annotation, row.names = annotation[,1])

pheatmap::pheatmap(gsva_immunome, color = colorRampPalette(c("navy", "white", "firebrick3"))(2345),
                   annotation_col = anno_df[, c("moffitt_tumor_type", "bailey_type","moffitt_stromal_type")],
                   main = "Immunome Gene Set Enrichment",
                   file = "results/2016-11-15-Immunome-GSVA.pdf",
                   width = 12, height = 8)

#-------------------------------------------------------------------------------------------------------------------------------
#
# GSVA enrichment with Broad datasets
#
#-------------------------------------------------------------------------------------------------------------------------------
biocart <- GSA.read.gmt("data/c2.cp.biocarta.v5.1.symbols.gmt.txt")
names(biocart$genesets) <- biocart$geneset.names

react <- GSA.read.gmt("data/c2.cp.reactome.v5.1.symbols.gmt.txt")
names(react$genesets) <- react$geneset.names

kegg <- GSA.read.gmt("data/c2.cp.kegg.v5.1.symbols.gmt.txt")
names(kegg$genesets) <- kegg$geneset.names

cgp <- GSA.read.gmt("data/c2.cgp.v5.1.symbols.gmt.txt")
names(cgp$genesets) <- cgp$geneset.names

immune <- GSA.read.gmt("data/c7.all.v5.2.symbols.gmt")
names(immune$genesets) <- immune$geneset.names

# convert MGI to paired hgnc_symbol and aggregate to gene level TPM values
gsva_tpm <- left_join(tpm, human_mouse_convesion[, c("Geneid", "hgnc_symbol"), with = FALSE], by = "Geneid") %>% data.table()
gsva_tpm %<>% dplyr::select(matches("PD|hgnc_symbol")) %>% aggregate(. ~ hgnc_symbol, data = ., mean) %>% data.table()

gsva_tpm %<>% data.frame()
row.names(gsva_tpm) <- gsva_tpm$hgnc_symbol
gsva_tpm$hgnc_symbol <- NULL

# filter out lowly expressed genes and select BULK samples
isexpr <- rowSums(cpm(gsva_tpm) >= 1) >= 3
gsva_tpm <- gsva_tpm[isexpr,] %>% dplyr::select(matches("BULK"))
dim(gsva_tpm)

#  GSVA versus biocarta, kegg, CGP, REACTOME, Immune-GSEA
biocart_gsva   <- gsva(as.matrix(gsva_tpm), biocart$genesets, rnaseq = TRUE, mx.diff = TRUE)$es.obs
kegg_gsva      <- gsva(as.matrix(gsva_tpm), kegg$genesets, rnaseq = TRUE, mx.diff = TRUE)$es.obs
cgp_gsva       <- gsva(as.matrix(gsva_tpm), cgp$genesets, rnaseq = TRUE, mx.diff = TRUE)$es.obs
react_gsva     <- gsva(as.matrix(gsva_tpm), react$genesets, rnaseq = TRUE, mx.diff = TRUE)$es.obs
immune_gsva    <- gsva(as.matrix(gsva_tpm), immune$genesets, rnaseq = TRUE, mx.diff = TRUE)$es.obs

# convert annotation data.table to data.frame for heatmap functions
annotation_df <- data.frame(annotation)
rownames(annotation_df) <- annotation_df$sample_id
annotation_df$sample_id <- NULL

pheatmap::pheatmap(biocart_gsva, scale = "row",
                   clustering_method = "complete",
                   annotation_col = annotation_df[,c(4:6)],
                   clustering_distance_rows = "correlation",
                   show_rownames = F,
                   file = paste0("Results/2015-11-08-GSVA-biocart-primary.pdf"),
                   width = 10, height = 8)

pheatmap::pheatmap(kegg_gsva, scale = "row",
                   clustering_method = "complete",
                   annotation_col = annotation_df[,c(4:6)],
                   clustering_distance_rows = "correlation",
                   show_rownames = F,
                   file = paste0("Results/2016-11-08-GSVA-keggprimary.pdf"),
                   width = 10, height = 8)

pheatmap::pheatmap(cgp_gsva, scale = "row",
                   clustering_method = "complete",
                   annotation_col = annotation_df[,c(4:6)],
                   clustering_distance_rows = "correlation",
                   show_rownames = F,
                   file = paste0("Results/2016-11-08-GSVA-cgp-primary.pdf"),
                   width = 10, height = 8)

pheatmap::pheatmap(react_gsva, scale = "row",
                   clustering_method = "complete",
                   annotation_col = annotation_df[,c(4:6)],
                   clustering_distance_rows = "correlation",
                   show_rownames = F,
                   file = paste0("Results/2016-11-08-GSVA-react-primary.pdf"),
                   width = 10, height = 8)

pheatmap::pheatmap(immune_gsva, scale = "row",
                   clustering_method = "complete",
                   annotation_col = annotation_df[,c(4:6)],
                   clustering_distance_rows = "correlation",
                   show_rownames = F,
                   file = paste0("Results/2016-11-08-GSVA-immune-primary.pdf"),
                   width = 10, height = 8)

#-------------------------------------------------------------------------------------------------------------------------------
#
# statistical analysis of gene set enrichment
#
# working through what comparisions to make
#
#-------------------------------------------------------------------------------------------------------------------------------

# create design matrix and contrast matrix using moffitt stromal type definitions
p_m <- paste(annotation[match(names(gsva_tpm), annotation$sample_id), ]$moffitt_stromal_type, sep = '')
p_m <- factor(p_m)
pm_design <- model.matrix(~0+p_m)
colnames(pm_design) <- c("Activated","Normal")
pm_design
ContrastMatrix <- limma::makeContrasts(Activated-Normal, levels = pm_design)

# function to run statistical tests of GSVA results against subtype level consensus clustering (found in src/functions.R)
# immune GSEA gene set
immune_res <- gene_set_statistic_test(gsva_result = immune_gsva, design_matrix = pm_design, contrast_matrix = ContrastMatrix)
immune_res %<>% arrange(desc(logFC)) %>% data.table()
immune_res
rio::export(immune_res, file = paste0("Results/", Sys.Date(), "-GSVA-immune-bulk-stroma.csv"))

# kegg GSEA gene set
kegg_res <- gene_set_statistic_test(gsva_result = kegg_gsva, design_matrix = pm_design, contrast_matrix = ContrastMatrix)
kegg_res %<>% arrange(desc(logFC)) %>% data.table()
kegg_res
rio::export(kegg_res, file = paste0("Results/", Sys.Date(), "-GSVA-kegg-bulk-stroma.csv"))

# biocarta GSEA gene set
biocart_res <- gene_set_statistic_test(gsva_result = biocart_gsva, design_matrix = pm_design, contrast_matrix = ContrastMatrix)
biocart_res %<>% arrange(desc(logFC)) %>% data.table()
biocart_res
rio::export(biocart_res, file = paste0("Results/", Sys.Date(), "-GSVA-biocarta-bulk-stroma.csv"))

# reactome GSEA gene set
reactome_res <- gene_set_statistic_test(gsva_result = react_gsva, design_matrix = pm_design, contrast_matrix = ContrastMatrix)
reactome_res %<>% arrange(desc(logFC)) %>% data.table()
reactome_res
rio::export(reactome_res, file = paste0("Results/", Sys.Date(), "-GSVA-reactome-bulk-stroma.csv"))

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
# data should be ready for further analysis
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------