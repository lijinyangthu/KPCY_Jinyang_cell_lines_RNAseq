# STAR/featureCount quantification of mm10 + ERCC spike in sequences
# Jinyang's KPCY cell lines sorted and bulk tumors

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



#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
## Immunome signature from Bindea et al. versus BULK tumor gene signatures
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------

signature <- read.table("data/bindea_genesig.txt", 
                        header = TRUE, sep = "\t", 
                        col.names = c("cell_type", "hgnc_symbol"))
names(signature)
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

head(human_mouse_convesion$hgnc_symbol)
human_mouse_convesion[grep("IL2RA", human_mouse_convesion$hgnc_symbol),]

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

immunome_tpm <- tpm %>% filter(Geneid %in% immunome$Geneid) %>% select(matches("Geneid|BULK")) %>% data.frame()
row.names(immunome_tpm) <- immunome_tpm$Geneid
immunome_tpm$Geneid <- NULL

immunome_tpm <- log2(immunome_tpm + 1)
# isexpr <- rowSums(immunome_tpm) >= 1 ) >= ncol(immunome_tpm) * 0.3
isexpr <- rowSums(immunome_tpm) >= 1

# gene set variation analysis 
gsva_immunome <- gsva(as.matrix(immunome_tpm[isexpr,]), Immune_PDA_Signature, 
                      rnaseq = TRUE, verbose = TRUE, mx.diff = TRUE)$es.obs

pheatmap::pheatmap(gsva_immunome, color = colorRampPalette(c("navy", "white", "firebrick3"))(2345),
                   annotation_col = anno_df[, c("yfp_bulk", "Moffitt_Tumor_type", "Moffitt_Stromal_type")])

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

# convert MGI to paired hgnc_symbol
gsva_tpm <- left_join(tpm, human_mouse_convesion[, c("Geneid", "hgnc_symbol"), with = FALSE], by = "Geneid") %>% data.table()
gsva_tpm %<>% dplyr::select(matches("PD|hgnc_symbol")) %>% aggregate(. ~ hgnc_symbol, data = ., mean) %>% data.table()
gsva_tpm %>% dplyr::select(matches("PD")) %>% boxplot(log2(.  + 1), las = 2, outline = F)

gsva_tpm %<>% data.frame()
row.names(gsva_tpm) <- gsva_tpm$hgnc_symbol
gsva_tpm$hgnc_symbol <- NULL

# filter out lowly expressed genes and select BULK samples 
isexpr <- (rowSums(cpm(gsva_tpm) >= 1) >= 3)
gsva_tpm <- gsva_tpm[isexpr,] %>% dplyr::select(matches("BULK"))

#  GSVA versus biocarta, kegg, CGP, REACTOME, Immune-GSEA 
biocart_gsva <- gsva(as.matrix(gsva_tpm), biocart$genesets, rnaseq = TRUE, mx.diff = TRUE)$es.obs
kegg_gsva    <- gsva(as.matrix(gsva_tpm), kegg$genesets, rnaseq = TRUE, mx.diff = TRUE)$es.obs
cgp_gsva     <- gsva(as.matrix(gsva_tpm), cgp$genesets, rnaseq = TRUE, mx.diff = TRUE)$es.obs
react_gsva   <- gsva(as.matrix(gsva_tpm), react$genesets, rnaseq = TRUE, mx.diff = TRUE)$es.obs
immune_gsva   <- gsva(as.matrix(gsva_tpm), immune$genesets, rnaseq = TRUE, mx.diff = TRUE)$es.obs

# convert annotation data.table to data.frame for heatmap functions
annotation_df <- data.frame(annotation)
rownames(annotation_df) <- annotation_df$sample_id
annotation_df$sample_id <- NULL
pheatmap::pheatmap(biocart_gsva, scale = "row", 
                   clustering_method = "complete", 
                   annotation_col = annotation_df[,c(1,2,4)], 
                   clustering_distance_rows = "correlation", 
                   show_rownames = F, 
                   # file = paste0("Results/",Sys.Date(),"-GSVA-biocart-primary.pdf"),
                   width = 10, height = 8)

pheatmap::pheatmap(kegg_gsva, scale = "row", 
                   clustering_method = "complete", 
                   annotation_col = annotation_df[,c(1,2,4)], 
                   clustering_distance_rows = "correlation", 
                   show_rownames = F,
                   # file = paste0("Results/",Sys.Date(),"-GSVA-keggprimary.pdf"),
                   width = 10, height = 8)

pheatmap::pheatmap(cgp_gsva, scale = "row", 
                   clustering_method = "complete",
                   annotation_col = annotation_df[,c(1,2,4)], 
                   clustering_distance_rows = "correlation", 
                   show_rownames = F,
                   # file = paste0("Results/",Sys.Date(),"-GSVA-cgp-primary.pdf"),
                   width = 10, height = 8)

pheatmap::pheatmap(react_gsva, scale = "row", 
                   clustering_method = "complete", 
                   annotation_col = annotation_df[,c(1,2,4)], 
                   clustering_distance_rows = "correlation", 
                   show_rownames = F, 
                   # file = paste0("Results/",Sys.Date(),"-GSVA-react-primary.pdf"), 
                   width = 10, height = 8)

pheatmap::pheatmap(immune_gsva, scale = "row", 
                   clustering_method = "complete", 
                   annotation_col = annotation_df[,c(1,2,4)], 
                   clustering_distance_rows = "correlation", 
                   show_rownames = F, 
                   # file = paste0("Results/",Sys.Date(),"-GSVA-immune-primary.pdf"), 
                   width = 10, height = 8)


#-------------------------------------------------------------------------------------------------------------------------------
#
# statistical analysis of gene set enrichment 
#
#-------------------------------------------------------------------------------------------------------------------------------

# create design matrix and contrast matrix using moffitt tumor type definitions
p_m <- paste(annotation[match(names(gsva_tpm), annotation$sample_id), ]$moffitt_tumor_type, sep = '')
p_m <- factor(p_m)
p_m

pm_design <- model.matrix(~0+p_m)
colnames(pm_design) <- c("Basal_like", "Classical")
pm_design
ContrastMatrix <- limma::makeContrasts(Basal_like-Classical, levels = pm_design)

# function to run statistical tests of GSVA results against subtype level consensus clustering (found in src/functions.R)
# immune GSEA gene set
immune_res <- gene_set_statistic_test(gsva_result = immune_gsva, design_matrix = pm_design, contrast_matrix = ContrastMatrix) 
immune_res %<>% arrange(desc(logFC)) %>% data.table()
immune_res
rio::export(immune_res, file = paste0("Results/", Sys.Date(), "-GSVA-immune-bulk.csv"))

# kegg GSEA gene set
kegg_res <- gene_set_statistic_test(gsva_result = kegg_gsva, design_matrix = pm_design, contrast_matrix = ContrastMatrix) 
kegg_res %<>% arrange(desc(logFC)) %>% data.table()
kegg_res
rio::export(kegg_res, file = paste0("Results/", Sys.Date(), "-GSVA-kegg-bulk.csv"))

# biocarta GSEA gene set
biocart_res <- gene_set_statistic_test(gsva_result = biocart_gsva, design_matrix = pm_design, contrast_matrix = ContrastMatrix)
biocart_res %<>% arrange(desc(logFC)) %>% data.table()
biocart_res
rio::export(biocart_res, file = paste0("Results/", Sys.Date(), "-GSVA-biocarta-bulk.csv"))

# reactome GSEA gene set
reactome_res <- gene_set_statistic_test(gsva_result = react_gsva, design_matrix = pm_design, contrast_matrix = ContrastMatrix)
reactome_res %<>% arrange(desc(logFC)) %>% data.table()
reactome_res
rio::export(reactome_res, file = paste0("Results/", Sys.Date(), "-GSVA-reactome-bulk.csv"))






###
###
### everything below is saving for later reference
###
###


#-------------------------------------------------------------------------------------------------------------------------------
# heatmaps/boxplots of genes of interest
#-------------------------------------------------------------------------------------------------------------------------------

primary_names <- names(dgelist_counts)
gsva_anno_merge
# convert TPM to gene_id
tpm_plotting <- tpm %>% rownames_to_column("refseq_mrna") %>% data.table()
tpm_plotting <- left_join(tpm_plotting, gene_tx_info[,c("refseq_mrna", "external_gene_name")])
tpm_plotting <- tpm_plotting %>% 
  dplyr::select(contains("PD"), external_gene_name) %>% 
  na.omit() %>% 
  aggregate(. ~ external_gene_name, ., sum)
dim(tpm_plotting)

# boxplot of gene of interest
tpm_plotting %>% data.table() %>% 
  dplyr::filter(external_gene_name == "Prf1" | external_gene_name == "Gzma") %>% 
  dplyr::select(contains("PD")) %>% 
  melt(variable.name = "sample_id", value.name = "TPM") %>% 
  filter(sample_id %in% primary_names) %>% 
  left_join(., gsva_anno_merge[,1:2], by = "sample_id") %>% 
  
  # group_by(Moffitt_tumor_type) %>% summarise(mean_TPM = mean(TPM), stdev_TPM = sd(TPM))
  boxplot(log2(TPM + 1) ~ Moffitt_tumor_type, data = ., 
          ylim = c(0, 15), 
          las = 2, 
          col = c("cyan", "blue"))


#-------------------------------------------------------------------------------------------------------------------------------
# heatmaps of genes in react/cgp gene sets of interest
#-------------------------------------------------------------------------------------------------------------------------------
gsva_normalized <- gsva_count_norm$pseudo.counts %>% as.data.frame() %>%  rownames_to_column("Geneid") %>% data.table()

goi <- c("REACTOME_FGFR4_LIGAND_BINDING_AND_ACTIVATION", "REACTOME_SIGNALING_BY_FGFR3_MUTANTS", 
         "REACTOME_SIGNALING_BY_ACTIVATED_POINT_MUTANTS_OF_FGFR1", "REACTOME_FGFR2C_LIGAND_BINDING_AND_ACTIVATION",
         "REACTOME_FGFR1_LIGAND_BINDING_AND_ACTIVATION", 'REACTOME_FGFR_LIGAND_BINDING_AND_ACTIVATION',
         'REACTOME_REGULATION_OF_THE_FANCONI_ANEMIA_PATHWAY', "REACTOME_PD1_SIGNALING",
         "REACTOME_PHOSPHORYLATION_OF_CD3_AND_TCR_ZETA_CHAINS", "REACTOME_DOUBLE_STRAND_BREAK_REPAIR",
         "REACTOME_FANCONI_ANEMIA_PATHWAY", "REACTOME_EXTENSION_OF_TELOMERES", "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
         "REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_ACTIVITY_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS",
         "REACTOME_ACTIVATED_POINT_MUTANTS_OF_FGFR2", "REACTOME_CELL_EXTRACELLULAR_MATRIX_INTERACTIONS", "REACTOME_REGULATION_OF_GENE_EXPRESSION_IN_BETA_CELLS")

fgfr_genes_df <- (react$genesets[goi]) %>% unlist() %>% data.frame(., stringsAsFactors = FALSE)
fgfr_genes <- (fgfr_genes_df$.)

rownames_temp <- gsub('[[:digit:]]+','', row.names(fgfr_genes_df))
fgfr_genes_df$REACTOME <- gsub("REACTOME_*", "\\1", rownames_temp)
fgfr_genes_df <- fgfr_genes_df[fgfr_genes_df$. %in% row.names(temp),]
ords <- match(row.names(temp), fgfr_genes_df$.)
fgfr_genes_df <- fgfr_genes_df[ords,]
row.names(fgfr_genes_df) <- fgfr_genes_df$.
fgfr_genes_df$. <- NULL

temp <- gsva_normalized %>% data.table() %>% 
  dplyr::filter(Geneid %in% fgfr_genes) 
dim(temp)
row.names(temp) <- temp$Geneid

pheatmap::pheatmap(temp[, 2:ncol(temp)], show_colnames = FALSE, annotation_names_row = FALSE,  fontsize = 8, 
                   scale = "row", clustering_distance_rows = "correlation", 
                   annotation_col = pheat_anno[,c(1,3,4)], 
                   annotation_row = fgfr_genes_df,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(2345),
                   file = "Results/test.pdf", width = 15, height = 8)

