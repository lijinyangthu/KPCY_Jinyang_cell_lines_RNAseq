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

## Immunome signature from Bindea et al. versus BULK tumor gene signatures

signature <- read.table("/Users/dballi/Projects/PDA_TCGA_UNC_RNAseq_Immune/data/bindea_genesig.txt", 
                        header = TRUE, sep = "\t", col.names = c("cell_type", "hgnc_symbol"))
table(signature$cell_type)

immunome <- right_join(human_mouse_convesion, signature, by = "hgnc_symbol") %>%
    dplyr::select(Geneid, hgnc_symbol, cell_type) %>%
    unique() %>% na.omit() %>% data.table()

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
                             Th17_cells = c("BTLA", "CD200", "CD99"),
                             Th2_cells = immunome[immunome$cell_type == "Th2_cells", ]$Geneid,
                             Treg = immunome[immunome$cell_type == "TReg", ]$Geneid)

immunome_tpm <- tpm %>% filter(Geneid %in% immunome$Geneid) %>% select(matches("Geneid|BULK")) %>% data.frame()
row.names(immunome_tpm) <- immunome_tpm$Geneid
immunome_tpm$Geneid <- NULL



immunome_tpm <- log2(immunome_tpm + 1)
# isexpr <- rowSums(immunome_tpm) >= 1 ) >= ncol(immunome_tpm) * 0.3
isexpr <- rowSums(immunome_tpm) >= 1

gsva_immunome <- gsva(as.matrix(immunome_tpm[isexpr,]), Immune_PDA_Signature, 
                      rnaseq = TRUE, verbose = TRUE, mx.diff = TRUE)$es.obs

pheatmap::pheatmap(gsva_immunome, 
                   annotation_col = anno_df[, c("yfp_bulk", "Moffitt_Tumor_type", "Moffitt_Stromal_type")])
