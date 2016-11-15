# Correlation of Cell type gene set enrichment versus flow immunophenotyping
# Jinyang's KPCY cell lines sorted and bulk tumors

library("tidyverse") # loads ggplot2, tibble, tidyr, readr, purr, dplyr
library("genefilter")
library("edgeR")
library("ggrepel")
library("gplots")
library("reshape2")
library("data.table")
library("stringr")
library("magrittr")
source("src/functions.R")
library("rio")

flowdata <- fread("data/20160430_Leukocyte_immune 2.csv")
corrected_names <- paste0("PD", names(flowdata[,2:ncol(flowdata), with = FALSE]))
corrected_names <- str_split_fixed(corrected_names, "c", 2)
corrected_names <- paste0(corrected_names[,1], "_C", corrected_names[,2], "_BULK")

# transpose, convert ID names, and merge with annotation file
flowdata <- data.table::transpose(flowdata)
flowdata <- flowdata[-1,]
setnames(flowdata, c("flow_CD4_T_cells", "flow_CD8_T_cells", "flow_gMDSCs", "flow_DCs", "flow_CD4_CD8_null_T_cells", "flow_F480_macrophages", "flow_B_cells"))
flowdata <- flowdata[, sample_id := corrected_names]
rio::export(flowdata, file = "data/2016-11-15-processed-flowdata.csv")

# transpose GSVA results, fix names and merge
gsva_names <- gsva_immunome %>% data.frame() %>% rownames_to_column("sample_id") %>% names()
gsva_for_flow <- gsva_immunome %>% data.frame() %>% rownames_to_column("") %>% data.table::transpose() %>% data.table()
gsva_for_flow[, sample_id := gsva_names]
colnames(gsva_for_flow) <- as.character(gsva_for_flow[1,])
colnames(gsva_for_flow) <- c(paste0("gsva_", colnames(gsva_for_flow)[1:27]), "sample_id")
gsva_for_flow <- gsva_for_flow[-1,]

flowdata_clean <- left_join(annotation, gsva_for_flow, by = "sample_id")
flowdata_clean <- left_join(flowdata_clean, flowdata, by = "sample_id")
left_join(annotation, flowdata, by = "sample_id") %>% filter(yfp_bulk == "BULK" & bailey_type == "Squamous") %>% dplyr::select(sample_id, flow_gMDSCs)
left_join(annotation, flowdata, by = "sample_id") %>% filter(yfp_bulk == "BULK" & bailey_type == "Pancreatic_Progenitor") %>% dplyr::select(sample_id, flow_gMDSCs)
# temp <- flowdata_clean %>% na.omit() %>% data.frame
flowdata_clean %>% filter(yfp_bulk == "YFP")
temp <- flowdata_clean %>% filter(yfp_bulk == "BULK")
flowdata_clean %>% dplyr::select(-matches("total_read_count")) %>%  export(., paste0("results/", Sys.Date(), "-annotation-gsva.csv"))
flowdata_clean %>% dplyr::select(-matches("total_read_count")) %>%  export(., paste0("results/", Sys.Date(), "-annotation-gsva.xlsx"))


plot(temp$gsva_CYT_index, temp$flow_CD8_T_cells)
cor.test(as.numeric(temp$gsva_CYT_index), as.numeric(temp$flow_CD8_T_cells))

plot(temp$flow_CD8_T_cells, temp$flow_CD4_T_cells)
cor.test(as.numeric(temp$flow_CD8_T_cells), as.numeric(temp$flow_CD4_T_cells))

plot(as.numeric(temp$gsva_DC), as.numeric(temp$flow_DCs))
cor.test(as.numeric(temp$gsva_aDC), as.numeric(temp$flow_DCs))
