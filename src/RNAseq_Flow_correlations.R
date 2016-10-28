# Correlation of Cell type gene set enrichment versus flow immunophenotyping
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

flowdata <- fread("data/20160430_Leukocyte_immune 2.csv")
corrected_names <- paste0("PD", names(flowdata[,2:ncol(flowdata), with = FALSE]))
corrected_names <- str_split_fixed(corrected_names, "c", 2)
corrected_names <- paste0(corrected_names[,1], "_C", corrected_names[,2], "_BULK")

# transpose, convert ID names, and merge with annotation file
flowdata <- data.table::transpose(flowdata)
flowdata <- flowdata[-1,]
setnames(flowdata, c("CD4_T_cells", "CD8_T_cels", "gMDSCs", "DCs", "CD4_CD8_null_T_cells", "F480_macrophages", "B_cells"))
flowdata <- flowdata[, sample_id := corrected_names]


left_join(annotation, flowdata, by = "sample_id") %>% filter(yfp_bulk == "BULK")
