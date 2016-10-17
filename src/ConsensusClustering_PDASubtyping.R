# Consensus clustering of KPCY cell lines versus Moffitt tumor signatures (Tumor and Stroma)

library("tidyverse") # loads ggplot2, tibble, tidyr, readr, purr, dplyr
library("rio")
library("ConsensusClusterPlus")
library("biomaRt")
library("gplots")
library("reshape2")
library("data.table")
library("magrittr")
source("src/functions.R")

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
# set graphics parameters
op <- par(mar = c(8, 5, 4, 2)+ 0.1)
options(op)
# set seed for reproducibility
set.seed(123)
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------

# need to map mouse MGI symbols to HUGO symbols
human <- useMart("ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "www.ensembl.org")
mouse <- useMart("ENSEMBL_MART_ENSEMBL", dataset = 'mmusculus_gene_ensembl', host = "www.ensembl.org")

human_mouse_convesion <- getLDS(attributes = c("external_gene_name", "ensembl_gene_id", "entrezgene", "description"),
                                filters = "external_gene_name", values = tpm$Geneid, mart = mouse,
                                attributesL = c('hgnc_symbol'), martL = human)

names(human_mouse_convesion) <- c("Geneid", "ensembl_gene_id", "entrezgene", "description","hgnc_symbol")
human_mouse_convesion %<>% data.table() %>%  dplyr::select(Geneid, hgnc_symbol, everything())
human_mouse_convesion

#------------------------------------------------------------------------------------------------------------------------------------------------------
# consensusclusteringplus definiation of KPCX with Moffitt Classical and Basal-like tumor types
#------------------------------------------------------------------------------------------------------------------------------------------------------
basal <- read.delim("~/Projects/KPCY_RNAseq_EMT/GeneSignature/Moffitt_Basal-like_Siganture.txt", col.names = "hgnc_symbol")
classic <- read.delim("~/Projects/KPCY_RNAseq_EMT/GeneSignature/Moffitt_Classical_Siganture.txt", col.names = "hgnc_symbol")

basal$Subtype <- c(rep("Basal-like", 25))
classic$Subtype <- c(rep("Classical", 25))

moffitt <- rbind(basal, classic)
moffitt <- right_join(human_mouse_convesion, moffitt, by = "hgnc_symbol") %>%
  dplyr::select(Geneid, hgnc_symbol, Subtype) %>%
  unique() %>% data.table()
moffitt

# only primary tumors
# clustering_tpm_moffitt <- tpm %>% filter(Geneid %in% moffitt$Geneid) %>% dplyr::select(matches("YFP")) %>% data.frame()

# all samples
clustering_tpm_moffitt <- tpm %>% filter(Geneid %in% moffitt$Geneid) %>% data.frame()

# convert to data.frame then matrix for CCP algorithm
row.names(clustering_tpm_moffitt) <- clustering_tpm_moffitt$Geneid
clustering_tpm_moffitt$Geneid <- NULL

# ConsensusClusterPlus function
title <- "results/consensuscluster-Moffitt/"
moffitt_res <- ConsensusClusterPlus(as.matrix(clustering_tpm_moffitt), maxK = 4, reps = 1000, pItem = 0.8,
                               pFeature = 1, title = title, clusterAlg = "hc", distance = "pearson",
                               seed = 42, plot = "pdf")

# two clusters with Moffitt PDA type signature
# extract item consensus
icl <- calcICL(moffitt_res, title = title, plot = "pdf")
icl[["clusterConsensus"]]
# 2 group is really strong concensus
dt_icl <- data.table(icl$itemConsensus)
dt_icl %>% dplyr::select(item) %>% summarise(n_dist = n_distinct(item))

# select sample in cluster 1 with consensus score greater than 0.9
clust1 <- dt_icl %>%
  filter(k == 2 & itemConsensus > 0.9 & cluster == 1) %>%
  dplyr::select(cluster, item, itemConsensus)
clust1 # 7

# select sample in cluster 2 with consensus score greater than 0.9
clust2 <- dt_icl %>%
  filter(k == 2 & itemConsensus > 0.9 & cluster == 2) %>%
  dplyr::select(cluster, item, itemConsensus)
clust2 # 3

table(clust2$item %in% clust1$item); table(clust1$item %in% clust2$item)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# consensusclusteringplus definiation of KPCX with Moffitt normal and activated stroma types only in BULK samples
#------------------------------------------------------------------------------------------------------------------------------------------------------

activated_stroma <- data.table(hgnc_symbol = factor(c("SPARC",'COL1A2','COL3A1','POSTN',
                                               'COL5A2',' THBS2','FN1','COL10A1',
                                               'COL5A1','SFRP2','CDH11', 'CTHRC1',
                                               'FNDC1','SULF1','FAP','LUM','COL11A1','ITGA11',
                                               'MMP11', 'INHBA', "ZNF469",
                                               'VCAN', 'GREM1', 'COMP')),
                               Subtype = "Activated_Stroma")

normal_stroma <- data.table(hgnc_symbol = factor(c("FABP4", 'ACTG2','DES','IFG1','PTX3','OGN',
                                            'ADAMTS1','GPM6B','ANGTL7', 'MEOX2', 'SYNM',
                                            'MYH11','ID4','RSPO3','LMOD1','RBPMS2','PLP1',
                                            'VIT', 'LPHN3','SCRG1','CDH19','RERGL', "ABCA8")),
                            Subtype = "Normal_stroma")

moffitt_stromal <- rbind(normal_stroma, activated_stroma)

moffitt_stromal <- right_join(human_mouse_convesion, moffitt_stromal, by = "hgnc_symbol") %>%
  dplyr::select(Geneid, hgnc_symbol, Subtype) %>%
  unique() %>% na.omit() %>% data.table()
moffitt_stromal

# only primary tumors
clustering_tpm_stromal <- tpm %>% filter(Geneid %in% moffitt_stromal$Geneid) %>% dplyr::select(matches("BULK")) %>% data.frame()

# convert to data.frame then matrix for CCP algorithm
row.names(clustering_tpm_stromal) <- clustering_tpm_stromal$Geneid
clustering_tpm_stromal$Geneid <- NULL

# ConsensusClusterPlus function
title <- "results/consensuscluster-Moffitt-stromal/"
moffitt_stromal_res <- ConsensusClusterPlus(as.matrix(clustering_tpm_stromal), maxK = 4, reps = 1000, pItem = 0.8,
                                    pFeature = 1, title = title, clusterAlg = "hc", distance = "pearson",
                                    seed = 42, plot = "pdf")

# two clusters with Moffitt PDA type signature
# extract item consensus
icl_strom <- calcICL(moffitt_stromal_res, title = title, plot = "pdf")
icl_strom[["clusterConsensus"]]
# group 3?
dt_strom <- data.table(icl_strom$itemConsensus)
dt_strom %>% dplyr::select(item) %>% summarise(n_dist = n_distinct(item))

dt_strom %>% filter(k == 3 & cluster == 1) %>% summarise(range = mean(itemConsensus))
dt_strom %>% filter(k == 3 & cluster == 2) %>% summarise(range = mean(itemConsensus))
dt_strom %>% filter(k == 3 & cluster == 3) %>% summarise(range = mean(itemConsensus))

# select sample in cluster 1 with consensus score greater than 0.79
(clust1_strom <- dt_strom %>%
  filter(k == 3 & itemConsensus > 0.6 & cluster == 1) %>%
  dplyr::select(cluster, item, itemConsensus))#  %>% summarise(n_dist = n_distinct(item))
#  ZERO

# select sample in cluster 2 with consensus score greater than 0.4
(clust2_strom <- dt_strom %>%
  filter(k == 3 & itemConsensus > 0.6 & cluster == 2) %>%
  dplyr::select(cluster, item, itemConsensus)) # %>% summarise(n_dist = n_distinct(item))
# 2

# select sample in cluster 3 with consensus score greater than 0.32
(clust3_strom <-  dt_strom %>%
  filter(k == 3 & itemConsensus > 0.6 & cluster == 3) %>%
  dplyr::select(cluster, item, itemConsensus)) #   %>% summarise(n_dist = n_distinct(item))
# 2

table(clust2_strom$item %in% clust3_strom$item)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Merged annnotation file for all samples versus selected tumor type definitions
#------------------------------------------------------------------------------------------------------------------------------------------------------
pda_typing <- data.table(sample_id = colnames(tpm[,grep("PD", colnames(tpm)), with = FALSE]),
                         moffitt_tumor_type = character(length(colnames(tpm[,grep("PD", colnames(tpm)), with = FALSE]))),
                         moffitt_stromal_type = character(length(colnames(tpm[,grep("PD", colnames(tpm)), with = FALSE]))))

# merge with annotation dt from preprocessing
annotation <- left_join(annotation, pda_typing, by = "sample_id")
annotation$moffitt_tumor_type[annotation$sample_id %in% clust1$item] <- "Clust 1"
annotation$moffitt_tumor_type[annotation$sample_id %in% clust2$item] <- "Clust 2"
annotation$moffitt_stromal_type[annotation$sample_id %in% clust1_strom$item] <- "unknown stroma"
annotation$moffitt_stromal_type[annotation$sample_id %in% clust2_strom$item] <- "Stromal 2"
annotation$moffitt_stromal_type[annotation$sample_id %in% clust3_strom$item] <- "Stromal 3"
annotation$moffitt_stromal_type[annotation$yfp_bulk == "YFP"] <- "NA"
annotation

export(annotation, file = "results/annotation.csv")

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
# data should be ready for further analysis
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------
