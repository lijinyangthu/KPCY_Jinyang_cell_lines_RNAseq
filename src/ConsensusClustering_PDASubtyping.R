#----------------------------------------------------------------------------------------------------------------------------------------------------------
# script 2 for analysis of KPCY cell lines
# Consensus clustering of KPCY cell lines versus Moffitt tumor signatures (Tumor and Stroma) and Bailey signatures
# references for both gene sets taken from publications:
# Moffitt manuscript - Nature Genetics, Virtual microdissection identifies distinct tumor- and stroma-specific subtypes of pancreatic ductal adenocarcinoma, 2015
# Bailey manuscript  - Nature, Genomic analyses identify molecular subtypes of pancreatic cancer, 2016
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------
library("tidyverse") # loads ggplot2, tibble, tidyr, readr, purr, dplyr
library("rio")
library("ConsensusClusterPlus")
library("gplots")
library("reshape2")
library("data.table")
library("magrittr")
library("biomaRt")
source("src/functions.R")

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
# set graphics parameters
# op <- par(mar = c(8, 5, 4, 2)+ 0.1)
# options(op)
# set seed for reproducibility
set.seed(123)
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------
library("biomaRt")
# need to map mouse MGI symbols to HUGO symbols
human <- useMart("ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "www.ensembl.org")
mouse <- useMart("ENSEMBL_MART_ENSEMBL", dataset = 'mmusculus_gene_ensembl', host = "www.ensembl.org")

sub_tpm$ensembl_transcript_id <- gsub(x = sub_tpm$Geneid, pattern = "*\\.[0-9]+", replacement = "\\1")
sub_tpm[, Geneid := NULL]

human_mouse_convesion <- getLDS(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", "entrezgene", "description"),
                                filters = "ensembl_transcript_id", values = sub_tpm$ensembl_transcript_id, mart = mouse,
                                attributesL = c('hgnc_symbol'), martL = human)

names(human_mouse_convesion) <- c("Geneid", "ensembl_gene_id", "ensembl_transcript_id", "entrezgene", "description","hgnc_symbol")
human_mouse_convesion %<>% data.table() %>%  dplyr::select(Geneid, hgnc_symbol, everything())
human_mouse_convesion

#------------------------------------------------------------------------------------------------------------------------------------------------------
# ConsensusClusterPlus definiation of KPCY with Moffitt (Nature Medicine 2015) Classical and Basal-like tumor types gene signatures
#------------------------------------------------------------------------------------------------------------------------------------------------------
basal <- read.delim("data/Moffitt_Basal-like_Siganture.txt", col.names = "hgnc_symbol")
classic <- read.delim("data/Moffitt_Classical_Siganture.txt", col.names = "hgnc_symbol")

basal$Subtype <- c(rep("Basal-like", 25))
classic$Subtype <- c(rep("Classical", 25))

moffitt <- rbind(basal, classic)
moffitt <- right_join(human_mouse_convesion, moffitt, by = "hgnc_symbol") %>%
  dplyr::select(ensembl_transcript_id, hgnc_symbol, Subtype) %>%
  unique() %>% data.table()
moffitt

#------------------------------------------------------------------------------------------------------------------------------------------------------
#  start here
#------------------------------------------------------------------------------------------------------------------------------------------------------
clustering_tpm_moffitt <- sub_tpm %>% filter(ensembl_transcript_id %in% moffitt$ensembl_transcript_id) %>% data.frame()
# convert to data.frame then matrix for CCP algorithm
row.names(clustering_tpm_moffitt) <- clustering_tpm_moffitt$ensembl_transcript_id
clustering_tpm_moffitt$ensembl_transcript_id <- NULL
dim(clustering_tpm_moffitt)
isexpr <- rowSums(clustering_tpm_moffitt >= 1) >=  4
clustering_tpm_moffitt <- clustering_tpm_moffitt[isexpr, ]
dim(clustering_tpm_moffitt)

log_moffitt <- log2(clustering_tpm_moffitt + 1)

pheatmap::pheatmap(as.matrix(log_moffitt), scale = "row")

# ConsensusClusterPlus function with K (clusters) set to 4
title <- "results/consensuscluster-Moffitt-sailfish/"
moffitt_res <- ConsensusClusterPlus(as.matrix(log_moffitt), maxK = 4, reps = 1000, pItem = 0.9,
                               pFeature = 1, title = title, clusterAlg = "hc", distance = "pearson",
                               seed = 42, plot = "pdf")

# two clusters with Moffitt PDA type signature
# extract item consensus
icl <- calcICL(moffitt_res, title = title, plot = "pdf")
icl[["clusterConsensus"]]
# only one group.  all Classical subtype

#------------------------------------------------------------------------------------------------------------------------------------------------------
# consensusclusteringplus definiation of KPCX with Moffitt normal and activated stroma types only in BULK samples
# these gene ids were manually pulled from Moffitt paper figure
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
clustering_tpm_stromal <- tpm %>% filter(Geneid %in% moffitt_stromal$Geneid) %>% dplyr::select(matches("Geneid|BULK")) %>% data.frame()

# convert to data.frame then matrix for CCP algorithm
row.names(clustering_tpm_stromal) <- clustering_tpm_stromal$Geneid
clustering_tpm_stromal$Geneid <- NULL

isexpr <- rowSums(clustering_tpm_stromal >= 1) >=  4
clustering_tpm_stromal <- clustering_tpm_stromal[isexpr,]
log_stromal <- log2(clustering_tpm_stromal + 1)
pheatmap::pheatmap(as.matrix(log_stromal), scale = "row")
###
##
# no clear groups in dataset - subcutaneous is potentitally different microenvironment
# running code, but not clustering samples as either groups
#
##
###
# ConsensusClusterPlus function
title <- "results/consensuscluster-Moffitt-stromal/"
moffitt_stromal_res <- ConsensusClusterPlus(as.matrix(log_stromal), maxK = 4, reps = 1000, pItem = 0.9,
                                    pFeature = 1, title = title, clusterAlg = "hc", distance = "pearson",
                                    seed = 42, plot = "pdf")

# two clusters with Moffitt PDA type signature
# extract item consensus
icl_strom <- calcICL(moffitt_stromal_res, title = title, plot = "pdf")
icl_strom[["clusterConsensus"]]
# group 2 but only 1 subtype - kind of a mess 
dt_strom <- data.table(icl_strom$itemConsensus)
dt_strom %>% dplyr::select(item) %>% summarise(n_dist = n_distinct(item))
# 17
dt_strom %>% filter(k == 2 & cluster == 1) %>% summarise(range = mean(itemConsensus))
dt_strom %>% filter(k == 2 & cluster == 2) %>% summarise(range = mean(itemConsensus))

# dt_strom %>% filter(k == 3 & item == "PD6419_C1_BULK")

# select sample in cluster 1 with consensus score greater than 0.79
(clust1_strom <- dt_strom %>%
  filter(k == 2 & itemConsensus > 0.5 & cluster == 1) %>%
  dplyr::select(cluster, item, itemConsensus))#  %>% summarise(n_dist = n_distinct(item))
#  13

# select sample in cluster 2 with consensus score greater than 0.4
(clust2_strom <- dt_strom %>%
  filter(k == 2 & itemConsensus > 0.5 & cluster == 2) %>%
  dplyr::select(cluster, item, itemConsensus)) # %>% summarise(n_dist = n_distinct(item))
# 4

table(clust2_strom$item %in% clust1_strom$item)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Bailey only Pancreatic Progenitor, Squamous GP3 signatures, and ADEX signatures from manuscript
#------------------------------------------------------------------------------------------------------------------------------------------------------
# GP1 = blue - pancreatic progenitor
blue <- read.table("data/blue_pancreatic_progenitor.txt")

# GP3 = Yellow - squamous
yellow <- read.table("data/yellow_squamous.txt")

Bailey_GP1_PP       <- tibble(hgnc_symbol = blue$V1,
                              Subtype = "GP1-Pancreatic Progenitor")
Bailey_GP3_Squamous <- tibble(hgnc_symbol = yellow$V1,
                              Subtype = "GP3-Squamous")
Bailey_GP9_ADEX     <- tibble(hgnc_symbol = c('REG3G','SYCN','REG1P','SERPINI2','CPB1',
                                   'RP11−331F4.4', 'AQP8','RBPJL','CUZD1','PLA2G1B',
                                   'CLPS','CEL','CELA3B','CELA3A','CTRB1','CTRC',
                                   'CTRB2','PNLIPRP1','CPA1','CPA2'),
                              Subtype = "GP9-ADEX")
Bailey_GP10_ADEX    <- tibble(hgnc_symbol = c('UNC79', 'GCK','GJD2','SYT4','KCNK16','NOL4',
                                    'SCGN','INS','CABP7','CHGB','BEX1','SVOP',
                                    'ABCC8','HMGCLL1','SLC30A8','SST','CELF3',
                                    'PCSK2','SCG3'),
                              Subtype = "GP10-ADEX")
bailey_merged <- rbind(Bailey_GP1_PP, Bailey_GP3_Squamous, Bailey_GP9_ADEX, Bailey_GP10_ADEX)

bailey_merged <- right_join(human_mouse_convesion, bailey_merged, by = "hgnc_symbol") %>%
  dplyr::select(Geneid, hgnc_symbol, Subtype) %>%
  unique() %>% na.omit() %>% data.table()
bailey_merged

# all samples
clustering_tpm_bail <- tpm %>% filter(Geneid %in% bailey_merged$Geneid) %>% data.frame()

# convert to data.frame then matrix for CCP algorithm
row.names(clustering_tpm_bail) <- clustering_tpm_bail$Geneid
clustering_tpm_bail$Geneid <- NULL

isexpr <- rowSums(clustering_tpm_bail >= 1) >=  4
clustering_tpm_bail <- clustering_tpm_bail[isexpr, ]

log_bail <- log2(clustering_tpm_bail + 1)

# ConsensusClusterPlus function with K set to 4
title <- "results/consensuscluster-Bailey/"
bail_res <- ConsensusClusterPlus(as.matrix(log_bail), maxK = 7, reps = 1000, pItem = 0.9,
                                    pFeature = 1, title = title, clusterAlg = "hc", distance = "pearson",
                                    seed = 42, plot = "pdf")

# 4 clusters with Bailey PDA type signature
# extract item consensus
icl <- calcICL(bail_res, title = title, plot = "pdf")
icl[["clusterConsensus"]]
# group 2 but only one cluster 

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Merged annnotation file for all samples versus selected tumor type definitions
#
# frozen on 12/02/2016
# uncomment to rerun
#------------------------------------------------------------------------------------------------------------------------------------------------------
names <- names(tpm)[grep("PD", names(tpm))]
(annotation <- data.table(sample_id = names(tpm)[grep("PD", names(tpm))],
                          primary = stringr::str_split_fixed(names, "_", 2)[,1],
                          yfp_bulk = regmatches(names(tpm), regexpr("(BULK|YFP)", names(tpm), perl = TRUE))))
annotation <- left_join(annotation, read_sum_info)
# merge with annotation dt from preprocessing
pda_typing <- data.table(sample_id = colnames(tpm[,grep("PD", colnames(tpm)), with = FALSE]),
                         moffitt_tumor_type = character(length(colnames(tpm[,grep("PD", colnames(tpm)), with = FALSE]))),
                         bailey_type = character(length(colnames(tpm[,grep("PD", colnames(tpm)), with = FALSE]))))
annotation <- left_join(annotation, pda_typing, by = "sample_id")

annotation$moffitt_tumor_type <- "Classical"
annotation$bailey_type <- "Pancreatic_progenitor"

# annotation$moffitt_stromal_type[annotation$yfp_bulk == "YFP"] <- "NA"

annotation %>% export(., file = paste0("results/", Sys.Date(), "-annotation.csv"))
annotation %>% export(., file = paste0("results/", Sys.Date(), "-annotation.xlsx"))
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
# data should be ready for further analysis
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------
