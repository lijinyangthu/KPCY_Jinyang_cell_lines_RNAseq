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
unloadNamespace("biomaRt")

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
# set graphics parameters 
op <- par(mar = c(8, 5, 4, 2)+ 0.1)
options(op)
# set seed for reproducibility
set.seed(123)
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------

# use data.table fread to load file.  will be loaded as a data.table.  ??data.table
counts <- fread("results/2016-09-29-star-counts.txt")
glimpse(counts)
summary(counts)

# reformat column names from featureCount output
names(counts) <- gsub("*-Aligned.sortedByCoord.out.bam", "\\1", names(counts))
names(counts) <- gsub("data.bam.", "\\1", names(counts))

# fix sample PD5422_C5_YFP 
names(counts)[names(counts) == "PD6322_C5_YFP"] <- "PD6422_C5_YFP"

# exclude Chr, Geneid, Start, End, Strand, Length
non_samples <- paste(names(counts)[1:6], collapse = "|")

read_sum_info <- counts %>% select(-matches(non_samples)) %>% 
  colSums() %>% 
  data.table(total_read_count = ., 
             sample_id = names(counts)[7:ncol(counts)]) %>% 
  arrange(desc(total_read_count))
  
read_sum_info

read_sum_info %>% summarize(read_mean = mean(total_read_count))
# 35.65e6

# convert to counts per million and plot distribution removing outlier dots
counts %>% 
  select(-matches(non_samples)) %>% 
  cpm() %>% 
  data.table() %>% 
  mutate_all(funs(log2(. + 1))) %>% 
  boxplot(., las = 2, outline = F, 
          main = "raw CPM for initial run",
          ylab = "Log2(CPM + 1)")
dev.copy(pdf, file = paste0("results/", Sys.Date(), "-KPCY-cell-lines.pdf"))
dev.off()

# drop samples with less than 1e6 total reads, no ERCCs, and keep Geneid
(poor_samples <- read_sum_info %>% 
  filter(read_sum < 1e6) %>% 
  select(sample_id) %>% 
  unlist() %>% 
  paste(., collapse = "|"))
non_geneid <- paste(names(counts)[2:5], collapse = "|")

sub_counts <- counts %>% 
  select(-matches(non_geneid))
sub_counts

# convert counts to TPM 
tpm <- STAR_to_TPM(sub_counts)
dim(tpm)
# 23513 11
# save to data/
export(tpm, file = paste0("data/", Sys.Date(), "-tpm.csv"))

# create annotation table
(annotation <- data.table(sample_id = names(tpm)[grep("PD", names(tpm))],
                          primary = stringr::str_split_fixed(names, "_", 2)[,1],
                          yfp_bulk = regmatches(names(tpm), regexpr("(BULK|YFP)", names(tpm), perl = TRUE))))

annotation <- left_join(annotation, read_sum_info)

# PCA analysis + plot against total read sum number
(pca_su <- PCA_analysis(expr_obj = tpm, annotation = annotation, colors = NULL, top_var = 1000) +
  geom_point(size = 4, aes(col = yfp_bulk)) + theme(legend.position = "bottom"))

# ggsave(pca_su, file = )

# heatmap of same genes as in PCA
anno_df <- data.frame(annotation, row.names = annotation[,1])

pheatmap::pheatmap(tpm_pca, scale = "row", fontsize_col = 6,
                   show_rownames = FALSE, 
                   annotation_names_col = FALSE,
                   clustering_distance_rows = "correlation", 
                   clustering_distance_cols = "correlation",
                   annotation_col = anno_df[, c("yfp_bulk", "moffitt_tumor_type", "moffitt_stromal_type")],
                   file = "results/PCA-heatmap.pdf")

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
# data should be ready for further analysis
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------
