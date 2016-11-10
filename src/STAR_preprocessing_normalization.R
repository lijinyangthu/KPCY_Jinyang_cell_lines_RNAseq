# STAR/featureCount quantification against mm10 reference

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
unloadNamespace("GenomicFeatures")

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
# set graphics parameters
op <- par(mar = c(10, 5, 2, 2)+ 0.1)
options(op)
# set seed for reproducibility
set.seed(123)
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------

# extract gene exon length info from gtf (function found in src/functions.r)
# gene_length_info <- gene_exon_length("data/genes.gtf")
# names(gene_length_info)[1] <- "Geneid"

# use data.table fread to load file.  will be loaded as a data.table.  ??data.table
counts <- fread("data/2016-11-10-star-counts.txt", skip = 1)
glimpse(counts)
summary(counts)
dim(counts)
# reformat column names from featureCount output
names(counts) <- gsub("*-Aligned.sortedByCoord.out.bam", "\\1", names(counts))
names(counts) <- gsub("data.bam.", "\\1", names(counts))

# fix sample PD5422_C5_YFP
names(counts)[names(counts) == "PD6322_C5_YFP"] <- "PD6422_C5_YFP"

# exclude Chr, Geneid, Start, End, Strand, Length
non_samples <- paste(names(counts)[1:6], collapse = "|")

read_sum_info <- counts %>% dplyr::select(-matches(non_samples)) %>%
  colSums() %>%
  data.table(total_read_count = .,
             sample_id = names(counts)[7:ncol(counts)]) %>%
  arrange(desc(total_read_count))

read_sum_info

read_sum_info %>% summarize(read_mean = mean(total_read_count))
# 33.5e6

# convert to counts per million and plot distribution removing outlier dots
counts %>%
  dplyr::select(-matches(non_samples)) %>%
  cpm() %>%
  data.table() %>%
  mutate_all(funs(log2(. + 1))) %>%
  boxplot(., las = 2, outline = F,
          main = "raw CPM for initial run",
          ylab = "Log2(CPM + 1)")
dev.copy(pdf, file = paste0("results/", Sys.Date(), "-KPCY-cell-lines.pdf"))
dev.off()

# drop samples with less than 1e6 total reads and keep Geneid
(poor_samples <- read_sum_info %>%
  filter(total_read_count < 1e6) %>%
  dplyr::select(sample_id) %>%
  unlist() %>%
  paste(., collapse = "|"))
non_geneid <- paste(names(counts)[2:5], collapse = "|")

sub_counts <- counts %>%
  dplyr::select(-matches(non_geneid))


# convert counts to TPM (function found in src/functions.R)
tpm <- STAR_to_TPM(sub_counts)
dim(tpm)
# 23513 25
# save to data/
export(tpm, file = paste0("data/", Sys.Date(), "-tpm.csv"))


# create annotation table
names <- names(tpm)[grep("PD", names(tpm))]
(annotation <- data.table(sample_id = names(tpm)[grep("PD", names(tpm))],
                          primary = stringr::str_split_fixed(names, "_", 2)[,1],
                          yfp_bulk = regmatches(names(tpm), regexpr("(BULK|YFP)", names(tpm), perl = TRUE))))

annotation <- left_join(annotation, read_sum_info)

# PCA analysis + plot against total read sum number (function found in src/functions.R)
(pca_su <- PCA_analysis(expr_obj = tpm, annotation = annotation, colors = NULL, top_var = 1000) +
  geom_point(size = 4, aes(col = yfp_bulk)) + theme(legend.position = "bottom"))

# ggsave(pca_su, file = )

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
# data should be ready for further analysis
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------
