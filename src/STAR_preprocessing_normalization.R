#----------------------------------------------------------------------------------------------------------------------------------------------------------
# script 1 for analysis of KPCY cell lines
# processing of STAR/featureCount output from Jinyang's KPCY cell lines sorted and bulk tumors
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------

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
gene_length_info <- gene_exon_length("data/genes.gtf")
names(gene_length_info)[1] <- "Geneid"

# use data.table fread to load file.  will be loaded as a data.table.
counts <- fread("data/2016-11-14-star-counts.txt", skip = 1)
glimpse(counts)
summary(counts)
dim(counts)
# reformat column names from featureCount output
names(counts) <- gsub("*-Aligned.sortedByCoord.out.bam", "\\1", names(counts))
names(counts) <- gsub("data.bam.", "\\1", names(counts))

# fix sample PD5422_C5_YFP id - incorrectly labelled as PD6322_C5_YFP
names(counts)[names(counts) == "PD6322_C5_YFP"] <- "PD6422_C5_YFP"

# exclude Chr, Geneid, Start, End, Strand, Length
non_samples <- paste(names(counts)[1:6], collapse = "|")

# calculate total read count metrics for each library
read_sum_info <- counts %>% dplyr::select(-matches(non_samples)) %>%
  colSums() %>%
  data.table(total_read_count = .,
             sample_id = names(counts)[7:ncol(counts)]) %>%
  arrange(desc(total_read_count))

read_sum_info
read_sum_info %>% summarize(read_mean = mean(total_read_count))
# 41e6 reads per sample

# convert to counts per million (edgeR cpm() function), log transform and plot distribution removing outlier dots
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
# all samples > 1e6 reads per
non_geneid <- paste(names(counts)[2:5], collapse = "|")

# removing PD6883_C1_BULK because poor library (probably over amplification or to much loaded on nextseq, e.g 100e6 total reads)
# removing PD6883_C4_BULK as it is a significant outlier for consensus clustering
sub_counts <- counts %>%
  dplyr::select(-matches(non_geneid), -matches("PD6883_C1_BULK"), -matches("PD6883_C4_BULK"))
names(sub_counts)
dim(sub_counts)

# convert featureCount gene-level counts to transcripts per million transcript TPM (function found in src/functions.R)
tpm <- STAR_to_TPM(sub_counts)
dim(tpm)
# 23513 39
# save to data/
export(tpm, file = paste0("data/", Sys.Date(), "-tpm.csv"))

# create annotation table.  Will fill this table with consensus clustering results, GSVA enrichment scores, and flow-based data per sample
names <- names(tpm)[grep("PD", names(tpm))]
(annotation <- data.table(sample_id = names(tpm)[grep("PD", names(tpm))],
                          primary = stringr::str_split_fixed(names, "_", 2)[,1],
                          yfp_bulk = regmatches(names(tpm), regexpr("(BULK|YFP)", names(tpm), perl = TRUE))))
annotation <- left_join(annotation, read_sum_info)

my_palette <- RColorBrewer::brewer.pal(8, "Set2")[3:4]
# PCA analysis + plot against total YFP or BULK tomorrow (function found in src/functions.R)
(pca_su <- PCA_analysis(expr_obj = tpm, annotation = annotation, colors = NULL, top_var = 2000) +
  geom_point(size = 4, aes(col = yfp_bulk)) + theme(legend.position = "bottom") +
    scale_color_manual(values = my_palette))
ggsave(pca_su, file = "results/2016-11-13-PCA.pdf", width = 6, height = 6)

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#
# data should be ready for further analysis
#
#----------------------------------------------------------------------------------------------------------------------------------------------------------


# # heatmap
# # filter out lowly expressed transcripts - more than 1 tpm in at least 30% of samples
# isexpr <- rowSums(tpm >= 1) >= (ncol(tpm) * 0.3)
# tpm_filt <- tpm[isexpr,]
#
# # running PCA on 1000 most variable transcripts
# rv <- tpm_filt %>%
#   dplyr::select(-Geneid) %>%
#   rowVars()
#
# sel <- order(rv, decreasing = TRUE)[1:top_var]
#
# # using data.table syntax to log2 transform only rows in sel and columns not Geneid
# exp_pca <- log2(tpm_filt[sel, !"Geneid", with = FALSE] + 1)
#
# names <- tpm_filt[sel,]
# names$Geneid
# pheatmap::pheatmap(exp_pca, scale = "row",
#                    clustering_method = "average",
#                    clustering_distance_rows = "correlation",
#                    color = colorRampPalette(c("navy", "white", "firebrick3"))(2345))
