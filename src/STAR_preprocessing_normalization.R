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

op <- par(mar = c(8, 5, 4, 2)+ 0.1)
options(op)

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

# filter out lowly expressed transcripts for PCA analysis
isexpr <- rowSums(tpm >= 1) >= 10
table(isexpr)
# FALSE  TRUE
# 10363 13150

tpm_filt <- tpm[isexpr,]

# running PCA on 1000 most variable transcripts
rv <- tpm_filt %>% select(-Geneid) %>% rowVars()
sel <- order(rv, decreasing = TRUE)[1:1000]

# using data.table syntax to log2 transform only rows in sel and columns not Geneid
tpm_pca <- log2(tpm_filt[sel, !"Geneid", with = FALSE] + 1)

# PCA
pc <- prcomp(t(tpm_pca))

# create annotation table
(annotation <- data.table(sample_id = names(tpm_filt)[grep("PD", names(tpm_filt))],
                          primary = stringr::str_split_fixed(names, "_", 2)[,1],
                          yfp_bulk = c("BULK", "YFP")))

# sample annotations and merging with read_sum_info
pca_dt <- data.table(sample_id = colnames(tpm_pca), pc$x)
pca_dt <- left_join(pca_dt, read_sum_info, by = "sample_id")
pca_dt <- left_join(pca_dt, annotation)

# colors_primary <- RColorBrewer::brewer.pal(10, "Paired")

# calculate variance explained by PC1 and PC2
varExp <- (pc$sdev^2/sum(pc$sdev^2) * 100)

# inital QC assesment using principal components 
(pcaplot <- ggplot(pca_dt, aes(x = PC1, y = PC2)) + 
  geom_text_repel(data = pca_dt, aes(label = sample_id)) + 
  geom_point(size = 4, aes(col = total_read_count, shape = primary)) + 
  ggtitle("Principal Components") +
  xlab(paste0("PC1\n(", round(varExp[1], 2), "% variance)")) +
  ylab(paste0("PC2\n(", round(varExp[2], 2), "% variance)")) +
  scale_color_gradient(low = "lightblue", high = "red") + 
  pca_theme() + theme(legend.position = "right"))
ggsave(pcaplot, file = "results/2016-09-29-PCA-noERCCnormalization-allsamples.pdf", width = 6, height = 6)

# heatmap of same genes as in PCA
anno_df <- data.frame(annotation, row.names = annotation[,1])
pheatmap::pheatmap(tpm_pca, scale = "row", fontsize_col = 6,
                   show_rownames = FALSE, 
                   annotation_names_col = FALSE,
                   clustering_distance_rows = "correlation", 
                   clustering_distance_cols = "correlation",
                   annotation_col = anno_df)

tpm
INDEX = c(immune_suppresion_index_genes)
# dplyr function chain to make a boxplot of log-average of genes in index 
# will need ConsensusCluster Annotations 
index_boxplot <- function(INDEX, TITLE){
  suppressWarnings(
    suppressMessages(
      tpm %>% 
        data.table() %>% 
        dplyr::filter(Geneid %in% INDEX) %>% 
        select(contains("BULK")) %>% 
        melt(variable.name = "sample_id", value.name = "TPM") %>% 
        group_by(sample_id) %>% 
        summarize(gm = geometric_mean(TPM)) %>% 
        # left_join(., gsva_anno_merge[,1:2], by = "sample_id") %>% 
        ggplot(., aes(x = Moffitt_tumor_type, y = log2(gm + 1), fill = Moffitt_tumor_type)) + geom_boxplot() + theme_bw() + 
        ggtitle(paste0(TITLE)) + xlab(" ") + ylab("log-average expression:\nlog2(TPM + 1)\n") +
        scale_fill_manual(values = colors_primary) + 
        ylim(c(0, 5)) +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "none"
        )
    )
  )
}

cyt_index_genes = c("Gzma", "Prf1")
immune_suppresion_index_genes = c("Ctla4", "Pdcd1", "Cd274", "Pdcd1lg2", "Ido1", "Ido2", "Lag3","Adora2a", "Havcr2", "Tigit", "Vtcn1", "Vsir")

moffitt_color <- RColorBrewer::brewer.pal(9, "Paired")[8:9]
cyt_index_boxplot <- index_boxplot(cyt_index_genes, "Cytolytic Index")
immune_suppresion_boxplot <- index_boxplot(immune_suppresion_index_genes, "Immune Suppresion Index")
p <- gridExtra::arrangeGrob(cyt_index_boxplot, immune_suppresion_boxplot, ncol = 2)
plot(p)
ggsave("Results/2016-09-28-Cyt-Suppresion-indexes.pdf", p, width = 6, height = 5)
