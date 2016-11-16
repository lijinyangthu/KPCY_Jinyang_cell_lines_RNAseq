# common functions

# function to combine multiple Sailfish output files into merged table and save to RDS
combine_sailfish <- function(files, target, target_file) {

  for(i in 1:length(files)){
    tab <- read.table(file = paste0(directory, files[i]), header = FALSE)
    names(tab) <- c("Transcript", "Length", "TPM", "RPKM", "KPKM", "EstimatedNumKmers", "EstimatedNumReads")
    assign(files[i], tab)
  }

  sample_list <- mget(files)
  ListNames <- Map(function(x,i)
    setNames(x, ifelse(names(x) %in% "Transcript", names(x), sprintf("%s.%d", names(x), i))),
    sample_list, seq_along(sample_list))

  # merge and subset for counts/TPM
  full_results <- Reduce(function(...) merge(..., by = "Transcript", all = TRUE), ListNames)
  cnt_tmp <- full_results[, grep(target, names(full_results))]
  row.names(cnt_tmp) <- full_results$Transcript
  colnames(cnt_tmp) <- gsub("\\-sailfish.txt", "", files)
  # save and return
  saveRDS(cnt_tmp, file = target_file)
  return(cnt_tmp)
}

#' function to convert STAR/featureCounts output to transcripts per million (TPM)
#'
#' @param counts data.table or data.frame of unprocessed STAR/featureCounts output.  columns are samples and rows are gene features. expected to contain featureCounts columns for "Geneid", "Strand", "Chr", "Start", "End", and"Lenght"
#' @param gene_exon_length data.table or data.frame of Geneid and effective gene length for TPM calculation. Must have same row length as counts. Optional
#' @usage
#' tpm <- STAR_to_TPM(counts = star_counts, gene_exon_length = gene_exon_length)
#' tpm <- STAR_to_TPM(counts = star_counts)
STAR_to_TPM <- function(counts, gene_exon_length = NULL) {
  packs <- c("dplyr", "magrittr","data.table", "stringr")
  missing_packs <- setdiff(packs, row.names(installed.packages()))
  if(length(missing_packs >= 1)){
    install.packages(missing_packs)
    lapply(packs, library, character = TRUE)
  } else {
    lapply(packs, library, character = TRUE)
  }

  counts_temp <- counts
  # if counts is data.table = convert to df and set row.names equal to Geneid
  if( class(counts)[1] == "data.table" )
    counts_temp <- counts %>%
    dplyr::select(-matches("Geneid|Strand|Chr|Start|End|Length")) %>%
    data.frame()
  row.names(counts_temp) <- counts$Geneid

  if(is.null(gene_exon_length))
    eff_len <- counts$Length
  else if (!is.null(gene_exon_length))
    eff_len <- gene_exon_length
  else
    message("need length info")

  # scale counts to transcript per million transcripts
  counts_temp[counts_temp == 0] <- NA
  rate <- log(counts_temp) - log(eff_len)
  denom <- log(colSums(exp(rate), na.rm = TRUE))
  out <- exp(t(t(rate) - denom) + log(1e6))
  out[is.na(out)] <- 0
  # convert matrix back to data.table
  out %<>%
    data.frame %>%
    rownames_to_column("Geneid") %>%
    data.table()

  # fix column names for downstream merges and joins
  names(out) <- str_replace_all(names(out), "\\.", "-")
  return(out)
}

#' function to run PCA analysis on a transcripts per million table
#' @param expr_obj data.table or data.frame from STAR_to_TPM output. should tpm table
#' @param annotation data.table or data.frame containing annotation information for samples in expr_obj. must contain column of sample IDs called "sample_id".
#' @param colors Optional vector of colors for PCA plot.  Defaults to RColorBrewer's Paired palette.
#' @param top_var Optional numeric value of the desired top X variable genes (e.g. top_var = 1000 will return top 1000 variable genes). defaults to 1000
#' @usage
#' pca_plot <- PCA_analysis(expr_obj = single_cell_tpm, annotation = single_cell_annotation)
#' pca_plot + geom_point(size = 4, aes(color = group))
#' ggsave(pca_plot, file = "PCA_plot.pdf", width = 6, height = 6)
PCA_analysis <- function(expr_obj, annotation, colors, top_var = 1000) {
  packs <- c("dplyr", "data.table", "genefilter", "RColorBrewer")
  missing_packs <- setdiff(packs, row.names(installed.packages()))
  if(length(missing_packs >= 1)){
    install.packages(missing_packs)
    lapply(packs, library, character = TRUE)
  } else {
    lapply(packs, library, character = TRUE)
  }

  # filter out lowly expressed transcripts - more than 1 tpm in at least 30% of samples
  isexpr <- rowSums(expr_obj >= 1) >= round(ncol(expr_obj) * 0.3)
  expr_filt <- expr_obj[isexpr,]
  message(paste0("Keeping ", table(isexpr)[2], " genes after filtering"))

  # running PCA on 1000 most variable transcripts
  rv <- expr_filt %>%
    dplyr::select(-Geneid) %>%
    rowVars()

  sel <- order(rv, decreasing = TRUE)[1:top_var]

  # using data.table syntax to log2 transform only rows in sel and columns not Geneid
  exp_pca <- log2(expr_filt[sel, !"Geneid", with = FALSE] + 1)

  # PCA
  pc <- prcomp(t(exp_pca))

  # sample annotations and merging with read_sum_info
  pca_dt <- data.table(sample_id = colnames(exp_pca),
                       pc$x)

  if (!is.null(annotation)) {
    pca_dt <- left_join(pca_dt, annotation, by = "sample_id")
  }

  if (!is.null(colors)){
    colors_primary <- colors
  } else {
    colors_primary <- brewer.pal(10, "Paired")
  }

  # calculate variance explained by PC1 and PC2
  varExp <- (pc$sdev^2/sum(pc$sdev^2) * 100)

  # generic PCA plot - coloring based on sample_id
  pca_p <- ggplot(pca_dt, aes(x = PC1, y = PC2)) +
    geom_text_repel(data = pca_dt, aes(label = sample_id)) +
    geom_point(size = 4) +
    ggtitle("Principal Components") +
    xlab(paste0("PC1\n(", round(varExp[1], 2), "% variance)")) +
    ylab(paste0("PC2\n(", round(varExp[2], 2), "% variance)")) +
    pca_theme() +
    theme(legend.position = "right")

  return(pca_p)
}

# PCA theme options used in PCA_analysis() function
pca_theme <- function(base_size = 12, base_family = "Helvetica") {
  theme(
    plot.title = element_text(size = 15, face = 'bold'),
    legend.key.size = unit(0.3, 'cm'),
    legend.position = "bottom",
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA))
}

# function to run statistical tests of GSVA results against subtype level consensus clustering
#' @param gsva_result matrix from GSVA output.
#' @param design_matrix matrix of column names of GSVA ouput with appropriate comparisons. See limma/edgeR design matrix.
#' @param contrast_matrix matrix of comparison names between design_matrix groups (e.g. Basal-Classical)
#' @usage
#' p_m <- paste(annotation[match(names(gsva_tpm), annotation$sample_id), ]$moffitt_tumor_type, sep = '')
#' p_m <- factor(p_m)
#' pm_design <- model.matrix(~0+p_m)
#' colnames(pm_design) <- c("Basal_like", "Classical")
#' pm_design
#' ContrastMatrix <- limma::makeContrasts(Basal_like-Classical, levels = pm_design)
#'
#' immune_res <- gene_set_statistic_test(gsva_result = immune_gsva, design_matrix = pm_design, contrast_matrix = ContrastMatrix)
#' immune_res %<>% arrange(desc(logFC)) %>% data.table()
#' immune_res
#' rio::export(immune_res, file = paste0("Results/", Sys.Date(), "-GSVA-immune-bulk.csv"))
#'
gene_set_statistic_test <- function(gsva_result, design_matrix, contrast_matrix){
  packs <- c("limma", "data.table", "edgeR", "tidyverse", "magrittr")
  missing_packs <- setdiff(packs, row.names(installed.packages()))
  if(length(missing_packs >= 1)){
    install.packages(missing_packs)
    lapply(packs, library, character = TRUE)
  } else {
    lapply(packs, library, character = TRUE)
  }

  fit_tmp <- lmFit(gsva_result, design = design_matrix)
  fit_tmp <- contrasts.fit(fit_tmp, contrast_matrix)
  fit_tmp <- eBayes(fit_tmp)
  print(summary(decideTests(fit_tmp)))

  res <- topTable(fit_tmp, adjust = "BH", sort = "l", n = Inf, p = 0.05)

  res <- rownames_to_column(res, "Gene_set")
  return(res)
}

#' function to calculate total exonic length of a gene. input for STAR_to_TPM
#' @param gtf GTF annotation file containing Geneid, exon length information
#' @examples
#'  gtf <- "data/mm10.gtf"
#'  mm10_gene_lenth <- gene_exon_length(gtf)
#'
gene_exon_length <- function(gtf) {
  packs <- c("GenomicFeatures")
  missing_packs <- setdiff(packs, row.names(installed.packages()))
  if(length(missing_packs >= 1)){
    install.packages(missing_packs)
    lapply(packs, library, character = TRUE)
  } else {
    lapply(packs, library, character = TRUE)
  }

  txdb <- makeTxDbFromGFF(gtf, format = "gtf")
  exons_per_gene <- exonsBy(txdb, by = "gene")
  gene_size <- lapply(exons_per_gene, function(x) sum(width(reduce(x))))
  gene_df <- data.frame(Gene_id = names(gene_size), Length = unlist(gene_size), row.names = NULL)
  return(gene_df)
}

# function to calculate geometric mean for normalized expression values)
geometric_mean <- function(x, na.rm = TRUE) {
  if (is.null(nrow(x))){
    2^(mean(log2(x + 1), na.rm = TRUE))
  } else {
    2^(apply(log2(x + 1), 2, mean, na.rm = na.rm))
  }
}