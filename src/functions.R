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



STAR_to_TPM <- function(counts, eff_length = NULL) {
  counts_temp <- counts
  # if counts is data.table = convert to df and set row.names equal to Geneid
  if( class(counts)[1] == "data.table" )
    counts_temp <- counts %>% 
      select(-matches("Geneid|Strand|Chr|Start|End|Length")) %>% 
      data.frame()
  row.names(counts_temp) <- counts$Geneid
  
  if( is.null(eff_length))
    eff_len <- counts$Length
  else if (!is.null(eff_length))
    eff_len <- eff_length
  else
    message("need length info")
  
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
  names(out) <- stringr::str_replace_all(names(out), "\\.", "-")
  return(out)
}



# PCA theme options 
pca_theme <- function(base_size = 12, base_family = "Helvetica") {
  theme(
    plot.title = element_text(size = 15, face = 'bold'),
    legend.key.size = unit(0.3, 'cm'),
    legend.position = "bottom", 
    legend.key = element_rect(fill = 'NA'),
    #legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(colour = "Black"), 
    axis.text.x = element_text(colour = "Black"),
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(), 
    panel.grid.minor.y = element_blank(), 
    panel.background = element_rect(color = 'black',fill = NA))
}

# function to calculate total exonic length of a gene
# example:
# gtf <- "data/mm10.gtf"
# mm10_gene_lenth <- gene_exon_length("mm10.gtf")
gene_exon_length <- function(gtf) {
  suppressPackageStartupMessages(library("GenomicFeatures"))
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