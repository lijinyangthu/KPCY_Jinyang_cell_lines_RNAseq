# --------------------------------------------------------------------------------------------------------------------------------
# heatmap of T cell related genes from recent Baily, Scientific reports
# --------------------------------------------------------------------------------------------------------------------------------
t_cell_goi <- grep("Lck|Cd8a|Cd8b1|Gzma|Prf1|Ifng|H2-|Gdpd5|Ccr2|Ccr5|Mefv|Fli1|Tbc1d5|Ddx17|Akt3|Ewsr1|Tbcd|Nfatc4", tpm$Geneid)

# filter TPM values for Geneid and samples labeled as BULK
goi_tpm <- tpm[t_cell_goi] %>% dplyr::select(matches("Geneid|BULK")) %>% data.frame()
row.names(goi_tpm) <- goi_tpm$Geneid
goi_tpm$Geneid <- NULL

isexpr <- rowSums(goi_tpm) >= 0.5
goi_tpm <- log2(goi_tpm[isexpr,] + 1)
pheatmap::pheatmap(goi_tpm, annotation_col = annotation_df[,c(4,5)],
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(2345), scale = "row",
                   main = "Antigen presentation, cytolytic and effector activity",
                   file = "results/2016-11-30-antigen-cytolytic.pdf",
                   width = 10, height = 8.5)

# --------------------------------------------------------------------------------------------------------------------------------
# heatmap of dendritic cell related genes 
# --------------------------------------------------------------------------------------------------------------------------------
dc_goi <- c(immunome[immunome$cell_type == "DC", ]$Geneid, immunome[immunome$cell_type == "aDC", ]$Geneid, immunome[immunome$cell_type == "iDC", ]$Geneid)

# filter TPM values for Geneid and samples labeled as BULK
dc_tpm <- tpm %>% filter(Geneid %in% dc_goi) %>% dplyr::select(matches("Geneid|BULK")) %>% data.frame()
row.names(dc_tpm) <- dc_tpm$Geneid
dc_tpm$Geneid <- NULL
isexpr <- rowSums(dc_tpm) >= 0.5
dc_tpm <- log2(dc_tpm[isexpr,] + 1)

pheatmap::pheatmap(dc_tpm, annotation_col = annotation_df[,c(4,5)],
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(2345), scale = "row",
                   main = "Immature and activated dendritic cell expression",
                   file = "results/2016-11-30-dendritic.pdf",
                   width = 10, height = 8)

# --------------------------------------------------------------------------------------------------------------------------------
# heatmap of MDSC/Treg related genes 
# --------------------------------------------------------------------------------------------------------------------------------
treg_goi <- c(immunome[immunome$cell_type == "MDSC", ]$Geneid, immunome[immunome$cell_type == "TReg", ]$Geneid)
treg_goi

# filter TPM values for Geneid and samples labeled as BULK
trg_tpm <- tpm %>% filter(Geneid %in% treg_goi) %>% dplyr::select(matches("Geneid|BULK")) %>% data.frame()
row.names(trg_tpm) <- trg_tpm$Geneid
trg_tpm$Geneid <- NULL
isexpr <- rowSums(trg_tpm) >= 0.5
trg_tpm <- log2(trg_tpm[isexpr,]+ 1)

pheatmap::pheatmap(trg_tpm, annotation_col = annotation_df[,c(4,5)],
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(2345), scale = "row",
                   main = "Treg and MDSC expression",
                   file = "results/2016-11-30-treg.pdf",
                   width = 10, height = 4)

# --------------------------------------------------------------------------------------------------------------------------------
# heatmap of Macrophage related genes 
# --------------------------------------------------------------------------------------------------------------------------------
m_goi <- immunome[immunome$cell_type == "Macrophages", ]$Geneid
m_goi

# filter TPM values for Geneid and samples labeled as BULK
m_tpm <- tpm %>% filter(Geneid %in% m_goi) %>% dplyr::select(matches("Geneid|BULK")) %>% data.frame()
row.names(m_tpm) <- m_tpm$Geneid
m_tpm$Geneid <- NULL
isexpr <- rowSums(m_tpm) >= 0.5
m_tpm <- log2(m_tpm[isexpr,]+ 1)

pheatmap::pheatmap(m_tpm, annotation_col = annotation_df[,c(4,5)],
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(2345), scale = "row",
                   main = "Macrophage expression",
                   file = "results/2016-11-30-mphage.pdf",
                   width = 10, height = 8)
