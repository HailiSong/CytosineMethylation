rm(list=ls())
library(tidyverse)
library(gplots)
library(amap)
library(RColorBrewer)
library(fastcluster)
library(ggplot2)
library(reshape2)
library(viridis)
library(ggpubr)
library(qvalue)


root_wd <- "/Users/haili/Documents/WorkUT/mc_florence/bismarkCX/1to1orthologous/correction/avgMsites"

setwd(paste0(root_wd, "/figures"))

# Read in Expression data
exp_D1 <- read_csv("/Users/haili/Documents/WorkUT/result/mc_florence/D1_Cont_Sig.csv")
names(exp_D1) <- c("Gene", names(exp_D1)[2:length(names(exp_D1))])

exp_D4 <- read_csv("/Users/haili/Documents/WorkUT/result/mc_florence/D4_Cont_Sig.csv")
names(exp_D4) <- c("Gene", names(exp_D4)[2:length(names(exp_D4))])

# stage = "P1"
# feature = "promoter"
# ctype = "CG"
# meth = "median"
# exp = exp_D1


addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

generateBoxplotOverlayViolin <- function(stage, feature, ctype, exp){
  # define directories of input files
  fil_dir <- paste0(root_wd, "/data/FIL")
  hal_dir <- paste0(root_wd, "/data/HAL")
  
  # Read in mc files
  mc_fil <- read_tsv(paste0(fil_dir, "/", stage, "_orthologous_joined_",ctype, "_", feature, ".txt"))
  mc_hal <- read_tsv(paste0(hal_dir, "/", stage, "_orthologous_joined_",ctype, "_", feature, ".txt"))
  
  # m_inner_joined <- mc_hal %>% inner_join(mc_fil, by=c("joined_ID"))
  # m_inner_joined$fold_change <- m_inner_joined$FIL_Mean / m_inner_joined$HAL_Mean
  # m_inner_joined[is.na(m_inner_joined)] = 0
  # fpval <- m_inner_joined %>%
  #   rowwise() %>%
  #   mutate(pval = t.test(c(HAL_r1, HAL_r2, HAL_r3),
  #                        c(FIL_r1, FIL_r2, FIL_r3))$p.value) %>%
  #   ungroup()
  # 
  # # Save p value
  # write.table(fpval,
  #             paste0(stage, "_", feature,"_", ctype, "_inner_noexp_avg_methy_pval.tsv"),
  #             sep = "\t", quote = F)
  
  m_full_join <- mc_hal %>% full_join(mc_fil, by=c("joined_ID"))
  
  # set NA to zero
  NAs <- m_full_join %>% filter(is.na(Gene.x) | is.na(Gene.y))
  
  # Create empty columns
  NAs$pval <- "NA"
  NAs$qvalue <- "NA"
  NAs$fdr <- "NA"
  NAs$fold_change <- "NAs"
  
  
  m_full_join <- m_full_join %>% filter(!is.na(Gene.x) & !is.na(Gene.y))
  m_full_join[is.na(m_full_join)] <- 0
  
  ffpval <- m_full_join %>%
    rowwise() %>%
    mutate(pval = t.test(c(HAL_r1, HAL_r2, HAL_r3),
                         c(FIL_r1, FIL_r2, FIL_r3))$p.value) %>%
    ungroup()
  
  qobj <- qvalue(p = ffpval$pval)
  
  # qvalues
  ffpval$qvalue <- qobj$qvalues
  
  #fdr
  ffpval$fdr <- qobj$lfdr
  
  # fold change
  ffpval$fold_change <- ffpval$FIL_Mean/ffpval$HAL_Mean
  
  ffpval <- rbind(ffpval, NAs)
  
  write.table(ffpval,
              paste0(stage, "_", feature,"_", ctype, "_full_noexp_avg_methy_pval.tsv"),
              sep = "\t", quote = F)
  
  # ----------------------------------
  # 
  # # Inner join
  # inner_joined <- mc_hal %>% inner_join(mc_fil, by=c("joined_ID"))
  # 
  # # selected columns
  # inner_selected <- inner_joined %>% dplyr::select(joined_ID, HAL_r1, HAL_r2, HAL_r3,HAL_Mean,
  #                                                  FIL_r1, FIL_r2, FIL_r3, FIL_Mean)
  # 
  # # Change to data frame
  # inner_selected <- as.data.frame(inner_selected)
  # rownames(inner_selected) <- inner_selected$joined_ID
  # inner_selected$joined_ID <- NULL
  # 
  # # Rename
  # names(inner_selected) <- c("HAL_r1", "HAL_r2", "HAL_r3", "HAL_Mean",
  #                            "FIL_r1", "FIL_r2", "FIL_r3", "FIL_Mean")
  # 
  # # Convert all columns to numeric
  # inner_selected[, c(1:8)] <- sapply(inner_selected[, c(1:8)], as.numeric)
  # 
  # # replace all NA with 0
  # inner_selected[is.na(inner_selected)] = 0
  # 
  # # remove rows having all zeros
  # inner_selected <- inner_selected[rowSums(inner_selected[])>0,]
  # 
  # ### Link with expression
  # inner_selected <- as_tibble(tibble::rownames_to_column(inner_selected, "joined_ID"))
  # inner_m <- inner_selected %>% separate(joined_ID, c("HAL_ID", "FIL_ID"), sep = "_")
  # 
  # exp_m <- exp %>% dplyr::inner_join(inner_m, by = c("Gene" = "HAL_ID") )
  # 
  # # Add regulation data
  # exp_m <- exp_m %>% mutate(regu =
  #                             case_when(log2FoldChange <=0 ~ "down", 
  #                                       log2FoldChange >0 ~ "up"))
  # 
  # exp_m <- exp_m %>% dplyr::select(1, 9:17)
  # exp_m$fold_change <- exp_m$FIL_Mean/exp_m$HAL_Mean
  # 
  # pval <- exp_m %>% 
  #   rowwise() %>%
  #   mutate(pval = t.test(c(HAL_r1, HAL_r2, HAL_r3),
  #                        c(FIL_r1, FIL_r2, FIL_r3))$p.value) %>%
  #   ungroup()
  # 
  # # Save p value
  # write.table(pval, 
  #             paste0(stage, "_", feature,"_", ctype, "_inner_exp_avg_methy_pval.tsv"),
  #             sep = "\t", quote = F)
  # 
  # # calculate average
  # avg_exp_m <- exp_m %>% dplyr::select(Gene, HAL_Mean, FIL_Mean, regu)
  # # avg_exp_m <- as.data.frame(cbind(Gene=exp_m$Gene, 
  # #                                  HAL=as.numeric(rowMeans(exp_m[2:4])), 
  # #                                  FIL=as.numeric(rowMeans(exp_m[5:7])), 
  # #                                  regu=exp_m$regu))
  # # 
  # write.table(avg_exp_m,
  #             paste0(stage, "_", feature,"_", ctype, "_inner_exp_avg_methy.tsv"),
  #             sep="\t", quote = F)
  # 
  # 
  # melt_exp_m <- melt(avg_exp_m, id=c("Gene", "regu"))
  # colnames(melt_exp_m) <- c("Gene", "Regulation", "Sample", "Avg_mLevel")
  # 
  # melt_exp_m$Avg_mLevel <- as.numeric(melt_exp_m$Avg_mLevel)
  # melt_exp_m$Regulation <- as.factor(melt_exp_m$Regulation)
  # melt_exp_m <- melt_exp_m %>% unite("Sample_Regualtion", Sample:Regulation, remove = FALSE)
  # melt_exp_m$Sample_Regualtion <- as.factor(melt_exp_m$Sample_Regualtion)
  # 
  # # p <- ggplot(melt_exp_m, aes(x=Sample_Regualtion, y=Avg_mLevel, fill=Regulation)) +
  # #   geom_violin(width=1.4, trim=FALSE) +
  # #   geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  # #   scale_fill_viridis(discrete = TRUE)  +
  # #   labs(title= paste0(ctype, " in ", feature), y="Avg. Methylation Level", x="") +
  # #   scale_x_discrete(breaks=c("FIL_down","FIL_up","HAL_down", "HAL_up"),
  # #                    labels=c("FIL", "FIL", "HAL", "HAL")) +
  # #   theme(legend.position="top", legend.title = element_text(size=5), 
  # #         legend.text = element_text(size=4),
  # #         axis.title.x = element_blank()) +
  # #   coord_cartesian(ylim = c(-0.1, 1.2)) 
  # my_comparisons <- list( c("FIL_Mean_down", "FIL_Mean_up"),
  #                         c("HAL_Mean_down", "HAL_Mean_up"),
  #                         c("FIL_Mean_down", "HAL_Mean_down"),
  #                         c("FIL_Mean_up", "HAL_Mean_up") )
  # # my_comparisons <- list( c("FIL_Mean_down", "HAL_Mean_down"), 
  # #                         c("FIL_Mean_up", "HAL_Mean_up") )
  # p <- ggboxplot(melt_exp_m, x="Sample_Regualtion", y="Avg_mLevel", 
  #                color="Regulation", palette = "jco")+ 
  #   stat_compare_means(comparisons = my_comparisons, label="p.signif", color="red") + 
  #   stat_compare_means(label.y = 3, method = "t.test") + 
  #   coord_cartesian(ylim = c(-0.1, 1.5)) + 
  #   scale_x_discrete(breaks=c("FIL_down","FIL_up","HAL_down", "HAL_up"),
  #                    labels=c("FIL", "FIL", "HAL", "HAL")) +
  #   theme(legend.position="top", legend.title = element_text(size=5), 
  #         legend.text = element_text(size=4),axis.title.x = element_blank()) +
  #   labs(title= paste0(ctype, " in ", feature), y="Avg. Methylation Level", x="")
  # 
  # # solution 3
  # # stat_pvalue <- melt_exp_m %>% 
  # #   rstatix::wilcox_test(Avg_mLevel ~ Sample_Regualtion) %>%
  # #   filter(p < 0.05) %>% 
  # #   rstatix::add_significance("p") %>% 
  # #   rstatix::add_y_position() %>% 
  # #   mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))
  # # 
  # # p<- ggplot(melt_exp_m, aes(x=Sample_Regualtion, y=Avg_mLevel)) + 
  # #   geom_boxplot() +
  # #   ggpubr::stat_pvalue_manual(stat_pvalue, label = "p.signif") +
  # #   theme_bw(base_size = 16)
  # 
  # #addSmallLegend(p)
  
  # return(p)
}

# Stage 1 figures
# p_avg_gene_CG_1 <- generateBoxplotOverlayViolin("P1", "gene", "CG", exp_D1)
# print(p_avg_gene_CG_1)
# p_avg_gene_CHG_1 <- generateBoxplotOverlayViolin("P1", "gene", "CHG", exp_D1)
# print(p_avg_gene_CHG_1)
# p_avg_gene_CHH_1 <- generateBoxplotOverlayViolin("P1", "gene", "CHH", exp_D1)
# print(p_avg_gene_CHH_1)


# Includes all genes even no match
# and includes fdr and qvalues

for (feature in c("CDS", "promoter", "five_prime_UTR", "introns", "three_prime_UTR")) {
  for (ctype in c("CG", "CHG", "CHH")) {
    generateBoxplotOverlayViolin("P1", feature, ctype, "exp_D1")
    generateBoxplotOverlayViolin("P4", feature, ctype, "exp_D4")
  }
}








#####################################################3
# compare P1 and P4 within the same ecotypes

setwd(paste0(root_wd, "/figures/intraCompare"))

genotype="HAL"
stage="P1"
feature="CDS"
ctype="CG"

getFDRintraHAL <- function(genotype="HAL", feature, ctype){
  
  file1 <- paste0(root_wd, "/data/",genotype,"/P1_orthologous_joined_", ctype, "_", feature, ".txt")
  file2 <- paste0(root_wd, "/data/",genotype,"/P4_orthologous_joined_", ctype, "_", feature, ".txt")
  
  mc1 <- read_tsv(file1)
  mc2 <- read_tsv(file2)
  
  m_full_join <- mc1 %>% full_join(mc2, by=c("joined_ID"))
  
  # set NA to zero
  NAs <- m_full_join %>% filter(is.na(Gene.x) | is.na(Gene.y))
  
  # Create empty columns
  NAs$pval <- "NA"
  NAs$qvalue <- "NA"
  NAs$fdr <- "NA"
  NAs$fold_change <- "NAs"
  
  
  m_full_join <- m_full_join %>% filter(!is.na(Gene.x) & !is.na(Gene.y))
  m_full_join[is.na(m_full_join)] <- 0
  
  nm <- names(m_full_join)
  
  ffpval <- m_full_join %>%
    rowwise() %>%
    mutate(pval = t.test(c(HAL_r1.x, HAL_r2.x, HAL_r3.x),
                         c(HAL_r1.y, HAL_r2.y, HAL_r3.y))$p.value) %>%
    ungroup()
  
  qobj <- qvalue(p = ffpval$pval)
  
  # qvalues
  ffpval$qvalue <- qobj$qvalues
  
  #fdr
  ffpval$fdr <- qobj$lfdr
  
  # fold change
  ffpval$fold_change <- ffpval$HAL_Mean.y/ffpval$HAL_Mean.x
  
  ffpval <- rbind(ffpval, NAs)
  
  write.table(ffpval,
              paste0("intra_", genotype, "_", feature,"_", ctype, "_full_noexp_avg_methy_pval.tsv"),
              sep = "\t", quote = F)
}

for (ctype in c("CG", "CHG", "CHH")) {
  for (feature in c("promoter", "five_prime_UTR", "CDS", "introns", "three_prime_UTR")) {
    getFDRintraHAL("HAL", feature = feature, ctype = ctype)
  }
}



getFDRintraFIL <- function(genotype="FIL", feature, ctype){
  
  file1 <- paste0(root_wd, "/data/", genotype, "/P1_orthologous_joined_", ctype, "_", feature, ".txt")
  file2 <- paste0(root_wd, "/data/", genotype, "/P4_orthologous_joined_", ctype, "_", feature, ".txt")
  
  mc1 <- read_tsv(file1)
  mc2 <- read_tsv(file2)
  
  m_full_join <- mc1 %>% full_join(mc2, by=c("joined_ID"))
  
  # set NA to zero
  NAs <- m_full_join %>% filter(is.na(Gene.x) | is.na(Gene.y))
  
  # Create empty columns
  NAs$pval <- "NA"
  NAs$qvalue <- "NA"
  NAs$fdr <- "NA"
  NAs$fold_change <- "NAs"
  
  
  m_full_join <- m_full_join %>% filter(!is.na(Gene.x) & !is.na(Gene.y))
  m_full_join[is.na(m_full_join)] <- 0
  #m_full_join[m_full_join$FIL_r1.x == m_full_join$FIL_r2.x && m_full_join$FIL_r2.x == m_full_join$FIL_r3.x]
  
  nm <- names(m_full_join)
  
  ffpval <- m_full_join %>%
    rowwise() %>%
    mutate(pval = t.test(c(FIL_r1.x, FIL_r2.x, FIL_r3.x),
                         c(FIL_r1.y, FIL_r2.y, FIL_r3.y))$p.value) %>%
    ungroup()
  
  qobj <- qvalue(p = ffpval$pval)
  
  # qvalues
  ffpval$qvalue <- qobj$qvalues
  
  #fdr
  ffpval$fdr <- qobj$lfdr
  
  # fold change
  ffpval$fold_change <- ffpval$FIL_Mean.x/ffpval$FIL_Mean.y
  
  ffpval <- rbind(ffpval, NAs)
  
  write.table(ffpval,
              paste0("intra_", genotype, "_", feature,"_", ctype, "_full_noexp_avg_methy_pval.tsv"),
              sep = "\t", quote = F)
}


for (ctype in c("CG", "CHG", "CHH")) {
  for (feature in c("promoter", "five_prime_UTR", "CDS", "introns", "three_prime_UTR")) {
    getFDRintraFIL("FIL", feature = feature, ctype = ctype)
  }
}




###########################################

# Below are the ones generating figures, and no fdr and qvalues
# 
# p_avg_prom_CG_1 <- generateBoxplotOverlayViolin("P1", "promoter", "CG", exp_D1)
# print(p_avg_prom_CG_1)
# p_avg_prom_CHG_1 <- generateBoxplotOverlayViolin("P1", "promoter", "CHG", exp_D1)
# print(p_avg_prom_CHG_1)
# p_avg_prom_CHH_1 <- generateBoxplotOverlayViolin("P1", "promoter", "CHH", exp_D1)
# print(p_avg_prom_CHH_1)
# 
# p_avg_cds_CG_1 <- generateBoxplotOverlayViolin("P1", "CDS", "CG", exp_D1)
# print(p_avg_cds_CG_1)
# p_avg_cds_CHG_1 <- generateBoxplotOverlayViolin("P1", "CDS", "CHG", exp_D1)
# print(p_avg_cds_CHG_1)
# p_avg_cds_CHH_1 <- generateBoxplotOverlayViolin("P1", "CDS", "CHH", exp_D1)
# print(p_avg_cds_CHH_1)
# 
# p_avg_5UTR_CG_1 <- generateBoxplotOverlayViolin("P1", "five_prime_UTR", "CG", exp_D1)
# print(p_avg_5UTR_CG_1)
# p_avg_5UTR_CHG_1 <- generateBoxplotOverlayViolin("P1", "five_prime_UTR", "CHG", exp_D1)
# print(p_avg_5UTR_CHG_1)
# p_avg_5UTR_CHH_1 <- generateBoxplotOverlayViolin("P1", "five_prime_UTR", "CHH", exp_D1)
# print(p_avg_5UTR_CHH_1)
# 
# p_avg_3UTR_CG_1 <- generateBoxplotOverlayViolin("P1", "five_prime_UTR", "CG", exp_D1)
# print(p_avg_3UTR_CG_1)
# p_avg_3UTR_CHG_1 <- generateBoxplotOverlayViolin("P1", "five_prime_UTR", "CHG", exp_D1)
# print(p_avg_3UTR_CHG_1)
# p_avg_3UTR_CHH_1 <- generateBoxplotOverlayViolin("P1", "five_prime_UTR", "CHH", exp_D1)
# print(p_avg_3UTR_CHH_1)
# 
# p_avg_mRNA_CG_1 <- generateBoxplotOverlayViolin("P1", "three_prime_UTR", "CG", exp_D1)
# print(p_avg_mRNA_CG_1)
# p_avg_mRNA_CHG_1 <- generateBoxplotOverlayViolin("P1", "three_prime_UTR", "CHG", exp_D1)
# print(p_avg_mRNA_CHG_1)
# p_avg_mRNA_CHH_1 <- generateBoxplotOverlayViolin("P1", "three_prime_UTR", "CHH", exp_D1)
# print(p_avg_mRNA_CHH_1)
# 
# library(gridExtra)
# p11 <- grid.arrange(p_avg_gene_CG_1, p_avg_gene_CHG_1, p_avg_gene_CHH_1, nrow=1, ncol=3)
# p12 <- grid.arrange(p_avg_prom_CG_1, p_avg_prom_CHG_1, p_avg_prom_CHH_1, nrow=1, ncol=3)
# p13 <- grid.arrange(p_avg_cds_CG_1, p_avg_cds_CHG_1, p_avg_cds_CHH_1, nrow=1, ncol=3)
# p14 <- grid.arrange(p_avg_5UTR_CG_1, p_avg_5UTR_CHG_1, p_avg_5UTR_CHH_1, nrow=1, ncol=3)
# p15 <- grid.arrange(p_avg_3UTR_CG_1, p_avg_3UTR_CHG_1, p_avg_3UTR_CHH_1, nrow=1, ncol=3)
# p16 <- grid.arrange(p_avg_mRNA_CG_1, p_avg_mRNA_CHG_1, p_avg_mRNA_CHH_1, nrow=1, ncol=3)
# ggsave("D1_gene_body_methylation.pdf", device="pdf", plot=p11, width = 21, height=6)
# ggsave("D1_promoter_methylation.pdf", device="pdf", plot=p12, width = 21, height=6)
# ggsave("D1_cds_methylation.pdf", device="pdf", plot=p13, width = 21, height=6)
# ggsave("D1_5UTR_methylation.pdf", device="pdf", plot=p14, width = 21, height=6)
# ggsave("D1_3UTR_methylation.pdf", device="pdf", plot=p15, width = 21, height=6)
# ggsave("D1_mRNA_methylation.pdf", device="pdf", plot=p16, width = 21, height=6)
# 
# 
# 
# # ------------------------- Stage 4 -------------------------
# 
# p_avg_gene_CG_4 <- generateBoxplotOverlayViolin("P4", "gene", "CG", exp_D4)
# print(p_avg_gene_CG_4)
# p_avg_gene_CHG_4 <- generateBoxplotOverlayViolin("P4", "gene", "CHG", exp_D4)
# print(p_avg_gene_CHG_4)
# p_avg_gene_CHH_4 <- generateBoxplotOverlayViolin("P4", "gene", "CHH", exp_D4)
# print(p_avg_gene_CHH_4)
# 
# p_avg_prom_CG_4 <- generateBoxplotOverlayViolin("P4", "promoter", "CG", exp_D4)
# print(p_avg_prom_CG_4)
# p_avg_prom_CHG_4 <- generateBoxplotOverlayViolin("P4", "promoter", "CHG", exp_D4)
# print(p_avg_prom_CHG_4)
# p_avg_prom_CHH_4 <- generateBoxplotOverlayViolin("P4", "promoter", "CHH", exp_D4)
# print(p_avg_prom_CHH_4)
# 
# p_avg_cds_CG_4 <- generateBoxplotOverlayViolin("P4", "CDS", "CG", exp_D4)
# print(p_avg_cds_CG_4)
# p_avg_cds_CHG_4 <- generateBoxplotOverlayViolin("P4", "CDS", "CHG", exp_D4)
# print(p_avg_cds_CHG_4)
# p_avg_cds_CHH_4 <- generateBoxplotOverlayViolin("P4", "CDS", "CHH", exp_D4)
# print(p_avg_cds_CHH_4)
# 
# p_avg_5UTR_CG_4 <- generateBoxplotOverlayViolin("P4", "five_prime_UTR", "CG", exp_D4)
# print(p_avg_5UTR_CG_4)
# p_avg_5UTR_CHG_4 <- generateBoxplotOverlayViolin("P4", "five_prime_UTR", "CHG", exp_D4)
# print(p_avg_5UTR_CHG_4)
# p_avg_5UTR_CHH_4 <- generateBoxplotOverlayViolin("P4", "five_prime_UTR", "CHH", exp_D4)
# print(p_avg_5UTR_CHH_4)
# 
# p_avg_3UTR_CG_4 <- generateBoxplotOverlayViolin("P4", "five_prime_UTR", "CG", exp_D4)
# print(p_avg_3UTR_CG_4)
# p_avg_3UTR_CHG_4 <- generateBoxplotOverlayViolin("P4", "five_prime_UTR", "CHG", exp_D4)
# print(p_avg_3UTR_CHG_4)
# p_avg_3UTR_CHH_4 <- generateBoxplotOverlayViolin("P4", "five_prime_UTR", "CHH", exp_D4)
# print(p_avg_3UTR_CHH_4)
# 
# p_avg_mRNA_CG_4 <- generateBoxplotOverlayViolin("P4", "three_prime_UTR", "CG", exp_D4)
# print(p_avg_mRNA_CG_4)
# p_avg_mRNA_CHG_4 <- generateBoxplotOverlayViolin("P4", "three_prime_UTR", "CHG", exp_D4)
# print(p_avg_mRNA_CHG_4)
# p_avg_mRNA_CHH_4 <- generateBoxplotOverlayViolin("P4", "three_prime_UTR", "CHH", exp_D4)
# print(p_avg_mRNA_CHH_4)
# 
# library(gridExtra)
# P41 <- grid.arrange(p_avg_gene_CG_4, p_avg_gene_CHG_4, p_avg_gene_CHH_4, nrow=1, ncol=3)
# P42 <- grid.arrange(p_avg_prom_CG_4, p_avg_prom_CHG_4, p_avg_prom_CHH_4, nrow=1, ncol=3)
# P43 <- grid.arrange(p_avg_cds_CG_4, p_avg_cds_CHG_4, p_avg_cds_CHH_4, nrow=1, ncol=3)
# P44 <- grid.arrange(p_avg_5UTR_CG_4, p_avg_5UTR_CHG_4, p_avg_5UTR_CHH_4, nrow=1, ncol=3)
# P45 <- grid.arrange(p_avg_3UTR_CG_4, p_avg_3UTR_CHG_4, p_avg_3UTR_CHH_4, nrow=1, ncol=3)
# P46 <- grid.arrange(p_avg_mRNA_CG_4, p_avg_mRNA_CHG_4, p_avg_mRNA_CHH_4, nrow=1, ncol=3)
# ggsave("D4_gene_body_methylation.pdf", device="pdf", plot=P41, width = 21, height=6)
# ggsave("D4_promoter_methylation.pdf", device="pdf", plot=P42, width = 21, height=6)
# ggsave("D4_cds_methylation.pdf", device="pdf", plot=P43, width = 21, height=6)
# ggsave("D4_5UTR_methylation.pdf", device="pdf", plot=P44, width = 21, height=6)
# ggsave("D4_3UTR_methylation.pdf", device="pdf", plot=P45, width = 21, height=6)
# ggsave("D4_mRNA_methylation.pdf", device="pdf", plot=P46, width = 21, height=6)
