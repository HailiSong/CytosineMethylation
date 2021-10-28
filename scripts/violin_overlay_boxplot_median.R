rm(list=ls())
library(tidyverse)
library(gplots)
library(amap)
library(RColorBrewer)
library(fastcluster)
library(ggplot2)
library(reshape2)
library(viridis)

# stage = "P1"
# feature = "gene"
# ctype = "CG"
meth = "median"

root_wd <- "/Users/haili/Documents/WorkUT/mc_florence/bismarkCX/1to1orthologous"

setwd(paste0(root_wd, "/figures/", meth, "/avg_merged"))

# Read in Expression data
exp_D1 <- read_csv("/Users/haili/Documents/WorkUT/result/mc_florence/D1_Cont_Sig.csv")
names(exp_D1) <- c("Gene", names(exp_D1)[2:length(names(exp_D1))])

exp_D4 <- read_csv("/Users/haili/Documents/WorkUT/result/mc_florence/D4_Cont_Sig.csv")
names(exp_D4) <- c("Gene", names(exp_D4)[2:length(names(exp_D4))])

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
  fil_dir <- paste0(root_wd, "/mc_by_feature/", meth, "/", stage, "/FIL")
  hal_dir <- paste0(root_wd, "/mc_by_feature/", meth, "/", stage, "/HAL")
  
  # Read in mc files
  mc_fil <- read_tsv(paste0(fil_dir, "/", "orthologous_joined_",feature, "_", ctype, ".txt"))
  mc_hal <- read_tsv(paste0(hal_dir, "/", "orthologous_joined_",feature, "_", ctype, ".txt"))
  
  # Inner join
  inner_joined <- mc_hal %>% inner_join(mc_fil, by=c("joined_ID"))
  
  # selected columns
  inner_selected <- inner_joined %>% dplyr::select(joined_ID, HAL_P1, HAL_mc_P2, HAL_mc_P3,
                                                   FIL_P1, FIL_mc_P2, FIL_mc_P3)
  
  # Change to data frame
  inner_selected <- as.data.frame(inner_selected)
  rownames(inner_selected) <- inner_selected$joined_ID
  inner_selected$joined_ID <- NULL
  
  # Rename
  names(inner_selected) <- c("HAL_r1", "HAL_r2", "HAL_r3",
                             "FIL_r1", "FIL_r2", "FIL_r3")
  
  # Convert all columns to numeric
  inner_selected[, c(1:6)] <- sapply(inner_selected[, c(1:6)], as.numeric)
  
  # replace all NA with 0
  inner_selected[is.na(inner_selected)] = 0
  
  # remove rows having all zeros
  inner_selected <- inner_selected[rowSums(inner_selected[])>0,]
  
  ### Link with expression
  inner_selected <- as_tibble(tibble::rownames_to_column(inner_selected, "joined_ID"))
  inner_m <- inner_selected %>% separate(joined_ID, c("HAL_ID", "FIL_ID"), sep = "_")
  
  exp_m <- exp %>% inner_join(inner_m, by = c("Gene" = "HAL_ID") )
  
  # Add regulation data
  exp_m <- exp_m %>% mutate(regu =
                              case_when(log2FoldChange <=0 ~ "down", 
                                        log2FoldChange >0 ~ "up"))
  
  exp_m <- exp_m %>% dplyr::select(1, 9:15)
  
  # calculate average
  avg_exp_m <- as.data.frame(cbind(Gene=exp_m$Gene, 
                                   HAL=as.numeric(rowMeans(exp_m[2:4])), 
                                   FIL=as.numeric(rowMeans(exp_m[5:7])), 
                                   regu=exp_m$regu))
  
  write.table(avg_exp_m,
              paste0(stage, "_", feature,"_", ctype, "_inner_exp_avg_methy.tsv"),
              sep="\t", quote = F)
  
  
  melt_exp_m <- melt(avg_exp_m, id=c("Gene", "regu"))
  colnames(melt_exp_m) <- c("Gene", "Regulation", "Sample", "Avg_mLevel")
  
  melt_exp_m$Avg_mLevel <- as.numeric(melt_exp_m$Avg_mLevel)
  melt_exp_m$Regulation <- as.factor(melt_exp_m$Regulation)
  melt_exp_m <- melt_exp_m %>% unite("Sample_Regualtion", Sample:Regulation, remove = FALSE)
  melt_exp_m$Sample_Regualtion <- as.factor(melt_exp_m$Sample_Regualtion)
  
  p <- ggplot(melt_exp_m, aes(x=Sample_Regualtion, y=Avg_mLevel, fill=Regulation)) +
    geom_violin(width=1.4, trim=FALSE) +
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    scale_fill_viridis(discrete = TRUE)  +
    labs(title= paste0(ctype, " in ", feature), y="Avg. Methylation Level", x="") +
    scale_x_discrete(breaks=c("FIL_down","FIL_up","HAL_down", "HAL_up"),
                     labels=c("FIL", "FIL", "HAL", "HAL")) +
    theme(legend.position="top", legend.title = element_text(size=5), 
          legend.text = element_text(size=4),
          axis.title.x = element_blank()) +
    coord_cartesian(ylim = c(-0.1, 1.2)) 
  
  addSmallLegend(p)
  
  return(p)
}

# Stage 1 figures
p_avg_gene_CG_1 <- generateBoxplotOverlayViolin("P1", "gene", "CG", exp_D1)
print(p_avg_gene_CG_1)

p_avg_gene_CHG_1 <- generateBoxplotOverlayViolin("P1", "gene", "CHG", exp_D1)
print(p_avg_gene_CHG_1)

p_avg_gene_CHH_1 <- generateBoxplotOverlayViolin("P1", "gene", "CHH", exp_D1)
print(p_avg_gene_CHH_1)

p_avg_prom_CG_1 <- generateBoxplotOverlayViolin("P1", "promoter", "CG", exp_D1)
print(p_avg_prom_CG_1)

p_avg_prom_CHG_1 <- generateBoxplotOverlayViolin("P1", "promoter", "CHG", exp_D1)
print(p_avg_prom_CHG_1)

p_avg_prom_CHH_1 <- generateBoxplotOverlayViolin("P1", "promoter", "CHH", exp_D1)
print(p_avg_prom_CHH_1)

library(gridExtra)
p1 <- grid.arrange(p_avg_gene_CG_1, p_avg_gene_CHG_1, p_avg_gene_CHH_1,
             p_avg_prom_CG_1, p_avg_prom_CHG_1, p_avg_prom_CHH_1,
             nrow = 2, ncol=3)
ggsave("gene_body_promoter_D1_methylation.pdf", device="pdf", plot=p1, width = 16, height=10)




# Stage 4 figures
p_avg_gene_CG_4 <- generateBoxplotOverlayViolin("P4", "gene", "CG", exp_D4)
print(p_avg_gene_CG_4)

p_avg_gene_CHG_4 <- generateBoxplotOverlayViolin("P4", "gene", "CHG", exp_D4)
print(p_avg_gene_CHG_4)

p_avg_gene_CHH_4 <- generateBoxplotOverlayViolin("P4", "gene", "CHH", exp_D4)
print(p_avg_gene_CHH_4)

p_avg_prom_CG_4 <- generateBoxplotOverlayViolin("P4", "promoter", "CG", exp_D4)
print(p_avg_prom_CG_4)

p_avg_prom_CHG_4 <- generateBoxplotOverlayViolin("P4", "promoter", "CHG", exp_D4)
print(p_avg_prom_CHG_4)

p_avg_prom_CHH_4 <- generateBoxplotOverlayViolin("P4", "promoter", "CHH", exp_D4)
print(p_avg_prom_CHH_4)

library(gridExtra)
p4 <- grid.arrange(p_avg_gene_CG_4, p_avg_gene_CHG_4, p_avg_gene_CHH_4,
                  p_avg_prom_CG_4, p_avg_prom_CHG_4, p_avg_prom_CHH_4,
                  nrow = 2, ncol=3)
ggsave("gene_body_promoter_D4_methylation.pdf", device="pdf", plot=p4, width = 16, height=10)
# p <- ggplot(melt_exp_m, aes(x=Sample_Regualtion, y=Avg_mLevel, fill=Regulation)) 
# p <- p + geom_violin(width=1.4, trim=FALSE) 
# p <- p + geom_boxplot(width=0.1, color="grey", alpha=0.2) 
# p <- p + scale_fill_viridis(discrete = TRUE) 
# p <- p + labs(title= paste0(ctype, " in ", feature),
#        y="Avg. Methylation Level", x = "Samples")
# p <- p + scale_x_discrete(labels=c("FIL_down" = "FIL", "FIL_up" = "FIL",
#                                    "HAL_down" = "HAL", "HAL_down" = "HAL"))
#   # theme(
#   #   axis.text.x = element_blank(),
#   #   axis.text.y = element_blank(),
#   #   axis.ticks = element_blank())
# print(p)












