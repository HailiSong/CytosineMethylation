rm(list=ls())
library(tidyverse)
library(gplots)
library(amap)
library(RColorBrewer)
library(fastcluster)
library(ggplot2)
library(reshape2)

cols <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(50))
root_wd <- "/Users/haili/Documents/WorkUT/mc_florence/bismarkCX/1to1orthologous"

exp <- read_csv("/Users/haili/Documents/WorkUT/result/mc_florence/D1_Cont_Sig.csv")
names(exp) <- c("Gene", names(exp)[2:length(names(exp))])

stage = "P1"

for(meth in c("median")){
  print(meth)
  for (feature in c("CDS", "five_prime_UTR", "mRNA", "three_prime_UTR", "gene", "promoter")) {
    print(feature)
    for (ctype in c("CG", "CHG", "CHH")) {
      print(ctype)
      # set working directory
      setwd(paste0(root_wd, "/figures/", meth, "/", stage))
      
      # define directories of input files
      fil_dir <- paste0(root_wd, "/mc_by_feature/", meth, "/", stage, "/FIL")
      hal_dir <- paste0(root_wd, "/mc_by_feature/", meth, "/", stage, "/HAL")
      
      # Read in mc files
      mc_fil <- read_tsv(paste0(fil_dir, "/", "orthologous_joined_",feature, "_", ctype, ".txt"))
      mc_hal <- read_tsv(paste0(hal_dir, "/", "orthologous_joined_",feature, "_", ctype, ".txt"))
      
      # Join mc data of HAL and FIL by joined ID
      full_joined <- mc_hal %>% full_join(mc_fil, by=c("joined_ID"))
      
      # Select the columns of mc level
      full_selected <- full_joined %>% dplyr::select(joined_ID, HAL_P1, HAL_mc_P2, HAL_mc_P3,
                                                     FIL_P1, FIL_mc_P2, FIL_mc_P3)
      
      # change to data frame and set joined_ID as row name
      full_selected <- as.data.frame(full_selected)
      rownames(full_selected) <- full_selected$joined_ID
      full_selected$joined_ID <- NULL
      
      # change column name
      names(full_selected) <- c("HAL_P1", "HAL_P2", "HAL_P3",
                                "FIL_P1", "FIL_P2", "FIL_P3")
      
      # Convert all columns to numeric
      full_selected[, c(1:6)] <- sapply(full_selected[, c(1:6)], as.numeric)
      
      # replace all NA with 0
      full_selected[is.na(full_selected)] = 0
      
      
      # remove rows having all zeros
      full_selected <- full_selected[rowSums(full_selected[])>0,]
      
      
      # Save picture as pdf
      pdf(paste0(feature,"_", ctype, "_full_NA2zero_nozero.pdf"))
      heatmap.2(as.matrix(full_selected) ,scale='none',  na.rm=TRUE,
                trace="none", 
                symbreaks=F,
                col=cols,
                breaks=seq(0,1,length.out=51))
      dev.off()
      
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
      
      # Save picture as pdf
      pdf(paste0(feature,"_", ctype, "_inner_NA2zero_nozero.pdf"))
      heatmap.2(as.matrix(inner_selected) ,scale='none',  na.rm=TRUE,
                trace="none",
                symbreaks=F,
                col=cols,
                breaks=seq(0,1,length.out=51))
      dev.off()
      
      ### Link with expression
      inner_selected <- as_tibble(tibble::rownames_to_column(inner_selected, "joined_ID"))
      inner_m <- inner_selected %>% separate(joined_ID, c("HAL_ID", "FIL_ID"), sep = "_")
      
      exp_m <- exp %>% inner_join(inner_m, by = c("Gene" = "HAL_ID") )
      
      write.table(exp_m, paste0(feature,"_", ctype, "_inner_exp_methy.tsv"), sep="\t", quote = F)
      
      # Add regulation data
      exp_m <- exp_m %>% mutate(regu =
                                  case_when(log2FoldChange <=0 ~ "down", 
                                            log2FoldChange >0 ~ "up"))
      
      exp_m <- exp_m %>% dplyr::select(1, 9:15)
      
      melt_exp_m <- melt(exp_m, id=c("Gene", "regu"))
      colnames(melt_exp_m) <- c("Gene", "regu", "Sample", "mValue")
      
      #pdf(paste0(feature,"_", ctype, "_inner_exp_methy_boxplot.pdf"))
      p <- ggplot(melt_exp_m, aes(x=Sample, y=mValue, color = regu)) + 
        geom_boxplot() +
        labs(title= paste0(ctype, " in ", feature),
             y="Methylation Level", x = "Samples")
      #  coord_cartesian(ylim = c(0, 0.5))
      #dev.off()
      ggsave(paste0(feature,"_", ctype, "_inner_exp_methy_boxplot.pdf"), plot = p)
      
      #### violin plot
      p <- ggplot(melt_exp_m, aes(x=Sample, y=mValue, color = regu)) + 
        geom_violin(trim=FALSE) +
        labs(title= paste0(ctype, " in ", feature),
             y="Methylation Level", x = "Samples")
      ggsave(paste0(feature,"_", ctype, "_inner_exp_methy_violin.pdf"), plot = p)
      
      
      #### Scatter plots
      
      # calculate average
      avg_exp <- as.data.frame(cbind(Gene=exp_m$Gene, 
                                     HAL=as.numeric(rowMeans(exp_m[2:4])), 
                                     FIL=as.numeric(rowMeans(exp_m[5:7])), 
                                     regu=exp_m$regu)
      )
      #rownames(avg_exp) <- avg_exp$Gene
      avg_exp$Gene <- NULL
      p1 <- ggplot(avg_exp, aes(x=HAL, y=FIL, col=regu)) + 
        geom_point() +
        theme(
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          axis.line = element_line(colour = "black", 
                                   size = 0.5, linetype = "solid")) +
        geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
        labs(title= paste0("HAL vs. HAL of", ctype, " in ", feature))
      
      ggsave(paste0(feature,"_", ctype, "_inner_avg_methy_scatter.pdf"), plot = p1)
    }
  }
}










