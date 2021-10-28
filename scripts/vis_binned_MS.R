library("tidyverse")
library(gplots)

setwd("/Users/haili/Documents/WorkUT/mc_florence/bismarkCX/FIL_vs_HAL/binned/P1")

orthologous <- read_tsv("HalvFil_1_to_1_orthologous")

fil_dir <- "./04.ext.dedup.FIL.mapping.default"
hal_dir <- "./04.ext.dedup.HAL.mapping.default"

features <- c("CDS")
ctype <- c("CG")

for (features in c("CDS", "five_prime_UTR", "mRNA", "three_prime_UTR")) {
  for (ctype in c("CG", "CHG", "CHH")) {
    print(paste(features, ctype))
    FIL <- read_tsv(paste0(fil_dir, "/", "orthologous_bygene_joined_", features, "_", ctype,".txt"))
    head(FIL)
    HAL <- read_tsv(paste0(hal_dir, "/", "orthologous_bygene_joined_", features, "_", ctype,".txt"))
    
    FIL$joined_ID <- paste(FIL$HAL_ID, FIL$FIL_ID)
    FIL$HAL_ID <- NULL
    FIL$FIL_ID <- NULL
    
    HAL$joined_ID <- paste(HAL$HAL_ID, HAL$FIL_ID)
    HAL$HAL_ID <- NULL
    HAL$FIL_ID <- NULL
    
    
    joined <- HAL %>% inner_join(FIL, by=c("joined_ID"))
    
    filtered <- joined %>% select(joined_ID, P1.x, P2.x, P3.x, P1.y, P2.y, P3.y)
    filtered <- filtered %>% rename(HAL.P1 = P1.x, HAL.P2 = P2.x, HAL.P3 = P3.x, 
                                    FIL.P1=P1.y, FIL.P2=P2.y, FIL.P3=P3.y)
    
    filtered <- as.data.frame(filtered)
    rownames(filtered) <- filtered$joined_ID
    filtered$joined_ID <- NULL
    
    pdf(paste0("heatmap_", features, "_", ctype, ".pdf"))
    heatmap.2(as.matrix(filtered) ,scale='none',  na.rm=TRUE,
              trace="none")
    dev.off()
  }
  
}


features <- c("CDS")

for (features in c("CDS", "five_prime_UTR", "mRNA", "three_prime_UTR")) {
  fil_CG <- read_tsv(paste0(fil_dir, "/", "orthologous_bygene_joined_", features, "_", "CG",".txt"))
  fil_CHG <- read_tsv(paste0(fil_dir, "/", "orthologous_bygene_joined_", features, "_", "CHG",".txt"))
  fil_CHH <- read_tsv(paste0(fil_dir, "/", "orthologous_bygene_joined_", features, "_", "CHH",".txt"))
  
  hal_CG <- read_tsv(paste0(hal_dir, "/", "orthologous_bygene_joined_", features, "_", "CG",".txt"))
  hal_CHG <- read_tsv(paste0(hal_dir, "/", "orthologous_bygene_joined_", features, "_", "CHG",".txt"))
  hal_CHH <- read_tsv(paste0(hal_dir, "/", "orthologous_bygene_joined_", features, "_", "CHH",".txt"))
  
  library(plyr)
  joined <- join_all(list(fil_CG, fil_CHG, fil_CHH, hal_CG, hal_CHG, hal_CHH), by="HAL_ID")
  
  # Rename columns
  names(joined) <- c("HAL_ID", "FIL_ID_1", "Fake_feature_ID_1", 
                     "FIL_CG_r1", "FIL_CG_r2", "FIL_CG_r3", 
                     "FIL_ID_2", "Fake_feature_ID_2", 
                     "FIL_CHG_r1", "FIL_CHG_r2","FIL_CHG_r3", 
                     "FIL_ID_3", "Fake_feature_ID_3",
                     "FIL_CHH_r1", "FIL_CHH_r2", "FIL_CHH_r3",
                     "FIL_ID_4", "Fake_feature_ID_4",
                     "HAL_CG_r1", "HAL_CG_r2", "HAL_CG_r3",
                     "FIL_ID_5", "Fake_feature_ID_5",
                     "HAL_CHG_r1", "HAL_CHG_r2", "HAL_CHG_r3",
                     "FIL_ID_6", "Fake_feature_ID_6",
                     "HAL_CHH_r1", "HAL_CHH_r2", "HAL_CHH_r3")
  
  # Set NA to Zero
  joined[is.na(joined)] <- 0
  mc <- joined %>% select(
    FIL_CG_r1, FIL_CG_r2, FIL_CG_r3, 
    FIL_CHG_r1, FIL_CHG_r2, FIL_CHG_r3, 
    FIL_CHH_r1, FIL_CHH_r2, FIL_CHH_r3,
    HAL_CG_r1, HAL_CG_r2, HAL_CG_r3,
    HAL_CHG_r1, HAL_CHG_r2, HAL_CHG_r3,
    HAL_CHH_r1, HAL_CHH_r2, HAL_CHH_r3)
  
  #mc$joined_ID <- paste0(mc$HAL_ID,"_", mc$FIL_ID_1)
  # rownames(mc) <- mc$joined_ID
  
  pdf(paste0("heatmap_joined_", features,"_HAL_FIL_CG_CHG_CHH.pdf"))
  heatmap.2(as.matrix(mc) ,scale='none', trace="none", 
            key.xlab="mCG level",
            density.info = 'none',key.title=NA,
            main=paste0("Binned(100bp) mC on HAL and FIL"),
            ylab = NULL)
  dev.off()
  
}





