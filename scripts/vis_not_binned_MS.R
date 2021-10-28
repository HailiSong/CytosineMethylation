library("tidyverse")


setwd("/Users/haili/Documents/WorkUT/mc_florence/bismarkCX/FIL_vs_HAL/P1")

fil_dir <- "/Users/haili/Documents/WorkUT/mc_florence/bismarkCX/FIL_vs_HAL/P1/FIL"
hal_dir <- "/Users/haili/Documents/WorkUT/mc_florence/bismarkCX/FIL_vs_HAL/P1/HAL"

features <- c("CDS")
ctype <- c("CG")

FIL <- read_tsv(paste0(fil_dir, "/", "orthologous_joined_", features, "_", ctype,".txt"))
head(FIL)
HAL <- read_tsv(paste0(hal_dir, "/", "orthologous_joined_", features, "_", ctype,".txt"))

joined <- FIL %>% inner_join(HAL, by=c("HAL_ID"))

write.table(joined, paste0("FIL_HAL_", features, "_", ctype, ".txt"), quote = FALSE)
filtered <- joined %>% select(FIL_ID.x:Fake_ID.x, strand.x, P1_median.x, P2_median.x, P3_median.x,
                              HAL_ID, Chr.y:Fake_ID.y, strand.y, P1_median.y, P2_median.y, P3_median.y)

# Join IDs of FIL and HAL
filtered$joined_ID <- paste0(filtered$FIL_ID.x,"_", filtered$HAL_ID)
filtered$FIL_ID.x <- NULL
filtered$HAL_ID <- NULL
median_mc <- filtered %>% select(P1_median.x, P2_median.x, P3_median.x, 
                                P1_median.y, P2_median.y, P3_median.y)

names(median_mc) <- c("P1_median.fil", "P2_median.fil", "P3_median.fil",
                      "P1_median.hal", "P2_median.hal", "P3_median.y")

df_median_mc <- as.data.frame(median_mc)
rownames(df_median_mc) <- df_median_mc$joined_ID

library(gplots)
heatmap.2(as.matrix(df_median_mc), na.rm=TRUE)












