
require(RIdeogram)
library("tidyverse")

######## FIL ###############
genome <- "/Users/haili/Documents/WorkUT/annotation/Phallii_495_v3.0.genome.size.reduced"
gff3 <- "/Users/haili/Documents/WorkUT/annotation/Phallii_FIL_495_v3.1.gene.gff3"
input_dir <- "/home/haili/results/mc_florence/mc_by_features/100kb/FIL"

fil_karyotype <- read.table(genome, sep = "\t", header = T, stringsAsFactors = F)

gene_density <- GFFex(input = gff3, 
                      karyotype = genome, feature = "gene", window = 500000)
write.table(format(gene_density, scientific = F), "/Users/haili/Documents/WorkUT/annotation/FIL_gene_density_500k.tsv",
            sep = "\t", quote = F, row.names = F)

for (fil in c("P1-1", "P1-2", "P1-3")) {
  print(fil)
  mc <- paste0(input_dir, "/FIL-", fil, "_CX_report.CGCHGCHH.sort_bin.noZero_RIdeogram.txt")
  # Draw file label with marker for all C type in one
  Cs <- read.table(mc, sep = "\t", header = T, stringsAsFactors = F)
  ideogram(karyotype = fil_karyotype, overlaid = gene_density, label = Cs, label_type = "marker")
  convertSVG("chromosome.svg", device = "png")
  file.rename("chromosome.png", paste0("FIL-", fil, ".chromosome_marker.png"))
  
  for (c in c("CG", "CHG", "CHH")) {
    print(c)
    print(paste0(input_dir, "/FIL-", fil, "_CX_report.", c, ".sort_bin.noZero_marked.bed"))
    # Draw figure for each C type separated, labelled by marker
    mc_marker <- read.table(paste0(input_dir, "/FIL-", fil, "_CX_report.", c, ".sort_bin.noZero_marked.bed"),
                            sep = "\t", header = T, stringsAsFactors = F)
    ideogram(karyotype = fil_karyotype, overlaid = gene_density, label = mc_marker, label_type = "marker")
    convertSVG("chromosome.svg", device = "png")
    file.rename("chromosome.png", paste0("FIL-", fil, "_", c, ".chromosome_marker.png"))
    
    # Draw figure for each C type, labelled by line
    mc_line <- read.table(paste0(input_dir, "/FIL-", fil, "_CX_report.", c, ".sort_bin.noZero_line.bed"),
                          sep = "\t", header = T, stringsAsFactors = F)
    
    ideogram(karyotype = fil_karyotype, overlaid = gene_density, 
             label = mc_line, 
             label_type = "line", 
             colorset1 = c("#e5f5f9", "#99d8c9", "#2ca25f")) 
    convertSVG("chromosome.svg", device = "png")
    file.rename("chromosome.png", paste0("FIL-", fil, "_", c, ".chromosome_line.png"))
    
    # one-polygon lable
    ideogram(karyotype = fil_karyotype, overlaid = gene_density, 
             label = mc_line, 
             label_type = "polygon", 
             colorset1 = c("#e5f5f9", "#99d8c9", "#2ca25f")) 
    
    convertSVG("chromosome.svg", device = "png")
    file.rename("chromosome.png", paste0("FIL-", fil, "_", c, ".chromosome_polygon.png"))
    
    }
  
}




# ====================================================
# ====================================================


genome <- "/Users/haili/Documents/WorkUT/annotation/PhalliiHAL_496_v2.0.genome.size.reduced"
gff3 <- "/Users/haili/Documents/WorkUT/annotation/Phallii_HAL_496_v2.1.gene.gff3"
input_dir <- "/home/haili/results/mc_florence/mc_by_features/100kb/HAL"

hal_karyotype <- read.table(genome, sep = "\t", header = T, stringsAsFactors = F)

gene_density <- GFFex(input = gff3, 
                      karyotype = genome, feature = "gene", window = 500000)
write.table(format(gene_density, scientific = F), "/Users/haili/Documents/WorkUT/annotation/HAL_gene_density_500k.tsv",
            sep = "\t", quote = F, row.names = F)


for (hal in c("P1-1", "P1-2", "P1-3")) {
#for (hal in c("P1-3")) {
  print(hal)
  mc <- paste0(input_dir, "/HAL-", hal, "_CX_report.CGCHGCHH.sort_bin.noZero_RIdeogram.txt")
  # Draw file label with marker for all C type in one
  Cs <- read.table(mc, sep = "\t", header = T, stringsAsFactors = F)
  ideogram(karyotype = hal_karyotype, overlaid = gene_density, label = Cs, label_type = "marker")
  convertSVG("chromosome.svg", device = "png")
  file.rename("chromosome.png", paste0("HAL-", hal, ".chromosome_marker.png"))
  
  for (c in c("CG", "CHG", "CHH")) {
    print(c)
    print(paste0(input_dir, "/HAL-", hal, "_CX_report.", c, ".sort_bin.noZero_marked.bed"))
    # Draw figure for each C type separated, labelled by marker
    mc_marker <- read.table(paste0(input_dir, "/HAL-", hal, "_CX_report.", c, ".sort_bin.noZero_marked.bed"),
                            sep = "\t", header = T, stringsAsFactors = F)
    ideogram(karyotype = hal_karyotype, overlaid = gene_density, label = mc_marker, label_type = "marker")
    convertSVG("chromosome.svg", device = "png")
    file.rename("chromosome.png", paste0("HAL-", hal, "_", c, ".chromosome_marker.png"))
    
    # Draw figure for each C type, labelled by line
    mc_line <- read.table(paste0(input_dir, "/HAL-", hal, "_CX_report.", c, ".sort_bin.noZero_line.bed"),
                          sep = "\t", header = T, stringsAsFactors = F)
    
    ideogram(karyotype = hal_karyotype, overlaid = gene_density, 
             label = mc_line, 
             label_type = "line", 
             colorset1 = c("#e5f5f9", "#99d8c9", "#2ca25f")) 
    convertSVG("chromosome.svg", device = "png")
    file.rename("chromosome.png", paste0("HAL-", hal, "_", c, ".chromosome_line.png"))
    
    # one-polygon lable
    ideogram(karyotype = hal_karyotype, overlaid = gene_density, 
             label = mc_line, 
             label_type = "polygon", 
             colorset1 = c("#e5f5f9", "#99d8c9", "#2ca25f")) 
    
    convertSVG("chromosome.svg", device = "png")
    file.rename("chromosome.png", paste0("HAL-", hal, "_", c, ".chromosome_polygon.png"))
    
  }
  
}

































