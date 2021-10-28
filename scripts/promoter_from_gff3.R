#### Convert gff3 to gRanges object and 
#### generate promoter regions

library(genomation)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

library(Homo.sapiens)

gff_FIL <- "/Users/haili/Documents/WorkUT/annotation/Phallii_FIL_495_v3.1.gene.gff3"
gff_HAL <- "/Users/haili/Documents/WorkUT/annotation/Phallii_HAL_496_v2.1.gene.gff3"

granges_FIL <- gffToGRanges(gff_FIL)
granges_HAL <- gffToGRanges(gff_HAL)

TxDb_FIL <- makeTxDbFromGFF(gff_FIL, format = "gff3")
TxDb_HAL <- makeTxDbFromGFF(gff_HAL, format = "gff3")

promoter_FIL <- promoters(genes(TxDb_FIL), upstream=2000, downstream=0)
promoter_HAL <- promoters(genes(TxDb_HAL), upstream=2000, downstream=0)

write.table( x = data.frame(promoter_FIL), 
             file = "/Users/haili/Documents/WorkUT/annotation/Phallii_FIL_495_v3.1.promoter.tsv",
             sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )

write.table( x = data.frame(promoter_HAL), 
             file = "/Users/haili/Documents/WorkUT/annotation/Phallii_HAL_496_v2.1.promoter.tsv",
             sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )



# Get introns from gff3 file ====================
# https://www.biostars.org/p/165226/

gtf <- makeTxDbFromGFF(gff_FIL) #change me!
exons <- exonsBy(gtf, by="gene")

#make introns
exons <- reduce(exons)
exons <- exons[sapply(exons, length) > 1]

introns <- lapply(exons, function(x) {
  #Make a "gene" GRange object
  gr = GRanges(seqnames=seqnames(x)[1], ranges=IRanges(start=min(start(x)),
                                                       end=max(end(x))), 
               strand=strand(x)[1])
  db = disjoin(c(x, gr))
  ints = db[countOverlaps(db, x) == 0]
  #Add an ID
  if(as.character(strand(ints)[1]) == "-") {
    ints$exon_id = c(length(ints):1)
  } else {
    ints$exon_id = c(1:length(ints))
  }
  ints
})
introns <- GRangesList(introns)
write.table(introns, "/Users/haili/Documents/WorkUT/annotation/Phallii_FIL_495_v3.1.gene.introns.txt",
            quote = F, sep = "\t")

# ---------------------------------------------

gtf <- makeTxDbFromGFF(gff_HAL) #change me!
exons <- exonsBy(gtf, by="gene")

#make introns
exons <- reduce(exons)
exons <- exons[sapply(exons, length) > 1]

introns <- lapply(exons, function(x) {
  #Make a "gene" GRange object
  gr = GRanges(seqnames=seqnames(x)[1], ranges=IRanges(start=min(start(x)),
                                                       end=max(end(x))), 
               strand=strand(x)[1])
  db = disjoin(c(x, gr))
  ints = db[countOverlaps(db, x) == 0]
  #Add an ID
  if(as.character(strand(ints)[1]) == "-") {
    ints$exon_id = c(length(ints):1)
  } else {
    ints$exon_id = c(1:length(ints))
  }
  ints
})
introns <- GRangesList(introns)
write.table(introns, "/Users/haili/Documents/WorkUT/annotation/Phallii_HAL_496_v2.1.gene.introns.txt",
            quote = F, sep = "\t")









