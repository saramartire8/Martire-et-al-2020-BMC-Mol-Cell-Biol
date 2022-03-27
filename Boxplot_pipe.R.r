#load required packages
require("groHMM")
library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicFeatures")
#options(mc.cores=getCores(4))
#install.packages("vioplot")
#library(vioplot)


a <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/p300CBP_project_Sara/KNOCKOUT/ChIP-Seq/run72_20_12_16_H3K27ac/nodup_norm/WT_H3K27ac_rd.norm.bam")
a1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/p300CBP_project_Sara/ChiP-Seq/run43_2019_06_17/Sara_alignment/nodup_norm/p300scr_H3K27ac_2rep_S1.norm.bam") 
b <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/p300CBP_project_Sara/KNOCKOUT/ChIP-Seq/run72_20_12_16_H3K27ac/nodup_norm/P300KO1_H3K27ac_rd.norm.bam")
b1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/p300CBP_project_Sara/KNOCKOUT/ChIP-Seq/run72_20_12_16_H3K27ac/nodup_norm/P300KO4_H3K27ac_rd.norm.bam")

a.gr <- granges(a)
a1.gr <- granges(a1)
b.gr <- granges(b)
b1.gr <- granges(b1)


##Assign to variable
a <- a.gr 
a1 <- a1.gr 
b <- b.gr
b1 <- b1.gr 

 
##Combine replicates
##If there are replicates, only then.
a <- c(a.gr, a1.gr)
b <- c(b.gr, b1.gr)

#take either enhancer or refseq.bed
#enhancers <- import("p300_narrowpeaks.bed", format = "BED")
enhancers <- import("p300allWT_enhancers.bed", format = "BED")
#superenhancers <- import("super_enhancers_Whyte_mm10.bed", format = "BED")
#genes <- import("mm10_RefSeqgenes.bed", format = "BED")

## Find average
library_a <- NROW(a)
library_b <- NROW(b)

# Calculate RPKM
#Under enhancers
## Import files RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
rpkm_a <- countOverlaps(enhancers, a) / (width(enhancers)/1000 * library_a/1000000)
rpkm_b <- countOverlaps(enhancers, b) / (width(enhancers)/1000 * library_b/1000000)

rpkm_enhancers_p300 <- data.frame(rpkm_a, rpkm_b)
head(rpkm_enhancers_p300)
write.table(rpkm_enhancers_p300, file="rpkm_p300KO-enhancers.xls", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
pdf('boxplot_p300-enhancers_p300wtall.pdf')
boxplot(rpkm_enhancers_p300, col=(c("blue","red")), main="p300KOenhancers", ylab="RPKM",outline=FALSE, notch=TRUE, names=c("Ctrl","p300KO"), yaxt="n", cex.axis=1,las=2,lwd=4,lty=1, cex.axis=1,las=2,lwd=4,lty=1, ylim=c(0,5))
axis(2, at=c(0,1,2,3,4,5))

dev.off()

wilcox.test(rpkm_enhancers_p300$rpkm_a, rpkm_enhancers_p300$rpkm_b,conf.int=TRUE)


##End of Script##
