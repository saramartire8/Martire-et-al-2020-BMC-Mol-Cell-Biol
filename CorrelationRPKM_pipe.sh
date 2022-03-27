#!/bin/tcsh

#SBATCH --job-name=bed_intersect
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --time=0-24:00:00
#SBATCH --output=bedintersect.%j.out
#SBATCH --error=bedintersect.%j.time
#SBATCH --mail-user=aishte@utsouthwestern.edu
#SBATCH --mail-type=ALL

module load deeptools/2.3.5  

###RUN boxplot (R)

require("groHMM")
library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicFeatures")

a <- readGAlignments("/your-path/129.merged.sorted.norm.bam")
b <- readGAlignments("/your-path/129.merged.sorted.norm.bam")


a.gr <- granges(a)
#a1.gr <- granges(a1)
b.gr <- granges(b)
#b1.gr <- granges(b1)

##Assign to variable
a <- a.gr 
#a1 <- a1.gr 
b <- b.gr
#b1 <- b1.gr 

 
##Combine replicates
##If there are replicates, only then.
#a <- c(a.gr, a1.gr)
#b <- c(b.gr, b1.gr)

#take either enhancer or refseq.bed

enhancers <- import("Activegenes_elongating_correct.bed", format = "BED")

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
write.table(rpkm_enhancers_p300, file="rpkm_GRO-Seq-Activegenes_elongating.xls", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
pdf('boxplot_GRO-Seq-Activegenes_elongating.pdf')
boxplot(rpkm_enhancers_p300, col=(c("blue","red")), main="GRO-Seq-Activegenes_elongating", ylab="RPKM",outline=FALSE, notch=TRUE, names=c("Ctrl","H33KO"), yaxt="n", cex.axis=1,las=2,lwd=4,lty=1, ylim=c(0,10))
axis(2, at=c(0,2,4,6,8,10))

dev.off()

wilcox.test(rpkm_enhancers_p300$rpkm_a, rpkm_enhancers_p300$rpkm_b,conf.int=TRUE)

##End of Script##

##RUN RPKM down up (R)
#This will take the output from boxplot script (rpkm.excel) and calculate fold change and log fold change


R/3.3

library(readxl)
library(dplyr)
library(tidyr)
library(gdata)
library(tidyverse)


RPKM.gene <- read.table("rpkm_Brd4-H3K27ac-enhancers.xls")
head(RPKM.gene)

colnames(RPKM.gene) <- c("WT_RPKM_Brd4", "KO_RPKM_Brd4", "WT_RPKM_Ac", "KO_RPKM_Ac")

"/" <- function(x,y) ifelse(y==0,x,base:::"/"(x,y))


RPKM.gene[,5] <- RPKM.gene[,2]/RPKM.gene[,1]
head(RPKM.gene)

RPKM.gene[,6] <- RPKM.gene[,4]/RPKM.gene[,3]
head(RPKM.gene)

RPKM.gene[,7] <- log2(RPKM.gene$V5)
head(RPKM.gene)
RPKM.gene[,8] <- log2(RPKM.gene$V6)
head(RPKM.gene)


#import promoters or enhancers file
enhancers<-read.table("p300merged_2500bp_distal_peaks.bed", header = FALSE, sep="\t",stringsAsFactors=FALSE)

#add  chr dan Start and End, merge with bed file/create new file

joint <- cbind(enhancers, RPKM.gene)
head(joint)

names(joint)[10] <- "log2Brd4"
names(joint)[11] <- "log2Ac"

write.table(joint, file="Brd4-H3K27ac-enhancers.bed", quote=F, sep="\t", row.names=F, col.names=F)

##RUN correlation (R) - either scatterplot/density or new contour plot

#Load the package
library(gdata)
library(gplots)
library(reshape)
library(R.utils)
library(dplyr)
library(readxl)
library(ggplot2)
library(plyr)
library(rowr)
library(devtools)
source_gist("524eade46135f6348140")
library(ggpmisc)
library(ggplot2)
library(tibble)
library(ggrepel)


#(Remove #NAME? MANUALLY FROM THE BED)
Joint <- read.table(file = "Brd4-H3K27ac-enhancers.bed", header= F, sep="\t",stringsAsFactors=FALSE)
head(Joint)

names(Joint) <- c("chr", "start", "end", "RPKM.WT.Brd4", "RPKM.KO.Brd4","RPKM.WT.Ac", "RPKM.KO.Ac", "fold.changeBrd4","fold.changeAc", "log2fold.changeBrd4","log2fold.changeAc")
head(Joint)


pdf("Ac-Brd4_scatterplot_enhancers-log2_nocolorscale.pdf")
ggplot(Joint, aes(x=log2fold.changeAc, y=log2fold.changeBrd4, color=log2fold.changeAc)) + geom_point(color='blue') + xlim(-6, 6) + ylim(-6, 6) + theme(aspect.ratio=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
cor.test(Joint$log2fold.changeAc, Joint$log2fold.changeBrd4, method = "spearman")
dev.off()

##Contour plot ##without geom_point()
pdf("Ac-Brd4_Contour plot_enhancers-log2_nocolorscale.pdf")
ggplot(data=Joint,aes(log2fold.changeAc,log2fold.changeBrd4)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..), geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") + stat_cor(method = "pearson") +
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + xlim(-2.5, 2.5) + ylim(-2.5, 2.5) + theme(aspect.ratio=1)
  guides(alpha="none")
dev.off()


