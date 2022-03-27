#!/bin/tcsh

#SBATCH --job-name=bed_intersect
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --time=0-24:00:00
#SBATCH --output=bedintersect.%j.out
#SBATCH --error=bedintersect.%j.time
#SBATCH --mail-user=aishte@utsouthwestern.edu
#SBATCH --mail-type=ALL


##RUN multiBamSummary with your BED of interest

module load deeptools/2.3.5  

##multiBamSummary bins --bamfiles file1.bam file2.bam -o results.npz
##multiBamSummary BED-file --BED selection.bed --bamfiles file1.bam file2.bam -o results.npz

multiBamSummary BED-file --BED p300merged_2500bp_distal_peaks.bed --bamfiles /your-path/129_BRD4_S7.norm.bam /your-path/282_BRD4_S8.norm.bam /your-path/129_H3K27ac_SI_rep3_norm.bam /your-path/282_H3K27ac_SI_rep3_norm.bam -o Brd4-Ac_multibamsummary_enhancers.npz --outRawCounts Brd4-Ac_multibamsummary_enhancers.txt

multiBamSummary BED-file --BED p300merged_2500bp_prom_proximal_peaks.bed --bamfiles /your-path/129_BRD4_S7.norm.bam /your-path/282_BRD4_S8.norm.bam /your-path/129_H3K27ac_SI_rep3_norm.bam /your-path/282_H3K27ac_SI_rep3_norm.bam -o Brd4-Ac_multibamsummary_promoters.npz --outRawCounts Brd4-Ac_multibamsummary_promoters.txt


##LOAD R and RUN this in R
###RUN RPKM down up (R) This will take the output from multiBamSummary (.txt) and calculate fold change and log fold change

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

Joint <- read.table(file = "Brd4-Ac_multibamsummary_promoters.txt", header= F, sep="\t",stringsAsFactors=FALSE)
head(Joint)

names(Joint) <- c("chr", "start", "end", "WT.Brd4", "KO.Brd4","WT.Ac", "KO.Ac")
head(Joint)


Joint[,"FCBrd4"] <- Joint[,5]/Joint[,4]
head(Joint)

Joint[,"FCGRO"] <- Joint[,7]/Joint[,6]
head(Joint)

Joint[,"log2FCBrd4"] <- log2(Joint[, "FCBrd4"])
head(Joint)
Joint[,"log2FCGRO"] <- log2(Joint[, "FCGRO"])
head(Joint)

#SCATTERPLOT
pdf("Brd4-GRO_mm10_RefSeq_promoter-scatterplot.pdf")
ggplot(Joint, aes(x=log2FCGRO, y=log2FCBrd4, color=log2FCGRO)) + geom_point(color='blue') + xlim(-6, 6) + ylim(-6, 6) + theme(aspect.ratio=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
cor.test(Joint$log2fold.changeAc, Joint$log2fold.changeBrd4, method = "spearman")
dev.off()

##Contour plot ##without geom_point()

pdf("Brd4-GRO_mm10_RefSeq_promoter_Contour_plot.pdf")
ggplot(data=Joint,aes(log2FCGRO,log2FCBrd4)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..), geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") + stat_cor(method = "pearson") +
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + xlim(-2.5, 2.5) + ylim(-2.5, 2.5) + theme(aspect.ratio=1)
  guides(alpha="none")
dev.off()


