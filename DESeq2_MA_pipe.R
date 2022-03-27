#!/usr/bin/env Rscript
library(reshape2)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(ggpubr)
library(DESeq2)
library(cummeRbund)

rm(list = ls())

# DE detecion ------------------------------------------------------------------


# import read count table


expr_readcount_combined <- read.table("combined_counts_Sara_WT-P300-CBPkd.txt", 
                                      header = T, row.names= 1, sep="\t")

# check the dimensions
dim(expr_readcount_combined)


# check the top 6 lines of this read count table
head(expr_readcount_combined)


# create sample information table
col_data <- data.frame(library = colnames(expr_readcount_combined),
                       sample = c(rep("KO",2),
                                  rep("WT",2)))

col_data$sample <- factor(col_data$sample, 
                           levels = c("KO", "WT"))


#check this sample information table
col_data


# create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = expr_readcount_combined,
                              colData = col_data,
                              design = ~ sample)

# perform DE detection
dds <- DESeq(dds)
res <- results(dds, 
               contrast = c("sample", "KO", "WT"),
               alpha = 0.05)


# check the result
head(res)


# print the column names of the result
mcols(res)$description


# print a quick summary
summary(res)

# add column with gene name

res1 <- read.table("gene_symbol_mapping.txt", 
                                      header = T, row.names= 1, sep="\t")
res$geneid<-res1$gene

# output results
write.table(cbind(gene_name = rownames(as.data.frame(res)), 
                  as.data.frame(res)), 
            file = "DESeq2_Sara_WT-P300-CBPkd.txt",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names =T)
          
            

pdf("WT-P300-CBPkd-MA.pdf")
# plot MA plot
#plotMA(res, ylim =c(-6, 6))
#with(subset(res, padj<.05 & (log2FoldChange)<(1*-1))
# plot MA plot with ggplot

#p <- MAplot(res,logMode=T,pseudocount=1,smooth=F)

ggmaplot(res, main = expression("WT" %->% "p300-CBP KD"),
   fdr = 0.05, fc = 2, size = 0.4,
   palette = c("#B31B21", "#1465AC", "darkgray"),
   genenames = as.vector(res$geneid),
   legend = "top", top = 20,
   font.label = c("bold", 11),
   font.legend = "bold",
   font.main = "bold",
   ggtheme = ggplot2::theme_minimal())
   
#p <- p + geom_point(aes(A, log2(M), colour = factor(ifelse(log2(M) < 2, 1,2))), size = 0.8) + geom_hline(yintercept = c(-2,2)) + theme(legend.position = "none") + scale_colour_manual(values = c("black","red"))


## add two lines
#abline(h = c(-1, 1), col = "steelblue", lty = "dashed")



