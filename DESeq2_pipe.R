#!/usr/bin/env Rscript
library(DESeq2)

rm(list = ls())


# DE detecion ------------------------------------------------------------------


# import read count table


expr_readcount_combined <- read.table("combined_counts_Sara_WT-CBPkd.txt", 
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


# output results
write.table(cbind(gene_name = rownames(as.data.frame(res)), 
                  as.data.frame(res)), 
            file = "DESeq2_Sara_WT-CBPkd.txt",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names =T)

pdf("WT-CBPkd.pdf")
# plot MA plot
plotMA(res, ylim =c(-6, 6))


## add two lines
abline(h = c(-1, 1), col = "steelblue", lty = "dashed")



# Hierarchical clustering ------------------------------------------------------


expr_readcount_norm <- 
    t(t(expr_readcount_combined) / 
          estimateSizeFactorsForMatrix(expr_readcount_combined))

sample_distance <- dist(t(log(expr_readcount_combined + 1, 10)),
                        method = "euclidean")

plot(hclust(sample_distance,
            method = "average"))


# Reduce dimensions by PCA ------------------------------------------------------


calc_cpm <- function(expr) {
    cpm <- t(t(expr) / colSums(expr)) * 1e+06
    
    if (is.data.frame(expr)) {
        cpm <- as.data.frame(cpm)
    }
    dimnames(cpm) <- dimnames(expr)
    return(cpm)
}


expr_cpm_combined <- calc_cpm(expr_readcount_combined)


pc_transcripts_var <- apply(expr_cpm_combined, 1, var)
transcripts_use <- rownames(expr_cpm_combined)[pc_transcripts_var > 0]
length(transcripts_use)


expr_use <- log10(expr_cpm_combined[transcripts_use, ] + 1)


pca_out <- prcomp(t(expr_use),
                  center = TRUE, scale = TRUE)


plot(pca_out$x[, "PC1"], pca_out$x[, "PC2"],
     xlab = "PC1",
     ylab = "PC2",
     main = "Principal component analysis", 
     col = "red")

text(pca_out$x, labels = rownames(pca_out$x))

#Volcano Plot
res <- read.table("DESeq2_Sara_WT-CBPkd.txt", header=TRUE)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-4,4)))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

library(calibrate)
library(MASS)
with(subset(res, padj<.05 & abs(log2FoldChange)>2), textxy(log2FoldChange, -log10(pvalue), labs=gene, cex=.8))
