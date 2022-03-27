  setwd("/Users/aishwaryasundaresan/Documents/Laura_research/ATRX_Project/run30_23_combined_peak_calls/intersect_peaks/")
library(ggplot2)
library(gridExtra)
library(dplyr)
library(VennDiagram)
  library(UpSetR)
############
# read in the venn text
venn_table_df<-read.table(file = "WT_KO_given.TXT",header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE)

# get the venn categories
venn_categories<-colnames(venn_table_df)[!colnames(venn_table_df) %in% c("Total","Name")] 
cat("Venn categories are:\n"); venn_categories

# venn_categories
num_categories<-length(venn_categories)
cat("Num categories are:\n"); num_categories

# make a summary table
venn_summary<-venn_table_df[!colnames(venn_table_df) %in% venn_categories]
cat("Venn summary table is categories are:\n"); venn_summary
# venn_summary

# write summary table
#write.table(venn_summary,file = "venn_summary.tsv",quote = FALSE,row.names = FALSE)

# PAIRWISE VENN
cat("CREATING PAIR-WISE VENN DIAGRAM\n")
# area1
area_n1<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
# area2
area_n2<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
# n12
area_n12<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
venn <-draw.pairwise.venn(area1=area_n1,
                            area2=area_n2,
                            cross.area=area_n12,
                            category=c("WT26","ATRX_KO"),
                            fill=c('red','blue'),
                          scaled=F,
                          lty = "blank",
                          alpha=0.4,
                          cat.pos = c(295, 120),
                          cat.dist = 0.09,
                          cex = 2,
                          cat.cex = 2,
                          margin = 0.1
                          )
###########
#upSET
##SET UP THE PLOT ~~~~~~~ #
  # convert the summary table to a numeric/int vector, with element names as the combination names
  # # swap the | character with &; for passing to UpSet fromExpression
upset_expression <- setNames(venn_summary[['Total']], gsub("|","&",venn_summary[['Name']], fixed = TRUE))


upset(fromExpression(upset_expression), 
      nsets = num_categories, 
      order.by = "freq", 
      decreasing = T, 
      point.size=5,
      mainbar.y.label = "Overlapping Peaks", 
      sets.x.label = "Peaks per Category",
      sets.bar.color = "#56B4E9") 

