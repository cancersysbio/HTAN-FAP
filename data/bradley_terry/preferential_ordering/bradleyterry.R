#!/usr/bin/env Rscript

library(BradleyTerry2)
library(tidyverse)
library(stringr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
infile = args[1]
wd = args[2]

# infile = "/Users/ryanschenck/Dropbox/Projects/FAP_Project/mutation_timing/tmp.csv"
# wd = "/Users/ryanschenck/Dropbox/Projects/FAP_Project/mutation_timing"

setwd(wd)
df <- read.csv(infile)

if (length(df$gene_1)>1){
  # Prepare data for BT by determining the references based on median number of wins
  refcats = as.character(unique(df[df$n_1_before_2>=median(df$n_1_before_2),]$gene_1))

  # Get the unique genes
  genes = unique(c(as.character(df$gene_1), as.character(df$gene_2)))
  
  # Set the appropriate levels and data types
  df$gene_1 = factor(df$gene_1, levels=genes)
  df$gene_2 = factor(df$gene_2, levels=genes)
  
  # Initialize the BT model
  init = BTm(cbind(n_1_before_2, n_2_before_1), gene_1, gene_2, data = df, id = "gene_")
  # # Extract coefficient DF
  # coefEstimates = as.data.frame(summary(init)$coefficients)
  # 
  # # Remove gene_ prefix
  # coefEstimates$gene = str_remove(rownames(coefEstimates),"gene_")
  # initRefcat = genes[genes %in% coefEstimates$gene==FALSE]
  
  
  # Run BT model for each reference
  ret = list()
  for (i in 1:length(refcats)){
    # Exe BT model for each reference gene
    bt_coef <- update(init, refcat = refcats[i])
    
    # Extract coefficient DF
    coefEstimates = as.data.frame(summary(bt_coef)$coefficients)
    
    # Remove gene_ prefix
    coefEstimates$gene = str_remove(rownames(coefEstimates),"gene_")
    
    # Add informaiton about what reference was used
    coefEstimates$refcat = refcats[i]
    
    # add this BT model to the output
    ret[[i]] = coefEstimates
  }
  
  # Combine the different BT models with different references
  bt_coef = do.call(rbind, ret)
  
  # Get median coefficient estimates for each gene
  # bt_summarized = aggregate(bt_coef$Estimate, by=list(bt_coef$gene), median)
  # colnames(bt_summarized) <- c("gene", "Estimate")
  
  # Add the median p-values
  # bt_summarized$`Pr(>|z|)` <- aggregate(bt_coef$`Pr(>|z|)`, by=list(bt_coef$gene), median)$x
  
  # Add FDR and fonttype columnn
  bt_coef$fdr = p.adjust(bt_coef$`Pr(>|z|)`, method = "fdr")
  bt_coef = bt_coef%>%
    mutate(FDR.symbol = if_else(fdr<0.001,"***",
                                if_else(fdr<0.01,"**",
                                        if_else(fdr<0.05,"*",""))))
  
  # Reorder the table and refactor it so it comes out in easy too understand order
  # bt_summarized = bt_summarized%>%arrange(Estimate,fdr)
  # bt_summarized$gene = factor(bt_summarized$gene, levels=bt_coef$gene)
  bt_coef = bt_coef%>%
    mutate(FDR.symbol.position = if_else(Estimate<0,Estimate-1,Estimate+1))
  
  # Add rank order
  # bt_summarized$rankOrder = rank(-bt_summarized$Estimate)
  
  # Output the model run to file
  write.csv(bt_coef, "./tmp_out.csv", quote = FALSE, row.names = FALSE)
} else {
  print( "No genes or mutations that won. Check input tmp.csv")
  print( length(df$gene_1) )
}

