#Integration of AMR++ output into phyloseq object for downstream analysis 
#Samantha Seibel- 3/15/23 

#install phyloseq and vegan----
BiocManager::install("phyloseq")
BiocManager::install("vegan", force = TRUE)

#load packges----
library(BiocManager)
library(stringr)
library(tidyverse)
library(tidylog)
library(phyloseq)
library(vegan) # version 2.6-4
library(ggplot2)
library(readxl)

#read in necessary files: count matrix, gene info, metadata---- 
countsDF <- read.delim("tables/countmatrix-cleanedall.txt", sep = "\t") 
met <- readxl::read_excel("data/20230315_Ransom_AMR_Metadata.xlsx")
genes <- read.delim("tables/geneinfo-all.txt", sep = "\t") %>%
  select(-c(MEG_ID)) %>%
  unique()

#need counts as a matrix to create phyloseq
counts <- as.matrix(countsDF)

#need to set rows as taxa
otutab <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)

#need gene info as matrix
tax <- as.matrix(genes)

#row names need to be gene name
rownames(tax) <- genes$Gene

#make a taxa table
taxtab <- phyloseq::tax_table(tax)
taxa_names(taxtab)

samp <- phyloseq::sample_data(met)
rownames(samp) <- samp$`Sample-ID`

ps <- phyloseq::phyloseq(otutab, taxtab, samp)

#save the ps object 
saveRDS(ps, "data/rawps.rds")
