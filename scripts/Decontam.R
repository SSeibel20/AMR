#decontam
#Samantha Seibel- 03/29/23

#install microViz and decontam----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

#load packages ----
#generally useful ones
library(BiocManager)
library(stringr)
library(tidyverse)
library(tidylog)

#ones we specifically need
library(phyloseq)
library(microViz)
library(vegan) # version 2.6-4
library(ggplot2)
library(decontam)

#load ps obj and other data---- 
countsDF <- read.delim("tables/countmatrix-cleanedall.txt", sep = "\t") 
met <- readxl::read_excel("data/20230315_Ransom_AMR_Metadata.xlsx")
genes <- read.delim("tables/geneinfo-all.txt", sep = "\t") %>%
  select(-c(MEG_ID)) %>%
  unique()

ps <- phyloseq::phyloseq(otutab, taxtab, samp)

#save the ps object 
saveRDS(ps, "data/rawps.rds")

#dummy code samples vs controls----
ps <- ps %>%
  ps_mutate(
    SampleBinary = if_else(str_detect(Sample.Type,"Control"), true = "Control", false = "Sample")
  )

#inspect # of counts per sample----
sam <- as.data.frame(sample_data(ps))
sam$GeneCounts <- sample_sums(ps)
sam <- sam[order(sam$GeneCounts),]
sam$Index <- seq(nrow(sam))

ggplot(sam, aes(x = Index, y = GeneCounts, color = Sample.Type)) +
  geom_point()

sample_data(ps)$is.neg <- sample_data(ps)$SampleBinary == "Control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")

#make a phyloseq object of only the contaminants in the controls
con <- rownames(contamdf.prev[contamdf.prev$contaminant == "TRUE",])
conps <- prune_taxa(rownames(ps@tax_table) %in% con, ps)
negps <- subset_samples(conps, is.neg == "TRUE")
con

#remove contaminant sequences
#looking in our ps object to look for the genes that are contaminated, and removing them exactly as they are
nocontam <- tax_select(ps, tax_list = con, strict_matches = TRUE, deselect = TRUE)

controls <- ps %>%
  ps_filter(SampleBinary == "Control", .keep_all_taxa = TRUE) 

controls <- sample_data(controls)$SampleID

#remove all controls
nocontrol <- ps_filter(nocontam, SampleBinary != "Control")

#removing SNP confirmation genes----
psfilt <- ps %>% tax_select(tax_list = "SNP", strict_matches = FALSE, deselect = TRUE)
taxa_names(psfilt)
taxa_names(nocontrol)

saveRDS(psfilt, "data/psfilt-presabs.rds")

#save work----
save.image("data/decontam.RData")
