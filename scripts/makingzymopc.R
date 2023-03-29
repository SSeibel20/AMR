#Zymo PC AMR Profile - Sophia Kenney 3/27/23

#load packages
library(BiocManager)
library(tidyverse)
library(tidylog)
library(phyloseq)
library(microViz)
library(stringr)

#inspect Zymo file----
rep1 <- readxl::read_excel("D6300 AMR Profile.xlsx", sheet = "ExtPositive1.1.summary")
rep2 <- readxl::read_excel("D6300 AMR Profile.xlsx", sheet = "ExtPositive1.2.summary")
rep3 <- readxl::read_excel("D6300 AMR Profile.xlsx", sheet = "ExtPositive1.3.summary")


#which bacteria are there - same 7 + NA
bacall <- data.frame(table(rep1$species)) %>%
  rbind(data.frame(table(rep2$species))) %>%
  rbind(data.frame(table(rep3$species))) 

#combine into one
all <- rbind(rep1, rep2, rep3)
all[is.na(all)] <- "unknown"

#get list of bacteria
bacsp <- all$species %>%
  unique()

#aggregate megares gene info ---- 
#read in megares annotated file
meg <- read.csv("megares_full_annotations_with_notes_v2.00.csv")

genes <- data.frame(table(all$protein))

#add group info from megares manually 
genes$group <- rep("", nrow(genes))
test <- genes %>%
  mutate(group = case_when(
    str_detect(Var1, "Erm") ~ "ERMA",
    str_detect(Var1, "Lsa") ~ "LSA",
    #str_detect(Var1, "") ~ "", #something with the ant3/ant9
    str_detect(Var1, "O-phosphotransferase") ~ "APH3-PRIME", #this should be just prime not dprime but thats not in megares so idk
    str_detect(Var1, "BlaR1") ~ "BLAR",
    #str_detect(Var1, "MecR1") ~ "MEC?", #maybe mecA or something idk since R isnt there
    str_detect(Var1, "bifunctional") ~ "APH2-DPRIME",
    str_detect(Var1, "BlaEC") ~ "BLAEC",
    str_detect(Var1, "BlaZ") ~ "BLAZ",
    str_detect(Var1, "Bla1") ~ "BLA1",
    str_detect(Var1, "BlaI") ~ "BLAI",
    str_detect(Var1, "PDC") ~ "PDC",
    str_detect(Var1, "BPU") ~ "BPU",
    str_detect(Var1, "FosA") ~ "FOSA",
    str_detect(Var1, "FosB") ~ "FOSB",
    str_detect(Var1, "FosD") ~ "FOSD",
    str_detect(Var1, "FosX") ~ "FOSX",
    str_detect(Var1, "Lmr") ~ "LMRB",
    #str_detect(Var1, "") ~ "", #lmr something 
    str_detect(Var1, "MecI") ~ "MECI",
    str_detect(Var1, "PBP2a") ~ "MECA",
    str_detect(Var1, "MexX") ~ "MEXX",
    str_detect(Var1, "OqxB") ~ "OQXB",
    str_detect(Var1, "OXA") ~ "OXA",
    str_detect(Var1, "oxytetracycline") ~ "TET34",
    str_detect(Var1, "38") ~ "TET38", #the parenthesis probs weird
    str_detect(Var1, "DfrC") ~ "DFRC",
    str_detect(Var1, "CatB") ~ "CATB",
    str_detect(Var1, "Rph") ~ "RPH"
  ))

#just manually do these
test$group[4] <- "ANT9"
test$group[65] <- "TETK"
test$group[66] <- "TETL"
test$group[67] <- "TETM"
test$group[68] <- "TETS"

test[is.na(test)] <- "UNK"

#remove the three that wont match to megares anyway 
test <- test %>%
  filter(group != "UNK")

genes <- merge(test, meg, by = "group") %>%
  select(-c(headers, Freq, notes, RequiresSNPConfirmation.))%>%
  unique()

genes$group[3] <- "APH2-DPRIME"
colnames(genes)[2] <- "protein"

all <- merge(genes, all, by = "protein")

#write.table(all, "tables/zymopc.txt", sep = "\t")

#make taxtable object ----

appgenes <- read.delim("tables/geneinfo-all.txt", sep = "\t") %>%
  select(-c(MEG_ID)) %>%
  unique()

tax <- genes %>%
  select(type, class, mechanism, group) %>% 
  mutate_if(is.character, str_replace_all, ' ', '_') %>%
  unique()

colnames(tax) <- c("Broadclass", "Class", "Mechanism", "Gene")

taxtab <- as.matrix(tax)
rownames(taxtab) <- tax$Gene
taxtab <- phyloseq::tax_table(taxtab)
taxa_names(taxtab)


#make otu table and sample data----

#get average by species + gene (as if avging across replicates provided by zymo)
avgs <- all %>%
  select(species, group, read_counts) %>%
  aggregate(.~species + group, FUN = mean) 

avgs$read_counts <- round(avgs$read_counts, digits = 0)

#now that okay that I'm not inflating everything by adding batches, aggregate now by gene as we did with shotgun output
otu <-avgs %>%
  select(group, read_counts) %>%
  aggregate(.~group, FUN = sum) 

rownames(otu) <- otu$group
colnames(otu) <- c("Gene", "Zymo")
otu <- otu %>% select(Zymo)

counts <- as.matrix(otu)
otutab <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)


sam <- sample_data(data.frame(
  row.names = "Zymo",
  fakedata = 1
))


#combine----
ps <- phyloseq::phyloseq(otutab, taxtab, sam)

saveRDS(ps, file = "rdata/zymo-amr.rds")

#save work----
save.image("rdata/makingzymopc.RData")