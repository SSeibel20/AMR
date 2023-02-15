#Exploratory plots
#Samantha Seibel - last updated 02/15/2023

#load packages
library(BiocManager)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

#load data from last time
load("data/wrangling_Sophia.RData")

#load Metadata 
meta <- read.delim("C:\\Users/saman/OneDrive - The Pennsylvania State University/Ganda Lab/Spring 2023/AMR++ Tutorial/AMR/data/rhamr-samples.txt")

#merge metadata with melted df
#first we need to correct the column names in meltdf

colnames(meltdf) <- c("Gene", "rh_ID", "counts")

meltdf2 <- merge(meltdf, meta, by = "rh_ID")

#add gene information, but need to remove duplicates of the genes, otherwise it gets confused because the rows are the same except for the MegID col
genefilt <- genes[2:5] %>%
  unique() #drop duplicate rows

meltdf2 <- merge(meltdf2, genefilt, by = "Gene")

#options for exploratory plots

#heatmap
ggplot(meltdf2 %>%
         filter(counts != "0"), aes(x=Gene, fill = Class, y = rh_ID)) +
  geom_tile() +
  theme_bw() +
  coord_fixed() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2))

#bar graph by Class
ggplot(meltdf2 %>%
         filter(counts != "0") %>%
         filter(Class == "Aminoglycosides"), aes(x=rh_ID, fill = rh_ID, y = counts)) +
  geom_bar(stat =  "identity") +
  scale_y_log10() +
  facet_wrap(~Gene) +
  theme(legend.position = "none")

#bar graoh by rh_ID
ggplot(meltdf2 %>%
        filter(counts != "0") %>% 
        filter(rh_ID == "rh02"), aes(x=rh_ID, fill = Class, y = counts)) +
  geom_bar(stat =  "identity") +
  scale_y_log10() +
  facet_wrap(~Broadclass) +
  theme(legend.position = "none")


#save work
save.image("data/preimplots.RData")
