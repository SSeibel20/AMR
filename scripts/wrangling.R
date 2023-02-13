#Wrangling AMR++ Output
#Samantha Seibel - last updated 20230213

#load packages
library(BiocManager)
library(dplyr)
library(stringr)


#read in data
out <- read.csv("amr-output/AMR_analytic_matrix.csv")

#drop extension from file name
colnames(out) <- str_split_i(colnames(out), "_", 1)

#break up gene information
genes <- as.data.frame(rownames(out))
genes <- as.data.frame(str_split_fixed(genes$`rownames(out)`, "\\|",5))
colnames(genes) <- c("MEG_ID", "Broadclass", "Class", "Mechanism", "Gene")

#merge gene info with counts
out2 <- cbind(out, genes[5])
out3 <- aggregate(.~Gene, data = out2, FUN = sum)
out4 <- out3[2:ncol(out3)]
rownames(out4) <-out3$Gene

#write.table(out4, file = "tables/countmatrix-cleaned.txt", sep = "\t")

#load packages for quick prelim plots
library(ggplot2)
library(data.table)

meltdf <- melt(as.data.table(out3))

ggplot(meltdf %>%
         filter(value = !="0"), aes(x=variable, fill = variable, y = value)) +
  geom_bar(stat =  "identity") +
  scale_y_log10() +
  facet_wrap(~Gene) +
  theme(legend.position = "none")

#heatmaps

ggplot(meltdf %>%
         filter(value = !="0"), aes(x=variable, fill = log(value), y = variable)) +
  geom_tile() +
  theme_bw() +
  coord_fixed() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
#save work
save.image("data/wrangling.RData")
