#Wrangling AMR++ Output
#Samantha Seibel - last updated 20230309

#load packages----
library(BiocManager)
library(plyr)
library(dplyr)
library(stringr)
library(rlist)


#read in and clean  data----
files <- Sys.glob("C:\\Users/saman/OneDrive - The Pennsylvania State University/Ganda Lab/Spring 2023/AMR++ Tutorial/AMR/amr-output/Ransom_practice/*.csv")

for (i in files) {
  filename <- str_split_fixed(str_split_fixed(paste0(i), "-", 3)[,3], ".csv", 2)[1]
  wd <- paste0(i)
  assign(filename, read.csv(wd))
}

dat<- lapply(Sys.glob("C:\\Users/saman/OneDrive - The Pennsylvania State University/Ganda Lab/Spring 2023/AMR++ Tutorial/AMR/amr-output/Ransom_practice/*.csv"), read.csv)
bovnames <- ls(pattern = "bov")
humnames <- ls(pattern = "hum")
batchnames <- list.append(bovnames, humnames)
filename2 <- paste(batchnames, ".2", sep = "")

for (i in 1:length(dat)) {
  assign(filename2[i],
         cbind(as.data.frame(dat[i]), as.data.frame(str_split_fixed(rownames(as.data.frame(dat[i])), "\\|", 5))[5])%>% #extract gene name from AMR++ output and add to count matrix
           aggregate(.~V5, FUN= sum))
}

#combine all matrices by gene ----
merged <- join_all(mget(grep(pattern = ".2", names(.GlobalEnv), value = TRUE, fixed = TRUE)), by = "V5", type = "full") %>% replace(., is.na(.), 0)
rownames(merged) <- merged$V5 #fix rownames so that they are the gene name
colnames(merged) <- trimws(colnames(merged), whitespace = "_.*") #fix colnames to drop file extension
merged_clean <- merged[2:ncol(merged)] #drop the gene name column so that its a count matrix with genes as rows and samples as columns

write.table(merged_clean, file = "C:\\Users/saman/OneDrive - The Pennsylvania State University/Ganda Lab/Spring 2023/AMR++ Tutorial/AMR/amr-output/Ransom_practice/tables/countmatrix-cleanedall.txt", sep = "\t") #save this matrix

#break up gene information for all batches (batch1 - batch n) - this will be useful later ----
genes <- data.frame()

for (i in 1:length(dat)) {
  genes <- rbind(genes, do.call(rbind, as.list(rownames(dat[[i]])))) %>% #add all rownames togther as a single column df
    unique() #remove redundant rows
}

colnames(genes)[1] <- "allmeg" #fix column name
genes <- as.data.frame(str_split_fixed(genes$allmeg, "\\|", 5)) #break up gene info
colnames(genes) <- c("MEG_ID", "Broadclass", "Class", "Mechanism", "Gene") #fix column names 

write.table(genes, file = "tables/geneinfo-all.txt", sep = "\t") #save this for later

#save work----
save.image("data/wrangling.RData")

#load packages for quick prelim plots----
library(ggplot2)
library(data.table)

meltdf <- melt(as.data.table(out3))

#bar plot
ggplot(meltdf %>%
         filter(value !="0"), aes(x=variable, fill = variable, y = value)) +
  geom_bar(stat =  "identity") +
  scale_y_log10() +
  facet_wrap(~Gene) +
  theme(legend.position = "none")

#heatmaps

ggplot(meltdf %>%
         filter(value !="0"), aes(x=Gene, fill = log(value), y = variable)) +
  geom_tile() +
  theme_bw() +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#save work----
save.image("data/wrangling.RData")



#earlier version: merge gene info with counts----
out2 <- cbind(out[2:ncol(out)], genes[5])
out3 <- aggregate(.~Gene, data = out2, FUN = sum)
out4 <- out3[2:ncol(out3)]
rownames(out4) <-out3$Gene

#write.table(genes, file = "tables/geneinfo-all.txt", sep = "\t")

#earlier version: read in and clean data batches----
out <- read.csv("amr-output/test_AMR_analytic_matrix.csv")
b1 <- read.csv("amr-output/Ransom_b1_AMR_analytic_matrix.csv")
b2 <- read.csv("amr-output/Ransom_b2_AMR_analytic_matrix.csv")
b3 <- read.csv("amr-output/Ransom_b3_AMR_analytic_matrix.csv")
b4 <- read.csv("amr-output/Ransom_b4_AMR_analytic_matrix.csv")

#do this for all batches, n being the batch n 
df1 <- cbind(b1, as.data.frame(str_split_fixed(rownames(b1), "\\|", 5))[5]) %>%
  unique() %>%
  aggregate(.~V5, FUN= sum)
df2 <- cbind(b2, as.data.frame(str_split_fixed(rownames(b2), "\\|", 5))[5]) %>%
  unique() %>%
  aggregate(.~V5, FUN= sum)
df3 <- cbind(b3, as.data.frame(str_split_fixed(rownames(b3), "\\|", 5))[5]) %>%
  unique() %>%
  aggregate(.~V5, FUN= sum)
df4 <- cbind(b4, as.data.frame(str_split_fixed(rownames(b4), "\\|", 5))[5]) %>%
  unique() %>%
  aggregate(.~V5, FUN= sum)

#install plyr
install.packages("plyr")
library(plyr)

#merge all batch gene info with counts
merged <- join_all(list(df1, df2, df3,dfn), by = "V5", type = "full")
merged[is.na(merged)] <- 0 #replace NAs 
rownames(merged) <- merged$V5 #fix rownames
colnames(merged) <- trimws(colnames(merged), whitespace = "_.*") #fix colnames
merged_clean <- merged[2:ncol(merged)]

write.table(merged_clean, file = "tables/countmatrix-cleanedall.txt", sep = "\t") #save this matrix

#break up gene information for alllll batches - this will be useful for later 
install.packages("rlist")
library("rlist")
genes <- as.data.frame(list.append(rownames(b1),
                                   rownames(b2),
                                   rownames(b3),
                                   rownames(b4))) %>%
  unique()

colnames(genes)[1] <- "allmeg"
genes <- as.data.frame(str_split_fixed(genes$allmeg, "\\|", 5))
colnames(genes) <- c("MEG_ID", "Broadclass", "Class", "Mechanism", "Gene")

write.table(genes, file = "tables/geneinfo-all.txt", sep = "\t") 
