#Wrangling AMR++ Output
#Samantha Seibel - last updated 20230823

#load packages
library(BiocManager)
library(dplyr)
library(stringr)


#read in data
out <- read.csv("amr-output/AMR_analytic_matrix.csv")

#drop extension from file name
str_split_i(colnames(out), "_", 1)


#save work
save.image("data/wrangling.RData")
