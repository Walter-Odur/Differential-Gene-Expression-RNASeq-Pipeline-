library(tximport)
library(DESeq2)
library(pheatmap)
library(enrichR)
library(tidyverse)
library(biomaRt)

setwd("LEARNING/PIPELINES/RNASeq/")

countdata <- read.csv("counts/feature_counts.txt", sep = "\t", row.names = 1, header = T, comment.char = "#")
View(countdata)
str(countdata)

# Remove unwanted columns
countdata[c('Chr','Start','End','Strand','Length')] <- NULL

View(countdata)

# Rename columns to remain with only sample ids
colnames(countdata)
#colnames(countdata) <- gsub("is this correct)
colnames(countdata) <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(countdata))
colnames(countdata) <- gsub("bams.", "", colnames(countdata))
colnames(countdata)


 