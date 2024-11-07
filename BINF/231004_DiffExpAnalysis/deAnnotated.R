#!/usr/bin/env Rscript
# de.R
library(tximport)
library(readr)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(knitr)

# TODO: update constants for your machine
# Define constants
TESTING <- FALSE # Change to FALSE if using entire Samples set
RESULTS_DIR <- "/home/wong.co/BINF6309/m04"
AIPTASIA_DIR <- "/work/courses/BINF6309/AiptasiaMiSeq"

# for testing purposes - alternative samples table
testing_samples <- data.frame(Sample = c("Aip02", "Aip02", "Aip02", "Aip02"),
                              Menthol = c("Control", "Control", "Menthol", "Menthol"),
                              Vibrio = c("Control", "Vibrio", "Control", "Vibrio"))
# head(testing_samples)

# True script begins
tx2gene <- read.csv(file.path(RESULTS_DIR, "tx2gene.csv"))
# head(tx2gene)

if (TESTING) {
  print("***Running test with Aip02 only***")
  samples <- testing_samples
} else {
  samples <- read.csv(file.path(AIPTASIA_DIR, "Samples.csv"), header=TRUE)
}
# head(samples)


files <- file.path(RESULTS_DIR, "quant", samples$Sample, "quant.sf")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

dds <- DESeqDataSetFromTximport(txi, colData = samples, 
                                design = ~ Menthol + Vibrio)

dds$Vibrio <- relevel(dds$Vibrio, ref = "Control")
dds$Menthol <- relevel(dds$Menthol, ref = "Control")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

padj <- .05
minLog2FoldChange <- .5
dfAll <- data.frame()
# Get all DE results except Intercept, and "flatten" into a single file.
for (result in resultsNames(dds)){
  if(result != 'Intercept'){
    res <- results(dds, alpha=.05, name=result)
    dfRes <- as.data.frame(res)
    dfRes <- subset(subset(dfRes, select=c(log2FoldChange, padj)))
    dfRes$Factor <- result
    dfRes$ko <- rownames(dfRes)
    dfAll <- rbind(dfAll, dfRes)
  }
}
rownames(dfAll) <- NULL
# head(dfAll)

write.csv(dfAll, file=file.path(RESULTS_DIR, "dfAll.csv"))
# end of de.R script


# TODO: update file to filter for adjusted p-value (padj < 0.05) 
# AND merge pathways and pathway names (use your tables from Annotation)
# with the results, writing them to deAnnotated.csv

# Define directory containing annotation files
KEGG_DIR <- "/work/courses/BINF6309/data_BINF6309/Module4/Annotation"

# Import relevant annotation files containing ko, pathway, and description
path_df <- read.delim(file.path(KEGG_DIR, "path.txt"), header=FALSE) %>%
  rename(ko = V1,
         pathway = V2)
desc_df <- read.delim(file.path(KEGG_DIR, "ko"), header=FALSE) %>%
  rename(pathway = V1,
         description = V2)

# Merge annotation files together
annotations <- merge(path_df, desc_df, by="pathway")

df_annotated <- dfAll %>%
  filter(padj < 0.05) %>% # filter for padj
  left_join(annotations, by="ko") %>% # merge in annotations
  filter(!is.na(pathway)) %>% # get rid of NA
  subset(select=c(ko, pathway, description, log2FoldChange, padj, Factor)) %>% # reorder columns to match assignment guidelines
  arrange(ko, pathway) # sort by ko then pathway for readability

# Export file as csv
write.csv(df_annotated, file=file.path(RESULTS_DIR, "deAnnotated.csv"))

kable(df_annotated)