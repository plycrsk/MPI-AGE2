## Analysis using DEseq2, controlling for AB diversity ##

library("DESeq2")
library("tidyverse")

## Loading Data ##

data <- read.csv(snakemake@input[[1]])
meta_data <- read.csv(snakemake@input[[2]])
#meta_data <- read.csv("meta_data_with_Qvalues_250721.csv")

output_path <- snakemake@output[[1]]

## Running Analysis ##

## Note: replace diversity order e.g. 'Q2.00' with desired diversity order ##
## in DESeqDataSetFromMatrix function, and name output file suitably ##

row.names(data) <- data$ensgene
data <- data[-1]

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = meta_data,
                              design = ~age +Q)

dds <- estimateSizeFactors(dds)
normalised <- counts(dds, normalized=TRUE)
write.csv(normalised, "data/Normalised_RNAseq_Counts_240821.csv")

dds <- DESeq(dds)
#contrast <- c('age', 'young', 'old')
#results <- results(dds, contrast = contrast)
results <- results(dds)
#results %>% data.frame() %>% View()
coeffs <- coef(dds)

write.csv(results, output_path)
#write.csv(coeffs, "DiffExp_age_plus_div_Q2.00_coefficients_150821.csv")

