## Analysis using DEseq2, controlling for AB diversity ##

library("DESeq2")
library("tidyverse")

##Parameters##

q_values_raw <- c(0, 1, 1.5, 2, 3, 4)
q_values_str <- format(q_values_raw, nsmall=2)

## Loading Data ##
for (q in q_values_str){
data <- read.csv("raw_data/rna_data_raw_counts.csv")
meta_data <- read.csv(paste0("raw_data/meta_data_Q", q, ".csv")))
#meta_data <- read.csv("meta_data_with_Qvalues_250721.csv")



output_path <- paste0("output/DESeq/Q", q, ".csv"

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
write.csv(normalised, "data/Normalised_RNAseq_Counts.csv")

dds <- DESeq(dds)
#contrast <- c('age', 'young', 'old')
#results <- results(dds, contrast = contrast)
results <- results(dds)
#results %>% data.frame() %>% View()
coeffs <- coef(dds)

write.csv(results, output_path)
#write.csv(coeffs, "DiffExp_age_plus_div_Q2.00_coefficients_150821.csv")
}

