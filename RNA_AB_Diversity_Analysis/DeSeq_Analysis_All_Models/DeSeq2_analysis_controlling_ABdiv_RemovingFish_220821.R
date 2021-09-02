## Analysis using DEseq2, controlling for AB diversity ##

library("DESeq2")
library("tidyverse")

## Loading Data ##

data <- read.csv("rna_data_raw_counts.csv", sep=',')
meta_data <- read.csv("meta_data_with_Qvalues_250721.csv")

data <- data[-12]
meta_data <- meta_data[-11,]
## Running Analysis ##

row.names(data) <- data$ensgene
data <- data[-1]

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = meta_data,
                              design = ~age +Q4.00)

dds <- DESeq(dds)
#contrast <- c('age', 'young', 'old')
#results <- results(dds, contrast = contrast)
results <- results(dds)
coeffs <- coef(dds)

write.csv(results, "Results/DiffExp_Analysis/RemovingFish/DiffExp_age_plus_div_Q4.00_NoFish3_220821.csv")
write.csv(coeffs, "Results/DiffExp_Analysis/RemovingFish/DiffExp_age_plus_div_Q4.00_coefficients_NoFish3_220821.csv")

## summarising results ##

#padj.cutoff <- 0.05

# Create a tibble of results
#results <- results %>%
#  data.frame() %>%
#  rownames_to_column(var="gene") %>% 
#  as_tibble()

#sigOE <- results %>%
#  filter(padj < padj.cutoff)

#write.csv(sigOE, "SignificantGenes_age_plus_div_Q4.00_170821.csv")
