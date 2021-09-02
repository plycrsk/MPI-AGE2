## Analysis using DEseq2, controlling for AB diversity ##

library("DESeq2")

## Loading Data ##

data <- read.csv("rna_data_raw_counts.csv", sep=',')
meta_data <- read.csv("meta_data_with_Qvalues_250721.csv")
Q_values <- read.csv("Hill_spectra_allQ_150721.csv")

## Adding Q values to meta data ##

#meta_data$Q4.00 <- plyr::mapvalues(meta_data$ID, from=c(Q_values$INDIVIDUAL), to=c(Q_values$X4.0))
#meta_data <- meta_data[-c(1,2,17,18),]

## Meta data row order needs to match count data column order ##

#meta_data <- meta_data[order(meta_data$sra),]

#write.csv(meta_data, "meta_data_with_Qvalues_250721.csv")

## 'Standard' Deseq2 analysis, comparing old v young ##

row.names(data) <- data$ensgene
data <- data[-1]

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = meta_data,
                              design = ~ age)

dds <- DESeq(dds)
contrast <- c('age', 'old', 'young')
res <- results(dds, contrast = contrast)
res

write.csv(res, "DEseq2_diffexp_old_v_young_260721.csv")

## Now running DiffExp whilst controlling for Q values ##

## Loading Data ##

data <- read.csv("rna_data_raw_counts.csv", sep=',')
meta_data <- read.csv("meta_data_with_Qvalues_250721.csv")
Q_values <- read.csv("Hill_spectra_allQ_150721.csv")

## Running Analysis ##

row.names(data) <- data$ensgene
data <- data[-1]

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = meta_data,
                              design = ~age + Q4.00)

dds <- DESeq(dds)
#contrast <- c('age', 'old', 'young')
#results <- results(dds, contrast = contrast)
results <- results(dds)
write.csv(results, "DiffExp_Q4.00_270721.csv")
