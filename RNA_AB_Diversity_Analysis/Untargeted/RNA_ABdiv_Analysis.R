## Add introductory explanation ##

## Import required libraries ##
library("DESeq2")
library("biomaRt")
library("plyr")
library("DOSE")
library("dplyr")
library("tidyr")
library("broom")
library("data.table")
library("clusterProfiler")
library("org.Dr.eg.db")
library("org.Hs.eg.db")
library("topGO")
library("Rgraphviz")
library("GOplot")

## Import read counts and meta data ## 

countData <- read.csv("rna_data_raw_counts.csv", header=TRUE, sep=",")

metaData <- read.csv("rna_seq_metadata.csv", header=TRUE, sep=",")

## Map to human orthologs so that I can summarise many:one mappings ##

## BioMart - Adding orthologue data - Mapping killifish -> human ##

ensembl <- useMart("ensembl")

gene_list = countData$ensgene

human = useDataset("hsapiens_gene_ensembl", mart=ensembl)
killifish = useMart("ensembl", dataset="nfurzeri_gene_ensembl")

human_ortho <- getLDS(attributes=c("ensembl_gene_id"), filters="ensembl_gene_id",
                      values=gene_list,
                      mart=killifish,attributesL=c("hgnc_symbol",
                                                   "ensembl_gene_id"), martL=human)

write.csv(human_ortho,"killifish_human_ortho.csv", row.names = FALSE)

## Code for mapping killifish_human_ortho onto dataset ##

human_orthos <- c(human_ortho$Gene.stable.ID.1)
names(human_orthos) <- human_ortho$Gene.stable.ID

countData$human_ortho <- countData$ensgene

countData$human_ortho <- plyr::mapvalues(countData$human_ortho, from=c(human_ortho$Gene.stable.ID), to=c(human_ortho$Gene.stable.ID.1))

## Summarise those transcripts that have duplicate matches ## 
# Some genes are duplicated, due to many:one mapping. Take mean of these and remove duplicates #

rownames(countData) <- countData$ensgene
countData <- countData[-1]
countData <- countData %>% group_by(human_ortho) %>% summarise_all(funs(mean))

countData <- as.data.frame(countData)

countData[, c(2:8)] <- sapply(countData[, c(2:8)], as.integer)

## Normalise with DeSeq2 ##

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~age, tidy = TRUE)

dds <- estimateSizeFactors(dds)
normalised_rna <- counts(dds, normalized=TRUE)

# log transform # 

normalised_rna <- log2(normalised_rna+1)

write.csv(normalised_rna, "DEseq2_normalised_rna_counts_170621.csv")

## Linear Regression lm(transcript ~ gene + age) ##

rna <- t(normalised_rna)
rna <- as.data.frame(rna)

# need to add diversity and age columns #

metadata <- read.csv("rna_seq_metadata.csv", header=TRUE, sep=",")
rna <- merge(rna, metadata, by.x=0, by.y='id')

# creating regression model controlling for age #

model_1 = list()
#names(rna) <- colnames(rna)

for (i in names(rna)[-1]){
  storage <- lm(get(i) ~ div + age, data = rna)
  tidy(storage)
  glance(storage)
  model_1[[i]] <- glance(storage)
}

# formatting output data from regression analysis #
store <- model_1
model_1 <- store
model_1 <- do.call(rbind, model_1)
model_1 <- cbind(GENE = rownames(model_1), model_1)
rownames(model_1) <- 1:nrow(model_1)
#model_1 <- t(model_1)
#model_1 <- lapply(model_1,as.numeric)
data <- model_1
#sorted_f <- model_1[order(model_1$r.squared, decreasing=TRUE),]
#test <- model_1[order(model_1$r.squared, decreasing=TRUE)]

# fdr correction for p-values, not required for this particular analysis #
# fdr correction done with gseGO later #

data$FDR_p_values = p.adjust(data$p.value, "fdr")

write.csv(data,"sorted_regression_table_model1_170621.csv", row.names = FALSE)

# Removing non-mapped genes #

data$GENE[!startsWith(data$GENE , "ENSG0")] <- NA
data <- na.omit(data)

write.csv(data,"killifish_human_ortho_regresson_table_170601.csv", row.names = FALSE)

## gseGO - enrichment analysis of transcripts correlated with AB diversity, based on kolmogorov smirnov test ##

# generating input required for gseGO - a list of genes, ranked, e.g by p.value #

data <- dplyr::select(data,c('GENE', 'p.value'))
data$p.value <- as.numeric(data$p.value)
data$p.value <- -log(data$p.value)

# Generate ranked geneList object for gseGO #

geneList = data[,2]
names(geneList) = as.character(data[,1])
geneList = sort(geneList, decreasing = TRUE)

# gseGO function #

set.seed(4321)
ego3 <-gseGO(
  geneList,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

ego_summary <- slot(ego3, "result")

## Visualise output of gseGO ##

GO_vis <- plotGOgraph(ego3)
GO_vis_plot <- plot(GO_vis$complete.dag)

write.csv(ego_summary, "gseGO_model1_humanOrtho_FDR_200620.csv")
