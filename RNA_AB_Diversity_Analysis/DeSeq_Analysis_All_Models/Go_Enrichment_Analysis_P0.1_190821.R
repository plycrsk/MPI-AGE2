# Load libraries
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)

data <- read.csv("Results/DiffExp_Analysis/DiffExp_age_plus_div_Q4.00_150821.csv", sep=',')

colnames(data)[1] <- "gene"

ortho <- read.csv("killifish_human_ortho.csv")

## mapping killifish transcripts to human orthologues ## 

data$gene <- plyr::mapvalues(data$gene, from=c(ortho$Gene.stable.ID), to=c(ortho$Gene.stable.ID.1))

data$gene[!startsWith(data$gene , "ENSG0")] <- NA
data <- na.omit(data)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
allOE_genes <- as.character(data$gene)

## Extract significant results
sigOE <- dplyr::filter(data, padj < 0.10)
sigOE_genes <- as.character(sigOE$gene)

## Run GO enrichment analysis 
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "fdr",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05, 
                readable = TRUE)

ego_summary <- slot(ego, "result")
save(ego, file="Results/GoEnrichment_Analysis/ego_GoEnrich_Simplified_age_plus_div_Q4.00_p0.1_190821.rda")

ego_simple <- clusterProfiler::simplify(
  ego,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)

ego_simple_summary  <- slot(ego_simple, "result")

write.csv(ego_simple_summary, "Results/GoEnrichment_Analysis/GoEnrich_Simplified_age_plus_div_Q4.00_p0.1_190821.csv")
