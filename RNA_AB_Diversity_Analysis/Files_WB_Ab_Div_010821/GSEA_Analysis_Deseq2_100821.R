## loading libraries ##

library("clusterProfiler")
library("org.Dr.eg.db")
library("org.Hs.eg.db")
library("DOSE")
library("dplyr")

data <- read.csv("DiffExp_Q4.00_100821.csv", sep=',')

colnames(data)[1] <- "gene"

ortho <- read.csv("killifish_human_ortho.csv")

## mapping killifish transcripts to human orthologues ## 

data$gene <- plyr::mapvalues(data$gene, from=c(ortho$Gene.stable.ID), to=c(ortho$Gene.stable.ID.1))

data$gene[!startsWith(data$gene , "ENSG0")] <- NA
data <- na.omit(data)

foldchanges <- data$log2FoldChange

names(foldchanges) <- data$gene

foldchanges <- sort(foldchanges, decreasing=TRUE)

head(foldchanges)

# GSEA using gene sets associated with BP Gene Ontology terms
ego <-gseGO(
  foldchanges,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)
gseaGO_results <- ego@result

ego_simple <- simplify(
  ego,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)

ego_simple_summary  <- slot(ego_simple, "result")

write.csv(ego_simple_summary, "gsea_GO_foldchange_Q4.00_100821.csv")
