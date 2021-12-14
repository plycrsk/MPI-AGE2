# Load libraries
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)
library(stringr)
library(dplyr)

## Note: replace diversity order, e.g. Q2.00, with desired diversity order ##
## throughout this R code                                                  ##

data <- read.csv("Output/DiffExp_age_plus_div_Q2.00.csv", sep=',')

colnames(data)[1] <- "gene"

ortho <- read.csv("killifish_human_ortho.csv")

## mapping killifish transcripts to human orthologues ## 

data$gene <- plyr::mapvalues(data$gene, from=c(ortho$Gene.stable.ID), to=c(ortho$Gene.stable.ID.1))

data$gene[!startsWith(data$gene , "ENSG0")] <- NA
data <- na.omit(data)
old_data <- data

## Summarise those transcripts that have duplicate matches ## 
# Some genes are duplicated, due to many:one mapping. Take mean of these and remove duplicates #

data$FC <- 2^data$log2FoldChange

data <- data %>% group_by(gene) %>% summarise_all(funs(mean))

data <- as.data.frame(data)

data$log2FoldChange <- log2(data$FC)

#data[-1] <- sapply(data[-1], as.integer)

##foldchanges for GSEA analysis##

foldchanges <- data$log2FoldChange

names(foldchanges) <- data$gene

## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)

head(foldchanges)

# GSEA using gene sets associated with BP Gene Ontology terms
gseaGO <-gseGO(
  foldchanges,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  exponent = 1,
  minGSSize = 10,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

gseaGO_results <- gseaGO@result
save(gseaGO, file="Output/GSEA_analysis/ego_GSEA_Simplified_age_plus_div_Q2.00.rda")

gsea_simple <- clusterProfiler::simplify(
  gseaGO,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

gsea_simple_results <- gsea_simple@result
write.csv(gsea_simple_results, "Output/GSEA_analysis/gsea_simplified_results_Q2.00.csv")

## Plotting ##

ego_summary <- read.csv("Output/GSEA_analysis/gsea_simplified_results_Q2.00.csv")

ego_summary <- ego_summary[order(-abs(ego_summary$NES)),]
ego_summary$type[ego_summary$NES < 0] = "Downregulated"
ego_summary$type[ego_summary$NES > 0] = "Upregulated"
#dot_df = ego_summary
dot_df_up <- ego_summary[ego_summary$type == "Upregulated", ][1:15,]
dot_df_down <- ego_summary[ego_summary$type == "Downregulated", ][1:15,]


gene_count_up <- dot_df_up %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/"))+ 1)
gene_count_down <- dot_df_down %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/"))+ 1)
## count the gene number

## merge with the original dataframe
dot_df_up <- left_join(dot_df_up, gene_count_up, by = "ID") %>% mutate(GeneRatio = count/setSize)
dot_df_down <- left_join(dot_df_down, gene_count_down, by = "ID") %>% mutate(GeneRatio = count/setSize)

X11()

## from Tommy's code ## 
p <- ggplot(dot_df_up, aes(x = NES, y = fct_reorder(Description, NES))) + 
  geom_point(aes(size = count, color = p.adjust)) +
  theme_bw(base_size = 11) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("GO pathway enrichment")

p + facet_grid(.~type)

dev.copy2pdf(file = "Output/GSEA_analysis/gseGO_dotplot_Enriched_Q2.00.pdf")

ggsave(
  "Output/GSEA_analysis/gseGO_dotplot_Enriched_Q2.00.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 14,
  height = 6.99,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)
