## Plots for gseGO analysis at all Q values ## 

library("ggplot2")
library("dplyr")
library("stringr")
library("forcats")


ego_summary <- read.csv("Results/GSEA_Analysis/gsea_simplified_results_RemovedRedundantGenes_Q0.00_6GOterms_060921.csv")
ego_summary <- ego_summary[order(-abs(ego_summary$NES)),]
ego_summary$Q = "0.00"
ego_summary$type[ego_summary$NES < 0] = "Downregulated"
ego_summary$type[ego_summary$NES > 0] = "Upregulated"
#dot_df = ego_summary
dot_df_up <- ego_summary[ego_summary$type == "Upregulated", ][1:10,]
dot_df_down <- ego_summary[ego_summary$type == "Downregulated", ][1:10,]


ego_summary2 <- read.csv("Results/GSEA_Analysis/gsea_simplified_results_RemovedRedundantGenes_Q1.00_6GOterms_060921.csv")
ego_summary2 <- ego_summary2[order(-abs(ego_summary2$NES)),]
ego_summary2$Q = "1.00"
ego_summary2$type[ego_summary2$NES < 0] = "Downregulated"
ego_summary2$type[ego_summary2$NES > 0] = "Upregulated"
#dot_df2 = ego_summary2
dot_df2_up <- ego_summary2[ego_summary2$type == "Upregulated", ][1:10,]
dot_df2_down <- ego_summary2[ego_summary2$type == "Downregulated", ][1:10,]

ego_summary3 <- read.csv("Results/GSEA_Analysis/gsea_simplified_results_RemovedRedundantGenes_Q1.50_6GOterms_060921.csv")
ego_summary3 <- ego_summary3[order(-abs(ego_summary3$NES)),]
ego_summary3$Q = "1.50"
ego_summary3$type[ego_summary3$NES < 0] = "Downregulated"
ego_summary3$type[ego_summary3$NES > 0] = "Upregulated"
#dot_df3 = ego_summary3
dot_df3_up <- ego_summary3[ego_summary3$type == "Upregulated", ][1:10,]
dot_df3_down <- ego_summary3[ego_summary3$type == "Downregulated", ][1:10,]

ego_summary4 <- read.csv("Results/GSEA_Analysis/gsea_simplified_results_RemovedRedundantGenes_Q2.00_6GOterms_060921.csv")
ego_summary4 <- ego_summary4[order(-abs(ego_summary4$NES)),]
ego_summary4$Q = "2.00"
ego_summary4$type[ego_summary4$NES < 0] = "Downregulated"
ego_summary4$type[ego_summary4$NES > 0] = "Upregulated"
#dot_df4 = ego_summary4
dot_df4_up <- ego_summary4[ego_summary4$type == "Upregulated", ][1:10,]
dot_df4_down <- ego_summary4[ego_summary4$type == "Downregulated", ][1:10,]

ego_summary5 <- read.csv("Results/GSEA_Analysis/gsea_simplified_results_RemovedRedundantGenes_Q3.00_6GOterms_060921.csv")
ego_summary5 <- ego_summary5[order(-abs(ego_summary5$NES)),]
ego_summary5$Q = "3.00"
ego_summary5$type[ego_summary5$NES < 0] = "Downregulated"
ego_summary5$type[ego_summary5$NES > 0] = "Upregulated"
#dot_df5 = ego_summary5
dot_df5_up <- ego_summary5[ego_summary5$type == "Upregulated", ][1:10,]
dot_df5_down <- ego_summary5[ego_summary5$type == "Downregulated", ][1:10,]

ego_summary6 <- read.csv("Results/GSEA_Analysis/gsea_simplified_results_RemovedRedundantGenes_Q4.00_6GOterms_060921.csv")
ego_summary6 <- ego_summary6[order(-abs(ego_summary6$NES)),]
ego_summary6$Q = "4.00"
ego_summary6$type[ego_summary6$NES < 0] = "Downregulated"
ego_summary6$type[ego_summary6$NES > 0] = "Upregulated"
#dot_df6 = ego_summary6
dot_df6_up <- ego_summary6[ego_summary6$type == "Upregulated", ][1:10,]
dot_df6_down <- ego_summary6[ego_summary6$type == "Downregulated", ][1:10,]

dot_df_up <- rbind(dot_df_up, dot_df2_up, dot_df3_up, dot_df4_up, dot_df5_up, dot_df6_up)
dot_df_down <- rbind(dot_df_down, dot_df2_down, dot_df3_down, dot_df4_down, dot_df5_down, dot_df6_down)

dot_df_up$ID <- paste(dot_df_up$ID,dot_df_up$Q)
dot_df_down$ID <- paste(dot_df_down$ID,dot_df_down$Q)

gene_count <- dot_df_up %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/"))+ 1)
dot_df_up <- left_join(dot_df_up, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)

gene_count <- dot_df_down %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/"))+ 1)
dot_df_down <- left_join(dot_df_down, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)

#dot_df <- dot_df[order(dot_df['Description']),]
#dot_df <- dot_df[order(dot_df['Description']),]

## adding immune term column if TRUE for any go category ##

dot_df_up <- dot_df_up %>% mutate(Immune = case_when(immune_system_processes == TRUE | cytokine_production == TRUE
                                        | immune_response_regulation == TRUE ~ "Immune"))
dot_df_up[is.na(dot_df_up)] <- "Non-Immune"

dot_df_down <- dot_df_down %>% mutate(Immune = case_when(immune_system_processes == TRUE | cytokine_production == TRUE
                                                         | immune_response_regulation == TRUE ~ "Immune"))
dot_df_down[is.na(dot_df_down)] <- "Non-Immune"

## Toll-like 9 receptor signalling pathway does not pull correctly as immune from GO, ##
## so adding manually ##

dot_df_up[dot_df_up$Description == "toll-like receptor 9 signaling pathway",]$Immune <- "Immune"
 
X11()
p <- ggplot(dot_df_up, aes(x = NES, y = fct_reorder(Description, NES))) + 
  geom_point(aes(color=Immune), size=5) +
  theme(legend.position="none") +
  theme_bw(base_size = 11) +
  ylab(NULL)

p + facet_grid(.~Q)

dev.copy2pdf(file = "gseGO_dotplot_QALL_coloured_top10_RemovedRedundantGenes_upregulated_Simplified_080921.pdf")

ggsave(
  "gseGO_dotplot_QALL_coloured_top10_upregulated_RemovedRedundantGenes_Simplified_080921.png",
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


## Downregulated Plot ## 

X11()
p <- ggplot(dot_df_down, aes(x = NES, y = fct_reorder(Description, NES))) + 
  geom_point(aes(color=Immune), size=5) +
  theme(legend.position="none") +
  theme_bw(base_size = 11) +
  ylab(NULL)

p + facet_grid(.~Q)



dev.copy2pdf(file = "gseGO_dotplot_QALL_coloured_top10_downregulated_RemovedRedundantGenes_Simplified_080921.pdf")

ggsave(
  "gseGO_dotplot_QALL_coloured_top10_downregulated_RemovedRedundantGenes_Simplified_080921.png",
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
