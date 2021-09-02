## Plots for gseGO analysis at all Q values ## 

library("ggplot2")
library("dplyr")
library("stringr")
library("forcats")

ego_summary <- read.csv("Results/GoEnrichment_Analysis/GoEnrich_Simplified_age_plus_div_Q0.00_150821.csv")
dot_df = ego_summary ## small dataset
dot_df$Q = "0.00"
#dot_df$type[dot_df$NES < 0] = "downregulated"

ego_summary2 <- read.csv("Results/GoEnrichment_Analysis/GoEnrich_Simplified_age_plus_div_Q1.00_150821.csv")
dot_df2 = ego_summary2 ## small dataset
dot_df2$Q = "1.00"

ego_summary3 <- read.csv("Results/GoEnrichment_Analysis/GoEnrich_Simplified_age_plus_div_Q1.50_150821.csv")
dot_df3 = ego_summary3 ## small dataset
dot_df3$Q = "1.50"

ego_summary4 <- read.csv("Results/GoEnrichment_Analysis/GoEnrich_Simplified_age_plus_div_Q2.00_150821.csv")
dot_df4 = ego_summary4 ## small dataset
dot_df4$Q = "2.00"

ego_summary5 <- read.csv("Results/GoEnrichment_Analysis/GoEnrich_Simplified_age_plus_div_Q3.00_150821.csv")
dot_df5 = ego_summary5 ## small dataset
dot_df5$Q = "3.00"

ego_summary6 <- read.csv("Results/GoEnrichment_Analysis/GoEnrich_Simplified_age_plus_div_Q4.00_150821.csv")
dot_df6 = ego_summary6 ## small dataset
dot_df6$Q = "4.00"

dot_df <- rbind(dot_df, dot_df2, dot_df3, dot_df4, dot_df5, dot_df6)

dot_df$core_enrichment <- dot_df$Count
dot_df$setSize <- dot_df$Count

dot_df$ID <- paste(dot_df$ID,dot_df$Q)

gene_count <- dot_df %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/"))+ 1)
dot_df <- left_join(dot_df, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)

dot_df <- dot_df[order(dot_df['Description']),]

dot_df

immune_list <- c(3,4,5,6,18,19,20)

immune_terms <- dot_df[c(immune_list),]
other_terms <-dot_df[-c(immune_list),]

X11()

## from Tommy's code
p <- ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  #geom_point(aes(size = count, color = p.adjust)) +
  geom_point(aes(color="red", size = 10)) +
  theme_bw(base_size = 11) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  geom_point(data = immune_terms, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio)),
             color="yellow", size = 3) +
  ylab(NULL) +
  ggtitle("GO pathway enrichment")

p + facet_grid(.~Q)

dev.copy2pdf(file = "gseGO_dotplot_QALL_coloured_Simplified_220721.pdf")

ggsave(
  "gseGO_dotplot_QALL_coloured_Simplified_220721.png",
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
