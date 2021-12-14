## Plots for gseGO analysis at all Q values ## 

library("ggplot2")
library("dplyr")
library("stringr")
library("forcats")


ego_summary <- read.csv("PhD/MPI-AGE2/RNA_AB_Diversity_Analysis/DeSeq_Analysis_All_Models/Results/Revigo/Revigo_ImmuneTagged_NES_Q0.00.csv")
ego_summary$Q = "0.00"
ego_summary$type[ego_summary$Value < 0] = "Downregulated"
ego_summary$type[ego_summary$Value > 0] = "Upregulated"
ego_summary <- ego_summary[ego_summary$Eliminated == " False",]
dot_df_up <- ego_summary[ego_summary$type == "Upregulated", ][1:10,]
dot_df_down <- ego_summary[ego_summary$type == "Downregulated", ][1:10,]


ego_summary2 <- read.csv("PhD/MPI-AGE2/RNA_AB_Diversity_Analysis/DeSeq_Analysis_All_Models/Results/Revigo/Revigo_ImmuneTagged_NES_Q1.00.csv")
ego_summary2$Q = "1.00"
ego_summary2$type[ego_summary2$Value < 0] = "Downregulated"
ego_summary2$type[ego_summary2$Value > 0] = "Upregulated"
ego_summary2 <- ego_summary2[ego_summary2$Eliminated == " False",]
dot_df2_up <- ego_summary2[ego_summary2$type == "Upregulated", ][1:10,]
dot_df2_down <- ego_summary2[ego_summary2$type == "Downregulated", ][1:10,]

ego_summary3 <- read.csv("PhD/MPI-AGE2/RNA_AB_Diversity_Analysis/DeSeq_Analysis_All_Models/Results/Revigo/Revigo_ImmuneTagged_NES_Q1.50.csv")
ego_summary3$Q = "1.50"
ego_summary3$type[ego_summary3$Value < 0] = "Downregulated"
ego_summary3$type[ego_summary3$Value > 0] = "Upregulated"
ego_summary3 <- ego_summary3[ego_summary3$Eliminated == " False",]
dot_df3_up <- ego_summary3[ego_summary3$type == "Upregulated", ][1:10,]
dot_df3_down <- ego_summary3[ego_summary3$type == "Downregulated", ][1:10,]

ego_summary4 <- read.csv("PhD/MPI-AGE2/RNA_AB_Diversity_Analysis/DeSeq_Analysis_All_Models/Results/Revigo/Revigo_ImmuneTagged_NES_Q2.00.csv")
ego_summary4$Q = "2.00"
ego_summary4$type[ego_summary4$Value < 0] = "Downregulated"
ego_summary4$type[ego_summary4$Value > 0] = "Upregulated"
ego_summary4 <- ego_summary[ego_summary4$Eliminated == " False",]
dot_df4_up <- ego_summary4[ego_summary4$type == "Upregulated", ][1:10,]
dot_df4_down <- ego_summary4[ego_summary4$type == "Downregulated", ][1:10,]

ego_summary5 <- read.csv("PhD/MPI-AGE2/RNA_AB_Diversity_Analysis/DeSeq_Analysis_All_Models/Results/Revigo/Revigo_ImmuneTagged_NES_Q3.00.csv")
ego_summary5$Q = "3.00"
ego_summary5$type[ego_summary5$Value < 0] = "Downregulated"
ego_summary5$type[ego_summary5$Value > 0] = "Upregulated"
ego_summary5 <- ego_summary5[ego_summary5$Eliminated == " False",]
dot_df5_up <- ego_summary5[ego_summary5$type == "Upregulated", ][1:10,]
dot_df5_down <- ego_summary5[ego_summary5$type == "Downregulated", ][1:10,]

ego_summary6 <- read.csv("PhD/MPI-AGE2/RNA_AB_Diversity_Analysis/DeSeq_Analysis_All_Models/Results/Revigo/Revigo_ImmuneTagged_NES_Q4.00.csv")
ego_summary6$Q = "4.00"
ego_summary6$type[ego_summary6$Value < 0] = "Downregulated"
ego_summary6$type[ego_summary6$Value > 0] = "Upregulated"
ego_summary6 <- ego_summary6[ego_summary6$Eliminated == " False",]
dot_df6_up <- ego_summary6[ego_summary6$type == "Upregulated", ][1:10,]
dot_df6_down <- ego_summary6[ego_summary6$type == "Downregulated", ][1:10,]

dot_df_up <- rbind(dot_df_up, dot_df2_up, dot_df3_up, dot_df4_up, dot_df5_up, dot_df6_up)
dot_df_down <- rbind(dot_df_down, dot_df2_down, dot_df3_down, dot_df4_down, dot_df5_down, dot_df6_down)

dot_df_up$ID <- paste(dot_df_up$ID,dot_df_up$Q)
dot_df_down$ID <- paste(dot_df_down$ID,dot_df_down$Q)

#gene_count <- dot_df_up %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/"))+ 1)
#gene_count <- dot_df_up$Frequency
#dot_df_up <- left_join(dot_df_up, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)

#gene_count <- dot_df_down %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/"))+ 1)
#dot_df_down <- left_join(dot_df_down, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)

#dot_df <- dot_df[order(dot_df['Name']),]
#dot_df <- dot_df[order(dot_df['Name']),]

## adding immune term column if TRUE for any go category ##

dot_df_up <- dot_df_up %>% mutate(Immune = case_when(immune_system_processes == TRUE | inflammatory_response == TRUE
                                        | regulation_of_immune_system_processes == TRUE |regulation_of_inflammatory_response 
                                        | immune_response_regulating_signaling_pathway | cytokine_production 
                                        | cytokine_mediated_signaling_pathway | defence_response_to_other_organism
                                        | leukocyte_proliferation ~ "Immune"))
dot_df_up[is.na(dot_df_up)] <- "Non-Immune"

dot_df_down <- dot_df_down %>% mutate(Immune = case_when(immune_system_processes == TRUE | inflammatory_response == TRUE
                                                          | regulation_of_immune_system_processes == TRUE |regulation_of_inflammatory_response 
                                                          | immune_response_regulating_signaling_pathway | cytokine_production 
                                                          | cytokine_mediated_signaling_pathway | defence_response_to_other_organism
                                                          | leukocyte_proliferation ~ "Immune"))
dot_df_down[is.na(dot_df_down)] <- "Non-Immune"

## Toll-like 9 receptor signalling pathway does not pull correctly as immune from GO, ##
## so adding manually ##

dot_df_up[dot_df_up$Name == "toll-like receptor 9 signaling pathway",]$Immune <- "Immune"
 
X11()
p <- ggplot(dot_df_up, aes(x = Value, y = fct_reorder(Name, Value))) + 
  geom_point(aes(color=Immune), size=5) +
  theme(legend.position="none") +
  theme_bw(base_size = 11) +
  ylab(NULL)

p + facet_grid(.~Q)

dev.copy2pdf(file = "gseGO_dotplot_QALL_coloured_top10_RemovedRedundantGeValue_upregulated_Revigo_080921.pdf")

ggsave(
  "gseGO_dotplot_QALL_coloured_top10_upregulated_RemovedRedundantGeValue_Revigo_080921.png",
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
p <- ggplot(dot_df_down, aes(x = Value, y = fct_reorder(Name, Value))) + 
  geom_point(aes(color=Immune), size=5) +
  theme(legend.position="none") +
  theme_bw(base_size = 11) +
  ylab(NULL)

p + facet_grid(.~Q)



dev.copy2pdf(file = "gseGO_dotplot_QALL_coloured_top10_downregulated_RemovedRedundantGeValue_Revigo_080921.pdf")

ggsave(
  "gseGO_dotplot_QALL_coloured_top10_downregulated_RemovedRedundantGeValue_Revigo_080921.png",
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
