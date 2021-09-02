library(org.Hs.eg.db)
library(GO.db)

## load enrichment output ## 

data <- read.csv("Results/GSEA_Analysis/gsea_simplified_results_Q4.00_160821.csv")

rownames(data) <- data$ID

as.list(GOBPANCESTOR["GO:0034162"])

# immune system processes = GO:0002376 #
# cytokine production GO:0001816 #
# immune response-regulating signaling pathway # GO:0002764 #

immune_trees <- c("GO:0002376", "GO:0001816", "GO:0002764")

children <- as.list(GOBPOFFSPRING[immune_trees])

immune_system_processes <- children$`GO:0002376`
cytokine_production <- children$`GO:0001816`  
immune_response_regulation <- children$`GO:0002764`

## Generating new columns ##

data$immune_system_processes <- data$X
data$cytokine_production <- data$X
data$immune_response_regulation <- data$X

for (ID in data$ID){
  data[ID,]$immune_system_processes <- ID%in% immune_system_processes 
}

for (ID in data$ID){
  data[ID,]$cytokine_production <- ID%in% cytokine_production  
}

for (ID in data$ID){
  data[ID,]$immune_response_regulation <- ID%in% immune_response_regulation
}

write.csv(data, "Results/GSEA_Analysis/sea_simplified_results_Q4.00_WithGOTree_3terms_240821.csv")

