library(org.Hs.eg.db)
library(GO.db)

## load enrichment output ## 

data <- read.csv("Results/GSEA_Analysis/gsea_simplified_results_RemovedRedundantGenes_Q0.00_060921.csv")

rownames(data) <- data$ID


# immune system processes = GO:0002376 #
# inflammatory response = # GO:0006954 #
# regulation of immune system processes = # GO:0002682 #
# regulation of inflammatory response = # GO:0050727 #
# immune response-regulating signaling pathway # GO:0002764 #
# cytokine production GO:0001816 #

immune_trees <- c("GO:0002376", "GO:0006954", "GO:0002682", "GO:0050727",
                  "GO:0002764", "GO:0001816")

children <- as.list(GOBPOFFSPRING[immune_trees])

immune_system_processes <- children$`GO:0002376`
inflammatory_response <- children$`GO:0006954`
regulation_of_immune_system_processes <- children$`GO:0002682`
regulation_of_inflammatory_response <- children$`GO:0050727`
immune_response_regulating_signaling_pathway <- children$`GO:0002764`
cytokine_production <- children$`GO:0001816`

## Generating new columns ##

data$immune_system_processes <- data$X
data$inflammatory_response <- data$X
data$regulation_of_immune_system_processes <- data$X
data$regulation_of_inflammatory_response <- data$X
data$immune_response_regulating_signaling_pathway <- data$X
data$cytokine_production <- data$X

for (ID in data$ID){
  data[ID,]$immune_system_processes <- ID%in% immune_system_processes 
}

for (ID in data$ID){
  data[ID,]$inflammatory_response <- ID%in% inflammatory_response  
}

for (ID in data$ID){
  data[ID,]$regulation_of_immune_system_processes <- ID%in% regulation_of_immune_system_processes
}

for (ID in data$ID){
  data[ID,]$regulation_of_inflammatory_response <- ID%in% regulation_of_inflammatory_response
}

for (ID in data$ID){
  data[ID,]$immune_response_regulating_signaling_pathway <- ID%in% immune_response_regulating_signaling_pathway
}

for (ID in data$ID){
  data[ID,]$cytokine_production <- ID%in% cytokine_production
}

write.csv(data, "Results/GSEA_Analysis/gsea_simplified_results_RemovedRedundantGenes_Q0.00_6GOterms_060921.csv")

