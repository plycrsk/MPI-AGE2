library(org.Hs.eg.db)
library(GO.db)

## load enrichment output ## 

data <- read.csv("Results/GSEA_Analysis/gsea_simplified_results_Q0.00_160821.csv")

rownames(data) <- data$ID

# immune system processes = GO:0002376 #
# immune response =  GO:0006955 #
# inflammatory response = # GO:0006954 #
#regulation_of_immune_system = # GO:0002682 #

immune_trees <- c("GO:0006955", "GO:0006954", "GO:0002376", 'GO:0002682')

children <- as.list(GOBPOFFSPRING[immune_trees])

inflammatory_response <- children$`GO:0006954`  
immune_response <- children$`GO:0006955`
immune_system_processes <- children$`GO:0002376`
regulation_of_immune_system <- children$`GO:0002682`

## Generating new columns ##

data$Immune_System_Processes <- data$X
data$Inflammatory_Response <- data$X
data$Immune_Response <- data$X
data$Regulation <- data$X

for (ID in data$ID){
  data[ID,]$Immune_System_Processes <- ID%in% immune_system_processes 
}

for (ID in data$ID){
  data[ID,]$Inflammatory_Response <- ID%in% inflammatory_response  
}

for (ID in data$ID){
  data[ID,]$Immune_Response <- ID%in% immune_response
}

for (ID in data$ID){
  data[ID,]$Regulation <- ID%in% regulation_of_immune_system
}

write.csv(data, "Results/GoEnrichment_Analysis/test.csv")
