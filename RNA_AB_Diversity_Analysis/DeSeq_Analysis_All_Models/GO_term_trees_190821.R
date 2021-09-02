library(org.Hs.eg.db)
library(GO.db)

## load enrichment output ## 

data <- read.csv("Results/RemovingFish/GSEA/gsea_simplified_results_Q1.00_NoFish3_220821.csv")

rownames(data) <- data$ID

# immune system processes = GO:0002376 #
# inflammatory response = # GO:0006954 #
# regulation of immune system processes = # GO:0002682 #
# regulation of inflammatory response = # GO:0050727 # 

immune_trees <- c("GO:0002376", "GO:0006954", "GO:0002682", "GO:0050727")

children <- as.list(GOBPOFFSPRING[immune_trees])

immune_system_processes <- children$`GO:0002376`
inflammatory_response <- children$`GO:0006954`  
regulation_of_immune_system_processes <- children$`GO:0002682`
regulation_of_inflammatory_response <- children$`GO:0050727`


## Generating new columns ##

data$immune_system_processes <- data$X
data$inflammatory_response <- data$X
data$regulation_of_immune_system_processes <- data$X
data$regulation_of_inflammatory_response <- data$X

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

write.csv(data, "Results/RemovingFish/GSEA//sea_simplified_results_Q1.00_WithGOTreeFinal_NoFish4_230821.csv")

