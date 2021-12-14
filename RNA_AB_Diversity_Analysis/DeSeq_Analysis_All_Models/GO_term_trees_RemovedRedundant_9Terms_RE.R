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

# cytokine-mediated signaling pathway #  GO:0019221 #
# defence response to other organism #  GO:0098542 #
# leukocyte proliferation # GO:0070661 #

immune_trees <- c("GO:0002376", "GO:0006954", "GO:0002682", "GO:0050727",
                  "GO:0002764", "GO:0001816", "GO:0019221", "GO:0098542",
                  "GO:0070661")

children <- as.list(GOBPOFFSPRING[immune_trees])

immune_system_processes <- c(children$`GO:0002376`, "GO:0002376")
inflammatory_response <- c(children$`GO:0006954`, "GO:0006954")
regulation_of_immune_system_processes <- c(children$`GO:0002682`, "GO:0002682")
regulation_of_inflammatory_response <- c(children$`GO:0050727`, "GO:00050727")
immune_response_regulating_signaling_pathway <- c(children$`GO:0002764`, "GO:0002764")
cytokine_production <- c(children$`GO:0001816`, "GO:0001816")
cytokine_mediated_signaling_pathway <- c(children$`GO:0019221`, "GO:0019221")
defence_response_to_other_organism <- c(children$`GO:0098542`, "GO:0098542")
leukocyte_proliferation <- c(children$`GO:0070661`, "GO:0070661")

## labeling additional immune terms, not captured as children ##
#cytokine-mediated signaling pathway GO:0019221 #
#defense response to bacterium GO:0042742 #
#leukocyte proliferation GO:0070661 #
#phagocytosis GO:0006909 #
#phagosome acidification GO:0090383 #
#phagosome maturation GO:0090382 #

#unlabeled_immune_terms <- c("GO:0019221", "GO:0042742", "GO:0070661", "GO006909",
                            #"GO:0090383", "GO:0090382")

## Generating new columns ##

data$immune_system_processes <- data$X
data$inflammatory_response <- data$X
data$regulation_of_immune_system_processes <- data$X
data$regulation_of_inflammatory_response <- data$X
data$immune_response_regulating_signaling_pathway <- data$X
data$cytokine_production <- data$X
data$cytokine_mediated_signaling_pathway <- data$X
data$defence_response_to_other_organism <- data$X
data$leukocyte_proliferation <- data$X
#data$unlabeled_immune_terms <- data$X

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

for (ID in data$ID){
  data[ID,]$cytokine_mediated_signaling_pathway <- ID%in% cytokine_mediated_signaling_pathway
}

for (ID in data$ID){
  data[ID,]$defence_response_to_other_organism <- ID%in% defence_response_to_other_organism
}

for (ID in data$ID){
  data[ID,]$leukocyte_proliferation <- ID%in% leukocyte_proliferation
}

#for (ID in data$ID){
#  data[ID,]$unlabeled_immune_terms <- ID%in% unlabeled_immune_terms
#}

write.csv(data, "Results/GSEA_Analysis/gsea_simplified_results_RemovedRedundantGenes_Q0.00_9GOterms_230921.csv")

