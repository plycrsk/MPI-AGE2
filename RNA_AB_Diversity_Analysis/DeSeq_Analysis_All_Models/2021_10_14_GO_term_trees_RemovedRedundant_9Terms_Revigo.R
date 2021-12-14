library(org.Hs.eg.db)
library(GO.db)

## load enrichment output ## 

data <- read.csv("PhD/MPI-AGE2/RNA_AB_Diversity_Analysis/DeSeq_Analysis_All_Models/Results/Revigo/Revigo_Q4.00_NES_Sorted.csv")

rownames(data) <- data$TermID


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
#phagosome acTermIDification GO:0090383 #
#phagosome maturation GO:0090382 #

#unlabeled_immune_terms <- c("GO:0019221", "GO:0042742", "GO:0070661", "GO006909",
                            #"GO:0090383", "GO:0090382")

## Generating new columns ##

data$immune_system_processes <- data$TermID
data$inflammatory_response <- data$TermID
data$regulation_of_immune_system_processes <- data$TermID
data$regulation_of_inflammatory_response <- data$TermID
data$immune_response_regulating_signaling_pathway <- data$TermID
data$cytokine_production <- data$TermID
data$cytokine_mediated_signaling_pathway <- data$TermID
data$defence_response_to_other_organism <- data$TermID
data$leukocyte_proliferation <- data$TermID
#data$unlabeled_immune_terms <- data$X

for (TermID in data$TermID){
  data[TermID,]$immune_system_processes <- TermID%in% immune_system_processes 
}

for (TermID in data$TermID){
  data[TermID,]$inflammatory_response <- TermID%in% inflammatory_response  
}

for (TermID in data$TermID){
  data[TermID,]$regulation_of_immune_system_processes <- TermID%in% regulation_of_immune_system_processes
}

for (TermID in data$TermID){
  data[TermID,]$regulation_of_inflammatory_response <- TermID%in% regulation_of_inflammatory_response
}

for (TermID in data$TermID){
  data[TermID,]$immune_response_regulating_signaling_pathway <- TermID%in% immune_response_regulating_signaling_pathway
}

for (TermID in data$TermID){
  data[TermID,]$cytokine_production <- TermID%in% cytokine_production
}

for (TermID in data$TermID){
  data[TermID,]$cytokine_mediated_signaling_pathway <- TermID%in% cytokine_mediated_signaling_pathway
}

for (TermID in data$TermID){
  data[TermID,]$defence_response_to_other_organism <- TermID%in% defence_response_to_other_organism
}

for (TermID in data$TermID){
  data[TermID,]$leukocyte_proliferation <- TermID%in% leukocyte_proliferation
}

#for (TermID in data$TermID){
#  data[TermID,]$unlabeled_immune_terms <- TermID%in% unlabeled_immune_terms
#}

write.csv(data, "PhD/MPI-AGE2/RNA_AB_Diversity_Analysis/DeSeq_Analysis_All_Models/Results/Revigo/Revigo_ImmuneTagged_NES_Q4.00.csv")

