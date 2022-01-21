install.packages("sf")

library(tidyverse)
library(ggmap)

location_df <- read.csv("At_natural_accessions_reference.csv")
location_df <- location_df %>% 
  select(3,5,6)

location_df <- location_df %>%
  drop_na()

locations <- as_tibble(location_df)

register_google(key = "AIzaSyAoS63_VPAXmyFELCj9sfeqEB7S_VkAnNo")

map <- get_googlemap(center = c(31, -7), zoom=6)


#library(sf)
#library(mapview)
