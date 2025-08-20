rm(list = ls())
library(dplyr)
#setwd("C:/Users/wonby/Desktop/Project/NENP")
df <- read.csv("data/data_extraction/2613_Pollinator dataset.csv")
bees_df <- subset(df, Broad_group == "Bees")
bumblebees_df <- subset(df, Broad_group == "Bees" & Functional_morphotypes == "Bumblebees")

##derive coordinates
df$Latitude <- as.numeric(df$Latitude)
df$Longitude <- as.numeric(df$Longitude)
summary_df <- bees_df %>%
  group_by(Site) %>%
  summarise(
    BeeAbundance = n(), 
    Site_type = first(na.omit(Site_type)),  
    Latitude = median(Latitude, na.rm = TRUE),
    Longitude = median(Longitude, na.rm = TRUE),
    .groups = "drop"
  )
#The coordinate for MS2 is missing from summary_df, so I used one of the available entries for MS2 in bees_df.
lat_value <- 50.4095092
lon_value <- -4.2228162
summary_df <- summary_df %>%
  mutate(
    Latitude = ifelse(Site == "MS2", lat_value, Latitude),
    Longitude = ifelse(Site == "MS2", lon_value, Longitude)
  )

##bees(main)
df_final <- summary_df %>%
  group_by(Site_type) %>%
  summarise(
    n_sites = 25,
    mean_abundance = sum(BeeAbundance, na.rm = TRUE) / 25,
    se_abundance = sd(BeeAbundance, na.rm = TRUE) / sqrt(25),
    .groups = "drop"
  )

##bumblebees
summary_bumb_df <- bumblebees_df %>%
  group_by(Site) %>%
  summarise(
    BumblebeeAbundance = n(), 
    Site_type = first(na.omit(Site_type)),  
    .groups = "drop"
  )

bumb_df_final <- summary_bumb_df %>%
  group_by(Site_type) %>%
  summarise(
    n_sites = 25,
    mean_bumbabundance = sum(BumblebeeAbundance, na.rm = TRUE) / 25,
    se_bumbabundance = sd(BumblebeeAbundance, na.rm = TRUE) / sqrt(25),
    .groups = "drop"
  )
