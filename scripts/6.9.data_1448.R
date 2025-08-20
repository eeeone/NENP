rm(list = ls())
library(dplyr)
library(readxl)
library(stringr)
#setwd("C:/Users/wonby/Desktop/Project/NENP")
df <- read_excel("data/data_extraction/1448_Collected_Bees_(raw_data_pan_traps).xlsx", sheet = 1)
names(df)

##bumblebees
species_cols <- colnames(df)[10:ncol(df)]
bumb_cols <- species_cols[grepl("^Bombus\\.", species_cols)]

df_cleaned <- df %>%
  mutate(
    bumb_abundance = rowSums(select(., all_of(bumb_cols)), na.rm = TRUE),
    bumb_richness = rowSums(select(., all_of(bumb_cols)) > 0, na.rm = TRUE)
  )

df_bumb_summary <- df_cleaned  %>%
  select(Year, Season, Neighborhood, Treatment,
         bumb_abundance, bumb_richness)

##bees(main)
df_summary <- df  %>%
  select(Year, Season, Neighborhood, Treatment,
         `Bee Richness`, `Bee Abundances`)
df_summary <- df_summary %>%
  rename(treat = Treatment, year = Year, season = Season,neighborhood=Neighborhood,Bee_Richness = `Bee Richness`,Bee_Abundances = `Bee Abundances`)

df_filtered <- df_summary %>%
  filter(treat %in% c(1, 4, 6, 8))

df_summary1 <- df_filtered %>%
  group_by(year, season, treat) %>%
  summarise(
    mean_abundance = mean(Bee_Abundances, na.rm = TRUE),
    se_abundance = sd(Bee_Abundances, na.rm = TRUE) / sqrt(n()),
    mean_richness = mean(Bee_Richness, na.rm = TRUE),
    se_richness = sd(Bee_Richness, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

group_all_years <- df_summary1 %>%
  filter(treat %in% c(1, 4, 6, 8)) %>%  
  group_by(treat) %>%  
  summarise(
    mean_abundance = mean(mean_abundance, na.rm = TRUE),
    se_abundance = sqrt(sum(se_abundance^2, na.rm = TRUE)) / n(),
    mean_richness = mean(mean_richness, na.rm = TRUE),
    se_richness = sqrt(sum(se_richness^2, na.rm = TRUE)) / n(),
    mean_semi1500 = mean(mean_semi1500, na.rm = TRUE),
    cv_semi1500 = mean(cv_semi1500, na.rm = TRUE),
    n_total = sum(n),
    .groups = "drop"
  )

##bumb_subset
df_bumb_summary <- df_bumb_summary %>%
  rename(treat = Treatment, year = Year, season = Season,neighborhood=Neighborhood)

df_bumb_filtered <- df_bumb_summary %>%
  filter(treat %in% c(1, 4, 6, 8))

df_bumb_summary1 <- df_bumb_filtered %>%
  group_by(year, season, treat) %>%
  summarise(
    mean_bumb_abundance = mean(bumb_abundance, na.rm = TRUE),
    se_bumb_abundance = sd(bumb_abundance, na.rm = TRUE) / sqrt(n()),
    mean_bumb_richness = mean(bumb_richness, na.rm = TRUE),
    se_bumb_richness = sd(bumb_richness, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )
