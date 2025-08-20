rm(list = ls())
library(dplyr)
library(readxl)
#setwd("C:/Users/wonby/Desktop/Project/NENP")
df <- read_excel("data/data_extraction/1028_Bumble_pooled_for_analysis.xlsx")

df_summary <- df %>%
  group_by(Treatment) %>%
  summarise(
    area_mean = 10000*mean(area_ha, na.rm = TRUE),
    area_cv = 100 * sd(area_ha, na.rm = TRUE) / mean(area_ha, na.rm = TRUE),
    
    bumble_richness_mean = mean(Bumble_taxa, na.rm = TRUE),
    bumble_richness_se = sd(Bumble_taxa, na.rm = TRUE) / sqrt(sum(!is.na(Bumble_taxa))),
    
    bumble_abundance_mean = mean(Bumble_abu, na.rm = TRUE),
    bumble_abundance_se = sd(Bumble_abu, na.rm = TRUE) / sqrt(sum(!is.na(Bumble_abu))),
    
    semi1000_mean = mean(`%semi1000`, na.rm = TRUE),
    semi1000_cv = 100 * sd(`%semi1000`, na.rm = TRUE) / mean(`%semi1000`, na.rm = TRUE)
  )
