rm(list = ls())
library(dplyr)
library(readxl)
library(tidyr)
#setwd("C:/Users/wonby/Desktop/Project/NENP")
df <- read_excel("data/data_extraction/516_Scheper et al_J Appl Ecol_bee, flower & landscape data.xlsx")

make_meta_data <- function(df, response_var) {
  
  meta_data <- df %>%
    group_by(country, year, treatment) %>%
    summarise(
      n = n(),
      mean = mean(.data[[response_var]], na.rm = TRUE),
      se = sd(.data[[response_var]], na.rm = TRUE) / sqrt(n),
      .groups = 'drop'
    ) %>%
    ungroup() %>%
    pivot_wider(
      names_from = treatment,
      values_from = c(n, mean, se)
    ) %>%
    rename(
      c_bd_mean = mean_field_boundary,
      c_bd_error_value = se_field_boundary,
      i_bd_mean = mean_wildflower_strip,
      i_bd_error_value = se_wildflower_strip
    )
  
  moderator_vars <- df %>%
    filter(treatment == "wildflower_strip") %>%
    group_by(country, year) %>%
    summarise(
      SNH_treat = mean(SNH_perc, na.rm = TRUE),
      SNH_treat_CV = 100*sd(SNH_perc, na.rm = TRUE) / SNH_treat,
      .groups = 'drop'
    )
  
  meta_data_final <- left_join(meta_data, moderator_vars, by = c("country", "year"))
  
  return(meta_data_final)
}

meta_data_bba <- make_meta_data(df, "bumblebee_abundance")
meta_data_bbr <- make_meta_data(df, "bumblebee_richness")
meta_data_solba <- make_meta_data(df, "solitary_bee_abundance")
meta_data_solbr <- make_meta_data(df, "solitary_bee_richness")

