rm(list = ls())
library(metafor)
library(dplyr)
#setwd("C:/Users/wonby/Desktop/Project/NENP")
df <- read.csv("data/data_extraction/data7.csv")
df_baci <- subset(df, study_ID == 516)
df_other  <- subset(df, study_ID != 516)
eps <- 1e-4

#BACI
df_before <- subset(df_baci, replicate_ID == "year_0")
df_after  <- subset(df_baci, replicate_ID != "year_0")
df_baci_combined <- merge(
  df_after,
  df_before,
  by = c("intervention_treatment_ID", "control_treatment_ID",
         "bd_metric_type", "country", "taxon_other"),
  suffixes = c("_after", "_before")
)

df_baci_combined <- df_baci_combined %>%
  mutate(
    i_bd_mean_after  = ifelse(i_bd_mean_after <= 0, eps, i_bd_mean_after),
    i_bd_mean_before = ifelse(i_bd_mean_before <= 0, eps, i_bd_mean_before),
    c_bd_mean_after  = ifelse(c_bd_mean_after <= 0, eps, c_bd_mean_after),
    c_bd_mean_before = ifelse(c_bd_mean_before <= 0, eps, c_bd_mean_before),

    yi = log(i_bd_mean_after) - log(i_bd_mean_before) -
      log(c_bd_mean_after) + log(c_bd_mean_before),
    
    vi = (i_bd_error_value_after  / i_bd_mean_after)^2 +
      (i_bd_error_value_before / i_bd_mean_before)^2 +
      (c_bd_error_value_after  / c_bd_mean_after)^2 +
      (c_bd_error_value_before / c_bd_mean_before)^2
  )

#other
df_other <- df_other %>%
  mutate(
    i_bd_mean = ifelse(i_bd_mean <= 0, eps, i_bd_mean),
    c_bd_mean = ifelse(c_bd_mean <= 0, eps, c_bd_mean),
    
    yi = log(i_bd_mean / c_bd_mean),
    
    vi = ifelse(
      bd_error_metric == "standard error",
      (i_bd_error_value / i_bd_mean)^2 + (c_bd_error_value / c_bd_mean)^2,
      (i_bd_error_value^2 / (i_sample_size * i_bd_mean^2)) +
        (c_bd_error_value^2 / (c_sample_size * c_bd_mean^2))
    )
  )

##combine
df <- df %>%
  left_join(
    df_baci_combined %>%
      select(study_ID = study_ID_after,
             effect_size_ID = effect_size_ID_after,
             yi, vi),
    by = c("study_ID", "effect_size_ID")
  )
df[df$study_ID %in% df_other$study_ID & df$effect_size_ID %in% df_other$effect_size_ID, "yi"] <-
  df_other$yi

df[df$study_ID %in% df_other$study_ID & df$effect_size_ID %in% df_other$effect_size_ID, "vi"] <-
  df_other$vi
df_final <- subset(df, !(study_ID == 516 & replicate_ID == "year_0"))

write.csv(df_final,"output/dataframes/effect_size.csv",row.names = FALSE)
