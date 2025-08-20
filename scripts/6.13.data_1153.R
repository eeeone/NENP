rm(list = ls())
library(dplyr)
library(tidyr)
#setwd("C:/Users/wonby/Desktop/Project/NENP")
df <- read.csv("data/data_extraction/1153_Simao_et_al_2018_Data.csv")

df <- df %>%
  mutate(TreatmentGroup = case_when(
    paste(Year, Trt) == "2015 0"  ~ "Control",
    paste(Year, Trt) == "2015 3"  ~ "Treatment_3",
    paste(Year, Trt) == "2015 6"  ~ "Treatment_6",
    paste(Year, Trt) == "2015 10" ~ "Treatment_10",
    paste(Year, Trt) == "2016 1"  ~ "Control",
    paste(Year, Trt) == "2016 3"  ~ "Treatment_3",
    paste(Year, Trt) == "2016 6"  ~ "Treatment_6",
    paste(Year, Trt) == "2016 10" ~ "Treatment_10",
    TRUE ~ NA_character_
  ))

summary_df <- df %>%
  group_by(Year, TreatmentGroup) %>%
  summarise(
    n_samples = 4,
    total_BeeAbun = sum(BeeAbun, na.rm = TRUE),
    mean_BeeAbun = mean(BeeAbun, na.rm = TRUE),
    se_BeeAbun = sd(BeeAbun, na.rm = TRUE) / sqrt(n_samples),
    .groups = "drop"
  )

