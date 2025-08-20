rm(list = ls())
library(dplyr)
library(readxl)
#setwd("C:/Users/wonby/Desktop/Project/NENP")
df1 <- read_excel("data/data_extraction/1815_bee_community_upload.xlsx", sheet = 4)

trt_size <- df1 %>%
  group_by(Year) %>%
  summarise(
    n = n(),
    mean_size = mean(`Size (m2)`, na.rm = TRUE),
    cv_size = 100 * sd(`Size (m2)`, na.rm = TRUE) / mean(`Size (m2)`, na.rm = TRUE)
  )

df2 <- read_excel("1815_data.xlsx", sheet = 2)
df2_summary <- df2 %>%
  mutate(
    abundance = rowSums(select(., 5:ncol(.)), na.rm = TRUE),
    richness = rowSums(select(., 5:ncol(.)) > 0, na.rm = TRUE)
  ) %>%
  select(1:4, abundance, richness)  

df2_final <- df2_summary %>%
  group_by(Year,Treatment) %>%
  summarise(
    n = 10,
    mean_abundance = mean(abundance, na.rm = TRUE),
    se_abundance = sd(abundance, na.rm = TRUE) / sqrt(n),
    mean_richness = mean(richness, na.rm = TRUE),
    se_richness = sd(richness, na.rm = TRUE) / sqrt(n)
  )

