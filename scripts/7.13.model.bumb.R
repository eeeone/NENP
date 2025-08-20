rm(list = ls())
library(metafor)
library(dplyr)
#setwd("C:/Users/wonby/Desktop/Project/NENP")
df <- read.csv("output/dataframes/effect_size.csv")
df$replicate_ID[is.na(df$replicate_ID)] <- "0"
df$taxon_scientific[is.na(df$taxon_scientific)] <- "NA"
df$land_use_context_ur <- relevel(factor(df$land_use_context_ur), ref = "rural")
df$bd_metric_type <- as.factor(df$bd_metric_type)
df$land_use_context_details <- df$land_use_context_details/100

df$cluster_ID <- interaction(
  df$study_ID,
  df$intervention_treatment_ID,
  df$control_treatment_ID,
  df$replicate_ID,
  df$taxon_scientific,
  df$bd_metric_type,
  drop = TRUE
)                                                      
df <- df %>% filter(tolower(taxon_other) == "bumblebees")
df_abund <- df %>% filter(bd_metric_type == "abundance")
df_rich <- df %>% filter(bd_metric_type == "richness")
table(df_abund$study_ID)
table(df_rich$study_ID)
table(df$study_ID,df$land_use_context_ur)

df_context <- df[!is.na(df$land_use_context_details), ]
df_size <- df[!is.na(df$intervention_size_sqm), ]
df_context_size <- df_context[!is.na(df_context$intervention_size_sqm), ]

##overall
overall <- rma.mv(yi, vi, random = ~1 | study_ID/cluster_ID, data = df, test = "t",method = "REML")
summary(overall)

##moderator_urban_rural
urban_rural <- rma.mv(
  yi, vi,
  mods = ~ land_use_context_ur,
  random = ~1 | study_ID/cluster_ID,
  data = df,
  test = "t",
  method = "REML"
)
summary(urban_rural)

#moderator_semi1000
total_semi1000 <- rma.mv(
  yi, vi,
  mods = ~ land_use_context_details,
  random = ~1 | study_ID/cluster_ID,
  data = df_context,
  test = "t",
  method = "REML"
)
summary(total_semi1000)

##moderator_semi1000_sq
df_context$snh1000_c <- scale(log(df_context$land_use_context_details), center = TRUE, scale = FALSE)
df_context$snh1000_c2 <- df_context$snh1000_c ^2
total_semi1000_sq <- rma.mv(
  yi, vi,
  mods = ~ snh1000_c + snh1000_c2,
  random = ~1 | study_ID/cluster_ID,
  data = df_context,
  test = "t",
  method = "REML"
)
mean(df_context$land_use_context_details)
summary(total_semi1000_sq)

##moderator_size
total_size <- rma.mv(
  yi, vi,
  mods = ~ intervention_size_sqm,
  random = ~1 | study_ID/cluster_ID,
  data = df_size,
  test = "t",
  method = "REML"
)
summary(total_size)

##moderator_semi1000*size
total_snh_size <- rma.mv(
  yi, vi,
  mods = ~ land_use_context_details*log(intervention_size_sqm),
  random = ~1 | study_ID/cluster_ID,
  data = df_context_size,
  test = "t",
  method = "REML"
)
summary(total_snh_size)

