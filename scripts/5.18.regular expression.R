library(readr)   
library(dplyr)   
library(stringr)
rm(list = ls())
#setwd("C:/Users/wonby/Desktop/Project/NENP")
df <- read.csv("data/data_extraction/relevant_titles_flowers_abstract_completed.csv")
raw_data_simplified <- raw_data %>% select(study_ID,title, abstract, doi,author,year,journal)
raw_data_simplified <- raw_data_simplified %>%
  mutate(full_text = paste(title, abstract, sep = " "))

pollinator_keywords <- c("pollinator","pollinators", "pollination","bee","bees","bumblebee","bumblebees","honeybee","honeybees")
flower_keywords <- c("flower", "wildflower", "wildflowers","planting","floral")
outcome_keywords <- c("abundance", "diversity", "richness")


pollinator_pattern <- paste0("\\b(", paste(pollinator_keywords, collapse = "|"), ")\\b")
flower_pattern     <- paste0("\\b(", paste(flower_keywords, collapse = "|"), ")\\b")
outcome_pattern    <- paste0("\\b(", paste(outcome_keywords, collapse = "|"), ")\\b")

raw_data_filtered <- raw_data_simplified %>%
  mutate(
    pollinator_match = str_detect(full_text, regex(pollinator_pattern, ignore_case = TRUE)),
    flower_match     = str_detect(full_text, regex(flower_pattern, ignore_case = TRUE)),
    outcome_match    = str_detect(full_text, regex(outcome_pattern, ignore_case = TRUE))
  )
raw_data_filtered <- raw_data_filtered  %>%
  filter(pollinator_match & flower_match & outcome_match)
write.csv(raw_data_filtered, "data/data_extraction/relevant_titles_filtered.csv", row.names = FALSE)


library(revtools)
screen_abstracts(raw_data_filtered)



