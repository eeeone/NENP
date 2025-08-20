rm(list = ls())
library(terra)
library(sf)
library(exactextractr)
library(dplyr)
#setwd("C:/Users/wonby/Desktop/Project/NENP")
location <- read.csv("data/data_extraction/rural_coordinates.csv")
location$trt <- as.character(location$trt)
location_eu <- location %>% 
  filter(!(study_ID %in% c(4,1155,1815)))
tif1 <- rast("data/data_extraction/rural_snh_extraction/CLCplus_2018_010m_1.tif")
tif2 <- rast("data/data_extraction/rural_snh_extraction/CLCplus_2018_010m_2.tif")
tif3 <- rast("data/data_extraction/rural_snh_extraction/CLCplus_2018_010m_3.tif")
tif4 <- rast("data/data_extraction/rural_snh_extraction/CLCplus_2018_010m_4.tif")
tif5 <- rast("data/data_extraction/rural_snh_extraction/CLCplus_2018_010m_5.tif")
tif6 <- rast("data/data_extraction/rural_snh_extraction/CLCplus_2018_010m_6.tif")
tif7 <- rast("data/data_extraction/rural_snh_extraction/CLCplus_2018_010m_7.tif")
tif8 <- rast("data/data_extraction/rural_snh_extraction/CLCplus_2018_010m_8.tif")
tif <- merge(tif1,tif2,tif3,tif4,tif5,tif6,tif7,tif8)

pollinator_snh_codes <- c(2, 3, 4, 5, 6)
points_sf <- st_as_sf(location_eu, coords = c("longitude", "latitude"), crs = 4326)
points_proj <- st_transform(points_sf, crs = crs(tif))
#st_crs(points_proj)$units
buffers <- st_buffer(points_proj, dist = 1000)

snh_proportions <- exact_extract(tif, buffers, function(vals, covs) {
  sum(vals %in% pollinator_snh_codes, na.rm = TRUE) / length(vals)
})

location_eu$snh_CLC <- round(100 * snh_proportions, 5)

final_location_with_snhCLC <- location_eu %>%
  group_by(trt) %>%
  summarise(
    snh_CLC_mean = round(mean(snh_CLC, na.rm = TRUE), 5),
    snh_CLC_cv = if (sum(!is.na(snh_CLC)) <= 1) {
      "0%"
    } else {
      paste0(
        round(100 * sd(snh_CLC, na.rm = TRUE) / mean(snh_CLC, na.rm = TRUE), 2), "%"
      )
    },
    .groups = "drop"
  )

