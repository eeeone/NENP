rm(list = ls())
library(terra)
library(sf)
library(exactextractr)
library(dplyr)

#setwd("C:/Users/wonby/Desktop/Project/NENP")
location <- read.csv("data/data_extraction/rural_coordinates.csv")
location$trt <- as.character(location$trt)
tif1 <- rast("data/data_extraction/rural_snh_extraction/ESA_WorldCover_10m_2020_v100_N36E021_Map.tif")
tif2 <- rast("data/data_extraction/rural_snh_extraction/ESA_WorldCover_10m_2020_v100_N45E009_Map.tif")
tif3 <- rast("data/data_extraction/rural_snh_extraction/ESA_WorldCover_10m_2020_v100_N51E009_Map.tif")
tif4 <- rast("data/data_extraction/rural_snh_extraction/ESA_WorldCover_10m_2020_v100_N51W003_Map.tif")
tif5 <- rast("data/data_extraction/rural_snh_extraction/ESA_WorldCover_10m_2020_v100_N54W003_Map.tif")
tif6 <- rast("data/data_extraction/rural_snh_extraction/ESA_WorldCover_10m_2020_v100_N36W078_Map.tif")
tif7 <- rast("data/data_extraction/rural_snh_extraction/ESA_WorldCover_10m_2020_v100_N36W123_Map.tif")
tif8 <- rast("data/data_extraction/rural_snh_extraction/ESA_WorldCover_10m_2020_v100_N42W087_Map.tif")
tif9 <- rast("data/data_extraction/rural_snh_extraction/ESA_WorldCover_10m_2020_v100_N27W084_Map.tif")

tif516_1 <- rast("data/data_extraction/rural_snh_extraction/ESA_WorldCover_10m_2020_v100_N51E006_Map.tif")
tif516_2 <- rast("data/data_extraction/rural_snh_extraction/ESA_WorldCover_10m_2020_v100_N54E012_Map.tif")
tif516_3 <- rast("data/data_extraction/rural_snh_extraction/ESA_WorldCover_10m_2020_v100_N48E009_Map.tif")
#ext(tif516_3)
pollinator_snh_codes <- c(10, 20, 30, 90)
extract_snh <- function(tile_code, location_df, raster_layer, buffer_radius) {
  
  tile_location <- subset(location_df, tif == tile_code)
  if (nrow(tile_location) == 0) {
    warning(paste0("No points found for tile: ", tile_code))
    return(NULL)
  }
  
  tile_sf <- st_as_sf(tile_location, coords = c("longitude", "latitude"), crs = 4326)
  tile_proj <- st_transform(tile_sf, crs = crs(raster_layer))
  tile_buffers <- st_buffer(tile_proj, dist = buffer_radius)
  
  semi_values <- exact_extract(raster_layer, tile_buffers, function(vals, covs) {
    sum(vals %in% pollinator_snh_codes, na.rm = TRUE) / length(vals)
  })
  
  snh_colname <- paste0("snh", buffer_radius)
  tile_location[[snh_colname]] <- round(100 * semi_values, 5)
  
  return(tile_location[, c("longitude", "latitude", "tif", snh_colname)])
}

summarise_snh_radius <- function(buffer_radius) {
  snh_col <- paste0("snh", buffer_radius)
  
  N36E021_snh <- extract_snh("N36E021", location, tif1, buffer_radius = buffer_radius)
  N45E009_snh <- extract_snh("N45E009", location, tif2, buffer_radius = buffer_radius)
  N51E009_snh <- extract_snh("N51E009", location, tif3, buffer_radius = buffer_radius)
  N51W003_snh <- extract_snh("N51W003", location, tif4, buffer_radius = buffer_radius)
  N54W003_snh <- extract_snh("N54W003", location, tif5, buffer_radius = buffer_radius)
  N36W078_snh <- extract_snh("N36W078", location, tif6, buffer_radius = buffer_radius)
  N36W123_snh <- extract_snh("N36W123", location, tif7, buffer_radius = buffer_radius)
  N42W087_snh <- extract_snh("N42W087", location, tif8, buffer_radius = buffer_radius)
  N27W084_snh <- extract_snh("N27W084", location, tif9, buffer_radius = buffer_radius)
  N51E006_snh <- extract_snh("N51E006", location, tif516_1, buffer_radius = buffer_radius)
  N54E012_snh <- extract_snh("N54E012", location, tif516_2, buffer_radius = buffer_radius)
  N48E009_snh <- extract_snh("N48E009", location, tif516_3, buffer_radius = buffer_radius)
  
  all_snh <- rbind(
    N36E021_snh,
    N45E009_snh,
    N51E009_snh,
    N51W003_snh,
    N54W003_snh,
    N36W078_snh,
    N36W123_snh,
    N42W087_snh,
    N27W084_snh,
    N51E006_snh,
    N54E012_snh,
    N48E009_snh
  )
  
  merged <- all_snh[, c("longitude", "latitude", snh_col)]
  
  location_with_snh <- merge(location, merged,
                             by = c("longitude", "latitude"),
                             all.x = TRUE, sort = FALSE)
  
  final_location_with_snh <- location_with_snh %>%
    group_by(trt) %>%
    summarise(
      !!paste0(snh_col, "_mean") := round(mean(.data[[snh_col]], na.rm = TRUE), 5),
      !!paste0(snh_col, "_cv") := paste0(round(100 * sd(.data[[snh_col]], na.rm = TRUE) / mean(.data[[snh_col]], na.rm = TRUE), 2), "%"),
      .groups = "drop"
    )
  
  return(final_location_with_snh)
}
snh1000_summary <- summarise_snh_radius(1000)
snh500_summary <- summarise_snh_radius(500)
snh250_summary <- summarise_snh_radius(250)
final_location_with_snh <- snh1000_summary %>%
  full_join(snh500_summary, by = "trt") %>%
  full_join(snh250_summary, by = "trt")

