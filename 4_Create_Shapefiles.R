## Create shapefiles of UAs which meet IPA thresholds ##
# 1. Load libraries and import data
# 2. Map Criterion A thresholds
# 3. Map Criterion B thresholds

# 1. load libraries ------------------------------------------------

library(tidyverse)
library(sf)
library(dplyr)
library(mapview)
library(tmap) # open-source R-library for drawing thematic maps

# import data

# list of qualifying grids and species triggered of conservation concern
critA <- read.csv("Data/IPACritA_grids.csv",
                  fileEncoding = "UTF-8")
colnames(critA)

# list of qualifying grids for species richness
critB <- read.csv("Data/IPACritB_grids.csv",
                  fileEncoding = "UTF-8")

# Colombia 10x10km grid squares shapefile
# NB. this is the shapefile of grids you produced in #3 of "2_Identify_CritA_Sites"
col_10 <- st_read(dsn="Data/shapefiles/",layer="COL_10_bioregion")
st_crs(col_10) = 4326

# 2. combine data for Crit A ------------------------------------------------------------

# summarise the number of trigger species and criteria per grid cell
colnames(critA)
critA_trigger_grids <- critA[, c(1, 53, 57:60)] # keep only relevant columns about which criteria are triggered
head(critA_trigger_grids)

# number of unique species per grid_id to add back after grouping by grids
spp_per_gridA <- critA_trigger_grids %>% 
  group_by(grid_id) %>% # one row per grid
  summarise(spp_number = n_distinct(Accepted_species_name), # number of unique trigger spp per grid
            # columns with number of unique species triggering each sub-criterion per grid
            A1_trigger = sum(A1_trigger), 
            A3_trigger = sum(A3_trigger),
            A4_trigger = sum(A4_trigger),
            A5_trigger = sum(A5_trigger)) %>%
        as.data.frame() %>% # convert to df
mutate(total_sub_criteria = rowSums(.[3:6] != 0)) # number of unique sub-criteria triggered
head(spp_per_gridA)

# add data to spatial grids
trigger_gridsA <- col_10 %>% #need to join to COL to keep sf object
  left_join(spp_per_gridA, by = "grid_id") %>%
  filter(!is.na(spp_number)) # filter to grids with species data

colnames(trigger_gridsA)

# explore data by ploting grids, colours dependent on number of species triggered
mapview(trigger_gridsA, zcol = "spp_number")

# export shapefile of grids meeting IPA criterion A
st_write(trigger_gridsA, "Data/shapefiles", "CritA_trigger_grids", driver = "ESRI Shapefile")

# 3. combine data for Crit B ------------------------------------------------------------

# keep only relevant columns of critB
colnames(critB)
critB_trigger_grids <- critB[, c(1:3, 12:13)]
head(critB_trigger_grids)

# add info to spatial grids
trigger_gridsB <- col_10 %>% #need to join to COL to keep sf object
  left_join(critB_trigger_grids, by = "grid_id") %>%
  filter(!is.na(S.obs)) # filter to grids which are B triggers

head(trigger_gridsB)

# explore projections by ploting grids, colours dependent on number of species triggered
mapview(trigger_gridsB, zcol = "bio_reg")

# export shapefile of grids meeting IPA criterion B
st_write(trigger_gridsB, "Data/shapefiles", "CritB_trigger_grids", driver = "ESRI Shapefile")

