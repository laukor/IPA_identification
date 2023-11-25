## Assess sites against IPA Criterion A ##
# 1. Load libraries and species list
# 2. Calculate threshold levels per species
# 3. Create country grid cells for analysis
# 4. Calculate occurrence of each species in each grid
# 5. Identify cells meeting Criterion A threshold

# 1. Load libraries and data --------------------------------------------------

library(tidyverse) # for data manipulation
library(sf) # simple features for spatial geometry
library(sp) # for spatial data analysis
library(mapview) # to see mapped outputs
library(doBy) # to work with grouped data
library(raster) # for working with raster data
library(rgeoboundaries) # country political administrative boundaries dataset
library(vegan)
library(SPECIES)
library(rgeos) # 
library(rgdal) # 

# cleaned occurrences for all IPA criterion A triggers
occ <- read.csv("Data/IPACritA_records_thresholds.csv", 
                as.is=TRUE, header=TRUE, encoding = "UTF-8")

# list of IPA criterion A trigger species
critA <- read.csv("Data/IPACritA.csv", 
                  encoding = "UTF-8")

# 2. Calculate Crit A thresholds ----------------------------------------------------

# count global occurrences
# remove entries with no occurrences so that they are not included in count 
occ_present <- occ %>% filter(decimalLongitude != "NA")
# count number of global occurrences of each species
count <- count(occ_present, Accepted_species_name)

# add column of global occurrence count to species list
colnames(critA)[colnames(critA) == "taxon_name"] <- "Accepted_species_name"
critA_glo <- left_join(critA, count, by = "Accepted_species_name")
head(critA_glo)
colnames(critA_glo)[colnames(critA_glo) == "n"] <- "glob_occurrences"
# add 0 to species with no occurrences
critA_glo$glob_occurrences[is.na(critA_glo$glob_occurrences)] <- 0

# count national occurrences
# import shapefile of country boundary
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") # change projection as relevant
country_map  <- readOGR(dsn="Data/shapefiles/",layer="COL_adm0") # change COL for relevant country boundary

# project records to sp object to clip to country
occ_present$decimalLatitude <- as.numeric(occ_present$decimalLatitude)
occ_present$decimalLongitude <- as.numeric(occ_present$decimalLongitude)
coordinates(occ_present) <- c("decimalLongitude", "decimalLatitude")
proj4string(occ_present) <- CRS("+init=epsg:4326") # this is WGS84
is.projected(occ_present)

# form matrix of points which intersect with Colombia (or other country) polygon
in_COL <- gIntersects(occ_present, country_map, byid = T)
# check output
class(in_COL)  # matrix

# return values that don't overlap with country
clipped <- apply(in_COL == F, MARGIN = 2, all) 
# select points that are in country
records_col <- occ_present[which(!clipped), ]  
records_col <- as.data.frame(records_col)
head(records_col)

# check by plotting
records_plot <- st_as_sf(records_col, coords = c("decimalLongitude", "decimalLatitude"),
                         stringsAsFactors=TRUE, crs = 4326) #4326 is the EPSG code for WGS84
records_plot %>% st_geometry() %>% plot 

# count number of national occurrences of each species
count_nat <- count(records_col, Accepted_species_name)

# add count column to species list
critA_counts <- left_join(critA_glo, count_nat, by = "Accepted_species_name")
head(critA_counts)
colnames(critA_counts)[colnames(critA_counts) == "n"] <- "nat_occurrences"
# add 0 to species with no occurences
critA_counts$nat_occurrences[is.na(critA_counts$nat_occurrences)] <- 0

# calculate thresholds 
# (NB this is based on IPA methodology for Colombia, amend values as required)
# calculate and add column for global threshold (1% of occurrences)
critA_counts <- critA_counts %>% mutate(threshold_1per = glob_occurrences * 0.01)
head(critA_glo)

# calculate and add column for global threshold (10% of occurrences)
critA_counts <- critA_counts %>% mutate(threshold_10per = nat_occurrences * 0.1)
head(critA_counts)

# export list of IPA trigger species with thresholds
write.table(critA_counts, file = "Data/IPACritA_records_thresholds.csv",
            sep=",", row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8")

# 3. Create grid cells for country ---------------------------------------------

# load country shapefile (you'll  need to source this for your area of analysis)
col_map  <- st_read(dsn="Data/shapefiles/",layer="COL_adm0")
mapview(col_map) # check mapping
st_crs(col_map) # check projections, etc

# form grid that equals extent of country (5x5km)
cell_size <- c(0.04503, 0.04503) # 0.045 deg ~= 5km
area_square_grid <- st_make_grid(col_map, cellsize = cell_size, what = "polygons", square = TRUE)
# convert to sf and add grid ID
sq_grid_sf <- st_sf(area_square_grid) %>%
  mutate(grid_id = 1:length(lengths(area_square_grid)))

# form grid cells 10x10km
cell_size10 <- c(0.09006, 0.09006)
area_square_grid10 <- st_make_grid(col_map, cellsize = cell_size10, what = "polygons", square = TRUE)
# convert to sf and add grid ID
sq_grid_sf10 <- st_sf(area_square_grid10) %>%
  mutate(grid_id = 1:length(lengths(area_square_grid10)))

# grid currently across extent of maximum lat and long of area of analysis
# clip to just the boundaries
col_map_10 <- st_intersects(sq_grid_sf10, col_map, sparse = FALSE)
class(col_map_10) # matrix of T or F on whether intersects
col_map_10_df <- as.data.frame(col_map_10)

# return values that don't overlap with COL
clipped <- apply(col_map_10 == F, MARGIN = 1, all) 

# select grids that are in COL
grids_col <- sq_grid_sf10[which(!clipped), ]  
grids_col_sf <- st_as_sf(grids_col,
                         stringsAsFactors=TRUE, crs = 4326)

grids_col_sf %>% st_geometry() %>% plot # check output

# export 5x5 grids as shapefile
st_write(sq_grid_sf, "Data/shapefiles", "COL_grid_5x5", driver = "ESRI Shapefile")

# export 10x10 grids clipped to Colombian extent as shapefile
st_write(grids_col_sf, "Data/shapefiles", "COL_grid_10x10", driver = "ESRI Shapefile")


# 4. Species occurrence per grid ---------------------------------------------
# make occurrences into sf object (remove NA values first)
occ <- occ %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)
# check mapping
mapview(occ, cex = 2, alpha = .4, popup = NULL)

# count points
# https://gis.stackexchange.com/questions/323698/counting-points-in-polygons-with-sf-package-of-r
# add column with number of points within each grid
col_10_sf$n_colli <- lengths(st_intersects(col_10_sf, occ))

# remove grid with value of 0 (i.e. no occurrences in that grid)
col_10_present <- filter(col_10_sf, n_colli > 0)
mapview(col_10_present)

# form matrix of which occurrences are within each grid
species_grid_length <- st_intersects(col_10_present, occ, sparse = TRUE)
species_grid <- st_intersects(col_10_present, occ, sparse = FALSE)
species_grid_df <- as.data.frame(species_grid)
# rename columns to reflect species names and rows to grid id
colnames(species_grid_df) <- occ$Accepted_species_name
rownames(species_grid_df) <- col_10_present$grid_id

# change TRUE/FALSE to 1 and 0 to allow summing
species_grid_presence <- species_grid_df*1

# condense columns by summing occurrences of each species for each grid cell
# and transpose so that species are rows and grid cells are columns
species_grid_sum <- as.data.frame(rowsum(t(species_grid_presence), 
                                         group = colnames(species_grid_presence), na.rm = T))
View(species_grid_sum)

# export matrix of species x grids
write.table(species_grid_sum, file = "Data/speciesxgrids.csv",
            sep=",", row.names = TRUE, col.names = TRUE, fileEncoding = "UTF-8")


# 5. UAs meeting threasholds

# Crit A species x grid squares with occurrences
matrix <- species_grid_sum

# list of CUPC of conservation concern
critA <- read.csv("Data/CUPC_IPACritA_221214.csv",
                  fileEncoding = "UTF-8")

# join occurrence matrix to list of Crit A species
# format and join data -------------------------------------------------------------
colnames(critA)
colnames(matrix)[colnames(matrix) == ""] <- "Accepted_species_name"

# reformat species x grid data from wide to long 
# ie. one row per unique species-grid combination
data_long <- matrix %>% 
  pivot_longer(cols = c(2:2119), # change numbers as appropriate
               names_to = "grid_id",
               values_to = "count") %>%
  filter(count > 0)

# form table of all species of conservation concern, 
# with the grids they occur in and number of occurrences per grid
critA_grids <- left_join(critA, data_long,
                         by = "Accepted_species_name")
colnames(critA_grids)

# NB if there are fewer rows than before joining, 
# these species probably did not have occurrence records
# can check which species these are by using
spp_no_grid <- critA_grids %>% filter(is.na(grid))

# add columns to show if row meets criteria
critA_grids_thresh <- critA_grids %>%
  mutate(meets_1per = if_else(count >= threshold_1per, "Y", "N")) %>%
  mutate(meets_10per = if_else(count >= threshold_10per, "Y", "N"))

# filter to only rows which meet either global or national criteria
critA_grids_trigger <- critA_grids_thresh %>%
  subset(grepl("Y", meets_1per) | grepl("Y", meets_10per))

# identify which sub-criteria are triggered
# NB. this is based on Colombia IPA methodology thresholds, change as appropriate
# set as 1 or 0 so that can be summed at later stage to determine # of sub-crit met
critA_IPA_trigger <- critA_grids_trigger %>%
  mutate(A1_trigger = if_else(IPA_A1 == "Y" & (meets_1per == "Y" | meets_10per == "Y"),
                              1, 0)) %>%
  mutate(A3_trigger = if_else(IPA_A3 == "Y" & (meets_10per == "Y"),
                              1, 0)) %>%
  mutate(A4_trigger = if_else(IPA_A4 == "Y" & (meets_10per == "Y"),
                              1, 0)) %>%
  mutate(A5_trigger = if_else(IPA_A5 == "Y" & (meets_10per == "Y"),
                              1, 0)) 

head(critA_IPA_trigger)

# remove NA values to 0
critA_IPA_trigger$A4_trigger[is.na(critA_IPA_trigger$IPA_A4)] <- 0
critA_IPA_trigger$A5_trigger[is.na(critA_IPA_trigger$IPA_A5)] <- 0
critA_IPA_trigger[is.na(critA_IPA_trigger)] <- "NA"

# filter to only species and grids which meet thresholds
critA_IPA_triggers <- critA_IPA_trigger %>% filter(A1_trigger == 1 | A3_trigger == 1 |
                                                     A4_trigger == 1 | A5_trigger == 1)
# to check how many cells meet each sub-criterion,
# can use following code, changing "A1" for the criterion in question
table(critA_IPA_trigger$A1_trigger)

# export table of all Criterion A species with information on grids 
# which meet thresholds that potentially trigger IPAs
write.table(critA_IPA_triggers, file = "Data/IPACritA_grids.csv",
            sep=",", row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8")

