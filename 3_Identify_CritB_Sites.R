## Assess sites against IPA Criterion B ##
# 1. Load libraries and species list
# 2. Create grid cells per ecosystem
# 3. Generate species list per bioregion
# 4. Estimate species richness
# 5. Identify qualifying UAs for Criterion B

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

# load occurrence points
occ <- st_read(dsn="Data/shapefiles/",layer="all_records")
# check mapping
mapview(occ, cex = 2, alpha = .4, popup = NULL)

# 4. Create grid cells per bioregion --------------------------------------

# load bioregions map (or ecosystems / habitat map, as relevant for area of analysis)
bioregion <- st_read(dsn="Data/shapefiles/",
                     layer="BiomasOlson2001_RegionesColombia_v2")
mapview(bioregion)

# load 10x10km grid squares shapefile
col_10  <- st_read(dsn="Data/shapefiles/",layer="COL_grid_10x10")
mapview(col_10) # check mapping
st_crs(col_10) # check projections, etc

# intersect bioregions and grids  --------------------------------------------------------
bioreg_grid <- st_intersection(col_10, bioregion) %>% 
  bioreg_grid %>% mutate(intersect_area = st_area(.)) %>%   # create new column with shape area
  dplyr::select(grid_id, bio_reg, intersect_area) %>%   # only select columns needed to merge
  
# export 10x10 grids with bioregion data as shapefile
st_write(bioreg_grid, "Data/shapefiles", "COL_10_bioregion", driver = "ESRI Shapefile")

# 3. species list for each bioregion ------------------------------------------------

# with reference to https://gis.stackexchange.com/questions/323698/counting-points-in-polygons-with-sf-package-of-r

# matrix of which occurrences are within each grid
# load bioregion grids as shapefile
col_10_bioreg  <- st_read(dsn="Data/shapefiles/",layer="COL_10_bioregion")

# intersect occurrences and grids
grid_species_length <- st_intersects(occ, col_10_bioreg, sparse = TRUE)
grid_species_length_df <- as.data.frame(grid_species_length)

# rename columns to reflect species names and rows to grid id
grid_species_length_df <- grid_species_length_df %>% 
  mutate(species = occ$Accpt__)
col_10_bioreg$col.id<-1:nrow(col_10_bioreg)
grid_species_length_df <- left_join(grid_species_length_df, col_10_bioreg, by = "col.id")

# keep relevant columns
# species, grid ID, bioregion
colnames(grid_species_length_df)
grid_species_length_df <- grid_species_length_df[,c(3, 4, 8)]

# split dataframe into bioregions
split_bioreg <- split(grid_species_length_df, grid_species_length_df$bio_reg)

# get a list of all species in each bioregion
spp_all_bioreg <- lapply(split_bioreg, 
                         function(x) unique(x$species))

# export species list for each bioregion as separate .csv file
mapply(
  write.table,
  x=spp_all_bioreg, file=paste(names(spp_all_bioreg), "csv", sep="."),
  MoreArgs=list(row.names=FALSE, sep=",")
)

# calculate number of occurrences and species in each bioregion
# dataframe with number of occurrences
occ_bioregion <- table(grid_species_length_df$bio_reg) %>% 
  as.data.frame() %>% occ_bioregion %>% 
  rename(bioregion = Var1, occurences = Freq)

# dataframe with number of species
sppnum_bioregion <- summary(spp_all_bioreg) %>% 
  as.data.frame() 
View(sppnum_bioregion)

sppnum_bioregion <- sppnum_bioregion %>% slice(1:13) %>% 
  rename(bioregion = Var1, species = Freq) %>%
  select(bioregion, species)

# calculate 10% of species for IPA criterion B
sppnum_bioregion$species <- as.numeric(sppnum_bioregion$species)
sppnum_bioregion <- sppnum_bioregion %>% 
  mutate(spe_10per = (species/100*10))

# join dataframes to include occ and species numbers
sppocc_bioregion <- left_join(occ_bioregion, sppnum_bioregion, by = "bioregion")
View(sppocc_bioregion)

# export csv showing number of occurrences and species in each bioregion
write.table(sppocc_bioregion, file = "Data/species_per_bioregion.csv",
            sep=",", row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8")

# 4. Estimate species richness ---------------------------------------------------

# sum occurrence of each species for each grid cell nationally
national_grid_spp <- grid_species_length_df %>% 
  count(species, grid_id) 

# transpose so that species are rows and grid cells are columns
COL_grid_matrix <- national_grid_spp %>% 
  pivot_wider(names_from = "species", values_from = "n")

# account for sampling bias
# rarefraction calculations using vegan package to calculate Chao and ACE for each grid cell
COL_grid_matrix <- as.data.frame(COL_grid_matrix)
rownames(COL_grid_matrix) <- COL_grid_matrix$grid_id
COL_grid_matrix <- COL_grid_matrix[, -1]
COL_grid_matrix[is.na(COL_grid_matrix)] <- 0
veganCOL <- estimateR(COL_grid_matrix, smallsample = TRUE) 
veganCOL_df <- veganCOL %>% as.data.frame() %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("grid_id")

# add bioregion information
colnames(col_10_bioreg)
grid_bioreg <- col_10_bioreg %>% as.data.frame() %>%
  select(grid_id, bio_reg)
veganCOL_df$grid_id <- as.integer(veganCOL_df$grid_id)
veganCOL_bioreg <- left_join(veganCOL_df, grid_bioreg, by = "grid_id")

# export csv showing number of species and estimated species per grid cell, per bioregion
write.table(veganCOL_bioreg, file = "Data/spp_estimates_bioreg.csv",
            sep=",", row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8")

# 5. Identify qualifying UAs ----------------------------------------------

# call relevant data
# table of total species number and estimated richness in each grid cell, with bioregion data
spp_est <- veganCOL_bioreg
# number of species per bioregion
spp_bioreg <- sppocc_bioregion
colnames(spp_bioreg)

# find qualifying grids (according to Colombian IPA methodology for Criterion B)
# grids with highest estimated richness per bioregion
spp_est <- rename(spp_est, bioregion = bio_reg)
spp_est_highest <- spp_est %>%
  arrange(desc(S.chao1)) %>% 
  group_by(bioregion) %>%
  slice(1:1) %>%
  as.data.frame() %>%
  mutate(IPA_B1 = 1)

spp_est_highest

# grids with higest observed species richness per bioregion
spp_obs_highest <- spp_est %>%
  arrange(desc(S.obs)) %>% 
  group_by(bioregion) %>%
  slice(1:1) %>%
  as.data.frame()
spp_obs_highest

# add threshold level to qualifying grids
spp_bioreg <- spp_bioreg %>% 
  rename(bioreg_occ = occurences, bioreg_spec = species)
spp_est_highest_thres <- left_join(spp_est_highest, spp_bioreg, by = "bioregion") %>%
  mutate(meet_thres = if_else(S.obs >= spe_10per, "Y", "N"))

# find the bioregions where the most species rich grid cell 
# does not meet the 10% threshold
bioreg_notthres <- spp_est_highest_thres %>% 
  filter(meet_thres == "N") %>% select(bioregion)

# pull these bioregions from the list of grids with highest observed species richness per bioregion
spp_obs_highest2 <- spp_obs_highest %>% 
  filter(bioregion %in% bioreg_notthres$bioregion) %>%
  left_join(spp_bioreg, by = "bioregion")%>%
  mutate(meet_thres = if_else(S.obs >= spe_10per, "Y", "N")) %>%
  mutate(IPA_B2 = 1)

spp_est_highest_thres2 <- bind_rows(spp_est_highest_thres, spp_obs_highest2)
View(spp_est_highest_thres2)

# check that this doesn't duplicate grids
distinct(spp_est_highest_thres2, grid_id)
# all 19 there, so we're good

spp_est_highest_thres2$IPA_B1[is.na(spp_est_highest_thres2$IPA_B1)] <- 0
spp_est_highest_thres2$IPA_B2[is.na(spp_est_highest_thres2$IPA_B2)] <- 0

colnames(spp_est_highest_thres2)
spp_est_highest_thres2 <- spp_est_highest_thres2 %>% select(-IPA_B1, IPA_B1,
                                                            -IPA_B2, IPA_B2)

# export table of grids that potentially trigger IPAs for Criterion B
write.table(spp_est_highest_thres2, file = "Data/CUPC_IPACritB_grids_230127.csv",
            sep=",", row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8")



