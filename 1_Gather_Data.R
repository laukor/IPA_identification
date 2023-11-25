## Download and collate occurence record data ##
# 1. Load libraries and species list
# 2. Get GBIF data
# 3. Get BIEN data
# 4. Merge and clean datasets
# 5. Identify potential HRE / RRE species
# 6. Export records

# 1. Load libraries and data --------------------------------------------------
library(tidyverse)  # for data manipulation
library(rgbif)  # for downloading GBIF data
library(taxize)  # for matching GBIF names with taxon keys
library(BIEN) # to download BIEN data
library(magrittr)  # for piping
library(readr)  # for reading CSV files
library(data.table)  # for efficient data manipulation
library(bit64)  # for large integers
library(stringdist)  # for calculating string distances
library(CoordinateCleaner) # to clean occurrence records
library(ConR) # to calculate EOO

# Load data
# CSV file containing a list of accepted plant species meeting Criterion A in your study area
cons_species <- read.csv("list_of_species.csv")

# NB to download occurrence records for all species in the study area (e.g. for criterion B analysis)
# upload list of all species names and skip step #5 

# 2. Get GBIF occurrence records -----------------------------------------------

# Enter GBIF credentials
user <- "gbif_username"
pwd <- "gbif_password"
email <- "email_address"

# Get GBIF species keys 
gbif_upfp_keys <- read_csv("list_of_species.csv") %>%
  pull("taxon_name") %>% # Match the species names to GBIF backbone to get taxon keys
  get_gbifid_(method = "backbone") %>%
  imap(~ .x %>% mutate(original_sciname = .y)) %>%
  bind_rows() %>%
  filter(rank == "species") %>% 
    # only keep records idenitified to species level 
  mutate(dist = stringdist(canonicalname, original_sciname, method = "dl")) %>% 
    # only keep records with a close match to the original name
  filter(dist <= 1) %>%
  pull(usagekey)

# Download occurrence records 
# Use the matched GBIF species keys to download occurrence data from GBIF
occ_download(
  pred_in("taxonKey", gbif_upfp_keys),
  pred("hasCoordinate", TRUE),
  # only include records with coordinates
  format = "SIMPLE_CSV",
  user = user, pwd = pwd, email = email
)

# Match GBIF synonyms back to accepted names inputted
# Load downloaded occurrence data
records <- fread("occurrence_data.csv", header = TRUE, encoding = "UTF-8")
# Remove 'absent' records
records <- filter(records, occurrenceStatus == "PRESENT")
# Rename the column "speciesKey" to "usagekey"
colnames(records)[colnames(records) == "speciesKey"] <- "usagekey"
# Load the list of inputted names from your data and the specieskeys used by GBIF
keys <- read_csv("list_of_species.csv")

# Merge the occurrence data with the list of inputted accepted species names
merged_data <- inner_join(records, keys, by = "usagekey") %>%
  select(taxon_name, occurrenceID, decimalLongitude, decimalLatitude) %>%
  # Keep only relevant columns
  distinct()

# Write the merged data to a CSV file
write_csv(merged_data, "data_gbif.csv")

# 3. Get BIEN occurrence records -----------------------------------------------------

# Extract accepted target species names
species <- cons_species$taxon_name

# Download global occurrences for all species
occ <- BIEN_occurrence_species(species, 
                               all.taxonomy = TRUE, 
                               native.status = TRUE,
                               natives.only = FALSE, 
                               political.boundaries = TRUE)

# Filter out records from GBIF
occ_bien <- subset(occ, datasource != "GBIF") %>%
  subset(occ_bien, scrubbed_species_binomial != "NA")
  # Only keep records at species level

# Get and export data citations
citations <- BIEN_metadata_citation(occ_bien_species,
                                    bibtex_file = "path/BIEN_citations")
  # saves as bib file

# Remove unneeded columns according to your needs
# Export records
write.table(occ_bien_species, 
            file = "path/data_BIEN.csv",
            sep = ",", 
            row.names = FALSE, 
            col.names = TRUE, 
            fileEncoding = "UTF-8")

# 4. Merge and clean data -----------------------------------------------------

# Load occurrence records data
load_data <- function(file_path) {
  read_csv(file_path, locale = locale(encoding = "UTF-8"))
}

gbif_records <- load_data(here::here("data_gbif.csv"))
bien_records <- load_data(here::here("data_BIEN.csv"))
other_records <- load_data(here::here("data_other.csv")) # any other local records

# Ensure all data sets have matching columns, as relevant for your analyses

# Merge data sources
all_occ <- bind_rows(gbif_records, bien_records, other_records)

# Remove records with no coordinates
all_occ_clean <- all_occ %>% 
  filter(decimalLatitude != "NA", decimalLongitude != "NA")

# Spatially clean records
clean_coordinates <- cc_equ(data, lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude",
          species = "taxon_name", additions = NULL, verbose = TRUE) %>%
  cc_gbif(lon = "decimalLongitude", lat = "decimalLatitude",
          species = "taxon_name") %>%
  cc_zero(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_cen(lon = "decimalLongitude", lat = "decimalLatitude", 
         species = "taxon_name") %>%
  cc_val(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_inst(lon = "decimalLongitude", lat = "decimalLatitude", 
          species = "taxon_name") %>%
  cc_gbif(lon = "decimalLongitude", lat = "decimalLatitude", 
          species = "taxon_name")

# Combine with other data if relevant 
# e.g. in Colombia, we added data on uses, IUCN and national threat levels, etc.

# Output all clean records to CSV
write.table(clean_coordinates, file = "allrecords.csv",
            sep=",", row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8")

# NB for Criterion B analysis, stop here

# 5. Identify potential HRE/RRE species --------------------------------------
# Pull out occurrence records for endemic species with no extinction risk assessment

# load all clean occurrence records 
occ <- read_csv("allrecords.csv", 
                as.is=TRUE, header=TRUE, encoding = "UTF-8")
# load list of Criterion A species
critA <- read_csv("list_of_species.csv", 
                  encoding = "UTF-8")

# filter data to species which are endemic and have no extinction risk assessment
occ_end <- occ %>% 
  subset(grepl("Y", Endemic) & (grepl("NA", IUCNstatus) | 
        grepl("DD", IUCNstatus) | grepl("NA", National_RL) | 
        grepl("DD", National_RL) | grepl("NE", National_RL)))

# check how many species are in this category
end_notassessed <- unique(occ_end$Accepted_species_name)

# remove records with missing occurrence data
occs_end <- drop_na(occ_end, decimalLongitude, decimalLatitude)

# export file of occurrences of endemic species with no extinction risk
write.table(occs_end, file = "endemic_noassessment.csv",
            sep=",", row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8")

# convert to shapefile object
records_sf <- occs_end %>% 
  dplyr::select(Accepted_species_name,decimalLatitude,
  decimalLongitude, datasource_id)

occs_end_sf <- st_as_sf(records_sf, coords = c("decimalLongitude", "decimalLatitude"),
                       stringsAsFactors=TRUE, crs = 4326) #the EPSG code for WGS84

# identify HRE and RRE species
# calculate EOO for each species
# use sink to capture outputs, including any warnings for species where EOO can't be calculated
sink(file="Data/code_outputs.txt")
EOO.computing(dataset.ex, write_results = TRUE, file.name = "Data/EOO")
sink()

# identify HRE and RRE 
range <-read.csv("Data/EOO.csv", 
                   encoding = "UTF-8")

# replace NA values with 0 (these are species with <3 occurrences)
range[is.na(range)] <- 0

# add columns for HRE and RRE based on estimated EOOs
range_HRE <- range %>% mutate(IPA_A4 = if_else(EOO < 100 | EOO == "NA", "Y", "N")) %>%
  mutate(IPA_A5 = if_else((EOO > 100 & EOO < 5000), "Y", "N"))

# 6. join with other datasets and export ---------------------------------------------------

# change column names to match other dataframes
colnames(range_HRE)[colnames(range_HRE) == "X"] <- "taxon_name"

# join with list of all potential Crit A species
critA_EOO <- left_join(critA, range_HRE, by = "taxon_name")
critA_EOO[is.na(critA_EOO)] <- "NA"

head(critA_EOO)

# filter to species meeting an IPA trigger criterion
critA_EOO_cons <- critA_EOO %>% subset(grepl("Y", IPA_A1) | grepl("Y", IPA_A3) |
                                         grepl("Y", IPA_A4) | grepl("Y", IPA_A5))

# join with occurrence records
colnames(occ_end_HRE)[colnames(occ_end_HRE) == "taxon_name"] <- "Accepted_species_name"
occ_EOO <- left_join(occ, occ_end_HRE, by = "Accepted_species_name")
occ_EOO[is.na(occ_EOO)] <- "NA"

# filter to only potential trigger species
occ_EOO_cons <- occ_EOO %>% subset(grepl("Y", IPA_A1) | grepl("Y", IPA_A3) |
                                       grepl("Y", IPA_A4) | grepl("Y", IPA_A5))

# export list of IPA trigger species
write.table(critA_EOO_cons, file = "Data/IPACritA.csv",
            sep=",", row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8")

# list of all occurrences with 
write.table(occ_EOO_cons, file = "Data/occurrences/IPACritA_records.csv",
            sep=",", row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8")

