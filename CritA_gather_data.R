## Criterion A. Gather data ##


# Load libraries and data --------------------------------------------------
library(tidyverse)  # for data manipulation
library(rgbif)  # for downloading GBIF data
library(taxize)  # for matching GBIF names with taxon keys
library(BIEN) # to download BIEN data
library(magrittr)  # for piping
library(readr)  # for reading CSV files
library(data.table)  # for efficient data manipulation
library(bit64)  # for large integers
library(stringdist)  # for calculating string distances

# Load data
# CSV file containing a list of accepted plant species meeting Criterion A in your study area
cons_species <- read.csv("list_of_species.csv")

# Pull data from GBIF -----------------------------------------------------

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
write_csv(merged_data, "merged_data_gbif.csv")

# Pull data from BIEN -----------------------------------------------------

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




