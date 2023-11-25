# IPA_identification
Steps for identifying potential Important Plant Areas (IPA).

**This Repository**

This github repository is based on the analysis undertaken to identify potential IPAs for useful plants in Colombia ([Kor and Diazgranados, 2023](https://www.researchgate.net/publication/372187721_Identifying_important_plant_areas_for_useful_plant_species_in_Colombia)). This was undertaken using existing occurrence records and applying IPA Criteria A and B. While there are many ways to approach IPA identification, we hope that this will prove useful for application in other countries and regions. 

We have split the analysis into four workbooks corresponding to: 1. downloading and cleaning data from GBIF and BIEN, 2. analysing data against IPA Criterion A thresholds, 3. analysing data against IPA Criterion B thresholds, 4. creating shapefiles of UAs which meet IPA thresholds to support mapping of outputs. 

The code has been somewhat simplified to summarise the steps taken and functions used. You can replace "path/to/data.csv", etc. with the actual file paths in your own system. Please note that thresholds applied were primarily based on the Colombian methodology for IPA identification and parameters may therefore need to be changed.

**Context**

The IPA programme was established in 2002 to identify and protect a network of best sites for plant conservation in the world. These are identified according to [globally consistent criteria](https://link.springer.com/article/10.1007/s10531-017-1336-6) based on the presence of A. threatened species, B. botanical richness and C. threatened habitats, usually identified at the national level. In this study, the criteria and thresholds applied were primarily based on the Colombian methodology for IPA identification.
