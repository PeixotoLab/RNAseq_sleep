# download-geo-data.R
# -----------------------------------------------------------------------------
# Author:             Stephanie Hicks
# Date last modified: Dec 5, 2020
#
# Download GEO data

library(here)
library(tidyverse)

# Create directories
if(!dir.exists(here("01_quantification", "log"))){
  dir.create(here("01_quantification", "log"))
}
if(!dir.exists(here("01_quantification", "bulk", "log"))){
  dir.create(here("01_quantification", "bulk", "log"))
}
if(!dir.exists(here("01_quantification", "snrnaseq", "log"))){
  dir.create(here("01_quantification", "snrnaseq", "log"))
}
if(!dir.exists(here("fastq_files"))){
  dir.create(here("fastq_files"))
}


# Download phenotype table
gse <- GEOquery::getGEO("GSE113754")
sapply(gse, dim) # Number of rows per entry
gse <- gse[[1]] # this contains all the samples
pdata <- Biobase::pData(gse) # Get phenotype table
pdata$Experiment <- sapply(stringr::str_split(pdata$relation.1, "="), tail, 1)

# To get SRR IDs, go to https://www.ncbi.nlm.nih.gov/bioproject/PRJNA453921. 
# Then click on `SRA Experiments`. Click `Send to`. Choose `File`. 
# Change Format to `RunInfo`. Click `Create File`. 
# This will download a file called `SraRunInfo.csv`. Read this file in. 
sra <- readr::read_csv(here("01_quantification", "bulk", "data", "SraRunInfo.csv"))
pdata <- dplyr::left_join(pdata, sra, by = "Experiment")

# Save phenotype table
readr::write_csv(pdata, file = here("01_quantification", "bulk", "data", "pdata.csv"))

# Shows there are 10 WT and 10 Shank3 mutant samples
table(pdata$`genotype:ch1`)
# > table(pdata$`genotype:ch1`)
#   Shank3 Mutant          WT 
#              10          10 


# Select subset of phenotype table
sra_meta <- pdata %>% 
  select(BioProject, BioSample, Submission, Run, Experiment, SRAStudy, Sample, 
         title, geo_accession,  
         source_name_ch1, organism_ch1, starts_with("characteristics"), 
         molecule_ch1, taxid_ch1, description,
         platform_id, instrument_model, library_selection, library_source, 
         library_strategy, LibraryLayout, 
         spots, bases, spots_with_mates, avgLength, size_MB, download_path,
         `tissue:ch1`, `age:ch1`, `Sex:ch1`, `genotype:ch1`)

# All SRA files
readr::write_csv(sra_meta, file = here("01_quantification", "bulk", "data", "SRA_geo_metadata.csv"))

# Write to files just the SRR IDs and paths
write.table(pdata$Run, file = here("01_quantification", "bulk", "data", "SRR_files.txt"), 
            quote= FALSE,row.names = FALSE, col.names = FALSE)
write.table(pdata$download_path, file = here("01_quantification", "bulk", "data", "SRR_paths.txt"), 
            quote= FALSE,row.names = FALSE, col.names = FALSE)