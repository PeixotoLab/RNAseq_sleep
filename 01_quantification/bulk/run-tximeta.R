# run-tximeta.R
# -----------------------------------------------------------------------------
# Author:             Stephanie Hicks
# Date last modified: Nov 28, 2020
#
# Use tximeta to read in quant files (output from salmon quant) into SummarizedExperiment files

suppressPackageStartupMessages({
  library(here)
  library(tximeta)
  library(fishpond)
  library(SummarizedExperiment)
  library(org.Mm.eg.db) # org package for mouse
})

# create new data folder to store R objects
if(!file.exists(here("data"))){
  dir.create(here("data"))
}

# create linkedTranscriptome for decoys pipeline
index_dir = here("salmon_index_files", "gencode.vM25-salmon-index-v1.0.0-mouse-withdecoys")
fasta_path = here("salmon_index_files", "gencode.vM25.transcripts.mouse.fa.gz")
gtf_path = here("salmon_index_files", "gencode.vM25.annotation.gtf.gz")
json_file = here("salmon_index_files", paste0(basename(index_dir), ".json"))
makeLinkedTxome(indexDir=index_dir, 
                source="GENCODE", organism="Mus musculus", 
                release="M25", genome="GRCm38", 
                fasta=fasta_path, gtf=gtf_path, 
                write=TRUE, jsonFile=json_file) # this command will add the index to the cache automatically

# Import with tximeta
# Note: salmon quant can import multiple files at a time, but salmon alevin is only one at a time

######### !! PICK ONE (need to ask Lucia and Davide about this)

# If you use SRA files, use this: 
pdata <- readr::read_csv(here("01_quantification", "bulk", "data", "SRA_geo_metadata.csv"))
file_paths = here("01_quantification", "bulk", "salmon_quants", 
                  paste0(pdata$Run, "_quant"), "quant.sf")

coldata <- data.frame(files=file_paths, names=pdata$Run,
                      genotype=factor(ifelse(pdata$`genotype:ch1`=="WT", "wt", "mutant"), levels = c("wt", "mutant")),
                      status=factor(c(rep("cage_control", 10), rep("sleep_dep",10))),
                      pdata,
                      stringsAsFactors=FALSE)
coldata

# If you use FASTQ files from Dario's server, use this: 
all_files <- list.files(here("01_quantification", "bulk", "salmon_quants"))
all_files <- stringr::str_subset(all_files, "^S3|^WT")
file_paths = here("01_quantification", "bulk", "salmon_quants", 
                  all_files, "quant.sf")

coldata <- data.frame(files=file_paths, names=stringr::str_sub(all_files, end=-7),
                      stringsAsFactors=FALSE)
coldata

######### 


# Import samples using tximeta into a SummarizedExperiment object
se <- tximeta(coldata, type = "salmon", txOut = TRUE)

# add gene IDs to se object
se <- addIds(se, "SYMBOL")
mcols(se)

# Check se object
colData(se)
assayNames(se)
rowRanges(se)

# Save as a SummarizedExperiment object
# saveRDS(se, file = here("data", "se_mouse_sleep.rds"))
saveRDS(se, file = here("data", "se_mouse_sleep_complete.rds"))


# Summarize to gene level counts
gse <- summarizeToGene(se)
gse <- addIds(gse, "SYMBOL", gene=TRUE)

colData(gse)
assayNames(gse)
rowRanges(gse)

# Save as a SummarizedExperiment object
# saveRDS(gse, file = here("data", "gse_mouse_sleep.rds"))
saveRDS(gse, file = here("data", "gse_mouse_sleep_complete.rds"))

