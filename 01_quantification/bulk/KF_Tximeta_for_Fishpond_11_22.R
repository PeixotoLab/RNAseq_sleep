# Katie Ford
# November 2022 (Reorganized folder on Dropbox and needed to run again. I used a new Mac and updated the packages as seen below)

# I adapted this code from Stephanie Hicks
# Use tximeta to read in quant files (output from salmon quant) into SummarizedExperiment files

# This is for downstream analysis with Fishpond

#### Set the working directory and load packages ####

suppressPackageStartupMessages({
  library(here) # Version 1.0.1
  library(tximeta) # Version 1.16.0
  library(fishpond) # Version 2.4.0
  library(SummarizedExperiment) # Version 1.28.0
  library(org.Mm.eg.db) # Version 3.16.0
  library(readr) # Version 2.1.3
})

#### Create a new data folder to store R objects ####

if(!file.exists(here("data"))){
  dir.create(here("data"))
}

#### Create linkedTranscriptome for decoys pipeline ####
index_dir = here("IndexFiles", "gencode.vM25-salmon-index-v1.0.0-mouse-withdecoys")
fasta_path = here("IndexFiles", "gencode.vM25.transcripts.mouse.fa.gz")
gtf_path = here("IndexFiles", "gencode.vM25.annotation.gtf.gz")
json_file = here("IndexFiles", paste0(basename(index_dir), ".json")) 
makeLinkedTxome(indexDir=index_dir, 
                source="GENCODE", organism="Mus musculus", 
                release="M25", genome="GRCm38", 
                fasta=fasta_path,
                gtf=gtf_path, 
                write=TRUE, jsonFile=json_file) # this command will add the index to the cache automatically

#### Import with tximeta ####

# See the path to the quantification files here, as they are stored on Dropbox:
all_files <- list.files(here("data", "Salmon_Quants"))
all_files <- stringr::str_subset(all_files, "^S3|^WT")
file_paths = here("data", "Salmon_Quants", 
                  all_files, "quant.sf")

coldata <- data.frame(files=file_paths, names=stringr::str_sub(all_files),
                      condition = as.factor(rep(c("S3HC5", "S3HC7", "S3RS2", "S3SD5", "WTHC5", "WTHC7", "WTRS2", "WTSD5"), each = 5)),
                      stringsAsFactors=FALSE)


# View coldata
coldata

##### Import samples using tximeta into a SummarizedExperiment object ####

# I initially added countsFromAbundance=scaledTPM based off of tutorials from 
# DRIMseq and Fishpond, however, this is redundant with UQ normalization, so I 
# have removed it from here as that will be included in downstream analysis
se <- tximeta(coldata, type = "salmon", txOut = TRUE, useHub = FALSE)

#### Isoform #### 
# add gene IDs to se object (Isoform)
se <- addIds(se, "SYMBOL")
mcols(se)

# Check se object
colData(se)
assayNames(se)
rowRanges(se)

#### Save as a SummarizedExperiment object ####
saveRDS(se, file = here("data", "transcript_summarized_experiment.rds"))

write.table(colData(se), file = "coldata.txt")

# Check to make sure these match
colnames(se) == coldata[,2]

#### Gene level for GEO submission, not used for Fishpond ####
# note: we drop inferential replicates for gene-level analysis
se_gene <- tximeta(coldata, type = "salmon", txOut = TRUE, dropInfReps = TRUE, useHub = FALSE)

# Add gene IDs for gene level matrix
se_gene <- addIds(se_gene, "SYMBOL")
mcols(se_gene)

# Check se_gene object
colData(se_gene)
assayNames(se_gene)
rowRanges(se_gene)

# Summarize to gene level counts
gse <- summarizeToGene(se_gene)
gse <- addIds(gse, "SYMBOL", gene = TRUE)

colData(gse)
assayNames(gse)
rowRanges(gse)

saveRDS(gse, file = here("data", "gene_summarized_experiment.rds"))

sink('112822_Fishpond_SessionInfo.txt')
sessionInfo()
sink() 
