# run-tximeta.R
# -----------------------------------------------------------------------------
# Author:             Stephanie Hicks
# Date last modified: Dec 5, 2020
#
# Use tximeta to read in quant files (output from salmon alevin) into SingleCellExperiment files

suppressPackageStartupMessages({
  library(here)
  library(tximeta)
  library(fishpond)
  library(alevinQC)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(org.Mm.eg.db) # org package for mouse
})

# create new data folder to store R objects
if(!file.exists(here("data"))){
  dir.create(here("data"))
}


# create linkedTranscriptome for decoys pipeline
index_dir = here("salmon_index_files", "gencode.vM25-salmon-index-v1.0.0-mouse-withdecoys-annotation-expanded")
fasta_path = here("salmon_index_files", "gencode.vM25.annotation.expanded.fa")
gtf_path = here("salmon_index_files", "gencode.vM25.annotation.expanded.gtf")
json_file = here("salmon_index_files", paste0(basename(index_dir), ".json"))
makeLinkedTxome(indexDir=index_dir, 
                source="GENCODE", organism="Mus musculus", 
                release="M25", genome="GRCm38", 
                fasta=fasta_path, gtf=gtf_path, 
                write=TRUE, jsonFile=json_file) # this command will add the index to the cache automatically

# Import with tximeta
# Note: salmon alevin can only import a single experiment/sample at a time; we have 6 samples here
sample_names = c("1C" , "2E", "3C", "4E" , "5C", "8E")
file_paths = here("01_quantification", "snrnaseq", "salmon_quants",  
                  paste0(sample_names, "_quant"), 
                  "alevin", "quants_mat.gz")
alevinQC_path = here("01_quantification", "snrnaseq", "salmon_quants",  
                     paste0(sample_names, "_quant"))

counts_mat = variance_mat = mean_mat = colnames_mat = coldata_mat <- NULL              
for(i in seq_len(length(sample_names))){
  print(file_paths[i])
  
  # looking at alevinQC output
  if(!file.exists(here("01_quantification", "snrnaseq", "alevinQC", 
                       paste0("alevinReport_", sample_names[i],".html")))){
    alevinQC::checkAlevinInputFiles(baseDir = alevinQC_path[i])
    alevinQCReport(baseDir = alevinQC_path[i], sampleId = sample_names[i], 
                   outputFile = paste0("alevinReport_", sample_names[i],".html"), 
                   outputFormat = "html_document",
                   outputDir = here("01_quantification", "snrnaseq", "alevinQC"), 
                   forceOverwrite = TRUE)
    # alevin_df <- readAlevinQC(baseDir = alevinQC_path[i])
    # head(alevin_df$cbTable)
    # dim(alevin_df$cbTable)
  }
  
  se = tximeta(coldata = data.frame(names = sample_names[i],
                                    files = file_paths[i],
                                    stringsAsFactors = FALSE),
               type = "alevin", 
               dropInfReps=TRUE)# , 
               # alevinArgs=list(filterBarcodes=TRUE))
  
  whitelist <- read.table(here("01_quantification", "snrnaseq", "salmon_quants",  
                               paste0(sample_names[i], "_quant"), 
                               "alevin", "whitelist.txt"),
                          stringsAsFactors = FALSE)
  colData(se)$whitelist <- colnames(se) %in% whitelist[,1]
  colData(se)$sample_id <- rep(sample_names[i], ncol(se))
  
  ##### Only need this if we want to have two separate assays of "spliced" and "intron" counts (e.g. useful for RNA velocity)
  # cg <- read.delim(here("salmon_index_files", 
  #                       "gencode.vM25.annotation.expanded.features.tsv"),
  #                  header = TRUE, as.is = TRUE)
  # ## Rename the 'intron' column 'unspliced' to make assay names compatible with scVelo
  # colnames(cg)[colnames(cg) == "intron"] <- "unspliced"
  # ses <- tximeta::splitSE(se, cg, assayName = "counts")
  
  # Save as a SummarizedExperiment object
  saveRDS(se, file = here("data", paste0("se_mouse_sleep_snrnaseq_", sample_names[i], ".rds")))
  
  counts_mat <- cbind(counts_mat, assay(se, "counts"))
  variance_mat <- cbind(variance_mat, assay(se, "variance"))
  mean_mat <- cbind(mean_mat, assay(se, "mean"))
  colnames_mat <- c(colnames_mat, colnames(se))
  coldata_mat <- rbind(coldata_mat, colData(se))
}

# Create new SE
se_new = SummarizedExperiment(assays = SimpleList(counts = counts_mat, 
                                                  variance = variance_mat, 
                                                  mean = mean_mat),
                              colData = coldata_mat,
                              rowRanges = rowRanges(se),
                              metadata = metadata(se))

# Convert to a SingleCellExperiment (sce) object
sce <- as(se_new, "SingleCellExperiment")

# add gene IDs to sce object
sce <- addIds(sce, "SYMBOL")
mcols(sce)

# Check sce object
colData(sce)
assayNames(sce)
rowRanges(sce)
table(sce$whitelist, sce$sample_id)

# > table(sce$whitelist, sce$sample_id)
# 
# 1C   2E   3C   4E
# FALSE 3650 4033 4393 3780
# TRUE  7350 6964 6606 7214

saveRDS(sce, file = here("data", paste0("sce_mouse_sleep_snrnaseq_complete.rds")))

