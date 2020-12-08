# download-gencode-files.R
# -----------------------------------------------------------------------------
# Author:             Stephanie Hicks
# Date last modified: Dec 5, 2020
#
# Download GENCODE files and generate files for the pipelines


suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(Biostrings)
  library(BSgenome)
  library(eisaR)
  library(GenomicFeatures)
  library(SummarizedExperiment)
})

# Create directories
if(!dir.exists(here("salmon_index_files"))){
  dir.create(here("salmon_index_files"))
}
if(!dir.exists(here("01_quantification", "bulk", "salmon_quants"))){
  dir.create(here("01_quantification", "bulk", "salmon_quants"))
}
if(!dir.exists(here("01_quantification", "snrnaseq", "salmon_quants"))){
  dir.create(here("01_quantification", "snrnaseq", "salmon_quants"))
}
if(!dir.exists(here("01_quantification", "snrnaseq", "alevinQC"))){
  dir.create(here("01_quantification", "snrnaseq", "alevinQC"))
}


# Download GENCODE Files
# download GENCODE primary assembly fasta file
if(!file.exists(here("salmon_index_files", "GRCm38.primary_assembly.genome.fa.gz"))){
  tar_gz_file <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz"
  download.file(tar_gz_file, 
                destfile = here("salmon_index_files", "GRCm38.primary_assembly.genome.fa.gz"), 
                method = "wget")
}

# download GENCODE transcripts fasta file
if(!file.exists(here("salmon_index_files", "gencode.vM25.transcripts.fa.gz"))){
  tar_gz_file <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz"
  download.file(tar_gz_file,
                destfile = here("salmon_index_files", "gencode.vM25.transcripts.fa.gz"),
                method = "wget")
}

# download GENCODE gtf file
if(!file.exists(here("salmon_index_files", "gencode.vM25.annotation.gtf.gz"))){
  tar_gz_file <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
  download.file(tar_gz_file, 
                destfile = here("salmon_index_files", "gencode.vM25.annotation.gtf.gz"), 
                method = "wget")
}


###############################################
### Files below this point are for snRNAseq ###
# Following https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/
# to generate reference fastq file
###############################################

# gunzip the GRCm38.primary_assembly.genome.fa.gz file
GEOquery::gunzip(filename = here("salmon_index_files", "GRCm38.primary_assembly.genome.fa.gz"), 
                 destname = here("salmon_index_files", "GRCm38.primary_assembly.genome.fa"))


# Next, we load the eisaR package and extract a GRanges object containing 
# the genomic coordinates of each annotated transcript and intron. 
# We use the ‘separate’ approach to define introns separately for each 
# transcript, and add a flank length of 90nt to each intron. We note that the 
# length of the flanking sequence should be chosen depending on the RNA read length 
# and the desired amount of overlap with an intron that is required to consider a 
# read potentially intronic. For more details on the different options, we refer to 
# the help of the getFeatureRanges() function.

gtf <- here("salmon_index_files", "gencode.vM25.annotation.gtf.gz")
grl <- eisaR::getFeatureRanges(
  gtf = gtf,
  featureType = c("spliced", "intron"), 
  intronType = "separate", 
  flankLength = 90L, 
  joinOverlappingIntrons = FALSE, 
  verbose = TRUE
)

grl[4:6]


# After defining the genomic positions of all features of interest, we can
# extract the sequences of these, and write to a fasta file for later 
# indexing with Salmon.

genome <- Biostrings::readDNAStringSet(
  here("salmon_index_files", "GRCm38.primary_assembly.genome.fa"))

names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = grl)

Biostrings::writeXStringSet(
  seqs, filepath = here("salmon_index_files", 
                        "gencode.vM25.annotation.expanded.fa"))


# To enable reading the estimated abundances with tximeta, automatically 
# recognizing the underlying transcriptome, we write the expanded annotation 
# to a GTF file. This will later be used to create a linked transcriptome for tximeta.
eisaR::exportToGtf(
  grl, 
  filepath = here("salmon_index_files", "gencode.vM25.annotation.expanded.gtf"))

# Since alevin quantifies spliced and unspliced features jointly, we will also 
# need to split the imported abundances by feature type. The splitting needs to 
# be done in such a way that we can still match up a spliced feature with the 
# corresponding unspliced feature. To help with this, the metadata of the 
# GRanges object contains a data frame with corresponding spliced and unspliced gene IDs.

head(metadata(grl)$corrgene)
##                spliced                 intron
## 1 ENSMUSG00000102693.1 ENSMUSG00000102693.1-I
## 2 ENSMUSG00000064842.1 ENSMUSG00000064842.1-I
## 3 ENSMUSG00000102851.1 ENSMUSG00000102851.1-I
## 4 ENSMUSG00000089699.1 ENSMUSG00000089699.1-I
## 5 ENSMUSG00000103147.1 ENSMUSG00000103147.1-I
## 6 ENSMUSG00000102348.1 ENSMUSG00000102348.1-I

write.table(
  metadata(grl)$corrgene, 
  file =  here("salmon_index_files",
               "gencode.vM25.annotation.expanded.features.tsv"),
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Finally, we create a text file mapping transcript and intron identifiers 
# to the corresponding gene identifiers.
df <- eisaR::getTx2Gene(
  grl, filepath = here("salmon_index_files",
                       "gencode.vM25.annotation.expanded.tx2gene.tsv")
)

