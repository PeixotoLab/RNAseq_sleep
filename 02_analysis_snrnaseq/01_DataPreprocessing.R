# Data setting and pre-processing ####
# The goals of this analysis are:
# 1) add the UMI counts of spliced mRNA and introns
#    sharing the same Ensembl Gene ID;
# 2) quality control and doublets removal for each mouse.

# Bioconductor (v3.16) and R (v4.2.0)

# Set up
library(SingleCellExperiment)
library(EnsDb.Mmusculus.v79)
library(scuttle)
library(scran)
library(scDblFinder)

# Load single-nuclear RNA-seq snrna_dataset ####
snrna_data <- readRDS("sce_mouse_sleep_snrnaseq_complete.rds")

# The dataset contains the UMI counts of spliced mRNA and introns
# sharing the same Ensembl ID. So, we added the sum of them.
# First, we selected exons and introns.
exons <- snrna_data[-which(grepl("-I", rownames(snrna_data))), ]
rownames(exons) <- substring(rownames(exons), 1, 18)

introns <- snrna_data[which(grepl("-I", rownames(snrna_data))), ]
rownames(introns) <- substring(rownames(introns), 1, 18)

# We identified introns and exons with the same Ensembl ID
same_ensembl <- intersect(rownames(exons), rownames(introns))
exons <- exons[same_ensembl, ]
exons <- exons[order(rownames(exons)), ]
introns <- introns[same_ensembl, ]
introns <- introns[order(rownames(introns)), ]

rownames(snrna_data) <- substring(rownames(snrna_data), 1, 18)
# Remove the introns and exons with the same Ensembl ID
# from the single-nuclear data.
snrna_data <- snrna_data[! rownames(snrna_data) %in% same_ensembl, ]

# A new SingleCellExperiment object was created, where the sum of the UMI counts
# of spliced mRNA and introns sharing the same Ensembl ID was added.
sce <- SingleCellExperiment(assays = list(
  counts = rbind(counts(snrna_data), counts(exons) + counts(introns))
))

colData(sce) <- colData(snrna_data)
sce$sample_id[sce$sample_id == "8E"] <- "6E"

# We created the sleep condition variable
sce$condition <- sce$sample_id
sce$condition[grep("C", sce$condition)] <- "HC" # Home Cage (HC)
sce$condition[grep("E", sce$condition)] <- "SD" # Sleep Deprivated (SD)

sce$sample_id <- paste(substr(sce$sample_id, 1, 1), sce$condition, sep = "")

# To identify mitochondrial genes, we retrieved the chromosome location of
# each Ensembl Gene ID.
ensids <- rownames(sce)
map <- mapIds(EnsDb.Mmusculus.v79, keys = ensids,
              column = "SEQNAME", keytype = "GENEID")

stopifnot(length(map) == nrow(sce))
rowData(sce)$CHR <- map

# We split the data into six SingleCellExperiment objects,
# one for each mouse.
scelist <- list(sce[, sce$sample_id == "1HC"], sce[, sce$sample_id == "2SD"],
                sce[, sce$sample_id == "3HC"], sce[, sce$sample_id == "4SD"],
                sce[, sce$sample_id == "5HC"], sce[, sce$sample_id == "6SD"])

# Pre-processing ####
# First, for each sample, we detected low-quality and damaged droplets.
# We computed per-cell QC metrics. In particular, the sum of UMI counts,
# the number of detected genes, and the percentage of mitochondrial counts.
scelist_filt <- list()
for (i in seq_along(scelist)) {
  stats <- perCellQCMetrics(scelist[[i]], 
                            subsets = list(Mito = which(rowData(scelist[[i]])$CHR == "MT")))
  high_mito <- isOutlier(stats$subsets_Mito_percent, type = "higher")

  colData(scelist[[i]]) <- cbind(colData(scelist[[i]]), stats)
  scelist[[i]]$high_mito <- high_mito

  qc_lib <- isOutlier(scelist[[i]]$sum, log = TRUE, type = "lower")
  qc_nexprs <- isOutlier(scelist[[i]]$detected, log = TRUE, type = "lower")
  discard <- qc_lib | qc_nexprs | scelist[[i]]$high_mito

  scelist_filt[[i]] <- scelist[[i]][, !discard]
}

# Log-normalized counts
scelist_filt <- lapply(scelist_filt, function(x) logNormCounts(x))

# Lastly, for each sample, we removed potential doublets.
# HVGs were calculated for each sample
topgs <- lapply(scelist_filt, function(x) getTopHVGs(x, prop = 0.1))

for (i in seq_along(scelist_filt)) {
  set.seed(422)
  # Function to calculate the scores
  scores <- computeDoubletDensity(scelist_filt[[i]], subset.row = topgs[[i]])
  # Function to set the doublet scores threshold
  dbl_calls <- doubletThresholding(data.frame(score = scores), method = "griffiths", returnType = "call")

  colData(scelist_filt[[i]]) <- cbind.DataFrame(colData(scelist_filt[[i]]), dbl_calls, scores)
}

scelist_sgl <- lapply(scelist_filt, function(u) u[, !u$dbl_calls == "doublet"])
names(scelist_sgl) <- levels(factor(sce$sample_id))

# Save the new SCE object in .rds file
saveRDS(scelist_sgl, file = "snrna_scelist_sgl.rds")
