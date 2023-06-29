library(SingleCellExperiment)
library(EnsDb.Mmusculus.v79)
library(scuttle)
library(scran)
library(scDblFinder)

# Bioconductor (v3.16) and R (v4.2.0)
# Load single-nuclear RNA-seq dataset ####
data <- readRDS("/home/zuin/sce_mouse_sleep_snrnaseq_complete.rds")

# We added the UMI counts of spliced mRNA and introns sharing the same ensembl ID
# Exons
exons <- data[-which(grepl("-I", rownames(data))),]
rownames(exons) <- substring(rownames(exons), 1, 18)
# Introns
introns <- data[which(grepl("-I", rownames(data))),]
rownames(introns) <- substring(rownames(introns), 1, 18)

# We identified introns and exons with the same ensembl
same.ensembl <- intersect(rownames(exons), rownames(introns))
exons <- exons[same.ensembl,]
exons <- exons[order(rownames(exons)),]
introns <- introns[same.ensembl,]
introns <- introns[order(rownames(introns)),]

rownames(data) <- substring(rownames(data), 1, 18)
data <- data[! rownames(data) %in% same.ensembl,] # remove the introns and exons with the same ensembl

# A new SingleCellExperiment object was created, where the UMI counts of spliced mRNA and introns sharing the same ensembl ID were added
sce <- SingleCellExperiment(assays=list(counts=rbind(counts(data), counts(exons)+counts(introns))))
colData(sce) <- colData(data)

sce$sample_id[sce$sample_id=="8E"] <- "6E"

# We created the sleep condition variable
sce$condition <- sce$sample_id
sce$condition[grep("C", sce$condition)] <- "HC" # Home Cage (HC)
sce$condition[grep("E", sce$condition)] <- "SD" # Sleep Deprivated (SD)

sce$sample_id <- paste(substr(sce$sample_id,1,1), sce$condition, sep = "")

# To identify mitochondrial genes, we retrieve the chromosome location of each ensembl gene
ensids <- rownames(sce)
map <- mapIds(EnsDb.Mmusculus.v79, keys = ensids, column = "SEQNAME", keytype = "GENEID")
stopifnot(length(map) == nrow(sce))
rowData(sce)$CHR <- map

# We split the data into six SingleCellExperiment objects, one for each mouse.
sces <- list(sce[,sce$sample_id=="1HC"], sce[,sce$sample_id=="2SD"], sce[,sce$sample_id=="3HC"],
             sce[,sce$sample_id=="4SD"], sce[,sce$sample_id=="5HC"], sce[,sce$sample_id=="6SD"])

# Pre-processing ####
# For each sample, we detected low quality and damaged droplets. 
# Particularly, we computed per-cell quality-control metrics; these metrics include the sum of UMI counts, 
# the number of detected genes, and the percentage of mitochondrial counts.
filter.sce <- list()
for(i in seq_along(sces)) {
  print(i)
  stats <- perCellQCMetrics(sces[[i]], subsets=list(Mito=which(rowData(sces[[i]])$CHR=="MT")))
  high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
  
  colData(sces[[i]]) <- cbind(colData(sces[[i]]), stats)
  sces[[i]]$high.mito <- high.mito
  
  qc.lib <- isOutlier(sces[[i]]$sum, log=TRUE, type="lower")
  qc.nexprs <- isOutlier(sces[[i]]$detected, log=TRUE, type="lower")
  discard <- qc.lib | qc.nexprs | sces[[i]]$high.mito
  
  filter.sce[[i]] <- sces[[i]][,!discard]
  filter.sce[[i]] <- logNormCounts(filter.sce[[i]]) # Normalization
}

# Remove doublets####
# HVGs were calculated for each sample
topgs <- lapply(filter.sce, function(x) getTopHVGs(x, prop=0.1))

for (i in seq_along(filter.sce)) {
  print(i)
  set.seed(422)
  scores <- computeDoubletDensity(filter.sce[[i]], subset.row=topgs[[i]]) # function to calculate the scores
  dbl.calls <- doubletThresholding(data.frame(score=scores), method="griffiths", returnType="call") # function to set the doublet scores threshold
  
  colData(filter.sce[[i]]) <- cbind.DataFrame(colData(filter.sce[[i]]), dbl.calls, scores)
}

sces.singlet <- lapply(filter.sce, function(u) u[,!u$dbl.calls=="doublet"])
names(sces.singlet) <- levels(factor(sce$sample_id))

# Save singlet cells in .rds file 
saveRDS(sces.singlet, file = "snRNA_DataPreprocessing.rds")
