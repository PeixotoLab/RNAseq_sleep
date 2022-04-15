library(SingleR)
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(scran)
library(pheatmap)
library(EnsDb.Mmusculus.v79)
library(biomaRt)
library(BiocNeighbors)
library(HDF5Array)

## all cortex MOp Smart-seq####
allen.20 <- loadHDF5SummarizedExperiment(dir="/mnt/callisto/Zuin", prefix="Allen_mm_20")
allen.20.cortex <- as(allen.20, "SingleCellExperiment")

# select cortex, neurons-GABAergic and non-neuronal labels
allen.20.cortex <- allen.20.cortex[,!is.na(allen.20.cortex$subclass_label) & allen.20.cortex$subclass_label!="" & allen.20.cortex$subclass_label!="Meis2"] 
allen.20.cortex <- allen.20.cortex[,!(allen.20.cortex$subclass_label=="CA1-ProS" | allen.20.cortex$subclass_label=="CA2" | allen.20.cortex$subclass_label=="CA3" | allen.20.cortex$subclass_label=="DG" |
                                        allen.20.cortex$subclass_label=="CT SUB" | allen.20.cortex$subclass_label=="L2 IT ENTl"| allen.20.cortex$subclass_label=="L2 IT RHP"|
                                        allen.20.cortex$subclass_label=="L2/3 IT ENTl"|allen.20.cortex$subclass_label=="L2/3 IT PPP"|allen.20.cortex$subclass_label=="L3 IT ENT"|
                                        allen.20.cortex$subclass_label=="L3 RSP-ACA"|allen.20.cortex$subclass_label=="L5 IT TPE-ENT"|allen.20.cortex$subclass_label=="L5 PPP"|allen.20.cortex$subclass_label=="L6 IT ENTl"|
                                        allen.20.cortex$subclass_label=="L6b/CT ENT"|allen.20.cortex$subclass_label=="NP PPP"|allen.20.cortex$subclass_label=="NP SUB"|
                                        allen.20.cortex$subclass_label=="SUB-ProS"|allen.20.cortex$subclass_label=="V3d"|
                                        allen.20.cortex$subclass_label=="CR")]

symbol_allen <- rownames(allen.20.cortex)

map <- mapIds(EnsDb.Mmusculus.v79, keys= symbol_allen, keytype = "SYMBOL", column = "GENEID")
stopifnot(length(map) == nrow(allen.20.cortex))
rowData(allen.20.cortex)$symbol <- symbol_allen
rowData(allen.20.cortex)$ens <- map
allen.20.cortex<- allen.20.cortex[!is.na(rowData(allen.20.cortex)$ens),]
rownames(allen.20.cortex) <- rowData(allen.20.cortex)$ens

allen.10xgenomics <- aggregateAcrossCells(allen.20.cortex, use.assay.type = "counts",
                                        id=DataFrame(label=allen.20.cortex$subclass_label))
allen.10xgenomics <- logNormCounts(allen.10xgenomics)

save(allen.10xgenomics, allen.20.cortex, file = "allen_pseudo_10xgenomics.RData")