library(SingleCellExperiment)
library(EnsDb.Mmusculus.v75)

#### dataset setting#####
sce <- readRDS("/home/zuin/sce_mouse_sleep_snrnaseq_complete.rds")
sce$sample_id[sce$sample_id=="8E"] <- "6E"

# create condition variable
sce$Sleep <- sce$sample_id
sce$Sleep[sce$Sleep=="1C"|sce$Sleep=="3C"|sce$Sleep=="5C"] <- "HC" 
sce$Sleep[sce$Sleep=="2E"|sce$Sleep=="4E"|sce$Sleep=="6E"] <- "SD"

sces <- list(
  sce[,sce$sample_id=="1C"],
  sce[,sce$sample_id=="2E"],
  sce[,sce$sample_id=="3C"],
  sce[,sce$sample_id=="4E"],
  sce[,sce$sample_id=="5C"],
  sce[,sce$sample_id=="6E"]
)

exons <- introns <- same.ensembl <- exons.ensembl <- introns.ensembl<- sce.sum <- list()
for (n in seq_along(sces)) {
  exons[[n]] <- sces[[n]][-which(grepl("-I", rownames(sces[[n]]))),]
  introns[[n]] <- sces[[n]][which(grepl("-I", rownames(sces[[n]]))),]
  rownames(exons[[n]]) <- substring(rownames(exons[[n]]), 1, 18)
  rownames(introns[[n]]) <- substring(rownames(introns[[n]]), 1, 18)
  
  same.ensembl <- intersect(rownames(exons[[n]]),rownames(introns[[n]]))
  exons.ensembl[[n]] <- exons[[n]][same.ensembl,]
  exons.ensembl[[n]] <- exons.ensembl[[n]][order(rownames(exons.ensembl[[n]])),]
  introns.ensembl[[n]] <- introns[[n]][same.ensembl,]
  introns.ensembl[[n]] <- introns.ensembl[[n]][order(rownames(introns.ensembl[[n]])),]
  
  exons[[n]] <- exons[[n]][! rownames(exons[[n]]) %in% same.ensembl,]
  introns[[n]] <- introns[[n]][! rownames(introns[[n]]) %in% same.ensembl,]
  
  
  sce.sum[[n]] <- SingleCellExperiment(assays=list(counts=rbind(counts(exons[[n]]),
                                                                counts(exons.ensembl[[n]])+counts(introns.ensembl[[n]]),
                                                                counts(introns[[n]]))))

  colData(sce.sum[[n]]) <- colData(sces[[n]])
}

for (n in seq_along(sce.sum)) {
  # convert mouse ensembl to chromosome identifier 
  ensids <- sapply(strsplit(rownames(sce.sum[[n]]), ".", fixed = TRUE), function(x) x[1])
  
  map <- mapIds(EnsDb.Mmusculus.v75, keys = ensids, column = "SEQNAME", keytype = "GENEID")
  stopifnot(length(map) == nrow(sce.sum[[n]]))
  rowData(sce.sum[[n]])$CHR <- map
  
  no.mito <- sce.sum[[n]][-which(rowData(sce.sum[[n]])$CHR=="MT"),]
  mito <- sce.sum[[n]][which(rowData(sce.sum[[n]])$CHR=="MT"),]
  mito <- SingleCellExperiment(assays=list(counts=counts(mito)+1), 
                               rowData=rowData(mito), colData=colData(mito))
  sce.sum[[n]] <- SingleCellExperiment(assays=list(counts=rbind(counts(no.mito),counts(mito))),
                                       rowData=rbind(rowData(no.mito),rowData(mito)), 
                                       colData=colData(sce.sum[[n]]))

}

save(sce.sum, file = "sce_sum.RData")