library(SingleCellExperiment)
library(scuttle)
library(HDF5Array)
library(dplyr)
library(Seurat)

#devtools::install_github("satijalab/seurat-data")
#devtools::install_github("satijalab/azimuth", ref = "release/0.4.6")

library(Azimuth)
library(SeuratData)
library(ggplot2)
library(biomaRt)
library(SingleR)

# Reference dataset ####
# We used Allen Whole Cortex & Hippocampus - 10x genomics (v 2021) as reference dataset for cell annotation
reference <- loadHDF5SummarizedExperiment(dir="/mnt/callisto/Zuin", prefix="Allen_mm_21")
reference <- as(reference, "SingleCellExperiment")
names(assays(reference)) <- c("counts")

# We selected Non-Neuronal, Neurons Glutamatergic and Neurons GABA-ergic, coming from Visual cortex (VIS, VISl, VISm, VISp) region
reference <- reference[,c(reference$region_label=="VIS"|reference$region_label=="VISl"|reference$region_label=="VISm"|
                            reference$region_label=="VISp")]
reference <- reference[,!is.na(reference$subclass_label) & reference$subclass_label!=""]
reference <- reference[,-which(grepl("ENT", reference$subclass_label))]
reference <- reference[,-which(grepl("PPP", reference$subclass_label))]
reference <- reference[,-which(grepl("CR", reference$subclass_label))]
reference <- reference[,-which(grepl("Meis", reference$subclass_label))]
reference <- reference[,-which(grepl("SUB", reference$subclass_label))]

# For computational issues, we decided to select 100000 cortical random cells.
# We selected cell type with less 100 cells. These weren't randomly selected.
no.random <- reference[,c(reference$subclass_label=="SMC-Peri"|reference$subclass_label=="VLMC")]

# Cell type with more 100 cells were selected
reference <- reference[,!c(reference$subclass_label=="SMC-Peri"|reference$subclass_label=="VLMC")]

df <- data.frame(colData(reference))
# We randomly selected 100 cells for each cell type
set.seed(23)
random <- df %>% group_by(subclass_label) %>% slice_sample(n=100)
random100 <- reference[,colnames(reference) %in% random$sample_name]

reference <- reference[,!(colnames(reference) %in% random$sample_name)]

# We randomly selected 98100 cells for each cell type
set.seed(23)
reference <- reference[,sample(colnames(reference), 98046)]

reference <- cbind(reference, random100, no.random)

# Cell Annotation with Azimuth ####
# The reference dataset was converted into Seurat object
counts <- as.matrix(counts(reference))
colData <- colData(reference)
reference.so <- CreateSeuratObject(counts = counts, meta.data = data.frame(colData))

# And we created a compatible object for Azimuth cell annotation
referenceSC <- SCTransform(reference.so,assay = "RNA", new.assay.name = "SCT", variable.features.n = 2000,
                           verbose = TRUE,conserve.memory=TRUE) # conserve.memory=TRUE for large dataset (no crash R)
referenceSC <- RunPCA(referenceSC, assay = "SCT", npcs = 50, verbose = FALSE, reduction.name = "PCA", return.model=TRUE) # Number of PCs must be 50 
referenceSC <- RunUMAP(referenceSC, assay = "SCT", reduction = "PCA", dims = seq_len(50), seed.use = 1, verbose = FALSE,
                       reduction.name = "umap",return.model=TRUE) # "umap" because not read
referenceSC$subclass_label <- as.factor(referenceSC$subclass_label)
Idents(object = referenceSC) <- "subclass_label"

reference.azimuth <- AzimuthReference(referenceSC, refUMAP = "umap", refDR = "PCA", refAssay = "SCT", dims = 1:50, 
                                      metadata = c("subclass_label"), verbose = TRUE)
# save reference in a folder called "reference"
ref.dir <- "reference/"
SaveAnnoyIndex(object = reference.azimuth[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
saveRDS(object = reference.azimuth, file = file.path(ref.dir, "ref.Rds"))

# Query dataset 
data <- readRDS("snRNA_DataPreprocessing.rds")
# The samples were merged and were converted into a Seurat object
counts <- cbind(counts(data[[1]]),counts(data[[2]]),counts(data[[3]]),counts(data[[4]]),counts(data[[5]]),counts(data[[6]]))

sample_id <- c(data[[1]]$sample_id, data[[2]]$sample_id, data[[3]]$sample_id, data[[4]]$sample_id, 
               data[[5]]$sample_id, data[[6]]$sample_id)

condition <- c(data[[1]]$condition, data[[2]]$condition, data[[3]]$condition, data[[4]]$condition, 
               data[[5]]$condition, data[[6]]$condition)

data <-  SingleCellExperiment(assays=list(counts=counts))
data$sample_id <- sample_id
data$condition <- condition
colnames(data) <- paste(colnames(data), data$sample_id, sep="_")

# Dataset was converted into a Seurat object
seurat.obj <- CreateSeuratObject(counts = counts(data), meta.data = data.frame(colData(data)), project = "SD")

# Cell annotation with RunAzimuth on all cells
seurat.obj <- RunAzimuth(seurat.obj, reference = "reference/")

# Projections 
seurat.obj <- RunTSNE(seurat.obj, reduction = "integrated_dr", dims = seq_len(20),
                      seed.use = 1, do.fast = TRUE, verbose = FALSE, reduction.name = "TSNE")

seurat.obj <- RunUMAP(seurat.obj, reduction = "integrated_dr", dims = seq_len(20),
                      seed.use = 1, verbose = FALSE, reduction.name = "UMAP")

# Cell Annotation with SingleR ####
# For SingleR cell annotation the reference dataset was aggregated across groups of cell type and was normalized.
aggregref <- aggregateAcrossCells(reference, use.assay.type = "counts",id=DataFrame(label=reference$subclass_label))
aggregref <- logNormCounts(aggregref)

# Symbol was converted to ensembl
ensembl.biomaRt <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
symbol.allen <- rownames(aggregref)

ensembl.id <- getBM(attributes=c('external_gene_name','ensembl_gene_id'), 
                    filters = 'external_gene_name', values = symbol.allen, mart = ensembl.biomaRt)
ensembl.id <- ensembl.id[!duplicated(ensembl.id$external_gene_name),] # remove duplicated symbol genes

inters <- intersect(rownames(aggregref), ensembl.id$external_gene_name)
ensembl.id <- ensembl.id[ensembl.id$external_gene_name %in% inters,]
ensembl.id <- ensembl.id[order(ensembl.id$external_gene_name),]

aggregref <- aggregref[rownames(aggregref) %in% inters,]
aggregref <- aggregref[order(rownames(aggregref)),]
rowData(aggregref)$ensembl_id <- ensembl.id$ensembl_gene_id
rownames(aggregref) <- ensembl.id$ensembl_gene_id

# Normalized data
data <- logNormCounts(data)

pred.SigleR <- SingleR(data, ref=aggregref, labels=aggregref$subclass_label)

# Save SingleR results into Seurat object for visualization
seurat.obj$predicted.SingleR <- pred.SingleR$labels

# Save cell-type annotation results into SingleCellExperiment object
data$Azimuth.labels <- seurat.obj$predicted.subclass_label
data$SingleR.labels <- pred.SigleR$labels

# Visualization with t-SNE plot and UMAP plot####
# Allen color labels
allen.color <- data.frame(aggregref$subclass_label, aggregref$subclass_color, aggregref$class_label)
colnames(allen.color) <- c("Label", "color","class")

# change Astro color
allen.color$color[allen.color$color=="#665C47"] <- "#957b46"
# change Oligo color
allen.color$color[allen.color$color=="#53776C"] <- "#744700"
# change SMC-Peri color
allen.color$color[allen.color$color=="#807059"] <- "#4c1130"
# change VLMC color
allen.color$color[allen.color$color=="#697255"] <- "#a9bd4f"
# change Endo color
allen.color$color[allen.color$color=="#8D6C62"] <- "#c95f3f"
# change Sncg color 
allen.color$color[allen.color$color=="#D3408D"] <- "#ffff00"

subclass.color <- allen.color$color
names(subclass.color) <- allen.color$Label

# save Allen color labels: 
save(subclass.color, file = "AllenColorLabel.RData")

# Azimuth
plot.so <- DimPlot(seurat.obj, reduction = "TSNE", group.by = "predicted.subclass_label", pt.size=0.5)+ 
  NoLegend()+ labs(x = "TSNE1", y="TSNE2") + ggtitle("") + scale_color_manual(values=subclass.color)
plot.so <- LabelClusters(plot.so, id = "predicted.subclass_label",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_TSNEplot_Azimuth_NoLegend.pdf",sep = ""), width = 20, height = 20, units = "cm")

plot.so <- DimPlot(seurat.obj, reduction = "UMAP", group.by = "predicted.subclass_label", pt.size=0.5)+
  NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") + scale_color_manual(values=subclass.color)
plot.so <- LabelClusters(plot.so, id = "predicted.subclass_label",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_UMAPplot_Azimuth_NoLegend.pdf",sep = ""), width = 20, height = 20, units = "cm")

# SingleR
plot.so <- DimPlot(seurat.obj, reduction = "TSNE", group.by = "predicted.SingleR", pt.size=0.5)+
  NoLegend()+ labs(x = "TSNE1", y="TSNE2") + ggtitle("") + scale_color_manual(values=subclass.color)
plot.so <- LabelClusters(plot.so, id = "predicted.SingleR",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_TSNEplot_SingleR_NoLegend.pdf",sep = ""), width = 20, height = 20, units = "cm")

plot.so <- DimPlot(seurat.obj, reduction = "UMAP", group.by = "predicted.SingleR", pt.size=0.5)+
  NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") + scale_color_manual(values=subclass.color)
plot.so <- LabelClusters(plot.so, id = "predicted.SingleR",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_UMAPplot_SingleR_NoLegend.pdf",sep = ""), width = 20, height = 20, units = "cm")

# Save objects ####
# We save a SingleCellExperimet object with all genes for differential expression analysis
saveRDS(data, file = "snRNA_SCE_CellAnnotation.rds")

# We save Seurat object for visualization
saveRDS(seurat.obj, file = "snRNA_Seurat_CellAnnotation.rds")

# We save reference with 100000 random cells from Visual Cortex and aggregated reference
save(aggregref, file = "Allen_aggregref.RData")
saveRDS(reference, file = "Allen_MMv21_VIS.rds")

# Cortex Markers ####
# UMAP plots of cell-type markers expression were visualized to see if the cell-type annotation made sense.

# Markers
markers <- c("Gfap","C1qa","Mbp","Sox9","Sox10","Glul","Mog","Camk2a","Gad1","Gad2",
             "Slc32a1","Sst","Nos1","Pvalb","Htr3a","Vip","Cck","Reln","Npy","Calb2",
             "Slc17a7","Nrgn","Foxp2","Otx1","Rorb","Etv1","Cux1","Tle4")

seurat.obj <- NormalizeData(seurat.obj)

Idents(object = seurat.obj) <-  seurat.obj$predicted.subclass_label

for (i in seq_along(markers)) {
  p <- FeaturePlot(seurat.obj, features = markers[i], reduction="UMAP", 
                   label=T, repel=T,label.size = 5) + labs(x = "UMAP1", y="UMAP2") +
    scale_color_continuous(low = "#fffbf7", high = "red") 
  ggsave(p, file=paste("snRNA_UMAPplot_",markers[i],".pdf",sep = ""),  width = 20, height = 20, units = "cm")
}

