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
library(pheatmap)

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
sce.obj <- readRDS("snRNA_DataPreprocessing.rds")
# The samples were merged and were converted into a Seurat object
counts <- cbind(counts(sce.obj[[1]]),counts(sce.obj[[2]]),counts(sce.obj[[3]]),counts(sce.obj[[4]]),counts(sce.obj[[5]]),counts(sce.obj[[6]]))

sample_id <- c(sce.obj[[1]]$sample_id, sce.obj[[2]]$sample_id, sce.obj[[3]]$sample_id, sce.obj[[4]]$sample_id, 
               sce.obj[[5]]$sample_id, sce.obj[[6]]$sample_id)

condition <- c(sce.obj[[1]]$condition, sce.obj[[2]]$condition, sce.obj[[3]]$condition, sce.obj[[4]]$condition, 
               sce.obj[[5]]$condition, sce.obj[[6]]$condition)

sce.obj <-  SingleCellExperiment(assays=list(counts=counts))
sce.obj$sample_id <- sample_id
sce.obj$condition <- condition
colnames(sce.obj) <- paste(colnames(sce.obj), sce.obj$sample_id, sep="_")

# Dataset was converted into a Seurat object
seurat.obj <- CreateSeuratObject(counts = counts(sce.obj), meta.data = data.frame(colData(sce.obj)), project = "SD")

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

# Normalized SCE data
sce.obj <- logNormCounts(sce.obj)

pred.SigleR <- SingleR(sce.obj, ref=aggregref, labels=aggregref$subclass_label)

# Save SingleR results into Seurat object for visualization
seurat.obj$predicted.SingleR <- pred.SingleR$labels

# Save cell-type annotation results into SingleCellExperiment object
sce.obj$Azimuth.labels <- seurat.obj$predicted.subclass_label
sce.obj$SingleR.labels <- pred.SigleR$labels

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
saveRDS(sce.obj, file = "snRNA_SCE_CellAnnotation.rds")

# We save Seurat object for visualization
saveRDS(seurat.obj, file = "snRNA_Seurat_CellAnnotation.rds")

# We save reference with 100000 random cells from Visual Cortex and aggregated reference
save(aggregref, file = "Allen_aggregref.RData")
saveRDS(reference, file = "Allen_MMv21_VIS.rds")

# Proportion of intron-containing reads ####
# We computed the proportion of intron-containing reads to identify cell-type with low proportion of pre-mRNA
# Load original single-nuclear RNA-seq Sleep Deprivation dataset
data <- readRDS("/home/zuin/sce_mouse_sleep_snrnaseq_complete.rds")
data$sample_id[data$sample_id=="8E"] <- "6E"
data$condition <- data$sample_id
data$condition[grep("C", data$condition)] <- "HC" # Home Cage (HC)
data$condition[grep("E", data$condition)] <- "SD" # Sleep Deprivated (SD)

data$sample_id <- paste(substr(data$sample_id,1,1), data$condition, sep = "")
colnames(data) <- paste(colnames(data), data$sample_id, sep = "_")

# Intersection between the columns of original and filtered dataset
inter <- intersect(colnames(data), colnames(sce.obj))
data <- data[,colnames(data) %in% inter]
data$Azimuth.labels <- sce.obj$Azimuth.labels

# Sum of gene expression for each column
sum.tot <- colSums(counts(data))

# Select only intron-containing reads
introns <- data[which(grepl("-I", rownames(data))),]
# Sum of intronic expression 
introns.sum <- colSums(counts(introns))
# Proportion of intron-containing reads for each column
prop.Intr <- introns.sum/sum.tot

# Violin plot
df <- as.data.frame(prop.Intr)
df$sample_id <- data$sample_id
df$CellType <- data$Azimuth.labels

p <- ggplot(df, aes(x=sample_id, y=prop.Intr)) + 
  geom_violin() + stat_summary(fun=median, geom="point", size=1, color="red") + # Add median 
  labs(x="Sample", y = "Proportion Intronic Reads")+theme_classic(base_size = 7)
ggsave(p, file="snRNA_ViolinPlot_ProportIntronicReads_Sample.jpg", width = 7, height = 7, units = "cm")
ggsave(p, file="snRNA_ViolinPlot_ProportIntronicReads_Sample.pdf", width = 7, height = 7, units = "cm")


for (i in 1:6) {
  df.sample <- df[df$sample_id==levels(factor(df$sample_id))[i],]
  p <- ggplot(df.sample, aes(x=CellType, y=prop.Intr)) + 
    geom_violin() + stat_summary(fun=median, geom="point", size=1, color="red") + # Add median
    labs(title=paste0(levels(factor(df$sample_id))[i]), x="Cell type", y = "Proportion Intronic Reads") + theme_classic(base_size = 7) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(p, file=paste("snRNA_ViolinPlot_ProportIntronicReads_Sample",i,".jpg",sep = ""), width = 7, height = 7, units = "cm")
  ggsave(p, file=paste("snRNA_ViolinPlot_ProportIntronicReads_Sample",i,".pdf",sep = ""), width = 7, height = 7, units = "cm")
}

# Heatmap of cell-type specific markers #####
markers <- read.table("markers.txt", header = TRUE)
markers<- markers[order(markers$Ensembl_ID),]

sce.markers <- sce[rownames(sce) %in% markers$Ensembl_ID,]
sce.markers <- sce.markers[order(rownames(sce.markers)),]
rownames(sce.markers) <- markers$Gene_Name

cell.type <- levels(factor(sce$Azimuth.labels))
average.celltype <- matrix(nrow = length(rownames(sce.markers)), ncol =length(levels(factor(sce.markers$Azimuth.labels))))
for (j in seq_along(cell.type)) {
  average.celltype[,j] <- rowMeans(logcounts(sce.markers[,sce.markers$Azimuth.labels==cell.type[j]]))
}
rownames(average.celltype) <- rownames(sce.markers)
colnames(average.celltype) <- cell.type

p.scale <- pheatmap(average.celltype, scale="row", fontsize = 7)
ggsave(p.scale, file=paste("Peixoto_Figure4_Supplement3_part1.jpg",sep = ""),  width = 10, height = 10, units = "cm")
ggsave(p.scale, file=paste("Peixoto_Figure4_Supplement3_part1.pdf",sep = ""),  width = 10, height = 10, units = "cm")

p <- pheatmap(average.celltype, fontsize = 7)
ggsave(p, file=paste("Peixoto_Figure4_Supplement3_part2.jpg",sep = ""),  width = 10, height = 10, units = "cm")
ggsave(p, file=paste("Peixoto_Figure4_Supplement3_part2.pdf",sep = ""),  width = 10, height = 10, units = "cm")
