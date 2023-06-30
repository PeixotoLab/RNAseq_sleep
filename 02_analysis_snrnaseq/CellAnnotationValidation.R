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

sce.obj <- readRDS("snRNA_SCE_CellAnnotation.rds")

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

sce.markers <- sce.obj[rownames(sce.obj) %in% markers$Ensembl_ID,]
sce.markers <- sce.markers[order(rownames(sce.markers)),]
rownames(sce.markers) <- markers$Gene_Name

cell.type <- levels(factor(sce.obj$Azimuth.labels))
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

# Best cell annotation methods #####
# To identify the best cell annotation method, we used two datasets from the BRAIN Initiative Cell Census atlas for which cells were already annotated
# 10x Nuclei v3 Broad ####
# Single-nuclei RNA-seq data from Primary Motor Area 
data <- readRDS("snRNA_v3Broad.rds")
data$sample_name <- substring(data$sample_name, 1,27)

data <- logNormCounts(data)

# Cell-type Annotation 
# Azimuth 
counts <- as.matrix(counts(data))
macosko.broad <- CreateSeuratObject(counts = counts, meta.data = data.frame(colData(data)))
macosko.broad <- RunAzimuth(macosko.broad, reference = "mousecortexref")

macosko.broad <- RunTSNE(macosko.broad, reduction = "integrated_dr", dims = seq_len(20), seed.use = 1, do.fast = TRUE, verbose = FALSE,
                         reduction.name = "TSNE")

macosko.broad <- RunUMAP(macosko.broad, reduction = "integrated_dr", dims = seq_len(20),seed.use = 1, verbose = FALSE, reduction.name = "UMAP")

# SingleR 
# We used Allen Whole Cortex & Hippocampus - 10x genomics (v 2021) as reference dataset for cell annotation
reference <- loadHDF5SummarizedExperiment(dir="/mnt/callisto/Zuin", prefix="Allen_mm_21")
reference <- as(reference, "SingleCellExperiment")
names(assays(reference)) <- c("counts")

# We selected Primary Motor Area 
reference <- reference[,c(reference$region_label=="MOp")]
reference <- reference[,-which(grepl("ENT", reference$subclass_label))]
reference <- reference[,-which(grepl("PPP", reference$subclass_label))]

aggregref <- aggregateAcrossCells(reference, use.assay.type = "counts",id=DataFrame(label=reference$subclass_label))
aggregref <- logNormCounts(aggregref)

save(aggregref, file = "AllenMOp_aggreg.rds")

pred.SingleR <- SingleR(data, ref=aggregref,labels=aggregref$subclass_label)
macosko.broad$SingleR <- pred.SingleR$labels

# save cell-type annotation results into a Seurat object
saveRDS(macosko.broad, file = "snRNA_v3Broad_Seurat.rds")

# Note: for computational issues SingleR on MNN method wasn't made.

# Visualization 
# Palette color
label.color <- c("#957b46","#5100FF","#c95f3f","#0BE652","#0D5B78","#50B2AD","#3E9E64",
                 "#2D8CB8","#A19922","#7044AA","#DA808C","#94AF97","#744700","#cd8109","#33632c","#5d554c",
                 "#D93137","#4d7647","#ffff00","#FF9900","#B864CC","#54a04d")
names(label.color) <- levels(factor(macosko.broad$Allen.subclass_label))

Azimuth.color <- c("#957b46","#c95f3f","#0BE652",
                   "#0D5B78","#50B2AD","#3E9E64",
                   "#2D8CB8","#A19922","#5100FF",
                   "#7044AA","#DA808C","#33632c",
                   "#94AF97",
                   "#744700","#cd8109","#4d7647",
                   "#D93137","#ffff00","#FF9900",
                   "#c5ce00","#B864CC","#54a04d")

names(Azimuth.color) <- levels(factor(macosko.broad$predicted.subclass))

Allen.color <- data.frame(aggregref$subclass_label, aggregref$subclass_color)
colnames(Allen.color) <- c("Label", "color")

# change Astro color
Allen.color$color[Allen.color$color=="#665C47"] <- "#957b46"
# change Oligo color
Allen.color$color[Allen.color$color=="#53776C"] <- "#744700"
# change SMC-Peri color
Allen.color$color[Allen.color$color=="#807059"] <- "#4d7647"
# change VLMC color
Allen.color$color[Allen.color$color=="#697255"] <- "#54a04d"
# change Endo color
Allen.color$color[Allen.color$color=="#8D6C62"] <- "#c95f3f"
# change Sncg color 
Allen.color$color[Allen.color$color=="#D3408D"] <- "#ffff00"

SingleR.color <- Allen.color$color
names(SingleR.color) <- Allen.color$Label

save(label.color, Azimuth.color, SingleR.color, file="MacoskoBroad_color.RData")

plot.so <- DimPlot(macosko.broad, reduction = "TSNE", group.by = "Allen.subclass_label") +
  scale_color_manual(values=label.color) + NoLegend()+ labs(x = "TSNE1", y="TSNE2") + ggtitle("") 
plot.so <- LabelClusters(plot.so, id = "Allen.subclass_label",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_atBroad_TSNEplot.pdf",sep = ""), width = 20, height = 20, units = "cm")

plot.so <- DimPlot(macosko.broad, reduction = "UMAP", group.by = "Allen.subclass_label") +
  scale_color_manual(values=label.color) + NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") 
plot.so <- LabelClusters(plot.so, id = "Allen.subclass_label",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_atBroad_UMAPplot.pdf",sep = ""), width = 20, height = 20, units = "cm")

plot.so <- DimPlot(macosko.broad, reduction = "TSNE", group.by = "predicted.subclass") +
  scale_color_manual(values=Azimuth.color) + NoLegend()+ labs(x = "TSNE1", y="TSNE2") + ggtitle("") 
plot.so <- LabelClusters(plot.so, id = "predicted.subclass",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_atBroad_TSNEplot_Azimuth.pdf",sep = ""), width = 20, height = 20, units = "cm")

plot.so <- DimPlot(macosko.broad, reduction = "UMAP", group.by = "predicted.subclass") +
  scale_color_manual(values=Azimuth.color) + NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") 
plot.so <- LabelClusters(plot.so, id = "predicted.subclass",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_atBroad_UMAPplot_Azimuth.pdf",sep = ""), width = 20, height = 20, units = "cm")

plot.so <- DimPlot(macosko.broad, reduction = "TSNE", group.by = "SingleR") + 
  scale_color_manual(values=SingleR.color)+ NoLegend()+ labs(x = "TSNE1", y="TSNE2") + ggtitle("") 
plot.so <- LabelClusters(plot.so, id = "SingleR",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_atBroad_TSNEplot_SingleR.pdf",sep = ""), width = 20, height = 20, units = "cm")

plot.so <- DimPlot(macosko.broad, reduction = "UMAP", group.by = "SingleR") +
  scale_color_manual(values=SingleR.color) + NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") 
plot.so <- LabelClusters(plot.so, id = "SingleR",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_atBroad_UMAPplot_Azimuth_SingleR.pdf",sep = ""), width = 20, height = 20, units = "cm")

# 10x Nuclei v2 AIBS ####
data.v2 <- readRDS("snRNA_v2AIBS.RDS")
data.v2 <- logNormCounts(data.v2)
colnames(data.v2) <- data.v2$X

# Cell Annotation 
## Azimuth 
macosko.nucleusv2 <- CreateSeuratObject(counts = counts(data.v2), meta.data = data.frame(colData(data.v2)))
macosko.nucleusv2 <- RunAzimuth(macosko.nucleusv2, reference = "mousecortexref")

macosko.nucleusv2 <- RunTSNE(macosko.nucleusv2, reduction = "integrated_dr", dims = seq_len(20), seed.use = 1, do.fast = TRUE, verbose = FALSE,
                             reduction.name = "TSNE")
macosko.nucleusv2 <- RunUMAP(macosko.nucleusv2, reduction = "integrated_dr", dims = seq_len(20), seed.use = 1, verbose = FALSE, reduction.name = "UMAP")

## SingleR 
pred.SingleR <- SingleR(data.v2, ref=aggregref,labels=aggregref$subclass_label)
macosko.nucleusv2$SingleR <- pred.SingleR$labels

## SingleR on MNN 
# In this case we corrected data from batch effect with Mutual Nearest Neighbor (MNN) method
MNN.correction <- fastMNN(data.v2, batch=data.v2$Donor, d=50, k=20, subset.row=rownames(data.v2), BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
pred.MNN <- SingleR(MNN.correction, ref=aggregref,labels=aggregref$subclass_label,assay.type.test = "reconstructed")

macosko.nucleusv2$SingleR_MNN <- pred.MNN$labels

saveRDS(macosko.nucleusv2, file = "snRNA_NucleusV2AIBS_Seurat.rds")

# Visualization 
# Palette color
label.color <- c("#957b46","#c95f3f","#0BE652","#0D5B78","#50B2AD","#3E9E64",
                 "#2D8CB8","#A19922","#5100FF","#7044AA","#DA808C",
                 "#94AF97","#744700","#000000","#D93137","#ffff00","#FF9900","#B864CC")
names(label.color) <- levels(factor(macosko.nucleusv2$subclass_label))

Azimuth.color <- c("#957b46","#c95f3f","#0BE652",
                   "#0D5B78","#50B2AD","#3E9E64",
                   "#2D8CB8","#A19922","#5100FF",
                   "#7044AA","#DA808C",
                   "#94AF97","#744700","#000000",
                   "#4d7647","#D93137","#ffff00",
                   "#FF9900","#B1B10C", "#B864CC","#54a04d")

names(Azimuth.color) <- levels(factor(macosko.nucleusv2$predicted.subclass))

Allen.color <- data.frame(aggregref$subclass_label, aggregref$subclass_color)
colnames(Allen.color) <- c("Label", "color")

# change Astro color
Allen.color$color[Allen.color$color=="#665C47"] <- "#957b46"
# change Oligo color
Allen.color$color[Allen.color$color=="#53776C"] <- "#744700"
# change SMC-Peri color
Allen.color$color[Allen.color$color=="#807059"] <- "#4d7647"
# change VLMC color
Allen.color$color[Allen.color$color=="#697255"] <- "#54a04d"
# change Endo color
Allen.color$color[Allen.color$color=="#8D6C62"] <- "#c95f3f"
# change Sncg color 
Allen.color$color[Allen.color$color=="#D3408D"] <- "#ffff00"

SingleR.color <- Allen.color$color
names(SingleR.color) <- Allen.color$Label

save(label.color, Azimuth.color, SingleR.color, file = "MacoskoNucleusV2_color.RData")

plot.so <- DimPlot(macosko.nucleusv2, reduction = "TSNE", group.by = "subclass_label") +
  scale_color_manual(values=label.color) + NoLegend()+ labs(x = "TSNE1", y="TSNE2") + ggtitle("") 
plot.so <- LabelClusters(plot.so, id = "subclass_label",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_nucleus_v2_TSNEplot.pdf",sep = ""), width = 20, height = 20, units = "cm")

plot.so <- DimPlot(macosko.nucleusv2, reduction = "UMAP", group.by = "subclass_label") +
  scale_color_manual(values=label.color) + NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") 
plot.so <- LabelClusters(plot.so, id = "subclass_label",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_nucleus_v2_UMAPplot.pdf",sep = ""), width = 20, height = 20, units = "cm")

plot.so <- DimPlot(macosko.nucleusv2, reduction = "TSNE", group.by = "predicted.subclass") +
  scale_color_manual(values=Azimuth.color) + NoLegend()+ labs(x = "TSNE1", y="TSNE2") + ggtitle("") 
plot.so <- LabelClusters(plot.so, id = "predicted.subclass",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_nucleus_v2_TSNEplot_Azimuth.pdf",sep = ""), width = 25, height = 25, units = "cm")

plot.so <- DimPlot(macosko.nucleusv2, reduction = "UMAP", group.by = "predicted.subclass") +
  scale_color_manual(values=Azimuth.color) + NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") 
plot.so <- LabelClusters(plot.so, id = "predicted.subclass",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_nucleus_v2_UMAPplot_Azimuth.pdf",sep = ""), width = 20, height = 20, units = "cm")

plot.so <- DimPlot(macosko.nucleusv2, reduction = "TSNE", group.by = "SingleR") + 
  scale_color_manual(values=SingleR.color)+ NoLegend()+ labs(x = "TSNE1", y="TSNE2") + ggtitle("") 
plot.so <- LabelClusters(plot.so, id = "SingleR",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_nucleus_v2_TSNEplot_SingleR.pdf",sep = ""), width = 20, height = 20, units = "cm")

plot.so <- DimPlot(macosko.nucleusv2, reduction = "UMAP", group.by = "SingleR") +
  scale_color_manual(values=SingleR.color) + NoLegend()+ labs(x = "UMAP1", y="UMAP2") + ggtitle("") 
plot.so <- LabelClusters(plot.so, id = "SingleR",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_nucleus_v2_UMAPplot_Azimuth.pdf",sep = ""), width = 20, height = 20, units = "cm")

plot.so <- DimPlot(macosko.nucleusv2, reduction = "TSNE", group.by = "SingleR_MNN") + 
  labs(x = "TSNE1", y="TSNE2")+ ggtitle("")+ NoLegend()+scale_color_manual(values=SingleR.color)
plot.so <- LabelClusters(plot.so, id = "SingleR_MNN",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_nucleus_v2_TSNEplot_SingleR_MNN.pdf",sep = ""), width = 20, height = 20, units = "cm")

plot.so <- DimPlot(macosko.nucleusv2, reduction = "UMAP", group.by = "SingleR_MNN") + 
  labs(x = "UMAP1", y="UMAP2")+ ggtitle("")+ NoLegend()+scale_color_manual(values=SingleR.color)
plot.so <- LabelClusters(plot.so, id = "SingleR_MNN",  fontface = "bold", color = "black", size=5)
ggsave(plot.so, file=paste("snRNA_nucleus_v2_UMAPplot_SingleR_MNN.pdf",sep = ""), width = 20, height = 20, units = "cm")

