library(scuttle)
library(HDF5Array)
library(dplyr)
library(Seurat)
library(Azimuth)
library(SeuratData)
library(ggplot2)
library(biomaRt)
library(SingleR)
library(pheatmap)

# Heatmap ####
# To evaluate the cell-type assignments, we visualized cell-type specific
# markers with a heatmap of the log-normalized count average in each group.
# Load single-nuclear dataset
sce_obj <- readRDS("snRNA_SCE_CellAnnotation.rds")
# Load cell-type specific markers
markers <- read.table("markers.txt", header = TRUE)
markers <- markers[order(markers$Ensembl_ID), ]
# Cell-type specific markers were selected inside the dataset
sce_markers <- sce_obj[rownames(sce_obj) %in% markers$Ensembl_ID, ]
sce_markers <- sce_markers[order(rownames(sce_markers)), ]
rownames(sce_markers) <- markers$Gene_Name

# The log-normalized count average was computed in each group
label <- levels(factor(sce_obj$Azimuth.labels))
average_celltype <- matrix(nrow = length(rownames(sce_markers)), ncol = length(levels(factor(sce_markers$Azimuth.labels))))
for (j in seq_along(label)) {
  average_celltype[, j] <- rowMeans(logcounts(sce_markers[, sce_markers$Azimuth.labels == label[j]]))
}
rownames(average_celltype) <- rownames(sce_markers)
colnames(average_celltype) <- label

p <- pheatmap(average_celltype, fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE)
ggsave("Peixoto_Figure4_Supplement3_part2.jpg",  width = 10, height = 10, units = "cm")
ggsave("Peixoto_Figure4_Supplement3_part2.pdf",  width = 10, height = 10, units = "cm")
# Heatmap scaled by row
p_scale <- pheatmap(average_celltype, scale = "row", fontsize = 7, cluster_rows = FALSE, cluster_cols = FALSE)
ggsave("Peixoto_Figure4_Supplement3_part1.jpg", width = 10, height = 10, units = "cm")
ggsave("Peixoto_Figure4_Supplement3_part1.pdf", width = 10, height = 10, units = "cm")

# Proportion of intron-containing reads ####
# As an additional quality control, we checked if there were cell types
# with a low proportion of intronic reads, as this could be
# a sign of cytoplasmic RNA and assigned incorrectly to nuclei.
# Load original single-nuclear dataset
data <- readRDS("/home/zuin/sce_mouse_sleep_snrnaseq_complete.rds")
data$sample_id[data$sample_id == "8E"] <- "6E"
data$condition <- data$sample_id
data$condition[grep("C", data$condition)] <- "HC" # Home Cage (HC)
data$condition[grep("E", data$condition)] <- "SD" # Sleep Deprivated (SD)
data$sample_id <- paste(substr(data$sample_id, 1, 1), data$condition, sep = "")
colnames(data) <- paste(colnames(data), data$sample_id, sep = "_")

# First, the low and damaged droplets and the potential doublets were removed
# from the original dataset.
inter <- intersect(colnames(data), colnames(sce_obj))
data <- data[, colnames(data) %in% inter]
data$Azimuth.labels <- sce_obj$Azimuth.labels

# The total sum of gene expression was computed for each column/sample.
tot_sum <- colSums(counts(data))
# The intron-containing reads were selected
introns <- data[which(grepl("-I", rownames(data))), ]
# The intronic sum was computed
intronic_sum <- colSums(counts(introns))
# The intronic proportion was computed for each sample
intronic_prop <- intronic_sum / tot_sum

## Violin plot ####
# Create dataframe
df <- as.data.frame(intronic_prop)
df$sample_id <- data$sample_id
df$CellType <- data$Azimuth.labels

violin_plot_sample <- ggplot(df, aes(x = sample_id, y = intronic_prop)) +
  # Add mean
  geom_violin() + stat_summary(fun = median, geom = "point", size = 1, color = "red") +
  labs(x = "Sample", y = "Proportion Intronic Reads") + theme_classic(base_size = 7)
ggsave("snRNA_ViolinPlot_ProportIntronicReads_Sample.jpg", width = 7, height = 7, units = "cm")
ggsave("snRNA_ViolinPlot_ProportIntronicReads_Sample.pdf", width = 7, height = 7, units = "cm")

# Violin plot for each cell-type in each sample
for (i in 1:6) {
  df_sample <- df[df$sample_id == levels(factor(df$sample_id))[i], ]
  
  violin_plot <- ggplot(df_sample, aes(x = CellType, y = intronic_prop)) +
    geom_violin() + stat_summary(fun = median, geom = "point", size = 1, color = "red") +
    labs(title = paste0(levels(factor(df$sample_id))[i]), x = "Cell type", y = "Proportion Intronic Reads") +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(violin_plot, file = paste("snRNA_ViolinPlot_ProportIntronicReads_Sample", i, ".jpg", sep = ""), width = 7, height = 7, units = "cm")
  ggsave(violin_plot, file = paste("snRNA_ViolinPlot_ProportIntronicReads_Sample", i, ".pdf", sep = ""), width = 7, height = 7, units = "cm")
}

# Identification of the best cell type assignment method ####
# To identify the best cell annotation method, we used two datasets from
# the BRAIN Initiative Cell Census atlas for which cells were already annotated.
# Two automatic and reference-based annotation methods, Azimuth and SinglR, were used.

## 10X v3 Broad ####
# The first is a single-nuclear dataset, obtained with a technique "10X v3".
# The primary motor cortex (MOp) was dissected in each mouse.
sn10x_v3 <- readRDS("snRNA_v3Broad.rds")
sn10x_v3$sample_name <- substring(sn10x_v3$sample_name, 1, 27)
sn10x_v3 <- logNormCounts(sn10x_v3)

### Azimuth ####
# SCE object was converted into a Seurat object
counts <- as.matrix(counts(sn10x_v3))
sn10x_v3_so <- CreateSeuratObject(counts = counts, meta.data = data.frame(colData(sn10x_v3)))
# We used the "mousecortexref" dataset as a reference.
# This is coming from the MOp region and is found inside the Azimuth function.
sn10x_v3_so <- RunAzimuth(sn10x_v3_so, reference = "mousecortexref")

sn10x_v3_so <- RunTSNE(sn10x_v3_so, reduction = "integrated_dr", seed.use = 1,
                       dims = seq_len(20), do.fast = TRUE, verbose = FALSE, reduction.name = "TSNE")

sn10x_v3_so <- RunUMAP(sn10x_v3_so, reduction = "integrated_dr", seed.use = 1,
                       dims = seq_len(20), verbose = FALSE, reduction.name = "UMAP")

### SingleR ####
# For SingleR, we used Allen Whole Cortex & Hippocampus - 10x genomics (v 2021)
# as reference dataset.
reference <- loadHDF5SummarizedExperiment(dir = "/mnt/callisto/Zuin", prefix = "Allen_mm_21")
reference <- as(reference, "SingleCellExperiment")
names(assays(reference)) <- c("counts")

# We selected cells from Primary Motor Area
reference <- reference[, c(reference$region_label == "MOp")]
reference <- reference[, -which(grepl("ENT", reference$subclass_label))]
reference <- reference[, -which(grepl("PPP", reference$subclass_label))]
# The reference was aggregated
aggregref <- aggregateAcrossCells(reference, use.assay.type = "counts", id = DataFrame(label = reference$subclass_label))
aggregref <- logNormCounts(aggregref)
save(aggregref, file = "AllenMOp_aggreg.rds")

pred_SingleR <- SingleR(sn10x_v3, ref = aggregref, labels = aggregref$subclass_label)
sn10x_v3_so$SingleR <- pred_SingleR$labels

# Seurat object was saved into a .rds file
saveRDS(sn10x_v3_so, file = "snRNA_v3Broad_Seurat.rds")

# Note: For computational issues (R crashed), SingleR on MNN method wasn't made.

### Cell-type assignment visualization ####
# Color palettes creation
label_color <- c("#957b46", "#5100FF", "#c95f3f", "#0BE652",
                 "#0D5B78", "#50B2AD", "#3E9E64", "#2D8CB8",
                 "#A19922", "#7044AA", "#DA808C", "#94AF97",
                 "#744700", "#cd8109", "#33632c", "#5d554c",
                 "#D93137", "#4d7647", "#ffff00", "#FF9900",
                 "#B864CC", "#54a04d")
names(label_color) <- levels(factor(sn10x_v3_so$Allen.subclass_label))

azimuth_color <- c("#957b46", "#c95f3f", "#0BE652",
                   "#0D5B78", "#50B2AD", "#3E9E64",
                   "#2D8CB8", "#A19922", "#5100FF",
                   "#7044AA", "#DA808C", "#33632c",
                   "#94AF97", "#744700", "#cd8109",
                   "#4d7647", "#D93137", "#ffff00",
                   "#FF9900", "#c5ce00", "#B864CC",
                   "#54a04d")
names(azimuth_color) <- levels(factor(sn10x_v3_so$predicted.subclass))

allen_color <- sn10x_v3.frame(aggregref$subclass_label, aggregref$subclass_color)
colnames(allen_color) <- c("Label", "color")
# some changes
allen_color$color[allen_color$color == "#665C47"] <- "#957b46"
allen_color$color[allen_color$color == "#53776C"] <- "#744700"
allen_color$color[allen_color$color == "#807059"] <- "#4d7647"
allen_color$color[allen_color$color == "#697255"] <- "#54a04d"
allen_color$color[allen_color$color == "#8D6C62"] <- "#c95f3f"
allen_color$color[allen_color$color == "#D3408D"] <- "#ffff00"

singler_color <- allen_color$color
names(singler_color) <- allen_color$Label


p <- DimPlot(sn10x_v3_so, reduction = "UMAP", group.by = "Allen.subclass_label") +
  scale_color_manual(values = label_color) + NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") + ggtitle("")
p <- LabelClusters(p, id = "Allen.subclass_label",  fontface = "bold", color = "black", size = 5)
ggsave(p, file = "snRNA_atBroad_UMAPplot.pdf", width = 20, height = 20, units = "cm")

p <- DimPlot(sn10x_v3_so, reduction = "UMAP", group.by = "predicted.subclass") +
  scale_color_manual(values = azimuth_color) + NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") + ggtitle("")
p <- LabelClusters(p, id = "predicted.subclass",  fontface = "bold", color = "black", size = 5)
ggsave(p, file = "snRNA_atBroad_UMAPplot_Azimuth.pdf", width = 20, height = 20, units = "cm")

p <- DimPlot(sn10x_v3_so, reduction = "UMAP", group.by = "SingleR") +
  scale_color_manual(values = singler_color) + NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") + ggtitle("")
p <- LabelClusters(p, id = "SingleR",  fontface = "bold", color = "black", size = 5)
ggsave(p, file = "snRNA_atBroad_UMAPplot_Azimuth_SingleR.pdf", width = 20, height = 20, units = "cm")

## 10X v2 AIBS ####
# The second is a single-nuclear dataset, obtained with a technique "10X v2".
# The primary motor cortex (MOp) was dissected in each mouse.
sn10x_v2 <- readRDS("snRNA_v2AIBS.RDS")
sn10x_v2 <- logNormCounts(sn10x_v2)
colnames(sn10x_v2) <- sn10x_v2$X

### Azimuth ####
# SCE object was converted into a Seurat object
sn10x_v2_so <- CreateSeuratObject(counts = counts(sn10x_v2), meta.data = data.frame(colData(sn10x_v2)))

sn10x_v2_so <- RunAzimuth(sn10x_v2_so, reference = "mousecortexref")

sn10x_v2_so <- RunTSNE(sn10x_v2_so, reduction = "integrated_dr", seed.use = 1,
                       dims = seq_len(20), do.fast = TRUE, verbose = FALSE, reduction.name = "TSNE")
sn10x_v2_so <- RunUMAP(sn10x_v2_so, reduction = "integrated_dr", seed.use = 1,
                       dims = seq_len(20), verbose = FALSE, reduction.name = "UMAP")

### SingleR ####
pred_SingleR <- SingleR(sn10x_v2, ref = aggregref, labels = aggregref$subclass_label)
sn10x_v2_so$SingleR <- pred_SingleR$labels

### SingleR on MNN ####
# Then SingleR on corrected data with the Mutual Nearest Neighbor (MNN) method
# from the batch effects was computed.
mnn <- fastMNN(sn10x_v2, batch = sn10x_v2$Donor, d = 50, k = 20, subset.row = rownames(sn10x_v2),
               BSPARAM = BiocSingular::RandomParam(deferred = TRUE))
# Annotation on reconstructed data
pred_mnn <- SingleR(mnn, ref = aggregref, labels = aggregref$subclass_label,
                    assay.type.test = "reconstructed")
sn10x_v2_so$SingleR_MNN <- pred_mnn$labels

saveRDS(sn10x_v2_so, file = "snRNA_NucleusV2AIBS_Seurat.rds")

### Cell-type assignment visualization ####
# Color palettes creation
label_color <- c("#957b46", "#c95f3f", "#0BE652", "#0D5B78",
                 "#50B2AD", "#3E9E64", "#2D8CB8", "#A19922",
                 "#5100FF", "#7044AA", "#DA808C", "#94AF97",
                 "#744700", "#000000", "#D93137", "#ffff00",
                 "#FF9900", "#B864CC")
names(label_color) <- levels(factor(sn10x_v2_so$subclass_label))

azimuth_color <- c("#957b46", "#c95f3f", "#0BE652",
                   "#0D5B78", "#50B2AD", "#3E9E64",
                   "#2D8CB8", "#A19922", "#5100FF",
                   "#7044AA", "#DA808C", "#94AF97",
                   "#744700", "#000000", "#4d7647",
                   "#D93137", "#ffff00", "#FF9900",
                   "#B1B10C", "#B864CC", "#54a04d")
names(azimuth_color) <- levels(factor(sn10x_v2_so$predicted.subclass))

allen_color <- sn10x_v3.frame(aggregref$subclass_label, aggregref$subclass_color)
colnames(allen_color) <- c("Label", "color")
# some changes
allen_color$color[allen_color$color == "#665C47"] <- "#957b46"
allen_color$color[allen_color$color == "#53776C"] <- "#744700"
allen_color$color[allen_color$color == "#807059"] <- "#4d7647"
allen_color$color[allen_color$color == "#697255"] <- "#54a04d"
allen_color$color[allen_color$color == "#8D6C62"] <- "#c95f3f"
allen_color$color[allen_color$color == "#D3408D"] <- "#ffff00"

singler_color <- allen_color$color
names(singler_color) <- allen_color$Label

p <- DimPlot(sn10x_v2_so, reduction = "UMAP", group.by = "subclass_label") +
  scale_color_manual(values = label_color) + NoLegend() +
  labs(x = "UMAP1", y ="UMAP2") + ggtitle("")
p <- LabelClusters(p, id = "subclass_label",  fontface = "bold", color = "black", size = 5)
ggsave(p, file = "snRNA_nucleus_v2_UMAPplot.pdf", width = 20, height = 20, units = "cm")

p <- DimPlot(sn10x_v2_so, reduction = "UMAP", group.by = "predicted.subclass") +
  scale_color_manual(values = azimuth_color) + NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") + ggtitle("") 
p <- LabelClusters(p, id = "predicted.subclass",  fontface = "bold", color = "black", size = 5)
ggsave(p, file = "snRNA_nucleus_v2_UMAPplot_Azimuth.pdf", width = 20, height = 20, units = "cm")

p <- DimPlot(sn10x_v2_so, reduction = "UMAP", group.by = "SingleR") +
  scale_color_manual(values = singler_color) + NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") + ggtitle("") 
p <- LabelClusters(p, id = "SingleR",  fontface = "bold", color = "black", size = 5)
ggsave(p, file = "snRNA_nucleus_v2_UMAPplot_Azimuth.pdf", width = 20, height = 20, units = "cm")

p <- DimPlot(sn10x_v2_so, reduction = "UMAP", group.by = "SingleR_MNN") + 
  labs(x = "UMAP1", y = "UMAP2") + ggtitle("") + NoLegend() +
  scale_color_manual(values = singler_color)
p <- LabelClusters(p, id = "SingleR_MNN",  fontface = "bold", color = "black", size = 5)
ggsave(p, file = "snRNA_nucleus_v2_UMAPplot_SingleR_MNN.pdf", width = 20, height = 20, units = "cm")
