# Cell-types annotation ####
# The goal of this analysis is to identify cell types,
# using an automantic and reference-based method, Azimuth.

# Set up
library(SingleCellExperiment)
library(scuttle)
library(HDF5Array)
library(dplyr)
library(Seurat)
library(Azimuth)
library(SeuratData)
library(ggplot2)
library(biomaRt)
library(SingleR)
library(heatmap)
library(AllenInstituteBrainData)

# Reference dataset ####
# We used Allen Whole Cortex & Hippocampus - 10x genomics (v 2021)
# as reference dataset for cell-type assignment.
reference <- AllenInstituteBrainData("Allen_Mouse_2021")
rownames(reference) <- rowData(reference)$X

# We selected Non-Neuronal, Neurons Glutamatergic and Neurons GABA-ergic,
# coming from Visual cortex (VIS, VISl, VISm, VISp) region
reference <- reference[, c(reference$region_label == "VIS" |
                             reference$region_label == "VISl" |
                             reference$region_label == "VISm" |
                             reference$region_label == "VISp")]
reference <- reference[, !is.na(reference$subclass_label) & reference$subclass_label != ""]
reference <- reference[, -which(grepl("ENT", reference$subclass_label))]
reference <- reference[, -which(grepl("PPP", reference$subclass_label))]
reference <- reference[, -which(grepl("CR", reference$subclass_label))]
reference <- reference[, -which(grepl("Meis", reference$subclass_label))]
reference <- reference[, -which(grepl("SUB", reference$subclass_label))]

# For computational issues, we decided to select 100,000 cortical random cells.
# First, all cell was selected from the cell types with less than 100 cells.
no_random <- reference[, c(reference$subclass_label == "SMC-Peri" | 
                             reference$subclass_label == "VLMC")]

# The resting cell types were selected.
reference <- reference[, !c(reference$subclass_label == "SMC-Peri" | 
                              reference$subclass_label == "VLMC")]

df <- data.frame(colData(reference))
# First, we randomly selected 100 cells for each cell type.
set.seed(23)
random <- df %>% group_by(subclass_label) %>% slice_sample(n = 100)
random100 <- reference[, colnames(reference) %in% random$sample_name]

reference <- reference[, !(colnames(reference) %in% random$sample_name)]

# Then, we randomly selected 98046 of the resting cells
set.seed(23)
reference <- reference[, sample(colnames(reference), 98046)]

reference <- cbind(reference, random100, no_random)

# Cell Annotation with Azimuth ####
# The reference dataset was converted into an Azimuth-compatible object
# for cell-type assignment.
# First, the reference dataset was converted into Seurat object
counts <- as.matrix(counts(reference))
coldata <- colData(reference)
reference_so <- CreateSeuratObject(counts = counts,
                                   meta.data = data.frame(coldata))

# And we created a compatible object for Azimuth cell annotation
reference_so <- SCTransform(reference_so, assay = "RNA", new.assay.name = "SCT",
                            variable.features.n = 2000, verbose = TRUE,
                            conserve.memory = TRUE)
reference_so <- RunPCA(reference_so, assay = "SCT", npcs = 50, verbose = FALSE,
                       reduction.name = "PCA", return.model = TRUE)
reference_so <- RunUMAP(reference_so, assay = "SCT", reduction = "PCA",
                        dims = seq_len(50), seed.use = 1, verbose = FALSE,
                        reduction.name = "umap", return.model = TRUE)

reference_so$subclass_label <- as.factor(reference_so$subclass_label)
Idents(object = reference_so) <- "subclass_label"

# Azimuth-compatible object
reference_azimuth <- AzimuthReference(reference_so, refUMAP = "umap",
                                      refDR = "PCA", refAssay = "SCT",
                                      dims = 1:50, metadata = c("subclass_label"),
                                      verbose = TRUE)

# save reference in a folder called "reference"
ref_dir <- "reference/"
SaveAnnoyIndex(object = reference_azimuth[["refdr.annoy.neighbors"]], 
               file = file.path(ref_dir, "idx.annoy"))
saveRDS(object = reference_azimuth, file = file.path(ref_dir, "ref.Rds"))

# Load single-nuclear dataset after pre-processing
scelist_sgl <- readRDS("snrna_scelist_sgl.rds")
# The six SingleCellExperiment were combined and converted into a Seurat object.
counts <- cbind(counts(scelist_sgl[[1]]), counts(scelist_sgl[[2]]),
                counts(scelist_sgl[[3]]), counts(scelist_sgl[[4]]), 
                counts(scelist_sgl[[5]]), counts(scelist_sgl[[6]]))

sample_id <- c(scelist_sgl[[1]]$sample_id, scelist_sgl[[2]]$sample_id, 
               scelist_sgl[[3]]$sample_id, scelist_sgl[[4]]$sample_id, 
               scelist_sgl[[5]]$sample_id, scelist_sgl[[6]]$sample_id)

condition <- c(scelist_sgl[[1]]$condition, scelist_sgl[[2]]$condition, 
               scelist_sgl[[3]]$condition, scelist_sgl[[4]]$condition, 
               scelist_sgl[[5]]$condition, scelist_sgl[[6]]$condition)

sce_obj <-  SingleCellExperiment(assays = list(counts = counts))
sce_obj$sample_id <- sample_id
sce_obj$condition <- condition
colnames(sce_obj) <- paste(colnames(sce_obj), sce_obj$sample_id, sep = "_")

# SCE object was converted into a Seurat object
seurat_obj <- CreateSeuratObject(counts = counts(sce_obj), 
                                 meta.data = data.frame(colData(sce_obj)))

# Cell-type annotation with Azimuth
seurat_obj <- RunAzimuth(seurat_obj, reference = "reference/")

# Dimensional riduction
seurat_obj <- RunTSNE(seurat_obj, reduction = "integrated_dr", 
                      dims = seq_len(20), seed.use = 1, do.fast = TRUE, 
                      verbose = FALSE, reduction.name = "TSNE")

seurat_obj <- RunUMAP(seurat_obj, reduction = "integrated_dr", 
                      dims = seq_len(20), seed.use = 1, verbose = FALSE,
                      reduction.name = "UMAP")

# Visualization
load("AllenColorLabel.RData")

p <- DimPlot(seurat_obj, reduction = "TSNE", group.by = "predicted.subclass_label") +
  NoLegend() + labs(x = "TSNE1", y = "TSNE2") + ggtitle("") +
  scale_color_manual(values = subclass_color)
p <- LabelClusters(p, id = "predicted.subclass_label",  fontface = "bold", color = "black", size = 5)
ggsave(p, file = "snRNA_TSNEplot_Azimuth_NoLegend.pdf", width = 20, height = 20, units = "cm")

p1 <- DimPlot(seurat.obj, reduction = "UMAP", group.by = "predicted.subclass_label")+
  NoLegend() + labs(x = "UMAP1", y="UMAP2") + ggtitle("") + 
  scale_color_manual(values=subclass.color)
p1 <- LabelClusters(p1, id = "predicted.subclass_label",  fontface = "bold", color = "black", size=2)
p1 <- p1 + theme(axis.text=element_text(size=7), axis.title=element_text(size=7))

ggsave("Peixoto_Figure4_part1.jpg",  width = 10, height = 10, units = "cm")
ggsave("Peixoto_Figure4_part1.pdf",  width = 10, height = 10, units = "cm")

# Save objects ####
# The SingleCellExperiment object with Azimuth labels was saved
sce_obj$azimuth_labels <- seurat_obj$predicted.subclass_label
saveRDS(sce_obj, file = "snrna_sce_annot.rds")

# The Seurat object with Azimuth was saved
saveRDS(seurat_obj, file = "snrna_seurat_annot.rds")

# The reference dataset, randomly selected was saved
saveRDS(reference, file = "Allen_MMv21_VIS.rds")
