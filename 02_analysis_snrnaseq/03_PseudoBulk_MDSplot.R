# Pseudo - bulk ####
# The goal in this analysis is to create pseudo-bulk

# Set up
library(SingleCellExperiment)
library(scuttle)
library(muscat)
library(edgeR)
library(ggplot2)
library(patchwork)

# Load SCE object
sce_obj <- readRDS("snrna_sce_annot.rds")
# The cell-types with more than 500 cells were selected
sce_obj <- sce_obj[, !c(sce_obj$azimuth_labels == "Car3" |
                          sce_obj$azimuth_labels == "Endo" |
                          sce_obj$azimuth_labels == "Lamp5" |
                          sce_obj$azimuth_labels == "Sncg" |
                          sce_obj$azimuth_labels == "Sst Chodl" |
                          sce_obj$azimuth_labels == "VLMC")]

# Note: if scuttle::logNormCounts doesn't work, remake SingleCellExpression object.

# Pseudo-bulk creation
snrna_pb <- aggregateAcrossCells(sce_obj, use.assay.type = "counts", 
                                 id = DataFrame(label = sce_obj$azimuth_labels, 
                                                sample = sce_obj$sample_id))
colnames(snrna_pb) <- paste(snrna_pb$azimuth_labels, snrna_pb$sample_id, sep = "_")
snrna_pb <- logNormCounts(snrna_pb)

# Cell-type class variable was created
snrna_pb$class <- snrna_pb$azimuth_labels
snrna_pb$class[which(grepl("CTX", snrna_pb$class))] <- "Glutamatergic"

snrna_pb$class[which(grepl("-", snrna_pb$class))] <- "Other"

snrna_pb$class[!(snrna_pb$class == "Astro" | snrna_pb$class == "Oligo" |
                   snrna_pb$class == "Glutamatergic" |
                   snrna_pb$class == "Other")] <- "GABAergic"

# Pseudo-bulk object was saved
saveRDS(snrna_pb, file = "snrna_pb.rds")

# MDS plot ####
# Load Allen color palette
load("AllenColorLabel.RData")

# Load the code inside pbMDSfunction.R file if this error appears:
# "Error in element_line(linewidth = 0.2, color = "lightgrey") :unused argument (linewidth = 0.2)"
source("pbMDSfunction.R")

neuronal_color <- subclass_color[-c(1:3, 12:14, 16, 17, 19, 21)]

# MDS plot according to the Glutamatergic labels
prep_sce <- prepSCE(snrna_pb[, snrna_pb$class == "Glutamatergic"],
                    kid = "azimuth_labels", # subpopulation assignments
                    gid = "condition",  # group IDs
                    sid = "sample_id",   # sample_id IDs
                    drop = TRUE)
pb <-  aggregateData(prep_sce, assay = "logcounts",
                     by = c("cluster_id", "sample_id"))

mds_gluta <- pbMDS(pb) + scale_color_manual(values = neuronal_color) +
  labs(col = "Cell-types", shape = "Condition") +
  theme(axis.text = element_text(size = 7), axis.title = element_text(size = 7),
        legend.position = "none")

# MDS plot according to the GABA-ergic labels
prep_sce <- prepSCE(snrna_pb[, snrna_pb$class == "GABAergic"],
                    kid = "azimuth_labels", gid = "condition",
                    sid = "sample_id", drop = TRUE)
pb <-  aggregateData(prep_sce, assay = "logcounts",
                     by = c("cluster_id", "sample_id"))

mds_gaba <- pbMDS(pb) + scale_color_manual(values = neuronal_color) +
  labs(col = "Cell-types", shape = "Condition") +
  theme(axis.text = element_text(size = 7), axis.title = element_text(size = 7),
        legend.position = "none")

patched <- (mds_gluta | mds_gaba)
ggsave("Peixoto_Figure4_part2.jpg",  width = 18, height = 9, units = "cm")
ggsave("Peixoto_Figure4_part2.pdf",  width = 18, height = 9, units = "cm")

# MDS plot according to the Non-Neuronal labels
prep_sce <- prepSCE(snrna_pb[, (snrna_pb$class == "Astro" |
                                  snrna_pb$class == "Oligo" |
                                  snrna_pb$class == "Other")],
                    kid = "azimuth_labels", gid = "condition",
                    sid = "sample_id", drop = TRUE)
pb <-  aggregateData(prep_sce, assay = "logcounts", 
                     by = c("cluster_id", "sample_id"))

mds_other <- pbMDS(pb) + scale_color_manual(values = subclass_color) +
  labs(col = "Cell-types", shape = "Condition") +
  theme(legend.spacing.y = unit(0, "cm"), axis.text = element_text(size = 7),
        axis.title = element_text(size = 7), legend.title = element_text(size = 7),
        legend.text = element_text(size = 7), legend.key.size = unit(0.4, "cm"))

# MDS plot according to the neuronal labels on negative control genes
# Load negative control genes
negctrl <- read.table("SD_Negative_Controls.txt")
# Negative control gene was selected inside the pseudo-bulk object
pb_negctrl <- snrna_pb[rownames(snrna_pb) %in% negctrl$x, ]

prep_sce <- prepSCE(pb_negctrl[, !(pb_negctrl$class == "Astro" |
                                     pb_negctrl$class == "Oligo" |
                                     pb_negctrl$class == "Other")],
                    kid = "azimuth_labels", gid = "condition",
                    sid = "sample", drop = TRUE)
pb <-  aggregateData(prep_sce, assay = "logcounts", 
                     by = c("cluster_id", "sample_id"))

mds_negctrl <- pbMDS(pb) + scale_color_manual(values = subclass_color) +
  labs(col = "Cell-types", shape = "Condition") +
  theme(legend.spacing.y = unit(0, "cm"), axis.text = element_text(size = 7),
        axis.title = element_text(size = 7), legend.title = element_text(size = 7),
        legend.text = element_text(size = 7), legend.key.size = unit(0.4, "cm"))

patched <- (mds_other | mds_negctrl) +  plot_layout(guides = "collect")
ggsave("Peixoto_Figure4_Supplement4.jpg",  width = 20, height = 10, units = "cm")
ggsave("Peixoto_Figure4_Supplement4.pdf",  width = 20, height = 10, units = "cm")
