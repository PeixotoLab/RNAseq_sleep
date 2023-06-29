library(SingleCellExperiment)
library(scuttle)
library(muscat)
library(edgeR)
library(ggplot2)

# Pseudo-bulk ####
sce.obj <- readRDS("snRNA_SCE_CellAnnotation.rds")
# We remove the cell-type with less than 500 cells
sce.obj <- sce.obj[,!c(sce.obj$Azimuth.labels=="Car3"|sce.obj$Azimuth.labels=="Endo"|sce.obj$Azimuth.labels=="Lamp5"|
                   sce.obj$Azimuth.labels=="Sncg"|sce.obj$Azimuth.labels=="Sst Chodl"|sce.obj$Azimuth.labels=="VLMC")]

pseudo.bulk <- aggregateAcrossCells(sce.obj, use.assay.type = "counts", id=DataFrame(label=sce.obj$Azimuth.labels, sample=sce.obj$sample_id))
colnames(pseudo.bulk) <- paste(pseudo.bulk$Azimuth.labels, pseudo.bulk$sample_id, sep="_")
pseudo.bulk <- logNormCounts(pseudo.bulk)

# We create class variable
pseudo.bulk$class <- pseudo.bulk$Azimuth.labels
pseudo.bulk$class[which(grepl("CTX", pseudo.bulk$class))] <- "Glutamatergic"

pseudo.bulk$class[which(grepl("-", pseudo.bulk$class))] <- "Other"
pseudo.bulk$class[!(pseudo.bulk$class=="Astro"|pseudo.bulk$class=="Oligo"|pseudo.bulk$class=="Glutamatergic"|pseudo.bulk$class=="Other")] <- "GABAergic"

# Save pseudo-bulk
saveRDS(pseudo.bulk, file = "snRNA_PseudoBulk.rds")

# MDS plot ####
# If this error appears: Error in element_line(linewidth = 0.2, color = "lightgrey") :unused argument (linewidth = 0.2)
source("pbMDSfunction.R")

load("AllenColorLabel.RData")
# by neuronal subclass
prep.SCE <- prepSCE(pseudo.bulk[,c(pseudo.bulk$class=="Glutamatergic"|pseudo.bulk$class=="GABAergic")],
                    kid = "Azimuth.labels", # subpopulation assignments
                    gid = "condition",  # group IDs 
                    sid = "sample_id",   # sample_id IDs 
                    drop = TRUE)

pb <- aggregatesce.obj(prep.SCE, assay="logcounts",by = c("cluster_id", "sample_id"))

neuronal.color <- subclass.color[-c(1:3,12:14,16,17,19,21)] # remove Non-Neuronal labels

mds_plot <- pbMDS(pb) + scale_color_manual(values=neuronal.color)
ggsave(mds_plot, file=paste("snRNA_MDSplot_neuronal_subclass.pdf",sep = ""), width = 20, height = 20, units = "cm")

# by GABA-ergic labels 
prep.SCE <- prepSCE(pseudo.bulk[,pseudo.bulk$class=="GABAergic"],
                    kid = "Azimuth.labels", # subpopulation assignments
                    gid = "condition",  # group IDs 
                    sid = "sample_id",   # sample_id IDs 
                    drop = TRUE)
pb <- aggregatesce.obj(prep.SCE, assay="logcounts",by = c("cluster_id", "sample_id"))

gaba.color <- neuronal.color[c(9:11)] 

mds_plot <- pbMDS(pb) + scale_color_manual(values=gaba.color)
ggsave(mds_plot, file=paste("snRNA_MDSplot_GABAergic.pdf",sep = ""), width = 20, height = 20, units = "cm")

# by Glutamatergic labels 
prep.SCE <- prepSCE(pseudo.bulk[,pseudo.bulk$class=="Glutamatergic"],
                    kid = "Azimuth.labels", # subpopulation assignments
                    gid = "condition",  # group IDs 
                    sid = "sample_id",   # sample_id IDs 
                    drop = TRUE)
pb <- aggregatesce.obj(prep.SCE, assay="logcounts",by = c("cluster_id", "sample_id"))

gluta.color <- neuronal.color[-c(9:11)] 

mds_plot <- pbMDS(pb) + scale_color_manual(values=gluta.color)
ggsave(mds_plot, file=paste("snRNA_MDSplot_Glutamatergic.pdf",sep = ""), width = 20, height = 20, units = "cm")

# by Non-Neuronal subclass
prep.SCE <- prepSCE(pseudo.bulk[,(pseudo.bulk$class=="Astro"|pseudo.bulk$class=="Oligo"|pseudo.bulk$class=="Other")],
                    kid = "Azimuth.labels", # subpopulation assignments
                    gid = "condition",  # group IDs 
                    sid = "sample_id",   # sample_id IDs 
                    drop = TRUE)
pb <- aggregatesce.obj(prep.SCE, assay="logcounts",by = c("cluster_id", "sample_id"))

non.neuronal.color <- subclass.color[-c(2:12,15,17:21)] # select Non-Neuronal labels

mds_plot <- pbMDS(pb) + scale_color_manual(values=non.neuronal.color)
ggsave(mds_plot, file=paste("snRNA_MDSplot_NonNeuronal_subclass.pdf",sep = ""), width = 20, height = 20, units = "cm")

# MDS plot on Negative Controls
NegCtrls <- read.table("SD_Negative_Controls.txt") #Negative Controls from Microarray

# We select negative controls 
pseudo.neg <- pseudo.bulk[rownames(pseudo.bulk) %in% NegCtrls$x,]

# by neuronal subclass
prep.SCE <- prepSCE(pseudo.neg[,!(pseudo.neg$class=="Astro"|pseudo.neg$class=="Oligo"|pseudo.neg$class=="Other")],
                    kid = "Azimuth.labels", # subpopulation assignments
                    gid = "condition",  # group IDs 
                    sid = "sample",   # sample IDs 
                    drop = TRUE)
pb <- aggregatesce.obj(prep.SCE, assay="logcounts",by = c("cluster_id", "sample_id"))

mds_plot <- pbMDS(pb) + scale_color_manual(values=neuronal.color)
ggsave(mds_plot, file=paste("snRNA_MDSplot_NegControls_neuronal_subclass.pdf",sep = ""), width = 20, height = 20, units = "cm")
