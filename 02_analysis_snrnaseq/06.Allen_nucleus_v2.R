library(SingleCellExperiment)
library(SingleR)
library(scran)
library(scater)
library(batchelor)
library(BiocNeighbors)


## Nucleus 10x v2 dataset#####
sce <- readRDS("nucleusV2_SE.RDS")

mouse <- list(
  sce[, sce$Donor=="375738"],
  sce[, sce$Donor=="376425"],
  sce[, sce$Donor=="383185"]
)

for (n in seq_along(mouse)) {
  mouse[[n]] <- logNormCounts(mouse[[n]])
}

same.genes <- Reduce(intersect, list(rownames(mouse[[1]]),rownames(mouse[[2]]),
                                     rownames(mouse[[3]])))

for (n in seq_along(mouse)) {
  set.seed(422)
  mouse[[n]] <- runPCA(mouse[[n]], ntop=length(same.genes))
  
  set.seed(422)
  mouse[[n]] <- runTSNE(mouse[[n]], dimred="PCA")
}

## Dataset labels ####
subclass.label.color <- c("#957b46","#c95f3f","#0BE652",
                          "#0D5B78","#50B2AD","#3E9E64",
                          "#2D8CB8","#A19922","#5100FF",
                          "#7044AA","#DA808C","#94AF97",
                          "#744700","#000000","#D93137",
                          "#ffff00","#FF9900","#B864CC")

names(subclass.label.color) <- levels(factor(mouse[[1]]$subclass_label))

for (n in seq_along(mouse)) {
  tsne <- plotTSNE(mouse[[n]],colour_by="subclass_label")+
    scale_color_manual(values=subclass.label.color)
  ggsave(tsne, file=paste("tSNEplot_nucleus_mouse",n,"_label.png",sep = ""), width = 14, height = 10, units = "cm")
}

## SingleR for each mouse####
## 10x genomics (2020) ####
load("allen_pseudo_20.RData")
rownames(allen.pseudo.20) <- rowData(allen.pseudo.20)$symbol

xgenomics.mouse <- list()
for (n in seq_along(mouse)) {
  xgenomics.mouse[[n]] <- SingleR(mouse[[n]], ref=allen.pseudo.20, 
                                  labels=allen.pseudo.20$subclass_label)
  mouse[[n]]$subclass.10xgenomics <- xgenomics.mouse[[n]]$labels
}

subclass.10xgenomics <- data.frame(allen.pseudo.20$subclass_label, allen.pseudo.20$subclass_color)
colnames(subclass.10xgenomics) <- c("Label", "color")

# change Astro color
subclass.10xgenomics$color[subclass.10xgenomics$color=="#665C47"] <- "#957b46"
# change Oligo color
subclass.10xgenomics$color[subclass.10xgenomics$color=="#53776C"] <- "#744700"
# change SMC-Peri color
subclass.10xgenomics$color[subclass.10xgenomics$color=="#807059"] <- "#4d7647"
# change VLMC color
subclass.10xgenomics$color[subclass.10xgenomics$color=="#697255"] <- "#54a04d"
# change Endo color
subclass.10xgenomics$color[subclass.10xgenomics$color=="#8D6C62"] <- "#c95f3f"
# change Sncg color 
subclass.10xgenomics$color[subclass.10xgenomics$color=="#D3408D"] <- "#ffff00"

subclass.10xgenomics.color <- subclass.10xgenomics$color
names(subclass.10xgenomics.color) <- subclass.10xgenomics$Label

for (n in seq_along(mouse)) {
  tsne <- plotTSNE(mouse[[n]],colour_by="subclass.10xgenomics")+
    scale_color_manual(values=subclass.10xgenomics.color)
  ggsave(tsne, file=paste("tSNEplot_nucleus_mouse",n,"_10xgenomics.png",sep = ""), width = 14, height = 10, units = "cm")
}

## MNN correction####
same.genes <- Reduce(intersect, list(rownames(mouse[[1]]),rownames(mouse[[2]]),
                                     rownames(mouse[[3]])))

mnn <- fastMNN(mouse[[1]], mouse[[2]], mouse[[3]], d=50, k=20, subset.row=same.genes,
               BSPARAM=BiocSingular::RandomParam(deferred=TRUE))

uncorrected <- cbind(mouse[[1]],mouse[[2]],mouse[[3]])
uncorrected <- uncorrected[same.genes,]

assays(mnn)$logcounts <- logcounts(uncorrected)
assays(mnn)$counts <- counts(uncorrected)
mnn$subclass_label <- uncorrected$subclass_label

set.seed(422)
mnn <- runTSNE(mnn, dimred="corrected", external_neighbors = TRUE, BNPARAM=AnnoyParam())

subclass.label.color <- c("#957b46","#c95f3f","#0BE652",
                          "#0D5B78","#50B2AD","#3E9E64",
                          "#2D8CB8","#A19922","#5100FF",
                          "#7044AA","#DA808C","#94AF97",
                          "#744700","#000000","#D93137",
                          "#ffff00","#FF9900","#B864CC")

names(subclass.label.color) <- levels(factor(mnn$subclass_label))

tsne <- plotTSNE(mnn,colour_by="subclass_label")+
  scale_color_manual(values=subclass.label.color)
ggsave(tsne, file=paste("tSNEplot_nucleus_mnn_label.png",sep = ""), width = 14, height = 10, units = "cm")

## SingleR on MNN #####
## 10x genomics (2020) ####
load("allen_pseudo_20.RData")
rownames(allen.pseudo.20) <- rowData(allen.pseudo.20)$symbol

xgenomics <- SingleR(mnn, ref=allen.pseudo.20, 
                     labels=allen.pseudo.20$subclass_label, assay.type.test = "reconstructed")

mnn$subclass.10xgenomics <- xgenomics$labels

subclass.10xgenomics <- data.frame(allen.pseudo.20$subclass_label, allen.pseudo.20$subclass_color)
colnames(subclass.10xgenomics) <- c("Label", "color")

# change Astro color
subclass.10xgenomics$color[subclass.10xgenomics$color=="#665C47"] <- "#957b46"
# change Oligo color
subclass.10xgenomics$color[subclass.10xgenomics$color=="#53776C"] <- "#744700"
# change SMC-Peri color
subclass.10xgenomics$color[subclass.10xgenomics$color=="#807059"] <- "#4d7647"
# change VLMC color
subclass.10xgenomics$color[subclass.10xgenomics$color=="#697255"] <- "#54a04d"
# change Endo color
subclass.10xgenomics$color[subclass.10xgenomics$color=="#8D6C62"] <- "#c95f3f"
# change Sncg color 
subclass.10xgenomics$color[subclass.10xgenomics$color=="#D3408D"] <- "#ffff00"

subclass.10xgenomics.color <- subclass.10xgenomics$color
names(subclass.10xgenomics.color) <- subclass.10xgenomics$Label

tsne <- plotTSNE(mnn,colour_by="subclass.10xgenomics")+
  scale_color_manual(values=subclass.10xgenomics.color)
ggsave(tsne, file=paste("tSNEplot_nucleus_mnn_10xgenomics.png",sep = ""), width = 14, height = 10, units = "cm")


