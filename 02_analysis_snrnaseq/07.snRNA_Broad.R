library(SingleCellExperiment)
library(SingleR)
library(scran)
library(scater)
library(batchelor)
library(BiocNeighbors)

## snRNA-seq at Broad dataset ####
sce <- readRDS("se.rds")

sce$sample_id <- sce$sample_name
sce$sample_id <- substring(sce$sample_id, 1,27)

sce <- logNormCounts(sce)

sample <- levels(factor(sce$sample_id))

mouse <- list()
for (n in seq_along(sample)) {
  mouse[[n]] <- sce[,sce$sample_id==sample[n]]
}
names(mouse) <- sample

for (n in seq_along(mouse)) {
  mouse[[n]] <- logNormCounts(mouse[[n]])
}

same.genes <- Reduce(intersect, list(rownames(mouse[[1]]),rownames(mouse[[2]]),
                                     rownames(mouse[[3]]),rownames(mouse[[4]]),
                                     rownames(mouse[[5]]),rownames(mouse[[6]]),
                                     rownames(mouse[[7]]),rownames(mouse[[8]]),
                                     rownames(mouse[[9]]),rownames(mouse[[10]]),
                                     rownames(mouse[[11]]),rownames(mouse[[12]]),
                                     rownames(mouse[[13]]),rownames(mouse[[14]]),
                                     rownames(mouse[[15]])))

for (n in seq_along(mouse)) {
  set.seed(422)
  mouse[[n]] <- runPCA(mouse[[n]], ntop=length(same.genes))
  
  set.seed(422)
  mouse[[n]] <- runTSNE(mouse[[n]], dimred="PCA")
}

## Dataset labels ####
subclass.label.color <- c("#957b46","#00FF66","#c95f3f","#0BE652",
                          "#0D5B78","#50B2AD","#3E9E64","#2D8CB8",
                          "#A19922","#7044AA","#DA808C","#94AF97",
                          "#744700","#000000","#FF0000","#5100FF","#D93137",
                          "#4d7647","#ffff00","#FF9900","#B864CC","#54a04d")

names(subclass.label.color) <- levels(factor(sce$label.subclass_label))

for (n in seq_along(mouse)) {
  tsne <- plotTSNE(mouse[[n]],colour_by="label.subclass_label")+
    scale_color_manual(values=subclass.label.color)
  ggsave(tsne, file=paste("tSNEplot_broad_mouse",n,"_label.png",sep = ""), width = 14, height = 10, units = "cm")
}

## SingleR for each mouse ####
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
  ggsave(tsne, file=paste("tSNEplot_broad_mouse",n,"_10xgenomics.png",sep = ""), width = 14, height = 10, units = "cm")
}

## MNN correction####
mnn <- fastMNN(sce, batch=sce$sample_id, d=50, k=20, subset.row=rownames(sce),
               BSPARAM=BiocSingular::RandomParam(deferred=TRUE))

## Note: R Failed!
