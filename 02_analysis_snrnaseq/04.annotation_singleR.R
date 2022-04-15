library(SingleCellExperiment)
library(scRNAseq)
library(scran)
library(scater)
library(flexmix)
library(BiocParallel)
library(miQC)
library(batchelor)
library(SingleR)
library(BiocNeighbors)
library(pheatmap)
library(EnsDb.Mmusculus.v79)

load("pre_processing.RData")

mouse <- list()
for (n in seq_along(filter.polyn)) {
  cells.5perc <- length(colnames(filter.polyn[[n]]))*0.05
  mouse[[n]] <- filter.polyn[[n]][rowSums(counts(filter.polyn[[n]])>0)>=cells.5perc,]
}

same.genes <- Reduce(intersect, list(rownames(mouse[[1]]),rownames(mouse[[2]]),
                                     rownames(mouse[[3]]),rownames(mouse[[4]]),
                                     rownames(mouse[[5]]),rownames(mouse[[6]])))

for (n in seq_along(filter.polyn)) {
  mouse[[n]] <- mouse[[n]][same.genes,]
}

for (n in seq_along(mouse)) {
  set.seed(422)
  mouse[[n]] <- runPCA(mouse[[n]], ntop=length(same.genes))
  
  set.seed(422)
  mouse[[n]] <- runTSNE(mouse[[n]], dimred="PCA")
}

## SingleR ####
## reference 10x genomics (2020)####
pred19.mouse <- list()
for (n in seq_along(mouse)) {
  pred19.mouse[[n]] <- SingleR(mouse[[n]], ref=allen.10xgenomics, 
                               labels=allen.10xgenomics$subclass_label)
  mouse[[n]]$subclass.10xgenomics <- pred19.mouse[[n]]$labels
}

subclass.10xgenomics <- data.frame(allen.10xgenomics$subclass_label, allen.10xgenomics$subclass_color)
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
  ggsave(tsne, file=paste("tSNEplot_snRNA_mouse",n,"_10xgenomics.png",sep = ""), width = 14, height = 10, units = "cm")
}

## MNN correction####
mnn <- fastMNN(mouse[[1]], mouse[[2]], mouse[[3]], mouse[[4]],
               mouse[[5]], mouse[[6]], d=50, k=20, subset.row=same.genes,
               BSPARAM=BiocSingular::RandomParam(deferred=TRUE))

uncorrected <- cbind(mouse[[1]],mouse[[2]],mouse[[3]],mouse[[4]],
                     mouse[[5]],mouse[[6]])
assays(mnn)$logcounts <- logcounts(uncorrected)
assays(mnn)$counts <- counts(uncorrected)

set.seed(422)
mnn <- runTSNE(mnn, dimred="corrected", external_neighbors = TRUE, BNPARAM=AnnoyParam())

mnn$subclass.10xgenomics <- uncorrected$subclass.10xgenomics

tsne_mnn <- plotTSNE(mnn,colour_by="subclass.10xgenomics")+
  scale_color_manual(values=subclass.10xgenomics.color)
ggsave(tsne_mnn, file=paste("tSNEplot_snRNA_mnn_10xgenomics.png",sep = ""), width = 14, height = 10, units = "cm")

save(mouse, mnn, file = "mouse_annotation.RData")

## heatmap ####
mnn <- mnn[,!(mnn$subclass.10xgenomics=="Endo"|mnn$subclass.10xgenomics=="SMC-Peri"|
                mnn$subclass.10xgenomics=="Micro-PVM"|mnn$subclass.10xgenomics=="VLMC"|
                mnn$subclass.10xgenomics=="Sst Chodl")]

mnn$subclass.10xgenomics[mnn$subclass.10xgenomics=="L2/3 IT CTX-1"|mnn$subclass.10xgenomics=="L2/3 IT CTX-2"] <- "L2/3 IT CTX"
#mnn <- mnn[,!(mnn$subclass.10xgenomics=="Astro"|mnn$subclass.10xgenomics=="Oligo")]

mean<- matrix(nrow = length(rownames(mnn)), ncol =length(levels(factor(mnn$subclass.10xgenomics))))
for (i in 1:length(levels(factor(mnn$subclass.10xgenomics)))) {
  message("Mean: ", i)
  mean[,i] <- rowMeans(logcounts(mnn[,mnn$subclass.10xgenomics==levels(factor(mnn$subclass.10xgenomics))[i]])) 
}

rownames(mean) <- rownames(mnn)
colnames(mean) <- levels(factor(mnn$subclass.10xgenomics))

markers.ensembl <- which(rownames(mean)=="ENSMUSG00000020932"|
                           rownames(mean)=="ENSMUSG00000036887"|
                           rownames(mean)=="ENSMUSG00000041607"|
                           rownames(mean)=="ENSMUSG00000000567"|
                           rownames(mean)=="ENSMUSG00000033006"|
                           rownames(mean)=="ENSMUSG00000026473"|
                           rownames(mean)=="ENSMUSG00000076439"|
                           rownames(mean)=="ENSMUSG00000024617"|
                           rownames(mean)=="ENSMUSG00000070880"|
                           rownames(mean)=="ENSMUSG00000026787"|
                           rownames(mean)=="ENSMUSG00000037771"|
                           rownames(mean)=="ENSMUSG00000004366"|
                           rownames(mean)=="ENSMUSG00000029361"|
                           rownames(mean)=="ENSMUSG00000005716"|
                           rownames(mean)=="ENSMUSG00000032269"|
                           rownames(mean)=="ENSMUSG00000019772"|
                           rownames(mean)=="ENSMUSG00000032532"|
                           rownames(mean)=="ENSMUSG00000042453"|
                           rownames(mean)=="ENSMUSG00000029819"|
                           rownames(mean)=="ENSMUSG00000003657"|
                           rownames(mean)=="ENSMUSG00000070570"|
                           rownames(mean)=="ENSMUSG00000053310"|
                           rownames(mean)=="ENSMUSG00000029563"|
                           rownames(mean)=="ENSMUSG00000005917"|
                           rownames(mean)=="ENSMUSG00000036192"|
                           rownames(mean)=="ENSMUSG00000004151"|
                           rownames(mean)=="ENSMUSG00000029705"|
                           rownames(mean)=="ENSMUSG00000024642")

mean_mark <- mean[markers.ensembl,]

symbol_allen <- rownames(mean_mark)
map <- mapIds(EnsDb.Mmusculus.v79, keys= symbol_allen, keytype = "GENEID", column = "SYMBOL")
stopifnot(length(map) == nrow(mean_mark))
rownames(mean_mark) <- map

pheatmap_mark <- pheatmap(log2(mean_mark+1), color=viridis::viridis(100))
ggsave(pheatmap_mark, file=paste("pheatmap_mark.png",sep = ""), width = 14, height = 14, units = "cm")

pheatmap_mark_scale <- pheatmap(log2(mean_mark+1), scale = "row", color=viridis::viridis(100))
ggsave(pheatmap_mark_scale, file=paste("pheatmap_mark_scale.png",sep = ""), width = 14, height = 14, units = "cm")
