library(SingleCellExperiment)
library(scRNAseq)
library(scran)
library(scater)
library(flexmix)
library(BiocParallel)
library(miQC)
library(batchelor)

load("sce_sum.RData")

## miQC for each mouse####
qc.sce <- model.linear <- model.polyn <- filter.linear <- filter.polyn <- list()
for (n in seq_along(sce.sum)) {
  is.mito <- which(rowData(sce.sum[[n]])$CHR=="MT")
  feature_ctrls <- list(mito = rownames(sce.sum[[n]])[is.mito])
  qc.sce[[n]] <- addPerCellQC(sce.sum[[n]], subsets = feature_ctrls, BPPARAM = MulticoreParam())
  
  # linear model
  model.linear[[n]] <- mixtureModel(qc.sce[[n]])
  filter.linear[[n]] <- filterCells(qc.sce[[n]], model.linear[[n]])
  
  # polynomial model
  model.polyn[[n]] <- mixtureModel(qc.sce[[n]], model_type = "polynomial")
  filter.polyn[[n]] <- filterCells(qc.sce[[n]], model.polyn[[n]])
}

for (n in seq_along(sce.sum)) {
  # linear
  qc_model <- plotModel(qc.sce[[n]], model.linear[[n]])+
    xlab("Uniquely Found Genes") + ylab("Mitochondrial reads percentage")
  ggsave(qc_model, file=paste("qc_mitureModel_mouse",n,"_linear.png",sep = ""), width = 14, height = 14, units = "cm")
  
  # polynomial
  qc_model <- plotModel(qc.sce[[n]], model.polyn[[n]])+
    xlab("Uniquely Found Genes") + ylab("Mitochondrial reads percentage")
  ggsave(qc_model, file=paste("qc_mitureModel_mouse",n,"_polynomial.png",sep = ""), width = 14, height = 14, units = "cm")
}

## normalization by deconvolution #####
# I used filter.polyn data
lib.size <- f <- list()
for (n in seq_along(sce.sum)) {
  #set.seed(100)
  #lib.size[[n]] <- librarySizeFactors(filter.polyn[[n]])
  #clust <- quickCluster(filter.polyn[[n]])
  #f[[n]] <- computeSumFactors(filter.polyn[[n]], min.mean=0.1, cluster=clust)
  
  # normalized dataset 
  filter.polyn[[n]] <- logNormCounts(filter.polyn[[n]])
}

# for (n in seq_along(sce.sum)) {
#  size.factors <- as.data.frame(cbind(lib.size[[n]], sizeFactors(f[[n]])))
#  colnames(size.factors) <- c("lib", "decov")
  
#  decov_plot <- ggplot(size.factors, aes(lib, decov)) + geom_point()+  geom_abline(col="red")+
#    labs(x = "Library size factor", y="Deconvolution size factor")
  
#  ggsave(decov_plot, file=paste("norm_deconv_mouse",n,".png",sep = ""), width = 14, height = 14, units = "cm")
# }

save(filter.polyn, file = "pre_processing.RData")