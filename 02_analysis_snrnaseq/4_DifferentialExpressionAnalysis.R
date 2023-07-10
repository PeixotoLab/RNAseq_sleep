library(SingleCellExperiment)
library(scuttle)
library(edgeR)
library(ggplot2)
library(RUVSeq)
library(ggrepel)
library(EDASeq)

# Differential Expression Analysis
# Load Pseudo-bulk 
pseudo.bulk <- readRDS("snRNA_PseudoBulk.rds")
# We selected only neuronal labels
pseudo.bulk <- pseudo.bulk[,pseudo.bulk$class=="Glutamatergic"|pseudo.bulk$class=="GABAergic"]

# Remove Unwanted Variation - RUVs k=1 ####
# Load negative control genes
NegCtrls <- read.table("SD_Negative_Controls.txt")

# We randomly selected 10% negative control genes, that we didn't use for fit RUVs
set.seed(23)
random.NegCtrls <- sample(NegCtrls$x,403)

NegCtrls <- NegCtrls[!(NegCtrls$x %in% random.NegCtrls),]
NegCtrls <- intersect(rownames(pseudo.bulk), NegCtrls)

# We computed RUVs for each cell-type
CellType.PseudoBulk <- RUVs1W <- ruvedSExprData <- list()
for (i in 1:length(levels(factor(pseudo.bulk$Azimuth.labels)))) {
  CellType.PseudoBulk[[i]] <- pseudo.bulk[,pseudo.bulk$Azimuth.labels==levels(factor(pseudo.bulk$Azimuth.labels))[i]]
  
  # A matrix specifying the replicates constructed
  groups <- makeGroups(CellType.PseudoBulk[[i]]$condition)
  
  ruvedSExprData[[i]] <- RUVs(as.matrix(round(counts(CellType.PseudoBulk[[i]]))), cIdx=NegCtrls, scIdx=groups, k=1)
  # We saved the factors of unwanted variation in a list object
  RUVs1W[[i]] <- ruvedSExprData[[i]]$W
} 

# Negative binomial generalized model (GLM) ####
resRUVs1 <- df.RUVs1 <- list()
for (i in 1:length(levels(factor(pseudo.bulk$Azimuth.labels)))) {
  y <- DGEList(counts(CellType.PseudoBulk[[i]]), samples=colData(CellType.PseudoBulk[[i]]))
  # The genes were filtered for each cell-type
  keep <- filterByExpr(y, group=CellType.PseudoBulk[[i]]$condition)
  y <- y[keep,]
  
  # Design Matrix
  design <- model.matrix(~0+y$samples$condition + RUVs1W[[i]], y$samples) 
  colnames(design) <- c("HC", "SD", "W_1")
  
  y <- estimateDisp(y, design)
  
  fit <- glmFit(y$counts, design, dispersion = y$trended.dispersion, robust = TRUE, offset = 0)
  
  # Contrast creation (SD vs HC)
  contrast <- makeContrasts(SD - HC,levels=design)
  
  resRUVs1[[i]] <- glmLRT(fit, contrast = contrast)
  
  df.RUVs1[[i]] <- as.data.frame(resRUVs1[[i]]$table)
  df.RUVs1[[i]] <- df.RUVs1[[i]][order(df.RUVs1[[i]]$PValue),]
  
  FDR <- as.data.frame(topTags(resRUVs1[[i]], n=length(rownames(pseudo.bulk)))[,5]) 
  
  df.RUVs1[[i]] <- cbind(df.RUVs1[[i]], FDR)
  df.RUVs1[[i]] <- df.RUVs1[[i]][order(df.RUVs1[[i]]$FDR),]
}  

# Save results
saveRDS(resRUVs1, file = "snRNA_edgeRResults_RUVs1.rds")
saveRDS(df.RUVs1, file = "snRNA_edgeR_DFResults_RUVs1.rds")

# Load positive control genes
PosCtrls <- read.table("Additional_File2_Positive_Controls.txt", header = T)
# Load DGE from Fishpond analysis
DGE_Fishpond <- read.csv("Fishpond_Annotated_DGE.csv")

# N. of genes detected
n.genes <- lapply(df.RUVs1, function(u) length(rownames(u)))

# Select DE genes
df.DEGs <- lapply(df.RUVs1, function(u) u[u$FDR<0.05,]) 
# N. of DEGs
n.DEGs <- lapply(df.DEGs, function(u) length(rownames(u))) 

df.UpDEGs <- lapply(df.DEGs, function(u) u[u$logFC>=0,]) 
# N. of DEGs Up
n.UpDEGs <-  lapply(df.UpDEGs, function(u) length(rownames(u))) 

df.DownDEGs <- lapply(df.DEGs, function(u) u[u$logFC<0,])
# N. of DEGs Down 
n.DownDEgs <-  lapply(df.DownDEGs, function(u) length(rownames(u))) 

# Intersection between single-nuclear and positive control genes
n.PosCtrls <- lapply(df.RUVs1, function(u) length(intersect(rownames(u), PosCtrls$Gene_ID))) 
# Intersection between DEGs and positive control genes
n.DEGs.PosCtrls <- lapply(df.DEGs, function(u) length(intersect(rownames(u), PosCtrls$Gene_ID)))
# Percetage of DGE positive controls/ Positive Controls detected
perc.DEGs.PosCtrls <- mapply(function(u,y) 100*round(u/y,4), n.DEGs.PosCtrls, n.PosCtrls) 

# Intersection between single-nuclear and DGE from Fishpond
n.FP <- lapply(df.RUVs1, function(u) length(intersect(rownames(u), DGE_Fishpond$Gene_Stable_ID))) 
# Intersection between dataset DGE and DGE from Fishpond
n.DGE.FP <- lapply(df.DEGs, function(u) length(intersect(rownames(u), DGE_Fishpond$Gene_Stable_ID)))
# Percetage of DGE Fishpond/ DGE Fishpond detected
perc.DGE.FP <- mapply(function(u,y) 100*round(u/y,4), n.DGE.FP, n.FP)

n.nuclei <- matrix(pseudo.bulk$ncells,11,6,byrow=TRUE) 
n.nuclei <- rowSums(n.nuclei) # N. of nuclei

tab <- cbind(n.nuclei, as.matrix(n.genes), as.matrix(n.UpDEGs),as.matrix(n.DownDEgs), as.matrix(perc.DEGs.PosCtrls), 
             as.matrix(perc.DGE.FP))
colnames(tab) <- c("N. of nuclei","N. of genes detected","N. of DEGs Up","N. of DEGs Down","% DE Pos Ctrls","% DEGs from bulk")
rownames(tab) <- levels(factor(pseudo.bulk$Azimuth.labels))

# Scatter plot ####
load("AllenColorLabel.RData")
neuronal.color <- subclass.color[-c(1:3,12:14,16,17,19,21)]

tab <- as.data.frame(cbind(n.nuclei, as.matrix(n.DEGs)))
colnames(tab) <- c("nuclei","DEG")
rownames(tab) <- levels(factor(pseudo.bulk$Azimuth.labels))
tab$Label <- factor(rownames(tab))

tab$DEG <- as.numeric(tab$DEG)
tab$nuclei <- as.numeric(tab$nuclei)

p1 <- ggplot(tab, aes(x=log(nuclei),y=log(DEG)))+ 
  geom_point(aes(color = Label), size = 3) + scale_color_manual(values=neuronal.color) + 
  theme_classic()+theme(legend.position="none", axis.text=element_text(size=7), axis.title=element_text(size=7))+
  stat_smooth(method = "lm", formula = y ~ poly(x, 2),se=FALSE,col="grey") +
  geom_label_repel(aes(label = Label), size=2,
                   box.padding   = 0.1, 
                   point.padding = 0,
                   segment.color = 'grey50', colour =neuronal.color)
ggsave("Peixoto_Figure5_part1.jpg",  width = 9, height = 9, units = "cm")
ggsave("Peixoto_Figure5_part1.pdf",  width = 9, height = 9, units = "cm")

# PCA plot ####
label <- levels(factor(pseudo.bulk$Azimuth.labels))

# Calculate the logarithm of normalized counts matrix
normalizedCounts <- lapply(ruvedSExprData, function(u) log1p(u$normalizedCounts))

# Create a SeqExpressionSet object 
seq.RUVs1 <- mapply(function(u,y) 
  newSeqExpressionSet(as.matrix(round(log1p(u))), 
                      phenoData = data.frame(colData(y),row.names=colnames(y))),
  normalizedCounts, CellType.PseudoBulk)

file.name <- paste("Peixoto_Figure5_Supplement1_part1.pdf",sep = "")
#jpeg(file = file.name, width = 794, height = 1134, units = "px")
pdf(file = file.name, width = 8.27, height = 11.81)
par(mfrow = c(4, 3)) 
for (i in seq_along(label)) {
  print(i)
  plotPCA(seq.RUVs1[[i]], xlim = c(-1, 1),label=FALSE,  pch=20, theme_size=7,
          col=as.numeric(as.factor(seq.RUVs1[[i]]$condition)), cex=2, 
          main=paste(label[i],sep=" ") )
  #legend("right", inset = c(-0.17,0), legend=levels(factor(seq.RUVs1[[i]]$condition)),  pch=20, cex=1,
  #       col= levels(as.factor(as.numeric(as.factor(seq.RUVs1[[i]]$condition)))), xpd = TRUE, 
  #       horiz = FALSE,  bty = "n",title="Condition")
}
dev.off()

## Volcano plot ####
# Read smaller list of positive control genes 
Highlight_Genes <- readxl::read_excel("Highlight_Genes_Subset.xlsx")
Highlight_Genes <- Highlight_Genes[!is.na(Highlight_Genes$Ensembl_ID),]
Highlight_Genes <- Highlight_Genes[!duplicated(Highlight_Genes$Gene_Name),]
Highlight_Genes <- Highlight_Genes[order(Highlight_Genes$Ensembl_ID),]

# Select L4/5 IT CTX and Sst DEA results
ll <- list(df.RUVs1[[2]], df.RUVs1[[10]])
ll <- lapply(ll, function(x) x[order(rownames(x)),])

Highlight_GenesLabel <- list(L45 = as.data.frame(Highlight_Genes[Highlight_Genes$Ensembl_ID %in% intersect(rownames(ll[[1]]), Highlight_Genes$Ensembl_ID),]),
                             Sst = as.data.frame(Highlight_Genes[Highlight_Genes$Ensembl_ID %in% intersect(rownames(ll[[2]]), Highlight_Genes$Ensembl_ID),]))

for (i in seq_along(ll)) {
  ll[[i]]$Symbol_HighlightGenes <- ""
  ll[[i]]$Symbol_HighlightGenes[rownames(ll[[i]]) %in% Highlight_GenesLabel[[i]]$Ensembl_ID] <- Highlight_GenesLabel[[i]]$Gene_Name
  ll[[i]]$Significance <- "No Significant"
  ll[[i]]$Significance[ll[[i]]$FDR<0.05] <- "Significant"
  inter <- intersect(rownames(ll[[i]]), PosCtrls$Gene_ID)
  ll[[i]]$Significance[rownames(ll[[i]]) %in% inter] <- "SignificantPos"
  ll[[i]]$HighlightGenes <- ll[[i]]$Symbol_HighlightGenes
  ll[[i]]$HighlightGenes[ll[[i]]$HighlightGenes==""] <- FALSE
  ll[[i]]$HighlightGenes[ll[[i]]$HighlightGenes!=FALSE] <- TRUE
}


df.DEG <- lapply(ll, function(x) x[x$FDR<0.05,])
n.DEG <- lapply(df.DEG, function(x) length(rownames(x)))
n.DEG.PosCtrls <- lapply(df.DEG, function(x) length(intersect(rownames(x), PosCtrls$Gene_ID)))

p1 <- ggplot(data=ll[[1]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  geom_point(size=1)+xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")+
  scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"))+
  ggtitle(paste("L4/5 IT CTX",sep = "")) +
  geom_label_repel(data=ll[[1]][ll[[1]]$Significance=="No Significant",],aes(label = ifelse(HighlightGenes, Symbol_HighlightGenes, "")),
                   fill = "white",max.overlaps = Inf, fontface=2 , direction="x", xlim=c(-3,3),show.legend = FALSE, size=2)+
  
  geom_label_repel(data=ll[[1]][(ll[[1]]$Significance=="Significant"&ll[[1]]$logFC<0),],aes(label = ifelse(HighlightGenes, Symbol_HighlightGenes, "")),
                   fill = "white",max.overlaps = Inf, fontface=2 , direction="both", hjust="left", ylim=c(3,15), show.legend = FALSE, size=2)+
  geom_label_repel(data=ll[[1]][(ll[[1]]$Significance=="Significant"&ll[[1]]$logFC>=0),],aes(label = ifelse(HighlightGenes, Symbol_HighlightGenes, "")),
                   fill = "white",max.overlaps = Inf, fontface=2 , direction="both", hjust="right", ylim=c(3,15), show.legend = FALSE, size=2)+
  geom_label_repel(data=ll[[1]][(ll[[1]]$Significance=="SignificantPos"),],aes(label = ifelse(HighlightGenes, Symbol_HighlightGenes, "")),
                   fill = "white",max.overlaps = Inf, fontface=2 , direction="both", ylim=c(10,30), show.legend = FALSE, size=2)

p2 <- ggplot(data=ll[[2]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  geom_point(size=1)+xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")+
  scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"))+
  ggtitle(paste("Sst",sep = "")) +
  geom_label_repel(data=ll[[2]][ll[[2]]$Significance=="No Significant",],aes(label = ifelse(HighlightGenes, Symbol_HighlightGenes, "")),
                   fill = "white",max.overlaps = Inf, fontface=2 , direction="x",show.legend = FALSE, size=2, xlim=c(-3,3))+
  geom_label_repel(data=ll[[2]][(ll[[2]]$Significance=="Significant"&ll[[2]]$logFC>0),],aes(label = ifelse(HighlightGenes, Symbol_HighlightGenes, "")),
                   fill = "white",max.overlaps = Inf, fontface=2 , direction="both", show.legend = FALSE, size=2)+
  geom_label_repel(data=ll[[2]][(ll[[2]]$Significance=="Significant"&ll[[2]]$logFC<0),],aes(label = ifelse(HighlightGenes, Symbol_HighlightGenes, "")),
                   fill = "white",max.overlaps = Inf, fontface=2 , direction="x", hjust="left", ylim=c(3,15), show.legend = FALSE, size=2)+
  
  geom_label_repel(data=ll[[2]][(ll[[2]]$Significance=="SignificantPos" &ll[[2]]$logFC<0),],aes(label = ifelse(HighlightGenes, Symbol_HighlightGenes, "")),
                   fill = "white",max.overlaps = Inf, fontface=2 , direction="both", ylim=c(2,10), show.legend = FALSE, size=2)+
  geom_label_repel(data=ll[[2]][(ll[[2]]$Significance=="SignificantPos" &ll[[2]]$logFC>0),],aes(label = ifelse(HighlightGenes, Symbol_HighlightGenes, "")),
                   fill = "white",max.overlaps = Inf, fontface=2 , direction="both", ylim=c(5,10), show.legend = FALSE, size=2)

patch <- p1|p2
ggsave("Peixoto_Figure5_part2.jpg",  width = 18, height = 9, units = "cm")
ggsave("Peixoto_Figure5_part2.pdf",  width = 18, height = 9, units = "cm")

### Volcano plot for each cell-type ####
for (i in seq_along(df.RUVs1)) {
  df.RUVs1[[i]]$Significance <- "No Significant"
  df.RUVs1[[i]]$Significance[df.RUVs1[[i]]$FDR<0.05] <- "Significant"
  inter <- intersect(rownames(df.RUVs1[[i]]), PosCtrls$Gene_ID)
  df.RUVs1[[i]]$Significance[rownames(df.RUVs1[[i]]) %in% inter] <- "SignificantPos"
}
label <- levels(factor(pseudo.bulk$Azimuth.labels))

df.DEG <- lapply(df.RUVs1, function(x) x[x$FDR<0.05,])
n.DEG <- lapply(df.DEG, function(x) length(rownames(x)))
n.DEG.PosCtrls <- lapply(df.DEG, function(x) length(intersect(rownames(x), PosCtrls$Gene_ID)))

p1 <- ggplot(data=df.RUVs1[[1]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  geom_point(size=1)+xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")+
  scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"))+
  ggtitle(paste(label[1],sep = ""))

p2 <- ggplot(data=df.RUVs1[[2]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  geom_point(size=1)+xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")+
  scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"))+
  ggtitle(paste(label[2],sep = ""))

p3 <- ggplot(data=df.RUVs1[[3]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  geom_point(size=1)+xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")+
  scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"))+
  ggtitle(paste(label[3],sep = ""))

p4 <- ggplot(data=df.RUVs1[[4]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  geom_point(size=1)+xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")+
  scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"))+
  ggtitle(paste(label[4],sep = ""))

p5 <- ggplot(data=df.RUVs1[[5]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  geom_point(size=1)+xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")+
  scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"))+
  ggtitle(paste(label[5],sep = ""))

p6 <- ggplot(data=df.RUVs1[[6]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  geom_point(size=1)+xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")+
  scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"))+
  ggtitle(paste(label[6],sep = ""))

p7 <- ggplot(data=df.RUVs1[[7]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  geom_point(size=1)+xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")+
  scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"))+
  ggtitle(paste(label[7],sep = ""))

p8 <- ggplot(data=df.RUVs1[[8]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  geom_point(size=1)+xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")+
  scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"))+
  ggtitle(paste(label[8],sep = ""))

p9 <- ggplot(data=df.RUVs1[[9]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  geom_point(size=1)+xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")+
  scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"))+
  ggtitle(paste(label[9],sep = ""))

p10 <- ggplot(data=df.RUVs1[[10]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  geom_point(size=1)+xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")+
  scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"))+
  ggtitle(paste(label[10],sep = ""))

p11 <- ggplot(data=df.RUVs1[[11]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  geom_point(size=1)+xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")+
  scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"))+
  ggtitle(paste(label[11],sep = ""))

p12 <- ggplot(data=df.RUVs1[[11]], aes(x=logFC, y=-log10(PValue),col=Significance))+
  xlim(-5,5) + theme_classic(base_size = 7)+theme(legend.position="none")

patch <- (p1|p2|p3)/(p4|p5|p6)/(p7|p8|p9)/(p10|p11|plot_spacer())
ggsave("Peixoto_Figure5_Supplement1_part2.jpg",  width = 21, height = 30, units = "cm")
ggsave("Peixoto_Figure5_Supplement1_part2.pdf",  width = 21, height = 30, units = "cm")

# Histogram of p-value distribution #####
file.name <- paste("Peixoto_Figure5_Supplement1_part3.pdf",sep = "")
#jpeg(file = file.name, width = 794, height = 1134, units = "px")
pdf(file = file.name, width = 8.27, height = 11.81)
par(mfrow = c(4, 3)) 
for (i in seq_along(label)) {
  print(i)
  hist(df.RUVs1[[i]]$PValue, xlab = "p-value", main=paste(label[i],sep=" "))
}
dev.off()
                         
# Negative binomial generalized model (GLM) on Pseudo-bulk sum ####
# Aggregate across sample
PseudoBulkSum <- aggregateAcrossCells(pseudo.bulk, use.assay.type = "counts", id=DataFrame(sample=pseudo.bulk$sample_id))
colnames(PseudoBulkSum) <- paste(substr(PseudoBulkSum$sample_id,1,1),PseudoBulkSum$condition, sep="")

# A matrix specifying the replicates constructed
groups <- makeGroups(PseudoBulkSum$condition)

ruvedSExprData <- RUVs(as.matrix(round(counts(PseudoBulkSum))), cIdx=NegCtrls, scIdx=groups, k=1)
# We saved the factor of unwanted variation
RUVs1W.Sum <- ruvedSExprData$W

y <- DGEList(counts(PseudoBulkSum), samples=colData(PseudoBulkSum))
keep <- filterByExpr(y, group=PseudoBulkSum$condition)
y <- y[keep,]

# Design Matrix
design <- model.matrix(~0+y$samples$condition + RUVs1W.Sum, y$samples) 
colnames(design) <- c("HC", "SD", "W_1")

y <- estimateDisp(y, design)

fit <- glmFit(y$counts, design, dispersion = y$trended.dispersion, robust = TRUE, offset = 0)

# Contrast creation (SD vs HC)
contrast <- makeContrasts(SD - HC,levels=design)

resRUVs1.Sum <- glmLRT(fit, contrast = contrast)

df.RUVs1.Sum <- as.data.frame(resRUVs1.Sum$table)
df.RUVs1.Sum <- df.RUVs1.Sum[order(df.RUVs1.Sum$PValue),]

FDR <- as.data.frame(topTags(resRUVs1.Sum, n=length(rownames(PseudoBulkSum)))[,5]) 

df.RUVs1.Sum <- cbind(df.RUVs1.Sum, FDR)
df.RUVs1.Sum <- df.RUVs1.Sum[order(df.RUVs1.Sum$FDR),]

n.genes <- length(rownames(df.RUVs1.Sum)) # N. of genes

df.DEGs <- df.RUVs1.Sum[df.RUVs1.Sum$FDR<0.05,] # data.frame of DEGs results 
n.DEGs <- length(rownames(df.DEGs)) # N. of DEGs

df.UpDEGs <- df.DEGs[df.DEGs$logFC>=0,]
n.UpDEGs <-  length(rownames(df.UpDEGs)) # N. of DEGs 

df.DownDEGs <- df.DEGs[df.DEGs$logFC<0,]
n.DownDEGs <-  length(rownames(df.DownDEGs)) # N. of DEGs 

n.PosCtrls <- length(intersect(rownames(df.RUVs1.Sum), PosCtrls$Gene_ID)) # Intersection between dataset genes and positive controls
n.DEGs.PosCtrls <- length(intersect(rownames(df.DEGs), PosCtrls$Gene_ID)) # Intersection between dataset DEGs and positive controls
perc.DEGs.PosCtrls <- 100*round(n.DEGs.PosCtrls/n.PosCtrls,4) # Percetage of DEGs positive controls/ Positive Controls detected

n.FP <- length(intersect(rownames(df.RUVs1.Sum), DGE_Fishpond$Gene_Stable_ID)) # Intersection between dataset genes and DEGs from Fishpond
n.DEGs.FP <- length(intersect(rownames(df.DEGs), DGE_Fishpond$Gene_Stable_ID)) # Intersection between dataset DEGs and DEGs from Fishpond
perc.DEGs.FP <- 100*round(n.DEGs.FP/n.FP,4) # Percetage of DEGs Fishpond/ DEGs FP detected
