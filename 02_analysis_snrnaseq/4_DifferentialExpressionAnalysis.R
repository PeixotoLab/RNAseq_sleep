library(SingleCellExperiment)
library(scuttle)
library(edgeR)
library(ggplot2)
library(RUVSeq)
library(ggrepel)

# Differential Expression Analysis with RUVs k=1 ####
# Load Pseudo-bulk 
pseudo.bulk <- readRDS("snRNA_PseudoBulk.rds")
# We selected only neuronal labels
pseudo.bulk <- pseudo.bulk[,pseudo.bulk$class=="Glutamatergic"|pseudo.bulk$class=="GABAergic"]

# Load negative control genes
NegCtrls <- read.table("SD_Negative_Controls.txt")

# We randomly selected 10% negative control genes, that we didn't use for fit RUVs
set.seed(23)
random.NegCtrls <- sample(NegCtrls$x,403)

NegCtrls <- NegCtrls[!(NegCtrls$x %in% random.NegCtrls),]
NegCtrls <- intersect(rownames(pseudo.bulk), NegCtrls)

# We computed RUVs for each cell-type
CellType.PseudoBulk <- RUVs1W <- list()
for (i in 1:length(levels(factor(pseudo.bulk$Azimuth.labels)))) {
  CellType.PseudoBulk[[i]] <- pseudo.bulk[,pseudo.bulk$Azimuth.labels==levels(factor(pseudo.bulk$Azimuth.labels))[i]]
  
  # A matrix specifying the replicates constructed
  groups <- makeGroups(CellType.PseudoBulk[[i]]$condition)
  
  ruvedSExprData <- RUVs(as.matrix(round(counts(CellType.PseudoBulk[[i]]))), cIdx=SD.neg.ctrls, scIdx=groups, k=1)
  # We saved the factors of unwanted variation in a list object
  RUVs1W[[i]] <- ruvedSExprData$W
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

# Sort by genes to associate correctly the Up- and Down-regulated DGE 
df.RUVs1 <- lapply(df.RUVs1, function(u) u[order(rownames(u)),])

# Identify the Up- and Down-regulated DGE 
up.down.regulation <- lapply(resRUVs1, function(u) as.data.frame(decideTests(u)))
for (i in seq_along(up.down.regulation)) {
  up.down.regulation[[i]]$gene <- rownames(up.down.regulation[[i]])
  colnames(up.down.regulation[[i]]) <- c("regulation","gene")
  up.down.regulation[[i]] <- up.down.regulation[[i]][order(rownames(up.down.regulation[[i]])),]
  df.RUVs1[[i]]$regulation <- up.down.regulation[[i]]$regulation
  df.RUVs1[[i]]$regulation[df.RUVs1[[i]]$regulation=="-1"]<- "DOWN"
  df.RUVs1[[i]]$regulation[df.RUVs1[[i]]$regulation=="1"]<- "UP"
  df.RUVs1[[i]]$regulation[df.RUVs1[[i]]$regulation=="0"]<- "NO"
}

# N. of genes detected
n.genes <- lapply(df.RUVs1, function(u) length(rownames(u)))

# data.frame of DGE results 
df.DGE <- lapply(df.RUVs1, function(u) u[u$FDR<0.05,]) 
# N. of DGE
n.DGE <- lapply(df.DGE, function(u) length(rownames(u))) 

df.UpDGE <- lapply(df.RUVs1, function(u) u[u$regulation=="UP",]) 
# N. of DGE Up
n.UpDGE <-  lapply(df.UpDGE, function(u) length(rownames(u))) 

df.DownDGE <- lapply(df.RUVs1, function(u) u[u$regulation=="DOWN",])
# N. of DGE Down 
n.DownDGE <-  lapply(df.DownDGE, function(u) length(rownames(u))) 

# Intersection between dataset genes and positive control genes
n.PosCtrls <- lapply(df.RUVs1, function(u) length(intersect(rownames(u), PosCtrls$Gene_ID))) 
# Intersection between dataset DGE and positive control genes
n.DGE.PosCtrls <- lapply(df.DGE, function(u) length(intersect(rownames(u), PosCtrls$Gene_ID)))
# Percetage of DGE positive controls/ Positive Controls detected
perc.DGE.PosCtrls <- mapply(function(u,y) 100*round(u/y,4), n.DGE.PosCtrls, n.PosCtrls) 

# Intersection between dataset genes and DGE from Fishpond
n.FP <- lapply(df.RUVs1, function(u) length(intersect(rownames(u), DGE_Fishpond$Gene_Stable_ID))) 
# Intersection between dataset DGE and DGE from Fishpond
n.DGE.FP <- lapply(df.DGE, function(u) length(intersect(rownames(u), DGE_Fishpond$Gene_Stable_ID)))
# Percetage of DGE Fishpond/ DGE Fishpond detected
perc.DGE.FP <- mapply(function(u,y) 100*round(u/y,4), n.DGE.FP, n.FP)

n.nuclei <- matrix(pseudo.bulk$ncells,11,6,byrow=TRUE) 
n.nuclei <- rowSums(n.nuclei) # N. of nuclei

# Scatter plot ####
load("AllenColorLabel.RData")
neuronal.color <- subclass.color[-c(1:3,12:14,16,17,19,21)]

tab <- as.data.frame(cbind(n.nuclei, as.matrix(n.DGE)))
colnames(tab) <- c("nuclei","DGE")
rownames(tab) <- levels(factor(pseudo.bulk$Azimuth.labels))
tab$Label <- factor(rownames(tab))
tab$DGE <- as.numeric(tab$DGE)
tab$nuclei <- as.numeric(tab$nuclei)

p1 <- ggplot(tab, aes(x=log(nuclei),y=log(DGE)))+ 
  geom_point(aes(color = Label), size = 7.5 )+ scale_color_manual(values=neuronal.color)+ 
  stat_smooth(method = "lm", formula = y ~ poly(x, 2),se=FALSE,col="grey") +
  theme_classic(base_size = 18)+theme(legend.position="none")+
  geom_label_repel(aes(label = Label), size=5,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50', colour =neuronal.color)
ggsave(p1, file="snRNA_Scatterplot.pdf",  width = 20, height = 20, units = "cm")

# PCA plot ####
label <- levels(factor(pseudo.bulk$Azimuth.labels))

label.save <- gsub("/","", label)
label.save <- gsub(" ","", label.save)
# Calculate the logarithm of normalized counts matrix
normalizedCounts <- lapply(ruvedSExprData, function(u) log1p(u$normalizedCounts))

# Create a SeqExpressionSet object 
seq.RUVs1 <- mapply(function(u,y) 
  newSeqExpressionSet(as.matrix(round(u)), phenoData = data.frame(colData(y),row.names=colnames(y))),
  normalizedCounts, CellType.PseudoBulk)

for (i in seq_along(label)) {
  file.name <- paste("snRNA_PCA_RUVs1_",label.save[i],".pdf",sep = "")
  pdf(file = file.name)
  par(mar =  c(5, 4, 4, 6) + 0.1)
  plotPCA(seq.RUVs1[[i]], xlim = c(-1, 1),label=FALSE,  pch=20, 
          col=as.numeric(as.factor(seq.RUVs1[[i]]$condition)), cex=2, 
          main=paste(label[i],sep=" ") )
  legend("right", inset = c(-0.17,0), legend=levels(factor(seq.RUVs1[[i]]$condition)),  pch=20, cex=1,
         col= levels(as.factor(as.numeric(as.factor(seq.RUVs1[[i]]$condition)))), xpd = TRUE, 
         horiz = FALSE,  bty = "n",title="Condition")
  dev.off()
}

# Volcano plot ####
for (i in 1:length(label)) {
  df.RUVs1[[i]] <- df.RUVs1[[i]][!is.na(df.RUVs1[[i]]$PValue),]
  df.RUVs1[[i]]$Significance <- "No Significant"
  df.RUVs1[[i]]$Significance[df.RUVs1[[i]]$FDR<0.05] <- "Significant"
  
  inter <- intersect(rownames(df.RUVs1[[i]]), PosCtrls$Gene_ID)
  df.RUVs1Pos <-  df.RUVs1[[i]][rownames(df.RUVs1[[i]]) %in% inter,]
  df.RUVs1Pos$Significance[df.RUVs1Pos$Significance=="Significant"] <- "SignificantPos"
  
  df.RUVs1[[i]] <- df.RUVs1[[i]][!rownames(df.RUVs1[[i]]) %in% inter,]
  
  df.RUVs1[[i]] <- rbind(df.RUVs1[[i]], df.RUVs1Pos)
  
  volcano_plot <- ggplot(data=df.RUVs1[[i]], aes(x=logFC, y=-log10(PValue),col=Significance))+
    geom_point(size=2)+
    scale_color_manual(values=c("Significant"="black","No Significant"="grey","SignificantPos"="red"),
                       labels = c(paste("Pval Adj < 0.05 (",n.DGE[[i]],")",sep = ""),"Not Sign.",
                                  paste("Pos. Ctrls < 0.05 (",n.DGE.PosCtrls[[i]],")",sep = "")))+
    ggtitle(paste(label[i],sep = "")) +theme_classic(base_size = 18) + guides(colour = guide_legend(override.aes = list(size=2)))
  ggsave(volcano_plot, file=paste("snRNA_edgeR_VolcanoPlot_",label.save[i],".pdf", sep = "" ) ,width = 20, height = 20, units = "cm")
}

# Histogram of p-value distribution ####
for (i in 1:length(label)) {
  file.name <- paste("snRNA_edgeR_Histogram_pvalue_", label.save[i], ".pdf", sep = "")
  pdf(file.name)
  hist(df.RUVs1[[i]]$PValue, xlab = "p-value", main=paste(label[i],sep=" "))
  dev.off()
}

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

# Sort by genes to associate correctly the Up- and Down-regulated DGE 
df.RUVs1.Sum <- lapply(df.RUVs1.Sum, function(u) u[order(rownames(u)),])

# Identify the Up- and Down-regulated DGE 
up.down.regulation <- as.data.frame(decideTests(resRUVs1.Sum))
up.down.regulation$gene <- rownames(up.down.regulation)
colnames(up.down.regulation) <- c("regulation","gene")
up.down.regulation <- up.down.regulation[order(rownames(up.down.regulation)),]

df.RUVs1.Sum$regulation <- up.down.regulation$regulation
df.RUVs1.Sum$regulation[df.RUVs1.Sum$regulation=="-1"]<- "DOWN"
df.RUVs1.Sum$regulation[df.RUVs1.Sum$regulation=="1"]<- "UP"
df.RUVs1.Sum$regulation[df.RUVs1.Sum$regulation=="0"]<- "NO"

n.genes <- length(rownames(df.RUVs1.Sum)) # N. of genes

df.DEGs <- df.RUVs1.Sum[df.RUVs1.Sum$FDR<0.05,] # data.frame of DEGs results 
n.DEGs <- length(rownames(df.DEGs)) # N. of DEGs

df.UpDEGs <- df.RUVs1.Sum[df.RUVs1.Sum$regulation=="UP",]
n.UpDEGs <-  length(rownames(df.UpDEGs)) # N. of DEGs 

df.DownDEGs <- df.RUVs1.Sum[df.RUVs1.Sum$regulation=="DOWN",]
n.DownDEGs <-  length(rownames(df.DownDEGs)) # N. of DEGs 

n.PosCtrls <- length(intersect(rownames(df.RUVs1.Sum), PosCtrls$Gene_ID)) # Intersection between dataset genes and positive controls
n.DEGs.PosCtrls <- length(intersect(rownames(df.DEGs), PosCtrls$Gene_ID)) # Intersection between dataset DEGs and positive controls
perc.DEGs.PosCtrls <- 100*round(n.DEGs.PosCtrls/n.PosCtrls,4) # Percetage of DEGs positive controls/ Positive Controls detected

n.FP <- length(intersect(rownames(df.RUVs1.Sum), DGE_Fishpond$Gene_Stable_ID)) # Intersection between dataset genes and DEGs from Fishpond
n.DEGs.FP <- length(intersect(rownames(df.DEGs), DGE_Fishpond$Gene_Stable_ID)) # Intersection between dataset DEGs and DEGs from Fishpond
perc.DEGs.FP <- 100*round(n.DEGs.FP/n.FP,4) # Percetage of DEGs Fishpond/ DEGs FP detected

n.nuclei <- 47649
