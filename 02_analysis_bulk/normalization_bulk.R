# This file explores the bulk analysis data from WT animals for the following conditions:
# Home Cage 5 Hours (HC5), Home Cage 7 Hours (HC7), Sleep Deprivation 5 Hours (SD5),
# Recovery Sleep 2 Hours (RS2)

# RStudio Version 1.3.1903


# Downloaded Bioconductor Packages:
library(limma)
library(edgeR)
library(EDASeq)
library(RUVSeq)
library(ffpe)
library(RColorBrewer)

# set the working directory
setwd("~/R/data_processing")

# First, read the SD counts in R:
sdrs <- read.table("SD_RS_WT_salmon_gene_level_data.txt", row.names = 1, header = TRUE)
head(sdrs)

# See how many genes we are starting with:
dim(sdrs)
# [1] 54347    20

# Next, read the the positive control file:
positive <- read.table("Peixoto_SD_RS_Positive_Controls.txt", as.is = TRUE, sep = "\t", header = TRUE )
x <- as.factor(rep(c("HC5", "HC7", "SD5", "RS2"), each= 5))
names(x) <- colnames(sdrs)

# Then, filter out non expressed genes:
filter <- apply(sdrs, 1, function(x) length(x[which(x > 10)]) > 5)
filtered <- as.matrix(sdrs)[filter, ]
dim(filtered)
# [1] 18084    20

# The positive control file contains SD and RS positive controls
# There are 558 SD positive controls
SDup <- intersect(positive[positive[, 3] == "UP", 1], rownames(filtered))
SDdown <- intersect(positive[positive[, 3] == "DOWN", 1], rownames(filtered))

# RS positive controls
RSup <- intersect(positive[positive[, 4] == "UP", 1], rownames(filtered))
RSdown <- intersect(positive[positive[, 4] == "DOWN", 1], rownames(filtered))

colors <- brewer.pal(9, "Set1")
colLib <- colors[x]

# Upper Quartile Normalization (divides the read count by the 75th percentile of the read counts in its sample)
uq <- betweenLaneNormalization(filtered, which = "upper")
dim(uq)
# [1] 18084    20

# For negative controls, read the rownames of the uq data:
rownames(uq)
negCon <- rownames(uq)

# Generating RLE Plots: RLE plots reveal confounders when the mean and the variance are not similar.
plotRLE(uq, col= colLib, outline = FALSE, las = 3, ylim = c(-0.5, 0.5), ylab = "Relative Log Expression", cex.axis = 1, cex.lab = 1)

# Generating PCA Plots: We need to make sure that replicates are more similar to each other than they are to samples from other conditions.
plotPCA(uq, col = colLib, cex = 1, cex.axis = 1, cex.lab = 1, xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75))

# RUV Normalization with k=5.
groups <- matrix(data = c(1:5, 6:10, 11:15, 16:20), nrow = 4, byrow = TRUE)
s <- RUVs(x=uq, cIdx=rownames(uq), scIdx=groups, k= 5)

# Following RUV Normalization, plot RLE and PCA: 
plotRLE(s$normalizedCounts, col = colLib, outline = FALSE, las = 3, ylim = c(-0.5, 0.5), ylab= "Relative Log Expression", cex.axis = 1, cex.lab = 1)
plotPCA(s$normalizedCounts, col = colLib, cex = 1, cex.axis = 1, cex.lab = 1, xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75))

## Differential Expression: UQ Normalization

design <- model.matrix(~x) 
design
##    (Intercept) xHC7 xRS2 xSD5
## 1            1    0    0    0
## 2            1    0    0    0
## 3            1    0    0    0
## 4            1    0    0    0
## 5            1    0    0    0
## 6            1    1    0    0
## 7            1    1    0    0
## 8            1    1    0    0
## 9            1    1    0    0
## 10           1    1    0    0
## 11           1    0    0    1
## 12           1    0    0    1
## 13           1    0    0    1
## 14           1    0    0    1
## 15           1    0    0    1
## 16           1    0    1    0
## 17           1    0    1    0
## 18           1    0    1    0
## 19           1    0    1    0
## 20           1    0    1    0
## attr(,"assign")
## [1] 0 1 1 1
## attr(,"contrasts")
## attr(,"contrasts")$x
## [1] "contr.treatment"

y <- DGEList(counts = filtered, group = x)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design, verbose = TRUE)
## Disp = 0.01037 , BCV = 0.1018 

y <- estimateGLMTagwiseDisp(y, design) 
fit <- glmFit(y, design)
dim(y)

# SD Data 
lrt <- glmLRT(fit, coef = 4) 
topUQSD <- topTags(lrt, n = Inf)$table
topTags(lrt) # the results of de
topTags <- topTags(lrt)
topTags
save(topTags, file= ("'edgeRresult.txt'"))
('edgeRresult.txt')

## Histogram
hist(topUQSD$PValue, main = "", xlab = "p-value", breaks = 100, ylim = c(0, 1400))

## Volcano Plot
plot(topUQSD[, 1], -log10(topUQSD$PValue), pch = 20, col = "gray", cex= 0.5, ylab = "-log10(p-value)", xlab= "log2(SD/HC5)", ylim = c(0, 15), xlim = c(-3, 3), cex.lab = 1, cex.axis = 1)

# de is blue
de <- rownames(topUQSD[topUQSD$FDR <= 0.05, ]) # In eLife it was less than 0.05
points(topUQSD[de, 1], -log10(topUQSD[de, "PValue"]), pch = 20, col = colors[2], cex = 1, lwd = 2)

# positive controls are red
points(topUQSD[SDup, 1], -log10(topUQSD[SDup, "PValue"]), pch = 1, col = colors[1], cex = 1, lwd = 2)
points(topUQSD[SDdown, 1], -log10(topUQSD[SDdown, "PValue"]), pch = 1, col = colors[1], cex = 1, lwd = 2)

# negative controls are green
points(topUQSD[negCon, 1], -log10(topUQSD[negCon, "PValue"]), pch = 1, col = colors[3], cex = 1, lwd = 2)


### Differential Expression: RUV Normalization
design <- model.matrix(~x + s$W)
# View the matrix for k=5
design 
## (Intercept) xHC7 xRS2 xSD5    s$WW_1     s$WW_2     s$WW_3
## 1            1    0    0    0 0.5619589 -1.1265810 -0.7213212
## 2            1    0    0    0 0.8892309 -1.4192080 -1.0486684
## 3            1    0    0    0 0.7838863 -1.3837141 -0.8837692
## 4            1    0    0    0 1.0436993 -1.5700778 -0.8458786
## 5            1    0    0    0 1.1782482 -1.3692613 -0.9283006
## 6            1    1    0    0 0.7870349 -1.2961416 -0.8898897
## 7            1    1    0    0 0.7963250 -1.3819239 -0.9456760
## 8            1    1    0    0 0.9492631 -1.3838981 -1.2088477
## 9            1    1    0    0 1.2265329 -1.6994939 -0.7327716
## 10           1    1    0    0 1.1159228 -1.1435857 -0.6415171
## 11           1    0    0    1 0.8595869 -1.0839891 -1.3696074
## 12           1    0    0    1 0.8194187 -1.4033900 -1.0873567
## 13           1    0    0    1 0.7273045 -1.7320783 -0.4465339
## 14           1    0    0    1 1.2578119 -1.2874312 -0.7240840
## 15           1    0    0    1 1.1762271 -1.1092751 -0.5186372
## 16           1    0    1    0 0.4192706 -1.1554371 -0.4934086
## 17           1    0    1    0 1.0420752 -0.9451563 -0.9536648
## 18           1    0    1    0 0.9090918 -1.3456205 -0.9235614
## 19           1    0    1    0 1.0327911 -1.8111225 -0.9882471
## 20           1    0    1    0 1.3069646 -0.9689444 -0.6364333
## s$WW_4    s$WW_5
## 1   0.2001352731 -2.475223
## 2   0.2135928721 -2.861806
## 3   0.1644966958 -2.639608
## 4   0.0673295216 -2.620855
## 5   0.0241621760 -2.615410
## 6   0.1766663771 -2.745836
## 7  -0.0545126885 -2.689278
## 8  -0.1331174113 -2.626316
## 9   0.3301063610 -2.485010
## 10  0.2440193105 -2.779906
## 11 -0.0213078369 -2.323891
## 12 -0.0849982781 -2.682953
## 13 -0.4738331442 -2.899057
## 14  0.5458763291 -2.494891
## 15  0.0280037612 -2.767487
## 16  0.3234078719 -2.189663
## 17  0.1036001600 -3.152865
## 18  0.0002425417 -2.592624
## 19  0.2177694952 -2.406248
## 20 -0.4745418181 -2.070735
## attr(,"assign")
## [1] 0 1 1 1 2 2 2 2 2
## attr(,"contrasts")
## attr(,"contrasts")$x
## [1] "contr.treatment"

y <- DGEList(counts = filtered, group = x)
y <- calcNormFactors(y, method = "upperquartile")

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 4)
topSD <- topTags(lrt, n = Inf)$table

## Histogram
hist(topSD$PValue, main = "", xlab= "p-value", breaks = 100, ylim = c(0, 1400))

## Volcano Plot
plot(topSD[, 1], -log10(topSD$PValue), pch = 20, col = "gray", cex = 0.5, ylab = "-log10(p-value)", xlab = "log2(SD/HC5)", ylim = c(0, 15), xlim = c(-4, 4), cex.lab = 1, cex.axis = 1)

# de genes are blue
de <- rownames(topSD[topSD$FDR <= 0.05, ])
points(topSD[de, 1], -log10(topSD[de, "PValue"]), pch = 20, col = colors[2], cex = 1)

# positive controls are red
points(topSD[SDup, 1], -log10(topSD[SDup, "PValue"]), pch = 20, col = colors[1], cex =1, lwd = 2)
points(topSD[SDdown, 1], -log10(topSD[SDdown, "PValue"]), pch = 1, col = colors[1], lwd = 2)

# negative controls are green
points(topSD[negCon, 1], -log10(topSD[negCon, "PValue"]), pch = 1, col = colors[3], cex = 1, lwd = 2)

## Determine the percentage of positive control genes detected
((sum(topSD[SDup, "FDR"] < .05) + sum(topSD[SDdown, "FDR"] < .05))/558) * 100
# [1] 49.10394

## The number of non de genes
ndifexp <- sum(topSD$FDR >= 0.05)
ndifexp
# [1] 13345

## The number of de genes
difexp <- sum(topSD$FDR < 0.05)
difexp
# [1] 4739
