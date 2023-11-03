# Fishpond Differential Gene Expression Pipeline

# Adapted from Michael Love's Lab and the Zhu et al 2019 tutorial linked here:
# https://bioconductor.org/packages/release/bioc/vignettes/fishpond/inst/doc/swish.html#Differential_transcript_usage

#### STEP ONE: Data Import ####

# Establish the working directory
setwd("~/Dropbox/Sleep_RNAseq_Splicing/BulkRNAseq_SD/Final_Tximeta")
dir <- setwd("~/Dropbox/Sleep_RNAseq_Splicing/BulkRNAseq_SD/Final_Tximeta")

# Import the colData (formerly the sample data that will now become the column 
# data), which contains the sample names and condition (the condition is 
# necessary for Fishpond)

# In short, the colData is a file that contains descriptions of the samples,
# and is created during tximeta (Love et al., 2020)

coldata <- read.table(file.path(dir, "coldata.txt"))
head(coldata)
#                   names                 condition
# S3HC5_PFC_1_quant S3HC5_PFC_1_quant     S3HC5
# S3HC5_PFC_2_quant S3HC5_PFC_2_quant     S3HC5
# S3HC5_PFC_3_quant S3HC5_PFC_3_quant     S3HC5
# S3HC5_PFC_4_quant S3HC5_PFC_4_quant     S3HC5
# S3HC5_PFC_5_quant S3HC5_PFC_5_quant     S3HC5
# S3HC7_PFC_1_quant S3HC7_PFC_1_quant     S3HC7

# Add a path to locate the files
coldata$files <- file.path(dir, "data/Salmon_Quants", coldata$names, "quant.sf")
coldata$files 

# Before continuing, make sure that all files in your colData exist in this 
# location
all(file.exists(coldata$files))
# [1] TRUE

# Load ‘SummarizedExperiment’ and ‘tximeta’ packages:
suppressPackageStartupMessages(library(SummarizedExperiment)) # Version 1.28.0
suppressPackageStartupMessages(library(tximeta)) # Version 1.16.1

# I will next load the quant data with tximeta:
se <- tximeta(coldata) 
# importing quantifications
#  reading in files with read_tsv
# 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 
# found matching linked transcriptome:
# [ GENCODE - Mus musculus - release M25 ]
# loading existing TxDb created: 2022-01-13 22:28:02
# Loading required package: GenomicFeatures
# Loading required package: AnnotationDbi
# loading existing transcript ranges created: 2022-01-13 22:28:03
# fetching genome info for GENCODE

# View the ‘countsFromAbundance’ (This should return ‘no’ as we did not use ‘ScaledTPM’ or ‘Length-
# ScaledTPM’, due to concern with how this would affect our downstream normalization with RUVs to correct
# for batch effects):
metadata(se)$countsFromAbundance 

# Make sure assays are loaded. Assays are experimental results and necessary
# for Fishpond.

# You can call assay() to retrieve experimental data

# So what are inferential replicates? During salmon quantification, you can 
# specify "--numBootstraps 30" to generate 30 inferential replicates which are
# effectively re-samplings of samples
assayNames(se)
# [1] "counts"    "abundance" "length"    "infRep1"   "infRep2"   "infRep3"   "infRep4"  
# [8] "infRep5"   "infRep6"   "infRep7"   "infRep8"   "infRep9"   "infRep10"  "infRep11" 
# [15] "infRep12"  "infRep13"  "infRep14"  "infRep15"  "infRep16"  "infRep17"  "infRep18" 
# [22] "infRep19"  "infRep20"  "infRep21"  "infRep22"  "infRep23"  "infRep24"  "infRep25" 
# [29] "infRep26"  "infRep27"  "infRep28"  "infRep29"  "infRep30" 

# View the rownames (These are transcript IDs):
head(rownames(se))
# [1] "ENSMUST00000193812.1" "ENSMUST00000082908.1" "ENSMUST00000162897.1" 
# "ENSMUST00000159265.1" [5] "ENSMUST00000070533.4" "ENSMUST00000192857.1"

#### Differential Gene Expression with RUVs ####
# First we summarize our transcript results to the gene level using the summarizeToGene 
# function from Tximeta
gse <- summarizeToGene(se)
gy <- gse

# For the transcripts we do not have positive or negative controls in response
# to sleep deprivation. While we could potentially annotate our gene lists to 
# include all of the transcripts for that gene, that would assume that every
# transcript responds the same way to sleep deprivation, which is unlikely.

# We do however have a combination of positive/negative controls from microarray 
# and public data that we use for the gene level analysis. These controls are from
# Gerstner et al., 2016.

# After downloading the controls, I then mapped them to the 
# version of gencode that we are using (vM25).

# The first step to intersect the positive controls with the expressed matrix 
# (after filtering) and ensure that they are in the matrix. Note that we need to
# remove the version IDs, as well.

# We will do the same two group comparison as demonstrated with differential 
# transcript expression
gy <- gy[,gy$condition %in% c("WTHC5", "WTSD5")]
gy$condition <- factor(gy$condition, levels=c("WTHC5", "WTSD5"))

# Call the fishpond package here 
suppressPackageStartupMessages(library(fishpond)) # Version 2.4.1

# Here we scale inferential replicates to the mean sequencing depth, filter
# so we are left with samples that are present at least 10 times across 3 or
# more samples
gy <- scaleInfReps(gy)
gy <- labelKeep(gy)

# After filtering, we have 18,334 genes vs the initial 54,347
gy <- gy[mcols(gy)$keep,]
dim(gy)
# [1] 18334    10

# Export DGE background for functional annotation
# write.table(rownames(after_filter), "041723_Gene_Background.txt", quote = FALSE, 
#             row.names = FALSE, col.names = FALSE)
            
# After filtering, we will need to determine the sources of unwanted variation
# We did to assemble a matrix that contains the groups
# This works as well: x <- gy$condition
groups <- matrix(data = c(1:5, 6:10), nrow = 2, byrow = TRUE)

# Check to make sure groups appears as expected
groups
#      [,1] [,2] [,3] [,4] [,5]
# [1,]    1    2    3    4    5
# [2,]    6    7    8    9   10

# Assign the names of the variables to match the column names of your data
# This is important for RColorBrewer, which is what I use to designate colors in
# figures later on in this analysis
x <- as.factor(rep(c("WTHC5", "WTSD5"), c(5,5)) )
x
#  [1] WTHC5 WTHC5 WTHC5 WTHC5 WTHC5 WTSD5 WTSD5 WTSD5 WTSD5 WTSD5
# Levels: WTHC5 WTSD5

# Here I shorten the names of the samples for simplicity while making the matrix
names(x) <- c("HC5_1", "HC5_2", "HC5_3", "HC5_4", "HC5_5", "SD5_1", "SD5_2", 
              "SD5_3", "SD5_4", "SD5_5")

# See what type of data class x is
data.class(x)
# [1] "factor"

# Turn x into a matrix
as.matrix(x)
#         [,1]   
# HC5_1 "WTHC5"
# HC5_2 "WTHC5"
# HC5_3 "WTHC5"
# HC5_4 "WTHC5"
# HC5_5 "WTHC5"
# SD5_1 "WTSD5"
# SD5_2 "WTSD5"
# SD5_3 "WTSD5"
# SD5_4 "WTSD5"
# SD5_5 "WTSD5"

# Load RColorBrewer here and set the color palate
suppressPackageStartupMessages(library(RColorBrewer)) # Version 1.1-3

# I will use this palette for PCA and RLE plots.
# To view other palette options, you can use brewer.pal.info 
colors <- brewer.pal(9, "Greys")[c(6, 7)]
colLib <- colors[x]
HC_color  <- brewer.pal(9, "Greys")[6] 
SD_color <- brewer.pal(9, "Greys")[7] 

# Here we remove the version #s (#s after .) for each gene

# It is helpful to note that the version #s correspond to the release version
rownames(gy) <- lapply(rownames(gy),  sub, pattern = "\\.\\d+$", 
                       replacement = "")
data.class(rownames(gy))
# [1] "character"

# Extract the counts from the SummarizedExperiment object for UQ normalization
# Note that we will estimate the factors of unwanted variation on the counts, and then
# we will apply that across the inferential replicates (see below using Limma)
counts <- as.matrix((assays(gy)[["counts"]]))
data.class(counts)
# [1] "matrix"

# For gene level analysis, we can use a set of positive controls from 
# (Gerstner et al. 2016) to assess the performance of our pipeline downstream. 
# Load the positive controls here:
Gene_Positive_Controls <- read.table("Positive_Controls.txt",
                                       header = TRUE)
dim(Gene_Positive_Controls)
## [1] 675 2

# In preparation for determining our positive control recovery, we will then 
# want to intersect the expressed genes with the positive controls.
Intersect_Positive_Expressed <- intersect(Gene_Positive_Controls[,1], 
                                          row.names(counts))
length(Intersect_Positive_Expressed)
## [1] 671

# Upper Quartile (UQ) Normalization for sequencing depth is implemented here 
# with the EDASeq package
# As the counts were not scaled (only the inferential replicates), we need
# to account for sequencing depth prior to moving on with RUVs normalization,
# which allows us to obtain a the estimated factors of unwanted variation (s$W):
suppressPackageStartupMessages(library("EDASeq"))
g_uq <- betweenLaneNormalization(counts, which = "upper")

# The size of the matrix should not change following normalization:
dim(g_uq)
# [1] 18334    10

# Shorten colnames of counts matrix for figures
colnames(g_uq)
colnames(g_uq) <- c("HC5_1", "HC5_2", "HC5_3", "HC5_4", "HC5_5", "SD5_1", 
                    "SD5_2", "SD5_3", "SD5_4", "SD5_5")

# Here we plot two quality control plots (RLE and PCA) following UQ normalization:
# PCH 17 (Triangles) --> Sleep Deprived Animals
# PCH 19 (Circles) --> Homecage Controls (Non-Sleep Deprived Animals) collected
# at the same circadian time
plotRLE(g_uq, col= colLib, outline = FALSE, las = 3, ylim = c(-1, 1), 
        ylab = "Relative Log Expression", cex.axis = 1.3, cex.lab = 1.3)
plotPCA(g_uq, labels = FALSE, col = colLib, cex = 1.5, cex.axis = 1.3, cex.lab = 1.3, 
        xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), pch = rep(c(19, 17), times = c(5, 5)))

# Read the negative controls, as they will be used with RUVs:
Gene_Negative_Controls <- read.table("Negative_Controls.txt", header = TRUE)
dim(Gene_Negative_Controls)
# [1] 4039    2

# Intersect the negative controls with the rownames:
Negative <- intersect(Gene_Negative_Controls[, 1], rownames(g_uq))
length(Negative)
# [1] 3024

# Load RUVSeq package (Version 1.32.0):
suppressPackageStartupMessages(library(RUVSeq)) 

# RUVseq: Here we will estimate a matrix that contains 
# estimated factors of unwanted factors after UQ normalization

# RUVs uses technical replicates or negative controls. For gene
# level analysis, we have negative controls.

# We have set k = 4, which optimizes the number of differentially expressed
# genes and transcripts detected, and optimizes our positive control recovery
# (at the gene level) without removing signal

# Note that k is the number of factors of unwanted variation that are being estimated from the data
gs <- RUVs(x = g_uq, cIdx = Negative, scIdx = groups, k = 4)

# Plot the RLE and PCA plots, our quality control figures, again following 
# RUVs normalization before applying the estimated factors of variation to all of the inferential replicates
plotRLE(gs$normalizedCounts, col= colLib, outline = FALSE, las = 3, 
        ylim = c(-1, 1), ylab = "Relative Log Expression", cex.axis = 1.2, 
        cex.lab = 1.2)
plotPCA(gs$normalizedCounts, labels = FALSE, col = colLib, cex = 1.3, 
        cex.axis = 1.2, cex.lab = 1.2, xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75),
        pch = rep(c(19, 17), times = c(5, 5)))

# Continue with the Fishpond protocol
infRepIdx <- grep("infRep",assayNames(gy),value=TRUE)
nreps <- length(infRepIdx)

# Account for continuous variables with removeBatchEffect from limma
# Our samples do not fall into discrete clusters so we will use the following
# approach recommended by the Fishpond developers

# In short, this is done by directly scaling the estimated counts across
# inferential replicates

# First inferential replicates are logged as limma requires log-expression
# values for a series of samples. This is done with the assay function from the 
# SummarizedExperiment package

# We use the gs$W as input for covariates (what needs to be adjusted for)
# Design is the conditions (experimental factor we want) other than batch effect
suppressPackageStartupMessages(library(limma))
model_matrix <- model.matrix(~condition, colData(gy))

pc <- .1 # This is added to avoid negative InfReps
for (k in seq_len(nreps)) {
  logInfRep <- log(assay(gy, infRepIdx[k]) + pc)
  logInfRep <- limma::removeBatchEffect(
    logInfRep,
    covariates=gs$W,
    design=model_matrix)
  assay(gy, infRepIdx[k]) <- exp(logInfRep)
}

# The Swish method is described in (Zhu et al. 2019).
# The set.seed function allows for reproducibility of exact results in the future.
set.seed(1)
gy <- swish(gy, x="condition")

# View differential expressed genes in the table below
table(mcols(gy)$qvalue < .05)
# k = 4
# FALSE  TRUE 
# 9829  8505 

with(mcols(gy),
     table(sig=qvalue < .05, sign.lfc=sign(log2FC))
)
#        sign.lfc
# sig       -1    1
# FALSE 5456 4373
# TRUE  4364 4141

significant_qvalue <- (mcols(gy)[mcols(gy)$qvalue < .05,])
dim(significant_qvalue)
# k = 4
# [1] 8505    9

# Intersect qvalue < .05 with intersect positive to figure out positive control recovery
Positive_Controls_Recovery <- length(intersect(row.names(significant_qvalue), Intersect_Positive_Expressed))
# [1] 558 for k = 4

(length(intersect(row.names(significant_qvalue), Intersect_Positive_Expressed))/length(Intersect_Positive_Expressed)) *100
# 83.16% for k = 4

# write.table(intersect(row.names(significant_qvalue), Intersect_Positive_Expressed), 
#           "081123_Positive_Control_Recovery_Genes.txt",
#           sep = "\t",
#           row.names = FALSE,
#           col.names = FALSE)

# Histogram here:
hist(mcols(gy)$pvalue, col="grey", ylim= c(0,10000), main = "", xlab = "Pvalue",
     cex.axis = 0.9)

# Add Gene IDs for MA plot:
suppressPackageStartupMessages(library(org.Mm.eg.db))
gy <- addIds(gy, "SYMBOL", gene=TRUE)

# Turn ‘gy’ into a dataframe here, and save as a different variable. 
# We do this to intersect with other lists when making plots, such as genes 
# to highlight on the MA plot downstream. Check the ‘dim’ to ensure that
# it has not changed:
gy_rowData <- as.data.frame(rowData(gy))
dim(gy_rowData)

# Select Colors for MA plot
Significant_Color <- brewer.pal(9, "Greys")[8]  # Genes with a qvalue < 0.05 will be dark grey
Highlight_Color <- "dodgerblue" # Genes from "Highlight Genes" that are significant
# will be bright blue. We will apply an abs(log2FC) threshold of 0.2 
Positive_Color <- brewer.pal(9, "Set1")[1] # Positive controls

# Read the list of genes to be highlighted (Arc, Bdnf, Cirbp, Eif4ebp1, Fos, Homer1, 
# Mef2c, Hdac7, Hspa5, Wnt9a, Sst)
Highlight_Genes <- readxl::read_excel("Peixoto_Figure_2_Supplementary_Table_2.xlsx", sheet = 2, col_names = FALSE)

Gene_Names <- as.matrix(Highlight_Genes)
dim(Gene_Names)
# [1] 11  1
Gene_Names <- unlist(strsplit(Gene_Names,","))
Gene_Names <- noquote(Gene_Names)
# Remove trailing spaces
Gene_Names <- trimws(Gene_Names)

# Use ‘plotMASwish’ to generate the MA plot. 
# Include labels and points for the genes in ‘Highlight_Genes’
# and points for positive controls using the code below 
# (Only genes that fit the qvalue and log2FC threshold will be labeled):
plotMASwish(gy, alpha=.05, xlim=c(.5,5.5), ylim=c(-7,5), sigcolor = Significant_Color,
            cex.axis = 1.3, cex.lab = 1.3)

with(
  subset(mcols(gy), SYMBOL %in% Gene_Names & qvalue < .05 & abs(log2FC) > 0.2),
  text(log10mean, log2FC, SYMBOL,
       col= Highlight_Color, pos=4, cex=0.8, font=2) #second was 2 before 1
)

# Add the positive controls to the plot (red)

Positive_Controls_Signifcant <- intersect(row.names(significant_qvalue), Intersect_Positive_Expressed)
Positive_Controls_Plot <- mcols(gy)[ rownames(mcols(gy)) %in% Positive_Controls_Signifcant, ]
dim(Positive_Controls_Plot)
# [1] 558  10
points(x= Positive_Controls_Plot$log10mean, y = Positive_Controls_Plot$log2FC, pch = 20, col = Positive_Color, cex = .2, lwd = 2)

# Plot subset of highlighted genes (such as Arc, Bdnf, etc) with a qvalue < 0.05 and log2FC > 0.2(names)

# Filter so anly highlighted genes with a qvalue less than 0.05 and fold change greater than 0.2 show up
Highlight_Genes_Intersection <- dplyr::filter(gy_rowData, SYMBOL %in% Gene_Names & qvalue < 0.05 & abs(log2FC) > 0.2)
dim(Highlight_Genes_Intersection)
# [1] 10 10
sort(Highlight_Genes_Intersection$SYMBOL)

# Add points for highlight genes
points(Highlight_Genes_Intersection$log10mean, Highlight_Genes_Intersection$log2FC, pch= 0, col = Highlight_Color, cex = 1, lwd = 1)
points(Highlight_Genes_Intersection$log10mean, Highlight_Genes_Intersection$log2FC, pch= 20, col = Highlight_Color, cex = 0.2, lwd = 2)

# Just include these columns to print out list of significant differential gene expression
cols <- c("log10mean","log2FC","pvalue","qvalue")
gy <- mcols(gy)
gy <- print(as.data.frame(gy)[,cols], digits=3)
gy <- dplyr::filter(gy, qvalue < .05)
dim(gy)
# [1] 8505    4

# write.table(gy, file = "DGE_Significant_k=4.txt")

# sink('072023_Fishpond_SessionInfo.txt')
# sessionInfo()
# sink() 



