# Fishpond Differential Transcript Usage Pipeline

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

#### Differential Transcript Usage with RUVs ####
# Clear environment and return to STEP ONE: Data Import before continuing
y <- se

# I will do a two group comparison between WILD TYPE home cage 5 hours and sleep 
# deprivation 5 hours

# The log fold changes will be made comparing the second group to the first 
# group listed in levels, in this case, sleep deprivation compared to home cage
y <- y[,y$condition %in% c("WTHC5", "WTSD5")]
y$condition <- factor(y$condition, levels=c("WTHC5", "WTSD5"))

# I call fishpond package (2.4.1)
suppressPackageStartupMessages(library(fishpond))

# Scale inferential replicates to the mean sequencing depth. Only the inferential
# replicates are scaled, not the counts.
y <- scaleInfReps(y)

# Before filtering, we have 140,992 transcripts
dim(y)
# [1] 140992     10  

# Filter here, keeping a minimum count of 10 transcripts across a minimum of 3 replicates
y <- labelKeep(y, minCount = 10, minN = 3) 

# Following filtering, we have 54,030 transcripts
y <- y[mcols(y)$keep,]
dim(y)
# [1] 54030    10

# Apply secondary log10mean filter here.
# We will only keep transcripts that have a log10mean > 1. We apply a 
# more stringent filter for DTU analysis, as we want to ensure that we
# are only discovering transcripts that can be reproduced in future analyses.
y <- y[rowData(y)$log10mean > 1]

dim(y)
# [1] 50227    10

# After filtering, we will need to determine the sources of unwanted variation
# We did to assemble a matrix that contains the groups
# This works as well: x <- y$condition
groups <- matrix(data = c(1:5, 6:10), nrow = 2, byrow = TRUE)

# Check to make sure groups appears as expected
groups
#      [,1] [,2] [,3] [,4] [,5]
# [1,]    1    2    3    4    5
# [2,]    6    7    8    9   10

# Assign the names of the variables to match the column names of your data
# This is important for RColorBrewer, which is used to designate colors in
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
# [1] "character"
# colnames(assays(y)[["counts"]])

# Load RColorBrewer here and set the color palate
suppressPackageStartupMessages(library(RColorBrewer)) # Version 1.1-3

# I will use this palette for PCA and RLE plots.
# To view other palette options, you can use brewer.pal.info 
colors <- brewer.pal(9, "Greys")[c(6, 7)]
colLib <- colors[x]
HC_color  <- brewer.pal(9, "Greys")[6] 
SD_color <- brewer.pal(9, "Greys")[7] 

# Here we remove the version #s (#s after .) for each transcript

# It is helpful to note that the version #s correspond to the release version, 
# they are not different transcripts of the same gene
rownames(y) <- lapply(rownames(y),  sub, pattern = "\\.\\d+$", replacement = "")
data.class(rownames(y))
# [1] "character"

# Extract the counts from the SummarizedExperiment object for UQ normalization
# Note that we will estimate the factors of unwanted variation on the counts, and then
# we will apply that across the inferential replicates (see below using Limma)
counts <- as.matrix((assays(y)[["counts"]]))

dim(counts)
# [1] 50227    10

data.class(counts)
# [1] "matrix"

# Upper Quartile (UQ) Normalization for sequencing depth is implemented here 
# with the EDASeq package
# As the counts were not scaled (only the inferential replicates), we need
# to account for sequencing depth prior to moving on with RUVs normalization,
# which allows us to obtain a the estimated factors of unwanted variation (s$W)
suppressPackageStartupMessages(library("EDASeq"))
uq <- betweenLaneNormalization(counts, which = "upper")
dim(uq)
# [1] 50227    10

# Shorten colnames of counts matrix for figures
colnames(uq)
colnames(uq) <- c("HC5_1", "HC5_2", "HC5_3", "HC5_4", "HC5_5", "SD5_1", "SD5_2", 
                  "SD5_3", "SD5_4", "SD5_5")

# Here we plot two quality control plots (RLE and PCA) following UQ normalization:
# PCH 17 (Triangles) --> Sleep Deprived Animals
# PCH 19 (Circles) --> Homecage Controls (Non-Sleep Deprived Animals) collected
# at the same circadian time
plotRLE(uq, col= colLib, outline = FALSE, las = 3, ylim = c(-1, 1), 
        ylab = "Relative Log Expression", cex.axis = 1.3, cex.lab = 1.3)
plotPCA(uq, labels=FALSE, col = colLib, cex = 1.5, cex.axis = 1.3, cex.lab = 1.3, 
        xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), pch = rep(c(19, 17), times = c(5, 5)))


# Load RUVSeq package (Version 1.32.0)
suppressPackageStartupMessages(library(RUVSeq)) 

# RUVseq: Here we will estimate a matrix that contains 
# estimated factors of unwanted factors after UQ normalization

# RUVs uses technical replicates or negative controls samples.
# In this instance, as we do not have negative controls for transcripts, we will
# use the rownames of the matrix y after filtering

# We have set k = 4, which optimizes the number of differentially expressed
# genes and transcripts detected, and optimizes our positive control recovery
# (at the gene level) without removing signal

# Note that k is the number of factors of unwanted variation that are being estimated from the data
s <- RUVs(x = uq, cIdx = rownames(y), scIdx = groups, k = 4)

# Plot the RLE and PCA plots, our quality control figures, again following 
# RUVs normalization before applying the estimated factors of variation to all of the inferential replicates
plotRLE(s$normalizedCounts, col= colLib, outline = FALSE, las = 3, 
        ylim = c(-1, 1), ylab = "Relative Log Expression", cex.axis = 1.3, 
        cex.lab = 1.3)
plotPCA(s$normalizedCounts, labels = FALSE, col = colLib, cex = 1.5, cex.axis = 1.3, 
        cex.lab = 1.3, xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), pch = rep(c(19, 17), times = c(5, 5)))

# As mentioned above will use s$W from RUVs as the estimation of batch factors that we will
# use as we continue with the Fishpond pipeline

# Continue with the Fishpond pipeline
infRepIdx <- grep("infRep",assayNames(y),value=TRUE)
nreps <- length(infRepIdx)

# Account for continuous variables with removeBatchEffect from limma
# Our samples do not fall into discrete clusters so we will use the following
# approach recommended by the Fishpond developers

# In short, this is done by directly scaling the estimated counts across
# inferential replicates

# First inferential replicates are logged as limma requires log-expression
# values for a series of samples. This is done with the assay function from the 
# SummarizedExperiment package

# Load limma (Version 3.54.2)
suppressPackageStartupMessages(library(limma))
model_matrix <- model.matrix(~condition, colData(y))

pc <- .1
for (k in seq_len(nreps)) {
  logInfRep <- log(assay(y, infRepIdx[k]) + pc)
  logInfRep <- limma::removeBatchEffect(
    logInfRep,
    covariates=s$W,
    design=model_matrix)
  assay(y, infRepIdx[k]) <- exp(logInfRep)
}

# isoformProportions takes the output of the scaled/filtered counts and returns proportions
iso <- isoformProportions(y)

# The Swish method is described in (Zhu et al. 2019). Note that the set.seed function allows for the reproducibility
# of exact results in the future.
set.seed(1)
y <- swish(iso, x="condition")

# You can view a table of differently expressed transcripts with a qvalue < 0.05
table(mcols(y)$qvalue < .05)
# k = 4, log10mean > 1
# FALSE  TRUE 
# 41383  2314 

# You can also view a table that shows which transcripts are up regulated in 
# response to sleep deprivation (1) and which are downregulated (-1)
with(mcols(y),
     table(sig=qvalue < .05, sign.lfc=sign(log2FC))
)

# k = 4, log10mean > 1
#       sign.lfc
# sig        -1     1
# FALSE 22829 18554
# TRUE    899  1415

# Here is how to make your own infReps matrix for plotting. As we will want to create extra figures,
# so we export it here: I checked this by reviewing this site:
# https://support.bioconductor.org/p/p134531/

# This is to view one isoform at once, I have selected Homer1a as an example below:
# Homer1a_DTU <- unlist( assays(y["ENSMUST00000102752",])[ grep("infRep", assayNames(y)) ] )

# This is to select all features
infReps <- assays(y)[ grep("infRep", assayNames(y)) ]

# abind combines multi-dimensional arrays
infArray <- abind::abind( as.list(infReps), along=3 )

# Take the dimensions of the array (rows (transcripts) x columns (replicates) 
# x heights (replicates))
dim(infArray)
# [1] 43697    10    30

# After selecting all features, you can take the median or the mean:
infMed <- apply(infArray, 1:2, median)
head(infMed)
data.class(infMed)
# [1] "matrix"

infMean <- apply(infArray, 1:2, mean)
head(infMean)
data.class(infMean)
# [1] "matrix"

# Write tables with the median and mean inferential replicate
# write.table(infMed, file = "072023_DTU_InfMed.txt", sep = "\t")
# write.table(infMean, file = "012623_DTU_InfMean.txt", sep = "\t")

# Histogram here:
hist(mcols(y)$pvalue, col="grey", ylim= c(0,10000))

# You can use the org.Mm.eg.db package to add symbol IDs
suppressPackageStartupMessages(library(org.Mm.eg.db))

#  Once symbol IDs are added you can annotate transcripts in your MA plot
y <- addIds(y, "SYMBOL", gene=TRUE)

# Beyond this is currently a work in progress
# MA Plot
Significant_Color <- brewer.pal(9, "Greys")[8] 
Homer1_Color <- "red" 
Bdnf_Color <- "magenta" 

plotMASwish(y, alpha=.05, xlim=c(.5,5.5), sigcolor = Significant_Color, 
            cex.axis = 1.3, cex.lab = 1.3, ylab = "log2FC (Proportion)")
with(
  subset(mcols(y), SYMBOL == "Homer1" & qvalue < .05),
  text(log10mean, log2FC, SYMBOL,
       col= Homer1_Color, pos=c(4,1,4), cex=1, font=2) #second was 2 before 1
)

with(
  subset(mcols(y), SYMBOL == "Bdnf" & qvalue < .05),
  text(log10mean, log2FC, SYMBOL,
       col= Bdnf_Color, pos=c(3,3), cex=1, font=2)
)

# Load dplyr as we will need this to filter our data to just keep Homer1
suppressPackageStartupMessages(library(dplyr))

# dplyr does not accept summarized experiments so we save the rowdata
# as a data frame
y_rowData <- as.data.frame(rowData(y))

# Here I  filter the data frame so that only Homer1/Bdnf isoforms remain, and
# specifically those  with a qvalue < 0.05
Homer1_DTU <- dplyr::filter(y_rowData, SYMBOL == "Homer1" & qvalue < 0.05)
Bdnf_DTU <- dplyr::filter(y_rowData, SYMBOL == "Bdnf" & qvalue < 0.05)

# I add points to the MA plot for Homer 1
points(Homer1_DTU$log10mean, Homer1_DTU$log2FC, pch= 20, col = Homer1_Color, cex = .2, lwd = 2) # Homer1
points(Homer1_DTU$log10mean, Homer1_DTU$log2FC, pch= 0, col = Homer1_Color, cex = 1.5, lwd = 1.7) # Homer1
points(Bdnf_DTU$log10mean, Bdnf_DTU$log2FC, pch= 20, col = Bdnf_Color, cex = .2, lwd = 2) # Bdnf
points(Bdnf_DTU$log10mean, Bdnf_DTU$log2FC, pch= 0, col = Bdnf_Color, cex = 1.5, lwd = 1.7) # Bdnf

# Limit columns to the following when creating a text file
cols <- c("log10mean","log2FC","pvalue","qvalue")
y_filter <- mcols(y)
y_filter <- print(as.data.frame(y_filter)[,cols], digits=3) # Limit to 3 digits
y_filter <- dplyr::filter(y_filter, qvalue < .05)
dim(y_filter)
# [1] 2314    4

# write.table(y_filter, 
#           file = "072123_Differential_Transcript_Usage_Significant_k=4_log10mean1.txt")

# sink('072023_Fishpond_SessionInfo.txt')
# sessionInfo()
# sink() 


