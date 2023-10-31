# Fishpond Differential Transcript Expression Pipeline

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

#### Differential Transcript Expression with RUVs ####

# ‘se’ will be referred to as ‘y’ for the remainder of the analysis.
y <- se

# I will do a two group comparison between wild-type animals.
# Those that were left undisturbed in their home cage (HC) for 5 hours and those
# who were sleep sleep deprived (SD) via gentle handling for 5 hours

# The log fold changes will be made comparing the second group to the first 
# group listed in levels, in this case, SD compared to HC controls.
y <- y[,y$condition %in% c("WTHC5", "WTSD5")]
y$condition <- factor(y$condition, levels=c("WTHC5", "WTSD5"))

# Call the fishpond package here 
suppressPackageStartupMessages(library(fishpond)) # Version 2.4.1

# Scale inferential replicates to the mean sequencing depth
# It is important to note that only the inferential replicates are scaled
# during this process, not the counts.
y <- scaleInfReps(y)

# Before filtering, we have 140,992 transcripts:
dim(y)
# [1] 140992     10 

# Filter here, keeping a minimum count of 10 transcripts across a minimum of 
# 3 replicates:
y <- labelKeep(y) 

# Following filtering, we have 54,030 transcripts:
y <- y[mcols(y)$keep,]
dim(y)
# [1] 54030    10

# After filtering, we will need to determine the sources of unwanted variation.

# Assemble a matrix that contains the groups
# This is an alternative way to create the matrix: x <- y$condition
groups <- matrix(data = c(1:5, 6:10), nrow = 2, byrow = TRUE)

# Check to make sure 'groups' appears as expected.
groups
#      [,1] [,2] [,3] [,4] [,5]
# [1,]    1    2    3    4    5
# [2,]    6    7    8    9   10

# Assign the names of the variables to match the column names of your data
# This is important for RColorBrewer, which is used to designate colors in
# figures later on in this analysis.
x <- as.factor(rep(c("WTHC5", "WTSD5"), c(5,5)) )
x
#  [1] WTHC5 WTHC5 WTHC5 WTHC5 WTHC5 WTSD5 WTSD5 WTSD5 WTSD5 WTSD5
# Levels: WTHC5 WTSD5

# Here I shorten the names of the samples for simplicity while making the matrix:
names(x) <- c("HC5_1", "HC5_2", "HC5_3", "HC5_4", "HC5_5", "SD5_1", "SD5_2", 
              "SD5_3", "SD5_4", "SD5_5")

# See what type of data class x is:
data.class(x)
# [1] "factor"

# Turn x into a matrix:
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

# Load RColorBrewer here and set the color palate:
suppressPackageStartupMessages(library(RColorBrewer)) # Version 1.1-3

# This palette will be used for PCA and RLE plots.
# To view other palette options, use brewer.pal.info.
colors <- brewer.pal(9, "Greys")[c(6, 7)]
colLib <- colors[x]
HC_color  <- brewer.pal(9, "Greys")[6] 
SD_color <- brewer.pal(9, "Greys")[7] 

# Remove the version #s for each transcript.
# It is helpful to note that the version #s correspond to the release version, 
# they are not different transcripts of the same gene.
rownames(y) <- lapply(rownames(y),  sub, pattern = "\\.\\d+$", replacement = "")
data.class(rownames(y))
# [1] "character"

# Extract background for Functional Annotation with DAVID
# write.table(sQuote(rownames(y)), "080223_DAVID_DTU/080223_Preparing_Background/080223_Transcript_Background.txt", quote = FALSE, 
#                        row.names = FALSE, col.names = FALSE)


# Extract the counts from the SummarizedExperiment object for normalization
# Note that we will estimate the factors of unwanted variation on the counts, 
# and then we will apply that across the inferential replicates further downstream.
counts <- as.matrix((assays(y)[["counts"]]))
data.class(counts)
# [1] "matrix"

# Upper Quartile (UQ) normalization is implemented here with the EDASeq package,
# and accounts for variation in sequencing depth

# As the counts were not scaled (only the inferential replicates), we need
# to account for sequencing depth prior to oroceeding with RUVs normalization,
# which allows us to obtain a the estimated factors of unwanted variation (s$W):
suppressPackageStartupMessages(library("EDASeq")) # Version 2.32.0
uq <- betweenLaneNormalization(counts, which = "upper")
dim(uq)
# [1] 54030    10

# Shorten the colnames for figures:
colnames(uq)
colnames(uq) <- c("HC5_1", "HC5_2", "HC5_3", "HC5_4", "HC5_5", "SD5_1", "SD5_2", 
                  "SD5_3", "SD5_4", "SD5_5")

# Here we plot two quality control plots (RLE and PCA) following UQ normalization:
# PCH 17 (Triangles) --> SD Animals
# PCH 19 (Circles) --> HC Animals (Controls) collected at the same circadian time
plotRLE(uq, col= colLib, outline = FALSE, las = 3, ylim = c(-1, 1), 
        ylab = "Relative Log Expression", cex.axis = 1.3, cex.lab = 1.3)
plotPCA(uq, labels = FALSE, col = colLib, cex = 1.5, cex.axis = 1.3, cex.lab = 1.3, 
        xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), pch = rep(c(19, 17), times = c(5, 5)))

# Load RUVSeq package 
suppressPackageStartupMessages(library(RUVSeq)) # Version 1.32.0

# RUVseq: Here we will estimate a matrix that contains 
# estimated factors of unwanted factors after UQ normalization

# RUVs uses technical replicates or negative controls samples.
# In this instance, as we do not have negative controls for transcripts, we will
# use the rownames of the matrix y after filtering (expressed). For more information
# regarding RUVs, see Risso et al., 2014.

# We have set the numbers of factors of unwanted variation (k) = 4, which 
# optimizes the number of differentially expressed genes and transcripts detected,
# and the positive control recovery (at the gene level) without removing signal

s <- RUVs(x = uq, cIdx = rownames(y), scIdx = groups, k = 4)

# Plot the RLE and PCA plots (quality control figures) again following RUVs 
# normalization 
plotRLE(s$normalizedCounts, col= colLib, outline = FALSE, las = 3, 
        ylim = c(-1, 1), ylab = "Relative Log Expression", cex.axis = 1.2, 
        cex.lab = 1.2)
plotPCA(s$normalizedCounts, labels = FALSE, col = colLib, cex = 1.3, 
        cex.axis = 1.2, cex.lab = 1.2, xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), 
        pch = rep(c(19, 17), times = c(5, 5)))

# Downstream we will use s$W (from RUVs) as the estimation of batch factors.

# infRepIdx lists the inferential replicates (30):
infRepIdx <- grep("infRep",assayNames(y),value=TRUE)

# Save the number of inferential replicates as nreps:
nreps <- length(infRepIdx)
# [1] 30

# Account for continuous variables with removeBatchEffect from limma
# Our samples do not fall into discrete clusters so we will use the following
# approach recommended by the Fishpond developers.

# In short, this is done by directly scaling the estimated counts across
# inferential replicates.

# First inferential replicates are logged as limma requires log-expression
# values for a series of samples. This is done with the assay function from the 
# SummarizedExperiment package.

# Load limma 
suppressPackageStartupMessages(library(limma)) # Version 3.54.2
model_matrix <- model.matrix(~condition, colData(y))

pc <- .1 # This is added to avoid negative InfReps
for (k in seq_len(nreps)) {
  logInfRep <- log(assay(y, infRepIdx[k]) + pc)
  logInfRep <- limma::removeBatchEffect(
    logInfRep,
    covariates=s$W,
    design=model_matrix)
  assay(y, infRepIdx[k]) <- exp(logInfRep)
}

# The Swish method is described in (Zhu et al. 2019).
# The set.seed function allows for reproducibility of exact results in the future.
set.seed(1)
y <- swish(y, x="condition")

# View a table of differently expressed transcripts with a qvalue < 0.05
table(mcols(y)$qvalue < .05)
# FALSE  TRUE 
# 38505 15525 

# You can also view a table that shows which transcripts are up regulated (1) in 
# response to sleep deprivation and which are downregulated (-1)
with(mcols(y),
     table(sig=qvalue < .05, sign.lfc=sign(log2FC))
)

# k = 4
# sign.lfc
# sig        -1     0     1
# FALSE 25479    11 13015
# TRUE   9709     0  5816


# Here is how to make your own infReps matrix for plotting. As we will want to 
# create extra figures comparing counts and proportions, we export counts here: 

# This is to select all features
infReps <- assays(y)[ grep("infRep", assayNames(y)) ]

# abind combines multi-dimensional arrays
infArray <- abind::abind( as.list(infReps), along=3 )

# Take the dimensions of the array (rows (transcripts) x columns (replicates) 
# x heights (replicates))
dim(infArray)
# [1] 54030    10    30

# After selecting all features, you can take the median (or mean):
infMed <- apply(infArray, 1:2, median)
head(infMed)
data.class(infMed)
# [1] "matrix"

dim(infMed)
# [1] 54030    10

# Write tables with the median inferential replicate
# write.table(infMed, file = "071823_DTE_InfMed.txt", sep = "\t")

# A histogram shows the distribution of pvalues:
hist(mcols(y)$pvalue, col="grey", ylim= c(0,25000), main = "", xlab = "Pvalue",
     cex.axis = 0.9)

# We can select just the transcripts with a qvalue < .05:
sig <- mcols(y)$qvalue < .05

# Add symbols to the Summarized Experiment object:
suppressPackageStartupMessages(library(org.Mm.eg.db))
y <- addIds(y, "SYMBOL", gene=TRUE)
rowData(y)

# We turn y into a dataframe here to intersect with lists when making plots:
y_rowData <- as.data.frame(rowData(y))
dim(y_rowData)
# [1] 54030    11

# MA Plot
# An MA plot shows log2FC vs the log10mean
# Transcripts above M=0 are upregulated, while transcripts below are 
# downregulated

# The transcripts that are farthest away from the M=0 axis are the ones most
# affected by the treatment

# Set colors for MA plot
# Transcripts that do not have a significant change in expression are light grey
# Transcripts that are significantly affected (qvalue < 0.05) by SD are dark grey
Significant_Color <- brewer.pal(9, "Greys")[8]

# A subset of transcripts that we have chosen to highlight that also have 
# significant changes in expression (qvalue < 0.05) will be bright blue. 
# Downstream, we will apply an abs(log2FC) threshold of 0.2 to these transcripts. 
Highlight_Color <- "dodgerblue" # Genes from "Highlight Genes" that are significant

# Read the list of genes
# If they have a transcript with a significant change in expression, it will 
# be highlighted
Highlight_Genes <- readxl::read_excel("Peixoto_Figure_2_Supplementary_Table_2.xlsx", sheet = 2,
                                      col_names = FALSE)
# # A tibble: 11 × 1
# ...1    
# <chr>   
# 1 Arc     
# 2 Bdnf    
# 3 Cirbp   
# 4 Eif4ebp1
# 5 Fos     
# 6 Homer1  
# 7 Mef2c   
# 8 Hdac7   
# 9 Hspa5   
# 10 Sst     
# 11 Wnt9a 

# Save as a matrix:
Gene_Names <- as.matrix(Highlight_Genes)
dim(Gene_Names)
# [1] 11  1

# Remove commas:
Gene_Names <- unlist(strsplit(Gene_Names,","))

# Remove quotes:
Gene_Names <- noquote(Gene_Names)

# Remove trailing spaces and save non-duplicated values:
Gene_Names <- trimws(Gene_Names)

# Use plotMASwish to generate the MA plot:
plotMASwish(y, alpha=.05, xlim=c(0,6), ylim=c(-8,5),
            sigcolor = Significant_Color, cex.axis = 1.3, cex.lab = 1.3)

# Include labels for the genes in 'Highlight_Genes'
# Only transcripts that fit the qvalue and log2FC threshold will be labeled
with(
  subset(mcols(y), SYMBOL %in% Gene_Names & qvalue < .05 & abs(log2FC) > 0.2),
  text(log10mean, log2FC, SYMBOL,
       col= Highlight_Color, pos=4, cex=0.8, font=2) 
)


# Labels were plotted above, but now we need to add points to the genes
# We will intersect y_rowData (the dataframe saved above) with the 'Gene Names' matrix
Data_Intersection <- dplyr::filter(y_rowData, SYMBOL %in% Gene_Names & qvalue < 0.05 & abs(log2FC) > 0.2)
dim(Data_Intersection)
# [1] 33 11

# Plot points and squares for 'Data Intersection'
points(Data_Intersection$log10mean, Data_Intersection$log2FC, pch= 0, col = Highlight_Color, cex = 1, lwd = 1)
points(Data_Intersection$log10mean, Data_Intersection$log2FC, pch= 20, col = Highlight_Color, cex = 0.2, lwd = 2)

# Write text file using the code below.
# Select columns to include in the text tile. 
# These lists will then be annotated in Perl with the gene name and description
cols <- c("log10mean","log2FC","pvalue","qvalue")
y_filter <- mcols(y)
y_filter <- print(as.data.frame(y_filter)[,cols], digits=3)
y_filter <- dplyr::filter(y_filter, qvalue < .05)
dim(y_filter)
# [1] 15525     4

# write.table(y_filter, 
#            file = "102722_Differential_Transcript_Expression_Significant_k=4.txt")

# sink('072023_Fishpond_SessionInfo.txt')
# sessionInfo()
# sink() 


