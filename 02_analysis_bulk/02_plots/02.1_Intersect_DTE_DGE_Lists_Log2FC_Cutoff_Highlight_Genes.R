#### Intersecting Highlight Genes with DTE and DGE Venn Diagrams ####

# Summary: I will intersect transcript and gene expression lists so that we can 
# determine the following:

## What genes can be detected with both gene and transcript expression analyses?

## What genes from the literature can only be detected at the transcript level or
## only at the gene level?

## Additionally, if we apply a log2FC filter, what is the recovery of positive and negative
## controls like? (Extended Data Figure 10)

# R Version 4.2.2

#### Load packages and set the working directory ####
library("readxl") # 1.4.2
library("dplyr") # 1.1.1

setwd("~/Dropbox/Sleep_RNAseq_Splicing/BulkRNAseq_SD/Final_Tximeta")

#### READ LISTS ####
# Read lists without a fold change cutoff
# Gene
Gene <- readxl::read_excel("060123_Peixoto_Figure1_Supplement2.xlsx", sheet = "DGE")

dim(Gene)
# [1] 8505    8

# Transcript Expression
Transcript <- readxl::read_excel("060123_Peixoto_Figure1_Supplement2.xlsx", sheet = "DTE")

dim(Transcript)
# [1] 15525    10

#### Gene Expression, Filter Up- and Down-Regulated ####
Gene_Up <- dplyr::filter(Gene, log2FC > 0)
dim(Gene_Up)
# [1] 4141    8

Gene_Down <- dplyr::filter(Gene, log2FC < 0)
dim(Gene_Down)
# [1] 4364    8

#### Transcript Expression, Filter Up- and Down-Regulated ####
Transcript_Up <- dplyr::filter(Transcript, log2FC > 0)
dim(Transcript_Up)
# [1] 5816   10

Transcript_Down <- dplyr::filter(Transcript, log2FC < 0)
dim(Transcript_Down)
# [1] 9709   10

#### Intersect gene and transcript lists before filtering log2FC ####
#### Transcript Up Only ####
# Transcript Up Only (What is upregulated at the transcript level, but is 
# not upregulated at the gene level?)

# Here we will use anti_join to return all rows from x (Transcript Up) 
# withOUT a match in y (Gene Up)

Transcript_Up_Only <- anti_join(Transcript_Up, Gene_Up, by = "Gene_Stable_ID")
# A tibble: 1,536 × 10

# One way to double check is to see that this matches the #s that we get with 
# the Venn Diagram code (We expect 1381 Gene IDs that are only detected at the
# transcript level)
Transcript_Up_Only_Genes <- Transcript_Up_Only$Gene_Stable_ID
Transcript_Up_Only_Genes <- Transcript_Up_Only_Genes[!duplicated(Transcript_Up_Only_Genes)]
length(Transcript_Up_Only_Genes)
# [1] 1381

# Now, here we will apply a log2FC cutoff, only keeping transcripts that
# have a log2FC > 0.2
Transcript_Up_Only_0.2 <- dplyr::filter(Transcript_Up_Only, log2FC > 0.2)
dim(Transcript_Up_Only_0.2)
# [1] 1025   10

# Now we will remove Gene IDs that are duplicated
Transcript_Up_Only_0.2_Genes <- dplyr::distinct(Transcript_Up_Only_0.2, Gene_Stable_ID, .keep_all = TRUE)
dim(Transcript_Up_Only_0.2_Genes)
# [1] 1381 (log2FC > 0)
# [1] 1272 (log2FC > 0.1)
# [1] 920 (log2FC > 0.2)
# [1] 638 (log2FC > 0.3)
# [1] 412 (log2FC > 0.4)
# [1] 272 (log2FC > 0.5)

#### Gene Up Only ####
# Gene Up # (What is upregulated at the gene level, but is 
# not upregulated at the transcript level?)

# Here we will use anti_join to return all rows from x (Gene Up) 
# withOUT a match in y (Transcript Up)

Gene_Up_Only <- anti_join(Gene_Up, Transcript_Up, by = "Gene_Stable_ID")
# A tibble: 872 × 8

# We don't need to remove duplicates from gene level results
Gene_Up_Only_0.2_Genes <- dplyr::filter(Gene_Up_Only, log2FC > 0.2)
dim(Gene_Up_Only_0.2_Genes)
# [1] 872   8 (log2FC > 0)
# [1] 691   8 (log2FC > 0.1)
# [1] 381   8 (log2FC > 0.2)
# [1] 219   8 (log2FC > 0.3)
# [1] 137   8 (log2FC > 0.4)
# [1] 86  8 (log2FC > 0.5)

#### Transcript and Gene Both Up ####
# Here we identify which genes are detected by both transcript and gene 
# expression analysis

# Inner join will keep all rows in x (Transcript_Up) and y (Gene_Up)
Transcript_Gene_Up <- inner_join(Transcript_Up, Gene_Up, by = "Gene_Stable_ID")
dim(Transcript_Gene_Up)
# [1] 4280   17

# View length after removing duplicates. Why? This is a sanity check to make
# sure what we detect here matches what we see with the Venn Diagram Code
Transcript_Gene_Up_Genes <- Transcript_Gene_Up$Gene_Stable_ID
Transcript_Gene_Up_Genes <- Transcript_Gene_Up_Genes[!duplicated(Transcript_Gene_Up_Genes)]
length(Transcript_Gene_Up_Genes)
# [1] 3269

# Apply a log2FC threshold
Transcript_Gene_Up_0.2 <- dplyr::filter(Transcript_Gene_Up, log2FC.x > 0.2 & log2FC.y > 0.2)

# View the number of non-duplicated Gene IDs
Transcript_Gene_Up_0.2_Genes <- dplyr::distinct(Transcript_Gene_Up_0.2, Gene_Stable_ID, .keep_all = TRUE)
dim(Transcript_Gene_Up_0.2_Genes)
# [1] 3269 (log2FC > 0)
# [1] 2899 (log2FC > 0.1)
# [1] 1819 (log2FC > 0.2)
# [1] 1061 (log2FC > 0.3)
# [1] 675 (log2FC > 0.4)
# [1] 445 (log2FC > 0.5)

#### Transcript Down Only ####
# Transcript Down Only (What is downregulated at the transcript level, but is 
# not downregulated at the gene level?)

# Here we will use anti_join to return all rows from x (Transcript Down) 
# withOUT a match in y (Gene Down)
Transcript_Down_Only <- anti_join(Transcript_Down, Gene_Down, by = "Gene_Stable_ID")
# A tibble: 3,815 × 10

# One way to double check is to see that this matches the #s that we get with 
# the Venn Diagram code (We expect ### Gene IDs that are only detected at the
# transcript level)
Transcript_Down_Only_Genes <- Transcript_Down_Only$Gene_Stable_ID
Transcript_Down_Only_Genes <- Transcript_Down_Only_Genes[!duplicated(Transcript_Down_Only_Genes)]
length(Transcript_Down_Only_Genes)
# [1] 3117

# Now, here we will apply a log2FC cutoff, only keeping transcripts that
# have a log2FC > 0.2
Transcript_Down_Only_0.2 <- dplyr::filter(Transcript_Down_Only, abs(log2FC) > 0.2)
dim(Transcript_Down_Only_0.2)
# [1] 3306   11

# Now we will remove Gene IDs that are duplicated
Transcript_Down_Only_0.2_Genes <- dplyr::distinct(Transcript_Down_Only_0.2, Gene_Stable_ID, .keep_all = TRUE)
dim(Transcript_Down_Only_0.2_Genes)
# [1] 3117 (log2FC > 0)
# [1] 3001 (log2FC > 0.1)
# [1] 2724 (log2FC > 0.2)
# [1] 2510 (log2FC > 0.3)
# [1] 2303 (log2FC > 0.4)
# [1] 2060 (log2FC > 0.5)

#### Gene Down ####
# Gene Down # (What is upregulated at the gene level, but is 
# not upregulated at the transcript level?)

# Here we will use anti_join to return all rows from x (Gene Down) 
# withOUT a match in y (Transcript Down)
Gene_Down_Only <- anti_join(Gene_Down, Transcript_Down, by = "Gene_Stable_ID")
# A tibble: 528 × 8

# We don't need to remove duplicates from gene level results
Gene_Down_Only_0.2_Genes <- dplyr::filter(Gene_Down_Only, abs(log2FC) > 0.2)
dim(Gene_Down_Only_0.2_Genes)
# [1] 528   8 (log2FC > 0)
# [1] 447   8 (log2FC > 0.1)
# [1] 290   8 (log2FC > 0.2)
# [1] 205   8 (log2FC > 0.3)
# [1] 145   8 (log2FC > 0.4)
# [1] 105   8 (log2FC > 0.5)

##### Transcript and Gene Both Down #####
# Here we identify which genes are detected by both transcript and gene 
# expression analysis

# Inner join will keep all rows in x (Transcript_Down) and y (Gene_Down)
Transcript_Gene_Down <- inner_join(Transcript_Down, Gene_Down, by = "Gene_Stable_ID")
dim(Transcript_Gene_Down)
# [1] 5894   17

# View length after removing duplicates. Why? This is a sanity check to make
# sure what we detect here matches what we see with the Venn Diagram Code
Transcript_Gene_Down_Genes <- Transcript_Gene_Down$Gene_Stable_ID
Transcript_Gene_Down_Genes <- Transcript_Gene_Down_Genes[!duplicated(Transcript_Gene_Down_Genes)]
length(Transcript_Gene_Down_Genes)
# [1] 3836

# Apply a log2FC threshold
Transcript_Gene_Down_0.2 <- dplyr::filter(Transcript_Gene_Down, abs(log2FC.x) > 0.2 & abs(log2FC.y) > 0.2)

# View the number of non-duplicated Gene IDs
Transcript_Gene_Down_0.2_Genes <- dplyr::distinct(Transcript_Gene_Down_0.2, Gene_Stable_ID, .keep_all = TRUE)
dim(Transcript_Gene_Down_0.2_Genes)
# [1] 3836 (log2FC > 0)
# [1] 3469 (log2FC > 0.1)
# [1] 2373 (log2FC > 0.2)
# [1] 1542 (log2FC > 0.3)
# [1] 1027 (log2FC > 0.4)
# [1] 672 (log2FC > 0.5)

#### Now we can determine positive and negative control recovery ####
# Sum the transcript total
Transcript_Total <- c(Transcript_Up_Only_0.2_Genes$Gene_Stable_ID, Transcript_Down_Only_0.2_Genes$Gene_Stable_ID, Transcript_Gene_Up_0.2_Genes$Gene_Stable_ID, Transcript_Gene_Down_0.2_Genes$Gene_Stable_ID)
length(Transcript_Total[!duplicated(Transcript_Total)])
# [1] 10439 for abs(log2FC) > 0
# [1] 9687 for abs(log2FC) > 0.1
# [1] 7341 for abs(log2FC) > 0.2
# [1] 5482 for abs(log2FC) > 0.3
# [1] 4270 for abs(log2FC) > 0.4
# [1] 3367 for abs(log2FC) > 0.5

# Sum the gene total
Gene_Total <- c(Gene_Up_Only_0.2_Genes$Gene_Stable_ID, Gene_Down_Only_0.2_Genes$Gene_Stable_ID, Transcript_Gene_Up_0.2_Genes$Gene_Stable_ID, Transcript_Gene_Down_0.2_Genes$Gene_Stable_ID)
length(Gene_Total[!duplicated(Gene_Total)])
# [1] 8505 for abs(log2FC) > 0
# [1] 7506 for abs(log2FC) > 0.1
# [1] 4863 for abs(log2FC) > 0.2
# [1] 3027 for abs(log2FC) > 0.3
# [1] 1984 for abs(log2FC) > 0.4
# [1] 1308 for abs(log2FC) > 0.5

# Intersect with positive controls
Gene_Positive_Controls <- read.table("Additional_File2_Positive_Controls.txt", header = TRUE)
Intersect_Positive <- intersect(Gene_Positive_Controls[, 1], Gene_Total)
length(Intersect_Positive)
# [1] 558 for abs(log2FC) > 0
# [1] 535 for abs(log2FC) > 0.1
# [1] 411 for abs(log2FC) > 0.2
# [1] 294 for abs(log2FC) > 0.3
# [1] 214 for abs(log2FC) > 0.4
# [1] 157 for abs(log2FC) > 0.5

# Intersect with negative controls
Gene_Negative_Controls <- read.table("Final_Negative_Controls_v25.txt", header = TRUE)
Intersect_Negative <- intersect(Gene_Negative_Controls[, 1], Gene_Total)
length(Intersect_Negative)
# [1] 1309 for abs(log2FC) > 0
# [1] 1105 for abs(log2FC) > 0.1
# [1] 592 for abs(log2FC) > 0.2
# [1] 318 for abs(log2FC) > 0.3
# [1] 188 for abs(log2FC) > 0.4
# [1] 124 for abs(log2FC) > 0.5

#### Now intersect Highlight Gene Lists ####
Highlight_Genes <- readxl::read_excel("Peixoto_Figure_2_Supplementary_Table_2.xlsx", sheet = 1, col_names = TRUE)
Gene_Names <- Highlight_Genes$Gene_Names
Gene_Names <- as.matrix(Gene_Names)
dim(Gene_Names)
# [1] 47  1
Gene_Names <- unlist(strsplit(Gene_Names,","))
Gene_Names <- noquote(Gene_Names)
# Remove trailing spaces and save non-duplicated values
Gene_Names <- trimws(Gene_Names)
length(Gene_Names)
# [1] 47

# Transcript Up Only
Transcript_Up_Only_Gene_Names <- intersect(Transcript_Up_Only_0.2_Genes$Gene_Name, Gene_Names)
Transcript_Up_Only_Eifs <- Transcript_Up_Only_0.2_Genes[grep("^Eif", Transcript_Up_Only_0.2_Genes$Gene_Name),]
Transcript_Up_Only_Ctnn <- Transcript_Up_Only_0.2_Genes[grep("^Ctnna", Transcript_Up_Only_0.2_Genes$Gene_Name),]
Transcript_Up_Only_Cntn <- Transcript_Up_Only_0.2_Genes[grep("^Cntn", Transcript_Up_Only_0.2_Genes$Gene_Name),]
Transcript_Up_Only_Hdac <- Transcript_Up_Only_0.2_Genes[grep("^Hdac", Transcript_Up_Only_0.2_Genes$Gene_Name),]
Transcript_Up_Only_Fos <- Transcript_Up_Only_0.2_Genes[grep("^Fos", Transcript_Up_Only_0.2_Genes$Gene_Name),]
Transcript_Up_Only_Grin <- Transcript_Up_Only_0.2_Genes[grep("^Grin", Transcript_Up_Only_0.2_Genes$Gene_Name),]
Transcript_Up_Only_Gria <- Transcript_Up_Only_0.2_Genes[grep("^Gria", Transcript_Up_Only_0.2_Genes$Gene_Name),]

Transcript_Up_Gene_Names <- c(Transcript_Up_Only_Gene_Names, Transcript_Up_Only_Eifs$Gene_Name, Transcript_Up_Only_Ctnn$Gene_Name, Transcript_Up_Only_Cntn$Gene_Name, Transcript_Up_Only_Hdac$Gene_Name, Transcript_Up_Only_Fos$Gene_Name, Transcript_Up_Only_Grin$Gene_Name, Transcript_Up_Only_Gria$Gene_Name)
sort(Transcript_Up_Gene_Names[!duplicated(Transcript_Up_Gene_Names)])
# [1] "Cdkn1b"   "Cntn4"    "Cntnap5c" "Eif2a"    "Eif2b4"   "Eif2d"    "Eif2s3y"  "Eif5"     "Hdac3"   

# Gene Up Only
Gene_Up_Only_Gene_Names <- intersect(Gene_Up_Only_0.2_Genes$Gene_Name, Gene_Names)
Gene_Up_Only_Eifs <- Gene_Up_Only_0.2_Genes[grep("^Eif", Gene_Up_Only_0.2_Genes$Gene_Name),]
Gene_Up_Only_Ctnn <- Gene_Up_Only_0.2_Genes[grep("^Ctnna", Gene_Up_Only_0.2_Genes$Gene_Name),]
Gene_Up_Only_Cntn <- Gene_Up_Only_0.2_Genes[grep("^Cntn", Gene_Up_Only_0.2_Genes$Gene_Name),]
Gene_Up_Only_Hdac <- Gene_Up_Only_0.2_Genes[grep("^Hdac", Gene_Up_Only_0.2_Genes$Gene_Name),]
Gene_Up_Only_Fos <- Gene_Up_Only_0.2_Genes[grep("^Fos", Gene_Up_Only_0.2_Genes$Gene_Name),]
Gene_Up_Only_Grin <- Gene_Up_Only_0.2_Genes[grep("^Grin", Gene_Up_Only_0.2_Genes$Gene_Name),]
Gene_Up_Only_Gria <- Gene_Up_Only_0.2_Genes[grep("^Gria", Gene_Up_Only_0.2_Genes$Gene_Name),]

Gene_Up_Gene_Names <- c(Gene_Up_Only_Gene_Names, Gene_Up_Only_Eifs$Gene_Name, Gene_Up_Only_Ctnn$Gene_Name, Gene_Up_Only_Cntn$Gene_Name, Gene_Up_Only_Hdac$Gene_Name, Gene_Up_Only_Fos$Gene_Name, Gene_Up_Only_Grin$Gene_Name, Gene_Up_Only_Gria$Gene_Name)
sort(Gene_Up_Gene_Names[!duplicated(Gene_Up_Gene_Names)])
# [1] "Eif5b"

# Transcript and Gene Up 
Transcript_Gene_Up_Gene_Names <- intersect(Transcript_Gene_Up_0.2_Genes$Gene_Name.x, Gene_Names)
Transcript_Gene_Up_Eifs <- Transcript_Gene_Up_0.2_Genes[grep("^Eif", Transcript_Gene_Up_0.2_Genes$Gene_Name.x),]
Transcript_Gene_Up_Ctnn <- Transcript_Gene_Up_0.2_Genes[grep("^Ctnna", Transcript_Gene_Up_0.2_Genes$Gene_Name.x),]
Transcript_Gene_Up_Cntn <- Transcript_Gene_Up_0.2_Genes[grep("^Cntn", Transcript_Gene_Up_0.2_Genes$Gene_Name.x),]
Transcript_Gene_Up_Hdac <- Transcript_Gene_Up_0.2_Genes[grep("^Hdac", Transcript_Gene_Up_0.2_Genes$Gene_Name.x),]
Transcript_Gene_Up_Fos <- Transcript_Gene_Up_0.2_Genes[grep("^Fos", Transcript_Gene_Up_0.2_Genes$Gene_Name.x),]
Transcript_Gene_Up_Grin <- Transcript_Gene_Up_0.2_Genes[grep("^Grin", Transcript_Gene_Up_0.2_Genes$Gene_Name.x),]
Transcript_Gene_Up_Gria <- Transcript_Gene_Up_0.2_Genes[grep("^Gria", Transcript_Gene_Up_0.2_Genes$Gene_Name.x),]

Transcript_Gene_Up_Gene_Names <- c(Transcript_Gene_Up_Gene_Names, Transcript_Gene_Up_Eifs$Gene_Name.x, Transcript_Gene_Up_Ctnn$Gene_Name.x, Transcript_Gene_Up_Cntn$Gene_Name.x, Transcript_Gene_Up_Hdac$Gene_Name.x, Transcript_Gene_Up_Fos$Gene_Name.x, Transcript_Gene_Up_Grin$Gene_Name.x, Transcript_Gene_Up_Gria$Gene_Name.x)
sort(Transcript_Gene_Up_Gene_Names[!duplicated(Transcript_Gene_Up_Gene_Names)])
# [1] "Arc"      "Bdnf"     "Bhlhe40"  "Camk1g"   "Crh"      "Egr1"     "Eif3j2"   "Eif4ebp2" "Fos"     
# [10] "Fosb"     "Fosl2"    "Gadd45a"  "Gadd45b"  "Gria1"    "Grin2a"   "Grin2b"   "Hdac8"    "Homer1"  
# [19] "Hspa5"    "Hspa8"    "Nr4a1"    "P4ha1"    "Per1"     "Per2"     "Ptgs2"    "Sgk1"     "Sik1"    
# [28] "Sult1a1"  "Vip"  

# Transcript Down Only
Transcript_Down_Only_Gene_Names <- intersect(Transcript_Down_Only_0.2_Genes$Gene_Name, Gene_Names)
Transcript_Down_Only_Eifs <- Transcript_Down_Only_0.2_Genes[grep("^Eif", Transcript_Down_Only_0.2_Genes$Gene_Name),]
Transcript_Down_Only_Ctnn <- Transcript_Down_Only_0.2_Genes[grep("^Ctnna", Transcript_Down_Only_0.2_Genes$Gene_Name),]
Transcript_Down_Only_Cntn <- Transcript_Down_Only_0.2_Genes[grep("^Cntn", Transcript_Down_Only_0.2_Genes$Gene_Name),]
Transcript_Down_Only_Hdac <- Transcript_Down_Only_0.2_Genes[grep("^Hdac", Transcript_Down_Only_0.2_Genes$Gene_Name),]
Transcript_Down_Only_Fos <- Transcript_Down_Only_0.2_Genes[grep("^Fos", Transcript_Down_Only_0.2_Genes$Gene_Name),]
Transcript_Down_Only_Grin <- Transcript_Down_Only_0.2_Genes[grep("^Grin", Transcript_Down_Only_0.2_Genes$Gene_Name),]
Transcript_Down_Only_Gria <- Transcript_Down_Only_0.2_Genes[grep("^Gria", Transcript_Down_Only_0.2_Genes$Gene_Name),]

Transcript_Down_Gene_Names <- c(Transcript_Down_Only_Gene_Names, Transcript_Down_Only_Eifs$Gene_Name, Transcript_Down_Only_Ctnn$Gene_Name, Transcript_Down_Only_Cntn$Gene_Name, Transcript_Down_Only_Hdac$Gene_Name, Transcript_Down_Only_Fos$Gene_Name, Transcript_Down_Only_Grin$Gene_Name, Transcript_Down_Only_Gria$Gene_Name)
sort(Transcript_Down_Gene_Names[!duplicated(Transcript_Down_Gene_Names)])
# [1] "Bdnf"      "Cntn1"     "Cntn4"     "Cntn5"     "Eif2ak3"   "Eif2b1"    "Eif2s2"    "Eif2s3y"  
# [9] "Eif3a"     "Eif3i"     "Eif3m"     "Eif4a1"    "Eif4a2"    "Eif4e2"    "Eif4enif1" "Eif4g1"   
# [17] "Eif4h"     "Nr4a1"  

# Gene Down Only
Gene_Down_Only_Gene_Names <- intersect(Gene_Down_Only_0.2_Genes$Gene_Name, Gene_Names)
Gene_Down_Only_Eifs <- Gene_Down_Only_0.2_Genes[grep("^Eif", Gene_Down_Only_0.2_Genes$Gene_Name),]
Gene_Down_Only_Ctnn <- Gene_Down_Only_0.2_Genes[grep("^Ctnna", Gene_Down_Only_0.2_Genes$Gene_Name),]
Gene_Down_Only_Cntn <- Gene_Down_Only_0.2_Genes[grep("^Cntn", Gene_Down_Only_0.2_Genes$Gene_Name),]
Gene_Down_Only_Hdac <- Gene_Down_Only_0.2_Genes[grep("^Hdac", Gene_Down_Only_0.2_Genes$Gene_Name),]
Gene_Down_Only_Fos <- Gene_Down_Only_0.2_Genes[grep("^Fos", Gene_Down_Only_0.2_Genes$Gene_Name),]
Gene_Down_Only_Grin <- Gene_Down_Only_0.2_Genes[grep("^Grin", Gene_Down_Only_0.2_Genes$Gene_Name),]
Gene_Down_Only_Gria <- Gene_Down_Only_0.2_Genes[grep("^Hdac", Gene_Down_Only_0.2_Genes$Gene_Name),]

Gene_Down_Gene_Names <- c(Gene_Down_Only_Gene_Names, Gene_Down_Only_Eifs$Gene_Name, Gene_Down_Only_Ctnn$Gene_Name, Gene_Down_Only_Cntn$Gene_Name, Gene_Down_Only_Hdac$Gene_Name, Gene_Down_Only_Fos$Gene_Name, Gene_Down_Only_Grin$Gene_Name, Gene_Down_Only_Gria$Gene_Name)
sort(Gene_Down_Gene_Names[!duplicated(Gene_Down_Gene_Names)])
# character(0)

# Transcript and Gene Down 
Transcript_Gene_Down_Gene_Names <- intersect(Transcript_Gene_Down_0.2_Genes$Gene_Name.x, Gene_Names)
Transcript_Gene_Down_Eifs <- Transcript_Gene_Down_0.2_Genes[grep("^Eif", Transcript_Gene_Down_0.2_Genes$Gene_Name.x),]
Transcript_Gene_Down_Ctnn <- Transcript_Gene_Down_0.2_Genes[grep("^Ctnna", Transcript_Gene_Down_0.2_Genes$Gene_Name.x),]
Transcript_Gene_Down_Cntn <- Transcript_Gene_Down_0.2_Genes[grep("^Cntn", Transcript_Gene_Down_0.2_Genes$Gene_Name.x),]
Transcript_Gene_Down_Hdac <- Transcript_Gene_Down_0.2_Genes[grep("^Hdac", Transcript_Gene_Down_0.2_Genes$Gene_Name.x),]
Transcript_Gene_Down_Fos <- Transcript_Gene_Down_0.2_Genes[grep("^Fos", Transcript_Gene_Down_0.2_Genes$Gene_Name.x),]
Transcript_Gene_Down_Grin <- Transcript_Gene_Down_0.2_Genes[grep("^Grin", Transcript_Gene_Down_0.2_Genes$Gene_Name.x),]
Transcript_Gene_Down_Gria <- Transcript_Gene_Down_0.2_Genes[grep("^Gria", Transcript_Gene_Down_0.2_Genes$Gene_Name.x),]

Transcript_Gene_Down_Gene_Names <- c(Transcript_Gene_Down_Gene_Names, Transcript_Gene_Down_Eifs$Gene_Name.x, Transcript_Gene_Down_Ctnn$Gene_Name.x, Transcript_Gene_Down_Cntn$Gene_Name.x, Transcript_Gene_Down_Hdac$Gene_Name.x, Transcript_Gene_Down_Fos$Gene_Name.x, Transcript_Gene_Down_Grin$Gene_Name.x, Transcript_Gene_Down_Gria$Gene_Name.x)
sort(Transcript_Gene_Down_Gene_Names[!duplicated(Transcript_Gene_Down_Gene_Names)])
# [1] "Cirbp"    "Cntn6"    "Dact2"    "Dbp"      "Eif2ak2"  "Eif4ebp1" "Hdac7"    "Hdac9"    "Mef2c"   
# [10] "Nfatc3"   "Npas1"    "Opalin"   "Tfrc"     "Tipin"    "Usp43"    "Wnt9a"   

