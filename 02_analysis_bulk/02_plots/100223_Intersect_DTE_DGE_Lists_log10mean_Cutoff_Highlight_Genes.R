#### Applying a log10mean filter to DTE and DGE analysis ####

# Summary: I will look at the number of significant genes and transcripts if we 
# apply a log10mean filter, as well as the positive and negative controls.

# See organized results in Extended Data Figure 10

# R Version 4.2.2

#### LOAD PACKAGES ####
library("VennDiagram") # 1.7.3
library("readxl") # 1.4.1

#### READ LISTS ####
# Read lists without a fold change cutoff
# Gene
Gene <- readxl::read_excel("060123_Peixoto_Figure1_Supplement2.xlsx", 
                           sheet = "DGE")
dim(Gene)
# [1] 8505    8

# Transcript Expression
Transcript <- readxl::read_excel("060123_Peixoto_Figure1_Supplement2.xlsx", 
                           sheet = "DTE")
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
Transcript_Up_Only_ID <- Transcript_Up_Only$Gene_Stable_ID
length(Transcript_Up_Only_ID[!duplicated(Transcript_Up_Only_ID)])
# [1] 1381

# Now, here we will apply a log2FC cutoff, only keeping transcripts that
# have a log2FC > 0.2
Transcript_Up_Only_1 <- dplyr::filter(Transcript_Up_Only, log10mean > 1)
dim(Transcript_Up_Only_1)
# [1] 1536   11 (log10mean > 0)
# [1] 1535   11 (log10mean > 0.5)
# [1] 1523   11 (log10mean > 1)
# [1] 1397   11 (log10mean > 1.5)
# [1] 1165   11 (log10mean > 2)
# [1] 868  11 (log10mean > 2.5)
# [1] 487  11 (log10mean > 3)
# [1] 145  11 (log10mean > 3.5)
# [1] 25 11 (log10mean > 4)

# Now we will remove Gene IDs that are duplicated
Transcript_Up_Only_ID <- Transcript_Up_Only_1$Gene_Stable_ID
length(Transcript_Up_Only_ID[!duplicated(Transcript_Up_Only_ID)])

# [1] 1381 (log10mean > 0)
# [1] 1380 (log10mean > 0.5)
# [1] 1374 (log10mean > 1)
# [1] 1271 (log10mean > 1.5)
# [1] 1083 (log10mean > 2)
# [1] 821 (log10mean > 2.5)
# [1] 475 (log10mean > 3)
# [1] 143 (log10mean > 3.5)
# [1] 24 (log10mean > 4)

#### Gene Up Only ####
# Gene Up # (What is upregulated at the gene level, but is 
# not upregulated at the transcript level?)

# Here we will use anti_join to return all rows from x (Gene Up) 
# withOUT a match in y (Transcript Up)

Gene_Up_Only <- anti_join(Gene_Up, Transcript_Up, by = "Gene_Stable_ID")
# A tibble: 872 × 8

Gene_Up_Only_1 <- dplyr::filter(Gene_Up_Only, log10mean > 1)
dim(Gene_Up_Only_1)
# [1] 872   8 (log10mean > 0)
# [1] 872   8 (log10mean > 0.5)
# [1] 869   8 (log10mean > 1)
# [1] 798   8 (log10mean > 1.5)
# [1] 703   8 (log10mean > 2)
# [1] 580   8 (log10mean > 2.5)
# [1] 383   8 (log10mean > 3)
# [1] 121   8 (log10mean > 3.5)
# [1] 18  8 (log10mean > 4)

#### Transcript and Gene Both Up ####
# Here we identify which genes are detected by both transcript and gene 
# expression analysis

# Inner join will keep all rows in x (Transcript_Up) and y (Gene_Up)
Transcript_Gene_Up <- inner_join(Transcript_Up, Gene_Up, by = "Gene_Stable_ID")
dim(Transcript_Gene_Up)
# [1] 4280   17

# View length after removing duplicates. Why? This is a sanity check to make
# sure what we detect here matches what we see with the Venn Diagram Code
length(Transcript_Gene_Up$Gene_Stable_ID[!duplicated(Transcript_Gene_Up$Gene_Stable_ID)])
# [1] 3269

# Apply a log2FC threshold
Transcript_Gene_Up_1 <- dplyr::filter(Transcript_Gene_Up, log10mean.x > 1 & log10mean.y > 1)

# View the number of non-duplicated Gene IDs
Transcript_Gene_Up_ID <- Transcript_Gene_Up_1$Gene_Stable_ID
length(Transcript_Gene_Up_ID[!duplicated(Transcript_Gene_Up_ID)])
# [1] 3269 (log10mean > 0)
# [1] 3269 (log10mean > 0.5)
# [1] 3266 (log10mean > 1)
# [1] 3128 (log10mean > 1.5)
# [1] 2877 (log10mean > 2)
# [1] 2433 (log10mean > 2.5)
# [1] 1649 (log10mean > 3)
# [1] 630 (log10mean > 3.5)
# [1] 130 (log10mean > 4)

#### Transcript Down Only ####
# Transcript Down Only (What is downregulated at the transcript level, but is 
# not downregulated at the gene level?)

# Here we will use anti_join to return all rows from x (Transcript Down) 
# withOUT a match in y (Gene Down)
Transcript_Down_Only <- anti_join(Transcript_Down, Gene_Down, by = "Gene_Stable_ID")
# A tibble: 3,815 × 10
Transcript_Down_Only_ID <- Transcript_Down_Only$Gene_Stable_ID
length(Transcript_Down_Only_ID[!duplicated(Transcript_Down_Only_ID)])
# [1] 3117

# Apply a log10mean cutoff
Transcript_Down_Only_1 <- dplyr::filter(Transcript_Down_Only, log10mean > 1)

Transcript_Down_Only_ID <- Transcript_Down_Only_1$Gene_Stable_ID
length(Transcript_Down_Only_ID[!duplicated(Transcript_Down_Only_ID)])
# [1] 3117 (log10mean > 0)
# [1] 3117 (log10mean > 0.5)
# [1] 2773 (log10mean > 1)
# [1] 2082 (log10mean > 1.5)
# [1] 1462 (log10mean > 2)
# [1] 883 (log10mean > 2.5)
# [1] 428 (log10mean > 3)
# [1] 157 (log10mean > 3.5)
# [1] 35 (log10mean > 4)

#### Gene Down ####
# Gene Down # (What is upregulated at the gene level, but is 
# not upregulated at the transcript level?)

# Here we will use anti_join to return all rows from x (Gene Down) 
# withOUT a match in y (Transcript Down)
Gene_Down_Only <- anti_join(Gene_Down, Transcript_Down, by = "Gene_Stable_ID")
# A tibble: 528 × 8

Gene_Down_Only_1 <- dplyr::filter(Gene_Down_Only, log10mean > 1)
dim(Gene_Down_Only_1)
# [1] 528   8 (log10mean > 0)
# [1] 528   8 (log10mean > 0.5)
# [1] 504   8 (log10mean > 1)
# [1] 415   8 (log10mean > 1.5)
# [1] 320   8 (log10mean > 2)
# [1] 232   8 (log10mean > 2.5)
# [1] 132   8 (log10mean > 3)
# [1] 49  8 (log10mean > 3.5)
# [1] 10  8 (log10mean > 4)

##### Transcript and Gene Both Down #####
# Here we identify which genes are detected by both transcript and gene 
# expression analysis

# Inner join will keep all rows in x (Transcript_Down) and y (Gene_Down)
Transcript_Gene_Down <- inner_join(Transcript_Down, Gene_Down, by = "Gene_Stable_ID")
dim(Transcript_Gene_Down)
# [1] 5894   17

# View length after removing duplicates. Why? This is a sanity check to make
# sure what we detect here matches what we see with the Venn Diagram Code
length(Transcript_Gene_Down$Gene_Stable_ID[!duplicated(Transcript_Gene_Down$Gene_Stable_ID)])
# [1] 3836

Transcript_Gene_Down_1 <- dplyr::filter(Transcript_Gene_Down, log10mean.x > 1 & log10mean.y > 1)

Transcript_Gene_Down_ID <- Transcript_Gene_Down_1$Gene_Stable_ID
length(Transcript_Gene_Down_ID[!duplicated(Transcript_Gene_Down_ID)])
# [1] 3836 (log10mean > 0)
# [1] 3836 (log10mean > 0.5)
# [1] 3701 (log10mean > 1)
# [1] 3239 (log10mean > 1.5)
# [1] 2640 (log10mean > 2)
# [1] 1882 (log10mean > 2.5)
# [1] 1033 (log10mean > 3)
# [1] 379 (log10mean > 3.5)
# [1] 97 (log10mean > 4)

#### Now we can determine positive and negative control recovery ####
# Sum the transcript total
Transcript_Total <- rbind(Transcript_Up_Only_1[,1],Transcript_Down_Only_1[,1], Transcript_Gene_Up_1[,1], Transcript_Gene_Down_1[,1])
length(Transcript_Total$Gene_Stable_ID[!duplicated(Transcript_Total$Gene_Stable_ID)])
# [1] 10439 (log10mean > 0)
# [1] 10439 (log10mean > 0.5)
# [1] 10008 (log10mean > 1)
# [1] 8859 (log10mean > 1.5)
# [1] 7470 (log10mean > 2)
# [1] 5699 (log10mean > 2.5)
# [1] 3485 (log10mean > 3)
# [1] 1287 (log10mean > 3.5)
# [1] 285 (log10mean > 4)

# Sum the gene total
Gene_Total <- rbind(Gene_Up_Only_1[,1], Gene_Down_Only_1[,1], Transcript_Gene_Up_1[,1], Transcript_Gene_Down_1[,1])
length(Gene_Total$Gene_Stable_ID[!duplicated(Gene_Total$Gene_Stable_ID)])
# [1] 8505 (log10mean > 0)
# [1] 8505 (log10mean > 0.5)
# [1] 8340 (log10mean > 1)
# [1] 7580 (log10mean > 1.5)
# [1] 6540 (log10mean > 2)
# [1] 5127 (log10mean > 2.5)
# [1] 3197 (log10mean > 3)
# [1] 1179 (log10mean > 3.5)
# [1] 255 (log10mean > 4)

# Intersect with positive controls
Gene_Positive_Controls <- read.table("Additional_File2_Positive_Controls.txt", header = TRUE)
Intersect_Positive <- intersect(Gene_Positive_Controls[, 1], Gene_Total$Gene_Stable_ID)
length(Intersect_Positive)
# [1] 558 (log10mean > 0.5)
# [1] 556 (log10mean > 1)
# [1] 541 (log10mean > 1.5)
# [1] 499 (log10mean > 2)
# [1] 396 (log10mean > 2.5)
# [1] 227 (log10mean > 3)
# [1] 67 (log10mean > 3.5)
# [1] 10 (log10mean > 4)

# Intersect with negative controls
Gene_Negative_Controls <- read.table("Final_Negative_Controls_v25.txt", header = TRUE)
Intersect_Negative <- intersect(Gene_Negative_Controls[, 1], Gene_Total$Gene_Stable_ID)
length(Intersect_Negative)
# [1] 1309 (log10mean > 0.5)
# [1] 1289 (log10mean > 1)
# [1] 1200 (log10mean > 1.5)
# [1] 1028 (log10mean > 2)
# [1] 811 (log10mean > 2.5)
# [1] 498 (log10mean > 3)
# [1] 178 (log10mean > 3.5)
# [1] 46 (log10mean > 4)

