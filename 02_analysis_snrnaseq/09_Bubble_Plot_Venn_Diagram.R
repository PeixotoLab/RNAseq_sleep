#### Bubble Plot Code for Glutamatergic/GABAergic Venn Diagram Functional Annotation ####

##### Summary ####

### For the Venn Diagram, we have too many terms to fit the space, so we will need 
# to reduce clusters to one point vs n points.

# Main points for how I will do this:
### I will plot the geometric mean of the fold enrichments within a cluster 
# on the x-axis. I chose the geometric mean as this is less likely to be 
# affected by extreme values, and it is implemented by DAVID for the 
# calculation of the Enrichment Score (the only difference there is that the 
# pvalues are used as input)

### Each cluster will then be on the y-axis, then unclustered terms will be listed.
# For unclustered terms, I will plot the fold enrichment on the x-axis.

### The size of the bubbles will be determined by the # of unique genes within the term
# or cluster

### The color for now will be one shade of red for upregulated and one shade of 
# blue for downregulated

##### Load packages and set the working directory ####

# Load packages
library(ggplot2) # Version 3.4.2
library(RColorBrewer) # Version 1.1-3
library(grid) # Version 4.2.2
library(dplyr) # Version 1.1.1
library(psych) # Version 2.3.6

# Set the working directory
setwd("~/Dropbox/Sleep_RNAseq_Splicing/Single_nuclear_RNAseq_SD/KF_Venn_Diagram/DAVID_Output")

#### Read the text files that were exported from DAVID ####
# As the text file is tab delimited, we define the field separator as "\t"
# "fill = TRUE" adds blank spaces to uneven rows

##### Shared between Glutamatergic and GABAergic #####
Glutamatergic_GABAergic_Up_Unclustered <- read.table("100523_DAVID_Glutamatergic_GABAergic_Shared_Upregulated_Unclustered.txt", sep = "\t", fill = TRUE)
Glutamatergic_GABAergic_Down_Unclustered <- read.table("100523_DAVID_Glutamatergic_GABAergic_Shared_Downregulated_Unclustered.txt", sep = "\t", fill = TRUE)

##### Glutamatergic Only #####
Glutamatergic_Only_Up_Clustered <- read.table("100523_DAVID_Glutamatergic_Only_Upregulated_Clustered.txt", sep = "\t", fill = TRUE)
Glutamatergic_Only_Up_Unclustered <- read.table("100523_DAVID_Glutamatergic_Only_Upregulated_Unclustered.txt", sep = "\t", fill = TRUE)
Glutamatergic_Only_Down_Clustered <- read.table("100523_DAVID_Glutamatergic_Only_Downregulated_Clustered.txt", sep = "\t", fill = TRUE)

##### GABAergic Only #####
GABAergic_Only_Up_Unclustered <- read.table("100523_DAVID_GABAergic_Only_Upregulated_Unclustered.txt", sep = "\t", fill = TRUE)
GABAergic_Only_Down_Unclustered <- read.table("100523_DAVID_GABAergic_Only_Downregulated_Unclustered.txt", sep = "\t", fill = TRUE)

#### Record dimensions of matrix ####
# After loading the text files into R, you can view the entire object by calling the name of the object (in this case "L23_Unclustered_Up", etc.)
# You can view the first few rows with head(object name)
# You can see the dimensions of the object with dim(object name): the first number is the number of rows and the second number is the number of columns.

##### Shared between Glutamatergic and GABAergic #####
# View the dimensions
dim(Glutamatergic_GABAergic_Up_Unclustered)
# [1]  6 13

# Subset to remove headers, and keep all terms
Glutamatergic_GABAergic_Up_Unclustered <- Glutamatergic_GABAergic_Up_Unclustered[2:6,]

# Ensure the dimensions following subsetting are correct
dim(Glutamatergic_GABAergic_Up_Unclustered)
# [1]  5 13

# Repeat with the downregulated file 
# View the dimensions
dim(Glutamatergic_GABAergic_Down_Unclustered)
# [1]  5 13

# Subset to remove headers, and keep all terms
Glutamatergic_GABAergic_Down_Unclustered <- Glutamatergic_GABAergic_Down_Unclustered[2:5,]
# Ensure the dimensions following subsetting are correct
dim(Glutamatergic_GABAergic_Down_Unclustered)
# [1]  4 13

##### Glutamatergic Only #####
# Clustered
# View the dimensions
dim(Glutamatergic_Only_Up_Clustered)
# [1] 52 13

# Subset to remove headers, and keep all terms
# Note that clustered terms have an extra header
Glutamatergic_Only_Up_Clustered <- Glutamatergic_Only_Up_Clustered[c(3:6, 9:16, 19:52),]

# Ensure the dimensions following subsetting are correct
dim(Glutamatergic_Only_Up_Clustered)
# [1] 46 13

# Unclustered
# View the dimensions
dim(Glutamatergic_Only_Up_Unclustered)
# [1] 14 13

# Subset to remove headers, and keep all terms
Glutamatergic_Only_Up_Unclustered <- Glutamatergic_Only_Up_Unclustered[c(2:14),]

# Ensure the dimensions following subsetting are correct
dim(Glutamatergic_Only_Up_Unclustered)
# [1] 13 13


# Repeat with the downregulated file 
# View the dimensions
dim(Glutamatergic_Only_Down_Clustered)
# [1] 12 13

# Subset to remove headers, and keep all terms
Glutamatergic_Only_Down_Clustered <- Glutamatergic_Only_Down_Clustered[c(3:6, 9:12),]
# Ensure the dimensions following subsetting are correct
dim(Glutamatergic_Only_Down_Clustered)
# [1]  8 13

##### GABAergic Only #####
# View the dimensions
dim(GABAergic_Only_Up_Unclustered)
# [1]  2 13

# Subset to remove headers, and keep all terms
GABAergic_Only_Up_Unclustered <- GABAergic_Only_Up_Unclustered[2,]

# Ensure the dimensions following subsetting are correct
dim(GABAergic_Only_Up_Unclustered)
# [1]  1 13

# Repeat with the downregulated file 
# View the dimensions
dim(GABAergic_Only_Down_Unclustered)
# [1]  3 13

# Subset to remove headers, and keep all terms
GABAergic_Only_Down_Unclustered <- GABAergic_Only_Down_Unclustered[2:3,]
# Ensure the dimensions following subsetting are correct
dim(GABAergic_Only_Down_Unclustered)
# [1]  2 13

#### Change the colnames ####

##### Shared between Glutamatergic and GABAergic #####
colnames(Glutamatergic_GABAergic_Up_Unclustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")
colnames(Glutamatergic_GABAergic_Down_Unclustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")

##### Glutamatergic Only #####
colnames(Glutamatergic_Only_Up_Clustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")
colnames(Glutamatergic_Only_Up_Unclustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")
colnames(Glutamatergic_Only_Down_Clustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")

##### GABAergic Only #####
colnames(GABAergic_Only_Up_Unclustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")
colnames(GABAergic_Only_Down_Unclustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")

#### Join both Clustered and Non-Clustered Terms Using rbind ####
# Note that both data frames you want to combine need to have the same column names.
# In this case, only the upregulated Glutamatergic genes have both clustered and unclustered terms.
Glutamatergic_Only_Up <- rbind(Glutamatergic_Only_Up_Clustered, Glutamatergic_Only_Up_Unclustered)
dim(Glutamatergic_Only_Up)
# [1] 59 13

#### Add a column (Order) to preserve the order of the terms with ggplot2 ####

##### Shared between Glutamatergic and GABAergic #####
Glutamatergic_GABAergic_Order_Up <- c(1:5)
Glutamatergic_GABAergic_Up_Unclustered$Order <- Glutamatergic_GABAergic_Order_Up
dim(Glutamatergic_GABAergic_Up_Unclustered)
# [1]  5 14

# Add "Cluster #" to the rows and/or "Unclustered"
Glutamatergic_GABAergic_Data_Up <- rep(c("Unclustered Up"), times = 5)
Glutamatergic_GABAergic_Groups_Up <- matrix(Glutamatergic_GABAergic_Data_Up, ncol = 1, byrow = TRUE)

Glutamatergic_GABAergic_Up_Unclustered$Groups <- Glutamatergic_GABAergic_Groups_Up

# Down
Glutamatergic_GABAergic_Order_Down <- c(1:4)
Glutamatergic_GABAergic_Down_Unclustered$Order <- Glutamatergic_GABAergic_Order_Down
dim(Glutamatergic_GABAergic_Down_Unclustered)
# [1]  4 14

# Add "Cluster #" to the rows and/or "Unclustered"
Glutamatergic_GABAergic_Data_Down <- rep(c("Unclustered Down"), times = 4)
Glutamatergic_GABAergic_Groups_Down <- matrix(Glutamatergic_GABAergic_Data_Down, ncol = 1, byrow = TRUE)

Glutamatergic_GABAergic_Down_Unclustered$Groups <- Glutamatergic_GABAergic_Groups_Down

##### Glutamatergic Only #####
Glutamatergic_Order_Up <- c(1:59)
Glutamatergic_Only_Up$Order <- Glutamatergic_Order_Up
dim(Glutamatergic_Only_Up)
# [1] 59 14

# Add "Cluster #" to the rows and/or "Unclustered"
Glutamatergic_Data_Up <- rep(c("Cluster 1 Up", "Cluster 2 Up", "Cluster 3 Up", "Unclustered Up"), times = c(4, 8, 34, 13))
Glutamatergic_Groups_Up <- matrix(Glutamatergic_Data_Up, ncol = 1, byrow = TRUE)

Glutamatergic_Only_Up$Groups <- Glutamatergic_Groups_Up

# Down
Glutamatergic_Order_Down <- c(1:8)
Glutamatergic_Only_Down_Clustered$Order <- Glutamatergic_Order_Down
dim(Glutamatergic_Only_Down_Clustered)
# [1]  8 14

# Add "Cluster #" to the rows and/or "Unclustered"
Glutamatergic_Data_Down <- rep(c("Cluster 1 Down", "Cluster 2 Down"), times = c(4, 4))
Glutamatergic_Groups_Down <- matrix(Glutamatergic_Data_Down, ncol = 1, byrow = TRUE)

Glutamatergic_Only_Down_Clustered$Groups <- Glutamatergic_Groups_Down

##### GABAergic Only #####
GABAergic_Order_Up <- 1
GABAergic_Only_Up_Unclustered$Order <- GABAergic_Order_Up
dim(GABAergic_Only_Up_Unclustered)
# [1]  1 14

# Add "Cluster #" to the rows and/or "Unclustered"
GABAergic_Data_Up <- rep(c("Unclustered Up"), times = 1)
GABAergic_Groups_Up <- matrix(GABAergic_Data_Up, ncol = 1, byrow = TRUE)

GABAergic_Only_Up_Unclustered$Groups <- GABAergic_Groups_Up

# Down
GABAergic_Order_Down <- c(1:2)
GABAergic_Only_Down_Unclustered$Order <- GABAergic_Order_Down
dim(GABAergic_Only_Down_Unclustered)
# [1]  2 14

# Add "Cluster #" to the rows and/or "Unclustered"
GABAergic_Data_Down <- rep(c("Unclustered Down"), times = 2)
GABAergic_Groups_Down <- matrix(GABAergic_Data_Down, ncol = 1, byrow = TRUE)

GABAergic_Only_Down_Unclustered$Groups <- GABAergic_Groups_Down

#### Determine the number of unique genes within clusters ####
# We are going to reduce the clusters to a single bubble for plotting.
# To do that, the size will be representative to the number of unique genes within that cluster.
# The fold enrichment (the x-axis in plotting) will be the geometric mean of the enrichments.

##### Glutamatergic Only #####
###### Up ######
# Subset clusters here
Glutamatergic_Cluster1_Up <- dplyr::filter(Glutamatergic_Only_Up, Groups == "Cluster 1 Up")
Glutamatergic_Cluster2_Up <- dplyr::filter(Glutamatergic_Only_Up, Groups == "Cluster 2 Up")
Glutamatergic_Cluster3_Up <- dplyr::filter(Glutamatergic_Only_Up, Groups == "Cluster 3 Up")

# Select only the genes from the data frame
Glutamatergic_Cluster1_Up_Genes <- Glutamatergic_Cluster1_Up$Genes
# Use the unlist function to remove commas between the gene IDs
Glutamatergic_Cluster1_Up_Genes <- unlist(strsplit(Glutamatergic_Cluster1_Up_Genes,","))
# Remove quotes using "noquote"
Glutamatergic_Cluster1_Up_Genes <- noquote(Glutamatergic_Cluster1_Up_Genes)
# Remove trailing spaces and save non-duplicated values
Glutamatergic_Cluster1_Up_Genes <- trimws(Glutamatergic_Cluster1_Up_Genes)
# Determine the number of non-duplicated genes
Glutamatergic_Cluster1_Up_Genes <- Glutamatergic_Cluster1_Up_Genes[!duplicated(Glutamatergic_Cluster1_Up_Genes)]
length(Glutamatergic_Cluster1_Up_Genes)
# [1] 83

# Select only the genes from the data frame
Glutamatergic_Cluster2_Up_Genes <- Glutamatergic_Cluster2_Up$Genes
# Use the unlist function to remove commas between the gene IDs
Glutamatergic_Cluster2_Up_Genes <- unlist(strsplit(Glutamatergic_Cluster2_Up_Genes,","))
# Remove quotes using "noquote"
Glutamatergic_Cluster2_Up_Genes <- noquote(Glutamatergic_Cluster2_Up_Genes)
# Remove trailing spaces and save non-duplicated values
Glutamatergic_Cluster2_Up_Genes <- trimws(Glutamatergic_Cluster2_Up_Genes)
# Determine the number of non-duplicated genes
Glutamatergic_Cluster2_Up_Genes <- Glutamatergic_Cluster2_Up_Genes[!duplicated(Glutamatergic_Cluster2_Up_Genes)]
length(Glutamatergic_Cluster2_Up_Genes)
# [1] 99

# Select only the genes from the data frame
Glutamatergic_Cluster3_Up_Genes <- Glutamatergic_Cluster3_Up$Genes
# Use the unlist function to remove commas between the gene IDs
Glutamatergic_Cluster3_Up_Genes <- unlist(strsplit(Glutamatergic_Cluster3_Up_Genes,","))
# Remove quotes using "noquote"
Glutamatergic_Cluster3_Up_Genes <- noquote(Glutamatergic_Cluster3_Up_Genes)
# Remove trailing spaces and save non-duplicated values
Glutamatergic_Cluster3_Up_Genes <- trimws(Glutamatergic_Cluster3_Up_Genes)
# Determine the number of non-duplicated genes
Glutamatergic_Cluster3_Up_Genes <- Glutamatergic_Cluster3_Up_Genes[!duplicated(Glutamatergic_Cluster3_Up_Genes)]
length(Glutamatergic_Cluster3_Up_Genes)
# [1] 177

###### Down ######
# Repeat the steps for upregulated clusters, but for the downregulated ones.
# Subset clusters here
Glutamatergic_Cluster1_Down <- dplyr::filter(Glutamatergic_Only_Down_Clustered, Groups == "Cluster 1 Down")
Glutamatergic_Cluster2_Down <- dplyr::filter(Glutamatergic_Only_Down_Clustered, Groups == "Cluster 2 Down")

# Select only the genes from the data frame
Glutamatergic_Cluster1_Down_Genes <- Glutamatergic_Cluster1_Down$Genes
# Use the unlist function to remove commas between the gene IDs
Glutamatergic_Cluster1_Down_Genes <- unlist(strsplit(Glutamatergic_Cluster1_Down_Genes,","))
# Remove quotes using "noquote"
Glutamatergic_Cluster1_Down_Genes <- noquote(Glutamatergic_Cluster1_Down_Genes)
# Remove trailing spaces and save non-duplicated values
Glutamatergic_Cluster1_Down_Genes <- trimws(Glutamatergic_Cluster1_Down_Genes)
# Determine the number of non-duplicated genes
Glutamatergic_Cluster1_Down_Genes <- Glutamatergic_Cluster1_Down_Genes[!duplicated(Glutamatergic_Cluster1_Down_Genes)]
length(Glutamatergic_Cluster1_Down_Genes)
# [1] 67

# Select only the genes from the data frame
Glutamatergic_Cluster2_Down_Genes <- noquote(Glutamatergic_Cluster2_Down$Genes)
# Use the unlist function to remove commas between the gene IDs
Glutamatergic_Cluster2_Down_Genes <- unlist(strsplit(Glutamatergic_Cluster2_Down_Genes,","))
# Remove quotes using "noquote"
Glutamatergic_Cluster2_Down_Genes <- noquote(Glutamatergic_Cluster2_Down_Genes)
# Remove trailing spaces and save non-duplicated values
Glutamatergic_Cluster2_Down_Genes <- trimws(Glutamatergic_Cluster2_Down_Genes)
# Determine the number of non-duplicated genes
Glutamatergic_Cluster2_Down_Genes <- Glutamatergic_Cluster2_Down_Genes[!duplicated(Glutamatergic_Cluster2_Down_Genes)]
length(Glutamatergic_Cluster2_Down_Genes)
# [1] 57

#### Determine the geometric mean of the fold enrichments for clustered terms ####
## This is done via the psych package that was loaded above

# First, change all of the "Fold_Enrichment" values to numerics
## Upregulated
Glutamatergic_Cluster1_Up_Fold_Enrichment <- as.numeric(Glutamatergic_Cluster1_Up$Fold_Enrichment)
# [1] 2.221583 1.649979 1.406473 1.406473

Glutamatergic_Cluster2_Up_Fold_Enrichment <- as.numeric(Glutamatergic_Cluster2_Up$Fold_Enrichment)
# [1] 2.151967 1.959432 2.027527 2.097655 2.473709 1.732845 2.312173 2.238604

Glutamatergic_Cluster3_Up_Fold_Enrichment <- as.numeric(Glutamatergic_Cluster3_Up$Fold_Enrichment)
# [1] 2.889292 2.151967 2.861580 2.027527 2.600623 2.097655 2.658865 1.508772 2.093370 2.571248 1.567201
# [12] 2.180307 2.052054 2.525605 2.397906 2.402404 1.824048 2.147944 3.482273 2.128056 2.440280 2.016053
# [23] 2.414181 2.164804 1.938948 1.748496 2.175672 1.994520 1.733412 1.989870 2.432064 2.104671 1.671134
# [34] 2.078029

## Downregulated
Glutamatergic_Cluster1_Down_Fold_Enrichment <- as.numeric(Glutamatergic_Cluster1_Down$Fold_Enrichment)
# [1] 1.704514 2.649345 2.086432 1.854400

Glutamatergic_Cluster2_Down_Fold_Enrichment <- as.numeric(Glutamatergic_Cluster2_Down$Fold_Enrichment)
# [1] 1.628191 2.009956 1.347810 1.347810

## Upregulated
Glutamatergic_Cluster1_Up_GM <- geometric.mean(Glutamatergic_Cluster1_Up_Fold_Enrichment)
# [1] 1.640971

Glutamatergic_Cluster2_Up_GM <- geometric.mean(Glutamatergic_Cluster2_Up_Fold_Enrichment)
# [1] 2.113351

Glutamatergic_Cluster3_Up_GM <- geometric.mean(Glutamatergic_Cluster3_Up_Fold_Enrichment)
# [1] 2.173767

## Downregulated
Glutamatergic_Cluster1_Down_GM <- geometric.mean(Glutamatergic_Cluster1_Down_Fold_Enrichment)
# [1] 2.044498

Glutamatergic_Cluster2_Down_GM <- geometric.mean(Glutamatergic_Cluster2_Down_Fold_Enrichment)
# [1] 1.561483

#### Make matrices for plotting ####
## Upregulated
# The order of the columns of the matrix will be "Order", "Count" then "Enrichment"
# The rows of the matrix will be in order of the clusters
Glutamatergic_Up_Clustered_Matrix <- matrix(data = c(1, 2, 3, length(Glutamatergic_Cluster1_Up_Genes), 
                                       length(Glutamatergic_Cluster2_Up_Genes), length(Glutamatergic_Cluster3_Up_Genes), 
                                       Glutamatergic_Cluster1_Up_GM, Glutamatergic_Cluster2_Up_GM, Glutamatergic_Cluster3_Up_GM), 
                                       nrow = 3, byrow = FALSE)

# Set the column names using "colnames"
colnames(Glutamatergic_Up_Clustered_Matrix) <- c("Order", "Count", "Enrichment")

# Set the rownames using "rownames"
rownames(Glutamatergic_Up_Clustered_Matrix) <- c("Cluster1", "Cluster2", "Cluster3")

# Turn the matrix to a data frame using "as.data.frame"
Glutamatergic_Up_Clustered_DF <- as.data.frame(Glutamatergic_Up_Clustered_Matrix)

# View "Glutamatergic_Up_Clustered_DF" before combining with Unclustered
Glutamatergic_Up_Clustered_DF

# Extract the unclustered terms
Glutamatergic_Up_Unclustered <- dplyr::filter(Glutamatergic_Only_Up, Groups == "Unclustered Up")

# Select only the counts
Glutamatergic_Up_Unclustered_Count <- paste0(Glutamatergic_Up_Unclustered$Count, collapse = ",")

# Turn the comma separated values to a list
Glutamatergic_Up_Unclustered_Count <- as.list(strsplit(Glutamatergic_Up_Unclustered_Count, ",")[[1]])

# Repeat for the enrichment scores
Glutamatergic_Up_Unclustered_Enrichment <- paste0(Glutamatergic_Up_Unclustered$Fold_Enrichment, collapse = ",")
Glutamatergic_Up_Unclustered_Enrichment <- as.list(strsplit(Glutamatergic_Up_Unclustered_Enrichment, ",")[[1]])

# Make a matrix of the values
Glutamatergic_Up_Unclustered_Matrix <- matrix(data = c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 
                                             Glutamatergic_Up_Unclustered_Count, 
                                             Glutamatergic_Up_Unclustered_Enrichment),
                                  nrow = 13, byrow = FALSE)

# Set the column names using "colnames"
colnames(Glutamatergic_Up_Unclustered_Matrix) <- c("Order", "Count", "Enrichment")

# Set the rownames using "rownames"
rownames(Glutamatergic_Up_Unclustered_Matrix) <- Glutamatergic_Up_Unclustered$Term

# Turn the matrix to a data frame using "as.data.frame"
Glutamatergic_Up_Unclustered_DF <- as.data.frame(Glutamatergic_Up_Unclustered_Matrix)

# View the data frame
Glutamatergic_Up_Unclustered_DF
#                                                      Order Count         Enrichment
# KW-0346~Stress response                                 4    16  4.021795003598232
# KW-0358~Heparin-binding                                 5    13 3.2553513798327076
# KW-0675~Receptor                                        6    65 1.4607345935146765
# KW-0090~Biological rhythms                              7    18  2.311085438218768
# KW-0037~Angiogenesis                                    8    15   2.52547858421705
# KW-0904~Protein phosphatase                             9    17 2.2537048014226437
# KW-0343~GTPase activation                              10    18 1.8805748674122724
# KW-0805~Transcription regulation                       11   100  1.200990264365219
# mmu04392:Hippo signaling pathway - multiple species    12     5  3.908673951166869
# KW-0276~Fatty acid metabolism                          13    13 1.9174156964468843
# KW-0143~Chaperone                                      14    20  1.618044472816257
# mmu04710:Circadian rhythm                              15     6   3.07808073654391
# KW-0804~Transcription                                  16   101 1.1704663949189498

# Combine clustered and unclustered terms
Glutamatergic_Up_Plot <- rbind(Glutamatergic_Up_Clustered_DF, Glutamatergic_Up_Unclustered_DF)

# View the data frame
Glutamatergic_Up_Plot
#                                                      Order Count         Enrichment
# Cluster1                                                1    83           1.640971
# Cluster2                                                2    99           2.113351
# Cluster3                                                3   177           2.173767
# KW-0346~Stress response                                 4    16  4.021795003598232
# KW-0358~Heparin-binding                                 5    13 3.2553513798327076
# KW-0675~Receptor                                        6    65 1.4607345935146765
# KW-0090~Biological rhythms                              7    18  2.311085438218768
# KW-0037~Angiogenesis                                    8    15   2.52547858421705
# KW-0904~Protein phosphatase                             9    17 2.2537048014226437
# KW-0343~GTPase activation                              10    18 1.8805748674122724
# KW-0805~Transcription regulation                       11   100  1.200990264365219
# mmu04392:Hippo signaling pathway - multiple species    12     5  3.908673951166869
# KW-0276~Fatty acid metabolism                          13    13 1.9174156964468843
# KW-0143~Chaperone                                      14    20  1.618044472816257
# mmu04710:Circadian rhythm                              15     6   3.07808073654391
# KW-0804~Transcription                                  16   101 1.1704663949189498

# Add a column to show what the cluster is or unclustered terms
Glutamatergic_Groups_Up <- rep(c("Cluster1", "Cluster2", "Cluster3", "Unclustered"), times = c(1,1,1,13))
# Add to data frame that will be used for the plot
Glutamatergic_Up_Plot$Groups <- Glutamatergic_Groups_Up
# Add rownames as the terms
Glutamatergic_Up_Terms <- rownames(Glutamatergic_Up_Plot)
# Add to data frame that will be used for the plot
Glutamatergic_Up_Plot$Terms <- Glutamatergic_Up_Terms
# Verify data class is a data frame
data.class(Glutamatergic_Up_Plot)

## Downregulated
# Now we will repeat the steps with the downregulated glutamatergic terms
# The order of the columns of the matrix will be "Order", "Count" then "Enrichment"
# The rows of the matrix will be in order of the clusters
Glutamatergic_Down_Clustered_Matrix <- matrix(data = c(1, 2, length(Glutamatergic_Cluster1_Down_Genes), 
                                         length(Glutamatergic_Cluster2_Down_Genes), Glutamatergic_Cluster1_Down_GM, 
                                         Glutamatergic_Cluster2_Down_GM), nrow = 2, byrow = FALSE)

# Set the column names using "colnames"
colnames(Glutamatergic_Down_Clustered_Matrix) <- c("Order", "Count", "Enrichment")

# Set the rownames using "rownames"
rownames(Glutamatergic_Down_Clustered_Matrix) <- c("Cluster1", "Cluster2")

# Turn the matrix to a data frame using "as.data.frame"
Glutamatergic_Down_Clustered_DF <- as.data.frame(Glutamatergic_Down_Clustered_Matrix)

# View "Down_Clustered_DF"
Glutamatergic_Down_Clustered_DF

# Add a column to show what the cluster is or unclustered terms
Glutamatergic_Groups_Down <- rep(c("Cluster1","Cluster2"), times = c(1,1))
# Add to data frame that will be used for the plot
Glutamatergic_Down_Clustered_DF$Groups <- Glutamatergic_Groups_Down
# Add rownames as the terms
Glutamatergic_Terms_Down <- rownames(Glutamatergic_Down_Clustered_DF)
# Add to data frame that will be used for the plot
Glutamatergic_Down_Clustered_DF$Terms <- Glutamatergic_Terms_Down
# Verify data class is a data frame
data.class(Glutamatergic_Down_Clustered_DF)
# [1] "data.frame"

#### Change necessary variables to integers and numerics for plotting ####
##### Glutamatergic Only #####
## Upregulated
# Order (needs to be an integer)
Glutamatergic_Up_Plot[, 1] <- as.integer(Glutamatergic_Up_Plot[, 1])
# Count (needs to be a numeric)
Glutamatergic_Up_Plot[, 2] <- as.numeric(Glutamatergic_Up_Plot[, 2])
# Enrichment (needs to be a numeric)
Glutamatergic_Up_Plot[, 3] <- as.numeric(Glutamatergic_Up_Plot[, 3])

## Downregulated
# Order (needs to be an integer)
Glutamatergic_Down_Clustered_DF[, 1] <- as.integer(Glutamatergic_Down_Clustered_DF[, 1])
# Count (needs to be a numeric)
Glutamatergic_Down_Clustered_DF[, 2] <- as.numeric(Glutamatergic_Down_Clustered_DF[, 2])
# Enrichment (needs to be a numeric)
Glutamatergic_Down_Clustered_DF[, 3] <- as.numeric(Glutamatergic_Down_Clustered_DF[, 3])

##### Shared between Glutamatergic and GABAergic #####
# Count (needs to be a numeric)
Glutamatergic_GABAergic_Up_Unclustered[, 3] <- as.numeric(Glutamatergic_GABAergic_Up_Unclustered[, 3])
# Fold Enrichment (needs to be a numeric)
Glutamatergic_GABAergic_Up_Unclustered[, 10] <- as.numeric(Glutamatergic_GABAergic_Up_Unclustered[, 10])

# Count (needs to be a numeric)
Glutamatergic_GABAergic_Down_Unclustered[, 3] <- as.numeric(Glutamatergic_GABAergic_Down_Unclustered[, 3])
# Enrichment (needs to be a numeric)
Glutamatergic_GABAergic_Down_Unclustered[, 10] <- as.numeric(Glutamatergic_GABAergic_Down_Unclustered[, 10])

##### GABAergic Only #####
# Count (needs to be a numeric)
GABAergic_Only_Up_Unclustered[, 3] <- as.numeric(GABAergic_Only_Up_Unclustered[, 3])
# Enrichment (needs to be a numeric)
GABAergic_Only_Up_Unclustered[, 10] <- as.numeric(GABAergic_Only_Up_Unclustered[, 10])

# Count (needs to be a numeric)
GABAergic_Only_Down_Unclustered[, 3] <- as.numeric(GABAergic_Only_Down_Unclustered[, 3])
# Enrichment (needs to be a numeric)
GABAergic_Only_Down_Unclustered[, 10] <- as.numeric(GABAergic_Only_Down_Unclustered[, 10])

# Can view individual structure of the R object
str(Glutamatergic_Up_Plot)
str(Glutamatergic_Down_Clustered_DF)

str(Glutamatergic_GABAergic_Up_Unclustered)
str(Glutamatergic_GABAergic_Down_Unclustered)

str(GABAergic_Only_Up_Unclustered)
str(GABAergic_Only_Down_Unclustered)

#### Make Bubble Plots ####
##### Glutamatergic Only #####
Glutamatergic_Up_Plot_Condensed <- ggplot2::ggplot(data = Glutamatergic_Up_Plot, mapping = aes(x= Enrichment, # x axis
                                                                   y= reorder(Terms, -Order), # y axis
                                                                   size = Count)) + # size of bubbles based on the number of genes in term
  geom_point(colour = "red3") + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,200),breaks=c(50, 100, 150, 200),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  coord_cartesian(xlim = c(1,18)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text.y = element_text(size=7), axis.text.x = element_blank(),
        legend.text = element_text(size=7), axis.ticks.x = element_blank(),
        axis.title = element_blank(), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

Glutamatergic_Up_Plot_Condensed

Glutamatergic_Down_Plot_Condensed <- ggplot2::ggplot(data = Glutamatergic_Down_Clustered_DF, mapping = aes(x= Enrichment, # x axis
                                                                                  y= reorder(Terms, -Order), # y axis
                                                                                  size = Count)) + # size of bubbles based on the number of genes in term
  geom_point(colour = "dodgerblue3") + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,200),breaks=c(50, 100, 150, 200),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  coord_cartesian(xlim = c(1,18)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text = element_text(size=7), legend.text = element_text(size=7), 
        axis.title = element_text(size=7), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

Glutamatergic_Down_Plot_Condensed

##### Shared between Glutamatergic and GABAergic #####
# Not currently included
Glutamatergic_GABAergic_Up_Plot <- ggplot2::ggplot(data = Glutamatergic_GABAergic_Up_Unclustered, mapping = aes(x= Fold_Enrichment, # x axis
                                                                                                           y= reorder(Term, -Order), # y axis
                                                                                                           size = Count)) + # size of bubbles based on the number of genes in term
  geom_point(colour = "red3") + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,200),breaks=c(50, 100, 150, 200),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  coord_cartesian(xlim = c(1,18)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text.y = element_text(size=7), axis.text.x = element_blank(),
        legend.text = element_text(size=7), axis.ticks.x = element_blank(),
        axis.title = element_blank(), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

Glutamatergic_GABAergic_Up_Plot

# Not currently included
Glutamatergic_GABAergic_Down_Plot <- ggplot2::ggplot(data = Glutamatergic_GABAergic_Down_Unclustered, mapping = aes(x= Fold_Enrichment, # x axis
                                                                                                                y= reorder(Term, -Order), # y axis
                                                                                                                size = Count)) + # size of bubbles based on the number of genes in term
  geom_point(colour = "dodgerblue3") + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,200),breaks=c(50, 100, 150, 200),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  coord_cartesian(xlim = c(1,18)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text = element_text(size=7), legend.text = element_text(size=7), 
        axis.title = element_text(size=7), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

Glutamatergic_GABAergic_Down_Plot

##### GABAergic Only #####
GABAergic_Only_Up_Plot <- ggplot2::ggplot(data = GABAergic_Only_Up_Unclustered, mapping = aes(x= Fold_Enrichment, # x axis
                                                                                                                y= reorder(Term, -Order), # y axis
                                                                                                                size = Count)) + # size of bubbles based on the number of genes in term
  geom_point(colour = "red3") + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,200),breaks=c(50, 100, 150, 200),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  coord_cartesian(xlim = c(1,18)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text.y = element_text(size=7), axis.text.x = element_blank(),
        legend.text = element_text(size=7), axis.ticks.x = element_blank(),
        axis.title = element_blank(), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

GABAergic_Only_Up_Plot

GABAergic_Only_Down_Plot <- ggplot2::ggplot(data = GABAergic_Only_Down_Unclustered, mapping = aes(x= Fold_Enrichment, # x axis
                                                                                                                    y= reorder(Term, -Order), # y axis
                                                                                                                    size = Count)) + # size of bubbles based on the number of genes in term
  geom_point(colour = "dodgerblue3") + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,200),breaks=c(50, 100, 150, 200),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  coord_cartesian(xlim = c(1,18)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text = element_text(size=7), legend.text = element_text(size=7), 
        axis.title = element_text(size=7), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

GABAergic_Only_Down_Plot

##### Ensure the plots have the same sized grids #####
grid.newpage()
grid.draw(rbind(ggplotGrob(GABAergic_Only_Up_Plot),
                ggplotGrob(GABAergic_Only_Down_Plot),
                ggplotGrob(Glutamatergic_Up_Plot_Condensed), 
                ggplotGrob(Glutamatergic_Down_Plot_Condensed),
                size="first"))


