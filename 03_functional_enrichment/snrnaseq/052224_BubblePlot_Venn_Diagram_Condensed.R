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

### The color will be a gradient and is based upon the p-values: shades of 
# red for upregulated and blue for downregulated.

##### Load packages and set the working directory ####

# Load packages
library(ggplot2) # Version 3.4.2
library(RColorBrewer) # Version 1.1-3
library(grid) # Version 4.2.2
library(dplyr) # Version 1.1.1
library(psych) # Version 2.3.6

# Set the working directory
setwd("~/Library/CloudStorage/Dropbox/Sleep_RNAseq_Splicing/Single_nuclear_RNAseq_SD/iScienceRevisionFigures/DAVID")

#### Read the text files that were exported from DAVID ####
# As the text file is tab delimited, we define the field separator as "\t"
# "fill = TRUE" adds blank spaces to uneven rows

##### Shared between Glutamatergic and GABAergic #####
Glutamatergic_GABAergic_Up_Unclustered <- read.table("052124_Shared_Up_Unclustered.txt", sep = "\t", fill = TRUE)
Glutamatergic_GABAergic_Down_Unclustered <- read.table("052124_Shared_Down_Unclustered.txt", sep = "\t", fill = TRUE)

##### Glutamatergic Only #####
Glutamatergic_Only_Up_Clustered <- read.table("052124_Glutamatergic_Only_Up_Clustered.txt", sep = "\t", fill = TRUE)
Glutamatergic_Only_Up_Unclustered <- read.table("052124_Glutamatergic_Only_Up_Unclustered.txt", sep = "\t", fill = TRUE)
Glutamatergic_Only_Down_Clustered <- read.table("052124_Glutamatergic_Only_Down_Clustered.txt", sep = "\t", fill = TRUE)
Glutamatergic_Only_Down_Unclustered <- read.table("052124_Glutamatergic_Only_Down_Unclustered.txt", sep = "\t", fill = TRUE)

##### GABAergic Only #####
GABAergic_Only_Up_Unclustered <- read.table("052124_GABAergic_Only_Up_Unclustered.txt", sep = "\t", fill = TRUE)
GABAergic_Only_Down_Unclustered <- read.table("052124_GABAergic_Only_Down_Unclustered.txt", sep = "\t", fill = TRUE)

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
# [1]  4 13

# Subset to remove headers, and keep all terms
Glutamatergic_GABAergic_Down_Unclustered <- Glutamatergic_GABAergic_Down_Unclustered[2:4,]
# Ensure the dimensions following subsetting are correct
dim(Glutamatergic_GABAergic_Down_Unclustered)
# [1]  3 13

##### Glutamatergic Only #####
# Clustered
# View the dimensions
dim(Glutamatergic_Only_Up_Clustered)
# [1] 53 13

# Subset to remove headers, and keep all terms
# Note that clustered terms have an extra header
Glutamatergic_Only_Up_Clustered <- Glutamatergic_Only_Up_Clustered[c(3:6, 9:16, 19:53),]

# Ensure the dimensions following subsetting are correct
dim(Glutamatergic_Only_Up_Clustered)
# [1] 47 13

# Unclustered
# View the dimensions
dim(Glutamatergic_Only_Up_Unclustered)
# [1] 15 13

# Subset to remove headers, and keep all terms
Glutamatergic_Only_Up_Unclustered <- Glutamatergic_Only_Up_Unclustered[c(2:15),]

# Ensure the dimensions following subsetting are correct
dim(Glutamatergic_Only_Up_Unclustered)
# [1] 14 13


# Repeat with the downregulated file 
# View the dimensions
dim(Glutamatergic_Only_Down_Clustered)
# [1] 12 13

# Subset to remove headers, and keep all terms
Glutamatergic_Only_Down_Clustered <- Glutamatergic_Only_Down_Clustered[c(3:6, 9:12),]
# Ensure the dimensions following subsetting are correct
dim(Glutamatergic_Only_Down_Clustered)
# [1]  8 13

# Unclustered
# View the dimensions
dim(Glutamatergic_Only_Down_Unclustered)
# [1] 13 13

# Subset to remove headers, and keep all terms
Glutamatergic_Only_Down_Unclustered <- Glutamatergic_Only_Down_Unclustered[c(2:13),]

# Ensure the dimensions following subsetting are correct
dim(Glutamatergic_Only_Down_Unclustered)
# [1] 12 13


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
GABAergic_Only_Down_Unclustered <- GABAergic_Only_Down_Unclustered[2,]
# Ensure the dimensions following subsetting are correct
dim(GABAergic_Only_Down_Unclustered)
# [1]  1 13

#### Change the colnames ####

##### Shared between Glutamatergic and GABAergic #####
colnames(Glutamatergic_GABAergic_Up_Unclustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")
colnames(Glutamatergic_GABAergic_Down_Unclustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")

##### Glutamatergic Only #####
colnames(Glutamatergic_Only_Up_Clustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")
colnames(Glutamatergic_Only_Up_Unclustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")
colnames(Glutamatergic_Only_Down_Clustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")
colnames(Glutamatergic_Only_Down_Unclustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")

##### GABAergic Only #####
colnames(GABAergic_Only_Up_Unclustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")
colnames(GABAergic_Only_Down_Unclustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")

#### Join both Clustered and Non-Clustered Terms Using rbind ####
# Note that both data frames you want to combine need to have the same column names.
# In this case, only the Glutamatergic intersections have both clustered and unclustered terms.
Glutamatergic_Only_Up <- rbind(Glutamatergic_Only_Up_Clustered, Glutamatergic_Only_Up_Unclustered)
dim(Glutamatergic_Only_Up)
# [1] 61 13

Glutamatergic_Only_Down <- rbind(Glutamatergic_Only_Down_Clustered, Glutamatergic_Only_Down_Unclustered)
dim(Glutamatergic_Only_Down)
# [1] 20 13

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
Glutamatergic_GABAergic_Order_Down <- c(1:3)
Glutamatergic_GABAergic_Down_Unclustered$Order <- Glutamatergic_GABAergic_Order_Down
dim(Glutamatergic_GABAergic_Down_Unclustered)
# [1]  3 14

# Add "Cluster #" to the rows and/or "Unclustered"
Glutamatergic_GABAergic_Data_Down <- rep(c("Unclustered Down"), times = 3)
Glutamatergic_GABAergic_Groups_Down <- matrix(Glutamatergic_GABAergic_Data_Down, ncol = 1, byrow = TRUE)

Glutamatergic_GABAergic_Down_Unclustered$Groups <- Glutamatergic_GABAergic_Groups_Down

##### Glutamatergic Only #####
Glutamatergic_Order_Up <- c(1:61)
Glutamatergic_Only_Up$Order <- Glutamatergic_Order_Up
dim(Glutamatergic_Only_Up)
# [1] 61 14

# Add "Cluster #" to the rows and/or "Unclustered"
Glutamatergic_Data_Up <- rep(c("Cluster 1 Up", "Cluster 2 Up", "Cluster 3 Up", "Unclustered Up"), times = c(4, 8, 35, 14))
Glutamatergic_Groups_Up <- matrix(Glutamatergic_Data_Up, ncol = 1, byrow = TRUE)

Glutamatergic_Only_Up$Groups <- Glutamatergic_Groups_Up

# Down
Glutamatergic_Order_Down <- c(1:20)
Glutamatergic_Only_Down$Order <- Glutamatergic_Order_Down
dim(Glutamatergic_Only_Down)
# [1]  20 14

# Add "Cluster #" to the rows and/or "Unclustered"
Glutamatergic_Data_Down <- rep(c("Cluster 1 Down", "Cluster 2 Down", "Unclustered Down"), times = c(4, 4, 12))
Glutamatergic_Groups_Down <- matrix(Glutamatergic_Data_Down, ncol = 1, byrow = TRUE)

Glutamatergic_Only_Down$Groups <- Glutamatergic_Groups_Down

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
GABAergic_Order_Down <- 1
GABAergic_Only_Down_Unclustered$Order <- GABAergic_Order_Down
dim(GABAergic_Only_Down_Unclustered)
# [1]  1 14

# Add "Cluster #" to the rows and/or "Unclustered"
GABAergic_Data_Down <- rep(c("Unclustered Down"), times = 1)
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
Glutamatergic_Cluster1_Down <- dplyr::filter(Glutamatergic_Only_Down, Groups == "Cluster 1 Down")
Glutamatergic_Cluster2_Down <- dplyr::filter(Glutamatergic_Only_Down, Groups == "Cluster 2 Down")

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
# [1] 59

#### Determine the geometric mean of the fold enrichments for clustered terms ####
## This is done via the psych package that was loaded above

# First, change all of the "Fold_Enrichment" values to numerics
## Upregulated
Glutamatergic_Cluster1_Up_Fold_Enrichment <- as.numeric(Glutamatergic_Cluster1_Up$Fold_Enrichment)
# [1] 2.248857 1.658392 1.398890 1.398890

Glutamatergic_Cluster2_Up_Fold_Enrichment <- as.numeric(Glutamatergic_Cluster2_Up$Fold_Enrichment)
# [1] 2.146840 1.974143 2.046000 2.027933 2.506329 2.323944 2.291667 1.659218

Glutamatergic_Cluster3_Up_Fold_Enrichment <- as.numeric(Glutamatergic_Cluster3_Up$Fold_Enrichment)
# [1] 2.904000 2.146840 2.876147 2.046000 2.613861 2.664178 1.490265 2.027933 2.104027 2.584337
# [11] 1.575179 2.062500 2.410112 2.078740 1.833333 2.414634 2.158879 2.452703 3.500000 2.357143
# [21] 2.426471 1.948819 1.757396 2.186747 2.004673 2.004673 1.742236 2.000000 1.898230 2.444444
# [31] 1.941176 2.016667 2.115385 1.804688 2.088608

## Downregulated
Glutamatergic_Cluster1_Down_Fold_Enrichment <- as.numeric(Glutamatergic_Cluster1_Down$Fold_Enrichment)
# [1] 2.681431 1.704701 2.080263 1.850660

Glutamatergic_Cluster2_Down_Fold_Enrichment <- as.numeric(Glutamatergic_Cluster2_Down$Fold_Enrichment)
# [1] 1.566745 1.396942 1.396942 1.920218

## Upregulated
Glutamatergic_Cluster1_Up_GM <- geometric.mean(Glutamatergic_Cluster1_Up_Fold_Enrichment)
# [1] 1.64363

Glutamatergic_Cluster2_Up_GM <- geometric.mean(Glutamatergic_Cluster2_Up_Fold_Enrichment)
# [1] 2.107691

Glutamatergic_Cluster3_Up_GM <- geometric.mean(Glutamatergic_Cluster3_Up_Fold_Enrichment)
# [1] 2.1574

## Downregulated
Glutamatergic_Cluster1_Down_GM <- geometric.mean(Glutamatergic_Cluster1_Down_Fold_Enrichment)
# [1] 2.048165

Glutamatergic_Cluster2_Down_GM <- geometric.mean(Glutamatergic_Cluster2_Down_Fold_Enrichment)
# [1] 1.556598


#### Determine the geometric mean of the pvalues for clustered terms ####
## This is done via the psych package that was loaded above
# First change the "P_Values" to numerics
## Upregulated
Glutamatergic_Cluster1_Up_Pvalue <- as.numeric(Glutamatergic_Cluster1_Up$P_Value)
# [1] 5.557905e-05 4.545762e-04 6.190082e-03 6.190082e-03

Glutamatergic_Cluster2_Up_Pvalue <- as.numeric(Glutamatergic_Cluster2_Up$P_Value)
# [1] 2.641545e-05 1.022559e-04 2.061992e-04 2.455299e-03 7.456094e-03 2.616586e-02 2.837519e-02
# [8] 4.170599e-02

Glutamatergic_Cluster3_Up_Pvalue <- as.numeric(Glutamatergic_Cluster3_Up$P_Value)
# [1] 1.528883e-05 2.641545e-05 7.867332e-05 2.061992e-04 9.836879e-04 1.261708e-03 2.223204e-03
# [8] 2.455299e-03 3.575598e-03 3.882135e-03 4.013194e-03 5.788112e-03 6.870253e-03 9.224919e-03
# [15] 9.771811e-03 9.830449e-03 1.186210e-02 1.288588e-02 1.289134e-02 1.677188e-02 2.027161e-02
# [22] 2.047037e-02 2.584813e-02 2.704764e-02 2.717184e-02 2.717184e-02 3.316363e-02 3.564846e-02
# [29] 3.928161e-02 4.242018e-02 4.296170e-02 4.395989e-02 4.442168e-02 4.461293e-02 4.758662e-02

## Downregulated
Glutamatergic_Cluster1_Down_Pvalue <- as.numeric(Glutamatergic_Cluster1_Down$P_Value)
# [1] 2.345231e-05 2.823659e-05 1.125197e-03 2.780638e-03

Glutamatergic_Cluster2_Down_Pvalue <- as.numeric(Glutamatergic_Cluster2_Down$P_Value)
# [1] 0.004646965 0.013999545 0.013999545 0.023785451

## Upregulated
Glutamatergic_Cluster1_Up_Pvalue_GM <- geometric.mean(Glutamatergic_Cluster1_Up_Pvalue)
# [1] 0.0009919223

Glutamatergic_Cluster2_Up_Pvalue_GM <- geometric.mean(Glutamatergic_Cluster2_Up_Pvalue)
# [1] 0.002053124

Glutamatergic_Cluster3_Up_Pvalue_GM <- geometric.mean(Glutamatergic_Cluster3_Up_Pvalue)
# [1] 0.006888918

## Downregulated
Glutamatergic_Cluster1_Down_Pvalue_GM <- geometric.mean(Glutamatergic_Cluster1_Down_Pvalue)
# [1] 0.0002133501

Glutamatergic_Cluster2_Down_Pvalue_GM <- geometric.mean(Glutamatergic_Cluster2_Down_Pvalue)
# [1] 0.01213185

#### Make matrices for plotting ####
## Upregulated
# The order of the columns of the matrix will be "Order", "Count" then "Enrichment"
# The rows of the matrix will be in order of the clusters
Glutamatergic_Up_Clustered_Matrix <- matrix(data = c(1, 2, 3, length(Glutamatergic_Cluster1_Up_Genes), 
                                       length(Glutamatergic_Cluster2_Up_Genes), length(Glutamatergic_Cluster3_Up_Genes), 
                                       Glutamatergic_Cluster1_Up_GM, Glutamatergic_Cluster2_Up_GM, Glutamatergic_Cluster3_Up_GM,
                                       Glutamatergic_Cluster1_Up_Pvalue_GM, Glutamatergic_Cluster2_Up_Pvalue_GM,
                                       Glutamatergic_Cluster3_Up_Pvalue_GM), 
                                       nrow = 3, byrow = FALSE)

# Set the column names using "colnames"
colnames(Glutamatergic_Up_Clustered_Matrix) <- c("Order", "Count", "Enrichment", "P_Value")

# Set the rownames using "rownames"
rownames(Glutamatergic_Up_Clustered_Matrix) <- c("Cluster1", "Cluster2", "Cluster3")

# Turn the matrix to a data frame using "as.data.frame"
Glutamatergic_Up_Clustered_DF <- as.data.frame(Glutamatergic_Up_Clustered_Matrix)

# View "Glutamatergic_Up_Clustered_DF" before combining with Unclustered
Glutamatergic_Up_Clustered_DF
#          Order Count Enrichment      P_Value
# Cluster1     1    83   1.643630 0.0009919223
# Cluster2     2    99   2.107691 0.0020531239
# Cluster3     3   177   2.157400 0.0068889184

# Extract the unclustered terms
Glutamatergic_Up_Unclustered <- dplyr::filter(Glutamatergic_Only_Up, Groups == "Unclustered Up")

# Select only the counts
Glutamatergic_Up_Unclustered_Count <- paste0(Glutamatergic_Up_Unclustered$Count, collapse = ",")

# Turn the comma separated values to a list
Glutamatergic_Up_Unclustered_Count <- as.list(strsplit(Glutamatergic_Up_Unclustered_Count, ",")[[1]])

# Repeat for the enrichment scores
Glutamatergic_Up_Unclustered_Enrichment <- paste0(Glutamatergic_Up_Unclustered$Fold_Enrichment, collapse = ",")
Glutamatergic_Up_Unclustered_Enrichment <- as.list(strsplit(Glutamatergic_Up_Unclustered_Enrichment, ",")[[1]])

# Repeat for Pvalues
Glutamatergic_Up_Unclustered_Pvalues <- paste0(Glutamatergic_Up_Unclustered$P_Value, collapse = ",")
Glutamatergic_Up_Unclustered_Pvalues <- as.list(strsplit(Glutamatergic_Up_Unclustered_Pvalues, ",")[[1]])

# Make a matrix of the values
Glutamatergic_Up_Unclustered_Matrix <- matrix(data = c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
                                             Glutamatergic_Up_Unclustered_Count, 
                                             Glutamatergic_Up_Unclustered_Enrichment,
                                             Glutamatergic_Up_Unclustered_Pvalues),
                                             nrow = 14, byrow = FALSE)

# Set the column names using "colnames"
colnames(Glutamatergic_Up_Unclustered_Matrix) <- c("Order", "Count", "Enrichment", "P_Value")

# Set the rownames using "rownames"
rownames(Glutamatergic_Up_Unclustered_Matrix) <- Glutamatergic_Up_Unclustered$Term

# Turn the matrix to a data frame using "as.data.frame"
Glutamatergic_Up_Unclustered_DF <- as.data.frame(Glutamatergic_Up_Unclustered_Matrix)

# View the data frame
Glutamatergic_Up_Unclustered_DF
#                                                     Order Count   Enrichment      P_Value
# KW-0346~Stress response                                 4    16 4.088062.... 5.265754....
# KW-0358~Heparin-binding                                 5    13 3.261855.... 4.839213....
# KW-0675~Receptor                                        6    65 1.461571.... 0.001494....
# KW-0090~Biological rhythms                              7    18 2.349165.... 0.001534....
# KW-0037~Angiogenesis                                    8    15 2.567091.... 0.001888....
# KW-0904~Protein phosphatase                             9    17 2.258207.... 0.003159....
# KW-0343~GTPase activation                              10    18 1.871935.... 0.015007....
# KW-0650~Protein phosphatase inhibitor                  11     6 3.951863.... 0.015347....
# mmu04392:Hippo signaling pathway - multiple species    12     5 3.928571.... 0.034648....
# KW-0805~Transcription regulation                       13   100 1.188001.... 0.035543....
# KW-0276~Fatty acid metabolism                          14    13 1.917317.... 0.037479....
# KW-0143~Chaperone                                      15    20 1.621277.... 0.038566....
# mmu04710:Circadian rhythm                              16     6      3.09375 0.041046....
# KW-0649~Protein kinase inhibitor                       17     4 4.863831.... 0.044495....

# Combine clustered and unclustered terms
Glutamatergic_Up_Plot <- rbind(Glutamatergic_Up_Clustered_DF, Glutamatergic_Up_Unclustered_DF)

# View the data frame
Glutamatergic_Up_Plot
#                                                      Order Count   Enrichment      P_Value
# Cluster1                                                1    83 1.643629.... 0.000991....
# Cluster2                                                2    99 2.107691.... 0.002053....
# Cluster3                                                3   177 2.157399.... 0.006888....
# KW-0346~Stress response                                 4    16 4.088062.... 5.265754....
# KW-0358~Heparin-binding                                 5    13 3.261855.... 4.839213....
# KW-0675~Receptor                                        6    65 1.461571.... 0.001494....
# KW-0090~Biological rhythms                              7    18 2.349165.... 0.001534....
# KW-0037~Angiogenesis                                    8    15 2.567091.... 0.001888....
# KW-0904~Protein phosphatase                             9    17 2.258207.... 0.003159....
# KW-0343~GTPase activation                              10    18 1.871935.... 0.015007....
# KW-0650~Protein phosphatase inhibitor                  11     6 3.951863.... 0.015347....
# mmu04392:Hippo signaling pathway - multiple species    12     5 3.928571.... 0.034648....
# KW-0805~Transcription regulation                       13   100 1.188001.... 0.035543....
# KW-0276~Fatty acid metabolism                          14    13 1.917317.... 0.037479....
# KW-0143~Chaperone                                      15    20 1.621277.... 0.038566....
# mmu04710:Circadian rhythm                              16     6      3.09375 0.041046....
# KW-0649~Protein kinase inhibitor                       17     4 4.863831.... 0.044495....

# Add a column to show what the cluster is or unclustered terms
Glutamatergic_Groups_Up <- rep(c("Cluster1", "Cluster2", "Cluster3", "Unclustered"), times = c(1,1,1,14))
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
                                         Glutamatergic_Cluster2_Down_GM, Glutamatergic_Cluster1_Down_Pvalue_GM,
                                         Glutamatergic_Cluster2_Down_Pvalue_GM), nrow = 2, byrow = FALSE)

# Set the column names using "colnames"
colnames(Glutamatergic_Down_Clustered_Matrix) <- c("Order", "Count", "Enrichment", "P_Value")

# Set the rownames using "rownames"
rownames(Glutamatergic_Down_Clustered_Matrix) <- c("Cluster1", "Cluster2")

# Turn the matrix to a data frame using "as.data.frame"
Glutamatergic_Down_Clustered_DF <- as.data.frame(Glutamatergic_Down_Clustered_Matrix)

# View "Down_Clustered_DF"
Glutamatergic_Down_Clustered_DF
#           Order Count Enrichment      P_Value
# Cluster1     1    67   2.048165 0.0002133501
# Cluster2     2    59   1.556598 0.0121318537

# Extract the unclustered terms
Glutamatergic_Down_Unclustered <- dplyr::filter(Glutamatergic_Only_Down, Groups == "Unclustered Down")

# Select only the counts
Glutamatergic_Down_Unclustered_Count <- paste0(Glutamatergic_Down_Unclustered$Count, collapse = ",")

# Turn the comma separated values to a list
Glutamatergic_Down_Unclustered_Count <- as.list(strsplit(Glutamatergic_Down_Unclustered_Count, ",")[[1]])

# Repeat for the enrichment scores
Glutamatergic_Down_Unclustered_Enrichment <- paste0(Glutamatergic_Down_Unclustered$Fold_Enrichment, collapse = ",")
Glutamatergic_Down_Unclustered_Enrichment <- as.list(strsplit(Glutamatergic_Down_Unclustered_Enrichment, ",")[[1]])

# Repeat for Pvalues
Glutamatergic_Down_Unclustered_Pvalues <- paste0(Glutamatergic_Down_Unclustered$P_Value, collapse = ",")
Glutamatergic_Down_Unclustered_Pvalues <- as.list(strsplit(Glutamatergic_Down_Unclustered_Pvalues, ",")[[1]])

# Make a matrix of the values
Glutamatergic_Down_Unclustered_Matrix <- matrix(data = c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
                                                       Glutamatergic_Down_Unclustered_Count, 
                                                       Glutamatergic_Down_Unclustered_Enrichment,
                                                       Glutamatergic_Down_Unclustered_Pvalues),
                                              nrow = 12, byrow = FALSE)

# Set the column names using "colnames"
colnames(Glutamatergic_Down_Unclustered_Matrix) <- c("Order", "Count", "Enrichment", "P_Value")

# Set the rownames using "rownames"
rownames(Glutamatergic_Down_Unclustered_Matrix) <- Glutamatergic_Down_Unclustered$Term

# Turn the matrix to a data frame using "as.data.frame"
Glutamatergic_Down_Unclustered_DF <- as.data.frame(Glutamatergic_Down_Unclustered_Matrix)

# View the data frame
Glutamatergic_Down_Unclustered_DF

#                                                       Order Count   Enrichment      P_Value
# KW-0119~Carbohydrate metabolism                          4    10 3.133489.... 0.004144....
# KW-0716~Sensory transduction                             5     9 2.862870.... 0.012373....
# mmu05033:Nicotine addiction                              6     6 4.043095.... 0.014673....
# KW-0469~Meiosis                                          7     9 2.699277.... 0.017277....
# KW-0511~Multifunctional enzyme                           8     9 2.596448.... 0.021235....
# KW-0720~Serine protease                                  9     8 2.712858.... 0.026366....
# mmu04740:Olfactory transduction                         10     7 3.018844.... 0.026553....
# mmu04913:Ovarian steroidogenesis                        11     6 3.404712.... 0.029246....
# mmu04974:Protein digestion and absorption               12     8 2.613718.... 0.031615....
# KW-0844~Vision                                          13     6 3.314902.... 0.032524....
# KW-0505~Motor protein                                   14    11 2.024954.... 0.043880....
# mmu00520:Amino sugar and nucleotide sugar metabolism    15     6 3.008815.... 0.046618....

# Combine clustered and unclustered terms
Glutamatergic_Down_Plot <- rbind(Glutamatergic_Down_Clustered_DF, Glutamatergic_Down_Unclustered_DF)

# View the data frame
Glutamatergic_Down_Plot
#                                                       Order Count   Enrichment      P_Value
# Cluster1                                                 1    67 2.048165.... 0.000213....
# Cluster2                                                 2    59 1.556597.... 0.012131....
# KW-0119~Carbohydrate metabolism                          4    10 3.133489.... 0.004144....
# KW-0716~Sensory transduction                             5     9 2.862870.... 0.012373....
# mmu05033:Nicotine addiction                              6     6 4.043095.... 0.014673....
# KW-0469~Meiosis                                          7     9 2.699277.... 0.017277....
# KW-0511~Multifunctional enzyme                           8     9 2.596448.... 0.021235....
# KW-0720~Serine protease                                  9     8 2.712858.... 0.026366....
# mmu04740:Olfactory transduction                         10     7 3.018844.... 0.026553....
# mmu04913:Ovarian steroidogenesis                        11     6 3.404712.... 0.029246....
# mmu04974:Protein digestion and absorption               12     8 2.613718.... 0.031615....
# KW-0844~Vision                                          13     6 3.314902.... 0.032524....
# KW-0505~Motor protein                                   14    11 2.024954.... 0.043880....
# mmu00520:Amino sugar and nucleotide sugar metabolism    15     6 3.008815.... 0.046618....

# Add a column to show what the cluster is or unclustered terms
Glutamatergic_Groups_Down <- rep(c("Cluster1","Cluster2","Unclustered"), times = c(1,1,12))
# Add to data frame that will be used for the plot
Glutamatergic_Down_Plot$Groups <- Glutamatergic_Groups_Down
# Add rownames as the terms
Glutamatergic_Terms_Down <- rownames(Glutamatergic_Down_Plot)
# Add to data frame that will be used for the plot
Glutamatergic_Down_Plot$Terms <- Glutamatergic_Terms_Down
# Verify data class is a data frame
data.class(Glutamatergic_Down_Plot)
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
# Pvalue (needs to be a numeric)
Glutamatergic_Up_Plot[, 4] <- as.numeric(Glutamatergic_Up_Plot[, 4])

## Downregulated
# Order (needs to be an integer)
Glutamatergic_Down_Plot[, 1] <- as.integer(Glutamatergic_Down_Plot[, 1])
# Count (needs to be a numeric)
Glutamatergic_Down_Plot[, 2] <- as.numeric(Glutamatergic_Down_Plot[, 2])
# Enrichment (needs to be a numeric)
Glutamatergic_Down_Plot[, 3] <- as.numeric(Glutamatergic_Down_Plot[, 3])
# Pvalue (needs to be a numeric)
Glutamatergic_Down_Plot[, 4] <- as.numeric(Glutamatergic_Down_Plot[, 4])

##### Shared between Glutamatergic and GABAergic #####
# Count (needs to be a numeric)
Glutamatergic_GABAergic_Up_Unclustered[, 3] <- as.numeric(Glutamatergic_GABAergic_Up_Unclustered[, 3])
# Pvalue (needs to be a numeric)
Glutamatergic_GABAergic_Up_Unclustered[, 5] <- as.numeric(Glutamatergic_GABAergic_Up_Unclustered[, 5])
# Fold Enrichment (needs to be a numeric)
Glutamatergic_GABAergic_Up_Unclustered[, 10] <- as.numeric(Glutamatergic_GABAergic_Up_Unclustered[, 10])

# Count (needs to be a numeric)
Glutamatergic_GABAergic_Down_Unclustered[, 3] <- as.numeric(Glutamatergic_GABAergic_Down_Unclustered[, 3])
# Pvalue (needs to be a numeric)
Glutamatergic_GABAergic_Down_Unclustered[, 5] <- as.numeric(Glutamatergic_GABAergic_Down_Unclustered[, 5])
# Enrichment (needs to be a numeric)
Glutamatergic_GABAergic_Down_Unclustered[, 10] <- as.numeric(Glutamatergic_GABAergic_Down_Unclustered[, 10])

##### GABAergic Only #####
# Count (needs to be a numeric)
GABAergic_Only_Up_Unclustered[, 3] <- as.numeric(GABAergic_Only_Up_Unclustered[, 3])
# Pvalue (needs to be a numeric)
GABAergic_Only_Up_Unclustered[, 5] <- as.numeric(GABAergic_Only_Up_Unclustered[, 5])
# Enrichment (needs to be a numeric)
GABAergic_Only_Up_Unclustered[, 10] <- as.numeric(GABAergic_Only_Up_Unclustered[, 10])

# Count (needs to be a numeric)
GABAergic_Only_Down_Unclustered[, 3] <- as.numeric(GABAergic_Only_Down_Unclustered[, 3])
# Pvalue (needs to be a numeric)
GABAergic_Only_Down_Unclustered[, 5] <- as.numeric(GABAergic_Only_Down_Unclustered[, 5])
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
                                                                   size = Count,
                                                                   color = P_Value)) + # size of bubbles based on the number of genes in term
  geom_point(alpha = 0.8) + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,200),breaks=c(50, 100, 150, 200),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  shades::lightness(scale_colour_distiller(palette = "Reds", limits = c(0, 0.05)), scalefac(0.8)) + # We incorporate a red color palate for up regulated terms. We darken the palate so lighter shades are visible against the background. 
  coord_cartesian(xlim = c(1,12)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text.y = element_text(size=7), axis.text.x = element_blank(),
        legend.text = element_text(size=7), axis.ticks.x = element_blank(),
        axis.title = element_blank(), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

# pdf("052224_Glutamatergic_Up_Plot_Condensed.pdf")
# grid::grid.newpage()
# grid::grid.draw(Glutamatergic_Up_Plot_Condensed)
# dev.off()


Glutamatergic_Down_Plot_Condensed <- ggplot2::ggplot(data = Glutamatergic_Down_Plot, mapping = aes(x= Enrichment, # x axis
                                                                                  y= reorder(Terms, -Order), # y axis
                                                                                  size = Count,
                                                                                  color = P_Value)) + # size of bubbles based on the number of genes in term
  geom_point(alpha = 0.8) + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,200),breaks=c(50, 100, 150, 200),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  shades::lightness(scale_colour_distiller(palette = "Blues", limits = c(0, 0.05)), scalefac(0.8)) + # We incorporate a blue color palate for down regulated terms. We darken the palate so lighter shades are visible against the background. 
  coord_cartesian(xlim = c(1,12)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text = element_text(size=7), legend.text = element_text(size=7), 
        axis.title = element_text(size=7), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

# pdf("052224_Glutamatergic_Down_Plot_Condensed.pdf")
# grid::grid.newpage()
# grid::grid.draw(Glutamatergic_Down_Plot_Condensed)
# dev.off()


##### Shared between Glutamatergic and GABAergic #####
# Not currently included
Glutamatergic_GABAergic_Up_Plot <- ggplot2::ggplot(data = Glutamatergic_GABAergic_Up_Unclustered, mapping = aes(x= Fold_Enrichment, # x axis
                                                                                                           y= reorder(Term, -Order), # y axis
                                                                                                           size = Count,
                                                                                                           color = P_Value)) + # size of bubbles based on the number of genes in term
  geom_point(alpha = 0.8) + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,200),breaks=c(50, 100, 150, 200),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  shades::lightness(scale_colour_distiller(palette = "Reds", limits = c(0, 0.05)), scalefac(0.8)) + # We incorporate a red color palate for up regulated terms. We darken the palate so lighter shades are visible against the background. 
  coord_cartesian(xlim = c(1,12)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text.y = element_text(size=7), axis.text.x = element_blank(),
        legend.text = element_text(size=7), axis.ticks.x = element_blank(),
        axis.title = element_blank(), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

# pdf("052224_Glutamatergic_GABAergic_Up_Plot.pdf")
# grid::grid.newpage()
# grid::grid.draw(Glutamatergic_GABAergic_Up_Plot)
# dev.off()

# Not currently included
Glutamatergic_GABAergic_Down_Plot <- ggplot2::ggplot(data = Glutamatergic_GABAergic_Down_Unclustered, mapping = aes(x= Fold_Enrichment, # x axis
                                                                                                                y= reorder(Term, -Order), # y axis
                                                                                                                size = Count,
                                                                                                                color = P_Value)) + # size of bubbles based on the number of genes in term
  geom_point(alpha = 0.8) + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,200),breaks=c(50, 100, 150, 200),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  shades::lightness(scale_colour_distiller(palette = "Blues", limits = c(0, 0.05)), scalefac(0.8)) + # We incorporate a blue color palate for down regulated terms. We darken the palate so lighter shades are visible against the background. 
  coord_cartesian(xlim = c(1,12)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text = element_text(size=7), legend.text = element_text(size=7), 
        axis.title = element_text(size=7), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

# pdf("052224_Glutamatergic_GABAergic_Down_Plot.pdf")
# grid::grid.newpage()
# grid::grid.draw(Glutamatergic_GABAergic_Down_Plot)
# dev.off()

##### GABAergic Only #####
GABAergic_Only_Up_Plot <- ggplot2::ggplot(data = GABAergic_Only_Up_Unclustered, mapping = aes(x= Fold_Enrichment, # x axis
                                                                                                                y= reorder(Term, -Order), # y axis
                                                                                                                size = Count,
                                                                                                                color = P_Value)) + # size of bubbles based on the number of genes in term
  geom_point(alpha = 0.8) + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,200),breaks=c(50, 100, 150, 200),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  shades::lightness(scale_colour_distiller(palette = "Reds", limits = c(0, 0.05)), scalefac(0.8)) + # We incorporate a red color palate for up regulated terms. We darken the palate so lighter shades are visible against the background. 
  coord_cartesian(xlim = c(1,12)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text.y = element_text(size=7), axis.text.x = element_blank(),
        legend.text = element_text(size=7), axis.ticks.x = element_blank(),
        axis.title = element_blank(), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

# pdf("052224_GABAergic_Only_Up_Plot.pdf")
# grid::grid.newpage()
# grid::grid.draw(GABAergic_Only_Up_Plot)
# dev.off()


GABAergic_Only_Down_Plot <- ggplot2::ggplot(data = GABAergic_Only_Down_Unclustered, mapping = aes(x= Fold_Enrichment, # x axis
                                                                                                                    y= reorder(Term, -Order), # y axis
                                                                                                                    size = Count,
                                                                                                                    color = P_Value)) + # size of bubbles based on the number of genes in term
  geom_point(alpha = 0.8) + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,200),breaks=c(50, 100, 150, 200),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  shades::lightness(scale_colour_distiller(palette = "Blues", limits = c(0, 0.05)), scalefac(0.8)) + # We incorporate a blue color palate for down regulated terms. We darken the palate so lighter shades are visible against the background. 
  coord_cartesian(xlim = c(1,12)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text = element_text(size=7), legend.text = element_text(size=7), 
        axis.title = element_text(size=7), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

# pdf("052224_GABAergic_Only_Down_Plot.pdf")
# grid::grid.newpage()
# grid::grid.draw(GABAergic_Only_Down_Plot)
# dev.off()


##### Ensure the plots have the same sized grids #####
pdf("052224_GABAergic_Glutamatergic.pdf", width = 7, height = 8.5)
grid.newpage()
grid.draw(rbind(ggplotGrob(GABAergic_Only_Up_Plot),
                ggplotGrob(GABAergic_Only_Down_Plot),
                ggplotGrob(Glutamatergic_Up_Plot_Condensed), 
                ggplotGrob(Glutamatergic_Down_Plot_Condensed),
                size="first"))
dev.off()
