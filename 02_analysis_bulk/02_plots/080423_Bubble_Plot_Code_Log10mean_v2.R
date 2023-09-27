#### Bubble Plot Code for DTU Analysis ####

# The DTU list inputted into DAVID contained the 1,575 genes that have significant DTU ("072423_DTU_Genes_for_DAVID.txt")
# Of these 1,575 genes, 1,553 were DAVID IDs

# For the background, I took all of the genes (43,697 transcripts --> 10,792 genes)
# remaining after the isoformProportions step during differential expression analysis ("072423_DTU_Background_for_DAVID.txt")

# For this analysis, I used R version 4.2.2

#### Load packages and set the working directory ####
library(ggplot2) # Version 3.4.0
library(RColorBrewer) # Version 1.1-3
library(shades) # Version 1.4.0
library(ggforce) # Version 0.4.1

# The working directory is the default location where files will be referenced from and saved to
setwd("~/Dropbox/Sleep_RNAseq_Splicing/BulkRNAseq_SD/Final_Tximeta")

#### Load R file with code (run it line by line) or copy and paste commands form this file ####

#### Read the text files that were exported from DAVID ####
# As the text file is tab delimited, we define the field separator as "\t"
# fill = TRUE adds blank spaces to uneven rows

# DTU, Up and Down
DTU_Clustered <- read.table("080223_DAVID_DTU/080223_Output_from_DAVID/080223_DTU_Clustered.txt", sep = "\t", fill = TRUE)
DTU_Not_Clustered <- read.table("080223_DAVID_DTU/080223_Output_from_DAVID/080223_DTU_Unclustered.txt", sep = "\t", fill = TRUE)

#### Record dimensions of matrix ####
# After loading the text files into R, you can view the entire object by calling the name of the object (in this case Up_Regulated_Clustered)
# You can view the first few rows with head(object name)
# You can see the dimensions of the object with dim(object name). The first number is the number of rows and the second number is the number of columns.

dim(DTU_Clustered)
# [1] 16 13

#### View all columns and the first few rows of the file ####

head(DTU_Clustered)

#### Subset Clustered Matrix ####
# Before each cluster of terms there is a row that contains the Annotation Cluster and the Enrichment Score. For plotting purposes, we want to remove these terms.
# There are also additional headers before each cluster. We want to remove these as well.

# We will subset the data frame by specifying the rows and columns that we want to keep (dataframe[rows, columns])
DTU_Clustered <- DTU_Clustered[c(3:6, 9:11, 14:16),]
dim(DTU_Clustered)
# [1] 10 13

# We will then replace the colnames that were previously V1, V2, V3, etc with their respective name
colnames(DTU_Clustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")

# Repeat the same steps but with the non-clustered terms
dim(DTU_Not_Clustered)
# [1]  7 13
head(DTU_Not_Clustered)

DTU_Not_Clustered <- DTU_Not_Clustered[2:7,]
dim(DTU_Not_Clustered)
# [1]  6 13

# Change the colnames
# DTU
colnames(DTU_Not_Clustered) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")

# Join both clustered and non-clustered up regulated terms using rbind. Note that both data frames you wat to combine need to have the same column names.
Combined_DTU <- rbind(DTU_Clustered, DTU_Not_Clustered)
dim(Combined_DTU)
# [1] 16 13

# Set the following columns to be numerics (rather than characters) for bubble plots
# Count
Combined_DTU[, 3] <- as.numeric(Combined_DTU[, 3])
# PValue
Combined_DTU[, 5] <- as.numeric(Combined_DTU[, 5])
# Fold Enrichment
Combined_DTU[, 10] <- as.numeric(Combined_DTU[, 10])

# You can view the structure of an R object using str (can double check that the Count, PValue and Fold Enrichment are all numerics)
str(Combined_DTU)

#### Finally, add a column (Order) to preserve the order of the clusters with ggplot2 ####

# DTU
Order_DTU <- c(1:16)
Combined_DTU$Order <- Order_DTU
dim(Combined_DTU)
# [1] 16 14

data_DTU <- rep(c("Cluster1","Cluster2", "Cluster3", "Unclustered"), times = c(4,3,3,6))
groups_DTU <- matrix(data_DTU, ncol = 1, byrow = TRUE)

Combined_DTU$groups <- groups_DTU

#### Bubble plots are generated with ggplot2: DTU ####
BubblePlot <- ggplot2::ggplot(data = Combined_DTU, mapping = aes(x=Fold_Enrichment, # x axis
                                                                     y= reorder(Term, -Order), # y axis
                                                                     size = Count, # size of bubbles based on the number of genes in term
                                                                     color= P_Value)) + # the color of the bubbles is determined by the pvalue
  geom_point(alpha = 0.8) + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,200),breaks=c(50, 100, 150, 200),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  shades::lightness(scale_colour_distiller(palette = "Greys"), scalefac(0.8)) + # We incorporate a red color palate for up regulated terms. We darken the palate so lighter shades are visible against the background. 
  coord_cartesian(xlim = c(1,2.5)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text = element_text(size=7), legend.text = element_text(size=7), 
        axis.title = element_text(size=7), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

BubblePlot

#### Save Session Info if Needed ####

# sink('050423_BubblePlotSessionInfo.txt')
# sessionInfo()
# sink() 


