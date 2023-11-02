#### Bubble Plot Code for L2/3 IT CTX and L4/5 IT CTX ####

# R version 4.2.2

#### Load packages and set working directory ####

# Load packages
library(ggplot2) # Version 3.4.2
library(RColorBrewer) # Version 1.1-3
library(shades) # Version 1.4.0
library(grid) # Version 4.2.2

# Set the working directory
setwd("~/Dropbox/Sleep_RNAseq_Splicing/Single_nuclear_RNAseq_SD/KF_Unique_Functional_Annotation/DAVID_Output")

#### Read the text files that were exported from DAVID ####
# As the text file is tab delimited, we define the field separator as "\t"
# "fill = TRUE" adds blank spaces to uneven rows

##### L2/3 IT CTX #####
L23_Unclustered_Up <- read.table("100523_DAVID_L23ITCTX_Unique_Upregulated_Unclustered.txt", sep = "\t", fill = TRUE)
L23_Clustered_Down <- read.table("100523_DAVID_L23ITCTX_Unique_Downregulated_Clustered.txt", sep = "\t", fill = TRUE)

##### L4/5 IT CTX #####
L45_Clustered_Up <- read.table("100523_DAVID_L45ITCTX_Unique_Upregulated_Clustered.txt", sep = "\t", fill = TRUE)
L45_Unclustered_Down <- read.table("100523_DAVID_L45ITCTX_Unique_Downregulated_Unclustered.txt", sep = "\t", fill = TRUE)

#### Record dimensions of matrix ####
# After loading the text files into R, you can view the entire object by calling the name of the object (in this case "L23_Unclustered_Up", etc.)
# You can view the first few rows with head(object name)
# You can see the dimensions of the object with dim(object name): the first number is the number of rows and the second number is the number of columns.

##### L2/3 IT CTX #####
# View the dimensions
dim(L23_Unclustered_Up)
# [1]  3 13

# Subset to remove headers
L23_Unclustered_Up <- L23_Unclustered_Up[2:3,]

# Ensure the dimensions following subsetting are correct
dim(L23_Unclustered_Up)
# [1]  2 13

# Repeat with the downregulated file 
# View the dimensions
dim(L23_Clustered_Down)
# [1]  5 13

# Subset to remove headers
L23_Clustered_Down <- L23_Clustered_Down[3:5,]
# Ensure the dimensions following subsetting are correct
dim(L23_Clustered_Down)
# [1]  3 13

##### L4/5 IT CTX #####
# View the dimensions
dim(L45_Clustered_Up)
# [1] 16 13

# Subset to remove headers
L45_Clustered_Up <- L45_Clustered_Up[3:16,]
# Ensure the dimensions following subsetting are correct
dim(L45_Clustered_Up)
# [1] 14 13

# Repeat with the downregulated file
# View the dimensions
dim(L45_Unclustered_Down)
# [1]  5 13

# Subset to remove headers
L45_Unclustered_Down <- L45_Unclustered_Down[2:5,]
# Ensure the dimensions following subsetting are correct
dim(L45_Unclustered_Down)
# [1]  4 13

#### Change the colnames ####

##### L2/3 IT CTX #####
colnames(L23_Unclustered_Up) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")
colnames(L23_Clustered_Down) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")

##### L4/5 IT CTX #####
colnames(L45_Clustered_Up) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")
colnames(L45_Unclustered_Down) <- c("Category",	"Term",	"Count",	"%",	"P_Value",	"Genes", "List_Total", "Pop_Hits", "Pop_Total",	"Fold_Enrichment",	"Bonferroni",	"Benjamini",	"FDR")

##### Set the following columns to be numerics (rather than characters) for bubble plots ####

##### L2/3 IT CTX #####
# Count
L23_Unclustered_Up[, 3] <- as.numeric(L23_Unclustered_Up[, 3])
# PValue
L23_Unclustered_Up[, 5] <- as.numeric(L23_Unclustered_Up[, 5])
# Fold Enrichment
L23_Unclustered_Up[, 10] <- as.numeric(L23_Unclustered_Up[, 10])

# Count
L23_Clustered_Down[, 3] <- as.numeric(L23_Clustered_Down[, 3])
# PValue
L23_Clustered_Down[, 5] <- as.numeric(L23_Clustered_Down[, 5])
# Fold Enrichment
L23_Clustered_Down[, 10] <- as.numeric(L23_Clustered_Down[, 10])

##### L4/5 IT CTX #####
# Count
L45_Clustered_Up[, 3] <- as.numeric(L45_Clustered_Up[, 3])
# PValue
L45_Clustered_Up[, 5] <- as.numeric(L45_Clustered_Up[, 5])
# Fold Enrichment
L45_Clustered_Up[, 10] <- as.numeric(L45_Clustered_Up[, 10])

# Count
L45_Unclustered_Down[, 3] <- as.numeric(L45_Unclustered_Down[, 3])
# PValue
L45_Unclustered_Down[, 5] <- as.numeric(L45_Unclustered_Down[, 5])
# Fold Enrichment
L45_Unclustered_Down[, 10] <- as.numeric(L45_Unclustered_Down[, 10])

# You can view the structure of an R object using str (can double check that the Count, PValue and Fold Enrichment are all numerics)
str(L23_Unclustered_Up)
str(L23_Clustered_Down)
str(L45_Clustered_Up)
str(L45_Unclustered_Down)

#### Finally, add a column (Order) to preserve the order of the clusters with ggplot2 ####

##### L2/3 IT CTX #####
Order_L23_Up <- c(1:2)
L23_Unclustered_Up$Order <- Order_L23_Up
dim(L23_Unclustered_Up)
# [1] 2 14

Terms_L23_Up <- rep(c("Unclustered"), times = c(2))
Groups_L23_Up <- matrix(Terms_L23_Up, ncol = 1, byrow = TRUE)

L23_Unclustered_Up$Groups <- Groups_L23_Up
dim(L23_Unclustered_Up)
# [1]  2 15

# Repeat with downregulated
Order_L23_Down <- c(1:3)
L23_Clustered_Down$Order <- Order_L23_Down
dim(L23_Clustered_Down)
# [1] 3 14

Terms_L23_Down <- rep(c("Cluster 1"), times = c(3))
Groups_L23_Down <- matrix(Terms_L23_Down, ncol = 1, byrow = TRUE)

L23_Clustered_Down$Groups <- Groups_L23_Down
dim(L23_Clustered_Down)
# [1]  3 15

##### L4/5 IT CTX #####
Order_L45_Up <- c(1:14)
L45_Clustered_Up$Order <- Order_L45_Up
dim(L45_Clustered_Up)
# [1] 14 14

Terms_L45_Up <- rep(c("Cluster1"), times = c(14))
Groups_L45_Up <- matrix(Terms_L45_Up, ncol = 1, byrow = TRUE)

L45_Clustered_Up$Groups <- Groups_L45_Up
dim(L45_Clustered_Up)
# [1] 14 15

# Repeat with downregulated
Order_L45_Down <- c(1:4)
L45_Unclustered_Down$Order <- Order_L45_Down
dim(L45_Unclustered_Down)
# [1]  4 14

Terms_L45_Down <- rep(c("Unclustered"), times = c(4))
Groups_L45_Down <- matrix(Terms_L45_Down, ncol = 1, byrow = TRUE)

L45_Unclustered_Down$Groups <- Groups_L45_Down
dim(L45_Unclustered_Down)
# [1]  4 15

#### Bubble plots are generated with ggplot2 ####

##### L2/3 IT CTX #####
BubblePlot_L23_Up <- ggplot2::ggplot(data = L23_Unclustered_Up, mapping = aes(x=Fold_Enrichment, # x axis
                                                                            y= reorder(Term, -Order), # y axis
                                                                            size = Count, # size of bubbles based on the number of genes in term
                                                                            color= P_Value)) + # the color of the bubbles is determined by the pvalue
  geom_point(alpha = 0.8) + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,20),breaks=c(5, 10, 15, 20),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  shades::lightness(scale_colour_distiller(palette = "Reds", limits = c(0, 0.05)), scalefac(0.8)) + # We incorporate a red color palate for up regulated terms. We darken the palate so lighter shades are visible against the background. 
  coord_cartesian(xlim = c(1,15)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text.y = element_text(size=7), axis.text.x = element_blank(),
        legend.text = element_text(size=7), axis.ticks.x = element_blank(),
        axis.title = element_blank(), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

BubblePlot_L23_Up

BubblePlot_L23_Down <- ggplot2::ggplot(data = L23_Clustered_Down, mapping = aes(x=Fold_Enrichment, # x axis
                                                                       y= reorder(Term, -Order), # y axis
                                                                       size = Count, # size of bubbles based on the number of genes in term
                                                                       color= P_Value)) + # the color of the bubbles is determined by the pvalue
  geom_point(alpha = 0.8) + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,20),breaks=c(5, 10, 15, 20),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  shades::lightness(scale_colour_distiller(palette = "Blues", limits = c(0, 0.05)), scalefac(0.8)) + # We incorporate a red color palate for up regulated terms. We darken the palate so lighter shades are visible against the background. 
  coord_cartesian(xlim = c(1,15)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text = element_text(size=7), legend.text = element_text(size=7), 
        axis.title = element_text(size=7), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

BubblePlot_L23_Down

##### L4/5 IT CTX #####
BubblePlot_L45_Up <- ggplot2::ggplot(data = L45_Clustered_Up, mapping = aes(x=Fold_Enrichment, # x axis
                                                                              y= reorder(Term, -Order), # y axis
                                                                              size = Count, # size of bubbles based on the number of genes in term
                                                                              color= P_Value)) + # the color of the bubbles is determined by the pvalue
  geom_point(alpha = 0.8) + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,20),breaks=c(5, 10, 15, 20),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  shades::lightness(scale_colour_distiller(palette = "Reds", limits = c(0, 0.05)), scalefac(0.8)) + # We incorporate a red color palate for up regulated terms. We darken the palate so lighter shades are visible against the background. 
  coord_cartesian(xlim = c(1,15)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text.y = element_text(size=7), axis.text.x = element_blank(),
        legend.text = element_text(size=7), axis.ticks.x = element_blank(),
        axis.title = element_blank(), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

BubblePlot_L45_Up

BubblePlot_L45_Down <- ggplot2::ggplot(data = L45_Unclustered_Down, mapping = aes(x=Fold_Enrichment, # x axis
                                                                                y= reorder(Term, -Order), # y axis
                                                                                size = Count, # size of bubbles based on the number of genes in term
                                                                                color= P_Value)) + # the color of the bubbles is determined by the pvalue
  geom_point(alpha = 0.8) + # alpha is the opacity of the bubbles
  scale_size_continuous(limits= c(1,20),breaks=c(5, 10, 15, 20),range = c(1,7)) + # limits sets the min and max count to be included, breaks defines the size of bubbles included on the legend, range sets the min and max size of bubbles after transformation
  shades::lightness(scale_colour_distiller(palette = "Blues", limits = c(0, 0.05)), scalefac(0.8)) + # We incorporate a red color palate for up regulated terms. We darken the palate so lighter shades are visible against the background. 
  coord_cartesian(xlim = c(1,15)) + # set min and max for x-axis
  theme_light() + # set the background color
  theme(axis.text = element_text(size=7), legend.text = element_text(size=7), 
        axis.title = element_text(size=7), legend.title = element_text(size=7)) + # text 7pt max to fit Nature Neuroscience submission guidelines
  facet_grid(Groups ~ ., scales = "free", space = "free") + 
  labs(x = "Fold Enrichment", y = NULL) # change the x label, leave y blank

BubblePlot_L45_Down

#### Organize plots into one pane ####

grid.newpage()
grid.draw(rbind(ggplotGrob(BubblePlot_L45_Up), 
                ggplotGrob(BubblePlot_L45_Down),
                ggplotGrob(BubblePlot_L23_Up), 
                ggplotGrob(BubblePlot_L23_Down),
                size="first"))

#### Save sessionInfo if needed ####
# sink('100923_L23ITCTX_L45ITCTX_BubblePlot_SessionInfo.txt')
# sessionInfo()
# sink() 
