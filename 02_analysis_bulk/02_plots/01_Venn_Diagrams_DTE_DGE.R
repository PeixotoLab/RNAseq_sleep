#### Venn Diagrams for Gene and Transcript Expression ####

# Summary: I will download and use the Venn Diagram package to separate out
# transcript and gene expression lists so that we can determine the following:

# 1) Which genes are detected at the transcript level, but not at the gene level?

# 2) Are there transcripts that are equally up- and down-regulated whose 
# effects are canceled out when summarized to the gene level?

# I would recommend this helpful website when making a Venn Diagram: 
# https://r-graph-gallery.com/14-venn-diagramm.html

# R Version 4.2.2

setwd("/Users/kaitlynford/Dropbox/Sleep_RNAseq_Splicing/BulkRNAseq_SD/Final_Tximeta")

#### LOAD PACKAGES ####
library("VennDiagram") # 1.7.3
library("readxl") # 1.4.1
library("RColorBrewer") # 1.1-3
library("dplyr") # 1.1.1
library("grid") # 4.2.2

#### READ LISTS ####
# Read lists without a fold change cutoff
# Gene
Gene <- readxl::read_excel("Peixoto_Figure_5_Supplementary_Table_1.xlsx", 
                           sheet = "DGE")
dim(Gene)
# [1] 8505    8

# Transcript Expression
Transcript <- readxl::read_excel("Peixoto_Figure_5_Supplementary_Table_1.xlsx", 
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

#### Set Colors for Venn Diagrams ####

# Red for up-regulated
Just_Up2 <- brewer.pal(9, "Greys")[c(5,6)]

# Blue for down-regulated
Just_Down2 <- brewer.pal(9, "Greys")[c(5,6)]

#### UP REGULATED DGE DTE ####
Venn2_Up <- VennDiagram::venn.diagram(
  x = list(Gene_Up$Gene_Stable_ID, Transcript_Up$Gene_Stable_ID), # Lists$columns of interest for intersection
  filename = NULL, # What  to save  file as
  resolution = 300, # Resolution in DPI
  category.names = c("DGE", "DTE"), # Names of categories
  lwd = 1, # width of each circles circumference
  lty = 1, # outline of each circle
  col =  Just_Up2, # Color palate for outline
  cex = 1, # font size for numbers
  fontfamily = "sans", # sans font for numbers
  fontface = "bold", # numbers are bolded
  cat.cex = 1, # font size for category names
  cat.pos = c(220, 140), # position of category name, 0 is 12 o'clock, 180 is 6 o'clock
  cat.fontfamily = "sans", # sans font for category names
  cat.fontface = "bold", # bold font of category names
  fill = Just_Up2, # Color palate for fill
  alpha = 0.30 # Transparency of each circle
)

grid::grid.newpage()
grid::grid.draw(Venn2_Up)

#### Down Regulated DGE, DTE ####
Venn_Down2 <- VennDiagram::venn.diagram(
  x = list(Gene_Down$Gene_Stable_ID, Transcript_Down$Gene_Stable_ID), # Lists$columns of interest for intersection
  filename = NULL, # What  to save  file as
  resolution = 300, # Resolution in DPI
  category.names = c("DGE", "DTE"), # Names of categories
  lwd = 1, # width of each circles circumference
  lty = 1, # outline of each circle
  col =  Just_Down2, # Color palate for outline
  cex = 1, # font size for numbers
  fontfamily = "sans", # sans font for numbers
  fontface = "bold", # numbers are bolded
  cat.cex = 1, # font size for category names
  cat.pos = c(220, 140), # position of category name, 0 is 12 o'clock, 180 is 6 o'clock
  cat.fontfamily = "sans", # sans font for category names
  cat.fontface = "bold", # bold font of category names
  fill = Just_Down2, # Color palate for fill
  alpha = 0.30 # Transparency of each circle
)

grid::grid.newpage()
grid::grid.draw(Venn_Down2)

#### Save Session Info if Needed ####
# sink('040423_VennDiagram_DTE_DGE.txt')
# sessionInfo()
# sink() 


