# Venn Diagram intersecting Glutamatergic and GABAergic lists

## We are waiting to hear back if L5/6 NP CTX will be incuded
## For now, I am keeping it

## I am going to compare the intersections we get when applying a 0.05 and 0.01 FDR threshold
## If we decide to include 0.10, I will need the full list from Elena

#### Set working directory and load packages ####
setwd("~/Dropbox/Sleep_RNAseq_Splicing/Single_nuclear_RNAseq_SD/iScienceRevisionFigures")

library(RColorBrewer)
library(dplyr)
library(readxl)

#### Read sheets of significant lists ####
L23_IT_CTX <- read_excel("Peixoto_Figure_2_Supplementary_Table_1.xlsx", sheet = 1)
L45_IT_CTX <- read_excel("Peixoto_Figure_2_Supplementary_Table_1.xlsx", sheet = 2)
L5_IT_CTX <- read_excel("Peixoto_Figure_2_Supplementary_Table_1.xlsx", sheet = 3)
L5_PT_CTX <- read_excel("Peixoto_Figure_2_Supplementary_Table_1.xlsx", sheet = 4)
L56_NP_CTX <- read_excel("Peixoto_Figure_2_Supplementary_Table_1.xlsx", sheet = 5)
L6_CT_CTX <- read_excel("Peixoto_Figure_2_Supplementary_Table_1.xlsx", sheet = 6)
L6_IT_CTX <- read_excel("Peixoto_Figure_2_Supplementary_Table_1.xlsx", sheet = 7)
L6b_CTX <- read_excel("Peixoto_Figure_2_Supplementary_Table_1.xlsx", sheet = 8)
Pvalb <- read_excel("Peixoto_Figure_2_Supplementary_Table_1.xlsx", sheet = 9)
Sst <- read_excel("Peixoto_Figure_2_Supplementary_Table_1.xlsx", sheet = 10)
Vip <- read_excel("Peixoto_Figure_2_Supplementary_Table_1.xlsx", sheet = 11)

###### Filter to only keep FDR < 0.01 #####
L23_IT_CTX_0.01 <- dplyr::filter(L23_IT_CTX, FDR < 0.01)
L45_IT_CTX_0.01 <- dplyr::filter(L45_IT_CTX, FDR < 0.01)
L5_IT_CTX_0.01 <- dplyr::filter(L5_IT_CTX, FDR < 0.01)
L5_PT_CTX_0.01 <- dplyr::filter(L5_PT_CTX, FDR < 0.01)
L56_NP_CTX_0.01 <- dplyr::filter(L56_NP_CTX, FDR < 0.01)
L6_CT_CTX_0.01 <- dplyr::filter(L6_CT_CTX, FDR < 0.01)
L6_IT_CTX_0.01 <- dplyr::filter(L6_IT_CTX, FDR < 0.01)
L6b_CTX_0.01 <- dplyr::filter(L6b_CTX, FDR < 0.01)
Pvalb_0.01 <- dplyr::filter(Pvalb, FDR < 0.01)
Sst_0.01 <- dplyr::filter(Sst, FDR < 0.01)
Vip_0.01 <- dplyr::filter(Vip, FDR < 0.01)

###### Read in lists for FDR < 0.10 #####
L23_IT_CTX_0.10 <- read_excel("2024_05_21_DEG_FDR_0.10.xlsx", sheet = 1)
L45_IT_CTX_0.10 <- read_excel("2024_05_21_DEG_FDR_0.10.xlsx", sheet = 2)
L5_IT_CTX_0.10 <- read_excel("2024_05_21_DEG_FDR_0.10.xlsx", sheet = 3)
L5_PT_CTX_0.10 <- read_excel("2024_05_21_DEG_FDR_0.10.xlsx", sheet = 4)
L56_NP_CTX_0.10 <- read_excel("2024_05_21_DEG_FDR_0.10.xlsx", sheet = 5)
L6_CT_CTX_0.10 <- read_excel("2024_05_21_DEG_FDR_0.10.xlsx", sheet = 6)
L6_IT_CTX_0.10 <- read_excel("2024_05_21_DEG_FDR_0.10.xlsx", sheet = 7)
L6b_CTX_0.10 <- read_excel("2024_05_21_DEG_FDR_0.10.xlsx", sheet = 8)
Pvalb_0.10 <- read_excel("2024_05_21_DEG_FDR_0.10.xlsx", sheet = 9)
Sst_0.10 <- read_excel("2024_05_21_DEG_FDR_0.10.xlsx", sheet = 10)
Vip_0.10 <- read_excel("2024_05_21_DEG_FDR_0.10.xlsx", sheet = 11)

##### Set colors #####
Glut_Color <- brewer.pal(9, "Greys")[6]
GABA_Color <- brewer.pal(9, "Greys")[5]


#### 0.05 ####
## Glutamatergic 
# Combine Glutamatergic DEG Lists
# There are 2652 unique genes in glutamatergic cell-types
Glutamatergic <- rbind(L23_IT_CTX, L45_IT_CTX, L5_IT_CTX, L5_PT_CTX, L56_NP_CTX,
                       L6_CT_CTX, L6_IT_CTX, L6b_CTX)

Glutamatergic <- unique(Glutamatergic$Gene_Stable_ID)
length(Glutamatergic)
# [1] 2656

# Turn this into a dataframe for intersections later
GABAergic <- rbind(Pvalb, Sst, Vip)
GABAergic <- unique(GABAergic$Gene_Stable_ID)
length(GABAergic)
# [1] 525

# Change to data frame

GABAergic <- as.data.frame(GABAergic)
Glutamatergic <- as.data.frame(Glutamatergic)

##### Plot Venn Diagram #####
Venn <- VennDiagram::venn.diagram(
  x = list(Glutamatergic$Glutamatergic, GABAergic$GABAergic), # Lists$columns of interest for intersection
  filename = NULL, # What  to save  file as
  resolution = 300, # Resolution in DPI
  category.names = c("Glutamatergic", "GABAergic"), # Names of categories
  lwd = 1, # width of each circles circumference
  lty = 1, # outline of each circle
  col =  c(Glut_Color, GABA_Color), # Color palate for outline
  cex = 1, # font size for numbers
  fontfamily = "sans", # sans font for numbers
  fontface = "bold", # numbers are bolded
  cat.cex = 1, # font size for category names
  cat.pos = 0, # position of category name, 0 is 12 o'clock, 180 is 6 o'clock
  cat.fontfamily = "sans", # sans font for category names
  cat.fontface = "bold", # bold font of category names
  fill = c(Glut_Color, GABA_Color), # Color palate for fill
  alpha = 0.40, # Transparency of each circle
  euler.d = FALSE, # Prevent scaling by size
  scaled = FALSE
)

pdf("051524_VennDiagram_Glutamatergic_GABAergic.pdf")
grid::grid.newpage()
grid::grid.draw(Venn)
dev.off()

# Shared between Glutamatergic and GABAergic cell types
Shared <- intersect(Glutamatergic$Glutamatergic, GABAergic$GABAergic)
length(Shared)
# [1] 417

# Glutamaterigc Only
Glutamatergic_Only <- setdiff(Glutamatergic$Glutamatergic, GABAergic$GABAergic)
length(Glutamatergic_Only)
# [1] 2239

# GABAergic Only
GABAergic_Only <- setdiff(GABAergic$GABAergic, Glutamatergic$Glutamatergic)
length(GABAergic_Only)
# [1] 108

##### Write tables ##### 
write.table(Shared, "051423_Shared_Glutamatergic_GABAergic.txt", sep = "/t",
            col.names = FALSE, row.names = FALSE)

write.table(Glutamatergic_Only, "051423_Glutamatergic_Only.txt", sep = "/t",
            col.names = FALSE, row.names = FALSE)

write.table(GABAergic_Only, "051423_GABAergic_Only.txt", sep = "/t",
            col.names = FALSE, row.names = FALSE)


#### 0.01 ####
## Glutamatergic 
# Combine Glutamatergic DEG Lists
# There are 2652 unique genes in glutamatergic cell-types
Glutamatergic_0.01 <- rbind(L23_IT_CTX_0.01, L45_IT_CTX_0.01, L5_IT_CTX_0.01, 
                       L5_PT_CTX_0.01, L56_NP_CTX_0.01,
                       L6_CT_CTX_0.01, L6_IT_CTX_0.01, L6b_CTX_0.01)

Glutamatergic_0.01 <- unique(Glutamatergic_0.01$Gene_Stable_ID)
length(Glutamatergic_0.01)
# [1] 1394

# Turn this into a dataframe for intersections later
GABAergic_0.01 <- rbind(Pvalb_0.01, Sst_0.01, Vip_0.01)
GABAergic_0.01 <- unique(GABAergic_0.01$Gene_Stable_ID)
length(GABAergic_0.01)
# [1] 222

# Change to data frame

GABAergic_0.01 <- as.data.frame(GABAergic_0.01)
Glutamatergic_0.01 <- as.data.frame(Glutamatergic_0.01)

##### Plot Venn Diagram #####
Venn_0.01 <- VennDiagram::venn.diagram(
  x = list(Glutamatergic_0.01$Glutamatergic_0.01, GABAergic_0.01$GABAergic_0.01), # Lists$columns of interest for intersection
  filename = NULL, # What  to save  file as
  resolution = 300, # Resolution in DPI
  category.names = c("Glutamatergic", "GABAergic"), # Names of categories
  lwd = 1, # width of each circles circumference
  lty = 1, # outline of each circle
  col =  c(Glut_Color, GABA_Color), # Color palate for outline
  cex = 1, # font size for numbers
  fontfamily = "sans", # sans font for numbers
  fontface = "bold", # numbers are bolded
  cat.cex = 1, # font size for category names
  cat.pos = 0, # position of category name, 0 is 12 o'clock, 180 is 6 o'clock
  cat.fontfamily = "sans", # sans font for category names
  cat.fontface = "bold", # bold font of category names
  fill = c(Glut_Color, GABA_Color), # Color palate for fill
  alpha = 0.40, # Transparency of each circle
  euler.d = FALSE, # Prevent scaling by size
  scaled = FALSE
)

pdf("052024_VennDiagram_Glutamatergic_GABAergic_0.01.pdf")
grid::grid.newpage()
grid::grid.draw(Venn_0.01)
dev.off()

# Shared between Glutamatergic and GABAergic cell types
Shared_0.01 <- intersect(Glutamatergic_0.01$Glutamatergic_0.01, GABAergic_0.01$GABAergic_0.01)
length(Shared_0.01)
# [1] 172

# Glutamaterigc Only
Glutamatergic_Only_0.01 <- setdiff(Glutamatergic_0.01$Glutamatergic_0.01, GABAergic_0.01$GABAergic_0.01)
length(Glutamatergic_Only_0.01)
# [1] 1222

# GABAergic Only
GABAergic_Only_0.01 <- setdiff(GABAergic_0.01$GABAergic_0.01, Glutamatergic_0.01$Glutamatergic_0.01)
length(GABAergic_Only_0.01)
# [1] 50

##### Write tables #####
write.table(Shared_0.01, "052023_Shared_Glutamatergic_GABAergic_0.01.txt", sep = "/t",
            col.names = FALSE, row.names = FALSE)

write.table(Glutamatergic_Only_0.01, "052023_Glutamatergic_Only_0.01.txt", sep = "/t",
            col.names = FALSE, row.names = FALSE)

write.table(GABAergic_Only_0.01, "052023_GABAergic_Only_0.01.txt", sep = "/t",
            col.names = FALSE, row.names = FALSE)

#### 0.10 ####
## Glutamatergic 
# Combine Glutamatergic DEG Lists
# There are 2652 unique genes in glutamatergic cell-types
Glutamatergic_0.10 <- rbind(L23_IT_CTX_0.10, L45_IT_CTX_0.10, L5_IT_CTX_0.10, 
                            L5_PT_CTX_0.10, L56_NP_CTX_0.10,
                            L6_CT_CTX_0.10, L6_IT_CTX_0.10, L6b_CTX_0.10)

Glutamatergic_0.10 <- unique(Glutamatergic_0.10$Gene_ID)
length(Glutamatergic_0.10)
# [1] 3618

# Turn this into a dataframe for intersections later
GABAergic_0.10 <- rbind(Pvalb_0.10, Sst_0.10, Vip_0.10)
GABAergic_0.10 <- unique(GABAergic_0.10$Gene_ID)
length(GABAergic_0.10)
# [1] 778

# Change to data frame

GABAergic_0.10 <- as.data.frame(GABAergic_0.10)
Glutamatergic_0.10 <- as.data.frame(Glutamatergic_0.10)

##### Plot Venn Diagram #####
Venn_0.10 <- VennDiagram::venn.diagram(
  x = list(Glutamatergic_0.10$Glutamatergic_0.10, GABAergic_0.10$GABAergic_0.10), # Lists$columns of interest for intersection
  filename = NULL, # What  to save  file as
  resolution = 300, # Resolution in DPI
  category.names = c("Glutamatergic", "GABAergic"), # Names of categories
  lwd = 1, # width of each circles circumference
  lty = 1, # outline of each circle
  col =  c(Glut_Color, GABA_Color), # Color palate for outline
  cex = 1, # font size for numbers
  fontfamily = "sans", # sans font for numbers
  fontface = "bold", # numbers are bolded
  cat.cex = 1, # font size for category names
  cat.pos = 0, # position of category name, 0 is 12 o'clock, 180 is 6 o'clock
  cat.fontfamily = "sans", # sans font for category names
  cat.fontface = "bold", # bold font of category names
  fill = c(Glut_Color, GABA_Color), # Color palate for fill
  alpha = 0.40, # Transparency of each circle
  euler.d = FALSE, # Prevent scaling by size
  scaled = FALSE
)

pdf("052024_VennDiagram_Glutamatergic_GABAergic_0.10.pdf")
grid::grid.newpage()
grid::grid.draw(Venn_0.10)
dev.off()

# Shared between Glutamatergic and GABAergic cell types
Shared_0.10 <- intersect(Glutamatergic_0.10$Glutamatergic_0.10, GABAergic_0.10$GABAergic_0.10)
length(Shared_0.10)
# [1] 635

# Glutamaterigc Only
Glutamatergic_Only_0.10 <- setdiff(Glutamatergic_0.10$Glutamatergic_0.10, GABAergic_0.10$GABAergic_0.10)
length(Glutamatergic_Only_0.10)
# [1] 2983

# GABAergic Only
GABAergic_Only_0.10 <- setdiff(GABAergic_0.10$GABAergic_0.10, Glutamatergic_0.10$Glutamatergic_0.10)
length(GABAergic_Only_0.10)
# [1] 143

##### Write tables #####
write.table(Shared_0.10, "052023_Shared_Glutamatergic_GABAergic_0.10.txt", sep = "/t",
            col.names = FALSE, row.names = FALSE)

write.table(Glutamatergic_Only_0.10, "052023_Glutamatergic_Only_0.10.txt", sep = "/t",
            col.names = FALSE, row.names = FALSE)

write.table(GABAergic_Only_0.10, "052023_GABAergic_Only_0.10.txt", sep = "/t",
            col.names = FALSE, row.names = FALSE)


