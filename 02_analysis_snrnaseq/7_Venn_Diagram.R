# Venn Diagram intersecting Glutamatergic and GABAergic lists

#### Set working directory and load packages ####
setwd("~/Dropbox/Sleep_RNAseq_Splicing/Single_nuclear_RNAseq_SD/KF_Venn_Diagram")

library(RColorBrewer)
library(dplyr)
library(readxl)

#### Read sheets of significant lists ####
Glutamatergic <- read_excel("2023_10_04_Peixoto_Figure6_background.xlsx", sheet = "Glutamatergic")
GABAergic <- read_excel("2023_10_04_Peixoto_Figure6_background.xlsx", sheet = "GABAergic")

## Glutamatergic
# There are 2652 unique genes in glutamatergic cell-types
dim(Glutamatergic)
# [1] 2652    1

# Turn this into a dataframe for intersections later
Glutamatergic <- as.data.frame(Glutamatergic)

## GABAergic
# There are 525 unique genes in GABAergic cell-types
dim(GABAergic)
# [1] 525   1

# Turn this into a dataframe for intersections later
GABAergic <- as.data.frame(GABAergic)

#### Set colors ####
Glut_Color <- brewer.pal(9, "Blues")[6]
GABA_Color <- brewer.pal(9, "Oranges")[6]


#### Plot Venn Diagram ####
Venn <- VennDiagram::venn.diagram(
  x = list(Glutamatergic[,1], GABAergic[,1]), # Lists$columns of interest for intersection
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

grid::grid.newpage()
grid::grid.draw(Venn)

# Shared between Glutamatergic and GABAergic cell types
Shared <- intersect(Glutamatergic[,1], GABAergic[,1])
length(Shared)
# [1] 416

# Glutamaterigc Only
Glutamatergic_Only <- anti_join(Glutamatergic, GABAergic, by = "Gene_Stable_ID")
dim(Glutamatergic_Only)
# [1] 2236    1

# GABAergic Only
GABAergic_Only <- anti_join(GABAergic, Glutamatergic, by = "Gene_Stable_ID")
dim(GABAergic_Only)
# [1] 109   1

#### Write tables ####
write.table(Shared, "100423_Shared_Glutamatergic_GABAergic.txt", sep = "/t",
            col.names = FALSE, row.names = FALSE)

write.table(Glutamatergic_Only, "100423_Glutamatergic_Only.txt", sep = "/t",
            col.names = FALSE, row.names = FALSE)

write.table(GABAergic_Only, "100423_GABAergic_Only.txt", sep = "/t",
            col.names = FALSE, row.names = FALSE)

