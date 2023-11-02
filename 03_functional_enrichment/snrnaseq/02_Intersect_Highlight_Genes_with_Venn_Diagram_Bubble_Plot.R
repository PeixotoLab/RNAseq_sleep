
# Set the working directory
setwd("~/Dropbox/Sleep_RNAseq_Splicing/Single_nuclear_RNAseq_SD/KF_Venn_Diagram/DAVID_Output")

#### Read the text files that were exported from DAVID ####
# As the text file is tab delimited, we define the field separator as "\t"
# "fill = TRUE" adds blank spaces to uneven rows

DAVID_VennDiagram_Shared <- readxl::read_excel("100623_DAVID_FunctionalAnnotation_VennDiagram_Output.xlsx", sheet = 1)
DAVID_VennDiagram_Glut <- readxl::read_excel("100623_DAVID_FunctionalAnnotation_VennDiagram_Output.xlsx", sheet = 2)
DAVID_VennDiagram_GABA <- readxl::read_excel("100623_DAVID_FunctionalAnnotation_VennDiagram_Output.xlsx", sheet = 3)

## Unclustered ##
Glutamatergic_GABAergic_Up_Unclustered_1 <- DAVID_VennDiagram_Shared[2,7]
Glutamatergic_GABAergic_Up_Unclustered_2 <- DAVID_VennDiagram_Shared[3,7]
Glutamatergic_GABAergic_Up_Unclustered_3 <- DAVID_VennDiagram_Shared[4,7]
Glutamatergic_GABAergic_Up_Unclustered_4 <- DAVID_VennDiagram_Shared[5,7]
Glutamatergic_GABAergic_Up_Unclustered_5 <- DAVID_VennDiagram_Shared[6,7]

Glutamatergic_GABAergic_Down_Unclustered_1 <- DAVID_VennDiagram_Shared[10,7]
Glutamatergic_GABAergic_Down_Unclustered_2 <- DAVID_VennDiagram_Shared[11,7]
Glutamatergic_GABAergic_Down_Unclustered_3 <- DAVID_VennDiagram_Shared[12,7]
Glutamatergic_GABAergic_Down_Unclustered_4 <- DAVID_VennDiagram_Shared[13,7]

Glutamatergic_Only_Up_Unclustered_1 <- DAVID_VennDiagram_Glut[57,7]
Glutamatergic_Only_Up_Unclustered_2 <- DAVID_VennDiagram_Glut[58,7]
Glutamatergic_Only_Up_Unclustered_3 <- DAVID_VennDiagram_Glut[59,7]
Glutamatergic_Only_Up_Unclustered_4 <- DAVID_VennDiagram_Glut[60,7]
Glutamatergic_Only_Up_Unclustered_5 <- DAVID_VennDiagram_Glut[61,7]
Glutamatergic_Only_Up_Unclustered_6 <- DAVID_VennDiagram_Glut[62,7]
Glutamatergic_Only_Up_Unclustered_7 <- DAVID_VennDiagram_Glut[63,7]
Glutamatergic_Only_Up_Unclustered_8 <- DAVID_VennDiagram_Glut[64,7]
Glutamatergic_Only_Up_Unclustered_9 <- DAVID_VennDiagram_Glut[65,7]
Glutamatergic_Only_Up_Unclustered_10 <- DAVID_VennDiagram_Glut[66,7]
Glutamatergic_Only_Up_Unclustered_11 <- DAVID_VennDiagram_Glut[67,7]
Glutamatergic_Only_Up_Unclustered_12 <- DAVID_VennDiagram_Glut[68,7]
Glutamatergic_Only_Up_Unclustered_13 <- DAVID_VennDiagram_Glut[69,7]

GABAergic_Only_Up_Unclustered_1 <- DAVID_VennDiagram_GABA[2,7]

GABAergic_Only_Down_Unclustered_1 <- DAVID_VennDiagram_GABA[6,7]
GABAergic_Only_Down_Unclustered_2 <- DAVID_VennDiagram_GABA[7,7]

z <- c(Glutamatergic_GABAergic_Up_Unclustered_1, Glutamatergic_GABAergic_Up_Unclustered_2,
      Glutamatergic_GABAergic_Up_Unclustered_3, Glutamatergic_GABAergic_Up_Unclustered_4, 
      Glutamatergic_GABAergic_Up_Unclustered_5, Glutamatergic_GABAergic_Down_Unclustered_1,
      Glutamatergic_GABAergic_Down_Unclustered_2, Glutamatergic_GABAergic_Down_Unclustered_3,
      Glutamatergic_GABAergic_Down_Unclustered_4, Glutamatergic_Only_Up_Unclustered_1,
      Glutamatergic_Only_Up_Unclustered_2, Glutamatergic_Only_Up_Unclustered_3, 
      Glutamatergic_Only_Up_Unclustered_4, Glutamatergic_Only_Up_Unclustered_5, 
      Glutamatergic_Only_Up_Unclustered_6, Glutamatergic_Only_Up_Unclustered_7, 
      Glutamatergic_Only_Up_Unclustered_8, Glutamatergic_Only_Up_Unclustered_9,
      Glutamatergic_Only_Up_Unclustered_10, Glutamatergic_Only_Up_Unclustered_11,
      Glutamatergic_Only_Up_Unclustered_12, Glutamatergic_Only_Up_Unclustered_13,
      GABAergic_Only_Up_Unclustered_1, GABAergic_Only_Down_Unclustered_1,
      GABAergic_Only_Down_Unclustered_2)

Highlight_Genes <- readxl::read_excel("Peixoto_Figure_2_Supplementary_Table_2.xlsx", sheet = 1, col_names = TRUE)
Highlight_Genes <- Highlight_Genes$Gene_Names
length(Highlight_Genes)
# [1] 47

for(i in z) {
  i <- unlist(strsplit(i,","))
  i <- noquote(i)
  # Remove trailing spaces and save non-duplicated values
  i <- trimws(i)
  
  i_Genes <- intersect(i, Highlight_Genes)
  i_Cntn <- grep("^Cntn", i, value = TRUE)
  i_Ctnna <- grep("^Ctnna", i, value = TRUE)
  i_Eif <- grep("^Eif", i, value = TRUE)
  i_Fos <- grep("^Fos", i, value = TRUE)
  i_Hdac <- grep("^Hdac", i, value = TRUE)
  i_Gria <- grep("^Gria", i, value = TRUE)
  i_Grin <- grep("^Grin", i, value = TRUE)
  
  i_Gene_Intersection <- c(i_Genes, i_Cntn, i_Ctnna, i_Eif, i_Fos, i_Hdac, i_Gria, i_Grin)
  print(sort(i_Gene_Intersection[!duplicated(i_Gene_Intersection)]))
  
}

 # Glutamatergic/GABAergic Up
# [1] "Hspa5" "Hspa8"
# [1] "Per2"
# [1] "Hspa8"
# [1] "Hspa8"
# [1] "Egr1" "Per1" "Per2"

 # Glutamatergic/GABAergic Down
# character(0)
# character(0)
# character(0)
# character(0)

 # Glutamatergic Only Up
# [1] "Sgk1"
# character(0)
# [1] "Gria1"  "Grin2a"
# [1] "Bhlhe40" "Sik1"   
# character(0)
# character(0)
# character(0)
# [1] "Bhlhe40"
# character(0)
# [1] "Ptgs2"
# character(0)
# [1] "Bhlhe40"
# [1] "Bhlhe40"

 # GABAergic Up Only
# [1] "Wnt4"

 # GABAergic Down Only
# [1] "Sst" "Vip"
# character(0)

## Clustered ##
Cluster1 <- DAVID_VennDiagram_Glut[3:6, 7]
Cluster2 <- DAVID_VennDiagram_Glut[10:17, 7]
Cluster3 <- DAVID_VennDiagram_Glut[21:54, 7]
Cluster4 <- DAVID_VennDiagram_Glut[74:77, 7]
Cluster5 <- DAVID_VennDiagram_Glut[81:84, 7]

y <- c(Cluster1, Cluster2, Cluster3, Cluster4, Cluster5)

for(i in y) {
  i <- as.list(i)
  i <- noquote(i)
  i <- as.character(i)
  i <- unlist(strsplit(i,","))
  i <- noquote(i)
  # Remove trailing spaces and save non-duplicated values
  i <- trimws(i)
  
  i_Genes <- intersect(i, Highlight_Genes)
  i_Cntn <- grep("^Cntn", i, value = TRUE)
  i_Ctnna <- grep("^Ctnna", i, value = TRUE)
  i_Eif <- grep("^Eif", i, value = TRUE)
  i_Fos <- grep("^Fos", i, value = TRUE)
  i_Hdac <- grep("^Hdac", i, value = TRUE)
  i_Gria <- grep("^Gria", i, value = TRUE)
  i_Grin <- grep("^Grin", i, value = TRUE)
  
  i_Gene_Intersection <- c(i_Genes, i_Cntn, i_Ctnna, i_Eif, i_Fos, i_Hdac, i_Gria, i_Grin)
  print(sort(i_Gene_Intersection[!duplicated(i_Gene_Intersection)]))
  
}

 # Glutamaergic Up
# [1] "Gadd45b" "Reln"    "Sik1"   

 # Glutamatergic Up
# [1] "Cdkn1b"   "Cntnap5a" "Fos"      "Gadd45a"  "Gadd45b"  "Reln"     "Sgk1" 

 # Glutamatergic Up
# [1] "Camk1g"  "Cdkn1b"  "Fos"     "Fosb"    "Gadd45a" "Gadd45b" "Gria1"   "Grin2a"  "Ptgs2"   "Reln"   
# [11] "Sgk1"    "Sik1"   

 # Glutamatergic Down
# [1] "Gria4"  "Grin2d"

 # Glutamatergic Down
# [1] "Dact2"  "Mef2c"  "Nfatc3" "Wnt1"   "Wnt9a" 
