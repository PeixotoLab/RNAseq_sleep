
# Set the working directory
setwd("~/Dropbox/Sleep_RNAseq_Splicing/Single_nuclear_RNAseq_SD/KF_Unique_Functional_Annotation/DAVID_Output")

#### Read the text files that were exported from DAVID ####
# As the text file is tab delimited, we define the field separator as "\t"
# "fill = TRUE" adds blank spaces to uneven rows

DAVID_L23 <- readxl::read_excel("100623_DAVID_Output_Unique_Lists_CellTypes.xlsx", sheet = 1)
DAVID_L45 <- readxl::read_excel("100623_DAVID_Output_Unique_Lists_CellTypes.xlsx", sheet = 2)

DAVID_L23_1 <- DAVID_L23[2, 7]
DAVID_L23_2 <- DAVID_L23[3, 7]

DAVID_L23_3 <- DAVID_L23[8, 7]
DAVID_L23_4 <- DAVID_L23[9, 7]
DAVID_L23_5 <- DAVID_L23[10, 7]

DAVID_L45_1 <- DAVID_L45[3, 7]
DAVID_L45_2 <- DAVID_L45[4, 7]
DAVID_L45_3 <- DAVID_L45[5, 7]
DAVID_L45_4 <- DAVID_L45[6, 7]
DAVID_L45_5 <- DAVID_L45[7, 7]
DAVID_L45_6 <- DAVID_L45[8, 7]
DAVID_L45_7 <- DAVID_L45[9, 7]
DAVID_L45_8 <- DAVID_L45[10, 7]
DAVID_L45_9 <- DAVID_L45[11, 7]
DAVID_L45_10 <- DAVID_L45[12, 7]
DAVID_L45_11 <- DAVID_L45[13, 7]
DAVID_L45_12 <- DAVID_L45[14, 7]
DAVID_L45_13 <- DAVID_L45[15, 7]
DAVID_L45_14 <- DAVID_L45[16, 7]

DAVID_L45_15 <- DAVID_L45[20, 7]
DAVID_L45_16 <- DAVID_L45[21, 7]
DAVID_L45_17 <- DAVID_L45[22, 7]
DAVID_L45_18 <- DAVID_L45[23, 7]


z <- c(DAVID_L23_1, DAVID_L23_2, DAVID_L23_3, DAVID_L23_4, DAVID_L23_5, 
       DAVID_L45_1, DAVID_L45_2, DAVID_L45_3, DAVID_L45_4, DAVID_L45_5, 
       DAVID_L45_6, DAVID_L45_7, DAVID_L45_8, DAVID_L45_9, DAVID_L45_10, 
       DAVID_L45_11, DAVID_L45_12, DAVID_L45_13, DAVID_L45_14, DAVID_L45_15,
       DAVID_L45_16, DAVID_L45_17, DAVID_L45_18)

Highlight_Genes <- readxl::read_excel("Highlight_Genes.xlsx")
Highlight_Genes <- Highlight_Genes$Gene_Names
length(Highlight_Genes)
# [1] 46

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

 # L2/3 UP
# character(0)
# character(0)

 # L2/3 DOWN
# character(0)
# character(0)
# character(0)

 # L4/5 UP
# [1] "Cdkn1b" "Reln"   "Sgk1"  
# [1] "Cdkn1b"
# [1] "Cdkn1b"  "Gadd45a"
# character(0)
# [1] "Gadd45a"
# [1] "Cdkn1b"  "Gadd45a" "Sgk1"   
# [1] "Gadd45a"
# [1] "Cdkn1b"  "Gadd45a"
# character(0)
# [1] "Cdkn1b" "Reln"  
# [1] "Gadd45a"
# character(0)
# character(0)
# [1] "Gadd45a"

 # L4/5 DOWN
# character(0)
# character(0)
# character(0)
# character(0)
