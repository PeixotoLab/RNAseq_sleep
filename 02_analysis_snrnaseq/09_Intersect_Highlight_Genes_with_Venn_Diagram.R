#### Intersecting Highlight Genes with DTE and DGE Venn Diagrams ####

# Summary: I will intersect transcript and gene expression lists so that we can 
# determine the following:

## What genes can be detected in both Glut/GABA cell-types?

## Which genes can only be detected in Glutamatergic or GABAergic cell-types?

# R Version 4.2.2

#### Load packages and set the working directory ####
library("readxl") # 1.4.2
library("dplyr") # 1.1.1

setwd("~/Dropbox/Sleep_RNAseq_Splicing/Single_nuclear_RNAseq_SD/KF_Venn_Diagram")

#### READ LISTS ####

Shared_Glutamatergic_GABAergic <- readxl::read_excel("100423_Venn_Diagram_Glutamatergic_GABAergic.xlsx", sheet = "Shared_Glutamatergic_GABAergic")
dim(Shared_Glutamatergic_GABAergic)
# [1] 416   4

Glutamatergic_Only <- readxl::read_excel("100423_Venn_Diagram_Glutamatergic_GABAergic.xlsx", sheet = "Glutamatergic_Only")
dim(Glutamatergic_Only)
# [1] 2236    4

GABAergic_Only <- readxl::read_excel("100423_Venn_Diagram_Glutamatergic_GABAergic.xlsx", sheet = "GABAergic_Only")
dim(GABAergic_Only)
# [1] 109   4

#### Now intersect Highlight Gene Lists ####
Highlight_Genes <- readxl::read_excel("Highlight_Genes.xlsx", sheet = 1, col_names = TRUE)
Gene_Names <- Highlight_Genes$Gene_Names
Gene_Names <- as.matrix(Gene_Names)
dim(Gene_Names)
# [1] 46  1
Gene_Names <- unlist(strsplit(Gene_Names,","))
Gene_Names <- noquote(Gene_Names)
# Remove trailing spaces and save non-duplicated values
Gene_Names <- trimws(Gene_Names)
length(Gene_Names)
# [1] 46

# Glut/GABA
Shared_Glutamatergic_GABAergic_Gene_Names <- intersect(Shared_Glutamatergic_GABAergic$Gene_Name, Gene_Names)
Shared_Glutamatergic_GABAergic_Eifs <- Shared_Glutamatergic_GABAergic[grep("^Eif", Shared_Glutamatergic_GABAergic$Gene_Name),]
Shared_Glutamatergic_GABAergic_Ctnn <- Shared_Glutamatergic_GABAergic[grep("^Ctnna", Shared_Glutamatergic_GABAergic$Gene_Name),]
Shared_Glutamatergic_GABAergic_Cntn <- Shared_Glutamatergic_GABAergic[grep("^Cntn", Shared_Glutamatergic_GABAergic$Gene_Name),]
Shared_Glutamatergic_GABAergic_Hdac <- Shared_Glutamatergic_GABAergic[grep("^Hdac", Shared_Glutamatergic_GABAergic$Gene_Name),]
Shared_Glutamatergic_GABAergic_Fos <- Shared_Glutamatergic_GABAergic[grep("^Fos", Shared_Glutamatergic_GABAergic$Gene_Name),]
Shared_Glutamatergic_GABAergic_Grin <- Shared_Glutamatergic_GABAergic[grep("^Grin", Shared_Glutamatergic_GABAergic$Gene_Name),]
Shared_Glutamatergic_GABAergic_Gria <- Shared_Glutamatergic_GABAergic[grep("^Gria", Shared_Glutamatergic_GABAergic$Gene_Name),]

Glut_GABA_Gene_Names <- c(Shared_Glutamatergic_GABAergic_Gene_Names, Shared_Glutamatergic_GABAergic_Eifs$Gene_Name, Shared_Glutamatergic_GABAergic_Ctnn$Gene_Name, Shared_Glutamatergic_GABAergic_Cntn$Gene_Name, Shared_Glutamatergic_GABAergic_Hdac$Gene_Name, Shared_Glutamatergic_GABAergic_Fos$Gene_Name, Shared_Glutamatergic_GABAergic_Grin$Gene_Name, Shared_Glutamatergic_GABAergic_Gria$Gene_Name)
sort(Glut_GABA_Gene_Names[!duplicated(Glut_GABA_Gene_Names)])
# [1] "Arc"     "Bdnf"    "Cirbp"   "Cntnap3" "Egr1"    "Fosl2"   "Homer1"  "Hspa5"   "Hspa8"   "Nr4a1"  
# [11] "P4ha1"   "Per1"    "Per2"    "Ube3a"   "Usp50"      

# Glut Only
Glutamatergic_Only_Gene_Names <- intersect(Glutamatergic_Only$Gene_Name, Gene_Names)
Glutamatergic_Only_Eifs <- Glutamatergic_Only[grep("^Eif", Glutamatergic_Only$Gene_Name),]
Glutamatergic_Only_Ctnn <- Glutamatergic_Only[grep("^Ctnna", Glutamatergic_Only$Gene_Name),]
Glutamatergic_Only_Cntn <- Glutamatergic_Only[grep("^Cntn", Glutamatergic_Only$Gene_Name),]
Glutamatergic_Only_Hdac <- Glutamatergic_Only[grep("^Hdac", Glutamatergic_Only$Gene_Name),]
Glutamatergic_Only_Fos <- Glutamatergic_Only[grep("^Fos", Glutamatergic_Only$Gene_Name),]
Glutamatergic_Only_Grin <- Glutamatergic_Only[grep("^Grin", Glutamatergic_Only$Gene_Name),]
Glutamatergic_Only_Gria <- Glutamatergic_Only[grep("^Gria", Glutamatergic_Only$Gene_Name),]

Glut_Only_Gene_Names <- c(Glutamatergic_Only_Gene_Names, Glutamatergic_Only_Eifs$Gene_Name, Glutamatergic_Only_Ctnn$Gene_Name, Glutamatergic_Only_Cntn$Gene_Name, Glutamatergic_Only_Hdac$Gene_Name, Glutamatergic_Only_Fos$Gene_Name, Glutamatergic_Only_Grin$Gene_Name, Glutamatergic_Only_Gria$Gene_Name)
sort(Glut_Only_Gene_Names[!duplicated(Glut_Only_Gene_Names)])
# [1] "Bhlhe40"  "Camk1g"   "Cdkn1b"   "Cntn5"    "Cntn6"    "Cntnap4"  "Cntnap5a" "Cntnap5b" "Ctnna3"  
# [10] "Dact2"    "Eif4ebp1" "Fos"      "Fosb"     "Gadd45a"  "Gadd45b"  "Gria1"    "Gria4"    "Grin1os" 
# [19] "Grin2a"   "Grin2d"   "Mef2c"    "Nfatc3"   "Ptgs2"    "Reln"     "Sgk1"     "Sik1"     "Sult1a1" 
# [28] "Tipin"    "Usp43"    "Usp46"    "Wnt1"     "Wnt9a"   

# GABA Only
GABAergic_Only_Gene_Names <- intersect(GABAergic_Only$Gene_Name, Gene_Names)
GABAergic_Only_Eifs <- GABAergic_Only[grep("^Eif", GABAergic_Only$Gene_Name),]
GABAergic_Only_Ctnn <- GABAergic_Only[grep("^Ctnna", GABAergic_Only$Gene_Name),]
GABAergic_Only_Cntn <- GABAergic_Only[grep("^Cntn", GABAergic_Only$Gene_Name),]
GABAergic_Only_Hdac <- GABAergic_Only[grep("^Hdac", GABAergic_Only$Gene_Name),]
GABAergic_Only_Fos <- GABAergic_Only[grep("^Fos", GABAergic_Only$Gene_Name),]
GABAergic_Only_Grin <- GABAergic_Only[grep("^Grin", GABAergic_Only$Gene_Name),]
GABAergic_Only_Gria <- GABAergic_Only[grep("^Gria", GABAergic_Only$Gene_Name),]

GABA_Only_Gene_Names <- c(GABAergic_Only_Gene_Names, GABAergic_Only_Eifs$Gene_Name, GABAergic_Only_Ctnn$Gene_Name, GABAergic_Only_Cntn$Gene_Name, GABAergic_Only_Hdac$Gene_Name, GABAergic_Only_Fos$Gene_Name, GABAergic_Only_Grin$Gene_Name, GABAergic_Only_Gria$Gene_Name)
sort(GABA_Only_Gene_Names[!duplicated(GABA_Only_Gene_Names)])
# [1] "Eif2s3y" "Npas1"   "Sst"     "Vip"     "Wnt4"   
