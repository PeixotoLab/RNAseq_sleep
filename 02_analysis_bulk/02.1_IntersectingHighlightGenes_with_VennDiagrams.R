#### Intersecting "Highlight Genes" with Venn Diagrams ####

#### Load Packages and set Working Directory ####
suppressPackageStartupMessages(library(dplyr)) # Version 1.1.1

setwd("~/Dropbox/Sleep_RNAseq_Splicing/BulkRNAseq_SD/Final_Tximeta")

#### Load Files ####
# Gene Expression
DGE <- readxl::read_excel("060123_Peixoto_Figure1_Supplement2.xlsx", sheet = "DGE")
# Transcript Expression
DTE <- readxl::read_excel("060123_Peixoto_Figure1_Supplement2.xlsx", sheet = "DTE")

# Highlight Genes 
# Read File
Highlight_Genes <- readxl::read_excel("Highlight_Genes.xlsx", sheet = 1, col_names = TRUE)
# Extract just the gene names
Gene_Names <- Highlight_Genes$Gene_Names
# Turn into a matrix
Gene_Names <- as.matrix(Gene_Names)
# Use dim to ensure that there are 39 genes present
dim(Gene_Names)
# [1] 39  1
# Use unlist function to flatten lists
Gene_Names <- unlist(strsplit(Gene_Names,","))
# Remove quotes
Gene_Names <- noquote(Gene_Names)
# Remove trailing spaces
Gene_Names <- trimws(Gene_Names)
# Ensure that there are 39 genes remaining
length(Gene_Names)
# [1] 39

#### Filter Gene and Transcript Expression Files by Log2FC ####

DGE_Up <- dplyr::filter(DGE, log2FC > 0)
DGE_Down <- dplyr::filter(DGE, log2FC < 0)

DTE_Up <- dplyr::filter(DTE, log2FC > 0)
DTE_Down <- dplyr::filter(DTE, log2FC < 0)


#### Get Intersections, Up ####
Both_DGE_DTE_Up <- semi_join(DTE_Up, DGE_Up, by = "Gene_Stable_ID")
length(Both_DGE_DTE_Up$Gene_Stable_ID[!duplicated(Both_DGE_DTE_Up$Gene_Stable_ID)])
# [1] 3269

Just_DTE_Up <- anti_join(DTE_Up, DGE_Up, by = "Gene_Stable_ID")
length(Just_DTE_Up$Gene_Stable_ID[!duplicated(Just_DTE_Up$Gene_Stable_ID)])
# [1] 1381

Just_DGE_Up <- anti_join(DGE_Up, DTE_Up, by = "Gene_Stable_ID")
length(Just_DGE_Up$Gene_Stable_ID[!duplicated(Just_DGE_Up$Gene_Stable_ID)])
# [1] 872

#### Get Intersections, Down ####
Both_DGE_DTE_Down <- semi_join(DTE_Down, DGE_Down, by = "Gene_Stable_ID")
length(Both_DGE_DTE_Down$Gene_Stable_ID[!duplicated(Both_DGE_DTE_Down$Gene_Stable_ID)])
# [1] 3836

Just_DTE_Down <- anti_join(DTE_Down, DGE_Down, by = "Gene_Stable_ID")
length(Just_DTE_Down$Gene_Stable_ID[!duplicated(Just_DTE_Down$Gene_Stable_ID)])
# [1] 3117

Just_DGE_Down <- anti_join(DGE_Down, DTE_Down, by = "Gene_Stable_ID")
length(Just_DGE_Down$Gene_Stable_ID[!duplicated(Just_DGE_Down$Gene_Stable_ID)])
# [1] 528

#### Take Intersection Lists above and Intersect with Highlight Genes, Up ####

# DTE, DGE Both Up
Both_DGE_DTE_Up_Highlight_Genes <- intersect(Both_DGE_DTE_Up$Gene_Name, Gene_Names)
# Determine how many genes show up in the intersection
length(Both_DGE_DTE_Up_Highlight_Genes[!duplicated(Both_DGE_DTE_Up_Highlight_Genes)])
# [1] 18

# View list of intersecting genes
sort(Both_DGE_DTE_Up_Highlight_Genes[!duplicated(Both_DGE_DTE_Up_Highlight_Genes)])
# [1] "Arc"     "Bdnf"    "Bhlhe40" "Camk1g"  "Crh"     "Egr1"    "Fos"     "Gadd45b" "Homer1" 
# [10] "Hspa5"   "Nr4a1"   "P4ha1"   "Per1"    "Per2"    "Ptgs2"   "Sgk1"    "Sik1"    "Sult1a1"

Both_DGE_DTE_Up_Eifs <- Both_DGE_DTE_Up[grep("^Eif", Both_DGE_DTE_Up$Gene_Name),]
Both_DGE_DTE_Up_Eifs$Gene_Name[!duplicated(Both_DGE_DTE_Up_Eifs$Gene_Name)]
# [1] "Eif4ebp2"  "Eif4g2"    "Eif1b"     "Eif3a"     "Eif3j2"    "Eif4e2"    "Eif1ax"    "Eif2s2"   
# [9] "Eif2ak1"   "Eif4a1"    "Eif1a"     "Eif4a2"    "Eif4enif1"
Both_DGE_DTE_Up_Ctnn <- Both_DGE_DTE_Up[grep("^Ctnna", Both_DGE_DTE_Up$Gene_Name),]
Both_DGE_DTE_Up_Ctnn$Gene_Name[!duplicated(Both_DGE_DTE_Up_Ctnn$Gene_Name)]
# [1] "Ctnna1"
Both_DGE_DTE_Up_Cntn <- Both_DGE_DTE_Up[grep("^Cntn", Both_DGE_DTE_Up$Gene_Name),]
Both_DGE_DTE_Up_Cntn$Gene_Name[!duplicated(Both_DGE_DTE_Up_Cntn$Gene_Name)]
# character(0)
Both_DGE_DTE_Up_Hdac <- Both_DGE_DTE_Up[grep("^Hdac", Both_DGE_DTE_Up$Gene_Name),]
Both_DGE_DTE_Up_Hdac$Gene_Name[!duplicated(Both_DGE_DTE_Up_Hdac$Gene_Name)]
# [1] "Hdac8" "Hdac2"

# DTE, Up
Just_DTE_Up_Highlight_Genes <- intersect(Just_DTE_Up$Gene_Name, Gene_Names)
# Determine how many genes show up in the intersection
length(Just_DTE_Up_Highlight_Genes[!duplicated(Just_DTE_Up_Highlight_Genes)])
# [1] 1

# View list of intersecting genes
sort(Just_DTE_Up_Highlight_Genes[!duplicated(Just_DTE_Up_Highlight_Genes)])
# [1] Ube3a"

Just_DTE_Up_Eifs <- Just_DTE_Up[grep("^Eif", Just_DTE_Up$Gene_Name),]
Just_DTE_Up_Eifs$Gene_Name[!duplicated(Just_DTE_Up_Eifs$Gene_Name)]
# [1] "Eif2d"   "Eif2s1"  "Eif2s3y" "Eif3b"   "Eif2a"   "Eif5"    "Eif2b4"  "Eif5a2"  "Eif1" 
Just_DTE_Up_Ctnn <- Just_DTE_Up[grep("^Ctnna", Just_DTE_Up$Gene_Name),]
Just_DTE_Up_Ctnn$Gene_Name[!duplicated(Just_DTE_Up_Ctnn$Gene_Name)]
# character(0)
Just_DTE_Up_Cntn <- Just_DTE_Up[grep("^Cntn", Just_DTE_Up$Gene_Name),]
Just_DTE_Up_Cntn$Gene_Name[!duplicated(Just_DTE_Up_Cntn$Gene_Name)]
# [1] "Cntnap5c" "Cntn4" 
Just_DTE_Up_Hdac <- Just_DTE_Up[grep("^Hdac", Just_DTE_Up$Gene_Name),]
Just_DTE_Up_Hdac$Gene_Name[!duplicated(Just_DTE_Up_Hdac$Gene_Name)]
# [1] "Hdac3" "Hdac4"

# DGE, Up
Just_DGE_Up_Highlight_Genes <- intersect(Just_DGE_Up$Gene_Name, Gene_Names)
# Determine how many genes show up in the intersection
length(Just_DGE_Up_Highlight_Genes[!duplicated(Just_DGE_Up_Highlight_Genes)])
# [1] 0

# View list of intersecting genes
sort(Just_DGE_Up_Highlight_Genes[!duplicated(Just_DGE_Up_Highlight_Genes)])
# character(0)

Just_DGE_Up_Eifs <- Just_DGE_Up[grep("^Eif", Just_DGE_Up$Gene_Name),]
Just_DGE_Up_Eifs$Gene_Name[!duplicated(Just_DGE_Up_Eifs$Gene_Name)]
# [1] "Eif5b"  "Eif3j1" "Eif3c"  "Eif4g1" "Eif4e"  "Eif4a3" "Eif4e3"

Just_DGE_Up_Ctnn <- Just_DGE_Up[grep("^Ctnna", Just_DGE_Up$Gene_Name),]
Just_DGE_Up_Ctnn$Gene_Name[!duplicated(Just_DGE_Up_Ctnn$Gene_Name)]
# character(0)

Just_DGE_Up_Cntn <- Just_DGE_Up[grep("^Cntn", Just_DGE_Up$Gene_Name),]
Just_DGE_Up_Cntn$Gene_Name[!duplicated(Just_DGE_Up_Cntn$Gene_Name)]
# character(0)

Just_DGE_Up_Hdac <- Just_DGE_Up[grep("^Hdac", Just_DGE_Up$Gene_Name),]
Just_DGE_Up_Hdac$Gene_Name[!duplicated(Just_DGE_Up_Hdac$Gene_Name)]
# character(0)

#### Take Intersection Lists above and Intersect with Add File 1, Down ####

# DTE, DGE Both Down
Both_DGE_DTE_Down_Highlight_Genes <- intersect(Both_DGE_DTE_Down$Gene_Name, Gene_Names)
# Determine how many genes show up in the intersection
length(Both_DGE_DTE_Down_Highlight_Genes[!duplicated(Both_DGE_DTE_Down_Highlight_Genes)])
# [1] 12

# View list of intersecting genes
sort(Both_DGE_DTE_Down_Highlight_Genes[!duplicated(Both_DGE_DTE_Down_Highlight_Genes)])
# [1] "Cirbp"  "Dact2"  "Dbp"    "Mef2c"  "Nfatc3" "Opalin" "Sst"    "Tfrc"   "Tipin"  "Usp43" 
# [11] "Usp46"  "Wnt9a" 

Both_DGE_DTE_Down_Eifs <- Both_DGE_DTE_Down[grep("^Eif", Both_DGE_DTE_Down$Gene_Name),]
Both_DGE_DTE_Down_Eifs$Gene_Name[!duplicated(Both_DGE_DTE_Down_Eifs$Gene_Name)]
# [1] "Eif3h"    "Eif3d"    "Eif4b"    "Eif5a"    "Eif4g3"   "Eif6"     "Eif2b3"   "Eif2ak2" 
# [9] "Eif2d"    "Eif4ebp1" "Eif2ak4"  "Eif2b4"

Both_DGE_DTE_Down_Ctnn <- Both_DGE_DTE_Down[grep("^Ctnna", Both_DGE_DTE_Down$Gene_Name),]
Both_DGE_DTE_Down_Ctnn$Gene_Name[!duplicated(Both_DGE_DTE_Down_Ctnn$Gene_Name)]
# [1] "Ctnnal1"

Both_DGE_DTE_Down_Cntn <- Both_DGE_DTE_Down[grep("^Cntn", Both_DGE_DTE_Down$Gene_Name),]
Both_DGE_DTE_Down_Cntn$Gene_Name[!duplicated(Both_DGE_DTE_Down_Cntn$Gene_Name)]
# [1] "Cntn3" "Cntn6"

Both_DGE_DTE_Down_Hdac <- Both_DGE_DTE_Down[grep("^Hdac", Both_DGE_DTE_Down$Gene_Name),]
Both_DGE_DTE_Down_Hdac$Gene_Name[!duplicated(Both_DGE_DTE_Down_Hdac$Gene_Name)]
# [1] "Hdac11" "Hdac1"  "Hdac7"  "Hdac3"  "Hdac9"  "Hdac6"  "Hdac5" 

# DTE, Down
Just_DTE_Down_Highlight_Genes <- intersect(Just_DTE_Down$Gene_Name, Gene_Names)
# Determine how many genes show up in the intersection
length(Just_DTE_Down_Highlight_Genes[!duplicated(Just_DTE_Down_Highlight_Genes)])
# [1] 2

# View list of intersecting genes
sort(Just_DTE_Down_Highlight_Genes[!duplicated(Just_DTE_Down_Highlight_Genes)])
# [1] "Bdnf" "Nr4a1" 

Just_DTE_Down_Eifs <- Just_DTE_Down[grep("^Eif", Just_DTE_Down$Gene_Name),]
Just_DTE_Down_Eifs$Gene_Name[!duplicated(Just_DTE_Down_Eifs$Gene_Name)]
# [1] "Eif3f"     "Eif2s3y"   "Eif4h"     "Eif4e2"    "Eif4a1"    "Eif2s2"    "Eif2b1"   
# [8] "Eif3m"     "Eif3i"     "Eif2ak3"   "Eif4g1"    "Eif4a2"    "Eif4enif1" "Eif3a"  

Just_DTE_Down_Ctnn <- Just_DTE_Down[grep("^Ctnna", Just_DTE_Down$Gene_Name),]
Just_DTE_Down_Ctnn$Gene_Name[!duplicated(Just_DTE_Down_Ctnn$Gene_Name)]
# character(0)

Just_DTE_Down_Cntn <- Just_DTE_Down[grep("^Cntn", Just_DTE_Down$Gene_Name),]
Just_DTE_Down_Cntn$Gene_Name[!duplicated(Just_DTE_Down_Cntn$Gene_Name)]
# [1] "Cntn5" "Cntn4" "Cntn1"

Just_DTE_Down_Hdac <- Just_DTE_Down[grep("^Hdac", Just_DTE_Down$Gene_Name),]
Just_DTE_Down_Hdac$Gene_Name[!duplicated(Just_DTE_Down_Hdac$Gene_Name)]
# character(0)

# DGE, Down
Just_DGE_Down_Highlight_Genes <- intersect(Just_DGE_Down$Gene_Name, Gene_Names)
# Determine how many genes show up in the intersection
length(Just_DGE_Down_Highlight_Genes[!duplicated(Just_DGE_Down_Highlight_Genes)])
# [1] 0

# View list of intersecting genes
sort(Just_DGE_Down_Highlight_Genes[!duplicated(Just_DGE_Down_Highlight_Genes)])

Just_DGE_Down_Eifs <- Just_DGE_Down[grep("^Eif", Just_DGE_Down$Gene_Name),]
Just_DGE_Down_Eifs$Gene_Name[!duplicated(Just_DGE_Down_Eifs$Gene_Name)]
# character(0)

Just_DGE_Down_Ctnn <- Just_DGE_Down[grep("^Ctnna", Just_DGE_Down$Gene_Name),]
Just_DGE_Down_Ctnn$Gene_Name[!duplicated(Just_DGE_Down_Ctnn$Gene_Name)]
# character(0)

Just_DGE_Down_Cntn <- Just_DGE_Down[grep("^Cntn", Just_DGE_Down$Gene_Name),]
Just_DGE_Down_Cntn$Gene_Name[!duplicated(Just_DGE_Down_Cntn$Gene_Name)]
# character(0)

Just_DGE_Down_Hdac <- Just_DGE_Down[grep("^Hdac", Just_DGE_Down$Gene_Name),]
Just_DGE_Down_Hdac$Gene_Name[!duplicated(Just_DGE_Down_Hdac$Gene_Name)]
# [1] "Hdac10"

#### Save Session Info if Needed ####
# sink('VennDiagram_Intersections_DTE_DGE.txt')
# sessionInfo()
# sink() 
