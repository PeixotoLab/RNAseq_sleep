#### Intersect Highlight Genes w/ DAVID Output ####

#### Load packages and set the working directory ####
library(readxl) # Version 1.4.2

setwd("~/Dropbox/Sleep_RNAseq_Splicing/BulkRNAseq_SD/Final_Tximeta")

#### Read the list of genes that we would like highlighted ####
# Note: We have changed the list of genes to highlight to the ones below.

# Read the file containing the subset of highlight genes that we are interested in including
# on the bubble plot:
# Below we will also select to include all Rbms and Eifs
Highlight_Genes <- c("Hdac3", "Ep400", "Chd9", "Camk1", "Prkaa1", "Dnmt1", "Gar1")
Highlight_Genes <- noquote(Highlight_Genes)
length(Highlight_Genes)
# [1] 7

# Read file containing the output of DAVID
DTU_Clusters <- readxl::read_xlsx("~/Dropbox/Sleep_RNAseq_Splicing/BulkRNAseq_SD/Final_Tximeta/080223_DAVID_DTU/080223_Supplementary_Table_4.xlsx")
head(DTU_Clusters)
# A tibble: 25 × 14
# `Annotation Cluster 1`   Enrichment Score: 3.020…¹ ...3  ...4  ...5  ...6  ...7  ...8  ...9 
# <chr>                    <chr>                     <chr> <chr> <chr> <chr> <chr> <chr> <chr>
# 1 Category                 Term                      Count %     PVal… Genes Gene… List… Pop …
# 2 UP_KW_BIOLOGICAL_PROCESS KW-0507~mRNA processing   56    3.60… 3.06… ENSM… Prpf… 804   342  
# 3 UP_KW_MOLECULAR_FUNCTION KW-0694~RNA-binding       89    5.73… 1.04… ENSM… Ybx3… 838   629  
# 4 KEGG_PATHWAY             mmu03040:Spliceosome      24    1.54… 1.45… ENSM… Rbmx… 589   128  
# 5 UP_KW_BIOLOGICAL_PROCESS KW-0508~mRNA splicing     44    2.83… 1.78… ENSM… Rbmx… 804   271  
# 6 NA                       NA                        NA    NA    NA    NA    NA    NA    NA   
# ℹ abbreviated name: ¹​`Enrichment Score: 3.0209437533365633`
# ℹ 5 more variables: ...10 <chr>, ...11 <chr>, ...12 <chr>, ...13 <chr>, ...14 <chr>
DTU_Clusters[1,]
# A tibble: 1 × 14
# `Annotation Cluster 1` Enrichment Score: 3.02094…¹ ...3  ...4  ...5  ...6  ...7  ...8  ...9 
# <chr>                  <chr>                       <chr> <chr> <chr> <chr> <chr> <chr> <chr>
# 1 Category               Term                        Count %     PVal… Genes Gene… List… Pop …
# ℹ abbreviated name: ¹​`Enrichment Score: 3.0209437533365633`
# ℹ 5 more variables: ...10 <chr>, ...11 <chr>, ...12 <chr>, ...13 <chr>, ...14 <chr>

DTU_Clusters_Header <- DTU_Clusters[1,]

# Subset so only the rows that contain the terms remain:
DTU_Clusters_Terms <- DTU_Clusters[c(2:5, 9:11, 15:17, 20:25), ]
colnames(DTU_Clusters_Terms) <- DTU_Clusters_Header

dim(DTU_Clusters_Terms)
# [1] 16 14

# Separate out each row for intersections:
mRNAprocessing <- DTU_Clusters_Terms[1,]
RNAbinding <- DTU_Clusters_Terms[2,]
Spliceosome <- DTU_Clusters_Terms[3,]
mRNAsplicing <- DTU_Clusters_Terms[4,]
ChromatinRegulator <- DTU_Clusters_Terms[5,]
Transcription <- DTU_Clusters_Terms[6,]
TranscriptionRegulation <- DTU_Clusters_Terms[7,]
Kinase <- DTU_Clusters_Terms[8,]
SerineThreonineProteinKinase <- DTU_Clusters_Terms[9,]
Transferase <- DTU_Clusters_Terms[10,]
RibosomeBiogenesis <- DTU_Clusters_Terms[11,]
CellCycle <- DTU_Clusters_Terms[12,]
GuanineNucleotideReleasingFactor <- DTU_Clusters_Terms[13,]
CellCycle2 <- DTU_Clusters_Terms[14,]
mRNAsurveillancePathway <- DTU_Clusters_Terms[15,]
DNArecombination <- DTU_Clusters_Terms[16,]

z <- c(mRNAprocessing$GeneNames, RNAbinding$GeneNames, Spliceosome$GeneNames,
       mRNAsplicing$GeneNames, ChromatinRegulator$GeneNames, Transcription$GeneNames,
       TranscriptionRegulation$GeneNames, Kinase$GeneNames, SerineThreonineProteinKinase$GeneNames,
       Transferase$GeneNames, RibosomeBiogenesis$GeneNames, CellCycle$GeneNames,
       GuanineNucleotideReleasingFactor$GeneNames, CellCycle2$GeneNames, 
       mRNAsurveillancePathway$GeneNames, DNArecombination$GeneNames)

for(i in z) {
  i <- unlist(strsplit(i,","))
  i <- noquote(i)
  # Remove trailing spaces and save non-duplicated values
  i <- trimws(i)
  
  i_Genes <- intersect(i, Highlight_Genes)
  i_Eifs <- grep("^Eif", i, value = TRUE)
  i_Rbms <- grep("^Rbm", i, value = TRUE)
  
  i_Gene_Intersection <- c(i_Genes, i_Eifs, i_Rbms)
  print(sort(i_Gene_Intersection[!duplicated(i_Gene_Intersection)]))

}
# [1] "Rbm4" "Rbmx"
# [1] "Eif2ak4" "Gar1"    "Rbm27"   "Rbm3"    "Rbm34"   "Rbm4"    "Rbm7"    "Rbmx"   
# [1] "Rbmx"
# [1] "Rbm4" "Rbmx"
# [1] "Chd9"   "Dnmt1"  "Ep400"  "Hdac3"  "Prkaa1"
# [1] "Chd9"   "Dnmt1"  "Hdac3"  "Prkaa1" "Rbmx"  
# [1] "Chd9"   "Dnmt1"  "Hdac3"  "Prkaa1"
# [1] "Camk1"   "Eif2ak4" "Prkaa1" 
# [1] "Camk1"   "Eif2ak4" "Prkaa1" 
# [1] "Camk1"   "Dnmt1"   "Eif2ak4" "Prkaa1" 
# [1] "Gar1"
# character(0)
# character(0)
# [1] "Camk1"   "Eif2ak4" "Prkaa1" 
# character(0)
# character(0)

# With the following code you can verify the loop above

# Both mRNA processing and RNA binding are terms- here we separate and intersect the gene names

# mRNA processing
mRNAprocessingGenes <- noquote(mRNAprocessing$GeneNames)
mRNAprocessingGenes <- unlist(strsplit(mRNAprocessingGenes,","))
mRNAprocessingGenes <- noquote(mRNAprocessingGenes)
# Remove trailing spaces and save non-duplicated values
mRNAprocessingGenes <- trimws(mRNAprocessingGenes)

mRNAprocessing_Gene_Names <- intersect(mRNAprocessingGenes, Highlight_Genes)
mRNAprocessing_Eifs <- grep("^Eif", mRNAprocessingGenes, value = TRUE)
mRNAprocessing_Rbm <- grep("^Rbm", mRNAprocessingGenes, value = TRUE)

mRNAprocessing_Gene_Intersection <- c(mRNAprocessing_Gene_Names, mRNAprocessing_Eifs, mRNAprocessing_Rbm)
sort(mRNAprocessing_Gene_Intersection[!duplicated(mRNAprocessing_Gene_Intersection)])
# [1] "Rbm4" "Rbmx"

# RNAbinding
RNAbindingGenes <- noquote(RNAbinding$GeneNames)
RNAbindingGenes <- unlist(strsplit(RNAbindingGenes,","))
RNAbindingGenes <- noquote(RNAbindingGenes)
# Remove trailing spaces and save non-duplicated values
RNAbindingGenes <- trimws(RNAbindingGenes)

RNAbinding_Gene_Names <- intersect(RNAbindingGenes, Highlight_Genes)
RNAbinding_Eifs <- grep("^Eif", RNAbindingGenes, value = TRUE)
RNAbinding_Rbm <- grep("^Rbm", RNAbindingGenes, value = TRUE)

RNAbinding_Gene_Intersection <- c(RNAbinding_Gene_Names, RNAbinding_Eifs, RNAbinding_Rbm)
sort(RNAbinding_Gene_Intersection[!duplicated(RNAbinding_Gene_Intersection)])
# [1] "Eif2ak4" "Gar1"    "Rbm27"   "Rbm3"    "Rbm34"   "Rbm4"    "Rbm7"    "Rbmx"  


