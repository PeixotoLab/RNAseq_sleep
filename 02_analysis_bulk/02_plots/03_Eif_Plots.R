#### Eif Plots ####
# Goal: View the distribution of Eifs at the gene and transcript level

#### Load packages and set the working directory ####
library(readxl) # Version 1.4.2 (For importing excel file)
library(ggplot2) # Version 3.4.2 (For making plots)

setwd("~/Dropbox/Sleep_RNAseq_Splicing/BulkRNAseq_SD/Final_Tximeta")

#### Read Gene/Transcript Lists ####
DGE <- readxl::read_excel("060123_Peixoto_Figure1_Supplement2.xlsx", sheet = "DGE")
DTE <- readxl::read_excel("060123_Peixoto_Figure1_Supplement2.xlsx", sheet = "DTE")

#### Subset matrices so that only Eifs remain ####
# DGE
DGE_Eifs <- DGE[grep("^Eif", DGE$Gene_Name),]
dim(DGE_Eifs)
# [1] 32  8

# DTE
DTE_Eifs <- DTE[grep("^Eif", DTE$Transcript_Name),]
dim(DTE_Eifs)
# [1] 60 11

#### Filter so that only genes and transcripts with an abs(log2FC) > 0.2 remain ####

# DGE
DGE_Eifs_0.2 <- dplyr::filter(DGE_Eifs, abs(log2FC) > 0.2)
dim(DGE_Eifs_0.2)
# [1] 5 8
DGE_Eifs_0.2$Gene_Name
# [1] "Eif3j2"   "Eif4ebp2" "Eif5b"    "Eif2ak2"  "Eif4ebp1"

# DTE
DTE_Eifs_0.2 <- dplyr::filter(DTE_Eifs, abs(log2FC) > 0.2)
dim(DTE_Eifs_0.2)
# [1] 41 11
DTE_Eifs_0.2$Gene_Name
# [1] "Eif2d"     "Eif2s2"    "Eif2s3y"   "Eif2b4"    "Eif3j2"    "Eif2a"     "Eif4ebp2" 
# [8] "Eif4e2"    "Eif4e2"    "Eif2s3y"   "Eif4g2"    "Eif5"      "Eif1a"     "Eif2b3"   
# [15] "Eif2s3y"   "Eif4h"     "Eif2ak2"   "Eif5a"     "Eif2d"     "Eif4e2"    "Eif4g3"   
# [22] "Eif2s3y"   "Eif4a1"    "Eif4ebp1"  "Eif2s2"    "Eif5a"     "Eif2b1"    "Eif3m"    
# [29] "Eif3i"     "Eif5a"     "Eif2ak3"   "Eif4g1"    "Eif4b"     "Eif4g3"    "Eif4a2"   
# [36] "Eif2ak4"   "Eif5a"     "Eif4a1"    "Eif2b4"    "Eif4enif1" "Eif3a" 

#### Add an "A" (up) or "B" (down) based on log2FC for plotting purposes and color ####
# DGE
DGE_Eifs_0.2$Group <- ifelse(DGE_Eifs_0.2$log2FC > 0, "A", "B") 
DGE_Eifs_0.2$Group <- as.factor(DGE_Eifs_0.2$Group)

# DTE
DTE_Eifs_0.2$Group <- ifelse(DTE_Eifs_0.2$log2FC > 0, "A", "B") 
DTE_Eifs_0.2$Group <- as.factor(DTE_Eifs_0.2$Group)

# Intersect DTE list with DGE list so that way 
# we can see which transcripts are also detected at the gene level to highlight in the figure later on
DTE_and_DGE <- semi_join(DTE_Eifs_0.2, DGE_Eifs_0.2, by = "Gene_Stable_ID")
DTE_and_DGE$Gene_Name
# [1] "Eif3j2"   "Eif4ebp2" "Eif2ak2"  "Eif4ebp1"

#### Bar Plot ####
# Add column to adjust text position
DTE_Eifs_0.2$Eif_Angles <- ifelse(DTE_Eifs_0.2$log2FC > 0, 270, 90) 

DTE_Eifs_Bar_Plot <- ggplot(DTE_Eifs_0.2, aes(x = reorder(Transcript_Name, -log2FC), y = log2FC, fill = Group, color = Group)) +
  geom_col() + 
  scale_fill_manual(values = c("red1", "dodgerblue")) +
  scale_color_manual(values = c("red4", "blue2")) +
  geom_text(aes(label = Transcript_Name, y = ifelse(log2FC > 0, -0.05, 0.05)), color = "black", fontface = "italic",
             angle = DTE_Eifs_0.2$Eif_Angles, position=position_dodge(width = 1), size = 7/.pt, vjust= 0.5, hjust = "bottom") + 
  theme_light() +
  theme(axis.text.x = element_blank(), axis.title.x.bottom = element_blank(), 
        axis.ticks.x = element_blank(), legend.position = "none", 
        panel.grid = element_blank(), axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7))

DTE_Eifs_Bar_Plot

#### Save Session Info if Needed ####

# sink('08_2023_Eif_Plot_Session_Info.txt')
# sessionInfo()
# sink() 

