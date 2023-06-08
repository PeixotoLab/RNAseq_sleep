#### Eif Plots ####
# Goal: View the distribution of Eifs at the gene and transcript level

#### Step 1: Read Gene/Transcript Lists ####
DGE <- readxl::read_excel("060123_Peixoto_Figure1_Supplement2.xlsx", sheet = "DGE")
DTE <- readxl::read_excel("060123_Peixoto_Figure1_Supplement2.xlsx", sheet = "DTE")

#### Step 2: Subset matrices so that only Eifs remain ####
# DGE
DGE_Eifs <- DGE[grep("^Eif", DGE$Gene_Name),]
dim(DGE_Eifs)
# [1] 32  8

# DTE
DTE_Eifs <- DTE[grep("^Eif", DTE$Transcript_Name),]
dim(DTE_Eifs)
# [1] 60 11

#### Step 3: Filter so that only genes and transcripts with an abs(log2FC) > 0.2 remain ####
library(dplyr)
# DGE
DGE_Eifs_0.2 <- dplyr::filter(DGE_Eifs, abs(log2FC) > 0.2)
dim(DGE_Eifs_0.2)
# [1] 5 8

#DTE
DTE_Eifs_0.2 <- dplyr::filter(DTE_Eifs, abs(log2FC) > 0.2)
dim(DTE_Eifs_0.2)
# [1] 41 11

#### Step 4: Add an "A" or "B" based on log2FC for plotting purposes and color ####
# DGE
DGE_Eifs_0.2$Group <- ifelse(DGE_Eifs_0.2$log2FC > 0, "A", "B" ) 
DGE_Eifs_0.2$Group <- as.factor(DGE_Eifs_0.2$Group)

# DTE
DTE_Eifs_0.2$Group <- ifelse(DTE_Eifs_0.2$log2FC > 0,"A","B") 
DTE_Eifs_0.2$Group <- as.factor(DTE_Eifs_0.2$Group)

# Intersect DTE list with DGE list so that way 
# we can see which transcripts are also detected at the gene level
DTE_and_DGE <- semi_join(DTE_Eifs_0.2, DGE_Eifs_0.2, by = "Gene_Stable_ID")

#### Step 5: Bar Plot ####
suppressPackageStartupMessages(library(ggplot2)) #Version 3.4.2

# Give genes with multiple transcripts unique names
DTE_Eif_Gene_Names <- DTE_Eifs_0.2$Gene_Name
DTE_Eif_Gene_Names_Unique <- make.unique(DTE_Eif_Gene_Names)
DTE_Eifs_0.2$Gene_Name_Unique <- DTE_Eif_Gene_Names_Unique

# Add column to adjust text position
DTE_Eifs_0.2$Text_Position <- ifelse(DTE_Eifs_0.2$log2FC > 0, -0.5, .5) 
DTE_Eifs_0.2$Eif_Angles <- ifelse(DTE_Eifs_0.2$log2FC > 0, 270, 90) 

DTE_Eifs_Bar_Plot <- ggplot(DTE_Eifs_0.2, aes(x = reorder(Gene_Name_Unique, -log2FC), y = log2FC, fill = Group, color = Group)) +
  geom_col() + 
  scale_fill_manual(values = c("red1", "dodgerblue")) +
  scale_color_manual(values = c("red4", "blue2")) +
  geom_text(aes(label = Gene_Name_Unique, y = 0.2), color = "black", angle = DTE_Eifs_0.2$Eif_Angles, position=position_dodge(width = 1), size = 4, vjust= 0.5, hjust = DTE_Eifs_0.2$Text_Position) + 
  theme_light() +
  theme(axis.text.x = element_blank(), axis.title.x.bottom = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", panel.grid = element_blank())
DTE_Eifs_Bar_Plot


# DTE and DGE plot
DTE_and_DGE$Text_Position <- ifelse(DTE_and_DGE$log2FC > 0, -2, 1) 
DTE_and_DGE$Eif_Angles <- ifelse(DTE_and_DGE$log2FC > 0, 270, 90)

DTE_and_DGE_Eifs_Bar_Plot <- ggplot(DTE_and_DGE, aes(x = reorder(Gene_Name, -log2FC), y = log2FC, fill = Group, color = Group)) +
  geom_col() + 
  scale_fill_manual(values = c("red1", "dodgerblue")) +
  scale_color_manual(values = c("red4", "blue2")) +
  geom_text(aes(label = Gene_Name, y = 0.2), color = "black", angle = DTE_and_DGE$Eif_Angles, position=position_dodge(width = 1), size = 4, vjust= 0.5, hjust = DTE_and_DGE$Text_Position) + 
  theme_light() +
  theme(axis.text.x = element_blank(), axis.title.x.bottom = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", panel.grid = element_blank())
DTE_and_DGE_Eifs_Bar_Plot
