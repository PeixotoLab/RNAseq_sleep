# Bar plots for DTU and DTE analysis

# What does this plot show? How the proportion and expression of individual 
# transcripts (isoforms) change between conditions

# Last updated: September 11, 2023

#### Establish the working directory ####
setwd("~/Dropbox/Sleep_RNAseq_Splicing/BulkRNAseq_SD/Final_Tximeta")

#### Read Packages ####
library(readxl) # Version 1.4.2 (For importing excel file)
library(ggplot2) # Version 3.4.2 (For making plots)
library(RColorBrewer) # Version 1.1-3 (Use for color palette)
library(colorspace) # Version 2.1-0 (To modify color palette)
library(patchwork) # Version 1.1.2 (For layout)

# When I was running the Fishpond DTE/DTU code, I saved the necessary text files for this analysis:

# A brief summary for the workflow of generating these files:

# DTU
# 1) Applied a secondary filter (log10mean > 1) to be extra stringent about the transcripts we are reporting.
# 2) I determined the proportion of a given transcript relative to other transcripts in the gene.
  # a) At this point, any transcripts from genes that only have one transcript are removed.
# 3) This was applied to all 30 inferential replicates.
# 4) I took the median proportions of the 30 inferential replicates.
# 5) I annotated this using Perl.
# 6) I imported this annotated file into excel.
# 7) I determined the mean of these median values in both HC and SD replicates 

# DTE
# 1) I took the median normalized counts of the 30 inferential replicates.
# 2) I annotated this using Perl.
# 3) I imported this annotated file into excel.
# 4) I determined the mean of these median values in both HC and SD replicates 

##### Read DTU and DTE InfMed Files, and Lists of Differentially Expressed Transcripts ####
# DTU
DTU_infMed <- readxl::read_excel("072023_DTU_InfMed_Annotated.xlsx", col_names = TRUE)
# Calculate the dimensions 
dim(DTU_infMed)
# [1] 43697    19

# DTE
DTE_infMed <- readxl::read_excel("071823_DTE_InfMed_Annotated.xlsx", col_names = TRUE)
# Calculate the dimensions
# There are extra filters applied to DTU analysis, which is why the dims for DTE and DTU don't match
dim(DTE_infMed)
# [1] 54030    19

#### Set Colors ####

# Homer1
# Set HC Colors
Homer1_204_Color_HC <- brewer.pal(12, "Set3")[1]
Homer1_203_Color_HC <- brewer.pal(12, "Set3")[3] 
Homer1_202_Color_HC <- brewer.pal(12, "Set3")[4] 

# Use the "darken" function to set the SD colors
Homer1_204_Color_SD <- darken(Homer1_204_Color_HC, 0.2)
Homer1_203_Color_SD <- darken(Homer1_203_Color_HC, 0.2)
Homer1_202_Color_SD <- darken(Homer1_202_Color_HC, 0.2)

# Save the home cage and SD colors
Homer1_Colors <- c(Homer1_204_Color_HC, Homer1_203_Color_HC, Homer1_202_Color_HC, 
                   Homer1_204_Color_SD, Homer1_203_Color_SD, Homer1_202_Color_SD)

# Again, apply the "darken" function to set the colors for the lines connecting
# the HC and SD points
Homer1_204_Color <- darken(Homer1_204_Color_HC, 0.1)
Homer1_203_Color <- darken(Homer1_203_Color_HC, 0.1)
Homer1_202_Color <- darken(Homer1_202_Color_HC, 0.1)

# Make a palette just for DTU colors
Homer1_Colors_DTU <- c(Homer1_202_Color, Homer1_202_Color_HC, Homer1_202_Color_SD, 
                         Homer1_203_Color, Homer1_203_Color_HC, Homer1_203_Color_SD, 
                         Homer1_204_Color, Homer1_204_Color_HC, Homer1_204_Color_SD)

# Bdnf
# Set HC Colors
Bdnf_201_Color_HC <- brewer.pal(12, "Set3")[5] 
Bdnf_205_Color_HC <- brewer.pal(12, "Set3")[6] 

# Use the "darken" function to set the SD colors
Bdnf_201_Color_SD <- darken(Bdnf_201_Color_HC, 0.2)
Bdnf_205_Color_SD <- darken(Bdnf_205_Color_HC, 0.2)

# Save the home cage and SD colors
Bdnf_Colors <- c(Bdnf_201_Color_HC, Bdnf_205_Color_HC, 
                 Bdnf_201_Color_SD, Bdnf_205_Color_SD)

# Again, apply the "darken" function to set the colors for the lines connecting
# the HC and SD points
Bdnf_201_Color <- darken(Bdnf_201_Color_HC, 0.1)
Bdnf_205_Color <- darken(Bdnf_205_Color_HC, 0.1)

# Make a palette just for DTU plots, put the colors in chronological order:
Bdnf_Colors_DTU <- c(Bdnf_201_Color, Bdnf_201_Color_HC, Bdnf_201_Color_SD, 
                     Bdnf_205_Color, Bdnf_205_Color_HC, Bdnf_205_Color_SD)


##### Plots for DTU ####
# Each plot will show the expressed transcripts (expressed: filtered)
# Reminder, a log10mean filter (> 1) was applied downstream for DTU analysis

###### Homer1 ####

# Transcripts (names identified in Shiraishi-Yamaguchi et al. 2007)
# "ENSMUST00000102752" = Homer1-204 (Homer1a)    
# "ENSMUST00000080127" = Homer1-203 (Homer1c - can't be distinguished from b)  
# "ENSMUST00000079086" = Homer1-202 (Homer1d)  

# Filter so only the transcripts of Homer1 remain
Homer1_DTU_Med <- dplyr::filter(DTU_infMed, Gene_Name == "Homer1")
# A tibble: 6 × 19
# Gene_ID    Gene_ID_Version Transcript_ID Transcript_ID_Version Gene_Description Gene_Name
# <chr>      <chr>           <chr>         <chr>                 <chr>            <chr>    
# 1 ENSMUSG00… ENSMUSG0000000… ENSMUST00000… ENSMUST00000109494.7  homer scaffoldi… Homer1   
# 2 ENSMUSG00… ENSMUSG0000000… ENSMUST00000… ENSMUST00000109493.8  homer scaffoldi… Homer1   
# 3 ENSMUSG00… ENSMUSG0000000… ENSMUST00000… ENSMUST00000102752.9  homer scaffoldi… Homer1   
# 4 ENSMUSG00… ENSMUSG0000000… ENSMUST00000… ENSMUST00000080127.11 homer scaffoldi… Homer1   
# 5 ENSMUSG00… ENSMUSG0000000… ENSMUST00000… ENSMUST00000109492.8  homer scaffoldi… Homer1   
# 6 ENSMUSG00… ENSMUSG0000000… ENSMUST00000… ENSMUST00000079086.7  homer scaffoldi… Homer1   
# ℹ 13 more variables: Transcript_Name <chr>, WTHC5_PFC_1_quant <dbl>,
#   WTHC5_PFC_2_quant <dbl>, WTHC5_PFC_3_quant <dbl>, WTHC5_PFC_4_quant <dbl>,
#   WTHC5_PFC_5_quant <dbl>, HC_mean <dbl>, WTSD5_PFC_1_quant <dbl>,
#   WTSD5_PFC_2_quant <dbl>, WTSD5_PFC_3_quant <dbl>, WTSD5_PFC_4_quant <dbl>,
#   WTSD5_PFC_5_quant <dbl>, SD_mean <dbl>

# I took the mean of these median values in excel (first 5 columns are control
# animals- homecage, second 5 columns are sleep deprived animals)
# HC: Homecage
# SD: Sleep Deprived
Homer1_DTU_Med$Transcript_Name
# [1] "Homer1-207" "Homer1-206" "Homer1-204" "Homer1-203" "Homer1-205" "Homer1-202"

Homer1_204_DTU <- dplyr::filter(Homer1_DTU_Med, Transcript_Name == "Homer1-204")
Homer1_204_DTU_HC <- Homer1_204_DTU$HC_mean
# [1] 0.2846462
Homer1_204_DTU_SD <- Homer1_204_DTU$SD_mean
# [1] 0.4950942

Homer1_203_DTU <- dplyr::filter(Homer1_DTU_Med, Transcript_Name == "Homer1-203")
Homer1_203_DTU_HC <- Homer1_203_DTU$HC_mean
# [1] 0.2863055
Homer1_203_DTU_SD <- Homer1_203_DTU$SD_mean
# [1] 0.1882096

Homer1_202_DTU <- dplyr::filter(Homer1_DTU_Med, Transcript_Name == "Homer1-202")
Homer1_202_DTU_HC <- Homer1_202_DTU$HC_mean
# [1] 0.3646076
Homer1_202_DTU_SD <- Homer1_202_DTU$SD_mean
# [1] 0.2505409

# Make a data frame for plotting, for plotting purposes, we will only include
# transcripts that have a proportion of at least 0.10
Homer1_DTU_DF <- data.frame(1:6, 1:6, 1:6)

# The colnames will contain the count, condition (HC/SD) and transcript
colnames(Homer1_DTU_DF) <- c("Count", "Condition", "Transcript")

# The transcript and condition (HC/SD) will be the rowname
row.names(Homer1_DTU_DF) <- c("Homer1_204_DTU_HC", "Homer1_203_DTU_HC", "Homer1_202_DTU_HC", 
                              "Homer1_204_DTU_SD", "Homer1_203_DTU_SD", "Homer1_202_DTU_SD")

Homer1_DTU_DF$Count <- c(Homer1_204_DTU_HC, Homer1_203_DTU_HC, Homer1_202_DTU_HC, 
                         Homer1_204_DTU_SD, Homer1_203_DTU_SD, Homer1_202_DTU_SD)

Homer1_DTU_DF$Condition <- rep(c("HC", "SD"), times = c(3,3))

Homer1_DTU_DF$Transcript <- c("Homer1-204", "Homer1-203", "Homer1-202")
 
Homer1_DTU_DF$Transcript_Condition <- c("Homer1_204_DTU_HC", "Homer1_203_DTU_HC", "Homer1_202_DTU_HC", 
                                        "Homer1_204_DTU_SD", "Homer1_203_DTU_SD", "Homer1_202_DTU_SD")
# Ensure the data.class is a data frame
data.class(Homer1_DTU_DF)
# [1] "data.frame"

## Interaction Plot
Homer1_DTU_Interaction <- ggplot(Homer1_DTU_DF) +
  aes(x = Condition, y = Count, color = Transcript_Condition, group = Transcript, 
      label = Transcript, fill = Transcript_Condition) +
  geom_text(nudge_x = 0.3, color = "black") + # If specify "Arial" family, plot will go blank during export
  geom_line(linewidth = 2, aes(group = Transcript, color = Transcript)) + # specify the linewidths (order is the numeric order of the transcripts)
  geom_point(shape = rep(c(19,17), times = c(3,3)), size = 3, 
             color = Homer1_Colors, fill = Homer1_Colors) + # geom point size is determined upon order of the dataframe
  scale_color_manual(values=c(Homer1_202_Color, Homer1_203_Color, Homer1_204_Color)) + # colors are the numeric order of the transcripts
  scale_fill_manual(values= Homer1_Colors_DTU) +
  ylab("Average Proportion") +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme_bw() +
  theme(text = element_text(size = 14), axis.text = element_text(size=10), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = NULL)

Homer1_DTU_Interaction


###### Bdnf ####

# Transcripts (Names determined from Aid et al. 2006)
# "ENSMUST00000053317" = Bdnf-201 (Bdnf I)
# "ENSMUST00000111045" = Bdnf-205 (Bdnf VI)

# Filter so only the transcripts of Bdnf remain
Bdnf_DTU_Med <- dplyr::filter(DTU_infMed, Gene_Name == "Bdnf")

# A tibble: 5 × 19
# Gene_ID    Gene_ID_Version Transcript_ID Transcript_ID_Version Gene_Description Gene_Name
# <chr>      <chr>           <chr>         <chr>                 <chr>            <chr>    
# 1 ENSMUSG00… ENSMUSG0000004… ENSMUST00000… ENSMUST00000053317.11 brain derived n… Bdnf     
# 2 ENSMUSG00… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111051.9  brain derived n… Bdnf     
# 3 ENSMUSG00… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111047.8  brain derived n… Bdnf     
# 4 ENSMUSG00… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111050.9  brain derived n… Bdnf     
# 5 ENSMUSG00… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111045.8  brain derived n… Bdnf     
# ℹ 13 more variables: Transcript_Name <chr>, WTHC5_PFC_1_quant <dbl>,
#   WTHC5_PFC_2_quant <dbl>, WTHC5_PFC_3_quant <dbl>, WTHC5_PFC_4_quant <dbl>,
#   WTHC5_PFC_5_quant <dbl>, HC_mean <dbl>, WTSD5_PFC_1_quant <dbl>,
#   WTSD5_PFC_2_quant <dbl>, WTSD5_PFC_3_quant <dbl>, WTSD5_PFC_4_quant <dbl>,
#   WTSD5_PFC_5_quant <dbl>, SD_mean <dbl>

Bdnf_DTU_Med$Transcript_Name
# [1] "Bdnf-201" "Bdnf-210" "Bdnf-207" "Bdnf-209" "Bdnf-205"

# I took the mean of these median values in excel (first 5 columns are control
# animals- homecage, second 5 columns are sleep deprived animals)
# HC: Homecage
# SD: Sleep Deprived
Bdnf_201_DTU <- dplyr::filter(Bdnf_DTU_Med, Transcript_Name == "Bdnf-201")
Bdnf_201_DTU_HC <- Bdnf_201_DTU$HC_mean
# [1] 0.1934274
Bdnf_201_DTU_SD <- Bdnf_201_DTU$SD_mean
# [1] 0.4765869

Bdnf_205_DTU <- dplyr::filter(Bdnf_DTU_Med, Transcript_Name == "Bdnf-205")
Bdnf_205_DTU_HC <- Bdnf_205_DTU$HC_mean
# [1] 0.3652426
Bdnf_205_DTU_SD <- Bdnf_205_DTU$SD_mean
# [1] 0.1887001

# Make a data frame for plotting
Bdnf_DTU_DF <- data.frame(1:4, 1:4, 1:4)

# The colnames will contain the count, condition (HC/SD) and transcript
colnames(Bdnf_DTU_DF) <- c("Count", "Condition", "Transcript")

# The transcript and condition (HC/SD) will be the rowname
row.names(Bdnf_DTU_DF) <- c("Bdnf_201_DTU_HC", "Bdnf_205_DTU_HC", 
                            "Bdnf_201_DTU_SD", "Bdnf_205_DTU_SD")

Bdnf_DTU_DF$Count <- c(Bdnf_201_DTU_HC, Bdnf_205_DTU_HC, 
                       Bdnf_201_DTU_SD, Bdnf_205_DTU_SD)

Bdnf_DTU_DF$Condition <- rep(c("HC", "SD"), times = c(2,2))

Bdnf_DTU_DF$Transcript <- c("Bdnf-201", "Bdnf-205")

Bdnf_DTU_DF$Transcript_Condition <- c("Bdnf_201_DTU_HC", "Bdnf_205_DTU_HC", 
                                      "Bdnf_201_DTU_SD", "Bdnf_205_DTU_SD")

# Ensure the data.class is a data frame
data.class(Bdnf_DTU_DF)
# [1] "data.frame"

## Interaction Plot 
Bdnf_DTU_Interaction <- ggplot(Bdnf_DTU_DF) +
  aes(x = Condition, y = Count, color = Transcript_Condition, group = Transcript,
      label = Transcript, fill = Transcript_Condition) +
  geom_text(nudge_x = 0.3, color = "black") +
  geom_line(linewidth = 2, aes(group = Transcript, color = Transcript)) + # specify the linewidths (order is the numeric order of the transcripts)
  geom_point(shape = rep(c(19,17), times = c(2,2)), size = 3, 
             color = Bdnf_Colors, fill = Bdnf_Colors) + # geom point size is determined upon order of the dataframe
  scale_color_manual(values=c(Bdnf_201_Color, Bdnf_205_Color)) + # colors are the numeric order of the transcripts
  scale_fill_manual(values= Bdnf_Colors_DTU) +
  ylab("Average Proportion") +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme_bw() +
  theme(text = element_text(size = 14), axis.text = element_text(size=10), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = NULL)

Bdnf_DTU_Interaction

##### Plots for DTE ####

# We will repeat the above plots, but for DTE analysis, we will only include
# transcripts that were included for DTU

###### Homer1 ####

# Transcripts (names identified in Shiraishi-Yamaguchi et al. 2007)
# "ENSMUST00000102752" = Homer1-204 (Homer1a)    
# "ENSMUST00000080127" = Homer1-203 (Homer1c - can't be distinguished from b)  
# "ENSMUST00000079086" = Homer1-202 (Homer1d)  

# Filter so only the transcripts of Homer1 remain
Homer1_DTE_Med <- dplyr::filter(DTE_infMed, Gene_Name == "Homer1")
# # A tibble: 6 × 19
# Gene_ID    Gene_ID_Version Transcript_ID Transcript_ID_Version Gene_Description Gene_Name
# <chr>      <chr>           <chr>         <chr>                 <chr>            <chr>    
# 1 ENSMUSG00… ENSMUSG0000000… ENSMUST00000… ENSMUST00000109494.7  homer scaffoldi… Homer1   
# 2 ENSMUSG00… ENSMUSG0000000… ENSMUST00000… ENSMUST00000109493.8  homer scaffoldi… Homer1   
# 3 ENSMUSG00… ENSMUSG0000000… ENSMUST00000… ENSMUST00000102752.9  homer scaffoldi… Homer1   
# 4 ENSMUSG00… ENSMUSG0000000… ENSMUST00000… ENSMUST00000080127.11 homer scaffoldi… Homer1   
# 5 ENSMUSG00… ENSMUSG0000000… ENSMUST00000… ENSMUST00000109492.8  homer scaffoldi… Homer1   
# 6 ENSMUSG00… ENSMUSG0000000… ENSMUST00000… ENSMUST00000079086.7  homer scaffoldi… Homer1   
# ℹ 13 more variables: Transcript_Name <chr>, WTHC5_PFC_1_quant <dbl>,
#   WTHC5_PFC_2_quant <dbl>, WTHC5_PFC_3_quant <dbl>, WTHC5_PFC_4_quant <dbl>,
#   WTHC5_PFC_5_quant <dbl>, HC_mean <dbl>, WTSD5_PFC_1_quant <dbl>,
#   WTSD5_PFC_2_quant <dbl>, WTSD5_PFC_3_quant <dbl>, WTSD5_PFC_4_quant <dbl>,
#   WTSD5_PFC_5_quant <dbl>, SD_mean <dbl>

# I took the mean of these median values in excel (first 5 columns are control
# animals- homecage, second 5 columns are sleep deprived animals):
# HC: Homecage
# SD: Sleep Deprived
Homer1_DTE_Med$Transcript_Name
# [1] "Homer1-207" "Homer1-206" "Homer1-204" "Homer1-203" "Homer1-205" "Homer1-202"

Homer1_204_DTE <- dplyr::filter(Homer1_DTE_Med, Transcript_Name == "Homer1-204")
Homer1_204_DTE_HC <- Homer1_204_DTE$HC_mean
Homer1_204_DTE_SD <- Homer1_204_DTE$SD_mean

Homer1_203_DTE <- dplyr::filter(Homer1_DTE_Med, Transcript_Name == "Homer1-203")
Homer1_203_DTE_HC <- Homer1_203_DTE$HC_mean
Homer1_203_DTE_SD <- Homer1_203_DTE$SD_mean

Homer1_202_DTE <- dplyr::filter(Homer1_DTE_Med, Transcript_Name == "Homer1-202")
Homer1_202_DTE_HC <- Homer1_202_DTE$HC_mean
Homer1_202_DTE_SD <- Homer1_202_DTE$SD_mean

# Make a data frame for plotting
Homer1_DTE_DF <- data.frame(1:6, 1:6, 1:6)

# The colnames will contain the count, condition (HC/SD) and transcript
colnames(Homer1_DTE_DF) <- c("Count", "Condition", "Transcript")

# The transcript and condition (HC/SD) will be the rowname
row.names(Homer1_DTE_DF) <- c("Homer1_204_DTE_HC", "Homer1_203_DTE_HC", "Homer1_202_DTE_HC",  
                                  "Homer1_204_DTE_SD", "Homer1_203_DTE_SD", "Homer1_202_DTE_SD")

Homer1_DTE_DF$Count <- c(Homer1_204_DTE_HC, Homer1_203_DTE_HC, Homer1_202_DTE_HC, 
                         Homer1_204_DTE_SD, Homer1_203_DTE_SD, Homer1_202_DTE_SD)

Homer1_DTE_DF$Condition <- rep(c("HC", "SD"), times = c(3,3))

Homer1_DTE_DF$Transcript <- c("Homer1-204", "Homer1-203","Homer1-202")

# Ensure the data.class is a data frame
data.class(Homer1_DTE_DF)
# [1] "data.frame"

# Bar plots for Normalized Counts
Homer1_DTE_Bar <- ggplot(Homer1_DTE_DF) +
  aes(x = Transcript, y = Count, color = Transcript, group = Condition, 
      label = Transcript, fill = Transcript) +
  geom_bar(stat = "identity", position = position_dodge(), fill = Homer1_Colors) + 
 # geom_text(color = "black", position = position_stack(vjust = .5)) + # If specify "Arial" family, plot will go blank during export
  scale_color_manual(values= rep("black", times = 12)) + 
  ylab("Normalized Counts") +
  scale_y_continuous() +
  theme_bw() +
  theme(text = element_text(size = 14), axis.text = element_text(size=10), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = NULL)

Homer1_DTE_Bar

###### Bdnf ####

# Transcripts (Names determined from Aid et al. 2006)
# "ENSMUST00000053317" = Bdnf-201 (Bdnf I)
# "ENSMUST00000111051" = Bdnf-210 (Bdnf IIC)
# "ENSMUST00000111047" = Bdnf-207 (Bdnf IIB)
# "ENSMUST00000111050" = Bdnf-209 (Bdnf IV)
# "ENSMUST00000111045" = Bdnf-205 (Bdnf VI)
# "ENSMUST00000111044" = Bdnf-204 (Bdnf VIII) (No longer included after log10mean filter)

# Filter so only the transcripts of Bdnf remain
Bdnf_DTE_Med <- dplyr::filter(DTE_infMed, Gene_Name == "Bdnf")
# A tibble: 5 × 19
# Gene_ID    Gene_ID_Version Transcript_ID Transcript_ID_Version Gene_Description Gene_Name
# <chr>      <chr>           <chr>         <chr>                 <chr>            <chr>    
# 1 ENSMUSG00… ENSMUSG0000004… ENSMUST00000… ENSMUST00000053317.11 brain derived n… Bdnf     
# 2 ENSMUSG00… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111051.9  brain derived n… Bdnf     
# 3 ENSMUSG00… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111047.8  brain derived n… Bdnf     
# 4 ENSMUSG00… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111050.9  brain derived n… Bdnf     
# 5 ENSMUSG00… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111045.8  brain derived n… Bdnf     
# ℹ 13 more variables: Transcript_Name <chr>, WTHC5_PFC_1_quant <dbl>,
#   WTHC5_PFC_2_quant <dbl>, WTHC5_PFC_3_quant <dbl>, WTHC5_PFC_4_quant <dbl>,
#   WTHC5_PFC_5_quant <dbl>, HC_mean <dbl>, WTSD5_PFC_1_quant <dbl>,
#   WTSD5_PFC_2_quant <dbl>, WTSD5_PFC_3_quant <dbl>, WTSD5_PFC_4_quant <dbl>,
#   WTSD5_PFC_5_quant <dbl>, SD_mean <dbl>

Bdnf_DTE_Med$Transcript_Name
# [1] "Bdnf-201" "Bdnf-210" "Bdnf-207" "Bdnf-209" "Bdnf-205"

# I took the mean of these median values in excel (first 5 columns are control
# animals- homecage, second 5 columns are sleep deprived animals):
# HC: Homecage
# SD: Sleep Deprived
Bdnf_201_DTE <- dplyr::filter(Bdnf_DTE_Med, Transcript_Name == "Bdnf-201")
Bdnf_201_DTE_HC <- Bdnf_201_DTE$HC_mean
Bdnf_201_DTE_SD <- Bdnf_201_DTE$SD_mean

Bdnf_205_DTE <- dplyr::filter(Bdnf_DTE_Med, Transcript_Name == "Bdnf-205")
Bdnf_205_DTE_HC <- Bdnf_205_DTE$HC_mean
Bdnf_205_DTE_SD <- Bdnf_205_DTE$SD_mean

# Make a data frame for plotting
Bdnf_DTE_DF <- data.frame(1:4, 1:4, 1:4)

# The colnames will contain the count, condition (HC/SD) and transcript
colnames(Bdnf_DTE_DF) <- c("Count", "Condition", "Transcript")

# The transcript and condition (HC/SD) will be the rowname
row.names(Bdnf_DTE_DF) <- c("Bdnf_201_DTE_HC", "Bdnf_205_DTE_HC", 
                            "Bdnf_201_DTE_SD", "Bdnf_205_DTE_SD")

Bdnf_DTE_DF$Count <- c(Bdnf_201_DTE_HC, Bdnf_205_DTE_HC, 
                           Bdnf_201_DTE_SD, Bdnf_205_DTE_SD)

Bdnf_DTE_DF$Condition <- rep(c("HC", "SD"), times = c(2,2))

Bdnf_DTE_DF$Transcript <- c("Bdnf-201", "Bdnf-205")

# Ensure the data.class is a data frame
data.class(Bdnf_DTE_DF)
# [1] "data.frame"

Bdnf_DTE_Bar <- ggplot(Bdnf_DTE_DF) +
  aes(x = Transcript, y = Count, color = Transcript, group = Condition, 
      label = Transcript, fill = Transcript) +
  geom_bar(stat = "identity", position = position_dodge(), fill = Bdnf_Colors) + 
#  geom_text(color = "black", position = position_stack(vjust = .5)) + # If specify "Arial" family, plot will go blank during export
  scale_color_manual(values= rep("black", times = 4)) + 
  ylab("Normalized Counts") +
  scale_y_continuous() +
  theme_bw() +
  theme(text = element_text(size = 14), axis.text = element_text(size=10), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = NULL)

Bdnf_DTE_Bar



# Now arrange all supplemental figures in a grid using patchwork
Homer1_DTU_Interaction + Homer1_DTE_Bar + Bdnf_DTU_Interaction + Bdnf_DTE_Bar + 
  patchwork::plot_layout(ncol = 2)

  