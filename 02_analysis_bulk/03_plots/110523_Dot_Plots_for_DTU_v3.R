# Dot plots for DTU and DTE analysis

# What does this plot show? How the proportion and expression of individual 
# transcripts (isoforms) change between conditions

#### Establish the working directory ####
setwd("~/Dropbox/Sleep_RNAseq_Splicing/BulkRNAseq_SD/Final_Tximeta")

#### Read Packages ####
library(readxl) # Version 1.4.2 (For importing excel file)
library(ggplot2) # Version 3.4.2 (For making plots)
library(RColorBrewer) # Version 1.1-3 (Use for color palette)
library(colorspace) # Version 2.1-0 (To modify color palette)
library(patchwork) # Version 1.1.2 (For layout) 
library(plotrix) # Version 3.8-3 (To calculate standard error)

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
# 7) I determined the mean of these median values in both HC and SD replicates (calculated again here, just to double check!)

# DTE
# 1) I took the median normalized counts of the 30 inferential replicates.
# 2) I annotated this using Perl.
# 3) I imported this annotated file into excel.
# 4) I determined the mean of these median values in both HC and SD replicates (calculated again here, just to double check!)

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
Homer1_Colors <- c(Homer1_202_Color_HC, Homer1_202_Color_SD, Homer1_203_Color_HC, 
                   Homer1_203_Color_SD, Homer1_204_Color_HC, Homer1_204_Color_SD)

# Bdnf
# Set HC Colors
Bdnf_201_Color_HC <- brewer.pal(12, "Set3")[5] 
Bdnf_205_Color_HC <- brewer.pal(12, "Set3")[6] 

# Use the "darken" function to set the SD colors
Bdnf_201_Color_SD <- darken(Bdnf_201_Color_HC, 0.2)
Bdnf_205_Color_SD <- darken(Bdnf_205_Color_HC, 0.2)

# Save the home cage and SD colors
Bdnf_Colors <- c(Bdnf_201_Color_HC, Bdnf_201_Color_SD, 
                 Bdnf_205_Color_HC, Bdnf_205_Color_SD)


##### Plots for DTU ####
# All of the transcripts shown are significant for both DTE and DTU analysis

###### Homer1 ####

# Transcripts (names identified in Shiraishi-Yamaguchi et al., 2007)
# "ENSMUST00000102752" = Homer1-204 (Homer1a)    
# "ENSMUST00000080127" = Homer1-203 (Homer1c - can't be distinguished from b from our analysis)  
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


# HC: Homecage
# SD: Sleep Deprived
Homer1_DTU_Med$Transcript_Name
# [1] "Homer1-207" "Homer1-206" "Homer1-204" "Homer1-203" "Homer1-205" "Homer1-202"

# Read the DTU proportions for the 5 replicates for both HC and SD conditions
### 202
Homer1_202_DTU <- dplyr::filter(Homer1_DTU_Med, Transcript_Name == "Homer1-202")
# HC
Homer1_202_DTU_HC <- Homer1_202_DTU[, 8:12]
# # A tibble: 1 × 5
# WTHC5_PFC_1_quant WTHC5_PFC_2_quant WTHC5_PFC_3_quant WTHC5_PFC_4_quant WTHC5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#   1             0.384             0.379             0.335             0.385             0.339

# SD
Homer1_202_DTU_SD <- Homer1_202_DTU[, 14:18]
# A tibble: 1 × 5
# WTSD5_PFC_1_quant WTSD5_PFC_2_quant WTSD5_PFC_3_quant WTSD5_PFC_4_quant WTSD5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#  1             0.244             0.249             0.251             0.255             0.254

### 203
Homer1_203_DTU <- dplyr::filter(Homer1_DTU_Med, Transcript_Name == "Homer1-203")
# HC
Homer1_203_DTU_HC <- Homer1_203_DTU[, 8:12]
# A tibble: 1 × 5
# WTHC5_PFC_1_quant WTHC5_PFC_2_quant WTHC5_PFC_3_quant WTHC5_PFC_4_quant WTHC5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#  1             0.276             0.293             0.324             0.285             0.253

# SD
Homer1_203_DTU_SD <- Homer1_203_DTU[, 14:18]
# A tibble: 1 × 5
# WTSD5_PFC_1_quant WTSD5_PFC_2_quant WTSD5_PFC_3_quant WTSD5_PFC_4_quant WTSD5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#   1             0.185             0.184             0.180             0.193             0.199

### 204
Homer1_204_DTU <- dplyr::filter(Homer1_DTU_Med, Transcript_Name == "Homer1-204")

# Read the DTU proportions for the 5 replicates for both HC and SD conditions
# HC
Homer1_204_DTU_HC <- Homer1_204_DTU[, 8:12]
# A tibble: 1 × 5
# WTHC5_PFC_1_quant WTHC5_PFC_2_quant WTHC5_PFC_3_quant WTHC5_PFC_4_quant WTHC5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#   1             0.276             0.256             0.268             0.280             0.342

# SD
Homer1_204_DTU_SD <- Homer1_204_DTU[, 14:18]
# A tibble: 1 × 5
# WTSD5_PFC_1_quant WTSD5_PFC_2_quant WTSD5_PFC_3_quant WTSD5_PFC_4_quant WTSD5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#  1             0.504             0.504             0.505             0.485             0.477

# Make a data frame for plotting, for plotting purposes, we will only include
# transcripts that have a proportion of at least 0.10
Homer1_DTU_DF <- data.frame(1:30, 1:30, 1:30, 1:30, 1:30, 1:30)

# The colnames will contain the count, condition (HC/SD) and transcript
colnames(Homer1_DTU_DF) <- c("Count", "Condition", "Transcript", "Transcript_Condition",
                             "Standard_Errors", "Mean")

# The transcript and condition (HC/SD) will be the rowname

row.names(Homer1_DTU_DF) <- c(1:30)
  
Homer1_DTU_DF$Count <- c(Homer1_202_DTU_HC, Homer1_202_DTU_SD, Homer1_203_DTU_HC, Homer1_203_DTU_SD,
                    Homer1_204_DTU_HC, Homer1_204_DTU_SD)

Homer1_DTU_DF$Condition <- rep(c("HC", "SD"), times = c(5, 5))

Homer1_DTU_DF$Transcript <- rep(c("1", "1.5", "2", "2.5", "3", "3.5"), times = c(5, 5, 5, 5, 5, 5))

Homer1_DTU_DF$Transcript_Condition <- rep(c("Homer1_202_DTU_HC", "Homer1_202_DTU_SD", "Homer1_203_DTU_HC", 
                                        "Homer1_203_DTU_SD", "Homer1_204_DTU_HC", "Homer1_204_DTU_SD"), times = c(5, 5, 5, 5, 5, 5))


# Calculate the standard error within home cage and sleep deprived animals for each significant transcript
Homer1_202_DTU_HC_SE <- std.error(as.numeric(Homer1_202_DTU_HC))
Homer1_202_DTU_SD_SE <- std.error(as.numeric(Homer1_202_DTU_SD))
Homer1_203_DTU_HC_SE <- std.error(as.numeric(Homer1_203_DTU_HC))
Homer1_203_DTU_SD_SE <- std.error(as.numeric(Homer1_203_DTU_SD))
Homer1_204_DTU_HC_SE <- std.error(as.numeric(Homer1_204_DTU_HC))
Homer1_204_DTU_SD_SE <- std.error(as.numeric(Homer1_204_DTU_SD))

Homer1_DTU_DF$Standard_Errors <- rep(c(Homer1_202_DTU_HC_SE, Homer1_202_DTU_SD_SE,
                     Homer1_203_DTU_HC_SE, Homer1_203_DTU_SD_SE,
                     Homer1_204_DTU_HC_SE, Homer1_204_DTU_SD_SE), times = c(5, 5, 5, 5, 5, 5))

# Calculate the mean for each transcript within home cage and sleep deprived animals
Homer1_202_DTU_HC_Mean <- mean(as.numeric(Homer1_202_DTU_HC))
Homer1_202_DTU_SD_Mean <- mean(as.numeric(Homer1_202_DTU_SD))
Homer1_203_DTU_HC_Mean <- mean(as.numeric(Homer1_203_DTU_HC))
Homer1_203_DTU_SD_Mean <- mean(as.numeric(Homer1_203_DTU_SD))
Homer1_204_DTU_HC_Mean <- mean(as.numeric(Homer1_204_DTU_HC))
Homer1_204_DTU_SD_Mean <- mean(as.numeric(Homer1_204_DTU_SD))

Homer1_DTU_DF$Mean <- rep(c(Homer1_202_DTU_HC_Mean, Homer1_202_DTU_SD_Mean,
                                       Homer1_203_DTU_HC_Mean, Homer1_203_DTU_SD_Mean,
                                       Homer1_204_DTU_HC_Mean, Homer1_204_DTU_SD_Mean), times = c(5, 5, 5, 5, 5, 5))


# Ensure the data.class is a data frame
data.class(Homer1_DTU_DF)
# [1] "data.frame"

Homer1_DTU_DF
# Count Condition Transcript Transcript_Condition Standard_Errors      Mean
# 1   0.384328        HC          1    Homer1_202_DTU_HC     0.011256203 0.3646076
# 2  0.3789204        HC          1    Homer1_202_DTU_HC     0.011256203 0.3646076
# 3  0.3350281        HC          1    Homer1_202_DTU_HC     0.011256203 0.3646076
# 4  0.3853528        HC          1    Homer1_202_DTU_HC     0.011256203 0.3646076
# 5  0.3394088        HC          1    Homer1_202_DTU_HC     0.011256203 0.3646076
# 6  0.2439448        SD        1.5    Homer1_202_DTU_SD     0.001952665 0.2505409
# 7  0.2492684        SD        1.5    Homer1_202_DTU_SD     0.001952665 0.2505409
# 8  0.2507163        SD        1.5    Homer1_202_DTU_SD     0.001952665 0.2505409
# 9  0.2552551        SD        1.5    Homer1_202_DTU_SD     0.001952665 0.2505409
# 10   0.25352        SD        1.5    Homer1_202_DTU_SD     0.001952665 0.2505409
# 11 0.2762195        HC          2    Homer1_203_DTU_HC     0.011524488 0.2863055
# 12  0.293149        HC          2    Homer1_203_DTU_HC     0.011524488 0.2863055
# 13 0.3236773        HC          2    Homer1_203_DTU_HC     0.011524488 0.2863055
# 14 0.2854993        HC          2    Homer1_203_DTU_HC     0.011524488 0.2863055
# 15 0.2529826        HC          2    Homer1_203_DTU_HC     0.011524488 0.2863055
# 16  0.185178        SD        2.5    Homer1_203_DTU_SD     0.003440206 0.1882096
# 17 0.1839112        SD        2.5    Homer1_203_DTU_SD     0.003440206 0.1882096
# 18 0.1800268        SD        2.5    Homer1_203_DTU_SD     0.003440206 0.1882096
# 19 0.1926679        SD        2.5    Homer1_203_DTU_SD     0.003440206 0.1882096
# 20 0.1992639        SD        2.5    Homer1_203_DTU_SD     0.003440206 0.1882096
# 21 0.2762969        HC          3    Homer1_204_DTU_HC     0.014902661 0.2846462
# 22 0.2564708        HC          3    Homer1_204_DTU_HC     0.014902661 0.2846462
# 23 0.2680429        HC          3    Homer1_204_DTU_HC     0.014902661 0.2846462
# 24 0.2804593        HC          3    Homer1_204_DTU_HC     0.014902661 0.2846462
# 25 0.3419613        HC          3    Homer1_204_DTU_HC     0.014902661 0.2846462
# 26 0.5038132        SD        3.5    Homer1_204_DTU_SD     0.005877623 0.4950942
# 27 0.5044076        SD        3.5    Homer1_204_DTU_SD     0.005877623 0.4950942
# 28  0.505281        SD        3.5    Homer1_204_DTU_SD     0.005877623 0.4950942
# 29 0.4846071        SD        3.5    Homer1_204_DTU_SD     0.005877623 0.4950942
# 30 0.4773624        SD        3.5    Homer1_204_DTU_SD     0.005877623 0.4950942

# View the structure of the data frame
# For plotting purposes, transcript and numeric will need to be counts
str(Homer1_DTU_DF)

# Dotplot

Homer1_DTU_Plot <- ggplot(Homer1_DTU_DF, aes(as.numeric(Transcript), as.numeric(Count), color= Transcript_Condition)) +
  geom_point(shape = rep(c(19,17,19,17,19,17), times = c(5,5,5,5,5,5)), size = 3) +
  labs(x = element_blank(), y = "Median Proportion", title = element_blank())  +
  scale_color_manual(values=c(Homer1_Colors)) +
  geom_errorbar(aes(ymin= as.numeric(Mean - Standard_Errors), ymax=as.numeric(Mean + Standard_Errors)), width=.2, color = "black") +
  geom_point(size = 5, stat = "summary", fun = mean, shape = "_", color = "black") +
  theme_bw() +
  scale_x_continuous(breaks = c(1.25, 2.25, 3.25)) +
  scale_y_continuous(limits = c(.15, 0.55)) +
  theme(text = element_text(size = 14), axis.text = element_text(size=10), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Homer1_DTU_Plot


###### Bdnf ####

# Transcripts (Names determined from Aid et al., 2006)
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

# HC: Homecage
# SD: Sleep Deprived

Bdnf_DTU_Med$Transcript_Name
# [1] "Bdnf-201" "Bdnf-210" "Bdnf-207" "Bdnf-209" "Bdnf-205"

# Read the DTU proportions for the 5 replicates for both HC and SD conditions
### 201
Bdnf_201_DTU <- dplyr::filter(Bdnf_DTU_Med, Transcript_Name == "Bdnf-201")
# HC
Bdnf_201_DTU_HC <- Bdnf_201_DTU[, 8:12]
# A tibble: 1 × 5
# WTHC5_PFC_1_quant WTHC5_PFC_2_quant WTHC5_PFC_3_quant WTHC5_PFC_4_quant WTHC5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#  1             0.200             0.210             0.218             0.179             0.160

# SD
Bdnf_201_DTU_SD <- Bdnf_201_DTU[, 14:18]
# A tibble: 1 × 5
# WTSD5_PFC_1_quant WTSD5_PFC_2_quant WTSD5_PFC_3_quant WTSD5_PFC_4_quant WTSD5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#  1             0.465             0.466             0.469             0.491             0.492

## 205
Bdnf_205_DTU <- dplyr::filter(Bdnf_DTU_Med, Transcript_Name == "Bdnf-205")
# HC
Bdnf_205_DTU_HC <- Bdnf_205_DTU[, 8:12]
# A tibble: 1 × 5
# WTHC5_PFC_1_quant WTHC5_PFC_2_quant WTHC5_PFC_3_quant WTHC5_PFC_4_quant WTHC5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#  1             0.441             0.308             0.350             0.356             0.371

# SD
Bdnf_205_DTU_SD <- Bdnf_205_DTU[, 14:18]
# A tibble: 1 × 5
# WTSD5_PFC_1_quant WTSD5_PFC_2_quant WTSD5_PFC_3_quant WTSD5_PFC_4_quant WTSD5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#   1             0.180             0.203             0.189             0.191             0.179

# Make a data frame for plotting
Bdnf_DTU_DF <- data.frame(1:20, 1:20, 1:20, 1:20, 1:20, 1:20)

# The colnames will contain the count, condition (HC/SD) and transcript
colnames(Bdnf_DTU_DF) <-  c("Count", "Condition", "Transcript", "Transcript_Condition",
                            "Standard_Errors", "Mean")

# The transcript and condition (HC/SD) will be the rowname
row.names(Bdnf_DTU_DF) <- c(1:20)

Bdnf_DTU_DF$Count <- c(Bdnf_201_DTU_HC, Bdnf_201_DTU_SD, Bdnf_205_DTU_HC, 
                       Bdnf_205_DTU_SD)

Bdnf_DTU_DF$Condition <- rep(c("HC", "SD"), times = c(5,5))

Bdnf_DTU_DF$Transcript <- rep(c("1", "1.5", "2", "2.5"), times = c(5, 5, 5, 5))

Bdnf_DTU_DF$Transcript_Condition <- rep(c("Bdnf_201_DTU_HC", "Bdnf_201_DTU_SD", "Bdnf_205_DTU_HC", 
                                      "Bdnf_205_DTU_SD"), times = c(5,5,5,5))

# Calculate the standard error within home cage and sleep deprived animals for each significant transcript
Bdnf_201_DTU_HC_SE <- std.error(as.numeric(Bdnf_201_DTU_HC))
Bdnf_201_DTU_SD_SE <- std.error(as.numeric(Bdnf_201_DTU_SD))
Bdnf_205_DTU_HC_SE <- std.error(as.numeric(Bdnf_205_DTU_HC))
Bdnf_205_DTU_SD_SE <- std.error(as.numeric(Bdnf_205_DTU_SD))

Bdnf_DTU_DF$Standard_Errors <- rep(c(Bdnf_201_DTU_HC_SE, Bdnf_201_DTU_SD_SE,
                                     Bdnf_205_DTU_HC_SE, Bdnf_205_DTU_SD_SE),
                                     times = c(5, 5, 5, 5))

# Calculate the mean for each transcript within home cage and sleep deprived animals
Bdnf_201_DTU_HC_Mean <- mean(as.numeric(Bdnf_201_DTU_HC))
Bdnf_201_DTU_SD_Mean <- mean(as.numeric(Bdnf_201_DTU_SD))
Bdnf_205_DTU_HC_Mean <- mean(as.numeric(Bdnf_205_DTU_HC))
Bdnf_205_DTU_SD_Mean <- mean(as.numeric(Bdnf_205_DTU_SD))


Bdnf_DTU_DF$Mean <- rep(c(Bdnf_201_DTU_HC_Mean, Bdnf_201_DTU_SD_Mean,
                            Bdnf_205_DTU_HC_Mean, Bdnf_205_DTU_SD_Mean),
                            times = c(5, 5, 5, 5))

# Ensure the data.class is a data frame
data.class(Bdnf_DTU_DF)
# [1] "data.frame"

Bdnf_DTU_DF
# Count Condition Transcript Transcript_Condition Standard_Errors      Mean
# 1   0.200326        HC          1      Bdnf_201_DTU_HC     0.010593640 0.1934274
# 2  0.2096222        HC          1      Bdnf_201_DTU_HC     0.010593640 0.1934274
# 3  0.2179973        HC          1      Bdnf_201_DTU_HC     0.010593640 0.1934274
# 4  0.1794715        HC          1      Bdnf_201_DTU_HC     0.010593640 0.1934274
# 5  0.1597199        HC          1      Bdnf_201_DTU_HC     0.010593640 0.1934274
# 6  0.4653976        SD        1.5      Bdnf_201_DTU_SD     0.006092141 0.4765869
# 7  0.4659736        SD        1.5      Bdnf_201_DTU_SD     0.006092141 0.4765869
# 8  0.4686864        SD        1.5      Bdnf_201_DTU_SD     0.006092141 0.4765869
# 9  0.4907703        SD        1.5      Bdnf_201_DTU_SD     0.006092141 0.4765869
# 10 0.4921064        SD        1.5      Bdnf_201_DTU_SD     0.006092141 0.4765869
# 11 0.4414682        HC          2      Bdnf_205_DTU_HC     0.021688295 0.3652426
# 12 0.3082558        HC          2      Bdnf_205_DTU_HC     0.021688295 0.3652426
# 13 0.3502204        HC          2      Bdnf_205_DTU_HC     0.021688295 0.3652426
# 14 0.3555433        HC          2      Bdnf_205_DTU_HC     0.021688295 0.3652426
# 15 0.3707252        HC          2      Bdnf_205_DTU_HC     0.021688295 0.3652426
# 16 0.1802555        SD        2.5      Bdnf_205_DTU_SD     0.004363939 0.1887001
# 17 0.2034005        SD        2.5      Bdnf_205_DTU_SD     0.004363939 0.1887001
# 18 0.1893209        SD        2.5      Bdnf_205_DTU_SD     0.004363939 0.1887001
# 19 0.1911541        SD        2.5      Bdnf_205_DTU_SD     0.004363939 0.1887001
# 20 0.1793695        SD        2.5      Bdnf_205_DTU_SD     0.004363939 0.1887001

# View the structure of the data frame with "str"
str(Bdnf_DTU_DF)

# Dotplot

Bdnf_DTU_Plot <- ggplot(Bdnf_DTU_DF, aes(as.numeric(Transcript), as.numeric(Count), color= Transcript_Condition)) +
  geom_point(shape = rep(c(19,17,19,17), times = c(5,5,5,5)), size = 3) +
  labs(x = element_blank(), y = "Median Proportion", title = element_blank())  +
  scale_color_manual(values=c(Bdnf_Colors)) +
  geom_errorbar(aes(ymin= as.numeric(Mean - Standard_Errors), ymax=as.numeric(Mean + Standard_Errors)), width=.13, color = "black") +
  geom_point(size = 5, stat = "summary", fun = mean, shape = "_", color = "black") +
  theme_bw() +
  scale_x_continuous(breaks = c(1.25, 2.25)) +
  scale_y_continuous(limits = c(0.15, 0.55)) +
  theme(text = element_text(size = 14), axis.text = element_text(size=10), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Bdnf_DTU_Plot


##### Plots for DTE ####
# We will repeat the above plots, but for DTE analysis

###### Homer1 ####

# Transcripts (names identified in Shiraishi-Yamaguchi et al., 2007)
# "ENSMUST00000102752" = Homer1-204 (Homer1a)    
# "ENSMUST00000080127" = Homer1-203 (Homer1c - can't be distinguished from b)  
# "ENSMUST00000079086" = Homer1-202 (Homer1d)  

# Filter so only the transcripts of Homer1 remain
Homer1_DTE_Med <- dplyr::filter(DTE_infMed, Gene_Name == "Homer1")
# # A tibble: 6 × 19
# Gene_ID          Gene_ID_Version Transcript_ID Transcript_ID_Version Gene_Description
# <chr>            <chr>           <chr>         <chr>                 <chr>           
# 1 ENSMUSG00000007… ENSMUSG0000000… ENSMUST00000… ENSMUST00000109494.7  homer scaffoldi…
# 2 ENSMUSG00000007… ENSMUSG0000000… ENSMUST00000… ENSMUST00000109493.8  homer scaffoldi…
# 3 ENSMUSG00000007… ENSMUSG0000000… ENSMUST00000… ENSMUST00000102752.9  homer scaffoldi…
# 4 ENSMUSG00000007… ENSMUSG0000000… ENSMUST00000… ENSMUST00000080127.11 homer scaffoldi…
# 5 ENSMUSG00000007… ENSMUSG0000000… ENSMUST00000… ENSMUST00000109492.8  homer scaffoldi…
# 6 ENSMUSG00000007… ENSMUSG0000000… ENSMUST00000… ENSMUST00000079086.7  homer scaffoldi…
# ℹ 14 more variables: Gene_Name <chr>, Transcript_Name <chr>,
#   WTHC5_PFC_1_quant <dbl>, WTHC5_PFC_2_quant <dbl>, WTHC5_PFC_3_quant <dbl>,
#   WTHC5_PFC_4_quant <dbl>, WTHC5_PFC_5_quant <dbl>, HC_mean <dbl>,
#   WTSD5_PFC_1_quant <dbl>, WTSD5_PFC_2_quant <dbl>, WTSD5_PFC_3_quant <dbl>,
#   WTSD5_PFC_4_quant <dbl>, WTSD5_PFC_5_quant <dbl>, SD_mean <dbl>

# I took the mean of these median values in excel (first 5 columns are control
# animals- homecage, second 5 columns are sleep deprived animals):
# HC: Homecage
# SD: Sleep Deprived
Homer1_DTE_Med$Transcript_Name
# [1] "Homer1-207" "Homer1-206" "Homer1-204" "Homer1-203" "Homer1-205" "Homer1-202"

# Read the DTE counts for the 5 replicates for both HC and SD conditions
### 202
Homer1_202_DTE <- dplyr::filter(Homer1_DTE_Med, Transcript_Name == "Homer1-202")

# HC
Homer1_202_DTE_HC <- Homer1_202_DTE[, 8:12]
# A tibble: 1 × 5
# WTHC5_PFC_1_quant WTHC5_PFC_2_quant WTHC5_PFC_3_quant WTHC5_PFC_4_quant WTHC5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#   1             2363.             2110.             1914.             2259.             2834.

# SD
Homer1_202_DTE_SD <- Homer1_202_DTE[, 14:18]
# A tibble: 1 × 5
# WTSD5_PFC_1_quant WTSD5_PFC_2_quant WTSD5_PFC_3_quant WTSD5_PFC_4_quant WTSD5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#  1             2959.             3001.             3093.             2854.             2726.

### 203
Homer1_203_DTE <- dplyr::filter(Homer1_DTE_Med, Transcript_Name == "Homer1-203")
# HC
Homer1_203_DTE_HC <- Homer1_203_DTE[, 8:12]
# A tibble: 1 × 5
# WTHC5_PFC_1_quant WTHC5_PFC_2_quant WTHC5_PFC_3_quant WTHC5_PFC_4_quant WTHC5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#   1             1710.             1655.             1868.             1677.             2128.

# SD
Homer1_203_DTE_SD <- Homer1_203_DTE[, 14:18]

### 204
Homer1_204_DTE <- dplyr::filter(Homer1_DTE_Med, Transcript_Name == "Homer1-204")
# HC
Homer1_204_DTE_HC <- Homer1_204_DTE[, 8:12]
# A tibble: 1 × 5
# WTHC5_PFC_1_quant WTHC5_PFC_2_quant WTHC5_PFC_3_quant WTHC5_PFC_4_quant WTHC5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#  1             1708.             1434.             1542.             1658.             2858.

# SD
Homer1_204_DTE_SD <- Homer1_204_DTE[, 14:18]
 # A tibble: 1 × 5
# WTSD5_PFC_1_quant WTSD5_PFC_2_quant WTSD5_PFC_3_quant WTSD5_PFC_4_quant WTSD5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#   1             6112.             6075.             6262.             5421.             5116.

# Make a data frame for plotting, for plotting purposes, we will only include
# transcripts that have a proportion of at least 0.10
Homer1_DTE_DF <- data.frame(1:30, 1:30, 1:30, 1:30, 1:30, 1:30)

# The colnames will contain the count, condition (HC/SD) and transcript
colnames(Homer1_DTE_DF) <- c("Count", "Condition", "Transcript", "Transcript_Condition",
                             "Standard_Errors", "Mean")

# The transcript and condition (HC/SD) will be the rowname

row.names(Homer1_DTE_DF) <- c(1:30)


Homer1_DTE_DF$Count <- c(Homer1_202_DTE_HC, Homer1_202_DTE_SD, Homer1_203_DTE_HC,  
                         Homer1_203_DTE_SD, Homer1_204_DTE_HC, Homer1_204_DTE_SD) 

Homer1_DTE_DF$Condition <- rep(c("HC", "SD"), times = c(5, 5))

Homer1_DTE_DF$Transcript <- rep(c("1", "1.5", "2", "2.5", "3", "3.5"), times = c(5, 5, 5, 5, 5, 5))

Homer1_DTE_DF$Transcript_Condition <- rep(c("Homer1_202_DTE_HC", "Homer1_202_DTE_SD", "Homer1_203_DTE_HC", 
                                            "Homer1_203_DTE_SD", "Homer1_204_DTE_HC", "Homer1_204_DTE_SD"), times = c(5, 5, 5, 5, 5, 5))


# Calculate the standard error within home cage and sleep deprived animals for each significant transcript
Homer1_202_DTE_HC_SE <- std.error(as.numeric(Homer1_202_DTE_HC))
Homer1_202_DTE_SD_SE <- std.error(as.numeric(Homer1_202_DTE_SD))
Homer1_203_DTE_HC_SE <- std.error(as.numeric(Homer1_203_DTE_HC))
Homer1_203_DTE_SD_SE <- std.error(as.numeric(Homer1_203_DTE_SD))
Homer1_204_DTE_HC_SE <- std.error(as.numeric(Homer1_204_DTE_HC))
Homer1_204_DTE_SD_SE <- std.error(as.numeric(Homer1_204_DTE_SD))

Homer1_DTE_DF$Standard_Errors <- rep(c(Homer1_202_DTE_HC_SE, Homer1_202_DTE_SD_SE,
                                       Homer1_203_DTE_HC_SE, Homer1_203_DTE_SD_SE,
                                       Homer1_204_DTE_HC_SE, Homer1_204_DTE_SD_SE), times = c(5, 5, 5, 5, 5, 5))

# Calculate the mean for each transcript within home cage and sleep deprived animals
Homer1_202_DTE_HC_Mean <- mean(as.numeric(Homer1_202_DTE_HC))
Homer1_202_DTE_SD_Mean <- mean(as.numeric(Homer1_202_DTE_SD))
Homer1_203_DTE_HC_Mean <- mean(as.numeric(Homer1_203_DTE_HC))
Homer1_203_DTE_SD_Mean <- mean(as.numeric(Homer1_203_DTE_SD))
Homer1_204_DTE_HC_Mean <- mean(as.numeric(Homer1_204_DTE_HC))
Homer1_204_DTE_SD_Mean <- mean(as.numeric(Homer1_204_DTE_SD))

Homer1_DTE_DF$Mean <- rep(c(Homer1_202_DTE_HC_Mean, Homer1_202_DTE_SD_Mean,
                            Homer1_203_DTE_HC_Mean, Homer1_203_DTE_SD_Mean,
                            Homer1_204_DTE_HC_Mean, Homer1_204_DTE_SD_Mean), times = c(5, 5, 5, 5, 5, 5))


# Ensure the data.class is a data frame
data.class(Homer1_DTE_DF)
# [1] "data.frame"

Homer1_DTE_DF
# Count Condition Transcript Transcript_Condition Standard_Errors     Mean
# 1   2363.35        HC          1    Homer1_202_DTE_HC       154.34616 2296.262
# 2  2110.334        HC          1    Homer1_202_DTE_HC       154.34616 2296.262
# 3  1913.682        HC          1    Homer1_202_DTE_HC       154.34616 2296.262
# 4  2259.476        HC          1    Homer1_202_DTE_HC       154.34616 2296.262
# 5  2834.466        HC          1    Homer1_202_DTE_HC       154.34616 2296.262
# 6  2958.678        SD        1.5    Homer1_202_DTE_SD        63.19912 2926.454
# 7  3001.434        SD        1.5    Homer1_202_DTE_SD        63.19912 2926.454
# 8  3092.683        SD        1.5    Homer1_202_DTE_SD        63.19912 2926.454
# 9  2853.795        SD        1.5    Homer1_202_DTE_SD        63.19912 2926.454
# 10 2725.678        SD        1.5    Homer1_202_DTE_SD        63.19912 2926.454
# 11 1710.412        HC          2    Homer1_203_DTE_HC        88.32544 1807.619
# 12  1655.02        HC          2    Homer1_203_DTE_HC        88.32544 1807.619
# 13 1867.512        HC          2    Homer1_203_DTE_HC        88.32544 1807.619
# 14 1677.108        HC          2    Homer1_203_DTE_HC        88.32544 1807.619
# 15 2128.043        HC          2    Homer1_203_DTE_HC        88.32544 1807.619
# 16 2267.954        SD        2.5    Homer1_203_DTE_SD        20.51829 2213.751
# 17 2216.407        SD        2.5    Homer1_203_DTE_SD        20.51829 2213.751
# 18 2248.941        SD        2.5    Homer1_203_DTE_SD        20.51829 2213.751
# 19  2168.05        SD        2.5    Homer1_203_DTE_SD        20.51829 2213.751
# 20 2167.404        SD        2.5    Homer1_203_DTE_SD        20.51829 2213.751
# 21 1707.663        HC          3    Homer1_204_DTE_HC       258.99109 1839.994
# 22 1434.356        HC          3    Homer1_204_DTE_HC       258.99109 1839.994
# 23 1541.912        HC          3    Homer1_204_DTE_HC       258.99109 1839.994
# 24  1657.55        HC          3    Homer1_204_DTE_HC       258.99109 1839.994
# 25 2858.491        HC          3    Homer1_204_DTE_HC       258.99109 1839.994
# 26 6112.002        SD        3.5    Homer1_204_DTE_SD       223.41958 5797.193
# 27 6074.696        SD        3.5    Homer1_204_DTE_SD       223.41958 5797.193
# 28  6262.46        SD        3.5    Homer1_204_DTE_SD       223.41958 5797.193
# 29 5420.895        SD        3.5    Homer1_204_DTE_SD       223.41958 5797.193
# 30 5115.915        SD        3.5    Homer1_204_DTE_SD       223.41958 5797.193

# View the structure of the data frame
# For plotting purposes, transcript and numeric will need to be counts
str(Homer1_DTE_DF)

# Dotplot

Homer1_DTE_Plot <- ggplot(Homer1_DTE_DF, aes(as.numeric(Transcript), as.numeric(Count), color= Transcript_Condition)) +
  geom_point(shape = rep(c(19,17,19,17,19,17), times = c(5,5,5,5,5,5)), size = 3) +
  labs(x = element_blank(), y = "Median Normalized Counts", title = element_blank())  +
  scale_color_manual(values=c(Homer1_Colors)) +
  geom_errorbar(aes(ymin= as.numeric(Mean - Standard_Errors), ymax=as.numeric(Mean + Standard_Errors)), width=.2, color = "black") +
  geom_point(size = 5, stat = "summary", fun = mean, shape = "_", color = "black") +
  theme_bw() +
  scale_x_continuous(breaks = c(1.25, 2.25, 3.25)) +
  scale_y_continuous(limits = c(1000, 6750)) +
  theme(text = element_text(size = 14), axis.text = element_text(size=10), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Homer1_DTE_Plot

###### Bdnf ####

# Transcripts (Names determined from Aid et al., 2006)
# "ENSMUST00000053317" = Bdnf-201 (Bdnf I)
# "ENSMUST00000111045" = Bdnf-205 (Bdnf VI)

# Filter so only the transcripts of Bdnf remain
Bdnf_DTE_Med <- dplyr::filter(DTE_infMed, Gene_Name == "Bdnf")
# A tibble: 6 × 19
# Gene_ID          Gene_ID_Version Transcript_ID Transcript_ID_Version Gene_Description
# <chr>            <chr>           <chr>         <chr>                 <chr>           
# 1 ENSMUSG00000048… ENSMUSG0000004… ENSMUST00000… ENSMUST00000053317.11 brain derived n…
# 2 ENSMUSG00000048… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111051.9  brain derived n…
# 3 ENSMUSG00000048… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111047.8  brain derived n…
# 4 ENSMUSG00000048… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111050.9  brain derived n…
# 5 ENSMUSG00000048… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111045.8  brain derived n…
# 6 ENSMUSG00000048… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111044.2  brain derived n…
# ℹ 14 more variables: Gene_Name <chr>, Transcript_Name <chr>,
#   WTHC5_PFC_1_quant <dbl>, WTHC5_PFC_2_quant <dbl>, WTHC5_PFC_3_quant <dbl>,
#   WTHC5_PFC_4_quant <dbl>, WTHC5_PFC_5_quant <dbl>, HC_mean <dbl>,
#   WTSD5_PFC_1_quant <dbl>, WTSD5_PFC_2_quant <dbl>, WTSD5_PFC_3_quant <dbl>,
#   WTSD5_PFC_4_quant <dbl>, WTSD5_PFC_5_quant <dbl>, SD_mean <dbl>

# HC: Homecage
# SD: Sleep Deprived

Bdnf_DTE_Med$Transcript_Name


# Read the Normalized Counts for the 5 replicates for both HC and SD conditions
### 201
Bdnf_201_DTE <- dplyr::filter(Bdnf_DTE_Med, Transcript_Name == "Bdnf-201")
# HC
Bdnf_201_DTE_HC <- Bdnf_201_DTE[, 8:12]
# A tibble: 1 × 5
# WTHC5_PFC_1_quant WTHC5_PFC_2_quant WTHC5_PFC_3_quant WTHC5_PFC_4_quant WTHC5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#   1              170.              189.              183.              175.              191.

# SD
Bdnf_201_DTE_SD <- Bdnf_201_DTE[, 14:18]
# A tibble: 1 × 5
# WTSD5_PFC_1_quant WTSD5_PFC_2_quant WTSD5_PFC_3_quant WTSD5_PFC_4_quant WTSD5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#   1              642.              612.              626.              612.              607.

## 205
Bdnf_205_DTE <- dplyr::filter(Bdnf_DTE_Med, Transcript_Name == "Bdnf-205")
# HC
Bdnf_205_DTE_HC <- Bdnf_205_DTE[, 8:12]
# # A tibble: 1 × 5
# WTHC5_PFC_1_quant WTHC5_PFC_2_quant WTHC5_PFC_3_quant WTHC5_PFC_4_quant WTHC5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#  1              387.              276.              297.              349.              455.

# SD
Bdnf_205_DTE_SD <- Bdnf_205_DTE[, 14:18]
# A tibble: 1 × 5
# WTSD5_PFC_1_quant WTSD5_PFC_2_quant WTSD5_PFC_3_quant WTSD5_PFC_4_quant WTSD5_PFC_5_quant
# <dbl>             <dbl>             <dbl>             <dbl>             <dbl>
#  1              246.              261.              259.              236.              220.

# Make a data frame for plotting
Bdnf_DTE_DF <- data.frame(1:20, 1:20, 1:20, 1:20, 1:20, 1:20)

# The colnames will contain the count, condition (HC/SD) and transcript
colnames(Bdnf_DTE_DF) <-  c("Count", "Condition", "Transcript", "Transcript_Condition",
                            "Standard_Errors", "Mean")

# The transcript and condition (HC/SD) will be the rowname
row.names(Bdnf_DTE_DF) <- c(1:20)

Bdnf_DTE_DF$Count <- c(Bdnf_201_DTE_HC, Bdnf_201_DTE_SD, Bdnf_205_DTE_HC, 
                       Bdnf_205_DTE_SD)

Bdnf_DTE_DF$Condition <- rep(c("HC", "SD"), times = c(5,5))

Bdnf_DTE_DF$Transcript <- rep(c("1", "1.5", "2", "2.5"), times = c(5, 5, 5, 5))

Bdnf_DTE_DF$Transcript_Condition <- rep(c("Bdnf_201_DTE_HC", "Bdnf_201_DTE_SD", "Bdnf_205_DTE_HC", 
                                          "Bdnf_205_DTE_SD"), times = c(5,5,5,5))

# Calculate the standard error within home cage and sleep deprived animals for each significant transcript
Bdnf_201_DTE_HC_SE <- std.error(as.numeric(Bdnf_201_DTE_HC))
Bdnf_201_DTE_SD_SE <- std.error(as.numeric(Bdnf_201_DTE_SD))
Bdnf_205_DTE_HC_SE <- std.error(as.numeric(Bdnf_205_DTE_HC))
Bdnf_205_DTE_SD_SE <- std.error(as.numeric(Bdnf_205_DTE_SD))

Bdnf_DTE_DF$Standard_Errors <- rep(c(Bdnf_201_DTE_HC_SE, Bdnf_201_DTE_SD_SE,
                                     Bdnf_205_DTE_HC_SE, Bdnf_205_DTE_SD_SE),
                                   times = c(5, 5, 5, 5))

# Calculate the mean for each transcript within home cage and sleep deprived animals
Bdnf_201_DTE_HC_Mean <- mean(as.numeric(Bdnf_201_DTE_HC))
Bdnf_201_DTE_SD_Mean <- mean(as.numeric(Bdnf_201_DTE_SD))
Bdnf_205_DTE_HC_Mean <- mean(as.numeric(Bdnf_205_DTE_HC))
Bdnf_205_DTE_SD_Mean <- mean(as.numeric(Bdnf_205_DTE_SD))


Bdnf_DTE_DF$Mean <- rep(c(Bdnf_201_DTE_HC_Mean, Bdnf_201_DTE_SD_Mean,
                          Bdnf_205_DTE_HC_Mean, Bdnf_205_DTE_SD_Mean),
                        times = c(5, 5, 5, 5))

# Ensure the data.class is a data frame
data.class(Bdnf_DTE_DF)
# [1] "data.frame"

Bdnf_DTE_DF
# Count Condition Transcript Transcript_Condition Standard_Errors     Mean
# 1  169.9225        HC          1      Bdnf_201_DTE_HC        4.007973 181.5972
# 2  188.7651        HC          1      Bdnf_201_DTE_HC        4.007973 181.5972
# 3  183.2897        HC          1      Bdnf_201_DTE_HC        4.007973 181.5972
# 4  175.0427        HC          1      Bdnf_201_DTE_HC        4.007973 181.5972
# 5   190.966        HC          1      Bdnf_201_DTE_HC        4.007973 181.5972
# 6  642.3897        SD        1.5      Bdnf_201_DTE_SD        6.484222 619.7273
# 7  611.9005        SD        1.5      Bdnf_201_DTE_SD        6.484222 619.7273
# 8  625.7589        SD        1.5      Bdnf_201_DTE_SD        6.484222 619.7273
# 9  611.7392        SD        1.5      Bdnf_201_DTE_SD        6.484222 619.7273
# 10 606.8483        SD        1.5      Bdnf_201_DTE_SD        6.484222 619.7273
# 11 386.9123        HC          2      Bdnf_205_DTE_HC       32.105975 352.8461
# 12 276.3434        HC          2      Bdnf_205_DTE_HC       32.105975 352.8461
# 13 296.8741        HC          2      Bdnf_205_DTE_HC       32.105975 352.8461
# 14 349.0044        HC          2      Bdnf_205_DTE_HC       32.105975 352.8461
# 15 455.0962        HC          2      Bdnf_205_DTE_HC       32.105975 352.8461
# 16 245.9287        SD        2.5      Bdnf_205_DTE_SD        7.665312 244.3552
# 17 261.0782        SD        2.5      Bdnf_205_DTE_SD        7.665312 244.3552
# 18 258.9428        SD        2.5      Bdnf_205_DTE_SD        7.665312 244.3552
# 19 236.2658        SD        2.5      Bdnf_205_DTE_SD        7.665312 244.3552
# 20 219.5607        SD        2.5      Bdnf_205_DTE_SD        7.665312 244.3552

# View the structure with "str"
str(Bdnf_DTE_DF)

# Dotplot

Bdnf_DTE_Plot <- ggplot(Bdnf_DTE_DF, aes(as.numeric(Transcript), as.numeric(Count), color= Transcript_Condition)) +
  geom_point(shape = rep(c(19,17,19,17), times = c(5,5,5,5)), size = 3) +
  labs(x = element_blank(), y = "Median Normalized Counts", title = element_blank())  +
  scale_color_manual(values=c(Bdnf_Colors)) +
  geom_errorbar(aes(ymin= as.numeric(Mean - Standard_Errors), ymax=as.numeric(Mean + Standard_Errors)), width=.13, color = "black") +
  geom_point(size = 5, stat = "summary", fun = mean, shape = "_", color = "black") +
  theme_bw() +
  scale_x_continuous(breaks = c(1.25, 2.25)) +
  scale_y_continuous(limits = c(100, 700)) +
  theme(text = element_text(size = 14), axis.text = element_text(size=10), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Bdnf_DTE_Plot



# Now arrange all supplemental figures in a grid using patchwork
Homer1_DTU_Plot + Homer1_DTE_Plot + Bdnf_DTU_Plot + Bdnf_DTE_Plot + 
  patchwork::plot_layout(ncol = 2)

#### Save Session Info if Needed ####

# sink('091223_BubblePlotSessionInfo.txt')
# sessionInfo()
# sink() 
  