# Using ggplot to make interaction plots 
# What does this plot show? How the counts of individual transcripts (isoforms) 
# change between conditions.
# June 1, 2023

# Establish the working directory
setwd("~/Dropbox/Sleep_RNAseq_Splicing/BulkRNAseq_SD/Final_Tximeta")

#### Read Packages ####
library(readxl) # Version 1.4.2 (For importing excel file)
library(ggplot2) # Version 3.4.2 (For making plots)
library(RColorBrewer) # Version 1.1-3 (For colors)
library(patchwork) # Version 1.1.2 (For layout)

# To obtain the scaled counts, read tables with the mean and median scaled counts
# Note that I will need to convert to matrices from data frames

# I copied the transcript IDs from the Fishpond code after viewing the symbols

##### Read InfMed File ####
# What does this file contain? It contains the median proportions of a transcript
# for a particular replicate (Median of all of the proportions of the inferential replicates)
# Each column is an animal
infMed <- readxl::read_excel("051723_DTU_InfMed_Annotated.xlsx", col_names = TRUE)
dim(infMed)
# [1] 47241    19

#### Plots ####
# Each plot will show the expressed transcripts for a given gene (such as Bdnf, Homer1),
# as long as that gene has more than one transcript
# Rather than plotting all of the transcripts for all of the replicates on one figure,
# I will take the median inferential replicates (which I exported in the Fishpond code)
# and take the average across the home cage and sleep deprived animals.

#### Homer1 ####
# Subset transcripts of Homer1
Homer1_Med <- dplyr::filter(infMed, Gene_Name == "Homer1")

# A tibble: 6 × 17
# Gene_ID           Gene_ID_Version Transcript_ID Transcript_ID_Version Gene_Description Gene_Name
# <chr>             <chr>           <chr>         <chr>                 <chr>            <chr>    
# 1 ENSMUSG000000076… ENSMUSG0000000… ENSMUST00000… ENSMUST00000109494.7  homer scaffoldi… Homer1   
# 2 ENSMUSG000000076… ENSMUSG0000000… ENSMUST00000… ENSMUST00000109493.8  homer scaffoldi… Homer1   
# 3 ENSMUSG000000076… ENSMUSG0000000… ENSMUST00000… ENSMUST00000102752.9  homer scaffoldi… Homer1   
# 4 ENSMUSG000000076… ENSMUSG0000000… ENSMUST00000… ENSMUST00000080127.11 homer scaffoldi… Homer1   
# 5 ENSMUSG000000076… ENSMUSG0000000… ENSMUST00000… ENSMUST00000109492.8  homer scaffoldi… Homer1   
# 6 ENSMUSG000000076… ENSMUSG0000000… ENSMUST00000… ENSMUST00000079086.7  homer scaffoldi… Homer1   
# ℹ 11 more variables: Transcript_Name <chr>, WTHC5_PFC_1_quant <dbl>, WTHC5_PFC_2_quant <dbl>,
#   WTHC5_PFC_3_quant <dbl>, WTHC5_PFC_4_quant <dbl>, WTHC5_PFC_5_quant <dbl>,
#   WTSD5_PFC_1_quant <dbl>, WTSD5_PFC_2_quant <dbl>, WTSD5_PFC_3_quant <dbl>,
#   WTSD5_PFC_4_quant <dbl>, WTSD5_PFC_5_quant <dbl>

# We will now take the mean of these median values (columns 1-5 are control
# animals- homecage, columns 6-10 are sleep deprived animals)
# HC: Homecage
# SD: Sleep Deprived
Homer1_Med$Transcript_Name
# [1] "Homer1-207" "Homer1-206" "Homer1-204" "Homer1-203" "Homer1-205" "Homer1-202"

Homer1_207 <- dplyr::filter(Homer1_Med, Transcript_Name == "Homer1-207")
Homer1_207_HC <- Homer1_207$HC_mean
Homer1_207_SD <- Homer1_207$SD_mean

Homer1_206 <- dplyr::filter(Homer1_Med, Transcript_Name == "Homer1-206")
Homer1_206_HC <- Homer1_206$HC_mean
Homer1_206_SD <- Homer1_206$SD_mean

Homer1_204 <- dplyr::filter(Homer1_Med, Transcript_Name == "Homer1-204")
Homer1_204_HC <- Homer1_204$HC_mean
Homer1_204_SD <- Homer1_204$SD_mean

Homer1_203 <- dplyr::filter(Homer1_Med, Transcript_Name == "Homer1-203")
Homer1_203_HC <- Homer1_203$HC_mean
Homer1_203_SD <- Homer1_203$SD_mean

Homer1_205 <- dplyr::filter(Homer1_Med, Transcript_Name == "Homer1-205")
Homer1_205_HC <- Homer1_205$HC_mean
Homer1_205_SD <- Homer1_205$SD_mean

Homer1_202 <- dplyr::filter(Homer1_Med, Transcript_Name == "Homer1-202")
Homer1_202_HC <- Homer1_202$HC_mean
Homer1_202_SD <- Homer1_202$SD_mean

Homer1_Matrix <- data.frame(1:12, 1:12, 1:12)
colnames(Homer1_Matrix) <- c("Count", "Condition", "Transcript")
row.names(Homer1_Matrix) <- c("Homer1_207_HC", "Homer1_206_HC", "Homer1_204_HC", "Homer1_203_HC", 
                              "Homer1_205_HC", "Homer1_202_HC", "Homer1_207_SD", "Homer1_206_SD", 
                              "Homer1_204_SD", "Homer1_203_SD", "Homer1_205_SD", "Homer1_202_SD")

Homer1_Matrix$Count <- c(Homer1_207_HC, Homer1_206_HC, Homer1_204_HC, Homer1_203_HC, Homer1_205_HC, 
                         Homer1_202_HC, Homer1_207_SD, Homer1_206_SD, Homer1_204_SD, Homer1_203_SD,
                         Homer1_205_SD, Homer1_202_SD)

Homer1_Matrix$Condition <- rep(c("HC", "SD"), times = c(6,6))

Homer1_Matrix$Transcript <- c("Homer1-207", "Homer1-206", "Homer1-204", "Homer1-203",
                              "Homer1-205", "Homer1-202")
 
Homer1_DF <- as.data.frame(Homer1_Matrix)

# Set Colors
Homer1_207_Color <- brewer.pal(9, "Greys")[4]
Homer1_206_Color <- brewer.pal(9, "Greys")[4] 
Homer1_204_Color <- brewer.pal(8, "Set2")[1] 
Homer1_203_Color <- brewer.pal(8, "Set2")[2] 
Homer1_205_Color <- brewer.pal(9, "Greys")[4] 
Homer1_202_Color <- brewer.pal(8, "Set2")[3] 

Homer1_Interaction <- ggplot(Homer1_DF) +
  aes(x = Condition, y = Count, color = Transcript, group = Transcript, 
      label = Transcript) +
  geom_text(nudge_x = 0.3, color = "black") + # If specify "Arial" family, plot will go blank during export FYI
  geom_line(linewidth = rep(c(2,1), times = c(6,6)), aes(group = Transcript, color = Transcript)) +
  scale_color_manual(values=c(Homer1_202_Color, Homer1_203_Color, Homer1_204_Color, 
                              Homer1_205_Color, Homer1_206_Color, Homer1_207_Color)) +
  geom_point(shape = rep(c(19,17), times = c(6,6)), size = rep(c(1,2,1,2,1,2,1,2), times = c(2,2,1,1,2,2,1,1)), 
             aes(color = Transcript)) +
  ylab("Average Proportion") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14), axis.text = element_text(size=14), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = NULL)

Homer1_Interaction

#### Bdnf ####

# Transcripts
# "ENSMUST00000053317" = Bdnf-201 (Bdnf I)
# "ENSMUST00000111051" = Bdnf-210 (Bdnf IIC)
# "ENSMUST00000111047" = Bdnf-207 (Bdnf IIB)
# "ENSMUST00000111050" = Bdnf-209 (Bdnf IV)
# "ENSMUST00000111045" = Bdnf-205 (Bdnf VI)
# "ENSMUST00000111044" = Bdnf-204 (Bdnf VIII)

# Subset transcripts of Bdnf
Bdnf_Med <- dplyr::filter(infMed, Gene_Name == "Bdnf")
# A tibble: 6 × 19
# Gene_ID           Gene_ID_Version Transcript_ID Transcript_ID_Version Gene_Description Gene_Name
# <chr>             <chr>           <chr>         <chr>                 <chr>            <chr>    
# 1 ENSMUSG000000484… ENSMUSG0000004… ENSMUST00000… ENSMUST00000053317.11 brain derived n… Bdnf     
# 2 ENSMUSG000000484… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111051.9  brain derived n… Bdnf     
# 3 ENSMUSG000000484… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111047.8  brain derived n… Bdnf     
# 4 ENSMUSG000000484… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111050.9  brain derived n… Bdnf     
# 5 ENSMUSG000000484… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111045.8  brain derived n… Bdnf     
# 6 ENSMUSG000000484… ENSMUSG0000004… ENSMUST00000… ENSMUST00000111044.2  brain derived n… Bdnf     
# ℹ 13 more variables: Transcript_Name <chr>, WTHC5_PFC_1_quant <dbl>, WTHC5_PFC_2_quant <dbl>,
#   WTHC5_PFC_3_quant <dbl>, WTHC5_PFC_4_quant <dbl>, WTHC5_PFC_5_quant <dbl>, HC_mean <dbl>,
#   WTSD5_PFC_1_quant <dbl>, WTSD5_PFC_2_quant <dbl>, WTSD5_PFC_3_quant <dbl>,
#   WTSD5_PFC_4_quant <dbl>, WTSD5_PFC_5_quant <dbl>, SD_mean <dbl>

Bdnf_Med$Transcript_Name
# [1] "Bdnf-201" "Bdnf-210" "Bdnf-207" "Bdnf-209" "Bdnf-205" "Bdnf-204"

# We will now take the mean of these median values (first 5 columns are control
# animals-- homecage, second 5 columns are sleep deprived animals)
# HC: Homecage
# SD: Sleep Deprived
Bdnf_201 <- dplyr::filter(Bdnf_Med, Transcript_Name == "Bdnf-201")
Bdnf_201_HC <- Bdnf_201$HC_mean
Bdnf_201_SD <- Bdnf_201$SD_mean

Bdnf_210 <- dplyr::filter(Bdnf_Med, Transcript_Name == "Bdnf-210")
Bdnf_210_HC <- Bdnf_210$HC_mean
Bdnf_210_SD <- Bdnf_210$SD_mean

Bdnf_207 <- dplyr::filter(Bdnf_Med, Transcript_Name == "Bdnf-207")
Bdnf_207_HC <- Bdnf_207$HC_mean
Bdnf_207_SD <- Bdnf_207$SD_mean

Bdnf_209 <- dplyr::filter(Bdnf_Med, Transcript_Name == "Bdnf-209")
Bdnf_209_HC <- Bdnf_209$HC_mean
Bdnf_209_SD <- Bdnf_209$SD_mean

Bdnf_205 <- dplyr::filter(Bdnf_Med, Transcript_Name == "Bdnf-205")
Bdnf_205_HC <- Bdnf_205$HC_mean
Bdnf_205_SD <- Bdnf_205$SD_mean

Bdnf_204 <- dplyr::filter(Bdnf_Med, Transcript_Name == "Bdnf-204")
Bdnf_204_HC <- Bdnf_204$HC_mean
Bdnf_204_SD <- Bdnf_204$SD_mean


Bdnf_Matrix <- data.frame(1:12, 1:12, 1:12)
colnames(Bdnf_Matrix) <- c("Count", "Condition", "Transcript")
row.names(Bdnf_Matrix) <- c("Bdnf_201_HC", "Bdnf_210_HC", "Bdnf_207_HC", "Bdnf_209_HC", 
                            "Bdnf_205_HC", "Bdnf_204_HC", "Bdnf_201_SD", "Bdnf_210_SD", "Bdnf_207_SD", "Bdnf_209_SD", 
                            "Bdnf_205_SD", "Bdnf_204_SD")

Bdnf_Matrix$Count <- c(Bdnf_201_HC, Bdnf_210_HC, Bdnf_207_HC, Bdnf_209_HC, Bdnf_205_HC, 
                       Bdnf_204_HC, Bdnf_201_SD, Bdnf_210_SD, Bdnf_207_SD, Bdnf_209_SD, 
                       Bdnf_205_SD, Bdnf_204_SD)

Bdnf_Matrix$Condition <- rep(c("HC", "SD"), times = c(6,6))

Bdnf_Matrix$Transcript <- c("Bdnf-201", "Bdnf-210", "Bdnf-207", "Bdnf-209", "Bdnf-205", "Bdnf-204")

Bdnf_DF <- as.data.frame(Bdnf_Matrix)

# Set Colors
Bdnf_201_Color <- brewer.pal(8, "Set2")[4] 
Bdnf_210_Color <- brewer.pal(9, "Greys")[4] 
Bdnf_207_Color <- brewer.pal(9, "Greys")[4] 
Bdnf_209_Color <- brewer.pal(9, "Greys")[4] 
Bdnf_205_Color  <- brewer.pal(8, "Set2")[5] 
Bdnf_204_Color  <- brewer.pal(9, "Greys")[4] 

Bdnf_Interaction <- ggplot(Bdnf_DF) +
  aes(x = Condition, y = Count, color = Transcript, group = Transcript,
      label = Transcript) +
  geom_text(nudge_x = 0.3, color = "black") +
  geom_line(linewidth = rep(c(2,1,2,1), times = c(2,2,2,6)), aes(group = Transcript, color = Transcript)) +
  scale_color_manual(values=c(Bdnf_201_Color, Bdnf_204_Color, Bdnf_205_Color, 
                              Bdnf_207_Color, Bdnf_209_Color, Bdnf_210_Color)) +
  geom_point(shape = rep(c(19,17), times = c(6,6)), size = rep(c(2,1,2,1,2,1,2,1), times = c(1,3,1,1,1,3,1,1)), 
             aes(color = Transcript)) +
  ylab(NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14), axis.text = element_text(size=14), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  labs(x = NULL)


Bdnf_Interaction



#### Cirbp ####

# Transcripts
# "ENSMUST00000154726" = Cirbp-207
# "ENSMUST00000105365" = Cirbp-202
# "ENSMUST00000054666" = Cirbp-201
# "ENSMUST00000140774" = Cirbp-204
# "ENSMUST00000151321" = Cirbp-206
# "ENSMUST00000147190" = Cirbp-205

# Subset transcripts of Cirbp
Cirbp_Med <- dplyr::filter(infMed, Gene_Name == "Cirbp")
# # A tibble: 6 × 19
# Gene_ID           Gene_ID_Version Transcript_ID Transcript_ID_Version Gene_Description Gene_Name
# <chr>             <chr>           <chr>         <chr>                 <chr>            <chr>    
# 1 ENSMUSG000000451… ENSMUSG0000004… ENSMUST00000… ENSMUST00000154726.7  cold inducible … Cirbp    
# 2 ENSMUSG000000451… ENSMUSG0000004… ENSMUST00000… ENSMUST00000105365.8  cold inducible … Cirbp    
# 3 ENSMUSG000000451… ENSMUSG0000004… ENSMUST00000… ENSMUST00000054666.6  cold inducible … Cirbp    
# 4 ENSMUSG000000451… ENSMUSG0000004… ENSMUST00000… ENSMUST00000140774.1  cold inducible … Cirbp    
# 5 ENSMUSG000000451… ENSMUSG0000004… ENSMUST00000… ENSMUST00000151321.1  cold inducible … Cirbp    
# 6 ENSMUSG000000451… ENSMUSG0000004… ENSMUST00000… ENSMUST00000147190.1  cold inducible … Cirbp    
# ℹ 13 more variables: Transcript_Name <chr>, WTHC5_PFC_1_quant <dbl>, WTHC5_PFC_2_quant <dbl>,
#   WTHC5_PFC_3_quant <dbl>, WTHC5_PFC_4_quant <dbl>, WTHC5_PFC_5_quant <dbl>, HC_mean <dbl>,
#   WTSD5_PFC_1_quant <dbl>, WTSD5_PFC_2_quant <dbl>, WTSD5_PFC_3_quant <dbl>,
#   WTSD5_PFC_4_quant <dbl>, WTSD5_PFC_5_quant <dbl>, SD_mean <dbl>

Cirbp_Med$Transcript_Name
# [1] "Cirbp-207" "Cirbp-202" "Cirbp-201" "Cirbp-204" "Cirbp-206" "Cirbp-205"

# We will now take the mean of these median values (columns 1-5 are control
# animals- homecage, columns 6-10 are sleep deprived animals)
# HC: Homecage
# SD: Sleep Deprived
Cirbp_207 <- dplyr::filter(Cirbp_Med, Transcript_Name == "Cirbp-207")
Cirbp_207_HC <- Cirbp_207$HC_mean
Cirbp_207_SD <- Cirbp_207$SD_mean

Cirbp_202 <- dplyr::filter(Cirbp_Med, Transcript_Name == "Cirbp-202")
Cirbp_202_HC <- Cirbp_202$HC_mean
Cirbp_202_SD <- Cirbp_202$SD_mean

Cirbp_201 <- dplyr::filter(Cirbp_Med, Transcript_Name == "Cirbp-201")
Cirbp_201_HC <- Cirbp_201$HC_mean
Cirbp_201_SD <- Cirbp_201$SD_mean

Cirbp_204 <- dplyr::filter(Cirbp_Med, Transcript_Name == "Cirbp-204")
Cirbp_204_HC <- Cirbp_204$HC_mean
Cirbp_204_SD <- Cirbp_204$SD_mean

Cirbp_206 <- dplyr::filter(Cirbp_Med, Transcript_Name == "Cirbp-206")
Cirbp_206_HC <- Cirbp_206$HC_mean
Cirbp_206_SD <- Cirbp_206$SD_mean

Cirbp_205 <- dplyr::filter(Cirbp_Med, Transcript_Name == "Cirbp-205")
Cirbp_205_HC <- Cirbp_205$HC_mean
Cirbp_205_SD <- Cirbp_205$SD_mean


Cirbp_Matrix <- data.frame(1:12, 1:12, 1:12)
colnames(Cirbp_Matrix) <- c("Count", "Condition", "Transcript")
row.names(Cirbp_Matrix) <- c("Cirbp_207_HC", "Cirbp_202_HC", "Cirbp_201_HC", 
                             "Cirbp_204_HC", "Cirbp_206_HC", "Cirbp_205_HC", 
                             "Cirbp_207_SD", "Cirbp_202_SD", "Cirbp_201_SD", 
                             "Cirbp_204_SD", "Cirbp_206_SD", "Cirbp_205_SD")

Cirbp_Matrix$Count <- c(Cirbp_207_HC, Cirbp_202_HC, Cirbp_201_HC, Cirbp_204_HC, 
                        Cirbp_206_HC, Cirbp_205_HC, Cirbp_207_SD, Cirbp_202_SD, 
                        Cirbp_201_SD, Cirbp_204_SD, Cirbp_206_SD, Cirbp_205_SD)

Cirbp_Matrix$Condition <- rep(c("HC", "SD"), times = c(6,6))

Cirbp_Matrix$Transcript <- c("Cirbp-207", "Cirbp-202", "Cirbp-201", "Cirbp-204", 
                             "Cirbp-206", "Cirbp-205")

Cirbp_DF <- as.data.frame(Cirbp_Matrix)

# Set Colors
Cirbp_201_Color <- brewer.pal(8, "Greys")[4]
Cirbp_202_Color <- brewer.pal(8, "Set2")[7] 
Cirbp_204_Color <- brewer.pal(8, "Greys")[4]
Cirbp_205_Color <- brewer.pal(8, "Greys")[4]
Cirbp_206_Color <- brewer.pal(8, "Greys")[4]
Cirbp_207_Color <- brewer.pal(8, "Greys")[4]

Cirbp_Interaction <- ggplot(Cirbp_DF) +
  aes(x = Condition, y = Count, color = Transcript, group = Transcript,
      label = Transcript) +
  geom_text(nudge_x = 0.3, color = "black") +
  geom_line(linewidth = rep(c(1,2,1), times = c(2,2,8)), aes(group = Transcript, color = Transcript)) +
  scale_color_manual(values=c(Cirbp_201_Color, Cirbp_202_Color, Cirbp_204_Color, 
                              Cirbp_205_Color, Cirbp_206_Color, Cirbp_207_Color)) + # List numbers in alphabetical order
  geom_point(shape = rep(c(19,17), times = c(6,6)), size = rep(c(1,2,1,2,1), times = c(1,1,5,1,4)), 
             aes(color = Transcript)) +
  ylab(NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14), axis.text = element_text(size=14), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  labs(x = NULL)

Cirbp_Interaction 


#### Timeless ####

# Transcripts
# "ENSMUST00000055539" = Timeless-201
# "ENSMUST00000145710" = Timeless-214
# "ENSMUST00000105244" = Timeless-205
# "ENSMUST00000105245" = Timeless-206
# "ENSMUST00000105240" = Timeless-202

# Subset transcripts of Timeless
Timeless_Med <- dplyr::filter(infMed, Gene_Name == "Timeless")
# A tibble: 5 × 19
# Gene_ID           Gene_ID_Version Transcript_ID Transcript_ID_Version Gene_Description Gene_Name
# <chr>             <chr>           <chr>         <chr>                 <chr>            <chr>    
# 1 ENSMUSG000000399… ENSMUSG0000003… ENSMUST00000… ENSMUST00000055539.10 timeless circad… Timeless 
# 2 ENSMUSG000000399… ENSMUSG0000003… ENSMUST00000… ENSMUST00000145710.7  timeless circad… Timeless 
# 3 ENSMUSG000000399… ENSMUSG0000003… ENSMUST00000… ENSMUST00000105244.7  timeless circad… Timeless 
# 4 ENSMUSG000000399… ENSMUSG0000003… ENSMUST00000… ENSMUST00000105245.2  timeless circad… Timeless 
# 5 ENSMUSG000000399… ENSMUSG0000003… ENSMUST00000… ENSMUST00000105240.7  timeless circad… Timeless 
# ℹ 13 more variables: Transcript_Name <chr>, WTHC5_PFC_1_quant <dbl>, WTHC5_PFC_2_quant <dbl>,
#   WTHC5_PFC_3_quant <dbl>, WTHC5_PFC_4_quant <dbl>, WTHC5_PFC_5_quant <dbl>, HC_mean <dbl>,
#   WTSD5_PFC_1_quant <dbl>, WTSD5_PFC_2_quant <dbl>, WTSD5_PFC_3_quant <dbl>,
#   WTSD5_PFC_4_quant <dbl>, WTSD5_PFC_5_quant <dbl>, SD_mean <dbl>

Timeless_Med$Transcript_Name
# [1] "Timeless-201" "Timeless-214" "Timeless-205" "Timeless-206" "Timeless-202"

# We will now take the mean of these median values (columns 1-5 are control
# animals- homecage, columns 6-10 are sleep deprived animals)
# HC: Homecage
# SD: Sleep Deprived
Timeless_201 <- dplyr::filter(Timeless_Med, Transcript_Name == "Timeless-201")
Timeless_201_HC <- Timeless_201$HC_mean
Timeless_201_SD <- Timeless_201$SD_mean

Timeless_214 <- dplyr::filter(Timeless_Med, Transcript_Name == "Timeless-214")
Timeless_214_HC <- Timeless_214$HC_mean
Timeless_214_SD <- Timeless_214$SD_mean

Timeless_205 <- dplyr::filter(Timeless_Med, Transcript_Name == "Timeless-205")
Timeless_205_HC <- Timeless_205$HC_mean
Timeless_205_SD <- Timeless_205$SD_mean

Timeless_206 <- dplyr::filter(Timeless_Med, Transcript_Name == "Timeless-206")
Timeless_206_HC <- Timeless_206$HC_mean
Timeless_206_SD <- Timeless_206$SD_mean

Timeless_202 <- dplyr::filter(Timeless_Med, Transcript_Name == "Timeless-202")
Timeless_202_HC <- Timeless_202$HC_mean
Timeless_202_SD <- Timeless_202$SD_mean


Timeless_Matrix <- data.frame(1:10, 1:10, 1:10)
colnames(Timeless_Matrix) <- c("Count", "Condition", "Transcript")
row.names(Timeless_Matrix) <- c("Timeless_201_HC", "Timeless_214_HC", "Timeless_205_HC", 
                                "Timeless_206_HC", "Timeless_202_HC", "Timeless_201_SD", 
                                "Timeless_214_SD", "Timeless_205_SD", "Timeless_206_SD", 
                                "Timeless_202_SD")

Timeless_Matrix$Count <- c(Timeless_201_HC, Timeless_214_HC, Timeless_205_HC, Timeless_206_HC, 
                           Timeless_202_HC,  Timeless_201_SD, Timeless_214_SD, Timeless_205_SD, 
                           Timeless_206_SD, Timeless_202_SD)

Timeless_Matrix$Condition <- rep(c("HC", "SD"), times = c(5,5))

Timeless_Matrix$Transcript <- c("Timeless-201", "Timeless-214", "Timeless-205", 
                                "Timeless-206", "Timeless-202")

Timeless_DF <- as.data.frame(Timeless_Matrix)

# Set Colors
Timeless_201_Color <- brewer.pal(9, "Set1")[2]
Timeless_202_Color <- brewer.pal(9, "Set1")[3]
Timeless_205_Color <- brewer.pal(9, "Greys")[4]
Timeless_206_Color <- brewer.pal(9, "Greys")[4]
Timeless_214_Color <- brewer.pal(9, "Greys")[4]


Timeless_Interaction <- ggplot(Timeless_DF) +
  aes(x = Condition, y = Count, color = Transcript, group = Transcript,
      label = Transcript) +
  geom_text(nudge_x = 0.3, color = "black") +
  geom_line(linewidth = rep(c(2,1), times = c(4,6)), aes(group = Transcript, color = Transcript)) +
  scale_color_manual(values=c(Timeless_201_Color, Timeless_202_Color, 
                              Timeless_205_Color, Timeless_206_Color, Timeless_214_Color)) + # List numbers in alphabetical order
  geom_point(shape = rep(c(19,17), times = c(5,5)), size = rep(c(2,1,2,1,2), times = c(1,3,2,3,1))
             , aes(color = Transcript)) +
  ylab("Average Proportion") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14), axis.text = element_text(size=14), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = NULL)

Timeless_Interaction

#### Dbp ####

# Subset transcripts of Dbp
Dbp_Med <- dplyr::filter(infMed, Gene_Name == "Dbp")
# A tibble: 5 × 19
# Gene_ID           Gene_ID_Version Transcript_ID Transcript_ID_Version Gene_Description Gene_Name
# <chr>             <chr>           <chr>         <chr>                 <chr>            <chr>    
# 1 ENSMUSG000000598… ENSMUSG0000005… ENSMUST00000… ENSMUST00000080885.11 D site albumin … Dbp      
# 2 ENSMUSG000000598… ENSMUSG0000005… ENSMUST00000… ENSMUST00000211513.1  D site albumin … Dbp      
# 3 ENSMUSG000000598… ENSMUSG0000005… ENSMUST00000… ENSMUST00000211357.1  D site albumin … Dbp      
# 4 ENSMUSG000000598… ENSMUSG0000005… ENSMUST00000… ENSMUST00000211748.1  D site albumin … Dbp      
# 5 ENSMUSG000000598… ENSMUSG0000005… ENSMUST00000… ENSMUST00000210120.1  D site albumin … Dbp      
# ℹ 13 more variables: Transcript_Name <chr>, WTHC5_PFC_1_quant <dbl>, WTHC5_PFC_2_quant <dbl>,
#   WTHC5_PFC_3_quant <dbl>, WTHC5_PFC_4_quant <dbl>, WTHC5_PFC_5_quant <dbl>, HC_mean <dbl>,
#   WTSD5_PFC_1_quant <dbl>, WTSD5_PFC_2_quant <dbl>, WTSD5_PFC_3_quant <dbl>,
#   WTSD5_PFC_4_quant <dbl>, WTSD5_PFC_5_quant <dbl>, SD_mean <dbl>

Dbp_Med$Transcript_Name
# [1] "Dbp-201" "Dbp-204" "Dbp-203" "Dbp-205" "Dbp-202"

# We will now take the mean of these median values (columns 1-5 are control
# animals- homecage, columns 6-10 are sleep deprived animals)
# HC: Homecage
# SD: Sleep Deprived
Dbp_201 <- dplyr::filter(Dbp_Med, Transcript_Name == "Dbp-201")
Dbp_201_HC <- Dbp_201$HC_mean
Dbp_201_SD <- Dbp_201$SD_mean

Dbp_204 <- dplyr::filter(Dbp_Med, Transcript_Name == "Dbp-204")
Dbp_204_HC <- Dbp_204$HC_mean
Dbp_204_SD <- Dbp_204$SD_mean

Dbp_203 <- dplyr::filter(Dbp_Med, Transcript_Name == "Dbp-203")
Dbp_203_HC <- Dbp_203$HC_mean
Dbp_203_SD <- Dbp_203$SD_mean

Dbp_205 <- dplyr::filter(Dbp_Med, Transcript_Name == "Dbp-205")
Dbp_205_HC <- Dbp_205$HC_mean
Dbp_205_SD <- Dbp_205$SD_mean

Dbp_202 <- dplyr::filter(Dbp_Med, Transcript_Name == "Dbp-202")
Dbp_202_HC <- Dbp_202$HC_mean
Dbp_202_SD <- Dbp_202$SD_mean


Dbp_Matrix <- data.frame(1:10, 1:10, 1:10)
colnames(Dbp_Matrix) <- c("Count", "Condition", "Transcript")
row.names(Dbp_Matrix) <- c("Dbp_201_HC", "Dbp_204_HC", "Dbp_203_HC", "Dbp_205_HC", 
                           "Dbp_202_HC", "Dbp_201_SD", "Dbp_204_SD", "Dbp_203_SD", 
                           "Dbp_205_SD", "Dbp_202_SD")


Dbp_Matrix$Count <- c(Dbp_201_HC, Dbp_204_HC, Dbp_203_HC, Dbp_205_HC, Dbp_202_HC,  
                      Dbp_201_SD, Dbp_204_SD, Dbp_203_SD, Dbp_205_SD, Dbp_202_SD)

Dbp_Matrix$Condition <- rep(c("HC", "SD"), times = c(5,5))

Dbp_Matrix$Transcript <- c("Dbp-201", "Dbp-204", "Dbp-203", "Dbp-205", "Dbp-202")

Dbp_DF <- as.data.frame(Dbp_Matrix)

# Set Colors
Dbp_201_Color <- brewer.pal(9, "Set1")[4]
Dbp_202_Color <- brewer.pal(9, "Greys")[4]
Dbp_203_Color <- brewer.pal(9, "Greys")[4]
Dbp_204_Color <- brewer.pal(9, "Greys")[4]
Dbp_205_Color <- brewer.pal(9, "Greys")[4]


Dbp_Interaction <- ggplot(Dbp_DF) +
  aes(x = Condition, y = Count, color = Transcript, group = Transcript,
      label = Transcript) +
  geom_text(nudge_x = 0.3, color = "black") +
  geom_line(linewidth = rep(c(2,1), times = c(2,8)), aes(group = Transcript, color = Transcript)) +
  scale_color_manual(values=c(Dbp_201_Color, Dbp_202_Color, Dbp_203_Color, 
                              Dbp_204_Color, Dbp_205_Color)) + # List numbers in alphabetical order
  geom_point(shape = rep(c(19,17), times = c(5,5)), size = rep(c(2,1,2,1), times = c(1,4,1,4)), 
             aes(color = Transcript)) +
  ylab("Average Proportion") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14), axis.text = element_text(size=14), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = NULL)

Dbp_Interaction

#### Arc ####
Arc_Med <- dplyr::filter(infMed, Gene_Name == "Arc")
# A tibble: 2 × 19
# Gene_ID           Gene_ID_Version Transcript_ID Transcript_ID_Version Gene_Description Gene_Name
# <chr>             <chr>           <chr>         <chr>                 <chr>            <chr>    
# 1 ENSMUSG000000226… ENSMUSG0000002… ENSMUST00000… ENSMUST00000023268.13 activity regula… Arc      
# 2 ENSMUSG000000226… ENSMUSG0000002… ENSMUST00000… ENSMUST00000110009.3  activity regula… Arc      
# ℹ 13 more variables: Transcript_Name <chr>, WTHC5_PFC_1_quant <dbl>, WTHC5_PFC_2_quant <dbl>,
#   WTHC5_PFC_3_quant <dbl>, WTHC5_PFC_4_quant <dbl>, WTHC5_PFC_5_quant <dbl>, HC_mean <dbl>,
#   WTSD5_PFC_1_quant <dbl>, WTSD5_PFC_2_quant <dbl>, WTSD5_PFC_3_quant <dbl>,
#   WTSD5_PFC_4_quant <dbl>, WTSD5_PFC_5_quant <dbl>, SD_mean <dbl>

Arc_Med$Transcript_Name
# [1] "Arc-201" "Arc-202"

# We will now take the mean of these median values (columns 1-5 are control
# animals- homecage, columns 6-10 are sleep deprived animals)
# HC: Homecage
# SD: Sleep Deprived

Arc_201 <- dplyr::filter(Arc_Med, Transcript_Name == "Arc-201")
Arc_201_HC <- Arc_201$HC_mean
Arc_201_SD <- Arc_201$SD_mean

Arc_202 <- dplyr::filter(Arc_Med, Transcript_Name == "Arc-202")
Arc_202_HC <- Arc_202$HC_mean
Arc_202_SD <- Arc_202$SD_mean

Arc_Matrix <- data.frame(1:4, 1:4, 1:4)
colnames(Arc_Matrix) <- c("Count", "Condition", "Transcript")
row.names(Arc_Matrix) <- c("Arc_201_HC", "Arc_202_HC", "Arc_201_SD", "Arc_202_SD")

Arc_Matrix$Count <- c(Arc_201_HC, Arc_202_HC, Arc_201_SD, Arc_202_SD)

Arc_Matrix$Condition <- rep(c("HC", "SD"), times = c(2,2))

Arc_Matrix$Transcript <- c("Arc-201", "Arc-202")

Arc_DF <- as.data.frame(Arc_Matrix)

# Set Colors
Arc201_Color <- brewer.pal(9, "Set1")[5]
Arc202_Color  <- brewer.pal(9, "Set1")[6]

Arc_Interaction <- ggplot(Arc_DF) +
  aes(x = Condition, y = Count, color = Transcript, group = Transcript,
      label = Transcript) +
  geom_text(nudge_x = 0.3, color = "black") +
  geom_line(linewidth = 2, aes(group = Transcript, color = Transcript)) +
  scale_color_manual(values=c(Arc201_Color, Arc202_Color)) +
  geom_point(shape = rep(c(19,17), times = c(2,2)), size = 2, 
             aes(color = Transcript)) +
  ylab(NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14), axis.text = element_text(size=14), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  labs(x = NULL)

Arc_Interaction

#### Eif2ak4 ####
Eif2ak4_Med <- dplyr::filter(infMed, Gene_Name == "Eif2ak4")
# A tibble: 4 × 19
# Gene_ID           Gene_ID_Version Transcript_ID Transcript_ID_Version Gene_Description Gene_Name
# <chr>             <chr>           <chr>         <chr>                 <chr>            <chr>    
# 1 ENSMUSG000000051… ENSMUSG0000000… ENSMUST00000… ENSMUST00000005233.11 eukaryotic tran… Eif2ak4  
# 2 ENSMUSG000000051… ENSMUSG0000000… ENSMUST00000… ENSMUST00000102527.9  eukaryotic tran… Eif2ak4  
# 3 ENSMUSG000000051… ENSMUSG0000000… ENSMUST00000… ENSMUST00000110870.7  eukaryotic tran… Eif2ak4  
# 4 ENSMUSG000000051… ENSMUSG0000000… ENSMUST00000… ENSMUST00000110869.1  eukaryotic tran… Eif2ak4  
# ℹ 13 more variables: Transcript_Name <chr>, WTHC5_PFC_1_quant <dbl>, WTHC5_PFC_2_quant <dbl>,
#   WTHC5_PFC_3_quant <dbl>, WTHC5_PFC_4_quant <dbl>, WTHC5_PFC_5_quant <dbl>, HC_mean <dbl>,
#   WTSD5_PFC_1_quant <dbl>, WTSD5_PFC_2_quant <dbl>, WTSD5_PFC_3_quant <dbl>,
#   WTSD5_PFC_4_quant <dbl>, WTSD5_PFC_5_quant <dbl>, SD_mean <dbl>
Eif2ak4_Med$Transcript_Name
# [1] "Eif2ak4-201" "Eif2ak4-202" "Eif2ak4-204" "Eif2ak4-203"

# We will now take the mean of these median values (columns 1-5 are control
# animals- homecage, columns 6-10 are sleep deprived animals)
# HC: Homecage
# SD: Sleep Deprived

Eif2ak4_201 <- dplyr::filter(Eif2ak4_Med, Transcript_Name == "Eif2ak4-201")
Eif2ak4_201_HC <- Eif2ak4_201$HC_mean
Eif2ak4_201_SD <- Eif2ak4_201$SD_mean

Eif2ak4_202 <- dplyr::filter(Eif2ak4_Med, Transcript_Name == "Eif2ak4-202")
Eif2ak4_202_HC <- Eif2ak4_202$HC_mean
Eif2ak4_202_SD <- Eif2ak4_202$SD_mean

Eif2ak4_203 <- dplyr::filter(Eif2ak4_Med, Transcript_Name == "Eif2ak4-203")
Eif2ak4_203_HC <- Eif2ak4_203$HC_mean
Eif2ak4_203_SD <- Eif2ak4_203$SD_mean

Eif2ak4_204 <- dplyr::filter(Eif2ak4_Med, Transcript_Name == "Eif2ak4-204")
Eif2ak4_204_HC <- Eif2ak4_204$HC_mean
Eif2ak4_204_SD <- Eif2ak4_204$SD_mean

Eif2ak4_Matrix <- data.frame(1:8, 1:8, 1:8)
colnames(Eif2ak4_Matrix) <- c("Count", "Condition", "Transcript")
row.names(Eif2ak4_Matrix) <- c("Eif2ak4_201_HC", "Eif2ak4_202_HC", "Eif2ak4_203_HC", 
                                "Eif2ak4_204_HC", "Eif2ak4_201_SD", "Eif2ak4_202_SD", 
                                "Eif2ak4_203_SD", "Eif2ak4_204_SD")

Eif2ak4_Matrix$Count <- c(Eif2ak4_201_HC, Eif2ak4_202_HC, Eif2ak4_203_HC, 
                           Eif2ak4_204_HC, Eif2ak4_201_SD, Eif2ak4_202_SD, 
                           Eif2ak4_203_SD, Eif2ak4_204_SD)

Eif2ak4_Matrix$Condition <- rep(c("HC", "SD"), times = c(4,4))

Eif2ak4_Matrix$Transcript <- c("Eif2ak4-201", "Eif2ak4-202", "Eif2ak4-203", "Eif2ak4-204")

Eif2ak4_DF <- as.data.frame(Eif2ak4_Matrix)

Eif2ak4_201_Color <- brewer.pal(9, "Set1")[1]
Eif2ak4_202_Color <- brewer.pal(9, "Set1")[8]
Eif2ak4_203_Color <- brewer.pal(9, "Greys")[4]
Eif2ak4_204_Color <- brewer.pal(9, "Greys")[4]

Eif2ak4_Interaction <- ggplot(Eif2ak4_DF) +
  aes(x = Condition, y = Count, color = Transcript, group = Transcript,
      label = Transcript) +
  geom_text(nudge_x = 0.3, color = "black") +
  geom_line(linewidth = rep(c(2,1), times = c(4,4)), aes(group = Transcript, color = Transcript)) +
  scale_color_manual(values=c(Eif2ak4_201_Color, Eif2ak4_202_Color, Eif2ak4_203_Color, 
                              Eif2ak4_204_Color)) + # List numbers in alphabetical order
  geom_point(shape = rep(c(19,17), times = c(4,4)), size = rep(c(2,1,2,1), times = c(2,2,2,2)), 
             aes(color = Transcript)) +
  ylab(NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 14), axis.text = element_text(size=14), 
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  labs(x = NULL)

Eif2ak4_Interaction

# Now arrange all supplemental figures in a grid
# For now, we are going to separate into upregulated and downregulated at the
# gene level (keep timeless separate as it is not detected), to minimize white
# space on the final figure

# Up
Homer1_Interaction + Bdnf_Interaction + Arc_Interaction + patchwork::plot_layout(ncol = 3)

# Down
Dbp_Interaction + Cirbp_Interaction + Eif2ak4_Interaction + patchwork::plot_layout(ncol = 3)

# Timeless
Timeless_Interaction


# Alternative for layout

design <- "
##1
234
567
"

Timeless_Interaction + Homer1_Interaction + Bdnf_Interaction + Arc_Interaction + 
Dbp_Interaction + Cirbp_Interaction + Eif2ak4_Interaction + patchwork::plot_layout(design = design)


