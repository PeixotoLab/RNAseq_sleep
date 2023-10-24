# Differential expression analysis ####
# The goal of this analysis is to identify
# the differentially expressed genes.
# The raw counts were normalized with upper-quartile and 
# RUV (k = 2) methods.

# Set up
library(SingleCellExperiment)
library(scuttle)
library(edgeR)
library(ggplot2)
library(RUVSeq)
library(ggrepel)
library(EDASeq)

# Load pseudo-bulk data
pb <- readRDS("snrna_pb.rds")
# For the differential expression analysis, we selected the neuronal labels
pb <- pb[, pb$class == "Glutamatergic" | pb$class == "GABAergic"]

# Create a list of pseudo-bulk for each cell-type
pb_ct <- lapply(unique(pb$azimuth_labels), function(ct) pb[, pb$azimuth_labels == ct])

# UQ + RUV normalization ####
# The 10% of meta-analysis negative control genes were randomly selected.
neg_ctrl <- read.table("/mnt/callisto/Zuin/SD_Negative_Controls.txt")
set.seed(23)
random_ng <- sample(neg_ctrl$x, round(length(neg_ctrl$x) * 0.1))
# The remaining control genes were used to fit RUV normalization.
neg_ctrl <- neg_ctrl[!(neg_ctrl$x %in% random_ng), ]
neg_ctrl <- intersect(rownames(pb), neg_ctrl)

# UQ normalization 
ct_counts <- lapply(pb_ct, function(x) as.matrix(counts(x)))
uq <- lapply(ct_counts, function(u) betweenLaneNormalization(u, which="upper"))

# RUV normalization
ruv2_expr_data <- ruv2_w <- list()
for (i in seq_along(unique(pb$azimuth_labels))) {
  # A matrix specifying the replicates constructed
  groups <- makeGroups(pb_ct[[i]]$condition)
  
  # The matrix of normalized counts was saved in a list object
  ruv2_expr_data[[i]] <- RUVs(uq[[i]], 
                              cIdx = neg_ctrl, scIdx = groups, k = 2)
  # The factors of unwanted variation were saved in a list object
  ruv2_w[[i]] <- ruv2_expr_data[[i]]$W
}
save(ruv2_expr_data, ruv2_w, file = "UQ_RUVk2_result.RData")

# DEA with UQ + RUVk2 ####
res_uq_ruv2 <- df_uq_ruv2 <- list()
for (i in seq_along(unique(pb$azimuth_labels))) {
  y <- DGEList(counts(pb_ct[[i]]), samples = colData(pb_ct[[i]]))
  
  # The genes were filtered for each cell-type
  keep <- filterByExpr(y, group = pb_ct[[i]]$condition)
  y <- y[keep, ]
  # Upper-quartile normalization
  y <- calcNormFactors(y, method = "upperquartile")
  
  # Design Matrix
  design <- model.matrix(~0 + y$samples$condition + ruv2_w[[i]], y$samples)
  colnames(design) <- c("HC", "SD", "W_1", "W_2")
  
  y <- estimateDisp(y, design)
  
  fit <- glmFit(y, design)
  
  # Contrast creation (SD vs HC)
  contrast <- makeContrasts(SD - HC, levels = design)
  
  res_uq_ruv2[[i]] <- glmLRT(fit, contrast = contrast)
  
  df_uq_ruv2[[i]] <- as.data.frame(res_uq_ruv2[[i]]$table)
  df_uq_ruv2[[i]] <- df_uq_ruv2[[i]][order(df_uq_ruv2[[i]]$PValue), ]
  
  FDR <- as.data.frame(topTags(res_uq_ruv2[[i]], n = length(rownames(pb)))[, 5])
  
  df_uq_ruv2[[i]] <- cbind(df_uq_ruv2[[i]], FDR)
  df_uq_ruv2[[i]] <- df_uq_ruv2[[i]][order(df_uq_ruv2[[i]]$FDR), ]
}
saveRDS(res_uq_ruv2, file = "snRNA_UQ_RUVs2_Results.rds")
saveRDS(df_uq_ruv2, file = "snRNA_UQ_RUVs2_DFResults.rds")

# N. of nuclei
n_nuclei <- matrix(pb$ncells, 11, 6, byrow = TRUE)
n_nuclei <- rowSums(n_nuclei)

n_genes <- lapply(df_uq_ruv2, function(u) length(rownames(u)))

# DEA results for differentially expressed genes were selected
df_degs <- lapply(df_uq_ruv2, function(u) u[u$FDR < 0.05, ])

# N. of DEGs
n_degs <- lapply(df_degs, function(u) length(rownames(u)))

df_up_degs <- lapply(df_degs, function(u) u[u$logFC >= 0, ])
# N. of up-regulated DEGs
n_up_degs <- lapply(df_up_degs, function(u) length(rownames(u)))

df_down_degs <- lapply(df_degs, function(u) u[u$logFC < 0, ])
# N. of down-regulated DEGs
n_down_degs <- lapply(df_down_degs, function(u) length(rownames(u)))

# Load meta-analysis positive control genes
posctrl <- read.table("/mnt/callisto/Zuin/Additional_File2_Positive_Controls.txt", header = TRUE)

# N. of meta-analysis positive control genes inside single-nuclear dataset
n_posctrl <- lapply(df_uq_ruv2, function(u) length(intersect(rownames(u), posctrl$Gene_ID)))
# N. of meta-analysis positive control DE genes
n_de_posctrl <- lapply(df_degs, function(u) length(intersect(rownames(u), posctrl$Gene_ID)))

tab <- as.data.frame(cbind(n_nuclei, as.matrix(n_up_degs), as.matrix(n_down_degs), as.matrix(n_de_posctrl)))
rownames(tab) <- unique(pb$azimuth_labels)

# Load Allen color palette
load("/mnt/callisto/Zuin/AllenColorLabel.RData")
subclass_color <- data.frame(subclass_color)
subclass_color$Label <- rownames(subclass_color)
colnames(subclass_color) <- c("color", "Label")
rm_labels <- subclass_color[rownames(subclass_color) %in% 
                              intersect(rownames(subclass_color), levels(factor(pb$azimuth_labels))), ]

neuronal_color <- rm_labels$color
names(neuronal_color) <- rm_labels$Label

tab <- as.data.frame(cbind(n_nuclei, as.matrix(n_degs)))
colnames(tab) <- c("nuclei", "DEG")
tab$Label <- levels(factor(pb$azimuth_labels))
tab$DEG <- as.numeric(tab$DEG)
tab$nuclei <- as.numeric(tab$nuclei)

# Scatter plot of log(nuclei) vs log(DEGs)
p1 <- ggplot(tab, aes(x = log(nuclei), y = log(DEG))) +
  geom_point(aes(color = Label), size = 3) +
  scale_color_manual(values = neuronal_color) + theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 7), axis.title = element_text(size = 7)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, col = "grey") +
  geom_label_repel(aes(label = Label), size = 2,
                   box.padding   = 0.5,
                   point.padding = 0,
                   segment.color = "grey50",min.segment.length = unit(1.5, 'lines'), colour = neuronal_color)
ggsave("Peixoto_Figure5_part1.jpg",  width = 9, height = 9, units = "cm")
ggsave("Peixoto_Figure5_part1.pdf",  width = 9, height = 9, units = "cm")

## PCA ####
label <- levels(factor(pb$azimuth_labels))
set_ct <- lapply(pb_ct, function(u) newSeqExpressionSet(as.matrix(round(counts(u))),
                                                        phenoData = data.frame(colData(u), 
                                                                               row.names = colnames(u))))

set_ct <- lapply(set_ct, function(u) betweenLaneNormalization(u, which="upper"))

set_ruv2 <- list()
for (i in seq_along(unique(pb$azimuth_labels))) {
  # A matrix specifying the replicates constructed
  groups <- makeGroups(pb_ct[[i]]$condition)
  
  # k = 2
  # The matrix of normalized counts was saved in a list object
  set_ruv2[[i]] <- RUVs(set_ct[[i]], cIdx = neg_ctrl, scIdx = groups, k = 2)
  
}
save(set_ruv2, file = "snrna_set_uq_ruv2.RData")

file_name <- paste("Peixoto_Figure5_Supplement1_part1.jpg", sep = "")
jpeg(file = file_name, width = 794, height = 1134, units = "px")
# pdf(file = file_name, width = 8.27, height = 11.81)
par(mfrow = c(4, 3))
for (i in seq_along(label)) {
  plotPCA(set_ruv2[[i]], xlim = c(-1, 1), label = FALSE, pch = 20, theme_size = 7,
          col = as.numeric(as.factor(set_ruv2[[i]]$condition)), cex = 2,
          main = paste(label[i], sep = " "))
}
dev.off()

### Volcano plot of L4/5 and Pvalb with highlighted genes ####
# Load a smaller list of meta-analysis positive control genes
highlight_genes <- readxl::read_excel("050123_Highlight_Genes_Subset.xlsx")
# highlight_genes <- highlight_genes[!is.na(highlight_genes$Ensembl_ID), ]
# highlight_genes <- highlight_genes[!duplicated(highlight_genes$Gene_Name), ]
highlight_genes <- highlight_genes[order(highlight_genes$Ensembl_ID), ]

# Select L4/5 IT CTX and Pvalb DEA results
ll <- list(df_uq_ruv2[[2]], df_uq_ruv2[[9]])
ll <- lapply(ll, function(x) x[order(rownames(x)), ])

highlight_ll <- list(L45 = as.data.frame(highlight_genes[highlight_genes$Ensembl_ID %in% intersect(rownames(ll[[1]]), highlight_genes$Ensembl_ID), ]),
                     Pvalb = as.data.frame(highlight_genes[highlight_genes$Ensembl_ID %in% intersect(rownames(ll[[2]]), highlight_genes$Ensembl_ID), ]))

for (i in seq_along(ll)) {
  ll[[i]] <- ll[[i]][order(rownames(ll[[i]])),]
  ll[[i]]$symbol_highlight <- ""
  ll[[i]]$symbol_highlight[rownames(ll[[i]]) %in% highlight_ll[[i]]$Ensembl_ID] <- highlight_ll[[i]]$Gene_Name
  ll[[i]]$Significance <- "No Significant"
  ll[[i]]$Significance[ll[[i]]$FDR < 0.05] <- "Significant"
  inter <- intersect(rownames(ll[[i]][ll[[i]]$FDR < 0.05, ]), posctrl$Gene_ID)
  ll[[i]]$Significance[rownames(ll[[i]]) %in% inter] <- "SignificantPos"
  ll[[i]]$highlight_genes <- ll[[i]]$symbol_highlight
  ll[[i]]$highlight_genes[ll[[i]]$highlight_genes == ""] <- FALSE
  ll[[i]]$highlight_genes[ll[[i]]$highlight_genes != FALSE] <- TRUE
}

p1 <- ggplot(data = ll[[1]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  geom_point(data= ll[[1]][(ll[[1]]$Significance == "Significant"|ll[[1]]$Significance == "No Significant"), ], size = 1) + 
  geom_point(data= ll[[1]][ll[[1]]$Significance == "SignificantPos", ], size = 1) + 
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none") +
  scale_color_manual(values = c("Significant" = "black", "No Significant" = "grey", "SignificantPos" = "red")) +
  ggtitle(paste("L4/5 IT CTX", sep = "")) +
  geom_label_repel(data = ll[[1]][(ll[[1]]$Significance == "Significant" & ll[[1]]$logFC < 0), ],
                   aes(label = ifelse(highlight_genes, symbol_highlight, "")),
                   fill = "white", max.overlaps = Inf, fontface = 2, direction = "both", hjust = "left",
                   ylim = c(3, 15), show.legend = FALSE, size = 2) +
  geom_label_repel(data = ll[[1]][(ll[[1]]$Significance == "SignificantPos" & ll[[1]]$logFC >= 0), ],
                   aes(label = ifelse(highlight_genes, symbol_highlight, "")),
                   fill = "white", max.overlaps = Inf, fontface = 2, direction = "both", ylim = c(10, 40), show.legend = FALSE, size = 2) +
  
  geom_label_repel(data = ll[[1]][(ll[[1]]$Significance == "Significant" & ll[[1]]$logFC >= 0), ],
                   aes(label = ifelse(highlight_genes, symbol_highlight, "")),
                   fill = "white", max.overlaps = Inf, fontface = 2, direction = "both", hjust = "right",
                   ylim = c(3, 20), show.legend = FALSE, size = 2) +
  geom_label_repel(data = ll[[1]][(ll[[1]]$Significance == "SignificantPos"& ll[[1]]$logFC < 0), ],
                   aes(label = ifelse(highlight_genes, symbol_highlight, "")),
                   fill = "white", max.overlaps = Inf, fontface = 2, direction = "both", ylim = c(5, 20), show.legend = FALSE, size = 2) +
  
  annotate("text", x = -4, y = 2.5, label = n_degs[[2]]) +
  annotate("text", x = -4, y = 0.5, label = n_de_posctrl[[2]], color = "red")

p2 <- ggplot(data = ll[[2]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  geom_point(data= ll[[2]][(ll[[2]]$Significance == "Significant"|ll[[2]]$Significance == "No Significant"), ], size = 1) + 
  geom_point(data= ll[[2]][ll[[2]]$Significance == "SignificantPos", ], size = 1) + 
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none") +
  scale_color_manual(values = c("Significant" = "black", "No Significant" = "grey", "SignificantPos" = "red")) +
  ggtitle(paste("Pvalb", sep = "")) +
  geom_label_repel(data = ll[[2]][(ll[[2]]$Significance == "Significant" & ll[[2]]$logFC > 0), ],
                   aes(label = ifelse(highlight_genes, symbol_highlight, "")),
                   fill = "white", max.overlaps = Inf, fontface = 2, direction = "both", show.legend = FALSE, size = 2) +
  geom_label_repel(data = ll[[2]][(ll[[2]]$Significance == "Significant" & ll[[2]]$logFC < 0), ],
                   aes(label = ifelse(highlight_genes, symbol_highlight, "")),
                   fill = "white", max.overlaps = Inf, fontface = 2, direction = "x", hjust = "left", ylim = c(3, 15), show.legend = FALSE, size = 2) +
  geom_label_repel(data = ll[[2]][(ll[[2]]$Significance == "SignificantPos" & ll[[2]]$logFC < 0), ],
                   aes(label = ifelse(highlight_genes, symbol_highlight, "")),
                   fill = "white", max.overlaps = Inf, fontface = 2, direction = "x", ylim = c(2, 10), show.legend = FALSE, size = 2) +
  geom_label_repel(data = ll[[2]][(ll[[2]]$Significance == "SignificantPos" & ll[[2]]$logFC > 0), ],
                   aes(label = ifelse(highlight_genes, symbol_highlight, "")),
                   fill = "white", max.overlaps = Inf, fontface = 2, direction = "both", ylim = c(5, 10), show.legend = FALSE, size = 2) +
  annotate("text", x = -4, y = 1.5, label = n_degs[[9]]) +
  annotate("text", x = -4, y = 0.5, label = n_de_posctrl[[9]], color = "red")

patch <- p1 | p2
ggsave("Peixoto_Figure5_part2.jpg",  width = 18, height = 9, units = "cm")
ggsave("Peixoto_Figure5_part2.pdf",  width = 18, height = 9, units = "cm")


## Volcano plot of all cell-types #####
for (i in seq_along(df_uq_ruv2)) {
  df_uq_ruv2[[i]]$Significance <- "No Significant"
  df_uq_ruv2[[i]]$Significance[df_uq_ruv2[[i]]$FDR < 0.05] <- "Significant"
  inter <- intersect(rownames(df_uq_ruv2[[i]][df_uq_ruv2[[i]]$FDR < 0.05, ]), posctrl$Gene_ID)
  df_uq_ruv2[[i]]$Significance[rownames(df_uq_ruv2[[i]]) %in% inter] <- "SignificantPos"
}
label <- levels(factor(pb$azimut_labels))

p1 <- ggplot(data = df_uq_ruv2[[1]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  geom_point(data= df_uq_ruv2[[1]][(df_uq_ruv2[[1]]$Significance == "Significant"| 
                                      df_uq_ruv2[[1]]$Significance == "No Significant") ,],size = 1) + 
  geom_point(data= df_uq_ruv2[[1]][df_uq_ruv2[[1]]$Significance == "SignificantPos",], size = 1) + 
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none") +
  scale_color_manual(values = c("Significant" = "black", "No Significant" = "grey", "SignificantPos" = "red")) +
  ggtitle(paste(label[1], sep = "")) + annotate("text", x = -4, y = 2.5, label = n_degs[[1]]) +
  annotate("text", x = -4, y = 0.5, label = n_de_posctrl[[1]], color = "red")

p2 <- ggplot(data = df_uq_ruv2[[2]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  geom_point(data= df_uq_ruv2[[2]][(df_uq_ruv2[[2]]$Significance == "Significant"|
                                      df_uq_ruv2[[2]]$Significance == "No Significant") ,],size = 1) + 
  geom_point(data= df_uq_ruv2[[2]][df_uq_ruv2[[2]]$Significance == "SignificantPos",], size = 1) + 
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none") +
  scale_color_manual(values = c("Significant" = "black", "No Significant" = "grey", "SignificantPos" = "red")) +
  ggtitle(paste(label[2], sep = "")) + 
  annotate("text", x = -4, y = 2.5, label = n_degs[[2]]) +
  annotate("text", x = -4, y = 0.5, label = n_de_posctrl[[2]], color = "red")

p3 <- ggplot(data = df_uq_ruv2[[3]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  geom_point(data= df_uq_ruv2[[3]][(df_uq_ruv2[[3]]$Significance == "Significant"|
                                      df_uq_ruv2[[3]]$Significance == "No Significant") ,],size = 1) + 
  geom_point(data= df_uq_ruv2[[3]][df_uq_ruv2[[3]]$Significance == "SignificantPos",], size = 1) + 
  
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none") +
  scale_color_manual(values = c("Significant" = "black", "No Significant" = "grey", "SignificantPos" = "red")) +
  ggtitle(paste(label[3], sep = "")) + annotate("text", x = -4, y = 2, label = n_degs[[3]]) +
  annotate("text", x = -4, y = 0.5, label = n_de_posctrl[[3]], color = "red")

p4 <- ggplot(data = df_uq_ruv2[[4]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  geom_point(data= df_uq_ruv2[[4]][(df_uq_ruv2[[4]]$Significance == "Significant"|
                                      df_uq_ruv2[[4]]$Significance == "No Significant") ,],size = 1) + 
  geom_point(data= df_uq_ruv2[[4]][df_uq_ruv2[[4]]$Significance == "SignificantPos",], size = 1) + 
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none") +
  scale_color_manual(values = c("Significant" = "black", "No Significant" = "grey", "SignificantPos" = "red")) +
  ggtitle(paste(label[4], sep = "")) + annotate("text", x = -4, y = 2, label = n_degs[[4]]) +
  annotate("text", x = -4, y = 0.5, label = n_de_posctrl[[4]], color = "red")

p5 <- ggplot(data = df_uq_ruv2[[5]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  geom_point(data= df_uq_ruv2[[5]][(df_uq_ruv2[[5]]$Significance == "Significant"|
                                      df_uq_ruv2[[5]]$Significance == "No Significant") ,],size = 1) + 
  geom_point(data= df_uq_ruv2[[5]][df_uq_ruv2[[5]]$Significance == "SignificantPos",], size = 1) + 
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none") +
  scale_color_manual(values = c("Significant" = "black", "No Significant" = "grey", "SignificantPos" = "red")) +
  ggtitle(paste(label[5], sep = "")) + annotate("text", x = -4, y = 1.5, label = n_degs[[5]]) +
  annotate("text", x = -4, y = 0.5, label = n_de_posctrl[[5]], color = "red")

p6 <- ggplot(data = df_uq_ruv2[[6]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  geom_point(data= df_uq_ruv2[[6]][(df_uq_ruv2[[6]]$Significance == "Significant"|
                                      df_uq_ruv2[[6]]$Significance == "No Significant") ,],size = 1) + 
  geom_point(data= df_uq_ruv2[[6]][df_uq_ruv2[[6]]$Significance == "SignificantPos",], size = 1) + 
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none") +
  scale_color_manual(values = c("Significant" = "black", "No Significant" = "grey", "SignificantPos" = "red")) +
  ggtitle(paste(label[6], sep = "")) + annotate("text", x = -4, y = 2, label = n_degs[[6]]) +
  annotate("text", x = -4, y = 0.5, label = n_de_posctrl[[6]], color = "red")

p7 <- ggplot(data = df_uq_ruv2[[7]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  geom_point(data= df_uq_ruv2[[7]][(df_uq_ruv2[[7]]$Significance == "Significant"|
                                      df_uq_ruv2[[7]]$Significance == "No Significant") ,],size = 1) + 
  geom_point(data= df_uq_ruv2[[7]][df_uq_ruv2[[7]]$Significance == "SignificantPos",], size = 1) + 
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none") +
  scale_color_manual(values = c("Significant" = "black", "No Significant" = "grey", "SignificantPos" = "red")) +
  ggtitle(paste(label[7], sep = "")) + annotate("text", x = -4, y = 2, label = n_degs[[7]]) +
  annotate("text", x = -4, y = 0.5, label = n_de_posctrl[[7]], color = "red")

p8 <- ggplot(data = df_uq_ruv2[[8]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  geom_point(data= df_uq_ruv2[[8]][(df_uq_ruv2[[8]]$Significance == "Significant"|
                                      df_uq_ruv2[[8]]$Significance == "No Significant") ,],size = 1) + 
  geom_point(data= df_uq_ruv2[[8]][df_uq_ruv2[[8]]$Significance == "SignificantPos",], size = 1) + 
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none") +
  scale_color_manual(values = c("Significant" = "black", "No Significant" = "grey", "SignificantPos" = "red")) +
  ggtitle(paste(label[8], sep = "")) + annotate("text", x = -4, y = 1.5, label = n_degs[[8]]) +
  annotate("text", x = -4, y = 0.5, label = n_de_posctrl[[8]], color = "red")

p9 <- ggplot(data = df_uq_ruv2[[9]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  geom_point(data= df_uq_ruv2[[9]][(df_uq_ruv2[[9]]$Significance == "Significant"|
                                      df_uq_ruv2[[9]]$Significance == "No Significant") ,],size = 1) + 
  geom_point(data= df_uq_ruv2[[9]][df_uq_ruv2[[9]]$Significance == "SignificantPos",], size = 1) + 
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none") +
  scale_color_manual(values = c("Significant" = "black", "No Significant" = "grey", "SignificantPos" = "red")) +
  ggtitle(paste(label[9], sep = "")) + annotate("text", x = -4, y = 1.5, label = n_degs[[9]]) +
  annotate("text", x = -4, y = 0.5, label = n_de_posctrl[[9]], color = "red")

p10 <- ggplot(data = df_uq_ruv2[[10]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  geom_point(data= df_uq_ruv2[[10]][(df_uq_ruv2[[10]]$Significance == "Significant"|
                                       df_uq_ruv2[[10]]$Significance == "No Significant") ,],size = 1) + 
  geom_point(data= df_uq_ruv2[[10]][df_uq_ruv2[[10]]$Significance == "SignificantPos",], size = 1) + 
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none") +
  scale_color_manual(values = c("Significant" = "black", "No Significant" = "grey", "SignificantPos" = "red")) +
  ggtitle(paste(label[10], sep = "")) + annotate("text", x = -4, y = 1.5, label = n_degs[[10]]) +
  annotate("text", x = -4, y = 0.5, label = n_de_posctrl[[10]], color = "red")

p11 <- ggplot(data = df_uq_ruv2[[11]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  geom_point(data= df_uq_ruv2[[11]][(df_uq_ruv2[[11]]$Significance == "Significant"|
                                       df_uq_ruv2[[11]]$Significance == "No Significant") ,],size = 1) + 
  geom_point(data= df_uq_ruv2[[11]][df_uq_ruv2[[11]]$Significance == "SignificantPos",], size = 1) + 
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none") +
  scale_color_manual(values=c("Significant" = "black", "No Significant" = "grey", "SignificantPos" = "red")) +
  ggtitle(paste(label[11], sep = "")) + annotate("text", x = -4, y = 1.5, label = n_degs[[11]]) +
  annotate("text", x = -4, y = 0.5, label = n_de_posctrl[[11]], color = "red")

p12 <- ggplot(data = df_uq_ruv2[[11]], aes(x = logFC, y = -log10(PValue), col = Significance)) +
  xlim(-5, 5) + theme_classic(base_size = 7) + theme(legend.position = "none")

patch <- (p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9) / (p10 | p11 | plot_spacer())
ggsave("Peixoto_Figure5_Supplement1_part2.jpg",  width = 21, height = 30, units = "cm")
ggsave("Peixoto_Figure5_Supplement1_part2.pdf",  width = 21, height = 30, units = "cm")

# Histogram of p-value distribution #####
label <- unique(pb$azimuth_labels)
file_name <- paste("Peixoto_Figure5_Supplement1_part3.pdf", sep = "")
# jpeg(file = file_name, width = 794, height = 1134, units = "px")
pdf(file = file_name, width = 8.27, height = 11.81)
par(mfrow = c(4, 3))
for (i in seq_along(label)) {
  hist(df_uq_ruv2[[i]]$PValue, xlab = "p-value", main = paste(label[i], sep = " "))
}
dev.off()
# Negative binomial generalized model (GLM) on Pseudo-bulk sum ####
# Pseudo-bulk was aggregate across samples
pb_sum <- aggregateAcrossCells(pb, use.assay.type = "counts",
                               id = DataFrame(sample = pb$sample_id))
colnames(pb_sum) <- paste(substr(pb_sum$sample_id, 1, 1), pb_sum$condition, sep = "")

# UQ normalization
sum_counts <- as.matrix(counts(pb_sum))
uq <- betweenLaneNormalization(sum_counts, which="upper")

# RUVs fitting
groups <- makeGroups(pb_sum$condition)

ruved_data_sum <- RUVs(uq, cIdx = neg_ctrl, scIdx = groups, k = 2)
ruv2_sum <- ruved_data_sum$W

# DEA
y <- DGEList(counts(pb_sum), samples = colData(pb_sum))
keep <- filterByExpr(y, group = pb_sum$condition)
y <- y[keep, ]
y <- calcNormFactors(y, method = "upperquartile")
# Design Matrix
design <- model.matrix(~ 0 + y$samples$condition + ruv2_sum, y$samples)
colnames(design) <- c("HC", "SD", "W_1", "W_2")

y <- estimateDisp(y, design)
fit <- glmFit(y, design)
# Contrast creation (SD vs HC)
contrast <- makeContrasts(SD - HC, levels = design)

res_sum <- glmLRT(fit, contrast = contrast)

df_sum <- as.data.frame(res_sum$table)
df_sum <- df_sum[order(df_sum$PValue), ]

FDR <- as.data.frame(topTags(res_sum, n = length(rownames(pb_sum)))[, 5])

df_sum <- cbind(df_sum, FDR)
df_sum <- df_sum[order(df_sum$FDR), ]

save(df_sum, file = "snrna_DEA_sum.RData")

# Total number of genes
n_gene_tot <- length(rownames(df_sum))

df_degs_sum <- df_sum[df_sum$FDR < 0.05, ]
# N. of up-regulated DEGs
n_up_de_sum <- length(rownames(df_degs_sum[df_degs_sum$logFC >= 0, ]))
# N. of down-regulated DEGs
n_down_de_sum <-  length(rownames(df_degs_sum[df_degs_sum$logFC < 0, ]))
# N. of meta-analysis DE positive control genes
n_de_posctrl_sum <- length(intersect(rownames(df_degs_sum), posctrl$Gene_ID))
n.nuclei <- 47649
