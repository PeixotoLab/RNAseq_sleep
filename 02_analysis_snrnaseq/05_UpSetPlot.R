# Upset plot ####
# The goal of this analysis is to identify the unique DE genes 
# for each cell-type, using the upset plot

library(patchwork)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(readxl)
library(UpSetR)
library(ComplexUpset)
library(dplyr)

# Load DEA results
df_uq_ruv2 <- readRDS("snRNA_UQ_RUVs2_DFResults.rds")
names(df_uq_ruv2) <- unique(pb$azimuth_labels)

# DEGs were selected
df_degs <- lapply(df_uq_ruv2, function(u) rownames(u[u$FDR < 0.05, ]))

# Since the p-value distribution wasn't good, the L6b CTX label was removed
ll <- list("L2/3 IT CTX" = df_degs[[1]], "L4/5 IT CTX" = df_degs[[2]],
           "L5 IT CTX" = df_degs[[3]], "L5 PT CTX" = df_degs[[4]],
           "L6 CT CTX" = df_degs[[6]], "L6 IT CTX" = df_degs[[7]],
           "L6b CTX" = df_degs[[8]],
           "Pvalb" = df_degs[[9]], "Sst" = df_degs[[10]], "Vip" = df_degs[[11]])

df.edgeRList <- fromList(ll)

p1 <- upset(df.edgeRList, colnames(df.edgeRList)[1:10],
            sort_intersections_by = c("degree", "cardinality"),
            sort_intersections = "ascending",
            name = "", width_ratio = 0.1, keep_empty_groups = TRUE,
            queries = list(
              upset_query(set = "L2/3 IT CTX", fill = "#0BE652"),
              upset_query(set = "L4/5 IT CTX", fill = "#00E5E5"),
              upset_query(set = "L5 IT CTX", fill = "#50B2AD"),
              upset_query(set = "L5 PT CTX", fill = "#0D5B78"),
              upset_query(set = "L6 CT CTX", fill = "#2D8CB8"),
              upset_query(set = "L6 IT CTX", fill = "#A19922"),
              upset_query(set = "L6b CTX", fill = "#7044AA"),
              upset_query(set = "Pvalb", fill = "#D93137"),
              upset_query(set = "Sst", fill = "#FF9900"),
              upset_query(set = "Vip", fill = "#B864CC")
            ),
            intersections = list(
              # Unique DEGs were visualized for each cell-type
              "L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX", "L5 PT CTX",
              "L6 CT CTX", "L6 IT CTX", "L6b CTX", "Pvalb", "Sst", "Vip"
              # Intersection of all cell-types was visualized
              # c("L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX", "L5 PT CTX",
              #  "L5/6 NP CTX", "L6 CT CTX", "L6 IT CTX", "Pvalb", "Sst", "Vip")
            ),
            base_annotations = list(
              "Intersection size" = (
                intersection_size(
                  bar_number_threshold = 1, # show all numbers on top of bars
                  width = 0.4,  # reduce width of the bars
                  text = list(size = 2)
                )
                # add some space on the top of the bars
                + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
                + theme(
                  # hide grid lines
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  # show axis lines
                  axis.line = element_line(colour = "black")
                )
              )
            ),
            stripes = upset_stripes(
              geom = geom_segment(size = 3),
              colors = c("grey95", "white")
            ),
            matrix = intersection_matrix(
              geom = geom_point(
                shape = "circle filled",
                size = 2,
                stroke = 0
              )
            ),
            set_sizes = (
              upset_set_size(geom = geom_bar(width = 0.5), filter_intersections = TRUE)
              + theme(
                axis.line.x = element_line(colour = "black"),
                axis.ticks.x = element_line()
              )
            ),
            themes = upset_default_themes(text = element_text(size = 7, color = "black"))) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

ggsave("Peixoto_Figure6_part1.jpg",  width = 18, height = 9, units = "cm")
ggsave("Peixoto_Figure6_part1.pdf",  width = 18, height = 9, units = "cm")

# Code for extracting the list of the unique differentially expressed genes
# and the intersections inside upset.
df2 <- data.frame(gene = unique(unlist(ll)))

df1 <- lapply(ll, function(x) {
  data.frame(gene = x)
}) %>%
  bind_rows(.id = "path")

df_int <- lapply(df2$gene, function(x) {
  # pull the name of the intersections
  intersection <- df1 %>%
    dplyr::filter(gene == x) %>%
    arrange(path) %>%
    pull("path") %>%
    paste0(collapse = "|")
  # build the dataframe
  data.frame(gene = x, int = intersection)
}) %>%
  bind_rows()

# Number of genes for each path
df_int %>%
  group_by(int) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

unique_ll <- list()
for (i in seq_along(ll)) {
  unique_ll[[i]] <- df_int[df_int$int==names(ll)[i],]
  unique_ll[[i]] <- unique(unique_ll[[i]]$gene)
  unique_ll[[i]] <- as.data.frame(unique_ll[[i]])
  colnames(unique_ll[[i]]) <- "Gene_Stable_ID"
}

unique_inter <- df_int[df_int$int=="L2/3 IT CTX|L4/5 IT CTX|L5 IT CTX|L5 PT CTX|L6 CT CTX|L6 IT CTX|L6b CTX|Pvalb|Sst|Vip",]
unique_inter <- unique(unique_inter$gene)
unique_inter <- as.data.frame(unique_inter)
colnames(unique_inter) <- "Gene_Stable_ID"

unique_ll[[11]] <- unique_inter

names(unique_ll) <- c(gsub("/","", names(ll)), "IntersectionOfAll")

writexl::write_xlsx(unique_ll, "Peixoto_Figure6_Supplement1.xlsx")

# Save combined list of DEGs for Gluta and GABA #####
degs_gluta <- c(df_degs[[1]], df_degs[[2]], df_degs[[3]], df_degs[[4]], 
                df_degs[[6]], df_degs[[7]], df_degs[[8]])
degs_gluta <- unique(degs_gluta)
degs_gluta <- as.data.frame(degs_gluta)
colnames(degs_gluta) <- "Gene_Stable_ID"

degs_gaba <- c(df_degs[[9]], df_degs[[10]], df_degs[[11]])
degs_gaba <- unique(degs_gaba)
degs_gaba <- as.data.frame(degs_gaba)
colnames(degs_gaba) <- "Gene_Stable_ID"

# Background ####
# Background for Glutamatergic 
bg_gluta <- c(rownames(df_uq_ruv2[[1]]), rownames(df_uq_ruv2[[2]]), rownames(df_uq_ruv2[[3]]),
              rownames(df_uq_ruv2[[4]]), rownames(df_uq_ruv2[[6]]),
              rownames(df_uq_ruv2[[7]]), rownames(df_uq_ruv2[[8]]))
bg_gluta <- unique(bg_gluta)
bg_gluta <- as.data.frame(bg_gluta)
colnames(bg_gluta) <- "Gene_Stable_ID"

# Background for GABAergic
bg_gaba <- c(rownames(df_uq_ruv2[[9]]), rownames(df_uq_ruv2[[10]]), rownames(df_uq_ruv2[[11]]))
bg_gaba <- unique(bg_gaba)
bg_gaba <- as.data.frame(bg_gaba)
colnames(bg_gaba) <- "Gene_Stable_ID"

comb_ll <- list(
  "Glutamatergic" = degs_gluta, "Background_Gluta" = bg_gluta, 
  "GABAergic" = degs_gaba, "Background_GABA" = bg_gaba
)
writexl::write_xlsx(comb_ll, "Peixoto_Figure6_background.xlsx")
