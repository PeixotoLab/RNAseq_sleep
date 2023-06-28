library(UpSetR)
library(ComplexUpset)
library(patchwork)

DGE_Fishpond <- read.csv("Fishpond_Annotated_DGE.csv")
DGE_Fishpond <- DGE_Fishpond$Gene_Stable_ID

DTE_Fishpond <- read.csv("Fishpond_Annotated_DTE.csv")
DTE_Fishpond <- DTE_Fishpond$Gene_Stable_ID

PosCtrls <- read.table("Additional_File2_Positive_Controls.txt", header = T)

df.RUVs1 <- readRDS("snRNA_edgeR_DFResults_RUVs1.rds")
# select DEGs
df.DEGs <- lapply(df.RUVs1, function(u) rownames(u[u$FDR<0.05,])) 

# UpSet plot #####
# Select intersection between DGE from bulk and snRNA DGE
Fishpond.DGE <- lapply(df.DEGs, function(u) intersect(u, DGE_Fishpond))

# For Glutamatergic labels
ll <- list("L2/3 IT CTX"=Fishpond.DGE[[1]], "L4/5 IT CTX"=Fishpond.DGE[[2]], "L5 IT CTX"=Fishpond.DGE[[3]],"L5 PT CTX"=Fishpond.DGE[[4]],
           "L5/6 NP CTX"=Fishpond.DGE[[5]], "L6 CT CTX" =Fishpond.DGE[[6]], "L6 IT CTX" =Fishpond.DGE[[7]],"L6b CTX" =Fishpond.DGE[[8]])

df.edgeRList <- fromList(ll)

p1 <- upset(
  df.edgeRList, colnames(df.edgeRList)[1:8],
  min_size=10, sort_intersections_by = c('degree', 'cardinality'),sort_intersections="ascending",
  name="",
  keep_empty_groups=TRUE,
  queries=list(
    upset_query(set='L2/3 IT CTX', fill='#0BE652'),
    upset_query(set='L4/5 IT CTX', fill='#00E5E5'),
    upset_query(set='L5 IT CTX', fill='#50B2AD'),
    upset_query(set='L5 PT CTX', fill='#0D5B78'),
    upset_query(set='L5/6 NP CTX', fill='#3E9E64'),
    upset_query(set='L6 CT CTX', fill='#2D8CB8'),
    upset_query(set='L6 IT CTX', fill='#A19922'),
    upset_query(set='L6b CTX', fill='#7044AA')
  ),
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        bar_number_threshold=1,  # show all numbers on top of bars
        width=0.5,   # reduce width of the bars
        text = list(size=6)
      )
      # add some space on the top of the bars
      + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line=element_line(colour='black')
      )
    )
  ),
  stripes=upset_stripes(
    geom=geom_segment(size=12),  # make the stripes larger
    colors=c('grey95', 'white')
  ),
  # to prevent connectors from getting the colorured
  # use `fill` instead of `color`, together with `shape='circle filled'`
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=5,
      stroke=0
    )
  ),
  set_sizes=(
    upset_set_size(geom=geom_bar(width=0.4),filter_intersections=TRUE)
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line()
    )
  ),
  themes=upset_default_themes(text=element_text(size=25, color="black"))
)  

# For GABAergic
ll <- list("Pvalb"=Fishpond.DGE[[9]], "Sst"=Fishpond.DGE[[10]],"Vip"=Fishpond.DGE[[11]])

df.edgeRList <- fromList(ll)

p2 <- upset(
  df.edgeRList, colnames(df.edgeRList)[1:3],
  min_size=4, sort_intersections_by = c('degree', 'cardinality'),sort_intersections="ascending",
  name="",
  keep_empty_groups=TRUE,
  queries=list(
    upset_query(set='Pvalb', fill='#D93137'),
    upset_query(set='Sst', fill='#FF9900'),
    upset_query(set='Vip', fill='#B864CC')
  ),
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        bar_number_threshold=1,  # show all numbers on top of bars
        width=0.5,   # reduce width of the bars
        text = list(size=6)
      )
      # add some space on the top of the bars
      + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line=element_line(colour='black')
      )
    )
  ),
  stripes=upset_stripes(
    geom=geom_segment(size=12),  # make the stripes larger
    colors=c('grey95', 'white')
  ),
  # to prevent connectors from getting the colorured
  # use `fill` instead of `color`, together with `shape='circle filled'`
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=5,
      stroke=0
    )
  ),
  set_sizes=(
    upset_set_size(geom=geom_bar(width=0.4),filter_intersections=TRUE)
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line()
    )
  ),
  themes=upset_default_themes(text=element_text(size=25, color="black"))
)

ggsave(p1, file="snRNA_UpSetPlot_Glutamatergic.pdf",  width = 30, height = 30, units = "cm")
ggsave(p2, file="snRNA_UpSetPlot_GABAergic.pdf",  width = 30, height = 30, units = "cm")

# Intersection with DGE and DTE from bulk and Microarray #####
# Up
df.UpDEGs <- lapply(df.DEGs, function(u) u[u$logFC>=0,]) 
# Down
df.DownDEGs <- lapply(df.DEGs, function(u) u[u$logFC<0,])

# DEGs single-nuclear
Up.DEG <- lapply(df.UpDEGs, function(x) as.data.frame(rownames(x)))
Down.DEG <- lapply(df.DownDEGs, function(x) as.data.frame(rownames(x)))

# Intersection with DGE from bulk 
Up.DGEbulk <- lapply(df.UpDEGs, function(x) as.data.frame(intersect(rownames(x), DGE_Fishpond)))
Down.DGEbulk <- lapply(df.DownDEGs, function(x) as.data.frame(intersect(rownames(x), DGE_Fishpond)))

# Intersection with DTE from bulk 
Up.DTEbulk <- lapply(df.UpDEGs, function(x) as.data.frame(intersect(rownames(x), DTE_Fishpond)))
Down.DTEbulk <- lapply(df.DownDEGs, function(x) as.data.frame(intersect(rownames(x), DTE_Fishpond)))

# Intersection with Microarray
Up.Pos <- lapply(df.UpDEGs, function(x) as.data.frame(intersect(rownames(x), PosCtrls$Gene_ID)))
Down.Pos <- lapply(df.DownDEGs, function(x) as.data.frame(intersect(rownames(x), PosCtrls$Gene_ID)))
