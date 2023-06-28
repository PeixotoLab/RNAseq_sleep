pbMDS <- function(x) {
  # check validity of input pseudobulk-SCE
  # (cells should have been aggregated by cluster-sample)
  .check_pbs(x, check_by = TRUE)
  
  y <- as.list(assays(x))
  y <- do.call("cbind", y)
  y <- y[, (j <- c(t(.n_cells(x))) != 0)]
  d <- DGEList(unname(y), remove.zeros = TRUE)
  d <- calcNormFactors(d)
  
  mds <- plotMDS.DGEList(d, plot = FALSE)
  nk <- length(kids <- assayNames(x))
  
  ss <- rep(colnames(x), nk)
  ks <- rep(kids, each = ncol(x))
  
  if (any(!j)) {
    txt <- paste(sQuote(ks[!j]), sQuote(ss[!j]), sep = "-")
    message("Removing cluster-sample instance(s) ", 
            paste(txt, collapse = ", "))
  }
  
  df <- data.frame(
    MDS1 = mds$x, MDS2 = mds$y, 
    cluster_id = factor(ks[j], levels = kids), 
    group_id = rep(x$group_id, nk)[j])
  
  cols <- .cluster_colors
  if (nk > length(cols)) 
    cols <- colorRampPalette(cols)(nk)
  
  ggplot(df, aes_string(x = "MDS1", y = "MDS2", 
                        col = "cluster_id", shape = "group_id")) +
    scale_color_manual(values = cols) +
    geom_point(size = 3, alpha = 0.8) + 
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme_bw() + theme(aspect.ratio = 1,
                       axis.text = element_text(color = "black"),
                       panel.grid.minor = element_blank())#,
  # panel.grid.major = element_line(linewidth = 0.2, color = "lightgrey"))
}

.n_cells <- function(x) {
  y <- int_colData(x)$n_cells
  if (is.null(y)) return(NULL)
  if (length(metadata(x)$agg_pars$by) == 2)
    y <- as.matrix(data.frame(y, check.names = FALSE))
  return(as.table(y))
}

.cluster_colors <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

.check_pbs <- function(pbs, sce = NULL, check_by = TRUE) {
  stopifnot(is(pbs, "SingleCellExperiment"),
            !is.null(ei <- metadata(pbs)$experiment_info),
            !is.null(agg_pars <- metadata(pbs)$agg_pars),
            !is.null(n_cells <- .n_cells(pbs)),
            identical(assayNames(pbs), rownames(n_cells)),
            identical(colnames(pbs), colnames(n_cells)))
  if (!is.null(sce)) {
    stopifnot(identical(ei, metadata(sce)$experiment_info),
              identical(assayNames(pbs), levels(sce[[agg_pars$by[1]]])),
              identical(rownames(pbs), rownames(sce)))
    if (length(agg_pars$by == 2))
      stopifnot(identical(colnames(pbs), levels(sce[[agg_pars$by[2]]])))
  }
  if (check_by)
    stopifnot(!is.null(pbs[["group_id"]]),
              identical(agg_pars$by, c("cluster_id", "sample_id")))
}