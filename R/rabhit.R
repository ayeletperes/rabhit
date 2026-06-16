# rabhit package documentation and import directives
"_PACKAGE"
#' The RAbHIT package
#'
#' The \code{rabhit} package provides a robust novel method for determining
#' antibody heavy and light chain haplotypes by adapting a Bayesian framework.
#' The key functions in \code{rabhit}, broken down by topic, are
#' described below.
#'
#'
#' @section  Haplotype and deletions inference:
#' \code{rabhit} provides tools to infer haplotypes based on given anchor genes,
#' deletion detection based on relative gene usage, pooling v genes, and a single anchor gene.
#'
#' \itemize{
#'   \item  \link{createFullHaplotype}:      Haplotypes inference and single chromosome deletions based on an anchor gene.
#'   \item  \link{deletionsByVpooled}:       Single chromosomal deletion detection by pooling V genes.
#'   \item  \link{deletionsByBinom}:         Double chromosomal deletion detection by relative gene usage.
#'   \item  \link{geneUsage}:                Relative gene usage.
#'   \item  \link{nonReliableVGenes}:        Non reliable gene assignment detection.
#' }
#'
#' @section  Haplotype and deletions visualization:
#' Functions for visualization of the inferred haplotypes and deletions
#'
#' \itemize{
#'   \item  \link{plotHaplotype}:            Haplotype inference map.
#'   \item  \link{deletionHeatmap}:          Single chromosome deletions heatmap.
#'   \item  \link{hapHeatmap}:               Chromosome comparison of multiple samples.
#'   \item  \link{hapDendo}:                 Hierarchical clustering of multiple haplotypes based on Jaccard distance.
#'   \item  \link{plotDeletionsByVpooled}:   V pooled based single chromosome deletions heatmap.
#'   \item  \link{plotDeletionsByBinom}:     Double chromosome deletions heatmap.
#' }
#'
#' @name     rabhit
#' @references
#' \enumerate{
#'   \item  Gidoni, M., Snir, O., Peres, A., Polak, P., Lindeman, I., Mikocziova, I., . . . Yaari, G. (2019).
#'   Mosaic deletion patterns of the human antibody heavy chain gene locus shown by Bayesian haplotyping.
#'   Nature Communications, 10(1). doi:10.1038/s41467-019-08489-3
#'  }
#'
#' @import   ggplot2
#' @import   utils
#' @importFrom  graphics         grid image axis points text par plot
#' @importFrom  cowplot          get_legend plot_grid ggdraw draw_label background_grid
#' @importFrom  dplyr            do n desc %>% distinct all_of across .data
#'                               bind_cols bind_rows rowwise slice
#'                               filter select arrange
#'                               group_by ungroup
#'                               mutate summarize
#'                               mutate_at summarize_at count
#'                               rename transmute pull ungroup row_number if_else
#' @importFrom  stats            hclust as.dendrogram as.dist binom.test p.adjust setNames weighted.mean
#' @importFrom  gtable           gtable_filter gtable_add_rows gtable_add_grob
#' @importFrom  grDevices        dev.off pdf recordPlot dev.control
#' @importFrom  tidyr            separate_rows separate pivot_longer expand_grid
#' @importFrom  grid             gpar textGrob
NULL

