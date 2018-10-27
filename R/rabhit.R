# rabhit package documentation and import directives

#' The RAbHIT package
#'
#' The \code{rabhit} package provides a robust novel method for determining
#' antibody heavy and light chain haplotypes by adapting a Bayesian framework.
#' The key functions in \code{rabhit}, broken down topic, are
#' described below.
#'
#'
#' @section  Haplotype and deletion infrence:
#' \code{rabhit} provides tools to infer haplotypes based on given anchor genes,
#' deletion detection based on relative gene usage, pooling v genes, and a single anchor gene.
#'
#' \itemize{
#'   \item  \link{createFullHaplotype}:      Infer haplotypes and single chromosome deletions based on an anchor gene.
#'   \item  \link{deletionsByVpooled}:       Single chromosomal deletion detection by pooling V genes.
#'   \item  \link{deletionsByBinom}:       Single chromosomal deletion detection by pooling V genes.
#' }
#'
#' @section  Haplotype visualization:
#' Functions for visualization of the inferred haplotypes and deletions
#'
#' \itemize{
#'   \item  \link{plotHaplotype}:            Haplotype inference plot.
#'   \item  \link{deletionHeatmap}:          Single chromosome deletion heatmap.
#'   \item  \link{hapHeatmap}:               Alleles per genes by chromosome heatmap.
#'   \item  \link{plotDeletionsByVpooled}:   V pooled based single chromosome deletion heatmap.
#'   \item  \link{plotDeletionsByBinom}:     Double chromosome deletion heatmap.
#' }
#'
#' @name     rabhit
#' @docType  package
#' @references
#' \enumerate{
#'   \item  our paper
#'  }
#'
#' @import   ggplot2
#' @import   graphics
#' @import   methods
#' @import   utils
#' @importFrom  cowplot     get_legend plot_grid ggdraw draw_label background_grid
#' @importFrom  plotly      ggplotly subplot
#' @importFrom  dplyr       do n desc funs %>%
#'                          as_data_frame data_frame data_frame_
#'                          bind_cols bind_rows combine rowwise slice
#'                          filter filter_ select select_ arrange arrange_
#'                          group_by group_by_ ungroup
#'                          mutate mutate_ summarize summarize_
#'                          mutate_at summarize_at count_
#'                          rename rename_ transmute transmute_ pull
#' @importFrom  reshape2    melt
#' @importFrom  mltools     bin_data
#' @importFrom  gtools      ddirichlet
NULL
