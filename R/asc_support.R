# Allele Similarity Cluster (ASC) support ---------------------------------

#' @include rabhit.R
NULL

##############################################################################################################
#' Prepare AIRR-seq data and germline for haplotype inference with allele similarity clusters (ASC)
#'
#' \code{convertToASC} converts IMGT-named V alleles to PIgLET Allele Similarity
#' Cluster (ASC) names and returns the pieces needed to run \code{\link{createFullHaplotype}}
#' on the clustered data. ASCs collapse near-identical and duplicated IGHV alleles
#' (for example \emph{IGHV1-69} / \emph{IGHV1-69D}) that are otherwise difficult to
#' assign unambiguously, and provide a naming scheme compatible with OGRDB reference
#' sets. It is a thin wrapper around the \pkg{piglet} package.
#'
#' @param  clip_db                AIRR-seq repertoire \code{data.frame}.
#' @param  germline               named vector of IMGT-gapped germline sequences for the
#'                                V genes (e.g. \code{HVGERM}).
#' @param  chain                  the chain. One of \code{"IGH"}, \code{"IGK"}, \code{"IGL"} or \code{"TRB"}.
#' @param  allele_cluster_table   optional precomputed ASC table (a \code{data.frame} with
#'                                \code{iuis_allele} and \code{new_allele} columns, as produced by
#'                                \code{piglet::inferAlleleClusters()} or downloaded from OGRDB with
#'                                \code{piglet::recentAlleleClusters()} / \code{piglet::extractASCTable()}).
#'                                If \code{NULL} it is inferred from \code{germline}.
#' @param  v_call                 the column holding the V allele calls to convert. Default \code{"v_call"}.
#' @param  ...                    further arguments passed to \code{piglet::inferAlleleClusters()}
#'                                when \code{allele_cluster_table} is inferred.
#'
#' @return a \code{list} with:
#' \itemize{
#'   \item \code{clip_db} - the input data with \code{v_call} converted to ASC names.
#'   \item \code{germline} - the ASC-named germline vector to pass as \code{toHap_GERM}.
#'   \item \code{genes_order} - an ASC gene order (derived from the chromosomal order in
#'         \code{GENE.loc}) to pass to the plotting functions.
#'   \item \code{allele_cluster_table} - the ASC table that was used.
#' }
#'
#' @details
#' The returned objects feed directly into the standard workflow:
#' \preformatted{
#' asc <- convertToASC(clip_db, HVGERM, chain = "IGH")
#' hap <- createFullHaplotype(asc$clip_db, toHap_col = "v_call",
#'                            hapBy_col = "j_call", hapBy = "IGHJ6",
#'                            toHap_GERM = asc$germline, chain = "IGH")
#' plotHaplotype(hap, genes_order = asc$genes_order)
#' }
#' Allele similarity clustering is applied to the V genes only; D and J calls are left unchanged.
#'
#' @examples
#' \dontrun{
#' data(samples_db, HVGERM)
#' clip_db <- samples_db[samples_db$subject == "I5", ]
#' asc <- convertToASC(clip_db, HVGERM, chain = "IGH")
#' hap <- createFullHaplotype(asc$clip_db, toHap_col = "v_call",
#'                            hapBy_col = "j_call", hapBy = "IGHJ6",
#'                            toHap_GERM = asc$germline, chain = "IGH")
#' plotHaplotype(hap, genes_order = asc$genes_order)
#' }
#'
#' @seealso \code{\link{createFullHaplotype}}; the \pkg{piglet} package for the underlying
#'   clustering (\code{inferAlleleClusters}, \code{assignAlleleClusters}, \code{germlineASC}).
#' @export
convertToASC <- function(clip_db,
                         germline,
                         chain = c("IGH", "IGK", "IGL", "TRB"),
                         allele_cluster_table = NULL,
                         v_call = "v_call",
                         ...) {
  if (!requireNamespace("piglet", quietly = TRUE)) {
    stop("convertToASC() requires the 'piglet' package. ",
         "Install it with install.packages('piglet').")
  }
  if (missing(chain)) chain <- "IGH"
  chain <- match.arg(chain)

  if (!v_call %in% names(clip_db)) {
    stop("Column '", v_call, "' not found in clip_db.")
  }

  # Obtain an ASC table, inferring one from the germline if none was supplied.
  if (is.null(allele_cluster_table)) {
    germline_cluster <- piglet::inferAlleleClusters(germline,
                                                    locus = paste0(chain, "V"),
                                                    ...)
    allele_cluster_table <- germline_cluster$alleleClusterTable
  }

  # Convert the V calls and the germline to ASC names.
  clip_db <- piglet::assignAlleleClusters(clip_db, allele_cluster_table, v_call = v_call)
  germline_asc <- piglet::germlineASC(allele_cluster_table, germline)
  germline_asc <- germline_asc[!is.na(names(germline_asc)) & names(germline_asc) != "NA"]

  # Derive an ASC gene order from the bundled chromosomal IMGT order.
  imgt_order <- GENE.loc[[chain]]
  imgt2asc <- stats::setNames(gsub("\\*.*", "", allele_cluster_table$new_allele),
                              gsub("\\*.*", "", allele_cluster_table$iuis_allele))
  ordered_asc <- unique(unname(imgt2asc[imgt_order[imgt_order %in% names(imgt2asc)]]))
  ordered_asc <- ordered_asc[!is.na(ordered_asc)]
  asc_genes <- unique(gsub("\\*.*", "", names(germline_asc)))
  genes_order <- c(ordered_asc, setdiff(asc_genes, ordered_asc))

  list(clip_db = clip_db,
       germline = germline_asc,
       genes_order = genes_order,
       allele_cluster_table = allele_cluster_table)
}
