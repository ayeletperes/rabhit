#' Example human naive b-cell data
#'
#' Example VDJ-rearranged immunoglobulin naive b-cell sequences annotation derived from a single
#' individual, sequenced on the HiSeq platform.
#'
#' @name sample_db
#' @docType data
#' @format A \code{data.frame} where rows correspond to unique VDJ sequences and
#' columns include:
#' \itemize{
#'   \item \code{"SUBJECT"}   containing the subject name
#'   \item \code{"V_CALL"}    containing the IMGT/V-QUEST V allele call(s)
#'   \item \code{"D_CALL"}    containing the IMGT/V-QUEST V allele call(s), for heavy chain
#'   \item \code{"J_CALL"}    containing the IMGT/V-QUEST J allele call(s)
#' }
#' 
#' @references our paper.
#' @keywords data
"sample_db"

#' Example of multiple human naive b-cell data
#'
#' Example VDJ-rearranged immunoglobulin naive b-cell sequences annotation derived from multiple
#' individuals, sequenced on the HiSeq platform.
#'
#' @name samples_db
#' @docType data
#' @format A \code{data.frame} where rows correspond to unique VDJ sequences and
#' columns include:
#' \itemize{
#'   \item \code{"SUBJECT"}   containing the subject name
#'   \item \code{"V_CALL"}    containing the IMGT/V-QUEST V allele call(s)
#'   \item \code{"D_CALL"}    containing the IMGT/V-QUEST V allele call(s), for heavy chain
#'   \item \code{"J_CALL"}    containing the IMGT/V-QUEST J allele call(s)
#' }
#' 
#' @references our paper.
#' @keywords data
"samples_db"

#' Human IGHV germlines
#'
#' A \code{character} vector of all 342 human IGHV germline gene segment alleles
#' in IMGT Gene-db release 201408-4.
#'
#' @name VGERM
#' @docType data
#' @format Values correspond to IMGT-gaped nuceltoide sequences (with
#' nucleotides capitalized and gaps represented by ".").
#' 
#' @references Xochelli \emph{et al}. (2014) Immunoglobulin heavy variable
#' (IGHV) genes and alleles: new entities, new names and implications for
#' research and prognostication in chronic lymphocytic leukaemia.
#' \emph{Immunogenetics}. 67(1):61-6.
#' @keywords data
"VGERM"

#' Human IGHD germlines
#'
#' A \code{character} vector of all 37 human IGHD germline gene segment alleles
#' in IMGT Gene-db release 201408-4.
#'
#' @name DGERM
#' @docType data
#' @format Values correspond to IMGT nuceltoide sequences.
#' 
#' @references Xochelli \emph{et al}. (2014) Immunoglobulin heavy variable
#' (IGHV) genes and alleles: new entities, new names and implications for
#' research and prognostication in chronic lymphocytic leukaemia.
#' \emph{Immunogenetics}. 67(1):61-6.
#' @keywords data
"DGERM"

#' Human IGHJ germlines
#'
#' A \code{character} vector of all 13 human IGHJ germline gene segment alleles
#' in IMGT Gene-db release 201408-4.
#'
#' @name JGERM
#' @docType data
#' @format Values correspond to IMGT nuceltoide sequences.
#' 
#' @references Xochelli \emph{et al}. (2014) Immunoglobulin heavy variable
#' (IGHV) genes and alleles: new entities, new names and implications for
#' research and prognostication in chronic lymphocytic leukaemia.
#' \emph{Immunogenetics}. 67(1):61-6.
#' @keywords data
"JGERM"

# Human IG germlines location
# @source Corey T watson + IMGT.
# @format A \code{character} list of all three human IG chains, each cahin is a vector of all germline gene segment location.
# @examples
# \dontrun{
#  GENE.loc
# }
NULL

# Human IG germlines pseudo genes
# @source Corey T watson + IMGT.
# @format A \code{character} list of all three human IG chains, each cahin is a vector of all pseudo germline gene segment.
# @examples
# \dontrun{
#  PSEUDO
# }
NULL

