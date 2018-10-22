# Documentation and definitions for data and constants

#### Sysdata ####

# Human IG germlines location
# GENE.loc

# Human IG germlines pseudo genes
#  PSEUDO

#### Data ####

#' Example of a IGH human naive b-cell data
#' 
#' Example IGH human naive b-cell data from a single
#' individual (see \cite{XXX}).
#'
#' @name sample_db
#' @docType data
#' @format A \code{data.frame} in Change-O format (\link{https://changeo.readthedocs.io/en/version-0.4.1---airr-standards/standard.html}) containing the following columns:
#' \itemize{
#'   \item \code{"SUBJECT"}: subject names
#'   \item \code{"V_CALL"}: V allele call(s) (in an IMGT format)
#'   \item \code{"D_CALL"}: D allele call(s) (in an IMGT format, only for heavy chains)
#'   \item \code{"J_CALL"}: J allele call(s) (in an IMGT format)
#' }
#' 
#' 
#' @references \cite{XXX}.
#' @keywords data antibody AIRR NGS
"sample_db"

#' Example of IGH human naive b-cell data from multiple
#' individuals
#'
#' Example IGH human naive b-cell data from multiple
#' individuals (see \cite{XXX}).
#'
#' @name samples_db
#' @docType data
#' @format A \code{data.frame} in Change-O format (\link{https://changeo.readthedocs.io/en/version-0.4.1---airr-standards/standard.html}) containing the following columns:
#' \itemize{
#'   \item \code{"SUBJECT"}: subject names
#'   \item \code{"V_CALL"}: V allele call(s) (in an IMGT format)
#'   \item \code{"D_CALL"}: D allele call(s) (in an IMGT format, only for heavy chains)
#'   \item \code{"J_CALL"}: J allele call(s) (in an IMGT format)
#' }
#' 
#' 
#' @references \cite{XXX}.
#' @keywords data antibody AIRR NGS
"samples_db"

#' Human IGHV germlines
#'
#' A \code{character} vector of all 342 human IGHV germline gene segment alleles
#' in IMGT Gene-db release 201408-4.
#'
#' @name HVGERM
#' @docType data
#' @format Values correspond to IMGT-gaped nuceltoide sequences (with
#' nucleotides capitalized and gaps represented by ".").
#' 
#' @references Xochelli \emph{et al}. (2014) Immunoglobulin heavy variable
#' (IGHV) genes and alleles: new entities, new names and implications for
#' research and prognostication in chronic lymphocytic leukaemia.
#' \emph{Immunogenetics}. 67(1):61-6.
#' @keywords data
"HVGERM"

#' Human IGHD germlines
#'
#' A \code{character} vector of all 37 human IGHD germline gene segment alleles
#' in IMGT Gene-db release 201408-4.
#'
#' @name HDGERM
#' @docType data
#' @format Values correspond to IMGT nuceltoide sequences.
#' 
#' @references Xochelli \emph{et al}. (2014) Immunoglobulin heavy variable
#' (IGHV) genes and alleles: new entities, new names and implications for
#' research and prognostication in chronic lymphocytic leukaemia.
#' \emph{Immunogenetics}. 67(1):61-6.
#' @keywords data
"HDGERM"

#' Human IGHJ germlines
#'
#' A \code{character} vector of all 13 human IGHJ germline gene segment alleles
#' in IMGT Gene-db release 201408-4.
#'
#' @name HJGERM
#' @docType data
#' @format Values correspond to IMGT nuceltoide sequences.
#' 
#' @references Xochelli \emph{et al}. (2014) Immunoglobulin heavy variable
#' (IGHV) genes and alleles: new entities, new names and implications for
#' research and prognostication in chronic lymphocytic leukaemia.
#' \emph{Immunogenetics}. 67(1):61-6.
#' @keywords data
"HJGERM"

#### Constants ####
## Defualt object