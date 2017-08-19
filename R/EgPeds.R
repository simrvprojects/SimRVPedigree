#' Example pedigrees
#'
#' A dataset containing five example pedigrees.
#'
#' @docType data
#'
#' @usage data(EgPeds)
#'
#' @format A data frame with 65 rows and 14 variables:
#' \describe{
#' \item{FamID}{Family identification number}
#' \item{ID}{Individual identification number}
#' \item{sex}{Gender identification variable: \code{sex = 0} for males, and \code{sex = 1} females. }
#' \item{dadID}{Identification number of father}
#' \item{momID}{Identification number of mother}
#' \item{affected}{Affection status: \code{affected = TRUE} if individual has developed lymphoid cancer, and \code{FALSE} otherwise.}
#' \item{DA1}{Paternally inherited allele at the assumed disease locus: \code{DA1 = 1} if rare variant is present, and 0 otherwise.}
#' \item{DA2}{Maternally inherited allele at the assumed disease locus: \code{DA2 = 1} if rare variant is present, and 0 otherwise.}
#' \item{birthYr}{Year of birth}
#' \item{onsetYr}{Year of disease onset, when applicable.}
#' \item{deathYr}{Year of death, when applicable.}
#' \item{RR}{The subject's relative-risk of disease}
#' \item{available}{Availability status, \code{available  = TRUE} if the individual is unavailable, and \code{FALSE} 1 otherwise}
#' \item{Gen}{The subject's generation number relative to the founder who introduced the rare variant.  That is, the founder who introduced the rare variant will have \code{Gen = 1}, his or her offspring will have \code{Gen = 2}, etc.}
#' \item{proband}{Proband identification variable, \code{proband = TRUE} if the individual is the proband, and \code{FALSE} otherwise.}
#' }
#'
"EgPeds"
