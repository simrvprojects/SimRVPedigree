#' Example pedigrees
#'
#' A dataset containing four example pedigrees.
#'
#' @docType data
#'
#' @usage data(exp_peds)
#'
#' @format A data frame with 50 rows and 14 variables:
#' \describe{
#' \item{FamID}{Family identification number}
#' \item{ID}{Individual identification number}
#' \item{gender}{Gender identification variable: `gender = 0` for males, and `gender = 1` females. }
#' \item{dad_id}{Identification number of father}
#' \item{mom_id}{Identification number of mother}
#' \item{affected}{Affection status: `affected = 1` if individual has developed lymphoid cancer, and 0 otherwise.}
#' \item{DA1}{Paternally inherited allele at the assumed disease locus: `DA1 = 1` if rare variant is present, and 0 otherwise.}
#' \item{DA2}{Maternally inherited allele at the assumed disease locus: `DA2 = 1` if rare variant is present, and 0 otherwise.}
#' \item{birth_year}{Subject's year of birth}
#' \item{onset_year}{Subject's year of disease onset, when applicable.}
#' \item{death_year}{Subject's year of death, when applicable.}
#' \item{RR}{Subject's relative risk of disease}
#' \item{available}{Availability status, 0 = unavailable, 1 = available}
#' \item{Gen}{The subject's generation number relative to the founder who introduced the rare variant.  That is, the founder who introduced the rare variant will have `Gen = 1`, his or her offspring will have `Gen = 2`, etc.}
#' }
#'
"exp_peds"
