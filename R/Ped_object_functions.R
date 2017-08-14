#' Create an object of class ped.
#'
#' Create a ped object, required input for \code{\link{reassignGen.ped}}, \code{\link{censor.ped}}, \code{\link{affInfo.ped}}, and \code{\link{trim.ped}} functions.
#'
#' The data frame supplied to \code{new.ped}, \code{ped_file}, \emph{must} contain the following columns:
#' \tabular{lll}{
#' \strong{name} \tab \strong{type} \tab \strong{description} \cr
#' \code{FamID} \tab numeric \tab family identification number \cr
#' \code{ID} \tab numeric \tab individual identification number \cr
#' \code{dadID} \tab numeric \tab identification number of father \cr
#' \code{momID} \tab numeric \tab identification number of mother \cr
#' \code{sex} \tab numeric \tab gender identification; if male \code{sex = 0}, if female \code{sex = 1} \cr
#' \code{affected  } \tab numeric \tab disease-affection status, if affected by disease \code{affected  = 1}, otherwise, \code{affected = 0} \cr
#' }
#'
#' Optionally, \code{ped_file} \emph{may} contain any of the following columns:
#' \tabular{lll}{
#' \strong{name} \tab \strong{type}\tab \strong{description} \cr
#' \code{available} \tab numeric \tab availibility status; \code{available} = 1 if available, and 0 otherwise.  If missing all pedigree members are assumed to be available \cr
#' \code{DA1} \tab numeric \tab paternally inherited allele at the assumed disease locus: \code{DA1} = 1 if rare variant is present, and 0 otherwise\cr
#' \code{DA2} \tab numeric \tab maternally inherited allele at the assumed disease locus: \code{DA2} = 1 if rare variant is present, and 0 otherwise\cr
#' \code{birthYr} \tab numeric \tab the individual's birth year \cr
#' \code{onsetYr} \tab numeric \tab the individual's year of disease onset, when applicable \cr
#' \code{deathYr} \tab numeric \tab the individual's year of death, when applicable \cr
#' \code{RR} \tab numeric \tab the individual's relative-risk of disease \cr
#' \code{Gen} \tab numeric \tab the individual's generation number relative to the founder who introduced the rare variant.  That is, the founder who introduced the rare variant will have \code{Gen} = 1, his or her offspring will have \code{Gen} = 2, etc. \cr
#' \code{proband} \tab numeric \tab a proband identifier: \code{proband} = 1 if the individual is the proband, and 0 otherwise.\cr
#' }
#'
#' \emph{We note that some of the optional fields above may be required for different ped functions}
#'
#' @param ped_file Data.frame. A pedigree, see details.
#'
#' @return An object of class ped.
#' @export
#'
#' @examples
#' data(EgPeds)
#' head(EgPeds)
#'
#' ped1 = new.ped(EgPeds[EgPeds$FamID == 1, ])
#' head(ped1, n = 3)
#' class(ped1)
#' summary(ped1)
#'
#' AllPeds = new.ped(EgPeds)
#' head(AllPeds)
#' class(AllPeds)
#' summary(AllPeds)
#'
new.ped <- function(ped_file) {
  n <- length(unique(ped_file$FamID))

  if (n > 1){
    lapply(seq(1, n, by = 1), function(x){
      check_ped(ped_file[ped_file$FamID == unique(ped_file$FamID)[x], ])})
  } else {
    check_ped(ped_file)
  }

  if (!"available" %in% colnames(ped_file)) ped_file$available <- 1

  return(ped(ped_file))

}

#' Constructor function for an object of class ped
#'
#' @inheritParams new.ped
#'
#' @return an object of class \code{ped}.
#' @export
ped <- function(ped_file) {
  class(ped_file) <- c("ped", class(ped_file))
  return(ped_file)
}

#' Check to see if object is of class ped
#'
#' @param x An R object.
#'
#' @return Logical. Indicates if \code{x} is of class \code{ped}.
#' @keywords internal
#'
is.ped <- function(x) {
  return(inherits(x, "ped"))
}

#' @export
summary.ped <- function(object, ...) {
  n <- length(unique(object$FamID))

  famdat <- lapply(unique(object$FamID), function(x){
    sumVars(object[object$FamID == x, ])
  })

  return(do.call(rbind, famdat))
}
