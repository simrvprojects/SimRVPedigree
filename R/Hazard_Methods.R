#' Create an object of class hazard.
#'
#' Create a hazard object, required input for \code{\link{sim_RVped}}, \code{\link{sim_ped}}, and \code{\link{sim_life}} functions.
#'
#' Users are permitted to specify harzard objects for two scenarios: (1) for a disease without subtypes or  (2) for a disease with multiple subtypes.
#'
#'
#' When simulating a disease \emph{without subtypes}, \code{hazardDF} must contain 3 columns that meet the following criteria:
#' \describe{
#'  \item{column 1:}{age-specific hazard rates of \emph{disease} for the population of interest}
#'  \item{column 2:}{age-specific hazard rates of \emph{death} for the \strong{unaffected} population.  If the disease of interest is sufficiently rare, so that death by the disease is rare, the user may choose to use the population, age-specific, hazard rates of death instead.}
#'  \item{column 3:}{age-specific hazard rates of \emph{death} for the \strong{affected} population.}
#' }
#'
#' When simulating a disease \emph{with n disease subtypes}, \code{hazardDF} must contain n + 2 columns that meet the following criteria:
#' \describe{
#'  \item{column 1:}{age-specific hazard rates of \emph{disease} for the first subtype of interest}
#'  \item{column 2:}{age-specific hazard rates of \emph{disease} for the second subtype of interest}
#'  \item{...}{}
#'  \item{column n:}{age-specific hazard rates of \emph{disease} for the \eqn{n^{th}} subtype of interest}
#'  \item{column n + 1:}{age-specific hazard rates of \emph{death} for the \strong{unaffected} population.  If the disease of interest is sufficiently rare, so that death by the disease is rare, the user may choose to use the population, age-specific, hazard rates of death instead.}
#'  \item{column n + 2:}{age-specific hazard rates of \emph{death} for an individual affectd by any of the \eqn{n} subtypes.}
#' }

#' Users must provide \code{partition} in years; e.g. a hazard rate for a baby between 6 months and 1 year of age should have lower bound 0.5 years and an upper bound 1 year.  Additionally, \code{partition} must apply to all of the age-specific hazard rates in \code{hazardDF}.
#'
#' @param hazardDF Data.frame. A data.frame contianing the age-specific hazard rate(s) of \emph{disease} in the population of interest, the age-specific hazard rate for \emph{death} in the \strong{unaffected} population, and the age-specific hazard rate for \emph{death} in the \strong{affected} population.  See details.
#' @param partition Numeric vector. The partition of ages, in years, over which to apply the age-specific hazard rates in \code{hazardDF}. If not supplied,  defaults to a partition that starts at 0 and increases in yearly increments.  See details.
#' @param subtype_ID List. If specifying a disease with multiple subtypes, a list of character subtype IDs.  By default, \code{subtype_ID = NULL}, i.e. no subtypes to simulate.
#'
#' @return An object of class hazard.
#' @export
#'
#' @examples
#' # Specifying the hazard rates for a disease
#' # with only one subtype (or grouped subtypes).
#' data(AgeSpecific_Hazards)
#'
#' head(AgeSpecific_Hazards)
#' nrow(AgeSpecific_Hazards)
#'
#' my_HR <- hazard(hazardDF = AgeSpecific_Hazards)
#' class(my_HR)
#' head(my_HR[[1]])
#'
#' #NOTE: since partition was not supplied, the partition has been assummed to
#' # start at 0 and increase in yearly increments.
#' my_HR[[2]]
#'
#' my_HR
#'
#'
#' # Specifying the hazard rates for a with teo disease subtypes
#' data(SubtypeHazards)
#'
#' head(SubtypeHazards)
#' nrow(SubtypeHazards)
#'
#' my_SHR <- hazard(hazardDF = SubtypeHazards,
#'                  subtype_ID = c("Hodgkin Lymphoma", "Non-Hodgkin Lymphoma"))
#' class(my_SHR)
#' head(my_SHR[[1]])
#' my_SHR
#'
#'
#'
hazard <- function(hazardDF, partition = NULL,
                   subtype_ID = NULL) {
  # Set up a new object of class Hazard

  if(is.null(partition)){
    partition = seq(0, nrow(hazardDF), by = 1)
  } else {

    if (any(is.na(partition))) {
      stop('partition contains missing values')
    }

    if (any(duplicated(partition))) {
      stop('partition contains duplicated entries')
    }

    if (min(partition) != 0) {
      stop('parition does not start at birth (age 0)')
    } else if (max(partition) < 65){
      warning ('For optimal results please specify age-specific hazard rates that begin at birth and end near the life expectancy of the population to which the age-specific hazards apply.')
    }

    if (length(partition) == 1 |
        length(partition) != (nrow(hazardDF) + 1) ) {
      stop ('Incorrect partition length. Please ensure: length(partition) == nrow(hazardDF) + 1')
    }

  }

  if (any(is.na(hazardDF))) {
    stop('hazardDF contains missing values')
  }

  if (any(hazardDF < 0)) {
    stop('hazardDF contains negative values')
  }

  # Choose Hazard Format
  # Execute if hazard_type and subtype_ID are NULL
  # otherwise check to see that both have been properly specified
  # and create the multi hazard object

  if (is.null(subtype_ID)) {
    #print('Creating the hazard object for a disease with no subtypes')
    if (class(hazardDF) != "data.frame") {
      stop("hazardDF must be a data frame with 3 columns:
           column 1 = population age-specific hazard rate of disease,
           column 2 = age-specific hazard rate of death in the unaffected population,
           column 3 = age-specific hazard rate of death in the affected population")
    } else if (ncol(hazardDF) != 3) {
      stop("Perhaps you are trying to specify the hazard rates for a disease with multiple subtypes and forgot to supply a list of subtype labels to subtype_ID?\n
           Otherwise, for diseases with only one subtype hazardDF must be a data frame with 3 columns:
           column 1 = population age-specific hazard rate of disease,
           column 2 = age-specific hazard rate of death in the unaffected population,
           column 3 = age-specific hazard rate of death in the affected population")
    } else if(sum(hazardDF[, 2] > hazardDF[, 3]) > nrow(hazardDF)/2 |
              sum(which(hazardDF[, 2] == hazardDF[, 3])) != 0){
      warning("Please check that you have specified hazardDF such that:
              column 2 = age-specific hazard rate of death in the UNAFFECTED population,
              column 3 = age-specific hazard rate of death in the AFFECTED population")
    }

    #create dummy subtype IDs
    s_ID <- c("no_subtypes")
    } else {
      num_subs <- length(subtype_ID)
      s_ID <- subtype_ID
      #print(paste0('Creating the hazard object for a disease with ', num_subs,' subtypes'))

      if (class(hazardDF) != "data.frame") {
        stop(paste0("hazardDF must be a data frame with ", num_subs + 2, " columns:
                    columns 1:", num_subs," = population age-specific hazard rates of disease for each subtype in subtype_ID,
                    column ", num_subs + 1, " = age-specific hazard rate of death in the UNAFFECTED population,
                    column ", num_subs + 2, "= age-specific hazard rate of death in the AFFECTED population"))
      } else if (ncol(hazardDF) != num_subs + 2) {
        stop(paste0("hazardDF must be a data frame with ", num_subs + 2, " columns:
                    columns 1:", num_subs," = population age-specific hazard rates of disease for each subtype in subtype_ID,
                    column ", num_subs + 1, " = age-specific hazard rate of death in the UNAFFECTED population,
                    column ", num_subs + 2, "= age-specific hazard rate of death in the AFFECTED population"))
      } else if(sum(hazardDF[, num_subs + 1] > hazardDF[, num_subs + 2]) > nrow(hazardDF)/2 |
                sum(which(hazardDF[, num_subs + 1] == hazardDF[, num_subs + 2])) != 0){
        warning(paste0("Please check that you have specified hazardDF such that:
                       column ", num_subs + 1, " = age-specific hazard rate of death in the UNAFFECTED population,
                       column ", num_subs + 2, " = age-specific hazard rate of death in the AFFECTED population"))
      }
      }


  obj <- list(hazardDF = hazardDF,
              partition = partition,
              subtype_ID = s_ID)
  class(obj) <- "hazard"
  return(obj)
      }

#' Check to see if object is of class hazard
#'
#' @param x An R object.
#'
#' @return Logical. Indicates if \code{x} is of class \code{hazard}.
#' @keywords internal
#'
is.hazard <- function(x) {
  return(inherits(x, "hazard"))
}

#' @export
print.hazard <- function(x, ...) {
  cat("Hazard object with age-specific hazard rates spanning from age",
      x$partition[1], "to age",  x$partition[length(x$partition)])
  cat("\n")
  if (length(x$subtype_ID) == 1){
    cat("Lifetime risk of disease = ", 1 - exp(-sum(diff(x$partition)*x$hazardDF[, 1])), "\n")
  } else {
    for (i in 1:length(x$subtype_ID)) {
      cat("Lifetime risk of ", x$subtype_ID[i], " = ",  1 - exp(-sum(diff(x$partition)*x$hazardDF[, i])), "\n")
    }
  }

}
