#' Create an object of class hazard.
#'
#' Create a new hazard object for input to \code{\link{sim_RVped}}, \code{\link{sim_ped}}, and \code{\link{sim_lifeEvents}}.
#'
#' The data frame supplied as \code{hazardDF} must contain 3 columns that meet the following criteria:
#' \describe{
#'  \item{column 1:}{age-specific hazard rates of \emph{disease} for the population of interest}
#'  \item{column 2:}{age-specific hazard rates of \emph{death} for the \strong{unaffected} population.  If the disease of interest is sufficiently rare, so that death by the disease is rare, the user may choose to use the population, age-specific, hazard rates of death instead.}
#'  \item{column 3:}{age-specific hazard rates of \emph{death} for the \strong{affected} population.}
#' }
#'
#' We expect that users provide \code{partition} in year units; e.g. a hazard rate for a baby between 6 months and 1 year of age should have lower bound 0.5 years and an upper bound 1 year.  Additionally, \code{partition} must be a valid parition for each of the age-specific hazard rates in \code{hazardDF}.
#'
#' @param hazardDF Data.frame. Column 1 specifies the age-specific hazard rate of \emph{disease} in the population of interest, column 2 specifies the age-specific hazard rate for \emph{death} in the \strong{unaffected} population, and column 3 specifies the age-specific hazard rate for \emph{death} in the \strong{affected} population.  See details.
#' @param partition Numeric vector. Optional. The partition of ages, in years, over which to apply the age-specific hazard rates in \code{hazardDF}. If missing, it is assumed that \code{partition} starts at 0 and increases in yearly increments.  See details.
#'
#' @return An object of class hazard.
#' @export
#'
#' @examples
#' data(AgeSpecific_Hazards)
#'
#' head(AgeSpecific_Hazards)
#' nrow(AgeSpecific_Hazards)
#'
#' my_part = seq(0, 100, by = 1)
#' my_HR <- new.hazard(hazardDF = AgeSpecific_Hazards,
#'                     partition = my_part)
#' class(my_HR)
#' head(my_HR[[1]])
#' my_HR[[2]]
#'
#'
#' my_HR <- new.hazard(hazardDF = AgeSpecific_Hazards)
#' class(my_HR)
#' my_HR[[2]]
#'
new.hazard <- function(hazardDF, partition) {
  # Set up a new object of class Hazard

  if(missing(partition)){
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


  if (class(hazardDF) != "data.frame") {
    stop("hazardDF must be a data frame with 3 columns:
         column 1 = population age-specific hazard rate of disease,
         column 2 = age-specific hazard rate of death in the unaffected population,
         column 3 = age-specific hazard rate of death in the affected population")
  } else if (ncol(hazardDF) != 3) {
    stop("hazardDF must be a data frame with 3 columns:
         column 1 = population age-specific hazard rate of disease,
         column 2 = age-specific hazard rate of death in the unaffected population,
         column 3 = age-specific hazard rate of death in the affected population")
  } else if(sum(hazardDF[, 2] > hazardDF[, 3]) > nrow(hazardDF)/2 |
            sum(which(hazardDF[, 2] == hazardDF[, 3])) != 0){
    warning("Please check that you have specified hazardDF such that:
            column 2 = age-specific hazard rate of death in the UNAFFECTED population,
            column 3 = age-specific hazard rate of death in the AFFECTED population")
  }


  # Return Hazard object with the user-supplied hazard rates and partition
  return(hazard(hazardDF, partition))
}

#' Constructor function for an object of class \code{hazard}
#'
#' @inheritParams new.hazard
#'
#' @return an object of class \code{hazard}.
#' @export
#'
#' @examples
#' data(AgeSpecific_Hazards)
#' my_HR <- hazard(partition = seq(0, 100, by = 1),
#'                 hazardDF = AgeSpecific_Hazards)
#' class(my_HR)
#'
hazard <- function(hazardDF, partition) {
  obj <- list(hazardDF = hazardDF,
              partition = partition)
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
