#' Create an object of class hazard.
#'
#' Create a new hazard object for input to \code{\link{sim_RVped}}, \code{\link{sim_ped}}, and \code{\link{sim_lifeEvents}}.
#'
#' @param hazardDF Data.frame. Column 1 should specify the age-specific hazard rate for individuals with no genetic predisposition to the disease, column 2 should specify the age-specific hazard rate for death in the unaffected population, and column 3 should specify the age-specific hazard rate for death in the affected population.
#' @param partition Numeric vector. Optional. The partition of ages over which to apply the age-specific hazard rates in \code{hazardDF}. If missing, it is assumed that \code{partition} starts at 0 and increases in yearly increments.
#'
#' @return An object of class hazard.
#' @export
#'
#' @examples
#' my_HR <- new.hazard(hazardDF = AgeSpecific_Hazards)
#'
#' class(my_HR)
#' head(my_HR[[2]])
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
