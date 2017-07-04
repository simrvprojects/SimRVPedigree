#' Create a hazard object for input to \code{sim_RVped}, \code{sim_ped}, and \code{sim_lifeEvents}.
#'
#' @param partition Numeric vector. The partition of ages over which to apply the age-specific hazard rates in \code{hazardDF}.
#' @param hazardDF Data.frame. Column 1 should specify the age-specific hazard rate for individuals with no genetic predisposition to the disease, column 2 should specify the age-specific hazard rate for death in the unaffected population, and column 3 should specify the age-specific hazard rate for death in the affected population.
#'
#' @return An object of class \code{hazard}.
#' @export
#'
#' @examples
#' my_HR <- new.hazard(partition = seq(0, 100, by = 1),
#'                     hazardDF = AgeSpecific_Hazards)
#'
#' class(my_HR)
#' head(my_HR[[2]])
#'
new.hazard <- function(partition, hazardDF) {
  # Set up a new object of class Hazard

  if (any(is.na(partition))) {
    stop('partition cannot contain missing values')
  }

  if (any(is.na(hazardDF))) {
    stop('hazardDF cannot contain missing values')
  }

  if (min(partition) != 0) {
    stop('age-specific hazards must begin at birth (age 0)')
  } else if (max(partition) < 65){
    warning ('For optimal results please specify age-specific hazard rates that begin at birth and end near the life expectancy of the population to which the age-specific hazards apply.')
  }

  if (length(partition) == 1 |
      length(partition) != (nrow(hazardDF) + 1) ) {
    stop ('please provide hazard rates, such that length(partition) == length(hazard) + 1')
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
  return(hazard(partition, hazardDF))
}

#' Constructor function of object of class \code{hazard}
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
hazard <- function(partition, hazardDF) {
  obj <- list(partition = partition,
              hazardDF = hazardDF)
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
