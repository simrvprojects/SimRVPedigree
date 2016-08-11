#' Check to see if hazard and partition vectors are specified correctly
#'
#'\code{check_hazpart} Check to see if hazard and partition vectors are specified correctly
#'
#'@param hazard A numeric vector.  A vector of age-specific hazards.
#'@param part A numeric vector. The partition over which to apply hazard, note length(part) == length(hazard) + 1 should return TRUE.
#'
#'@examples
#'check_hazpart(hazard = c(0.1, 0.4, 0.7, 0.8, 0.9), part = c(0, 20, 40, 60, 80, 100))
#'check_hazpart(hazard = c(0.1, 0.4, 0.7, 0.8, 0.9), part = c(0, 20, 40, 60, 80))

check_hazpart = function(hazard, part){
  check1 <- (class(hazard) != "numeric")
  check2 <- (length(hazard) < 1)
  check3 <- (length(part) != (length(hazard) + 1))
  if ( check1 | check2 | check3 ) {
    stop ('please provide numeric hazard and part vectors,
          such that length(part) == length(hazard) + 1 returns TRUE')
  }
}



#' define a function that issues a warning if age-specific hazards are notsupplied for an appropriate span of ages
#'
#'\code{check_part} Check to see if the partition vector meets some basic properties
#'
#'@param part A numeric vector. Partiton over which to apply hazards
#'@examples
#'check_part(part = c("a"))
#'check_part(part = c(-1, 0, 10))

check_part = function(part){
  check1 <- min(part)
  check2 <- max(part)
  if ( check1 > 20 | check2 < 60 ) {
    warning ('For optimal results please specify age-specific hazards which begin near birth and end near the life expectancy of the population to which the age specific hazards should be applied.')
  }
}


#check_recall_probs
#check_spans
