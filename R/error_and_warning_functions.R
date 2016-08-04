##-----------------##
##  check_hazpart  ##
##-----------------##
## define a function that checks to make sure all the appropriate arguments 
##  have been provided in NHPP functions

## Arguments--------------------------------------------------------------------
## hazard    - numeric vector of length n
## part - numeric vector of length n + 1

## Function Requirements--------------------------------------------------------
## NONE

## Package Requirements---------------------------------------------------------
## NONE

check_hazpart = function(hazard, part){
  check1 <- (class(hazard) != "numeric")
  check2 <- (length(hazard) < 1)
  check3 <- (length(part) != (length(hazard) + 1))
  if ( check1 | check2 | check3 ) {
    stop ('please provide numeric hazard and part vectors, with length(part) = length(hazard) + 1')
  }
}


##--------------##
##  check_part  ##
##--------------##
## define a function that issues a warning if the partition does not apply to 
## a wide enough span of ages

## Arguments--------------------------------------------------------------------
## part - numeric vector of length n + 1

## Function Requirements--------------------------------------------------------
## NONE

## Package Requirements---------------------------------------------------------
## NONE

check_part = function(part){
  check1 <- min(part)
  check2 <- max(part)
  check3 <- length(part)
  if ( check1 > 20 | check2 < 60 | check3 < 5 ) {
    stop ('please provide numeric hazard and part vectors, with length(part) = length(hazard) + 1')
  }
}
