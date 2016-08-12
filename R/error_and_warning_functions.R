check_hazpart = function(hazard, part){
  check1 <- (class(hazard) != "numeric")
  check2 <- (length(hazard) < 1)
  check3 <- (length(part) != (length(hazard) + 1))
  if ( check1 | check2 | check3 ) {
    stop ('please provide numeric hazard and part vectors,
          such that length(part) == length(hazard) + 1 returns TRUE')
  }
}

check_part = function(part){
  check1 <- min(part)
  check2 <- max(part)
  if ( check1 > 20 | check2 < 60 ) {
    warning ('For optimal results please specify age-specific hazards which begin near birth and end near the life expectancy of the population to which the age specific hazards should be applied.')
  }
}


check_rprobs = function(recall_probs){
  if (sum(recall_probs > 1) > 0 ){
    stop ('recall probabilities must be less than or equal to 1')
  }
}

check_spans = function(span){
  if (length(span) != 2 | span[1] >= span[2]){
    stop ('please provide an appropriate time span')
  }
}
