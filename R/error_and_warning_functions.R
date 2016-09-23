#' Check to see if age-specific hazard and partition of ages are specified correctly
#'
#' @inheritParams get_WaitTime
#' @export
#'
#' @examples
#' #the length of part should be equal to the length of hazard + 1
#' check_hazpart(hazard = seq(0, 1, by = 0.05), part = seq(0, 21, by = 1))
#' \dontrun{check_hazpart(hazard = seq(0, 1, by = 0.05),
#'                        part = seq(0, 20, by = 1))}
#' \dontrun{check_hazpart(hazard = c(1, NA, 2, 4),
#'                        part = c(0, 10, 20, 30, 40))}
#'
check_hazpart = function(hazard, part){
  if ( class(hazard) != "numeric" | class(part) != "numeric" ){
    stop ('please provide numeric hazard and partition vectors')
  } else if (length(part) == 1 | length(part) != (length(hazard) + 1)) {
    stop ('please provide numeric hazard and partition vectors, such that length(part) == length(hazard) + 1')
  } else if (any(is.na(hazard)) | any(is.na(part))) {
    stop('hazard and partition vectors cannot contain missing values')
  }
}

#' Check to see if hazards span an appropriate range of ages
#'
#' @inheritParams get_WaitTime
#'
#' @export
#'
#' @examples
#' \dontrun{check_part(part = seq(10, 20, by = 1))}
#' \dontrun{check_part(part = seq(-10, 10, by = 1))}
#' check_part(part = seq(0, 85, by = 1))
check_part = function(part){
  if (min(part) != 0) {
    stop('age-specific hazards must begin at birth')
  } else if (max(part) < 65){
    warning ('For optimal results please specify age-specific hazards that begin near birth and end near the life expectancy of the population to which the age-specific hazards apply.')
    }
}


#' Check to see if recall probabilies are correctly specified
#'
#' @inheritParams trim_pedigree
#'
#' @export
#'
#' @examples
#' \dontrun{check_rprobs(c(2, 0.3, 0.4))
#'          check_rprobs(NA)
#'          check_rprobs("a")}
#' check_rprobs(recall_probs = c(1))
check_rprobs = function(recall_probs){
  if (any(recall_probs > 1) | any(recall_probs < 0) ){
    stop ('recall probabilities must be between 0 and 1')
  }
}

#' Check to see if time spans are appropriately specified
#'
#' @param span numeric.  A span of years or ages.
#'
#' @export
#'
#' @examples
#' check_spans(c(1975, 1982))
#' \dontrun{
#' check_spans(c(1, 3, 6))
#' check_spans(c(1982, 1975))
#' check_spans(c(1975, 1975))
#' check_spans(c(18))}
#'
check_spans = function(span){
  if (length(span) != 2 | span[1] >= span[2]){
    stop ('please provide appropriate time/age spans')
  }
}

#' Check to see if death_hazard is appropriately specified
#'
#' @param death_hazard data.frame. The age-specific death hazards for affected and unaffected individuals.
#'
#' @export
check_dhaz = function(death_hazard){
  if(class(death_hazard) != "data.frame"  ){
    stop("death_hazard must be a data frame with 2 columns:,
         column 1 = unaffected age-specific death hazard,
         column 2 = affected age-specific death hazard")
  }else if(ncol(death_hazard) != 2){
    stop("death_hazard must be a data frame with 2 columns:,
         column 1 = unaffected age-specific death hazard,
         column 2 = affected age-specific death hazard")
  }
}
