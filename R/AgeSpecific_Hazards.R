#' Example age-specific hazards
#'
#' A dataset containing 3 columns of age-specific hazards.  Each age-specific hazard column lists the age-specific hazards in yearly increments, beginning at age 0 ending with age 100.  That is, values in the first row are the age-specific hazards for individuals whose ages fall in the interval [0, 1), the values in the first row are the age-specific hazards for individuals whose ages fall in the interval [1, 2), etc.
#'
#' The hazards in this dataset were simulated randomly and do not represent the hazards for a real disease.  These simulated hazards have been included in the `SimRVPedigree` package to simplify illustrating usage of package functions.
#'
#' @docType data
#'
#' @usage data(AgeSpecific_Hazards)
#'
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#' \item{pop_onset_hazard}{The age-specific population onset hazard}
#' \item{unaffected_death_hazard}{The age-specific death hazard for the UNAFFECTED population}
#' \item{affected_death_hazard}{The age-specific death hazard for the AFFECTED population}
#' }
#'
"AgeSpecific_Hazards"
