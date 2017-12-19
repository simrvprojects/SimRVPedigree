#' Age-specific hazard rate dataset
#'
#' A dataset that contains age-specific hazard rates to \strong{roughly mimic}: (1) the age-specific hazard rates for lymphoid cancer in the United States, (2) the age-specific hazard rates for death in the United States, and (3) the age-specific hazard rates for death for individuals, living in the United States, who have been diagnosed with a lymphoid cancer.
#'
#'  The \code{AgeSpecific_Hazards} dataset contains age-specific hazard rates which \strong{roughly mimic}: (1) the age-specific hazard rates for lymphoid cancer onset in the United States, (2) the age-specific hazard rates for death in the United States, and (3) the age-specific hazard rates for death for individuals, living in the United States, who have been diagnosed with a lymphoid cancer.  The age-specific hazard rates of lymphoid cancer onset and death in the affected population may be estimated by a program such as the Surveillance, Epidemiology, and End Results Program (SEER), and the age-specific hazard rates of death in the United States may be estimated from actuarial life tables provided by the Social Security Administration.
#'
#'  The three columns in the \code{AgeSpecific_Hazards} dataset provide age-specific hazard rates, in yearly increments, beginning at age 0 and ending with age 100.  That is, the values in the first row describe the hazard rates for an individual whose age is contained in the interval [0, 1), while the values in the second row describe the hazard rates for an individual whose age is contained in the interval [1, 2), and so on.
#'
#'
#'
#' @docType data
#'
#' @references The Surveillance, Epidemiology, and End Results (SEER) Program. \url{https://seer.cancer.gov/}
#' @references Bell, F. C., Miller, M. L. (2005). \emph{Life Tables for the United States Social Security Area, 1900-2100}. Baltimore, Md.: Social Security Administration, Office of the Chief Actuary.
#'
#' @usage data(AgeSpecific_Hazards)
#'
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#' \item{pop_onset_hazard}{The age-specific population hazard rate for lymphoid cancer}
#' \item{unaffected_death_hazard}{The age-specific hazard rate for death in the \strong{unaffected} population}
#' \item{affected_death_hazard}{The age-specific hazard rate for death in the \strong{affected} population}
#' }
#'
"AgeSpecific_Hazards"
