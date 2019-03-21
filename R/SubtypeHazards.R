#' Age-specific hazard rate dataset with subtypes
#'
#' A dataset that contains the following age-specific hazard rates: (1) the age-specific hazard rates for Hodgkin Lymphoma and Non-Hodgkin Lymphoma in the United States, (2) the age-specific hazard rates for death in the United States, and (3) the age-specific hazard rates for death for individuals, living in the United States, who have been diagnosed with either Hodgkin Lymphoma and Non-Hodgkin Lymphoma.
#'
#'  The \code{SubtypeHazards} dataset contains the following age-specific hazard rates: (1) the age-specific hazard rates for Hodgkin Lymphoma and Non-Hodgkin Lymphoma in the United States, (2) the age-specific hazard rates for death in the United States, and (3) the age-specific hazard rates for death for individuals, living in the United States, who have been diagnosed with either Hodgkin Lymphoma and Non-Hodgkin Lymphoma.  The age-specific hazard rates of disease onset and death in the affected population were estimated by the Surveillance, Epidemiology, and End Results Program (SEER) SEER*Stat Software, and the age-specific hazard rates of death in the United States may be estimated from actuarial life tables provided by the Social Security Administration.
#'
#'  The three columns in the \code{SubtypeHazards} dataset provide age-specific hazard rates, in yearly increments, beginning at age 0 and ending with age 100.  That is, the values in the first row describe the hazard rates for an individual whose age is contained in the interval [0, 1), while the values in the second row describe the hazard rates for an individual whose age is contained in the interval [1, 2), and so on.
#'
#' @docType data
#'
#' @references The Surveillance, Epidemiology, and End Results (SEER) Program. \url{https://seer.cancer.gov/}
#' @references Bell, F. C., Miller, M. L. (2005). \emph{Life Tables for the United States Social Security Area, 1900-2100}. Baltimore, Md.: Social Security Administration, Office of the Chief Actuary.
#'
#' @usage data(SubtypeHazards)
#'
#' @format A data frame with 100 rows and 4 variables:
#' \describe{
#' \item{pop_HL_hazard}{The population, age-specific hazard rate for Hodgkin Lymphoma}
#' \item{pop_NHL_hazard}{The population, age-specific hazard rate for Non-Hodgkin Lymphoma}
#' \item{unaffected_death_hazard}{The age-specific hazard rate for death in the \strong{unaffected} population}
#' \item{affected_death_hazard}{The age-specific hazard rate for death in the \strong{affected} population}
#' }
#'
"SubtypeHazards"
