#' Determine if a pedigree is ascertained
#'
#' Intended priamrily as an internal function, \code{is_ascertained} checks to see if a pedigree returned by \code{\link{sim_ped}} is ascertained.
#'
#' @inheritParams trim.ped
#' @inheritParams sim_RVped
#'
#' @return  A list containing the following data frames:
#' @return \code{ascertained} Logical.  Indicates if pedigree is ascertained.
#' @export
#'
#' @examples
#' data(AgeSpecific_Hazards)
#' my_HR <- hazard(hazardDF = AgeSpecific_Hazards)
#'
#' #Simulate a random pedigree
#' set.seed(2)
#' ex_ped <- sim_ped(hazard_rates = my_HR,
#'                   GRR = 50, carrier_prob = 0.002,
#'                   FamID = 1,
#'                   RVfounder = TRUE,
#'                   founder_byears = c(1900, 1910),
#'                   stop_year = 2015)
#'
#' ex_ped
#' is_ascertained(ped_file = ex_ped,
#'                num_affected = 2,
#'                ascertain_span = c(2000, 2015),
#'                recall_probs = c(1, 1, 0.5, 0.25))[[1]]
#'
#' data(EgPeds)
#' EgPeds[EgPeds$FamID == 1, ]
#' is_ascertained(ped_file = new.ped(EgPeds[EgPeds$FamID == 1, -15]),
#'                num_affected = 2,
#'                ascertain_span = c(2000, 2015),
#'                recall_probs = c(1))[[1]]
#'
is_ascertained <- function(ped_file, num_affected, ascertain_span, recall_probs){

  # prior to sending the simulated pedigree to the trim function,
  # we check to see if it meets the required criteria for number of
  # affected.  If it does, we choose a proband from the available
  # candidates prior to sending it to the trim.ped function.
  if (disqualify_ped(ped_file, num_affected, ascertain_span)) {
    ascertained <- FALSE
    return_ped = ped_file
  } else {
    #choose a proband
    pro_ped <- choose_proband(ped_file, num_affected, ascertain_span)

    # Now that we have a full pedigree that meets our conditions, we trim the
    # pedigree and check to see that the trimmed pedigree STILL meets our
    # conditions, we then update ascertained appropriately.
    if (missing(recall_probs)) {
      ascertained_ped <- trim.ped(ped_file = pro_ped)
    } else {
      ascertained_ped <- trim.ped(ped_file = pro_ped, recall_probs)
    }

    ascertained <- ascertainTrim_ped(ped_file = ascertained_ped, num_affected)
    return_ped = ascertained_ped
  }

  return(list(ascertain = ascertained,
              ped_file = return_ped))
}
