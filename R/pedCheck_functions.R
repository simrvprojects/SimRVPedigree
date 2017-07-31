#' Initial checks to disqualify a pedigree from ascertainment.
#'
#' @inheritParams choose_proband
#'
#' @return Logical. If TRUE, pedigree will not be ascertained.
#' @keywords internal
disqualify_ped <- function(ped_file, num_affected, ascertain_span){
  nrow(ped_file) < num_affected | sum(ped_file$affected) < num_affected |
    length(ped_file$ID[which(ped_file$onsetYr <= ascertain_span[2])]) < num_affected |
    length(ped_file$ID[which(ped_file$onsetYr %in%
                               ascertain_span[1]:ascertain_span[2])]) < 1
}

#' Check to see if a trimmed pedigree is ascertained.
#'
#' @inheritParams choose_proband
#'
#' @return Logical. If TRUE, pedigree is ascertained.
#' @keywords internal
ascertain_trimmedPed <- function(ped_file, num_affected){
  #Gather the onset years for all remaining affecteds, and for the proband.
  #We need to ensure that at least num_affected - 1 were affected prior to
  #the proband for the pedigree to be ascertained.
  POyear <- ped_file$onsetYr[ped_file$proband == 1]

  Oyears <- ped_file$onsetYr[which(ped_file$affected == 1 &
                                     ped_file$available == 1 &
                                     ped_file$proband == 0)]


  #determine the number of available affected individuals
  ascertained <- sum(Oyears <= POyear) >= (num_affected - 1)

  return(ascertained)
}

#' Determine if a pedigree is ascertained
#'
#' Intended priamrily as an internal function, \code{is_ascertained} checks to see if a pedigree returned by \code{\link{sim_ped}} is ascertained.
#'
#' @inheritParams trim_ped
#' @inheritParams sim_RVped
#'
#' @return  A list containing the following data frames:
#' @return \code{ascertained} Logical.  Indicates if pedigree is ascertained.
#' @export
#'
#' @examples
#' data(AgeSpecific_Hazards)
#' my_HR <- new.hazard(hazardDF = AgeSpecific_Hazards)
#'
#' #Simulate a random pedigree
#' set.seed(37)
#' ex_ped <- sim_ped(hazard_rates = my_HR,
#'                   GRR = 50, carrier_prob = 0.002,
#'                   FamID = 1,
#'                   RVfounder = "first",
#'                   founder_byears = c(1900, 1910),
#'                   stop_year = 2015)
#'
#' ex_ped
#' is_ascertained(ped_file = ex_ped,
#'                num_affected = 2,
#'                ascertain_span = c(2000, 2015),
#'                recall_probs = c(1, 1, 0.5, 0.25))
#'
#' data(EgPeds)
#' EgPeds[EgPeds$FamID == 1, ]
#' is_ascertained(ped_file = EgPeds[EgPeds$FamID == 1, -15],
#'                num_affected = 2,
#'                ascertain_span = c(2000, 2015),
#'                recall_probs = c(1))
#'
is_ascertained <- function(ped_file, num_affected, ascertain_span, recall_probs){

  # prior to sending the simulated pedigree to the trim function,
  # we check to see if it meets the required criteria for number of
  # affected.  If it does, we choose a proband from the available
  # candidates prior to sending it to the trim_ped function.
  if( disqualify_ped(ped_file, num_affected, ascertain_span) ){
    ascertained <- FALSE
    return_ped = ped_file
  } else {
    #choose a proband
    pro_ped <- choose_proband(ped_file, num_affected, ascertain_span)

    # Now that we have a full pedigree that meets our conditions, we trim the
    # pedigree and check to see that the trimmed pedigree STILL meets our
    # conditions, we then update ascertained appropriately.
    if (missing(recall_probs)) {
      ascertained_ped <- trim_ped(ped_file = pro_ped)
    } else {
      ascertained_ped <- trim_ped(ped_file = pro_ped, recall_probs)
    }

    ascertained <- ascertain_trimmedPed(ped_file = ascertained_ped, num_affected)
    return_ped = ascertained_ped
  }

  return(list(ascertain = ascertained,
              ped_file = return_ped))
}


#' Checks ped_file input for expected information.
#'
#' @inheritParams trim_ped
#'
#' @keywords internal
#'
check_ped <- function(ped_file){
  if (class(ped_file) != "data.frame") {
    stop("please provide a data.frame with the following variables: FamID, ID, dadID, momID, sex, affected")
  }


  if (!"FamID" %in% colnames(ped_file) |
      !"ID" %in% colnames(ped_file) |
      !"dadID" %in% colnames(ped_file) |
      !"momID" %in% colnames(ped_file) |
      !"sex" %in% colnames(ped_file) |
      !"affected" %in% colnames(ped_file)) {
    stop('please provide a data.frame with the following variables: FamID, ID, dadID, momID, sex, affected')
  }

  if(any(is.na(ped_file$ID))) {
    stop('ID contains missing values.  Please ensure all individuals have a valid ID.')
  }

  moms <- unique(ped_file$momID[!is.na(ped_file$momID)])
  dads <- unique(ped_file$dadID[!is.na(ped_file$dadID)])

  if (any(ped_file$sex[which(ped_file$ID %in% moms)] != 1) |
      any(ped_file$sex[which(ped_file$ID %in% dads)] != 0)){

    wrong_sex <- c(ped_file$ID[which(ped_file$sex[which(ped_file$ID %in% dads)] != 0)],
                   ped_file$ID[which(ped_file$sex[which(ped_file$ID %in% moms)] != 1)])

    stop(paste0('Sex improperly specifed ID: ', sep = '', wrong_sex, '.  Please ensure that for males: sex = 0; and for females: sex = 1.'))
  }

  if (any(!moms %in% ped_file$ID) | any(!dads %in% ped_file$ID)) {

    wrong_par <- c(ped_file$ID[which(ped_file$momID == moms[which(!moms %in% ped_file$ID)])],
                   ped_file$ID[which(ped_file$dadID == dads[which(!dads %in% ped_file$ID)])])

    stop(paste0('ID: ', sep = '', wrong_par, '.  Non-founders must have both a mother and a father, while founders have neither.'))
  }

  if (any(!is.na(ped_file$momID[is.na(ped_file$dadID)])) |
      any(!is.na(ped_file$dadID[is.na(ped_file$momID)]))) {
    stop("Non-founders must have both a mother and a father, while founders have missing momID and dadID.")
  }
}
