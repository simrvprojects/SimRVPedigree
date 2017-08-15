#' Initial checks to disqualify a pedigree from ascertainment.
#'
#' @inheritParams choose_proband
#'
#' @return Logical. If TRUE, pedigree is discarded.
#' @keywords internal
disqualify_ped <- function(ped_file, num_affected, ascertain_span){
  length(which(ped_file$onsetYr <= ascertain_span[2])) < num_affected |
    length(which(ped_file$onsetYr %in% ascertain_span[1]:ascertain_span[2])) < 1
}

#' Check to see if a trimmed pedigree is ascertained.
#'
#' @inheritParams choose_proband
#'
#' @return Logical. If TRUE, pedigree is ascertained.
#' @keywords internal
ascertainTrim_ped <- function(ped_file, num_affected){

  #Gather the onset years for all affecteds, and for the proband.
  #We need to ensure that at least num_affected - 1 were affected prior to
  #the proband for the pedigree to be ascertained.
  POyear <- ped_file$onsetYr[ped_file$proband == 1]

  Oyears <- ped_file$onsetYr[ped_file$affected == 1
                             & ped_file$available == 1
                             & ped_file$proband == 0]

  #determine the number of available affected individuals
  ascertained <- sum(Oyears <= POyear) >= (num_affected - 1)

  return(ascertained)
}


#' Determine if a previously ascertained pedigree meets new criteria.
#'
#' @inheritParams simRVped
#' @param asc_ped Data.frame.  The previously ascertained pedigree. Optionally, asc_ped may contain ped_file information for the affected relatives only.
#' @param ref_year A numeric constant.  The reference year used to determine current age for pedigree members.
#' @return Logical. If TRUE, pedigree is re-ascertained.
#' @export
#' @examples
#' #Read in age-specific hazards
#' data(AgeSpecific_Hazards)
#'
#' my_HR <- hazard(hazardDF = AgeSpecific_Hazards)
#'
#' #Simulate pedigree ascertained for multiple affected individuals
#' set.seed(8008135)
#' ex_ped <- simPed(hazard_rates = my_HR,
#'                   GRR = 50, carrier_prob = 0.002,
#'                   RVfounder = TRUE,
#'                   FamID = 1,
#'                   founder_byears = c(1900, 1910))
#'
#' head(ex_ped, n = 3)
#' Aped <- is_ascertained(ped_file = ex_ped,
#'                        num_affected = 2,
#'                        ascertain_span = c(2000, 2015),
#'                        recall_probs = c(1, 1, 0.5, 0.25))
#'
#' Aped[[1]]
#'
#' Aped[[2]][(Aped[[2]]$affected == 1 & Aped[[2]]$available == 1), ]
#'
#' is_reAscertained(asc_ped = Aped[[2]],
#'                  num_affected = 3,
#'                  ref_year = 2012)
#'
#' is_reAscertained(asc_ped = Aped[[2]],
#'                  num_affected = 3,
#'                  ref_year = 2016)
#'
is_reAscertained <- function(asc_ped, num_affected, ref_year){

  if (!is.ped(asc_ped)) {
    stop("\n \n Expecting a ped object. \n Please use new.ped to create an object of class ped.")
  }

  if (!("proband" %in% colnames(asc_ped))) {
    warning("Proband missing: is_reAscertained is intended for ascertained pedigrees.")
  }
  # Gather the onset years for all affecteds.
  # Check to see if the pedigree had at least num_affected
  # individuals at the selected ref_year.
  Oyears <- asc_ped$onsetYr[asc_ped$affected == 1 & asc_ped$available == 1]


  #determine the number of available affected individuals
  ascertained <- sum(Oyears <= ref_year) >= num_affected

  return(ascertained)
}


#' Determine if a pedigree is ascertained
#'
#' Intended priamrily as an internal function, \code{is_ascertained} checks to see if a pedigree returned by \code{\link{simPed}} is ascertained.
#'
#' @inheritParams trim.ped
#' @inheritParams simRVped
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
#' ex_ped <- simPed(hazard_rates = my_HR,
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
#' is_ascertained(ped_file = ped(EgPeds[EgPeds$FamID == 1, -15]),
#'                num_affected = 2,
#'                ascertain_span = c(2000, 2015),
#'                recall_probs = c(1))[[1]]
#'
is_ascertained <- function(ped_file, num_affected, ascertain_span, recall_probs){

  # prior to sending the simulated pedigree to the trim function,
  # we check to see if it meets the required criteria for number of
  # affected.  If it does, we choose a proband from the available
  # candidates prior to sending it to the trim.ped function.
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


#' Checks ped_file for expected information, used before converting to ped object.
#'
#' @param ped_file data.frame The pedigree in data frame format.
#'
#' @keywords internal
#'
check_ped <- function(ped_file){
  if (class(ped_file) != "data.frame") {
    stop("ped_file must be a data.frame with the following variables: FamID, ID, dadID, momID, sex, affected")
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
