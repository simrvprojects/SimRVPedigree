#' Determine founder genotype at the disease locus and determime their relative-risk of disease
#'
#' @inheritParams sim_RVped
#' @param intro_RV Logical. If \code{intro_RV = TRUE} a founder has introduced a rare, causal variant to the pedigree, otherwise no RV has been introduced to the pedigree.
#'
#' @return A list containing: 1. the founder's genotype at the disease locus, their relative-risk of disease, and an updated value for intro_RV.
#'
#' @keywords internal
#'
sim_founderRVstatus <- function(GRR, allele_freq, RVfounder, intro_RV){
  if (GRR == 1 |
      (RVfounder == "single" & intro_RV == TRUE) |
      (RVfounder == "first" & intro_RV == TRUE)) {
    # If GRR (genetic relative risk) = 1, the variant is not associated with
    # the disease; hence we do not allow an RV to segregate in the pedigree
    # Additionally, we set intro_RV to TRUE, so that we do not needlessly search
    # for a founder to introduce the RV
    d_locus <- c(0, 0)
    fRR <- 1
  } else if (RVfounder == "single" & intro_RV == FALSE){
    # If a single founder has not been located we simulate allele status at the
    # disease locus given allele_freq, set RR and update intro_RV appropriately
    carrier_prob <- 1 - (1 - allele_freq)^2
    d_locus <- sample(x = c(0, ifelse(runif(1) <= (1 - (1 - carrier_prob)^2), 1, 0)),
                      size = 2, replace = F)
    fRR <- ifelse(any(d_locus == 1), GRR, 1)
    intro_RV <- ifelse(any(d_locus == 1), T, F)
  } else if (RVfounder == "multiple"){
    # If multiple has been specified we allow fouders to introduce the RV with
    # probability proprotional to its allele frequency in the population, and
    # set RR and update intro_RV appropriately
    d_locus <- ifelse(runif(2) <= allele_freq, 1, 0)
    fRR <- ifelse(any(d_locus == 1), GRR, 1)
    intro_RV <- ifelse(any(d_locus == 1), T, F)
  } else {
    d_locus <- sample(x = c(0, 1), size = 2, replace = F)
    fRR <- GRR
    intro_RV <- TRUE
  }
  founder_dat <- list(d_locus, fRR, intro_RV)
  return(founder_dat)
}

#' Choose a proband from the disease-affected relatives in a pedigree
#'
#' @param ped_file Pedigree simulated by \code{sim_ped}.
#' @inheritParams sim_RVped
#'
#' @return Pedigree with proband selected.
#' @keywords internal
#'
choose_proband = function(ped_file, num_affected, ascertain_span){
  #initialize proband ID variable
  ped_file$proband <- 0

  #Gather info on affecteds
  A_ID <- ped_file[which(ped_file$affected == 1),
              which(colnames(ped_file) %in% c("onsetYr", "ID", "proband"))]
  A_ID <- A_ID[order(A_ID$onsetYr), ]
  A_ID <- A_ID[which(A_ID$onsetYr <= ascertain_span[2]), ]
  A_ID$proband <- ifelse(A_ID$onsetYr %in% ascertain_span[1]:ascertain_span[2], 1, 0)

  if (sum(A_ID$proband) == 1) {
    #In this scenario we have only 1 candidate proband
    #NOTE: sim_RVped has already checked to make sure
    #that there was another affected prior to this one

    ped_file$proband[which(ped_file$ID == A_ID$ID[which(A_ID$proband == 1)])] <- 1

  } else if (sum(abs(A_ID$proband - 1)) > (num_affected - 1)) {

    #multiple available probands and the n-1 affected condition has already
    #been met by start of ascertainment period, so simply choose randomly
    #amongst available probands
    probandID <- sample(size = 1,
                        x = A_ID$ID[which(A_ID$proband == 1)])
    ped_file$proband[which(ped_file$ID == probandID)] <- 1

  } else {

    #no affecteds before ascertainment period, must choose from among
    #the nth or greater to experience onset
    A_ID$proband[1:(num_affected - 1)] <- 0
    #must write additional if statement here because of R's interesting
    #take on how sample should work when there is only one 1 to sample from....
    if (sum(A_ID$proband) == 1) {
      ped_file$proband[which(ped_file$ID == A_ID$ID[which(A_ID$proband == 1)])] <- 1
    } else {
      probandID <- sample(size = 1,
                          x = A_ID$ID[which(A_ID$proband == 1)])
      ped_file$proband[which(ped_file$ID == probandID)] <- 1
    }
  }

  return(ped_file)
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
