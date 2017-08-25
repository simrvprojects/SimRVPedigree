#' Reassign generation number based on affected status
#'
#' The \code{reassign_gen} function assigns generation numbers among affected family members so that generation 1 represents the most recent generation that a putative disease variant shared identical by descent (IBD), as defined in Thompson (2000), by affected members could have been introduced into the pedigree.
#'
#' \emph{The \code{reassign_gen} function is primarily intended for pedigrees simulated by the \code{simRVPedigree} package.  Results will be inaccurate if inbreeding or loops are present.}
#'
#' The \code{reassign_gen} function accepts a pedigree and reassigns generation numbers among disease-affected relatives so that generation 1 represents the generation of the most recent common ancestor of all disease-affected relatives.  We note that the individual in generation 1 could themselves be disease-affected, i.e. an individual can be considered their own ancestor.
#'
#' For example, consider a family with 2 affected members.  If the disease-affected relatives are a parent and a child, the affected parent would be assigned generation 1, and the affected child generation 2.  However, if the disease-affected relatives are a pair of siblings, each is be assigned generation 2 since a common parent of the two is assumed to be a carrier of a latent susceptibility variant.  Similarly, if the disease-affected relatives are a pair of cousins, is assigned generation 3, since a common grandparent is the most recent common ancestor from whom they could have inherited a shared variant associated with the disease.
#'
#' Users who wish to assign generation number based on affection status in pedigrees that have not been simulated with the \code{SimRVpedigree} package must create a ped object using \code{new.ped}.  This \code{ped} object \emph{must} contain the following variables for each pedigree member:
#'
#' \tabular{lll}{
#' \strong{name} \tab \strong{type} \tab \strong{description} \cr
#' \code{FamID} \tab numeric \tab family identification number \cr
#' \code{ID} \tab numeric \tab individual identification number \cr
#' \code{dadID} \tab numeric \tab identification number of father \cr
#' \code{momID} \tab numeric \tab identification number of mother \cr
#' \code{sex} \tab numeric \tab gender identification; if male \code{sex = 0}, if female \code{sex = 1} \cr
#' \code{affected} \tab logical \tab disease-affection status, \code{affected  = TRUE} if affected by disease , and \code{FALSE} otherwise, \cr
#' \code{Gen} \tab numeric \tab the individual's generation number relative to the eldest founder. \cr
#' \tab \tab That is, for the eldest founder \code{Gen} = 1, for his or her offspring \code{Gen} = 2, etc. \cr
#' }
#'
#' @inheritParams censor_ped
#' @return \item{\code{reGen_ped} }{A pedigree containing only affected members, obligate carriers, and founders with generation number based on the most recent common ancestor of affected members as described in details.}
#' @export
#' @seealso \code{\link{new.ped}}
#'
#' @importFrom kinship2 kinship
#' @references OUR MANUSCRIPT
#' @references Thompson, E. (2000). \emph{Statistical Inference from Genetic Data on Pedigrees.} NSF-CBMS Regional Conference Series in Probability and Statistics, 6, I-169. Retrieved from http://www.jstor.org.proxy.lib.sfu.ca/stable/4153187
#'
#' @examples
#' # Read in example pedigrees
#' data(EgPeds)
#' class(EgPeds)
#'
#' # Create ped object
#' Bpeds <- new.ped(EgPeds)
#' summary(Bpeds)
#'
#' # Reassign generation numbers in the first four pedigrees in EgPeds
#' Apeds <- lapply(seq_len(4), function(x){
#'                  reassign_gen(Bpeds[Bpeds$FamID == x, ])})
#' Apeds <- do.call(rbind, Apeds)
#'
#' # Create kinship2 pedigree objects so that we can plot pedigrees
#' kin2ped_before <- ped2pedigree(Bpeds)
#' kin2ped_after  <- ped2pedigree(Apeds)
#'
#' # Compare pedigrees before and after reassigning
#' # generation number based on affected status
#' par(mfrow = c(1, 2))
#' for (k in 1:4) {
#'   ID1 = paste0("ID", sep = ":",
#'                Bpeds$ID[Bpeds$FamID == k],
#'                sep = "\n Gen:", Bpeds$Gen[Bpeds$FamID == k])
#'
#'   ID2 = paste0("ID", sep = ":",
#'                Apeds$ID[Apeds$FamID == k],
#'                sep = "\n Gen:", Apeds$Gen[Apeds$FamID == k])
#'
#'   plot(kin2ped_before[paste0(k)], id = ID1)
#'   mtext(paste0("Ped", k, ": before generation reassignment", sep = ""),
#'         side = 3)
#'
#'   plot(kin2ped_after[paste0(k)], id = ID2)
#'   mtext(paste0("Ped", k, ": after generation reassignment", sep = ""),
#'         side = 3)
#' }
#' par(mfrow = c(1, 1))
#'
#'
reassign_gen = function(ped_file){

  if (!is.ped(ped_file)) {
    stop("\n \n Expecting a ped object. \n Please use new.ped to create an object of class ped.")
  }

  if(!"Gen" %in% colnames(ped_file)){
    ped_file$Gen <- assign_gen(ped_file)
  }

  #create new ped file with available affecteds only
  reGen_ped <- ped_file[ped_file$affected & ped_file$available, ]

  if (nrow(reGen_ped) == 0) {
    warning("No affecteds to assign affected generation")
    return(reGen_ped)
  } else {
    d <- 0
    while (d == 0) {
      #find the dad IDs that are required but have been removed
      readd_dad <- find_missing_parent(reGen_ped)

      #find the mom IDs that are required but have been removed
      readd_mom <- find_missing_parent(reGen_ped, dad = FALSE)

      #check to see if we need to readd anyone
      if (length(c(readd_dad, readd_mom)) == 0) {
        d <- 1
      } else {
        #Now pull the rows containing the required parents
        # from the original ped_file
        readd <- ped_file[which(ped_file$ID %in% c(readd_dad, readd_mom)), ]

        #combine with affected ped file
        reGen_ped <- rbind(reGen_ped, readd)
      }
    }

    #Change Generation number so that only affecteds have a gen number
    reGen_ped$Gen <- ifelse(reGen_ped$affected, reGen_ped$Gen, NA)

    #table affected generation number
    Gen_tab <- table(reGen_ped$Gen)
    #minimum (earliest) generation number (i.e. 1)
    min_gen <- as.numeric(names(Gen_tab[1]))
    #second smallest generation number
    min_gen2 <- as.numeric(names(Gen_tab[2]))
    #difference between two earliest generation numbers
    gen_diff <- min_gen2-min_gen
    #number of affecteds in the earliest generation
    num_in_min_gen <- as.numeric(Gen_tab[1])


    if (min_gen == 1 | (min_gen == 2 & num_in_min_gen >= 2)) {
      #For this condition we do not need to reassign generation number
      # Covers case when the founder who introduced the RV is affected
      # and available, and case when RV founder not affected but 2 or more
      # of his or her children are affected.
      return(reGen_ped)
    } else if (num_in_min_gen >= 2) {
      #compute and store kinship matrix
      kin_mat <- kinship(reGen_ped,
                         id = reGen_ped$ID,
                         dadid = reGen_ped$dadID,
                         momid = reGen_ped$momID)

      #find the distance between only those in the lowest generation
      kin_distance <- -log(kin_mat[which(reGen_ped$Gen == min_gen),
                                  which(reGen_ped$Gen == min_gen)])/log(2)
      #find the new gen difference
      new_gen_diff <- min_gen - (max(kin_distance)/2 + 1)
      #assign new generation
      reGen_ped$Gen[!is.na(reGen_ped$Gen)] <- reGen_ped$Gen[!is.na(reGen_ped$Gen)] - new_gen_diff
      return(reGen_ped)
    } else if (num_in_min_gen == 1) {
      #compute and store kinship matrix
      kin_mat <- kinship(reGen_ped,
                         id = reGen_ped$ID,
                         dadid = reGen_ped$dadID,
                         momid = reGen_ped$momID)


      #find the distance between only those in the lowest 2 generations
      kin_distance <- -log(kin_mat[which(reGen_ped$Gen == min_gen),
                                  which(reGen_ped$Gen %in% c(min_gen,min_gen2))])/log(2)

      #find the difference between the maximum distance in kin_distance and
      # the value of the second smallest generation - this is the value we will
      # use to adjust everyone at or below the second smallest generation
      new_gen_diff <- min_gen2 - max(kin_distance)

      #find the difference between the maximum distance in kin_distance
      # and the difference between the smallest 2 generations
      # This will be the new gen no for the 1 individual in the lowest generation
      new_gen_oldest <- max(kin_distance) - gen_diff

      #assign new generation
      reGen_ped$Gen[!is.na(reGen_ped$Gen)] <- ifelse(reGen_ped$Gen[!is.na(reGen_ped$Gen)] == min_gen,
                                                    new_gen_oldest,
                                                    reGen_ped$Gen[!is.na(reGen_ped$Gen)] - new_gen_diff)

      return(reGen_ped)
    }
  }
}


#' Censor pedigree data
#'
#' \code{censor_ped} censors a pedigree of any information that occurs after a specified year.
#'
#' Upon supplying a pedigree and a censor year the \code{censor_ped} function will remove all individuals born after \code{censor_year} and censor all disease onset and death events after the \code{censor_year}.
#'
#' Users who wish to use \code{censor_ped} for pedigrees not generated by \code{\link{sim_ped}} or \code{\link{sim_RVped}} must use \code{\link{new.ped}} to create an object of class \code{ped}.  When creating the \code{ped} object please provide as much relevant date information as possible, i.e. years of birth, onset, and death.  When present please specify a proband as described in \code{\link{new.ped}}.
#'
#' By default, \code{censor_year} is set to the year that the pedigree is ascertained, i.e. the year the proband experienced disease onset. However, if \code{ped_file} does not contain the proband identification variable the user must supply a value for \code{censor_year}.
#'
#' @param ped_file An object of class \code{ped}. A pedigree generated by \code{sim_ped} or \code{sim_RVped}, or an object created by the function \code{\link{new.ped}}.  See details.
#' @param censor_year Numeric. The censor year. If not supplied, defaults to the year the pedigree was ascertained, i.e. the proband's onset year.  See details.
#'
#' @return \item{\code{censored_ped}}{The censored pedigree.}
#' @export
#' @seealso \code{\link{new.ped}}
#'
#' @examples
#' #Read in age-specific harard data and create hazard object.
#' data(AgeSpecific_Hazards)
#' haz_obj <- hazard(hazardDF = AgeSpecific_Hazards)
#'
#' #Simulate a pedigree ascertained for multiple affecteds
#' set.seed(3)
#' RVped2015 <- sim_RVped(hazard_rates = haz_obj,
#'                        num_affected = 2,
#'                        ascertain_span = c(1900, 2015),
#'                        GRR = 30, carrier_prob = 0.002,
#'                        RVfounder = TRUE,
#'                        stop_year = 2015,
#'                        recall_probs = c(1),
#'                        founder_byears = c(1900, 1905),
#'                        FamID = 1)[[2]]
#'
#' # Plot the 2015 pedigree
#' plot(RVped2015)
#' mtext(side = 3, line = 2, "Reference Year: 2017")
#'
#' # Censor RVped2015 after 1960
#' RVped1960 <- censor_ped(ped_file = RVped2015, censor_year = 1960)
#'
#' # Plot the 1960 pedigree
#' plot(RVped1960)
#' mtext(side = 3, line = 2, "Reference Year: 1960")
#'
censor_ped = function(ped_file, censor_year = NULL){

  if (!is.ped(ped_file)) {
    stop("\n \n Expecting a ped object. \n Please use new.ped to create an object of class ped.")
  }

  if (any(is.na(match(c("birthYr", "onsetYr", "deathYr"),
                      colnames(ped_file))))) {
    stop("\n \n Missing date data. \n Please ensure that ped_file includes the following variables:\n birthYr, onsetYr, deathYr" )
  }

  if (all(is.na(ped_file$birthYr))
      & all(is.na(ped_file$onsetYr))
      & all(is.na(ped_file$deathYr))) {
    stop("\n \n Nothing to censor, all date data is missing.")
  }

  if (is.null(censor_year)) {
    if ("proband" %in% colnames(ped_file)) {
      if(sum(ped_file$proband) == 1){
        if (is.na(ped_file$onsetYr[ped_file$proband])) {
          stop("\n \n Proband's onset year is missing. \n Specify the proband's onset year or specify censor_year.")
        } else {
          censor_year <- ped_file$onsetYr[ped_file$proband]
        }
      } else {
        stop("\n \n Proband cannot be uniquely identified.\n  Please identify a single proband or specify censor_year. \n ")
      }
    } else {
      stop("\n \n Proband cannot be identified. \n  Please identify a proband or specify censor_year.")
    }
  }


  #censor any onset or death info before censor year
  ped_file$affected <- ifelse(is.na(ped_file$onsetYr), FALSE,
                              ifelse(ped_file$onsetYr <= censor_year, T, F))
  ped_file$onsetYr <- ifelse(is.na(ped_file$onsetYr), NA,
                             ifelse(ped_file$onsetYr <= censor_year,
                                    ped_file$onsetYr, NA))
  ped_file$deathYr <- ifelse(is.na(ped_file$deathYr), NA,
                             ifelse(ped_file$deathYr <= censor_year,
                                    ped_file$deathYr, NA))

  # if proband exists remove proband status if the
  # proband does not experience onset before censor_year
  if ("proband" %in% colnames(ped_file)) {
    ped_file$proband[ped_file$proband] <- ifelse(is.na(ped_file$onsetYr[ped_file$proband]),
                                                 FALSE, ped_file$proband[ped_file$proband])
  }

  if (all(is.na(ped_file$birthYr))) {
    censored_ped <- ped_file
    warning("\n \n Birth data not detected. \n  Censoring onset and death data only")
  } else {

    #create new ped file containing only individuals born before the censor year
    censored_ped <- ped_file[which(ped_file$birthYr <= censor_year), ]

    if (nrow(censored_ped) == 0) {
      warning("\n \n No data recorded prior to the specified year.")
      return(censored_ped)
    } else {

      d <- 0
      while (d == 0) {
        #find the dad IDs that are required but have been removed
        readd_dad <- find_missing_parent(censored_ped)

        #find the mom IDs that are required but have been removed
        readd_mom <- find_missing_parent(censored_ped, dad = FALSE)

        #check to see if we need to readd anyone
        if (length(c(readd_dad, readd_mom)) == 0) {
          d <- 1
        } else {
          #Now pull the rows containing the required parents
          # from the original ped_file
          readd <- ped_file[which(ped_file$ID %in% c(readd_dad, readd_mom)), ]

          #combine with censored ped file
          censored_ped <- rbind(censored_ped, readd)
        }
      }
    }
  }

  return(censored_ped)
}
#' Gather information for the affected relatives
#'
#' @param ped_file ped object
#'
#' @return \item{\code{affVars}}{Information for the affected relatives.}
#' @keywords internal
#'
get_affectedInfo <- function(ped_file){

  aLoc <- match(c("available"), colnames(ped_file))

  if (is.na(aLoc)) {
    ped_file$available <- T
  }

  dA_loc <- match(c("DA1", "DA2"), colnames(ped_file))

  if (length(dA_loc[!is.na(dA_loc)]) == 2) {
    ped_file$RVstatus <- ped_file$DA1 + ped_file$DA2
  }

  keep_cols <- match(c("FamID", "ID",
                       "birthYr", "onsetYr", "deathYr",
                       "RR", "proband", "RVstatus"), colnames(ped_file))

  affected_info <- ped_file[ped_file$affected & ped_file$available,
                            keep_cols[!is.na(keep_cols)]]

  rownames(affected_info) <- NULL
  class(affected_info)    <- "data.frame"

  return(affected_info)
}

#' Compute kinship matrix for the affected relatives
#'
#' @param ped_file ped object
#'
#' @return \item{\code{affKin}}{The kinship matix for the affected relatives only.}
#' @keywords internal
#' @importFrom kinship2 kinship
#'
get_kinship <- function(ped_file){

  aLoc <- match(c("available"), colnames(ped_file))

  if (is.na(aLoc)) {
    ped_file$available <- T
  }

  kin_ped <- ped2pedigree(ped_file)

  kinMat <- kinship(kin_ped)[ped_file$affected & ped_file$available,
                             ped_file$affected & ped_file$available]

  return(kinMat)
}

#' Obtain summary information
#'
#' @param ped_file
#'
#' @return A data frame containing the relevant summary information
#' @keywords internal
get_famInfo <- function(ped_file){

  AV <- get_affectedInfo(ped_file)

  SRV <- ifelse(any(is.na(match(c("DA1", "DA2"), colnames(ped_file)))), NA,
                ifelse(any(ped_file$DA1 == 1) | any(ped_file$DA2 == 1), TRUE, FALSE))

  YRcols <- match(c("birthYr", "onsetYr"), colnames(ped_file))
  AOO <- ifelse(nrow(AV) > 0 & length(YRcols[!is.na(YRcols)]) == 2,
                  mean(AV$onsetYr - AV$birthYr, na.rm = T), NA)

  AY <- ifelse(nrow(AV) > 0
                 & !is.na(match(c("proband"), colnames(ped_file))) & sum(AV$proband) == 1,
                 AV$onsetYr[AV$proband], NA)

  AK <- get_kinship(ped_file)
  AIBD <- ifelse(nrow(AV) > 0, 2*mean(AK[upper.tri(AK)]), NA)

  return(data.frame(FamID = ped_file$FamID[1],
                    totalRelatives = nrow(ped_file),
                    numAffected = nrow(AV),
                    aveOnsetAge = AOO,
                    aveIBD = AIBD,
                    ascertainYear = AY,
                    segRV = SRV))
}
