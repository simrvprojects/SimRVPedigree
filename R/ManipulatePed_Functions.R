#' Assign affected generation number.
#'
#' \code{assign_affectedGen} reassigns the generation number based on affected status and reduces pedigree to affecteds and obligate carriers.
#'
#' The \code{assign_affectedGen} function accepts a ped file simulated by \code{sim_RVpedigree} and reassigns generation number based on the affecteds in the pedigree.  Specifically, given a pedigree it reassigns generation numbers of the affected memebers so that generation 1 represents the first generation containing an affected member or a suspected carrier of a genetic factor assumed to increase disease susceptibility.
#'
#' #' For user who are reassign generation number based on affection status in pedigrees that have not been simulated by \code{sim_RVpedigree} or \code{sim_ped}, the \code{ped_file} supplied to \code{assign_affectedGen} must contain the following variables for each pedigree member:
#'
#' \enumerate{
#' \item \code{ID}: an identification number.
#' \item \code{dad_id}: identification number of father.
#' \item \code{mom_id}: identification number of mother.
#' \item \code{gender}: gender identification; if male \code{gender = 0}, if female \code{gender = 1}.
#' \item \code{affected}: affection status, if affected by disease \code{affected  = 1}, otherwise, \code{affected = 0}.
#' \item \code{birth_year}: the individual's birth year.
#' \item \code{onset_year}: the individual's disease onset year, when applicable.
#' \item \code{death_year}: the individual's death year, when applicable.
#' \item \code{Gen}: the individual's current generation number.
#' }
#'
#'
#' @param ped_file data.frame. A pedigree to reassign generation number based on affection status, see details.
#'
#' @return reGen_ped A ped file with only affecteds and obligate carriers, with generation number assigned based on affected status.
#' @export
#'
#' @importFrom kinship2 kinship
#'
#' @examples
#' #Read in example pedigree to trim
#' data(exp_peds)
#'
#' library(kinship2)
#' #assign to pedigree object to show before and after behavior of
#' #the assign_affectedGen function
#' ex_pedigree <- pedigree(id = exp_peds$ID,
#'                         dadid = exp_peds$dad_id,
#'                         momid = exp_peds$mom_id,
#'                         sex = (exp_peds$gender + 1),
#'                         affected = exp_peds$affected,
#'                         famid = exp_peds$FamID)
#'
#'
#' #create df to store peds with reassigned generation number
#' RAG_peds <- exp_peds[1,]
#' RAG_peds <- RAG_peds[-1,]
#'
#' for(i in 1:4){
#'   RAG_peds <- rbind(RAG_peds,
#'                     assign_affectedGen(exp_peds[which(exp_peds$FamID == i), ]))
#' }
#'
#' RAG_pedigrees <-  pedigree(id = RAG_peds$ID,
#'                            dadid = RAG_peds$dad_id,
#'                            momid = RAG_peds$mom_id,
#'                            sex = (RAG_peds$gender + 1),
#'                            affected = RAG_peds$affected,
#'                            famid = RAG_peds$FamID)
#'
#' #Affecteds in Pedigrees 1 and 4 will keep their original generation
#' # assignment, while affecteds in pedigrees 2 and 3 will be  given new
#' # generation numbers so that generation 1 represents the generation of
#' # first affected or first obligate carrier
#' par(mfrow = c(1, 2))
#' for (k in 1:4) {
#'   ID1 = paste0("ID", sep = ":",
#'                exp_peds[which(exp_peds$FamID == k), 2],
#'                sep = "\n Gen:", exp_peds[which(exp_peds$FamID == k), 14])
#'   ID2 = paste0("ID", sep = ":",
#'                RAG_peds[which(RAG_peds$FamID == k), 2],
#'                sep = "\n Gen:", RAG_peds[which(RAG_peds$FamID == k), 14])
#'   plot(ex_pedigree[paste0(k)], id = ID1)
#'   mtext(paste0("Ped", k, ": before generation reassignment", sep = ""),
#'         side = 3)
#'   plot(RAG_pedigrees[paste0(k)], id = ID2)
#'   mtext(paste0("Ped", k, ": after generation reassignment", sep = ""),
#'         side = 3)
#' }
#'
#' #FamID2
#' #since there are at least two affected in generation 2, one of the
#' #founders is an obligate carrier and generation numbers for affecteds
#' #do not need to be re-assigned
#'
#' #NOTE: this pedigree now only includes affected individuals and
#' #individuals who are required to create a pedigree.
#'
#' #FamID3
#' #assign_affectedGen will reassign the generation for individual 5 to 2 and
#' #the generation for individual 9 to 3.
#'
#'
#' #FamID4
#' # reassign generation numbers for individuals 3, 5, and 10,
#' # to 1, 2, and, 2 respectively.
#'
#' #FamID5
#' #No reassignment of generation number required for these affecteds
#'
assign_affectedGen = function(ped_file){

  #create new ped file with affecteds only
  reGen_ped <- ped_file[which(ped_file$affected == 1), ]

  if (nrow(reGen_ped) == 0) {
    warning("No affecteds to assign affected generation")
    return(reGen_ped)
  } else {
    d <- 0
    while (d == 0) {
      #find the dad IDs that are required but have been removed
      miss_dad  <- !is.element(reGen_ped$dad_id,
                               reGen_ped$ID[which(reGen_ped$gender == 0)])
      readd_dad <- reGen_ped$dad_id[miss_dad]
      readd_dad <- unique(readd_dad[!is.na(readd_dad)])

      #find the mom IDs that are required but have been removed
      miss_mom  <- !is.element(reGen_ped$mom_id,
                               reGen_ped$ID[which(reGen_ped$gender == 1)])
      readd_mom <- reGen_ped$mom_id[miss_mom]
      readd_mom <- unique(readd_mom[!is.na(readd_mom)])

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
    reGen_ped$Gen <- ifelse((reGen_ped$affected == 1),
                            reGen_ped$Gen, NA)

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
                         dadid = reGen_ped$dad_id,
                         momid = reGen_ped$mom_id)

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
                         dadid = reGen_ped$dad_id,
                         momid = reGen_ped$mom_id)


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


#' Censor pedigree information after a specified year.
#'
#' The \code{censor_ped} function censors a pedigree relative to a specified year.
#'
#' Upon supplying a pedigree and a censor year the \code{censor_ped} function will remove all individuals born after the censor year and censor all disease onset and death events after the censor year.
#'
#' For user who are censoring pedigrees that have not been simulated by \code{sim_RVpedigree} or \code{sim_ped}, the \code{ped_file} supplied to \code{censor_ped} must contain the following variables for each pedigree member:
#'
#' \enumerate{
#' \item \code{ID}: an identification number.
#' \item \code{dad_id}: identification number of father.
#' \item \code{mom_id}: identification number of mother.
#' \item \code{gender}: gender identification; if male \code{gender = 0}, if female \code{gender = 1}.
#' \item \code{affected}: affection status, if affected by disease \code{affected  = 1}, otherwise, \code{affected = 0}.
#' \item \code{birth_year}: the individual's birth year.
#' \item \code{onset_year}: the individual's disease onset year, when applicable.
#' \item \code{death_year}: the individual's death year, when applicable.
#' }
#'
#' If an individual has not experienced disease onset and/or death, then \code{onset_year = NA} and/or \code{death_year = NA}.
#'
#' If \code{censor_year} is missing, when pedigree contains a proband, \code{censor_year} will assume the value of the proband's onset year. However, if \code{ped_file} does not contain a proband identification variable the user must supply a value for \code{censor_year}.
#'
#'
#'
#' @param ped_file data.frame. The pedigree to censor, see details.
#' @param censor_year Numeric. The censor year. If missing, when pedigree contains a proband, \code{censor_year} will assume the value of the proband's onset year. See details.
#'
#' @return censored_ped The censored pedigree.
#' @export
#'
#' @examples
#' #Read in example pedigree to trim
#' data(AgeSpecific_Hazards)
#'
#'
#' #Specify onset hazard and death hazard
#' my_onset_hazard <- AgeSpecific_Hazards[,1]
#' my_death_hazard <- AgeSpecific_Hazards[,c(2,3)]
#'
#' #Specify age partition.
#' age_part <- seq(0, 100, by = 1)
#'
#'
#' #Simulate pedigree
#' set.seed(3)
#' ex_RVped <- sim_RVpedigree(onset_hazard = my_onset_hazard,
#'                            death_hazard = my_death_hazard,
#'                            part = age_part,
#'                            num_affected = 2,
#'                            ascertain_span = c(1900, 2015),
#'                            RR = 10, stop_year = 2015,
#'                            recall_probs = c(1),
#'                            founder_byears = c(1900, 1980),
#'                            FamID = 1)
#'
#' library(kinship2)
#' Original_ped <- with(ex_RVped[[2]], pedigree(id = ID,
#'                                              dadid = dad_id,
#'                                              momid = mom_id,
#'                                              sex = gender + 1,
#'                                              affected = affected))
#' plot(Original_ped)
#'
#' Censor_ped1 <- with(censor_ped(ex_RVped[[2]]),
#'                     pedigree(id = ID,
#'                              dadid = dad_id,
#'                              momid = mom_id,
#'                              sex = gender + 1,
#'                              affected = affected))
#' plot(Censor_ped1)
#'
#' Censor_ped2 <- with(censor_ped(ex_RVped[[2]], 2000),
#'                     pedigree(id = ID,
#'                              dadid = dad_id,
#'                              momid = mom_id,
#'                              sex = gender + 1,
#'                              affected = affected))
#' plot(Censor_ped2)
#'
censor_ped = function(ped_file, censor_year){

  if (missing(censor_year)) {
    if ("proband" %in% colnames(ped_file)) {
      censor_year <- ped_file$onset_year[which(ped_file$proband == 1)]
    } else {
      stop("Ped file must contain a proband or user must supply censor_year")
    }
  }

  #create new ped file containing only individuals born before the censor year
  censored_ped <- ped_file[which(ped_file$birth_year <= censor_year), ]

  if (nrow(censored_ped) == 0) {
    warning("Please check censor_year, no pedigree information prior to censor_year")
    return(censored_ped)
  } else {
    #censor onset and death events prior to censor year
    censored_ped$affected <- ifelse(is.na(censored_ped$onset_year), 0,
                                    ifelse(censored_ped$onset_year <= censor_year,
                                           censored_ped$affected, 0))
    censored_ped$onset_year <- ifelse(censored_ped$onset_year <= censor_year,
                                      censored_ped$onset_year, NA)
    censored_ped$death_year <- ifelse(censored_ped$death_year <= censor_year,
                                      censored_ped$death_year, NA)

    d <- 0
    while (d == 0) {
      #find the dad IDs that are required but have been removed
      miss_dad  <- !is.element(censored_ped$dad_id,
                               censored_ped$ID[which(censored_ped$gender == 0)])
      readd_dad <- censored_ped$dad_id[miss_dad]
      readd_dad <- unique(readd_dad[!is.na(readd_dad)])

      #find the mom IDs that are required but have been removed
      miss_mom  <- !is.element(censored_ped$mom_id,
                               censored_ped$ID[which(censored_ped$gender == 1)])
      readd_mom <- censored_ped$mom_id[miss_mom]
      readd_mom <- unique(readd_mom[!is.na(readd_mom)])

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

  return(censored_ped)
}
