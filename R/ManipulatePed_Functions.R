#' Assign Generation Number Based on Affected Status
#'
#' \code{assign_affectedGen} assigns generation numbers among affected family members so that generation 1 represents the most recent generation that a putative disease variant shared identical by descent (IBD) [add citation] by affected members could have been introduced into the pedigree.
#'
#' The \code{assign_affectedGen} function accepts a pedigree simulated by \code{sim_RVped} and reassigns generation numbers among affected family members in the pedigree.  Specifically, given a pedigree this function reassigns the generation numbers of affected members so that generation 1 is assigned to the most recent common ancestor of all affected members.  We note that the individual in generation 1 could themselves be affected, i.e. an individual can be considered their own ancestor.
#'
#' For example, consider a family with 2 affected members.  If the two affected members are a parent and his or her offspring, the affected parent would be assigned generation 1, and the affected child generation 2.  However, if the two affected members are a pair of siblings, each affected sibling would be assigned generation 2 since a common parent of the two affected siblings is assumed to be a carrier of a latent susceptibility variant.  Similarly, if the two affected members are a pair of cousins, each affected cousin is assigned generation 3, since a common grandparent of the two affected cousins is the most recent common ancestor from whom they could have inherited a shared variant associated with the disease. [cite our paper for this example]
#'
#' Users who wish to assign generation number based on affection status in pedigrees that have not been simulated with the \code{SimRVpedigree} package must ensure that the pedigree, \code{ped_file}, supplied to \code{assign_affectedGen} contains the following variables for each pedigree member:
#'
#' \enumerate{
#' \item \code{ID}: an identification number.
#' \item \code{dad_id}: identification number of father.
#' \item \code{mom_id}: identification number of mother.
#' \item \code{gender}: gender identification; if male \code{gender = 0}, if female \code{gender = 1}.
#' \item \code{affected}: affection status, if affected by disease \code{affected  = 1}, otherwise, \code{affected = 0}.
#' \item \code{birth_year}: the individual's year of birth.
#' \item \code{onset_year}: the individual's disease onset year, when applicable.
#' \item \code{death_year}: the individual's death year, when applicable.
#' \item \code{Gen}: the individual's generation number relative to the eldest founder.  For the eldest founder \code{Gen = 1}, for his or her offspring \code{Gen = 2}, etc.
#' }
#'
#'
#' @param ped_file data.frame. A pedigree to reassign generation number based on affection status, see details.
#'
#' @return \code{reGen_ped} A pedigree containing only affected members, obligate carriers, and founders with generation number based on the most recent common ancestor of affected members as described in details.
#' @export
#'
#' @importFrom kinship2 kinship
#'
#' @examples
#' #Read in example pedigrees
#' data(ExamplePedigrees)
#'
#' library(kinship2)
#' #assign to pedigree object to show before and after behavior of
#' #the assign_affectedGen function
#' ex_pedigree <- pedigree(id = ExamplePedigrees$ID,
#'                         dadid = ExamplePedigrees$dad_id,
#'                         momid = ExamplePedigrees$mom_id,
#'                         sex = (ExamplePedigrees$gender + 1),
#'                         affected = ExamplePedigrees$affected,
#'                         famid = ExamplePedigrees$FamID)
#'
#'
#' #create df to store pedigrees with reassigned generation number
#' RAG_peds <- ExamplePedigrees[1,]
#' RAG_peds <- RAG_peds[-1,]
#'
#' for(i in 1:4){
#'   RAG_peds <- rbind(RAG_peds,
#'                     assign_affectedGen(ExamplePedigrees[which(ExamplePedigrees$FamID == i), ]))
#' }
#'
#' RAG_pedigrees <-  pedigree(id = RAG_peds$ID,
#'                            dadid = RAG_peds$dad_id,
#'                            momid = RAG_peds$mom_id,
#'                            sex = (RAG_peds$gender + 1),
#'                            affected = RAG_peds$affected,
#'                            famid = RAG_peds$FamID)
#'
#' # Compare pedigrees before and after reassigning
#' # generation number based on affcted status
#' par(mfrow = c(1, 2))
#' for (k in 1:4) {
#'   ID1 = paste0("ID", sep = ":",
#'                ExamplePedigrees[which(ExamplePedigrees$FamID == k), 2],
#'                sep = "\n Gen:", ExamplePedigrees[which(ExamplePedigrees$FamID == k), 14])
#'
#'   ID2 = paste0("ID", sep = ":",
#'                RAG_peds[which(RAG_peds$FamID == k), 2],
#'                sep = "\n Gen:", RAG_peds[which(RAG_peds$FamID == k), 14])
#'
#'   plot(ex_pedigree[paste0(k)], id = ID1)
#'   mtext(paste0("Ped", k, ": before generation reassignment", sep = ""),
#'         side = 3)
#'
#'   plot(RAG_pedigrees[paste0(k)], id = ID2)
#'   mtext(paste0("Ped", k, ": after generation reassignment", sep = ""),
#'         side = 3)
#' }
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


#' Censor Pedigree After a Specified Year
#'
#' \code{censor_ped} censors a pedigree of any information that occurs after a specified year.
#'
#' Upon supplying a pedigree and a censor year the \code{censor_ped} function will remove all individuals born after \code{censor_year} and censor all disease onset and death events after the \code{censor_year}.
#'
#' Users who wish to censor pedigrees which have not been simulated by \code{\link{sim_RVped}} or \code{\link{sim_ped}} must ensure that the pedigree, \code{ped_file}, supplied to \code{censor_ped} contains the following variables for each pedigree member:
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
#' \item \code{proband}: (Optional) Proband identification variable, \code{proband = 1} if the individual is the proband, and 0 otherwise.
#' }
#'
#' If an individual has not experienced disease onset and/or death, then \code{onset_year = NA} and/or \code{death_year = NA}.
#'
#' If \code{censor_year} is missing, when the pedigree contains a proband, \code{censor_year} is set, internally, to the year that the proband experienced disease onset. However, if \code{ped_file} does not contain the proband identification variable the user must supply a value for \code{censor_year}.
#'
#' For a detailed example please refer to vignette.
#'
#' @param ped_file Data.frame. The pedigree to censor, see details.
#' @param censor_year Numeric. The censor year. If missing, when pedigree contains a proband, \code{censor_year} will assume the value of the proband's onset year. See details.
#'
#' @return censored_ped The censored pedigree.
#' @export
#'
#' @examples
#' #Read in example pedigree to trim
#' data(AgeSpecific_Hazards)
#'
#' #Simulate a pedigree ascertained for multiple affecteds
#' set.seed(3)
#' ex_RVped <- sim_RVped(onset_hazard = AgeSpecific_Hazards[,1],
#'                       death_hazard = AgeSpecific_Hazards[,c(2,3)],
#'                       part = seq(0, 100, by = 1),
#'                       num_affected = 2,
#'                       ascertain_span = c(1900, 2015),
#'                       RR = 10, stop_year = 2015,
#'                       recall_probs = c(1),
#'                       founder_byears = c(1900, 1980),
#'                       FamID = 1)[[2]]
#'
#' #To plot pedigrees, use the kinship2 package
#' library(kinship2)
#'
#' #assign simulated pedigree to pedigree object, then pass to plot function
#' Original_ped <- pedigree(id = ex_RVped$ID,
#'                          dadid = ex_RVped$dad_id,
#'                          momid = ex_RVped$mom_id,
#'                          sex = ex_RVped$gender + 1,
#'                          affected = ex_RVped$affected)
#' plot(Original_ped)
#'
#' #Now censor ex_RVped after 1999
#' Cped <- censor_ped(ped_file = ex_RVped, censor_year = 1999)
#'
#' #Assign the the censored pedigree to a pedigree object and plot.
#' Censor_ped <- pedigree(id = Cped$ID,
#'                        dadid = Cped$dad_id,
#'                        momid = Cped$mom_id,
#'                        sex = Cped$gender + 1,
#'                        affected = Cped$affected)
#' plot(Censor_ped)
#'
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
