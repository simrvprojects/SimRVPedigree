#' Assign minimum carrier generation to a simulated pedigree
#'
#' \code{assign_affectedGen} reassigns the generation number based on affected status and trims the pedigree down to affecteds and obligate carriers.
#'
#' Add details concerning reassignment...
#'
#' @param ped_file data.frame. A ped file with the same format as one simulated with sim_RVpedigree
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

        #remove all of their simulated data, and mark unavailable
        readd$birth_year <- NA
        readd$onset_year <- NA
        readd$death_year <- NA
        readd$available  <- 0

        #combine with affected ped file
        reGen_ped <- rbind(reGen_ped, readd)
      }
    }

    #Change Generation number so that only affecteds have a gen number
    reGen_ped$Gen <- ifelse((reGen_ped$affected == 1 &
                               reGen_ped$available == 1),
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



