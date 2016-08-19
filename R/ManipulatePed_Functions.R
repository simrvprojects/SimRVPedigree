#' Assign minimum carrier generation to a simulated pedigree
#'
#' \code{AssignGeneration} reassigns the generation number based on affected status and trims the pedigree down to affecteds and obligate carriers.
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
#' part_vec <- seq(0, 100, by = 1)
#' unaffected_mort <- 0.00001 + pgamma(seq(0.16, 16, by = .16),
#'                                     shape = 9.5, scale = 1)/350
#' affected_mort <- c(0.55, 0.48, 0.37, 0.23, 0.15,
#'                    pgamma(seq(0.96, 16, by = .16),
#'                    shape = 4, scale = 1.5))/300
#' Dhaz_df  <- (as.data.frame(cbind(unaffected_mort, affected_mort)))
#' Ohaz_vec <- (dgamma(seq(0.1, 10, by = .1), shape = 8, scale = 0.75))/60
#'
#' set.seed(1)
#' ex_RVped <- sim_RVpedigree(onset_hazard = Ohaz_vec,
#'                            death_hazard = Dhaz_df,
#'                            part = part_vec, RR = 5, FamID = 1,
#'                            founder_byears = c(1900, 1910),
#'                            ascertain_span = c(1900, 2015),
#'                            num_affected = 2)
#'
#' plot_RVpedigree(ex_RVped[[1]])
#' plot_RVpedigree(ex_RVped[[2]])

assignGeneration = function(ped_file){

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



