#' Find parents that have been removed but are required for plotting.
#'
#' @param x An object of class ped.
#' @param dad Logical. If TRUE looking for missing dad, if FALSE looking for missing mom.
#'
#' @return The IDs of the missing parents.
#' @keywords internal
find_missing_parent <- function(x, dad = TRUE){
  if (dad) {
    miss_parent  <- !is.element(x$dadID, x$ID[which(x$sex == 0)])
    readd_parent <- x$dadID[miss_parent]
  } else {
    miss_parent  <- !is.element(x$momID, x$ID[which(x$sex == 1)])
    readd_parent <- x$momID[miss_parent]
  }
  return(unique(readd_parent[!is.na(readd_parent)]))
}

#' Assign generation number based on oldest founder
#'
#' @param x an object of class ped
#'
#' @return a list of generation numbers for pedigree members, in the order listed in \code{x}.
#' @importFrom kinship2 kindepth
#' @keywords internal
assign_gen <- function(x){
  Gen <- NA
  mates <- cbind(x$dadID, x$momID)
  #remove all rows with only zeros, these are founders
  mates <- unique(mates)
  mates <- mates[!is.na(mates[, 1]), ]

  #get kindepth and set Gen to kindepth when kindepth is non-zero
  kd <- kindepth(ped2pedigree(x))
  Gen[kd != 0] <- kd[kd != 0]

  if(class(mates) == "matrix"){
    for(i in 1:nrow(mates)){
      mate_gens <-  kd[x$ID %in% mates[i, ]]
      Gen[x$ID %in% mates[i, ]] <- max(mate_gens)
    }
  } else {
    Gen[x$ID %in% mates] <- 0
  }

  return(Gen + 1)
}


#' Get generation lables
#'
#' @param x An object of class ped.
#' @param nlevel The total number of levels in the pedigree.
#'
#' @return A character list of either blanks or roman numerals which wil be plotted along side the pedigree.
#' @keywords internal
get_gen_labs <- function(x, nlevel){
  max_gen <- max(x$Gen, na.rm = T)
  if (nlevel == max_gen) {
    return(as.character(as.roman(seq(1:max_gen))))
  } else {
    return(c(rep(" ", nlevel - max_gen),
             as.character(as.roman(seq(1:max_gen)))))
  }
}


#' Find available parent
#'
#' @param ped the ped object
#' @param ID the Id of the person whose availble parent we would like to find
#'
#' @return the ID of the availble parent
#' @keywords internal
#'
find_available_parent <- function(ped, ID){
  if (is.na(ID)) {
    #parent of missing person is also missing
    return(c())
  } else {
    dad <- ped$dadID[ped$ID == ID]
    mom <- ped$momID[ped$ID == ID]
    if (is.na(dad) & is.na(mom)) {
      #if mom and dad are both missing, then this person is a founder so return nothing
      available_parent <- c()
    } else {
      available_parent <- c(mom, dad)
    }

    return(available_parent)
  }

}


#' Find the most recent common ancestor of two pedigree members
#'
#' @param ped A ped object
#' @param ID1 The ID of the first relative
#' @param ID2 The ID of the second relative
#' @importFrom kinship2 kinship
#'
#' @return The ID of the common ancestor
#' @export
#'
#' @examples
#' library(SimRVPedigree)
#' data(AgeSpecific_Hazards)
#'
#'
#' set.seed(5)
#' ex_ped <- sim_ped(hazard_rates = hazard(hazardDF = AgeSpecific_Hazards),
#'                   GRR = 10, FamID = 1,
#'                   founder_byears = c(1800, 1900),
#'                   random_BR = TRUE,
#'                   stop_year = 2020)
#'
#' plot(ex_ped)
#'
#' # Find most recent common ancestor of individuals with IDs 19 and 21
#' find_mrca(ped = ex_ped, ID1 = 19, ID2 = 21)
#'
#' # Note that someone can be their own most recent common ancsetor.
#' # In the following example, since the individual with ID 8 is the grandmother
#' # of the individual with ID 21, the find_mrca function returns 8.
#' find_mrca(ped = ex_ped, ID1 = 8, ID2 = 21)
#'
#' # For unrelated individuals, the find_mcra function returns NA
#' find_mrca(ped = ex_ped, ID1 = 8, ID2 = 15)
#' find_mrca(ped = ex_ped, ID1 = 5, ID2 = 4)
#'
find_mrca <- function(ped, ID1, ID2){

  if (!is.ped(ped)) {
    stop("\n \n Expecting a ped object. \n Please use new.ped to create an object of class ped.")
  }

  #check to see if the individuals are actually related
  if (kinship(ped2pedigree(ped))[ped$ID == ID1, ped$ID == ID2] == 0) {
    common_ancestor <- NA
  } else {
    #the following lists will store all of the ancestors of the individuals
    # with ID1 and ID2, **including themselves**.
    ancestors_ID1 <- c(ID1)
    ancestors_ID2 <- c(ID2)

    # the following variables are the people for whom we want to
    # find parents.  We update these in the following while loop.
    find_parent_of_1 <- ID1
    find_parent_of_2 <- ID2
    mrca_found <- FALSE
    while (mrca_found == FALSE) {
      if (length(find_parent_of_1) == 1) {
        ancestors_ID1 <- c(ancestors_ID1, find_available_parent(ped, find_parent_of_1))
      } else {
        ancestors_ID1 <- c(ancestors_ID1,
                           unlist(lapply(find_parent_of_1, function(x){find_available_parent(ped, x)})))
      }

      if (length(find_parent_of_2) == 1) {
        ancestors_ID2 <- c(ancestors_ID2, find_available_parent(ped, find_parent_of_2))
      } else {
        ancestors_ID2 <- c(ancestors_ID2,
                           unlist(lapply(find_parent_of_2, function(x){find_available_parent(ped, x)})))
      }

      if (length(intersect(ancestors_ID1, ancestors_ID2)) > 0) {
        #common ancestor is found!
        common_ancestor <- intersect(ancestors_ID1, ancestors_ID2)
        mrca_found <- TRUE
      } else {
        #no common ancestor, update variables to
        #the individuals we just added to ancestor lists]
        if (length(find_parent_of_1) == 1) {
          find_parent_of_1 <- find_available_parent(ped, find_parent_of_1)
        } else {
          find_parent_of_1 <- unlist(lapply(find_parent_of_1, function(x){find_available_parent(ped, x)}))
        }

        if (length(find_parent_of_2) == 1) {
          find_parent_of_2 <- find_available_parent(ped, find_parent_of_2)
        } else {
          find_parent_of_2 <- unlist(lapply(find_parent_of_2, function(x){find_available_parent(ped, x)}))
        }

      }

    }
  }

  return(common_ancestor)
}
