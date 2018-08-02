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
#' @export
#'
find_available_parent <- function(ped, ID){
  dad <- ped$dadID[ped$ID == ID]
  mom <- ped$momID[ped$ID == ID]
  available_parent <- ifelse(is.na(ped$dadID[ped$ID == dad]), mom, dad)
  return(available_parent)
}


#' Find the most recent common ancestor of two related individuals
#'
#' @param ped The ped object
#' @param ID1 The ID of the first relative
#' @param ID2 The ID of the second relative
#'
#' @return The ID of the common ancestor
#' @export
#'
find_mrca <- function(ped, ID1, ID2){
  al1 <- c(ID1)
  al2 <- c(ID2)

  na1 <- ID1
  na2 <- ID2
  mrca_found <- FALSE
  while (mrca_found == FALSE) {
    al1 <- c(al1, find_available_parent(ped, na1))
    al2 <- c(al2, find_available_parent(ped, na2))

    if (length(intersect(al1, al2)) > 0) {
      common_ancestor <- intersect(al1, al2)
      mrca_found <- TRUE
    } else {
      na1 <- find_available_parent(ped, na1)
      na2 <- find_available_parent(ped, na2)

    }

  }

  return(common_ancestor)
}
