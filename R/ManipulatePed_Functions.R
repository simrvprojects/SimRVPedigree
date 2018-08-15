#' Reassign generation number based on affected status
#'
#' The \code{reassign_gen} function assigns generation numbers among affected family members so that generation 1 represents the most recent generation that a putative disease variant shared identical by descent (IBD), as defined in Thompson (2000), by affected members could have been introduced into the pedigree.
#'
#' \emph{**\code{reassign_gen} cannot be applied to pedigrees that contain loops or inbreeding.**}
#'
#' The \code{reassign_gen} function accepts a pedigree and reassigns generation numbers among disease-affected relatives so that generation 1 represents the generation of the most recent common ancestor of all disease-affected relatives.  We note that the individual in generation 1 could themselves be disease-affected, i.e. an individual can be considered their own ancestor.
#'
#' For example, consider a family with 2 affected members.  If the disease-affected relatives are a parent and a child, the affected parent would be assigned generation 1, and the affected child generation 2.  However, if the disease-affected relatives are a pair of siblings, each is be assigned generation 2 since a common parent of the two is assumed to be a carrier of a latent susceptibility variant.  Similarly, if the disease-affected relatives are a pair of cousins, is assigned generation 3, since a common grandparent is the most recent common ancestor from whom they could have inherited a shared variant associated with the disease.
#'
#' Users who wish to assign generation number based on affection status in pedigrees that have not been simulated with the \code{SimRVpedigree} package must create a ped object using \code{\link{new.ped}}.
#'
#' @inheritParams censor_ped
#' @return A \code{ped} object containing only affected members, obligate carriers, and founders with generation numbers reassigned among disease-affected relatives based on their most recent common ancestor, as described in details.
#' @export
#' @seealso \code{\link{new.ped}}
#'
#' @importFrom kinship2 kinship
#' @importFrom kinship2 align.pedigree
#' @importFrom utils combn
#'
#' @references Nieuwoudt, Christina and Jones, Samantha J and Brooks-Wilson, Angela and Graham, Jinko. (14 December 2017) \emph{Simulating Pedigrees Ascertained for Multiple Disease-Affected Relatives}. bioRxiv 234153.
#' @references Thompson, E. (2000). \emph{Statistical Inference from Genetic Data on Pedigrees.} NSF-CBMS Regional Conference Series in Probability and Statistics, 6, I-169.
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
#' Apeds <- lapply(seq_len(5), function(x){
#'                  reassign_gen(Bpeds[Bpeds$FamID == x, ])})
#' Apeds <- do.call(rbind, Apeds)
#'
#' # Compare pedigrees before and after reassigning
#' # generation number based on affected status
#' par(mfrow = c(1, 2))
#' for (k in 1:5) {
#'   plot(subset(Bpeds, FamID == k), gen_lab = TRUE, plot_legend = FALSE)
#'   mtext(paste0("Ped", k, ": before generation reassignment", sep = ""),
#'         side = 3, line = 1.5)
#'
#'   plot(subset(Apeds, FamID == k), gen_lab = TRUE, plot_legend = FALSE)
#'   mtext(paste0("Ped", k, ": after generation reassignment", sep = ""),
#'         side = 3, line = 1.5)
#' }
#' par(mfrow = c(1, 1))
reassign_gen = function(ped_file){

  if (!is.ped(ped_file)) {
    stop("\n \n Expecting a ped object. \n Please use new.ped to create an object of class ped.")
  }

  if(!"Gen" %in% colnames(ped_file)){
    ped_file$Gen <- assign_gen(ped_file)
  }

  if (sum(ped_file$affected[ped_file$available]) == 0) {
    stop("\n \n No disease-affected relatives present.  \n Cannot reassign generations.")
  }

  if (sum(ped_file$affected[ped_file$available]) == 1) {
    gped <- ped_file[ped_file$affected & ped_file$available, ]
    gped$Gen <- 1
    warning(paste0("\n \n Family ", gped$FamID, " only contains one disease-affected relative."))
    return(gped)
  }

  #create new ped file with available affecteds only
  gped <- ped_file[ped_file$affected & ped_file$available, ]

  #add back any unaffected relatives needed to create a full pedigree
  d <- 0
  while (d == 0) {
    #find the dad IDs that are required but have been removed
    readd_dad <- find_missing_parent(gped)

    #find the mom IDs that are required but have been removed
    readd_mom <- find_missing_parent(gped, dad = FALSE)

    #check to see if we need to readd anyone
    if (length(c(readd_dad, readd_mom)) == 0) {
      d <- 1
    } else {
      #Now pull the rows containing the required parents
      # from the original ped_file

      readd <- ped_file[which(ped_file$ID %in% c(readd_dad, readd_mom)), ]

      #combine with affected ped file
      gped <- rbind(gped, readd)
    }
  }

  pedgre <- ped2pedigree(gped)
  if (any(align.pedigree(pedgre)$spouse == 2)) {
    stop("\n \n Inbreeding detected. \n reassign_gen is not intended for pedigrees that contain loops or inbreeding.")
  }

  id_array <- as.numeric(align.pedigree(pedgre)$nid)
  id_array <- id_array[id_array != 0]
  if (any(duplicated(id_array))) {
    stop("\n \n Loop detected. \n reassign_gen is not intended for pedigrees that contain loops or inbreeding.")
  }

  #compute and store kinship matrix
  kin_mat <- kinship(pedgre)

  #check to make sure all available affecteds are related
  if (any(kin_mat[gped$affected & gped$available,
                  gped$affected & gped$available] == 0)) {

    #Reduce to kinship mat for affecteds only
    ak_mat <- kin_mat[gped$affected & gped$available,
                      gped$affected & gped$available]

    url1 <- colnames(ak_mat)[which(ak_mat == 0, arr.ind = T)[1, 1]]
    url2 <- colnames(ak_mat)[which(ak_mat == 0, arr.ind = T)[1, 2]]

    stop(paste0("\n \n The disease-affected relatives with ID numbers ", url1, " and ", url2, " are not related. \n Cannot determine most recent common ancestor for unrelated affecteds."))
  }

  #store the old generation numbers
  old_gen <- gped$Gen

  #Change Generation number so that only affecteds have a gen number
  gped$Gen <- ifelse(gped$affected, gped$Gen, NA)

  #table affected generation number
  gen_tab <- table(gped$Gen)

  #first two smallest generation numbers (i.e. 1 and 2)
  min_gens <- as.numeric(names(gen_tab[c(1, 2)]))

  #number of affecteds in the earliest generation
  num_in_min_gen <- as.numeric(gen_tab[1])


  if (min_gens[1] == 1 | (min_gens[1] == 2 & num_in_min_gen >= 2)) {
    # For this condition we do not need to reassign generation number
    # Covers case when the founder who introduced the RV is affected
    # and available, and case when RV founder not affected but 2 or more
    # of his or her children are affected.
    return(gped)
  } else {
    #find all pairwaise combinations of the disease-affected relatives
    pair_mat <- combn(x = gped$ID[gped$affected & gped$available], m = 2)

    #find mrca for each pair of disease affected relatives, reduce to unique mrcas
    mrca <- unique(unlist(lapply(1:ncol(pair_mat), function(x){
      find_mrca(gped, pair_mat[1, x], pair_mat[2, x])
    })))

    #find the amount by which we must shift the generation number
    #if mrca contains multiple mrcas, we take the one with the smallest
    #generation number.
    if (length(mrca) > 1){
      new_gen_diff <- min(old_gen[which(gped$ID %in% mrca)]) - 1
    } else {
      new_gen_diff <- old_gen[which(gped$ID == mrca)] - 1
    }

    #assign new generation
    gped$Gen[!is.na(gped$Gen)] <- gped$Gen[!is.na(gped$Gen)] - new_gen_diff

    return(gped)
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
#' @return The censored pedigree.
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
