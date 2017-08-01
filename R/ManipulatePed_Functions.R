#' Assign generation number based on affected status
#'
#' The \code{assign_affectedGen} function assigns generation numbers among affected family members so that generation 1 represents the most recent generation that a putative disease variant shared identical by descent (IBD), as defined in Thompson (2000), by affected members could have been introduced into the pedigree.
#'
#' The \code{assign_affectedGen} function accepts a pedigree simulated by \code{sim_RVped} and reassigns generation numbers among affected family members in the pedigree.  Specifically, given a pedigree this function reassigns the generation numbers of affected members so that generation 1 is assigned to the most recent common ancestor of all affected members.  We note that the individual in generation 1 could themselves be affected, i.e. an individual can be considered their own ancestor.
#'
#' For example, consider a family with 2 affected members.  If the two affected members are a parent and his or her offspring, the affected parent would be assigned generation 1, and the affected child generation 2.  However, if the two affected members are a pair of siblings, each affected sibling would be assigned generation 2 since a common parent of the two affected siblings is assumed to be a carrier of a latent susceptibility variant.  Similarly, if the two affected members are a pair of cousins, each affected cousin is assigned generation 3, since a common grandparent of the two affected cousins is the most recent common ancestor from whom they could have inherited a shared variant associated with the disease.
#'
#' Users who wish to assign generation number based on affection status in pedigrees that have not been simulated with the \code{SimRVpedigree} package must ensure that the pedigree, \code{ped_file}, supplied to \code{assign_affectedGen} contains the following variables for each pedigree member:
#'
#' \enumerate{
#' \item \code{ID}: an identification number.
#' \item \code{dadID}: identification number of father.
#' \item \code{momID}: identification number of mother.
#' \item \code{sex}: sex identification; if male \code{sex = 0}, if female \code{sex = 1}.
#' \item \code{affected}: affection status, if affected by disease \code{affected  = 1}, otherwise, \code{affected = 0}.
#' \item \code{birthYr}: the individual's year of birth.
#' \item \code{onsetYr}: the individual's disease onset year, when applicable.
#' \item \code{deathYr}: the individual's death year, when applicable.
#' \item \code{Gen}: the individual's generation number relative to the eldest founder.  For the eldest founder \code{Gen = 1}, for his or her offspring \code{Gen = 2}, etc.
#' }
#'
#'
#' @param ped_file data.frame. A pedigree to reassign generation number based on affection status, see details.
#'
#' @return \code{reGen_ped} A pedigree containing only affected members, obligate carriers, and founders with generation number based on the most recent common ancestor of affected members as, described in details.
#' @export
#'
#' @importFrom kinship2 kinship
#' @references OUR MANUSCRIPT
#' @references Thompson, E. (2000). \emph{Statistical Inference from Genetic Data on Pedigrees.} NSF-CBMS Regional Conference Series in Probability and Statistics, 6, I-169. Retrieved from http://www.jstor.org.proxy.lib.sfu.ca/stable/4153187
#'
#' @examples
#' #Read in example pedigrees
#' data(EgPeds)
#'
#' library(kinship2)
#' #assign to pedigree object to show before and after behavior of
#' #the assign_affectedGen function
#' ex_pedigree <- pedigree(id = EgPeds$ID,
#'                         dadid = EgPeds$dadID,
#'                         momid = EgPeds$momID,
#'                         sex = (EgPeds$sex + 1),
#'                         affected = EgPeds$affected,
#'                         famid = EgPeds$FamID)
#'
#'
#' #create df to store pedigrees with reassigned generation number
#' RAG_peds <- EgPeds[1,]
#' RAG_peds <- RAG_peds[-1,]
#'
#' for(i in 1:4){
#'   RAG_peds <- rbind(RAG_peds,
#'                     assign_affectedGen(EgPeds[which(EgPeds$FamID == i), ]))
#' }
#'
#' RAG_pedigrees <-  pedigree(id = RAG_peds$ID,
#'                            dadid = RAG_peds$dadID,
#'                            momid = RAG_peds$momID,
#'                            sex = (RAG_peds$sex + 1),
#'                            affected = RAG_peds$affected,
#'                            famid = RAG_peds$FamID)
#'
#' # Compare pedigrees before and after reassigning
#' # generation number based on affcted status
#' par(mfrow = c(1, 2))
#' for (k in 1:4) {
#'   ID1 = paste0("ID", sep = ":",
#'                EgPeds[which(EgPeds$FamID == k), 2],
#'                sep = "\n Gen:", EgPeds[which(EgPeds$FamID == k), 14])
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
#'         side = 3)}
#'
#'
assign_affectedGen = function(ped_file){

  check_ped(ped_file)

  #create new ped file with affecteds only
  reGen_ped <- ped_file[which(ped_file$affected == 1), ]

  if (nrow(reGen_ped) == 0) {
    warning("No affecteds to assign affected generation")
    return(reGen_ped)
  } else {
    d <- 0
    while (d == 0) {
      #find the dad IDs that are required but have been removed
      miss_dad  <- !is.element(reGen_ped$dadID,
                               reGen_ped$ID[which(reGen_ped$sex == 0)])
      readd_dad <- reGen_ped$dadID[miss_dad]
      readd_dad <- unique(readd_dad[!is.na(readd_dad)])

      #find the mom IDs that are required but have been removed
      miss_mom  <- !is.element(reGen_ped$momID,
                               reGen_ped$ID[which(reGen_ped$sex == 1)])
      readd_mom <- reGen_ped$momID[miss_mom]
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


#' Censor pedigree after a specified year
#'
#' \code{censor_ped} censors a pedigree of any information that occurs after a specified year.
#'
#' Upon supplying a pedigree and a censor year the \code{censor_ped} function will remove all individuals born after \code{censor_year} and censor all disease onset and death events after the \code{censor_year}.
#'
#' Users who wish to censor pedigrees which have not been simulated by \code{\link{sim_RVped}} or \code{\link{sim_ped}} must ensure that the pedigree, \code{ped_file}, supplied to \code{censor_ped} contains the following variables for each pedigree member:
#'
#' \enumerate{
#' \item \code{ID}: an identification number.
#' \item \code{dadID}: identification number of father.
#' \item \code{momID}: identification number of mother.
#' \item \code{sex}: sex identification; if male \code{sex = 0}, if female \code{sex = 1}.
#' \item \code{affected}: affection status, if affected by disease \code{affected  = 1}, otherwise, \code{affected = 0}.
#' \item \code{birthYr}: the individual's birth year.
#' \item \code{onsetYr}: the individual's disease onset year, when applicable.
#' \item \code{deathYr}: the individual's death year, when applicable.
#' \item \code{proband}: (Optional) Proband identification variable, \code{proband = 1} if the individual is the proband, and 0 otherwise.
#' }
#'
#' If an individual has not experienced disease onset and/or death, then \code{onsetYr = NA} and/or \code{deathYr = NA}.
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
#' haz_obj <- new.hazard(partition = seq(0, 100, by = 1),
#'                       hazardDF = AgeSpecific_Hazards)
#'
#' #Simulate a pedigree ascertained for multiple affecteds
#' set.seed(3)
#' ex_RVped <- sim_RVped(hazard_rates = haz_obj,
#'                       num_affected = 2,
#'                       ascertain_span = c(1900, 2015),
#'                       GRR = 30, carrier_prob = 0.002,
#'                       RVfounder = "first",
#'                       stop_year = 2015,
#'                       recall_probs = c(1),
#'                       founder_byears = c(1900, 1905),
#'                       FamID = 1)[[2]]
#'
#' #To plot pedigrees, use the kinship2 package
#' library(kinship2)
#'
#' #assign simulated pedigree to pedigree object, then pass to plot function
#' Original_ped <- pedigree(id = ex_RVped$ID,
#'                          dadid = ex_RVped$dadID,
#'                          momid = ex_RVped$momID,
#'                          sex = ex_RVped$sex + 1,
#'                          affected = ex_RVped$affected)
#' plot(Original_ped)
#'
#' #Now censor ex_RVped after 1999
#' Cped <- censor_ped(ped_file = ex_RVped, censor_year = 1999)
#'
#' #Assign the the censored pedigree to a pedigree object and plot.
#' Censor_ped <- pedigree(id = Cped$ID,
#'                        dadid = Cped$dadID,
#'                        momid = Cped$momID,
#'                        sex = Cped$sex + 1,
#'                        affected = Cped$affected)
#' plot(Censor_ped)
#'
#'
censor_ped = function(ped_file, censor_year){

  check_ped(ped_file)

  if (any(is.na(match(c("birthYr", "onsetYr", "deathYr"), colnames(ped_file))))) {
    stop("ped_file does not contain one or more of the following variables: birthYr, onsetYr, deathYr" )
  }

  if (missing(censor_year)) {
    if ("proband" %in% colnames(ped_file)) {
      if(sum(ped_file$proband) == 1){
        censor_year <- ped_file$onsetYr[which(ped_file$proband == 1)]
      } else {
        stop("the proband variable cannot identify more than proband")
      }
    } else {
      stop("ped_file must contain a proband variable or the user must supply censor_year")
    }
  }

  #create new ped file containing only individuals born before the censor year
  censored_ped <- ped_file[which(ped_file$birthYr <= censor_year), ]

  if (nrow(censored_ped) == 0) {
    warning("Please check censor_year, no pedigree information prior to censor_year")
    return(censored_ped)
  } else {
    #censor onset and death events prior to censor year
    censored_ped$affected <- ifelse(is.na(censored_ped$onsetYr), 0,
                                    ifelse(censored_ped$onsetYr <= censor_year,
                                           censored_ped$affected, 0))
    censored_ped$onsetYr <- ifelse(censored_ped$onsetYr <= censor_year,
                                      censored_ped$onsetYr, NA)
    censored_ped$deathYr <- ifelse(censored_ped$deathYr <= censor_year,
                                      censored_ped$deathYr, NA)

    d <- 0
    while (d == 0) {
      #find the dad IDs that are required but have been removed
      miss_dad  <- !is.element(censored_ped$dadID,
                               censored_ped$ID[which(censored_ped$sex == 0)])
      readd_dad <- censored_ped$dadID[miss_dad]
      readd_dad <- unique(readd_dad[!is.na(readd_dad)])

      #find the mom IDs that are required but have been removed
      miss_mom  <- !is.element(censored_ped$momID,
                               censored_ped$ID[which(censored_ped$sex == 1)])
      readd_mom <- censored_ped$momID[miss_mom]
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


#' Obtain specialized pedigree information
#'
#' Obtain specialized pedigree information
#'
#' Users who wish to use \code{ped_info} for pedigrees not simulated with the \code{SImRVPedigree} package must ensure that the pedigree, \code{ped_file}, supplied to \code{ped_info} contains the following variables for each pedigree member:
#' \enumerate{
#' \item \code{ID}: an identification number.
#' \item \code{dadID}: identification number of father.
#' \item \code{momID}: identification number of mother.
#' \item \code{sex}: sex identification; if male \code{sex = 0}, if female \code{sex = 1}.
#' \item \code{affected}: affection status, if affected by disease \code{affected  = 1}, otherwise, \code{affected = 0}.}
#' Optionally, \code{ped_file} may contain any of the additional variables contained in pedigrees simulated by \code{\link{sim_RVped}}.
#'
#' @inheritParams trim_ped
#' @param ref_year A numeric constant.  The reference year used to determine current age for pedigree members.
#'
#' @return  A list containing the following:
#' @return \code{link_format} data.frame. The pedigree in linkage format; i.e. containing only the fields: FamID, ID, dadID, momID, affected.
#' @return \code{affected_info} data.frame.  Relevant information for the affected realtives only.
#' @return \code{kinshipMat} The kinship matrix for the pedigree (see kinship2 package).
#' @return \code{kinshipPedigree} An object of class pedigree (see kinship2 package).  A pedigree object which can be plotted using R's plot function.  See example.
#' @return \code{pedLabs} ID labels which can be fed to kinship2 plotting function.  See example.
#'
#' @export
#' @importFrom kinship2 kinship
#' @importFrom kinship2 pedigree
#'
#' @references Terry M Therneau and Jason Sinnwell (2015). **kinship2: Pedigree Functions.** *R package version 1.6.4.* https://CRAN.R-project.org/package=kinship2
#'
#' @examples
#' data(EgPeds)
#'
#' Fam1 <- ped_info(EgPeds[EgPeds$FamID == 1, ], ref_year = 2015)
#' summary(Fam1)
#'
#' head(Fam1$link_format)
#' Fam1$affected_info
#' Fam1$kinshipMat
#'
#' library(kinship2)
#' plot(Fam1$kinshipPedigree)
#' pedigree.legend(Fam1$kinshipPedigree, location = "topleft",  radius = 0.25)
#'
#' plot(Fam1$kinshipPedigree, id = Fam1$pedLabs)
#' pedigree.legend(Fam1$kinshipPedigree, location = "topleft",  radius = 0.25)
#'
ped_info <- function(ped_file, ref_year){

  check_ped(ped_file)

  link_format <- ped_file[, match(c("FamID", "ID", "dadID",
                                    "momID", "sex", "affected"),
                                  colnames(ped_file))]

  dA_loc <- match(c("DA1", "DA2"), colnames(ped_file))

  if (length(dA_loc[!is.na(dA_loc)]) == 2) {
    ped_file$RVstatus <- ped_file$DA1 + ped_file$DA2
  }

  keep_cols <- match(c("FamID", "ID",
                       "birthYr", "onsetYr", "deathYr",
                       "RR", "proband", "RVstatus"), colnames(ped_file))

  affected_info <- ped_file[ped_file$affected == 1 & ped_file$available == 1,
                            keep_cols[!is.na(keep_cols)]]

  if (any(!is.na(match(c("affected", "proband", "RVstatus"), colnames(ped_file))))) {
    affected_vrs = cbind(Affected = ped_file$affected,
                         Proband = ped_file$proband,
                         RV_Status = ped_file$RVstatus)
  } else {
    affected_vrs = ped_file$affected
  }

  if (!is.na(match("deathYr", colnames(ped_file)))) {
    d_status = 1*!is.na(ped_file$deathYr)

  } else {
    d_status = rep(0, nrow(ped_file))
  }

  kin_ped <- with(ped_file, pedigree(id = ID,
                                     dadid = dadID,
                                     momid = momID,
                                     sex = sex + 1,
                                     affected = affected_vrs,
                                     status = d_status))

  kinMat <- kinship(kin_ped)

  #create pedigree labels that reflect age data, when appropriate
  if(!missing(ref_year) & !is.na(match(c("birthYr"), colnames(ped_file)))){

    if (!is.na(match("deathYr", colnames(ped_file)))) {
      #create age lable, if death has not ocurred
      age_lab <- ifelse(is.na(ped_file$birthYr) | !is.na(ped_file$deathYr),
                        " ", paste0("\n age: ",
                                    ref_year - ped_file$birthYr))

      # Create a death age label for individuals who have died.
      Dage_lab <- ifelse(is.na(ped_file$deathYr),
                         " ", paste0("\n death age: ",
                                     ped_file$deathYr - ped_file$birthYr))
    } else {
      warning('Death year information not provided. Creating age lables under the assumption that all pedigree members are still alive.')
      age_lab <- ifelse(is.na(ped_file$birthYr),
                        " ", paste0("\n age: ",
                                    ref_year - ped_file$birthYr))

      Dage_lab <- rep(" ", nrow(ped_file))
    }

    if (!is.na(match("onsetYr", colnames(ped_file)))) {
      # Create a death age label for individuals who have died.
      Oage_lab <- ifelse(is.na(ped_file$onsetYr),
                         " ", paste0("\n onset age: ",
                                     ped_file$onsetYr - ped_file$birthYr))
    } else {
      Oage_lab <- rep(" ", nrow(ped_file))
    }

    age_id = paste0("ID: ", sep = "", ped_file$ID,
                    age_lab, Oage_lab, Dage_lab)
  } else {
    age_id = ped_file$ID
  }

  ped_return = list(link_format = link_format,
                    affected_info = affected_info,
                    kinshipMat = kinMat,
                    kinshipPedigree = kin_ped,
                    pedLabs = age_id)
  return(ped_return)
}
