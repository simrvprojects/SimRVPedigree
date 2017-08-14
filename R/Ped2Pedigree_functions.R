#' Convert a ped object to a kinship2 pedigree object
#'
#' Convert a SimRVPedigree ped object to a kinship 2 pedigree object
#'
#' @param x A ped object.
#'
#' @return A pedigree object.  See ?kinship2::pedigree for details.
#' @importFrom kinship2 pedigree
#' @export
#'
#' @examples
#' data(EgPeds)
#'
#' ped1 = new.ped(EgPeds[EgPeds$FamID == 1, ])
#' head(ped1, n = 3)
#' class(ped1)
#' summary(ped1)
#'
#' kinPed <- ped2pedigree(ped1)
#' kinPed
#'
#' library(kinship2)
#' plot(kinPed)
#' pedigree.legend(kinPed, location="topleft", radius=.25)
#'
#'
#'
#' AllPeds = new.ped(EgPeds)
#' summary(AllPeds)
#'
#' kinPed_multi <- ped2pedigree(AllPeds)
#' kinPed_multi
#'
#' unique(AllPeds$FamID)
#'
#' # If pedigree object contains more than 1 family
#' # take care to specify the family ID when plotting
#' plot(kinPed_multi['2'])
#' pedigree.legend(kinPed_multi['2'], location="topleft", radius=.25)
#'
#' plot(kinPed_multi['4'])
#' pedigree.legend(kinPed_multi['4'], location="topleft", radius=.25)
#'
#'
#'
ped2pedigree <- function(x){

  if (!is.ped(x)) stop("please supply an object of class ped")

  m <- length(unique(x$FamID))

  dA_loc <- match(c("DA1", "DA2"), colnames(x))

  if (length(dA_loc[!is.na(dA_loc)]) == 2) {
    x$RVstatus <- x$DA1 + x$DA2
  }

  if (any(!is.na(match(c("affected", "proband", "RVstatus"), colnames(x))))) {
    affected_vrs = cbind(Affected = x$affected,
                         Proband = x$proband,
                         RV_Status = x$RVstatus)
  } else {
    affected_vrs = x$affected
  }

  if (!is.na(match("deathYr", colnames(x)))) {
    d_status <- 1*!is.na(x$deathYr)
  } else {
    d_status <- rep(0, nrow(x))
  }

  if (m == 1){
    kin_ped <- with(x, pedigree(id = ID,
                                dadid = dadID,
                                momid = momID,
                                sex = sex + 1,
                                affected = affected_vrs,
                                status = d_status))
  } else {
    kin_ped <- with(x, pedigree(famid = FamID,
                                id = ID,
                                dadid = dadID,
                                momid = momID,
                                sex = sex + 1,
                                affected = affected_vrs,
                                status = d_status))
  }

  return(kin_ped)
}



#' Create pedigree labels
#'
#' @param x An object of class ped.
#' @param ref_year A numeric constant.  The reference year used to determine current age for pedigree members.
#'
#' @return A list of labels that can be used with kinship2's plot function.
#' @export
#'
#' @examples
#' data(EgPeds)
#'
#' ped1 = new.ped(EgPeds[EgPeds$FamID == 2, ])
#' head(ped1, n = 3)
#' class(ped1)
#' summary(ped1)
#'
#' kinPed <- ped2pedigree(ped1)
#' kinPed
#' kinLabs <- pedigreeLabels(ped1, ref_year = 2016)
#'
#' library(kinship2)
#' plot(kinPed,
#'      id = kinLabs,
#'      mar = c(3, 3, 4, 3) + 0.1)
#'
#'
#' AllPeds = new.ped(EgPeds)
#' summary(AllPeds)
#'
#' kinPed_multi <- ped2pedigree(AllPeds)
#' kinPed_multi
#' kinLabs_multi <- pedigreeLabels(AllPeds, ref_year = 2016)
#'
#' unique(AllPeds$FamID)
#'
#' # If pedigree object contains more than 1 family
#' # take care to specify the family ID when plotting
#' plot(kinPed_multi['4'],
#'      id = kinLabs_multi[AllPeds$FamID == 4],
#'      mar = c(3, 3, 4, 3) + 0.1)
#'
#'
pedigreeLabels <- function(x, ref_year){

  if (!is.ped(x)) stop("please supply an object of class ped")

  m <- length(unique(x$FamID))

  #create pedigree labels that reflect age data, when appropriate
  if(!missing(ref_year) & !is.na(match(c("birthYr"), colnames(x)))){

    if (!is.na(match("deathYr", colnames(x)))) {
      #create age lable, if death has not ocurred
      age_lab <- ifelse(is.na(x$birthYr) | !is.na(x$deathYr),
                        " ", paste0("\n age: ",
                                    ref_year - x$birthYr))

      # Create a death age label for individuals who have died.
      Dage_lab <- ifelse(is.na(x$deathYr),
                         " ", paste0("\n death age: ",
                                     x$deathYr - x$birthYr))
    } else {
      warning('Death year information not provided. Creating age lables under the assumption that all pedigree members are still alive.')
      age_lab <- ifelse(is.na(x$birthYr),
                        " ", paste0("\n age: ",
                                    ref_year - x$birthYr))

      Dage_lab <- rep(" ", nrow(x))
    }

    if (!is.na(match("onsetYr", colnames(x)))) {
      # Create a death age label for individuals who have died.
      Oage_lab <- ifelse(is.na(x$onsetYr),
                         " ", paste0("\n onset age: ",
                                     x$onsetYr - x$birthYr))
    } else {
      Oage_lab <- rep(" ", nrow(x))
    }

    ped_labs = paste0("ID: ", sep = "", x$ID,
                    age_lab, Oage_lab, Dage_lab)
  } else {
    ped_labs = x$ID
  }

  return(ped_labs)

}

