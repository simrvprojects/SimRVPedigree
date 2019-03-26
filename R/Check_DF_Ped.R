#' Checks ped_file for expected information, used before converting to ped object.
#'
#' @param ped_file data.frame The pedigree in data frame format.
#'
#' @keywords internal
#'
check_ped <- function(ped_file){
  if (class(ped_file) != "data.frame") {
    stop("ped_file must be a data.frame with the following variables:\n FamID, ID, dadID, momID, sex, affected")
  }


  if (!"FamID" %in% colnames(ped_file) |
      !"ID" %in% colnames(ped_file) |
      !"dadID" %in% colnames(ped_file) |
      !"momID" %in% colnames(ped_file) |
      !"sex" %in% colnames(ped_file) |
      !"affected" %in% colnames(ped_file)) {
    stop('please provide a data.frame with the following variables:\n FamID, ID, dadID, momID, sex, affected')
  }

  if(any(is.na(ped_file$ID))) {
    stop('ID contains missing values.\n  Please ensure all individuals have a valid ID.')
  }

  if (any(!ped_file$affected %in% c(TRUE, FALSE, NA))) {
    stop('For the variable "affected" please use the following convention
         TRUE = affected by disease
         FALSE = unaffected
         NA = unknown disease-affection status.\n')
  }

  moms <- unique(ped_file$momID[!is.na(ped_file$momID)])
  dads <- unique(ped_file$dadID[!is.na(ped_file$dadID)])

  if (any(ped_file$sex[which(ped_file$ID %in% moms)] != 1) |
      any(ped_file$sex[which(ped_file$ID %in% dads)] != 0)){

    wrong_sex <- c(ped_file$ID[which(ped_file$sex[which(ped_file$ID %in% dads)] != 0)],
                   ped_file$ID[which(ped_file$sex[which(ped_file$ID %in% moms)] != 1)])

    stop(paste0('Sex improperly specifed ID: ', sep = '', wrong_sex, '.  Please ensure that for males: sex = 0; and for females: sex = 1.'))
  }

  if (any(!moms %in% ped_file$ID) | any(!dads %in% ped_file$ID)) {

    wrong_par <- c(ped_file$ID[which(ped_file$momID == moms[which(!moms %in% ped_file$ID)])],
                   ped_file$ID[which(ped_file$dadID == dads[which(!dads %in% ped_file$ID)])])

    stop(paste0('ID: ', sep = '', wrong_par, '.  Non-founders must have a mother and a father. Founders have neither.'))
  }

  if (any(!is.na(ped_file$momID[is.na(ped_file$dadID)])) |
      any(!is.na(ped_file$dadID[is.na(ped_file$momID)]))) {
    stop("Non-founders must have both a mother and a father, while founders have missing momID and dadID.")
  }
}
