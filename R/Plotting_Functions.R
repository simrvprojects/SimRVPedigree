#' Plot simulated pedigree
#'
#' Plot ped files generated using \code{SimRVPedigree} using the pedigree plotting tools provided by the \link{kinship2} package.
#'
#' @param ped_file A single pedigree simulated with either \code{sim_RVpedigree} or \code{ped_step}
#'
#' @return Values returned from \link[kinship2]{plot.pedigree}
#' @export
#'
#' @importFrom kinship2 pedigree
#' @importFrom kinship2 plot.pedigree
#' @importFrom graphics plot
#'
#' @examples
#' part_vec <- seq(0, 100, by = 1)
#' unaffected_mort <- 0.00001 + pgamma(seq(0.16, 16, by = .16),
#'                                     shape = 9.5, scale = 1)/350
#' affected_mort <- c(0.55, 0.48, 0.37, 0.23, 0.15,
#'                    pgamma(seq(0.96, 16, by = .16), shape = 4, scale = 1.5))/300
#' Dhaz_df  <- (as.data.frame(cbind(unaffected_mort, affected_mort)))
#' Ohaz_vec <- (dgamma(seq(0.1, 10, by = .1), shape = 8, scale = 0.75))
#' set.seed(22)
#' ex_ped <- ped_step(onset_hazard = Ohaz_vec, death_hazard = Dhaz_df,
#'                    part = part_vec, RR = 5, founder_byears = c(1900, 1910))
#' ex_ped$FamID = 1
#' plot_RVpedigree(ex_ped)
plot_RVpedigree = function(ped_file){
  RV_ped = pedigree(ped_file$ID, ped_file$dad_id,
                    ped_file$mom_id, (ped_file$gender + 1),
                    affected = cbind(ped_file$affected,
                                     ped_file$is_proband,
                                     ped_file$available),
                    famid = ped_file$FamID)
  rvPed = RV_ped[paste0(ped_file$FamID[1])]
  plot(rvPed)
}
