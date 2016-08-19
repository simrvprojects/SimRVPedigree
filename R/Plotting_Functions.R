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
#' @importFrom kinship2 pedigree.legend
#' @importFrom graphics plot
#' @importFrom graphics legend
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
plot_RVpedigree = function(ped_file, legend_location, A_colors){

  RV_status <- ped_file$DA1 + ped_file$DA2
  Affected  <- ped_file$affected
  Proband   <- ped_file$is_proband
  Available <- ped_file$available

  if (missing(A_colors)){
    A_colors <- c(rgb(red = 0/225, green = 112/225, blue = 150/225),  #SFUblue,
                 rgb(red = 166/225, green = 25/225, blue = 46/225))  #SFUred)
  }

  Ped_cols <- ifelse((Available == 1), A_colors[1], A_colors[2])

  RV_ped <- pedigree(ped_file$ID, ped_file$dad_id,
                    ped_file$mom_id, (ped_file$gender + 1),
                    affected = cbind(Affected,
                                     Proband,
                                     RV_status),
                    famid = ped_file$FamID)
  rvPed = RV_ped[paste0(ped_file$FamID[1])]
  plot(rvPed, col = Ped_cols)

  if (missing(legend_location)){
    pedigree.legend(rvPed, location="topleft", radius=.35)
  } else {
    pedigree.legend(rvPed, location= legend_location, radius=.35)
  }
  legend("topright", title = "Availability Status",
         legend = c("available", "unavailable"),
         col = A_colors,
         lwd = 4)

}
