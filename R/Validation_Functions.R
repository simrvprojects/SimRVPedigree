#' Simulate survival data for a population sample
#'
#' @inheritParams sim_RVped
#' @param nlives Numeric.  The number of individuals to simulate
#' @param YOB Numeric. The year of birth for all individuals in the sample.
#' @param stop_year Numeric. The last year of the study, i.e. the last year we can observe data.
#' @param RV_status Numeric.  The rare variant status for all individuals in the study. \code{RV_status = 1} if carrier of the cRV, otherwise \code{RV_status = 0}.  By default, \code{RV_status = NULL} so that the cRV status is simulated based on \code{carrier_prob}.
#'
#' @return A data frame of event times and classifiers for the simulated individuals.
#'
#' @export
#'
#' @examples
#' library(SimRVPedigree)
#' data(SubtypeHazards)
#'
#' a = Sys.time()
#' s_dat = sim_pop(nlives = 100,
#'                 hazard_rates = hazard(SubtypeHazards,
#'                                       subtype_ID = c("HL", "NHL")),
#'                 GRR = c(10, 1),
#'                 YOB = 1900,
#'                 stop_year = 2000)
#' b = Sys.time()
#' difftime(b, a, units = "secs")
#'
#' head(s_dat)
#' table(s_dat$RV_status)
#' table(s_dat$affected)
#' summary(s_dat$onset_age)
#' table(s_dat$subtype)
#' table(s_dat$RV_status, s_dat$subtype)
#'
#'
sim_pop <- function(nlives, hazard_rates, GRR,
                      YOB, stop_year,
                      RV_status = NULL,
                      carrier_prob = 0.002){

  if (is.null(RV_status)) {
    rvs <- sample(size = nlives, replace = TRUE,
                  x = c(1, 0), prob = c(carrier_prob, 1 - carrier_prob))
  } else {
    rvs <- rep(RV_status, nlives)
  }


  study_data <- do.call(rbind, lapply(1:nlives, function(x){
    get_ind_data(hazard_rates, GRR,
                 YOB, stop_year,
                 RVstat = rvs[x],
                 carrier_prob)
  }))
  # colnames(study_data) = c("RV_status", "affected",
  #                          "subtype", "onset_age",
  #                          "death_age", "censor_age",
  #                          "nchild", "first_birth")
  #
  # study_data[, c(1, 2, 4:8)] <- as.numeric(study_data[, c(1, 2, 4:8)])

  return(study_data)
}


#' Get individual data for sim_pop function
#'
#' This is an internal function
#'
#' @inheritParams sim_pop
#' @param RVstat The individual's cRV status
#'
#' @return A data frame of survival data
#' @keywords internal
get_ind_data <- function(hazard_rates, GRR,
                         YOB, stop_year,
                         RVstat,
                         carrier_prob){

  life_dat <- sim_life(hazard_rates, GRR, carrier_prob,
                       RV_status = RVstat, YOB, stop_year)

  # c(RV_status,
  #   !is.na(life_dat$onset_event),
  #   life_dat$subtype,
  #   ifelse(is.na(life_dat$onset_event), NA,
  #          life_dat$onset_event - as.numeric(life_dat$life_events[1])),
  #   ifelse(is.na(life_dat$death_event), NA,
  #          life_dat$death_event - as.numeric(life_dat$life_events[1])),
  #   ifelse(is.na(life_dat$death_event),
  #          life_dat$censor_year - as.numeric(life_dat$life_events[1]), NA),
  #   length(life_dat$repro_events),
  #   ifelse(length(life_dat$repro_events) > 0,
  #          life_dat$repro_events[1] - as.numeric(life_dat$life_events[1]),
  #          NA))

#
#   data.frame(RVstatus = RV_status,
#              affected = !is.na(life_dat$onset_event),
#              subtype = life_dat$subtype,
#              onset_age = ifelse(is.na(life_dat$onset_event), NA, life_dat$onset_event - as.numeric(life_dat$life_events[1])),
#              death_age = ifelse(is.na(life_dat$death_event), NA,
#                                 life_dat$death_event - as.numeric(life_dat$life_events[1])),
#              censor_age = ifelse(is.na(life_dat$death_event),
#                                  life_dat$censor_year - as.numeric(life_dat$life_events[1]), NA),
#              nchild = length(life_dat$repro_events),
#              first_birth = ifelse(length(life_dat$repro_events) > 0,
#                                   life_dat$repro_events[1] - as.numeric(life_dat$life_events[1]),
#                                   NA))

  return(data.frame(RV_status = RVstat,
                    affected = !is.na(life_dat$onset_event),
                    subtype = life_dat$subtype,
                    onset_age = ifelse(is.na(life_dat$onset_event), NA, life_dat$onset_event - as.numeric(life_dat$life_events[1])),
                    death_age = ifelse(is.na(life_dat$death_event), NA,
                                       life_dat$death_event - as.numeric(life_dat$life_events[1])),
                    censor_age = ifelse(is.na(life_dat$death_event),
                                        life_dat$censor_year - as.numeric(life_dat$life_events[1]), NA),
                    nchild = length(life_dat$repro_events),
                    first_birth = ifelse(length(life_dat$repro_events) > 0,
                                         life_dat$repro_events[1] - as.numeric(life_dat$life_events[1]),
                                         NA)))

}
