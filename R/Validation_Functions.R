#' Simulate a sample of individuals
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
#' data(AgeSpecific_Hazards)
#'
#' a = Sys.time()
#' s_dat = sim_pop(nlives = 1000,
#'                 hazard_rates = hazard(AgeSpecific_Hazards[, c(1, 1, 2, 3)],
#'                                       subtype_ID = c("HL", "NHL")),
#'                 GRR = 10,
#'                 YOB = 1900,
#'                 stop_year = 2000)
#' b = Sys.time()
#' difftime(b, a, units = "secs")
#'
#' head(s_dat)
#' table(s_dat$RV_status)
#' table(s_dat$affected)
#' table(s_dat$onset_age)
#' table(s_dat$subtype)
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


  study_data <- data.frame(RV_status   = rvs,
                           affected    = rep(NA, nlives),
                           subtype     = rep(NA, nlives),
                           onset_age   = rep(NA, nlives),
                           death_age   = rep(NA, nlives),
                           nchild      = rep(NA, nlives),
                           first_birth = rep(NA, nlives))

  for (k in 1:nlives) {
    life_dat <- sim_life2(hazard_rates, GRR, carrier_prob,
                          RV_status = study_data$RV_status[k], YOB, stop_year)

    death_yr <- as.numeric(life_dat[[1]][1,])[which(colnames(life_dat[[1]])=="Death")]
    onset_yr <- as.numeric(life_dat[[1]][1,])[which(colnames(life_dat[[1]])=="Onset")]

    #tabulate life events so that we can count number of offspring
    Life.counts <- table(colnames(life_dat[[1]]))

    #count the number of children for each individual
    study_data$nchild[k] <- ifelse(is.element("Child", names(Life.counts)),
                                   Life.counts[which(names(Life.counts)=="Child")],
                                   0)

    study_data$first_birth[k] <- ifelse(is.element("Child", names(Life.counts)),
                                        min(as.numeric(life_dat[[1]][1,])[which(colnames(life_dat[[1]])=="Child")]) - YOB,
                                        NA)

    #set Is.Onset = 1, if person develops cancer, otherwise set to zero.
    study_data$affected[k] <- ifelse(is.element("Onset", names(Life.counts)), 1, 0)
    study_data$subtype[k] <- life_dat[[2]]
    study_data$onset_age[k] <- ifelse(is.element("Onset", names(Life.counts)),
                                      onset_yr - YOB,
                                      NA)
    study_data$death_age[k] <- ifelse(is.element("Death", names(Life.counts)),
                                      death_yr - YOB,
                                      NA)
  }

  return(study_data)
}
