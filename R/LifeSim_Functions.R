#' Simulate next life event.
#'
#' \code{get_nextEvent} randomly simulates the next life event for an individual
#'  given their current age, disease status, and rare variant status.
#'
#' \code{get_nextEvent} randomly simulates the next life event for an individual by
#' generating the waiting times, via \link{get_WaitTime}, to reproduction, onset,
#' and death given the individuals current age.  The event with the shortest
#' waiting time is chosen as the next life event.  If get_nextEvent returns a value
#' named "Birth", then the next life event is reproduction, if get_nextEvent
#' returns a value named "Onset" then the next life event is onset of disease,
#' if get_nextEvent returns a value named "Death" then the next life event is death.
#'
#' @param current_age The individuals current age.
#' @param disease_status The disease status, disease_status = 1 if individual
#' has experienced onset, and 0 otherwise.
#' @param lambda_birth Numeric constant. The individuals birth rate.
#' @param onset_hazard Numeric vector. The population age-specific onset hazard.
#' @param death_hazard data.frame. Column 1 should specify the age specific
#' mortality rates in the unaffected population, while column 2 should provide
#' the age specific morality rates in the affected population.
#' @param part A numeric vector.  Partition of ages over which to apply the
#' age-specific hazards
#' @param birth_range A numeric vector of length 2.  The minimum and maximum
#' allowable birth ages in simulation.
#' @param RR A numeric constant. The relative risk of developing
#' disease for individuals who have inherited the rare variant.
#'
#' @return Named matrix. The number of years until the next life event,
#' named by event type.
#' @export
#' @importFrom stats runif
#' @importFrom stats rexp
#'
#' @examples
#' part_vec <- seq(0, 100, by = 1)
#' Birth_rate <- 0.05
#' unaffected_mort <- 0.00001 + pgamma(seq(0.16, 16, by = .16),
#'                                     shape = 9.5, scale = 1)
#' affected_mort <- c(0.55, 0.48, 0.37, 0.23, 0.15,
#'                    pgamma(seq(0.96, 16, by = .16), shape = 4, scale = 1.5))
#' Dhaz_df  <- as.data.frame(cbind(unaffected_mort, affected_mort))
#' Ohaz_vec <- dgamma(seq(0.1, 10, by = .1), shape = 8, scale = 0.75)
#' set.seed(17)
#' get_nextEvent(current_age = 23, disease_status = 0,
#'            lambda_birth = Birth_rate, onset_hazard = Ohaz_vec,
#'            death_hazard = Dhaz_df, part = part_vec,
#'            birth_range = c(17,45), RR = 5)
#'
#' set.seed(17)
#' get_nextEvent(current_age = 23, disease_status = 1,
#'            lambda_birth = Birth_rate, onset_hazard = Ohaz_vec,
#'            death_hazard = Dhaz_df, part = part_vec,
#'            birth_range = c(17,45), RR = 5)
#'
get_nextEvent = function(current_age, disease_status,
                         lambda_birth, onset_hazard, death_hazard, part,
                         birth_range, RR){

  # Assuming that the person is not yet affected, simulate the waiting time
  # until onset given current age
  t_onset <- ifelse(disease_status == 0,
                    get_WaitTime(p = runif(1), last_event = current_age,
                                 hazard = onset_hazard*RR,
                                 part), NA)

  #simulate the waiting time until death given current age.
  # NOTE: choosing rate = FASLE implies that we are assuming person will die
  t_death <- get_WaitTime(p = runif(1),
                          last_event = current_age,
                          hazard = death_hazard[, (1 + disease_status)],
                          part, scale = TRUE)

  # Want to adjust the waiting time until birth based on current age
  # and also ensure that birth cannot occur after the maximum birth age

  # First condition handles the case when person is less than the minimum
  # birth age when this occurs we return current age + waiting time
  # so long as this value is less than the maximum birth age
  # Second condition deals with people who have reached the minimum birth age,
  # if neiter is met we return NA
  nyear_birth <- rexp(1, lambda_birth)
  t_birth <- ifelse((current_age < birth_range[1] &
                       nyear_birth + birth_range[1] <= birth_range[2]),
                    nyear_birth + birth_range[1],
                    ifelse((current_age >= birth_range[1] &
                              (nyear_birth + current_age) <= birth_range[2]),
                           nyear_birth,
                           NA))

  #create list of simulated times
  times <- as.list(c(Birth = t_birth,  Onset = t_onset, Death = t_death))

  #find the smallest waiting time
  if (sum(which.min(times)) != 0) {
    nyears <- as.matrix(times[[which.min(times)]], ncol = 1)
    colnames(nyears) <- c(paste(names(which.min(times))))
  }else if (sum(which.min(times)) == 0) {
    nyears <- as.matrix(0, ncol = 1)
    colnames(nyears) <- c("retry")
  }

  return(nyears)
}

#' Simulate all life events.
#'
#' \code{get_lifeEvents} simulates all life events for an individual starting at age
#' 0 and ending with death by applying the \link{get_nextEvent} function until death
#' occurs.
#'
#' @param NB_params A numeric vector of length 2. The size and probabiliy
#' parameters of the negative binomial distribution that describes the number of
#' children per household in the population.
#' @param YOB numeric. The indivdiual's year of birth.
#' @inheritParams get_nextEvent
#'
#' @return Named matrix. The waiting times between all life events simulated for an individual, named according to which life event has occurred.
#' @export
#' @importFrom stats rgamma
#'
#' @examples
#' part_vec <- seq(0, 100, by = 1)
#' unaffected_mort <- 0.00001 + pgamma(seq(0.16, 16, by = .16),
#'                                     shape = 9.5, scale = 1)
#' affected_mort <- c(0.55, 0.48, 0.37, 0.23, 0.15,
#'                    pgamma(seq(0.96, 16, by = .16), shape = 4, scale = 1.5))
#' Dhaz_df  <- (as.data.frame(cbind(unaffected_mort, affected_mort)))
#' Ohaz_vec <- (dgamma(seq(0.1, 10, by = .1), shape = 8, scale = 0.75))
#' set.seed(17)
#' get_lifeEvents(onset_hazard = Ohaz_vec,
#'                death_hazard = Dhaz_df,
#'                part = part_vec,
#'                birth_range = c(17,45),
#'                NB_params = c(2, 4/7), RR = 1,
#'                YOB = 1900)
#'
#' set.seed(17)
#' get_lifeEvents(onset_hazard = Ohaz_vec,
#'                death_hazard = Dhaz_df,
#'                part = part_vec,
#'                birth_range = c(17,45),
#'                NB_params = c(2, 4/7), RR = 15,
#'                YOB = 1900)
#'
get_lifeEvents = function(onset_hazard, death_hazard, part,
                          birth_range, NB_params, RR, YOB){

  #initialize data frame to hold life events
  R.life  <- data.frame(Start = 0)
  min.age <- min(part)
  max.age <- max(part)
  #initialize disease status, start age at minumum age permissable under part
  DS <- 0; t <- min.age

  #generate and store the birth rate for this individual
  B.lambda <- rgamma(1, shape = NB_params[1],
                       scale = (1-NB_params[2])/NB_params[2])/(birth_range[2] -
                                                                 birth_range[1])
  while(t < max.age){
    #generate next event
    my.step <- get_nextEvent(current_age = t, disease_status = DS,
                             lambda_birth = B.lambda, onset_hazard, death_hazard,
                             part, birth_range, RR)

    #add to previous life events
    R.life <- cbind(R.life, my.step)

    if (colnames(my.step) == "Death") {
      #if death occurs stop simulation by setting t = 100
      t <- max.age
    } else if (colnames(my.step) == "Onset") {
      #if onset occurs change disease status and continue
      t  <- t + as.numeric(my.step[1,1])
      DS <- 1
    } else {
      #otherwise, increment counter, handles both birth event and case when no
      # event occurs (i.e. retry)
      t <- t + as.numeric(my.step[1,1])
    }
  }

  life_events <- round(cumsum(as.numeric(R.life[1,]))) + YOB
  names(life_events) <- names(R.life)
  return(life_events)

}
