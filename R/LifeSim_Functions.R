#' Simulate a random birth rate
#'
#' @param NB_size,NB_prob size and probability parameters of negative binomial distribution.
#' @param min_birth_age,max_birth_age minimum and maximum allowable birth ages
#'
#' @return birth_rate numeric, the randomly generated birth rate
#'
#' @examples
#' set.seed(17)
#' rbirth_rate(NB_size = 2, NB_prob = 4/7, min_birth_age = 17, max_birth_age = 45)
#'
rbirth_rate = function(NB_size, NB_prob, min_birth_age, max_birth_age){
  birth_rate = rgamma(1, shape = NB_size,
                      scale = (1-NB_prob)/NB_prob ) / (max_birth_age-min_birth_age)
}

#' Simulate the next life event for an individual
#'
#' \code{event_step} Randomly simulates the next life event for an individual by
#' generating the waiting times, via \link{findWaitTime}, to reproduction, onset, and death given the
#' individuals current age and then choosing the event with the shortest waiting
#' time as the winner.
#'
#' @param current_age The individuals current age.
#' @param disease_status The disease status, disease_status = 1 if individual
#' has experienced onset, and 0 otherwise.
#' @param RV_status Rare variant status, RV_status  = 1 if
#' individual has inherited the rare variant, and 0 otherwise.
#' @param lambda_birth Numeric constant. The individuals birth rate.
#' @param onset_hazard Numeric vector. The population age-specific onset hazard.
#' @param death_hazard data.frame. Column 1 should specify the age specific
#' mortality rates in the unaffected population, while column 2 should provide
#' the age specific morality rates in the affected population.
#' @param part A numeric vector.  Partition of ages over which to apply the
#' age-specific hazards
#' @param birth_range A numeric vector of length 2.  The minimum and maximum
#' allowable birth ages in simulation.
#' @param RR A numeric constant. The rare variant relative risk of developing
#' disease.
#'
#' @return Named matrix. The number of years until the next life event,
#' named by event type.
#'
#' @examples
#' part_vec <- seq(0, 100, by = 1)
#' Brate <- rbirth_rate(NB_size = 2, NB_prob = 4/7,
#'                      min_birth_age = 17, max_birth_age = 45)
#' unaffected_mort <- 0.00001 + pgamma(seq(0.16, 16, by = .16),
#'                                     shape = 9.5, scale = 1)
#' affected_mort <- c(0.55, 0.48, 0.37, 0.23, 0.15,
#'                    pgamma(seq(0.96, 16, by = .16), shape = 4, scale = 1.5))
#' Dhaz_df  <- as.data.frame(cbind(unaffected_mort, affected_mort))
#' Ohaz_vec <- dgamma(seq(0.1, 10, by = .1), shape = 8, scale = 0.75)
#' set.seed(17)
#' event_step(current_age = 23, disease_status = 0, RV_status = 0,
#'            lambda_birth = Brate, onset_hazard = Ohaz_vec,
#'            death_hazard = Dhaz_df, part = part_vec,
#'            birth_range = c(17,45), RR = 15)
#'
#' set.seed(17)
#' event_step(current_age = 23, disease_status = 1, RV_status = 1,
#'            lambda_birth = Brate, onset_hazard = Ohaz_vec,
#'            death_hazard = Dhaz_df, part = part_vec,
#'            birth_range = c(17,45), RR = 15)
#'
event_step = function(current_age, disease_status, RV_status,
                      lambda_birth, onset_hazard, death_hazard, part,
                      birth_range, RR){

  #first pull appropriate columns from death_hazard and onset_hazard for given
  # RV status and disease status
  risk.lambda  <- ifelse(RV_status == 1, RR, 1)
  onset.lambda <- onset_hazard*risk.lambda
  death.lambda <- death_hazard[, (1 + disease_status)]

  #Assuming that the person is not yet affected, simulate the waiting time until onset given current age
  t.onset <- ifelse(disease_status == 0,
                    findWaitTime(u = runif(1), last_event = current_age,
                                 hazard = onset.lambda, part), NA)

  #simulate the waiting time until death given current age.
  # NOTE: choosing rate = FASLE implies that we are assuming person will die
  t.death <- findWaitTime(u = runif(1),
                          last_event = current_age,
                          hazard = death.lambda,
                          part,
                          scale = TRUE)

  # Want to adjust the waiting time until birth based on current age
  # and also ensure that birth cannot occur after the maximum birth age
  nyear.birth <- rexp(1, lambda_birth)
  t.birth <- ifelse(
    # this condition handles the case when person is less than the minimum
    # birth age when this occurs we return current age + waiting time
    # so long as this value is less than the maximum birth age
    ((current_age < birth_range[1]) & (nyear.birth + birth_range[1] <= birth_range[2])),
    nyear.birth + birth_range[1],
    # next we deal with people who have reached the minimum birth age,
    ifelse(((current_age >= birth_range[1]) & (nyear.birth + current_age <= birth_range[2])),
           nyear.birth,
           #if neither condition is met we return NA
           NA))

  #create list of simulated times
  times <- as.list(c(Birth = t.birth,  Onset = t.onset, Death = t.death))

  #find the smallest waiting time
  if (sum(which.min(times)) == 0) {
    nyears = as.matrix(0, ncol = 1)
    colnames(nyears) = c("retry")
  }else if (sum(which.min(times)) != 0) {
    nyears = as.matrix(times[[which.min(times)]], ncol = 1)
    colnames(nyears) = c(paste(names(which.min(times))))
  }
  return(nyears)
}

#' Simulate all life events, starting at birth and ending with death.
#'
#' @param RV_status Rare variant status, RV_status  = 1 if
#' individual has inherited the rare variant, and 0 otherwise.
#' @param onset_hazard Numeric vector. The population age-specific onset hazard.
#' @param death_hazard data.frame. Column 1 should specify the age specific
#' mortality rates in the unaffected population, while column 2 should provide
#' the age specific morality rates in the affected population.
#' @param part A numeric vector.  Partition of ages over which to apply the
#' age-specific hazards
#' @param birth_range A numeric vector of length 2.  The minimum and maximum
#' allowable birth ages in simulation.
#' @param NB_params A numeric vector of length 2. The size and probabiliy parameters of the negative binomial distribution that describes the number of children per household in the population.
#' @param RR A numeric constant. The rare variant relative risk of developing
#' disease.
#'
#' @return Named matrix. All life events randomly simulated for an individual.
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
#' life_step(RV_status = 0, onset_hazard = Ohaz_vec,
#'            death_hazard = Dhaz_df, part = part_vec,
#'            birth_range = c(17,45), NB_params = c(2, 4/7), RR = 15)
#'
#' set.seed(17)
#' life_step(RV_status = 1, onset_hazard = Ohaz_vec,
#'            death_hazard = Dhaz_df, part = part_vec,
#'            birth_range = c(17,45), NB_params = c(2, 4/7), RR = 15)
#'

life_step = function(RV_status, onset_hazard, death_hazard, part,
                     birth_range, NB_params, RR){

  #initialize data frame to hold life events
  R.life  <- data.frame(Start = 0)
  min.age <- min(part)
  max.age <- max(part)
  #initialize disease status, start age at minumum age permissable under part
  DS <- 0; t <- min.age

  #generate and store the birth rate for this individual
  B.lambda <- rbirth_rate( NB_params[1], NB_params[2] ,
                           birth_range[1], birth_range[2])

  while(t < max.age){
    #generate next event
    my.step <- event_step(current_age = t, disease_status = DS, RV_status,
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

  return(R.life)
}
