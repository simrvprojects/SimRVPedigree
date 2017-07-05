#' Simulate the next life event
#'
#' Primarily intended as an internal function, \code{get_nextEvent} randomly simulates an individual's next life event given their current age, disease status, and relative-risk of disease.
#'
#' Given their current age, \code{get_nextEvent} randomly simulates an individual's next life event by generating waiting times to reproduction, onset, and death.  The event with the shortest waiting time is chosen as the next life event.
#'
#'  We assume that, given an individual's current age, their time to disease onset is the waiting time in a non-homogeneous Poisson process with an age-specific hazard rate that follows a proportional hazards model.  In this model, individuals who have NOT inherited the rare variant experience disease onset according to the baseline (or population) hazard rate of disease.  On the other hand, individuals who have inherited the rare variant are assumed to have an increased risk of disease onset relative to those who have inherited it.  The user is expected to supply the baseline hazard rate of disease, as well as the relative-risk of disease for genetic cases. Additionally, we impose the restriction that individuals may only experience disease onset once, and remain affected from that point on.
#'
#' We assume that, given an individual's current age, their time to death is the waiting time in a non-homogeneous Poisson process with age-specific hazard rate determined by their affection status.  We assume that unaffected individuals experience death according to the age-specific hazard rate for death in the unaffected population.  If the disease of interest is sufficiently rare, the user may instead choose to substitute the population age-specific hazard rate for death in the general population.  We assume that affected individuals experience death according to the age-specific hazard rate for death in the affected population.  The user is expected to supply both of these age-specific hazard rates.
#'
#' We assume that, given an individual's current age, their time to reproduction is the waiting time in a homogeneous Poisson process.  That is, we assume that individuals reproduce at uniform rate during their reproductive years.  For example, one's reproductive years may span from age 18 to age 45.  We do not allow for offspring to be produced outside of an individual's \code{birth_range}.
#'
#'
#' If get_nextEvent returns the waiting time to the next life event, named for event type.  The possible event types are as follows:
#' \itemize{
#'  \item "Child" a reproductive event, i.e. creation of offspring
#'  \item "Onset" disease onset event,
#'  \item "Death" death event
#' }
#'
#' @param current_age Numeric. The individual's current age.
#' @param disease_status Numeric. The individual's disease status, where \code{disease_status = 1} if individual has experienced disease onset, otherwise \code{disease_status = 0}.
#' @param lambda_birth Numeric. The individual's birth rate.
#' @param RR Numeric. The individual's relative-risk of disease.
#' @inheritParams sim_RVped
#'
#' @return A named matrix. The number of years until the next life event,
#' named by event type.  See Details.
#'
#' @keywords internal
#' @importFrom stats runif
#' @importFrom stats rexp
#'
#'
get_nextEvent = function(current_age, disease_status,
                         lambda_birth, hazard_rates,
                         birth_range, RR){

  # Assuming that the person is not yet affected, simulate the waiting time
  # until onset given current age
  t_onset <- ifelse(disease_status == 0,
                    get_WaitTime(p = runif(1), last_event = current_age,
                                 hazard = hazard_rates[[1]][, 1]*RR,
                                 part = hazard_rates[[2]]), NA)

  #simulate the waiting time until death given current age.
  # NOTE: choosing rate = FASLE implies that we are assuming person will die
  t_death <- get_WaitTime(p = runif(1),
                          last_event = current_age,
                          hazard = hazard_rates[[1]][, (2 + disease_status)],
                          part = hazard_rates[[2]], scale = TRUE)

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
  times <- as.list(c(Child = t_birth,  Onset = t_onset, Death = t_death))

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

#' Simulate all life events
#'
#' Primarily intended as an internal function, \code{sim_lifeEvents} simulates all life events for an individual starting at birth, age 0, and ending with death or the end of the study.
#'
#' Starting at birth, age 0, the \code{sim_lifeEvents} generates waiting times to reproduction, onset, and death. The event with the shortest waiting time is chosen as the next life event, and individual's age is updated by the waiting time of the winning event.  Conditioned on the individual's new age, this process is applied recursively, until death or until the end of the study is reached.
#'
#'  We make the following assumptions regarding the simulation of waiting times:
#'  \enumerate{
#'  \item We assume that, given an individual's current age, their time to disease onset is the waiting time in a non-homogeneous Poisson process with an age-specific hazard rate that follows a proportional hazards model.  In this model, individuals who have NOT inherited the rare variant experience disease onset according to the baseline (or population) hazard rate of disease.  On the other hand, individuals who have inherited the rare variant are assumed to have an increased risk of disease onset relative to those who have inherited it.  The user is expected to supply the baseline hazard rate of disease, as well as the relative-risk of disease for genetic cases. Additionally, we impose the restriction that individuals may only experience disease onset once, and remain affected from that point on.
#'
#'  \item We assume that, given an individual's current age, their time to death is the waiting time in a non-homogeneous Poisson process with age-specific hazard rate determined by their affection status.  We assume that unaffected individuals experience death according to the age-specific hazard rate for death in the unaffected population; if the disease of interest is sufficiently rare, the user may instead choose to substitute the population age-specific hazard rate for death in the general population.  We assume that affected individuals experience death according to the age-specific hazard rate for death in the affected population.  The user is expected to supply both of these age-specific hazard rates.
#'
#'  \item We assume that, given an individual's current age, their time to reproduction is the waiting time in a homogeneous Poisson process.  That is, we assume that individuals reproduce at uniform rate during their reproductive years.  For example, one's reproductive years may span from age 18 to age 45 years.  We do not allow for offspring to be produced outside of an individual's \code{birth_range}.
#'
#'  }
#'
#'
#' \code{sim_lifeEvents} will return a named matrix, which contains the years of the simulated life events, named by event type.  The possible event types are as follows:
#' \itemize{
#'  \item "Child" a reproductive event, i.e. creation of offspring
#'  \item "Onset" disease onset event,
#'  \item "Death" death event
#' }
#'
#' @param YOB A positive number. The indivdiual's year of birth.
#' @param RR Numeric. The individual's relative-risk of disease.
#' @inheritParams sim_RVped
#'
#' @return A named matrix containing the years of an individual's simulated life events, named by event type, see details.
#'
#' @references OUR MANUSCRIPT
#' @references Ken-Ichi Kojima, Therese M. Kelleher. (1962), \emph{Survival of Mutant Genes}. The American Naturalist 96, 329-346.
#' @export
#' @importFrom stats rgamma
#'
#' @examples
#' data(AgeSpecific_Hazards)
#'
#' my_HR <- new.hazard(hazardDF = AgeSpecific_Hazards)
#'
#' # The following commands simulate all life events for an individual, whose
#' # relative-risk of disease is 1, born in 1900.  From the output, this
#' # individual has 1 child in 1944, and then dies in 1961.
#' set.seed(1234)
#' sim_lifeEvents(hazard_rates = my_HR,
#'                birth_range = c(17,45),
#'                NB_params = c(2, 4/7), RR = 1,
#'                YOB = 1900, stop_year = 2000)
#'
#' # Using the same random seed, notice how the life events can vary for
#' # someone with an increased relative-risk of disease, say 25.
#' # From the output, this individual also has 1 child in 1944, then
#' # experiences disease onset in 1949, and dies in 1957.
#' set.seed(1234)
#' sim_lifeEvents(hazard_rates = my_HR,
#'                birth_range = c(17,45),
#'                NB_params = c(2, 4/7), RR = 25,
#'                YOB = 1900, stop_year = 2000)
#'
sim_lifeEvents = function(hazard_rates, birth_range, NB_params,
                          RR, YOB, stop_year){

  if(!is.hazard(hazard_rates)) {
    stop("hazard_rates must be an object of class hazard")
  }

  #initialize data frame to hold life events
  R_life  <- data.frame(Start = 0)
  min_age <- min(hazard_rates[[2]])
  max_age <- max(hazard_rates[[2]])
  #initialize disease status, start age at minumum age permissable under part
  DS <- 0; t <- min_age; yr <- YOB

  #generate and store the birth rate for this individual
  B_lambda <- rgamma(1, shape = NB_params[1],
                       scale = (1-NB_params[2])/NB_params[2])/(birth_range[2] -
                                                                 birth_range[1])
  while(t < max_age & yr <= stop_year){
    #generate next event
    l_event <- get_nextEvent(current_age = t, disease_status = DS,
                             lambda_birth = B_lambda, hazard_rates,
                             birth_range, RR)

    if(yr + as.numeric(l_event[1,1]) <= stop_year){
      #add to previous life events
      R_life <- cbind(R_life, l_event)

      if (colnames(l_event) == "Death") {
        #if death occurs stop simulation by setting t = 100
        t <- max_age
      } else if (colnames(l_event) == "Onset") {
        #if onset occurs change disease status and continue
        t  <- t + as.numeric(l_event[1,1])
        yr <- yr + as.numeric(l_event[1,1])
        DS <- 1
      } else {
        #otherwise, increment counter, handles both birth event and case when no
        # event occurs (i.e. retry)
        t <- t + as.numeric(l_event[1,1])
        yr <- yr + as.numeric(l_event[1,1])
      }
    } else {
      yr <- yr + as.numeric(l_event[1,1])
    }

  }

  life_events <- round(cumsum(as.numeric(R_life[1,]))) + YOB
  names(life_events) <- names(R_life)
  return(life_events)

}
