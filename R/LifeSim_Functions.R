#' Simulate the next life event
#'
#' Primarily intended as an internal function, \code{get_nextEvent} randomly simulates an individual's next life event given their current age, disease status, and relative-risk of disease.
#'
#' Given their current age, \code{get_nextEvent} randomly simulates an individual's next life event by generating waiting times to reproduction, onset, and death.  The event with the shortest waiting time is chosen as the next life event.
#'
#'  We assume that, given an individual's current age, their time to disease onset is the waiting time in a non-homogeneous Poisson process with an age-specific hazard rate that follows a proportional hazards model.  In this model, individuals who have NOT inherited the rare variant experience disease onset according to the baseline (or population) hazard rate of disease.  On the other hand, individuals who have inherited the rare variant are assumed to have an increased risk of disease onset relative to those who have inherited it.  The user is expected to supply the population age-specific hazard rate of disease, the relative-risk of disease for genetic cases, and the population allele frequency of the rare variant.  We impose the restriction that individuals may only experience disease onset once, and remain affected from that point on.
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
#' @param birth_range Numeric vector of length 2. The minimum and maximum allowable ages, in years, between which this individual may reproduce.
#' @inheritParams sim_RVped
#' @inheritParams sim_life
#'
#' @return A named matrix. The number of years until the next life event,
#' named by event type.  See Details.
#'
#' @export
#' @importFrom stats runif
#' @importFrom stats rexp
#'
#'
get_nextEvent = function(current_age, disease_status, RV_status,
                         hazard_rates, GRR, carrier_prob,
                         lambda_birth, birth_range, fert = 1){

  RR <- ifelse(RV_status, GRR, 1)

  #set the number of subtypes to simulate
  num_subs <- length(hazard_rates$subtype_ID)

  event_type <- "retry"
  while (event_type == "retry"){
    #Multi-Subtype option
    t_onset <- ifelse(rep(disease_status, num_subs), rep(NA, num_subs),
                      sapply(1:num_subs, function(x){
                        get_wait_time(p = runif(1), last_event = current_age,
                                      hazard = hazard_rates[[1]][, x]*RR/(1 + carrier_prob*(GRR[x] - 1)),
                                      part = hazard_rates[[2]])
                      }))

    #simulate the waiting time until death given current age.
    t_death <- get_wait_time(p = runif(1),
                             last_event = current_age,
                             hazard = hazard_rates[[1]][, (2 + 1*disease_status)],
                             part = hazard_rates[[2]], scale = TRUE)

    # Want to adjust the waiting time until birth based on current age
    # and also ensure that birth cannot occur after the maximum birth age

    # First condition handles the case when person is less than the minimum
    # birth age when this occurs we return current age + waiting time
    # so long as this value is less than the maximum birth age
    # Second condition deals with people who have reached the minimum birth age,
    # if neiter is met we return NA
    nyear_birth <- rexp(1, lambda_birth*fert*disease_status + lambda_birth*(1 - disease_status))
    t_birth <- ifelse((current_age < birth_range[1] &
                         nyear_birth + birth_range[1] <= birth_range[2]),
                      nyear_birth + birth_range[1] - current_age,
                      ifelse((current_age >= birth_range[1] &
                                (nyear_birth + current_age) <= birth_range[2]),
                             nyear_birth,
                             NA))

    #create list of simulated times
    times <- c(t_birth, t_onset, t_death)

    #create event vector and corresponding subtype vector
    events <- c("Child", rep("Onset", num_subs), "Death")
    sub_lab <- c(NA, hazard_rates$subtype_ID, NA)



    #find the smallest waiting time
    min_time <- which.min(times)
    duplicate_times <- anyDuplicated(times, incomparables = NA)

    # check to see if there are any ties,
    # if not store winning event and, if
    # applicable, the subtype
    if (duplicate_times != 0) {
      event_type <- "retry"
    } else {
      t_event    <- times[min_time]
      event_type <- events[min_time]
      sub_ID     <- sub_lab[min_time]
    }
  }

  return(list(t_event, event_type, sub_ID))
}

#' Simulate all life events
#'
#' Primarily intended as an internal function, \code{sim_life} simulates all life events for an individual starting at birth, age 0, and ending with death or the end of the study.
#'
#' Starting at birth, age 0, \code{sim_life} generates waiting times to reproduction, onset, and death. The event with the shortest waiting time is chosen as the next life event, and the individual's age is updated by the waiting time of the winning event.  Conditioned on the individual's new age, this process is applied recursively, until death or until the end of the study is reached.
#'
#'  We make the following assumptions regarding the simulation of waiting times:
#'  \enumerate{
#'  \item We assume that, given an individual's current age, their time to disease onset is the waiting time in a non-homogeneous Poisson process with an age-specific hazard rate that follows a proportional hazards model.  In this model, individuals who have NOT inherited the rare variant experience disease onset according to the baseline (or population) hazard rate of disease.  On the other hand, individuals who have inherited the rare variant are assumed to have an increased risk of disease onset relative to those who have inherited it.  The user is expected to supply the baseline hazard rate of disease, as well as the relative-risk of disease for genetic cases. Additionally, we impose the restriction that individuals may only experience disease onset once, and remain affected from that point on.
#'
#'  \item We assume that, given an individual's current age, their time to death is the waiting time in a non-homogeneous Poisson process with age-specific hazard rate determined by their affection status.  We assume that disease-affected individuals experience death according to the age-specific hazard rate for death in the \emph{affected} population.  On the other hand, we assume that \emph{unaffected} individuals experience death according to the age-specific hazard rate for death in the \emph{unaffected} population.  If the disease of interest is sufficiently rare, the user may choose to substitute the \emph{population} age-specific hazard rate for death for the aforementioned age-specific hazard rate for death in the \emph{unaffected} population.  The user is expected to supply age-specific hazard rates of death for both the \emph{affected} and \emph{unaffected} populations.
#'
#'  \item We assume that, given an individual's current age, their time to reproduction is the waiting time in a homogeneous Poisson process.  That is, we assume that individuals reproduce at uniform rate during their reproductive years.  For example, one's reproductive years may span from age 20 to age 35 years.  To mimic observed age-specific fertility data, the birth range for an individual is simulated as follows: first we sample the lower bound uniformly from ages 16 to 27, next we sample the range of the birth span uniformly from 10 to 18 years and add this value to the lower bound to determine the upper bound of the birth range.  We do not allow for offspring to be produced outside of an individual's simulated reproductive birth span.
#'  }
#'
#'
#' The events simulated by \code{sim_life} are labelled as follows:
#' \itemize{
#'  \item "Start" the individual's year of birth.
#'  \item "Child" a reproductive event, i.e. creation of offspring
#'  \item "Onset" disease onset event,
#'  \item "Death" death event
#' }
#'
#'
#' @param YOB A positive number. The indivdiual's year of birth.
#' @param RV_status Numeric. \code{RV_status = TRUE} if the individual is a carrier of a rare variant that increases disease suseptibility, and \code{FALSE} otherwise.
#' @inheritParams sim_RVped
#'
#' @return an object of class \code{event}.  An object of class \code{event} is a list that contains the following items.
#' @return \item{\code{life_events} }{A named numeric vector of life events, see details.}
#' @return \item{\code{repro_events} }{A vector of reproduction years, that is the year(s) that the individual produces offspring.  When the individual does not reproduce \code{repro_events = NULL}.}
#' @return \item{\code{onset_event} }{Numeric. When applicable the year of disease-onset. When onset does not occur \code{onset_event = NA}}
#' @return \item{\code{death_event} }{Numeric. When applicable the year of death. When death is censored \code{death_event = NA}}
#' @return \item{\code{subtype} }{Character. When applicable, the disease subtype.}
#' @return \item{\code{censor_year} }{Numeric. When applicable the last year that data was observed. Note: after death \code{censor_year = NA}}

#'
#' @references Nieuwoudt, Christina and Jones, Samantha J and Brooks-Wilson, Angela and Graham, Jinko. (14 December 2017) \emph{Simulating Pedigrees Ascertained for Multiple Disease-Affected Relatives}. bioRxiv 234153.
#' @references Ken-Ichi Kojima, Therese M. Kelleher. (1962), \emph{Survival of Mutant Genes}. The American Naturalist 96, 329-346.
#' @export
#' @importFrom stats rgamma
#'
#' @examples
#' data(AgeSpecific_Hazards)
#' my_HR <- hazard(hazardDF = AgeSpecific_Hazards)
#'
#' # The following commands simulate all life events for an individual, who
#' # has NOT inherited a causal variant, born in 1900.  From the output, this
#' # individual has two children, one in 1921 and another in 1923, and then
#' # dies in 1987.
#' set.seed(135)
#' sim_life(hazard_rates = my_HR, GRR = 10,
#'          carrier_prob = 0.002,
#'          RV_status = FALSE,
#'          YOB = 1900, stop_year = 2000)
#'
#' # Using the same random seed, notice how life events can vary for
#' # someone who has inherited the causal variant, which carries a
#' # relative-risk of 10. From the output, this individual also has
#' # two children, but then experiences disease onset in 1974,
#' # and dies in 1976.
#' set.seed(135)
#' sim_life(hazard_rates = my_HR, GRR = 10,
#'                carrier_prob = 0.002,
#'                RV_status = TRUE,
#'                YOB = 1900, stop_year = 2000)
#'
#' set.seed(1)
#' sim_life(hazard_rates = hazard(AgeSpecific_Hazards[, c(1, 1, 2, 3)], subtype_ID = c("HL", "NHL")),
#'             GRR = 10,
#'                carrier_prob = 0.002,
#'                RV_status = TRUE,
#'                YOB = 1900, stop_year = 2000)
#'
#'
sim_life = function(hazard_rates, GRR, carrier_prob,
                    RV_status, YOB, stop_year,
                    NB_params = c(2, 4/7),
                    fert = 1, birth_range = NULL){

  if(!is.hazard(hazard_rates)) {
    stop("hazard_rates must be an object of class hazard")
  }

  if (!is.null(birth_range)) {
    warning("The argument birth_range has been deprecated. Execute help(sim_life) for details.")
  }

  if (any(GRR <= 0)) {
    stop ('GRR must be greater than 0')
  }

  if (fert <= 0) {
    stop ('fert must be greater than 0')
  }

  if (fert > 1) {
    warning ('fert > 1 detected. \n Are you sure you want to increase fertility after disease-onset?')
  }


  if (length(hazard_rates$subtype_ID) > 1 & length(GRR) == 1) {
    GRR = rep(GRR, length(hazard_rates$subtype_ID))
  } else if (length(GRR) != length(hazard_rates$subtype_ID)) {
    stop('length(GRR) != length(subtype_counts)\n please ensure that GRR contains one entry for each subtype')
  }


  #initialize lists to store event times and types
  R_life <- c(0)
  R_life_names  <- c("Start")
  min_age <- min(hazard_rates[[2]])
  max_age <- max(hazard_rates[[2]])
  #initialize disease status, start age at minumum age permissable under part
  DS <- F; t <- min_age; yr <- YOB

  #initialize events
  repro_events <- NULL
  death_event  <- NA
  onset_event  <- NA
  subtype      <- NA
  #onset_age    <- NA
  #death_age    <- NA

  #initialize event locations
  death_ID <- NA
  onset_ID <- NA
  repro_IDs <- c()


  #sample the minimum and maximum reproductive ages.
  B_range <- c(NA, NA)
  B_range[1] <- round(runif(1, min = 16, max = 27))
  B_range[2] <- B_range[1] + round(runif(1, min = 10, max = 18))

  #generate and store the birth rate for this individual
  B_lambda <- rgamma(1, shape = NB_params[1],
                     scale = (1-NB_params[2])/NB_params[2])/(B_range[2] -
                                                               B_range[1])

  counter <- 1
  while(t < max_age & yr <= stop_year){
    #generate next event
    l_event <- get_nextEvent(current_age = t, disease_status = DS, RV_status,
                             hazard_rates, GRR, carrier_prob,
                             lambda_birth = B_lambda, birth_range = B_range,
                             fert)
    counter <- counter + 1

    if(yr + l_event[[1]] <= stop_year){
      #add to previous life events
      R_life <- c(R_life, l_event[[1]])
      R_life_names <- c(R_life_names, l_event[[2]])

      if (l_event[[2]] == "Death") {
        #if death occurs stop simulation by setting t = 100
        t <- max_age
        death_ID <- counter
      } else if (l_event[[2]] == "Onset") {
        #if onset occurs change disease status and continue
        t  <- t + l_event[[1]]
        yr <- yr + l_event[[1]]
        DS <- T
        subtype <- l_event[[3]]
        onset_ID <- counter
      } else {
        #otherwise, handle birth events
        t <- t + l_event[[1]]
        yr <- yr + l_event[[1]]
        repro_IDs <- repro_IDs <- c(repro_IDs, counter)
      }
    } else {
      #if event is after stop date, increment yr to stop while loop
      yr <- yr + l_event[[1]]
    }

  }

  life_events <- round(cumsum(as.numeric(R_life)) + YOB)

  if (length(repro_IDs) > 0) repro_events <- life_events[repro_IDs]
  if (!is.na(onset_ID)) {
    onset_event <- life_events[onset_ID]
    onset_age <- onset_event - YOB
  }

  if (!is.na(death_ID)) {
    death_event <- life_events[death_ID]
    death_age <- death_event - YOB
  }

  names(life_events) <- R_life_names

  evt_obj <- list(life_events = life_events,
                  repro_events = repro_events,
                  onset_event = onset_event,
                  death_event = death_event,
                  subtype = subtype,
                  censor_year = stop_year)
  #onset_age = onset_age,
  #death_age = death_age)

  return(events(evt_obj))
}
