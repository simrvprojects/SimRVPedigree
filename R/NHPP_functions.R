#' Determine the cumulative probability of the waiting time, given the age at last event
#'
#' @param last_event A numeric constant.  The age at last event.
#' @param wait_time A numeric constant.  The waiting time to next_event.
#' @param hazard A numeric vector.  A vector of age-specific hazards.
#' @param part A numeric vector.  Partition of ages over which to apply the
#' age-specific hazards.
#' @param scale Logical. By default scale = FALSE.  Specifying scale = TRUE ensures that W is a proper random variable, i.e. that this function is a proper CDF with upper limit 1.
#'
#' @return wait_prob The probability that the waiting time to next event is at least wait_time.
#' @export
#'
#' @importFrom stats approxfun
#'
#' @examples
#' haz_vec <- c(seq(0, 0.5, by = 0.05), rev(seq(0.46, 0.5, by = 0.01)))
#' part_vec <- seq(0, 80, by = 5)
#'
#' get_WaitProb(last_event = 0, wait_time = 18,
#'              hazard = haz_vec, part = part_vec,
#'              scale = FALSE)
#' ##Tests:
#' ##wait_prob in [0, 1]
#' ##when scale = TRUE, maxProb = 1
get_WaitProb = function(last_event, wait_time,
                        hazard, part, scale = FALSE) {

  approxCH <- approxfun(x = part, y = c(0, cumsum(hazard)), rule = 2)

  CumProb <- 1 - exp(approxCH(last_event)
                     - approxCH(last_event + wait_time))

  uplimit <- 1 - exp(approxCH(last_event)
                     -approxCH(part[length(part)]))

  if(scale == FALSE){
    wait_prob <- CumProb
  } else {
    wait_prob <- CumProb/uplimit
  }

  return(wait_prob)
}


#' Simulate the waiting time to next event for a non-homogeneous Poisson process
#'
#' \code{get_WaitTime} simulates the waiting time to next event for a non-homogeneous Poisson process
#'
#' \code{get_WaitTime} simulates the waiting time to next event for a non-homogeneous Poisson process.  The units of the simulated waiting time are the units specified in \code{part}, i.e. if \code{part} is specified in years, the simulated waiting time is in years.
#'
#' @param p A numeric constant. Argument of non-homogeneous poisson process quantile function
#' @param last_event A numeric constant.  The age at last event.

#' @param hazard A numeric vector.  A vector of age-specific hazards.
#' @param part A numeric vector.  Partition of ages over which to apply the
#' age-specific hazards.
#' @param scale Logical. By default scale = FALSE.  Specifying scale = TRUE ensures that W is a proper random variable, i.e. that this function is a proper CDF with upper limit 1.
#'
#' @return The waiting time to next event, units same as those in \code{part}.
#' @export
#'
#' @examples
#' haz_vec <- c(seq(0, 0.5, by = 0.05), rev(seq(0.46, 0.5, by = 0.01)))
#' part_vec <- seq(0, 80, by = 5)
#'
#' get_WaitTime(p = 0.05, last_event = 0,
#'              hazard = haz_vec, part = part_vec,
#'              scale = FALSE)
#'
#' get_WaitTime(p = 0.05, last_event = 0,
#'              hazard = haz_vec, part = part_vec,
#'              scale = TRUE)
#' ##Tests:
#' ##wait_time in [0, max(part) - last_event]
#' ##when scale = TRUE, wait_time != NA
get_WaitTime = function(p, last_event, hazard, part,
                         scale = FALSE){

  #Find the maximum value of findWaitProb so that we
  # can quickly return NA when u is greater than this value

  MaxProb <- get_WaitProb(last_event, wait_time = max(part),
                          hazard, part, scale)

  if (p > MaxProb) {
    w <- NA
  } else if (p < MaxProb) {
    wts   <- seq(0, (max(part) - last_event), by = 0.5)
    probs <- get_WaitProb(wait_time = wts, last_event, hazard, part, scale)

    a <- wts[(which.max((p - probs) < 0)-1)]
    b <- wts[which.max((p - probs) < 0)]
    w <- (a + b)/2
  } else if (p == MaxProb) {
    w <- max(part) - last_event
  }
  return(w)
}

