#' Obtain the Cumulative Probability of the Waiting Time to the Next Event
#'
#' Determine the cumulative probability of the waiting time associated with a non-homogenous Poisson process, given the time of the last event.
#'
#' @param last_event A positive number.  The time at last event.
#' @param wait_time A positive number.  The waiting time to the next event.
#' @param hazard A vector of positive numbers.  The time-specific hazard rate.
#' @param part A vector of positive numbers.  Time partition over which to apply the time-specific hazard rate.
#' @param scale Logical. By default \code{scale = FALSE}.  Specifying \code{scale = TRUE} ensures that this function is the CDF of a proper random variable.
#'
#' @return \code{wait_prob} The probability that the waiting time to next event is at least \code{wait_time}.
#' @export
#'
#' @importFrom stats approxfun
#'
#' @examples
#' data(AgeSpecific_Hazards)
#' haz_vec <- AgeSpecific_Hazards[,2]
#' part_vec <- seq(0, 100, by = 1)
#'
#' get_WaitProb(last_event = 0, wait_time = 18,
#'              hazard = haz_vec, part = part_vec,
#'              scale = FALSE)
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


#' Obtain the Waiting Time to the Next Event
#'
#' \code{get_WaitTime} approximates the result of the inverse cumulative distribution function of the waiting time to the next event associated with a non-homogeneous Poisson process conditioned of the time of the last event.
#'
#' \code{get_WaitTime} obtains the waiting time to next event associated with a non-homogeneous Poisson process.  The units of the simulated waiting time are the units specified in \code{part}, i.e. if \code{part} is specified in years, the simulated waiting time is in years.
#'
#' @inheritParams get_WaitProb
#' @param p Numeric.  A probability; i.e. the argument of the quantile function of the waiting time associated with a non-homogeneous Poisson process with rate \code{hazard}.
#' @param scale Logical. By default \code{scale = FALSE}.  Specifying \code{scale = TRUE} ensures that this function is the quantile function for a proper random variable.
#'
#' @return The waiting time to next event, units same as those in \code{part}.
#' @export
#'
#' @examples
#' data(AgeSpecific_Hazards)
#' haz_vec <- AgeSpecific_Hazards[,1]
#' part_vec <- seq(0, 100, by = 1)
#'
#' get_WaitTime(p = 0.05, last_event = 0,
#'              hazard = haz_vec, part = part_vec,
#'              scale = FALSE)
#'
#' get_WaitTime(p = 0.05, last_event = 0,
#'              hazard = haz_vec, part = part_vec,
#'              scale = TRUE)
get_WaitTime = function(p, last_event, hazard, part,
                         scale = FALSE){

  #Find the maximum value of findWaitProb so that we
  # can quickly return NA when u is greater than this value

  MaxProb <- get_WaitProb(last_event, wait_time = max(part),
                          hazard, part, scale)

  if (p > MaxProb) {
    w <- NA
  } else if (p < MaxProb) {
    wts <- seq(0, (max(part) - last_event), by = 0.5)
    if (max(wts) < (max(part) - last_event)){
      wts <- c(wts, (max(part) - last_event))
    }
    probs <- get_WaitProb(wait_time = wts, last_event, hazard, part, scale)

    a <- wts[(which.max((p - probs) < 0)-1)]
    b <- wts[which.max((p - probs) < 0)]
    w <- (a + b)/2
  } else if (p == MaxProb) {
    w <- max(part) - last_event
  }
  return(w)
}

