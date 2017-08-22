#' Obtain the cumulative probability of the waiting time to the next event
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
#'
#' @keywords internal
#'
#' @importFrom stats approxfun
get_wait_prob = function(last_event, wait_time,
                         hazard, part, scale = FALSE) {

  approxCH <- approxfun(x = part, y = c(0, cumsum(hazard)), rule = 2)

  cum_prob <- 1 - exp(approxCH(last_event)
                     - approxCH(last_event + wait_time))

  upper_limit <- 1 - exp(approxCH(last_event)
                     -approxCH(part[length(part)]))

  if(scale){
    wait_prob <- cum_prob/upper_limit
  } else {
    wait_prob <- cum_prob
  }

  return(wait_prob)
}


#' obtain the waiting time to the next event
#'
#' \code{get_wait_time} approximates the result of the inverse cumulative distribution function of the waiting time to the next event associated with a non-homogeneous Poisson process conditioned of the time of the last event.
#'
#' \code{get_wait_time} obtains the waiting time to next event associated with a non-homogeneous Poisson process.  The units of the simulated waiting time are the units specified in \code{part}, i.e. if \code{part} is specified in years, the simulated waiting time is in years.
#'
#' @inheritParams get_wait_prob
#' @param p Numeric.  A probability; i.e. the argument of the quantile function of the waiting time associated with a non-homogeneous Poisson process with rate \code{hazard}.
#' @param scale Logical. By default \code{scale = FALSE}.  Specifying \code{scale = TRUE} ensures that this function is the quantile function for a proper random variable.
#'
#' @return The waiting time to next event, units same as those in \code{part}.
#'
#' @keywords internal
get_wait_time = function(p, last_event, hazard, part,
                         scale = FALSE){

  #Find the maximum value of findWaitProb so that we
  # can quickly return NA when u is greater than this value
  max_age <- max(part)

  max_prob <- get_wait_prob(last_event, wait_time = max_age,
                            hazard, part, scale)

  if (p > max_prob) {
    w <- NA
  } else if (p < max_prob) {
    wts <- seq(0, (max_age - last_event), by = 0.5)
    if (max(wts) < (max_age - last_event)){
      wts <- c(wts, (max_age - last_event))
    }
    probs <- get_wait_prob(wait_time = wts, last_event, hazard, part, scale)

    a <- wts[(which.max((p - probs) < 0)-1)]
    b <- wts[which.max((p - probs) < 0)]
    w <- (a + b)/2
  } else if (p == max_prob) {
    w <- max_age - last_event
  }
  return(w)
}

