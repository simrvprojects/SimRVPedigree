#' Determine cumulative risk at a given time
#'
#'\code{approxCumHaz} returns the cumulative risk at time t, given a vector of age-specific hazards and a partition of ages over which to apply the age-specific hazards.
#'
#'@param t A numeric constant. The age at which to approximate the cumulative hazard
#'@param hazard A numeric vector.  A vector of age-specific hazards.
#'@param part A numeric vector.  Partition of ages over which to apply the age-specific hazards.
#'
#' @return cum_hazard numeric. The approximate cumulative hazard at time t.
#'
#' @examples
#' approxCumHaz(t = 5, hazard = seq(1, 2.5, by = 0.1), part = seq(0, 80, by = 5))
#' approxCumHaz(t = 5.5, hazard = seq(1, 2.5, by = 0.1), part = seq(0, 80, by = 5))
#' approxCumHaz(t = 85.5, hazard = seq(1, 2.5, by = 0.1), part = seq(0, 80, by = 5))
#'
approxCumHaz = function(t, hazard, part) {

  check_hazpart(hazard, part)

  CHazFun = approxfun(x = part, y = c(0, cumsum(hazard)))

  cum_hazard = ifelse((t >= min(part) & t <= max(part)), CHazFun(t),
                   ifelse(t > max(part), CHazFun(max(part)), 0))

  return(cum_hazard)
}

#' computes the probability that a waiting time, W, is less than or equal to some value s, given that the last event occured at s_0
#'
#' @param last_event A numeric constant.  The age at last event.
#' @param wait_time A numeric constant. The wiating time, in years, to next event.
#'@param hazard A numeric vector.  A vector of age-specific hazards.
#'@param part A numeric vector.  Partition of ages over which to apply the age-specific hazards.
#' @param scale Logical. By default scale = FALSE.  NOTE: Choosing scale = TRUE ensures that W is a proper random variable, i.e. that this function is a proper CDF with upper limit 1.
#'
#' @return wait_prob numeric. The probability that the waiting time is at most wait_time given that the last event occured at last_event
#'
#' @examples
#' findWaitProb(last_event = 2, wait_time = 18, hazard = seq(1, 2.5, by = 0.1), part = seq(0, 80, by = 5))
#' findWaitProb(last_event = 2, wait_time = 18, hazard = seq(0, 1.5, by = 0.1), part = seq(0, 80, by = 5), scale = TRUE)
#' findWaitProb(last_event = 2, wait_time = 18, hazard = seq(0, 1.5, by = 0.1), part = seq(0, 80, by = 5))
#'
findWaitProb = function(last_event, wait_time,
                        hazard, part, scale = FALSE) {

  prob <- 1 - exp(approxCumHaz(t = last_event, hazard, part)
                  - approxCumHaz(t = last_event + wait_time, hazard, part))

  #uplimit <- 1 - exp(-approxCumHaz(t = part[length(part)], hazard, part))

  uplimit <- 1 - exp(approxCumHaz(t = last_event, hazard, part)
                     -approxCumHaz(t = part[length(part)], hazard, part))

  wait_prob = ifelse(scale == FALSE, prob, prob/uplimit)

  return(wait_prob)
}

##----------------##
##  findWaitTime  ##
##----------------##
## deine a function that finds the waiting time for a given probability
## i.e. create a function that obtains solutions to the inverse CDF numerically

## Arguments--------------------------------------------------------------------
## u          - constant
## last_event - constant
## hazard     - numeric vector of length n
## part  - numeric vector of length n + 1
## scale      - logical (T/F)

## Function Requirements--------------------------------------------------------
## approxCumHaz

## Package Requirements---------------------------------------------------------
## NONE

#Define a function that will obtain the solution to the inverse CDF numerically
findWaitTime = function(u, last_event,
                        hazard, part, scale = FALSE){

  #Find the maximum value of findWaitProb so that we
  # can quickly return NA when u is greater than this value

  MaxProb <- findWaitProb(last_event, wait_time = max(part), hazard, part, scale)

  if (u <= MaxProb) {
    a    <- min(part)
    b    <- max(part) - last_event
    diff <- 1

    while (diff > 0.001) {
      w <- (a + b)/2

      if (findWaitProb(wait_time = w, last_event, hazard, part, scale) <= u){
        a <- w
      } else if (findWaitProb(wait_time = w, last_event, hazard, part, scale) > u){
        b <- w
      }

      diff <- abs(w - (a + b)/2)
    }

  } else {
    w <- NA
  }

  return(w)

}
