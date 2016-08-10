#' Determine cumulative risk at a given time
#'
#' \code{approxCumHaz} returns the cumulative risk at time t, given a vector of
#' age-specific hazards and a partition of ages over which to apply the
#' age-specific hazards.
#'
#' @param t A numeric constant. The age at which to approximate the cumulative
#' hazard
#' @param hazard A numeric vector.  A vector of age-specific hazards.
#' @param part A numeric vector.  Partition of ages over which to apply the
#' age-specific hazards.
#'
#' @return cum_hazard numeric. The approximate cumulative hazard at time t.
#'
#' @examples
#' haz_vec <- seq(1, 2.5, by = 0.1)
#' part_vec <- seq(0, 80, by = 5)
#'
#' approxCumHaz(t = 0, hazard = haz_vec, part = part_vec)
#' approxCumHaz(t = 5, hazard = haz_vec, part = part_vec)
#' approxCumHaz(t = 17, hazard = haz_vec, part = part_vec)
#' approxCumHaz(t = 80, hazard = haz_vec, part = part_vec)
#' approxCumHaz(t = 100, hazard = haz_vec, part = part_vec)
#'
approxCumHaz = function(t, hazard, part) {
  check_hazpart(hazard, part)

  CHazFun = approxfun(x = part, y = c(0, cumsum(hazard)), rule = 2)

  return(cum_hazard = CHazFun(t))
}


#' computes the probability that the waiting time, is at least some value given the time of the last event
#'
#' @param last_event A numeric constant.  The age at last event.
#' @param wait_time A numeric constant. The wiating time, in years, to next event.
#' @param scale Logical. By default scale = FALSE.  Specifying scale = TRUE ensures that W is a proper random variable, i.e. that this function is a proper CDF with upper limit 1.
#' @inheritParams approxCumHaz
#'
#' @return wait_prob numeric. The probability that the waiting time is at most wait_time given that the last event occured at last_event
#'
#' @examples
#' haz_vec <- seq(1, 2.5, by = 0.1)
#' part_vec <- seq(0, 80, by = 5)
#'
#' findWaitProb(last_event = 2, wait_time = 18,
#'              hazard = haz_vec, part = part_vec)
#'
#' findWaitProb(last_event = 10, wait_time = 18,
#'              hazard = haz_vec, part = part_vec)
#'
#' findWaitProb(last_event = 10, wait_time = 18,
#'              hazard = haz_vec, part = part_vec,
#'              scale = TRUE)

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

#' Randomly generate the waiting time to next event, given the last event time
#'
#' @param u A numeric constant. argument of inverse CDF
#' @inheritParams findWaitProb
#'
#' @return A numeric constant.  The waiting time to next event.
#'
#' @examples
#' haz_vec <- c(seq(0, 0.5, by = 0.05), rev(seq(0.46, 0.5, by = 0.01)))
#' part_vec <- seq(0, 80, by = 5)
#'
#' set.seed(17)
#' u_val = runif(1)
#' findWaitTime(u = u_val, last_event = 2,
#'              hazard = haz_vec, part = part_vec)
#'
#' findWaitTime(u = u_val, last_event = 2,
#'              hazard = haz_vec, part = part_vec,
#'              scale = TRUE)
#'
#' findWaitTime(u = u_val, last_event = 10,
#'              hazard = haz_vec, part = part_vec,
#'              scale = TRUE)
#'
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

#' Simulate the waiting time to next event for a non-homogeneous Poisson process
#'
#' \code{get_WaitTime} simulates the waiting time to next event for a non-homogeneous Poisson process
#'
#' \code{get_WaitTime} simulates the waiting time to next event for a non-homogeneous Poisson process.  The units of the simulated waiting time are the units specified in \code{part}, i.e. if \code{part} is specified in years, the simulated waiting time is in years.  Makes use of \link{approxfun} and \link{splinefun}
#'
#' @param p A numeric constant. Argument of non-homogeneous poisson process quantile function
#' @param last_event A numeric constant.  The age at last event.

#' @param hazard A numeric vector.  A vector of age-specific hazards.
#' @param part A numeric vector.  Partition of ages over which to apply the
#' age-specific hazards.
#' @param scale Logical. By default scale = FALSE.  Specifying scale = TRUE ensures that W is a proper random variable, i.e. that this function is a proper CDF with upper limit 1.
#'
#' @return The waiting time to next event, units same as those in \code{part}.
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
#'

get_WaitTime = function(p, last_event, hazard, part, scale = FALSE){
  check_hazpart(hazard, part)
  check_part(part)

  #create a function to approximate cumulative hazard using approxfun()
  approxCumHaz <- approxfun(x = part, y = c(0, cumsum(hazard)), rule = 2)

  #find F_T(x) = P(wait time <= x| last event occurred at time T), denote "CumProb"
  CumProb <- 1 - exp(approxCumHaz(last_event) - approxCumHaz(last_event + part))

  #Find maximum value of CDF,
  uplimit <- 1 - exp(approxCumHaz(last_event) - approxCumHaz(part[length(part)]))

  #rescale the CDF so that it's upper limit is 1 when scale==TRUE
  if(scale == TRUE){
    CumProb <- CumProb/uplimit
  }

  #use approxfun to approximate the inverse CDF
  approxInvCDF = splinefun(x = CumProb, y = part)

  wait_time = ifelse( approxInvCDF(p) > (part[length(part)]-last_event),
                      (part[length(part)]-last_event),
                      approxInvCDF(p))

  #check to see if u is greater than the maximum value of CumProb, if so
  #return NA
  wait_time = ifelse(p <= max(CumProb), wait_time, NA)

  return(wait_time)
}


