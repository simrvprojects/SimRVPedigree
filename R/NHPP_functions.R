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
#' @export
#'
#' @importFrom stats splinefun
#' @importFrom stats approxfun
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

get_WaitTime = function(p, last_event, hazard, part,
                        scale = FALSE){

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


