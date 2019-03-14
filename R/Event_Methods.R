#' Constructor function for an object of class events
#'
#' @param life_events The list of items returned by the \code{sim_life} function.
#'
#' @return an object of class \code{events}.
#' @keywords internal
events <- function(life_events) {
  class(life_events) <- c("events", class(life_events))
  return(life_events)
}

#' Check to see if object is of class ped
#'
#' @param x An R object.
#'
#' @return Logical. Indicates if \code{x} is of class \code{ped}.
#'
#' @keywords internal
is.events <- function(x) {
  return(inherits(x, "events"))
}
