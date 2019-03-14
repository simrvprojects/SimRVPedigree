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


#' @export
print.events <- function(x, ...) {
  cat("Life events starting at year-of-birth: ", x$life_events[1], "\n")
  cat("\n")
  print(x$life_events)
  cat("\n")
  if ( !is.na(x$onset_event) ) {
    if ( x$subtype != "no_subtypes" ) {
      cat("Onset of disease-subtype", x$subtype, "at age", x$onset_event - as.numeric(x$life_events[1]), "\n")
    } else {
      cat("Onset of disease at age", x$onset_event - as.numeric(x$life_events[1]), "\n")
    }

  }
}
