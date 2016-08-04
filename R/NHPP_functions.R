##-----------------##
##  check_hazpart  ##
##-----------------##
## define a function that checks to make sure all the appropriate arguments 
##  have been provided in NHPP functions

## Arguments--------------------------------------------------------------------
## hazard    - numeric vector of length n
## part - numeric vector of length n + 1

## Function Requirements--------------------------------------------------------
## NONE

## Package Requirements---------------------------------------------------------
## NONE

check_hazpart = function(hazard, part){
  check1 <- (class(hazard) != "numeric")
  check2 <- (length(hazard) < 1)
  check3 <- (length(part) != (length(hazard) + 1))
  if ( check1 | check2 | check3 ) {
    stop ('please provide numeric hazard and part vectors, with length(part) = length(hazard) + 1')
  }
}

##-----------------##
## approxHazFun ##
##-----------------##
## define a functon that creates a cumulative hazard function  
## that linearly interpolates the cumulative risk within a time interval  

## Arguments--------------------------------------------------------------------
## hazard    - numeric vector of length n
## part - numeric vector of length n + 1

## Function Requirements--------------------------------------------------------
## check_hazpart

## Package Requirements---------------------------------------------------------
## NONE
approxHazFun = function(hazard, part){
  check_hazpart(hazard, part)
  return (approxfun(x = part, y = c(0, hazard)))
}


##-----------------##
## approxCumHazFun ##
##-----------------##
## define a functon that creates a cumulative hazard function  
## that linearly interpolates the cumulative risk within a time interval  

## Arguments--------------------------------------------------------------------
## hazard    - numeric vector of length n
## part - numeric vector of length n + 1

## Function Requirements--------------------------------------------------------
## check_hazpart

## Package Requirements---------------------------------------------------------
## NONE
approxCumHazFun = function(hazard, part){
  check_hazpart(hazard, part)
  return (approxfun(x = part, y = c(0, cumsum(hazard))))
}


##----------------##
##  approxCumHaz  ##
##----------------##
## define a function returns the cumulative risk at time t, given
## a vector of hazards and a partition over which to apply hazards.  
## NOTE: the difference between approxCumHaz and approxCumHazFun is that 
##       approxCumHaz returns a constant while approxCumHazFun returns a function

## Arguments--------------------------------------------------------------------
## t         - constant 
## hazard    - numeric vector of length n
## part - numeric vector of length n + 1

## Function Requirements--------------------------------------------------------
## approxCumHazFun

## Package Requirements---------------------------------------------------------
## NONE
approxCumHaz = function(t, hazard, part) {
  
  CHazFun = approxCumHazFun(hazard, part)
  
  my.CHaz = ifelse((t >= min(part) & t <= max(part)), CHazFun(t),
                   ifelse(t > max(part), CHazFun(max(part)), 0))
  
  return(my.CHaz)
}


##----------------##
##  findWaitProb  ##
##----------------##
## deine a function that to computes the probability that a waiting time, W, 
## is less than or equal to some value s, given that the last event occured at 
## s_0
## NOTE: Choosing scale = TRUE ensures that W is a proper random variable,
## i.e. that this function is a proper CDF with upper limit 1.  

## Arguments--------------------------------------------------------------------
## last_event - constant 
## next_event - constant
## hazard     - numeric vector of length n
## part  - numeric vector of length n + 1
## scale      - logical (T/F)

## Function Requirements--------------------------------------------------------
## approxCumHaz

## Package Requirements---------------------------------------------------------
## NONE

findWaitProb = function(last_event, next_event, 
                        hazard, part, scale = FALSE) {
  
  prob <- 1 - exp(approxCumHaz(t = last_event, hazard, part)
                  - approxCumHaz(t = last_event + next_event, hazard, part))
  
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
  
  MaxProb <- findWaitProb(last_event, next_event = max(part), hazard, part, scale)
  
  if (u <= MaxProb) {
    a    <- min(part)
    b    <- max(part) - last_event
    diff <- 1
    
    while (diff > 0.001) {
      w <- (a + b)/2
      
      if (findWaitProb(next_event = w, last_event, hazard, part, scale) <= u){
        a <- w
      } else if (findWaitProb(next_event = w, last_event, hazard, part, scale) > u){
        b <- w
      }
      
      diff <- abs(w - (a + b)/2)
    }
    
  } else {
    w <- NA
  }
  
  return(w)

}
