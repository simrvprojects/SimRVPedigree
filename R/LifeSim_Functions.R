##---------------##
##  rbirth_rate  ##
##---------------##
## create function to simulate a random birth rate for individual who will
## reproduce from min_birth_age to max_birth_age
##
## Arguments--------------------------------------------------------------------
## NB_size        - constant; size parameter of NB distrbution
## NB_prob        - constant; probability parameter of NB distrbution 
## min_birth_age  - constant
## max_birth_age  - constant

## Function Requirements--------------------------------------------------------
## NONE

## Package Requirements---------------------------------------------------------
## NONE
rbirth_rate = function(NB_size, NB_prob, min_birth_age, max_birth_age){
  birth_rate = rgamma(1, shape = NB_size, 
         scale = (1-NB_prob)/NB_prob ) / (max_birth_age-min_birth_age)
}

##---------------##
##  event_step   ##
##---------------##
## create a function that determines the next life event: onset, death, or birth
## by treating these 3 events as competing life events
##
## Arguments--------------------------------------------------------------------
## current_age    - constant
## disease_status - constant; 1 if individuals has ever been affected, 0 otherwise
## RV_status      - constant; 1 if individual has RV, 0 otherwise
## lambda_birth   - constant; the individuals simulated birth rate
## onset_hazard   - vector; age specific incidence rates for the disease of interest
## death_hazard   - data.frame; age specific mortality rates
##                              column 1 = unaffected mortality rates
##                              column 2 = affected mortality rates
## part           - vector; partition over which to apply onset and death rates
## birth_range    - vector; max and min birth ages
## risk           - constant: the RR for developing disease
## Function Requirements--------------------------------------------------------
## findWaitTime

## Package Requirements---------------------------------------------------------
## NONE

event_step = function(current_age, disease_status, RV_status, 
                      lambda_birth, onset_hazard, death_hazard, part, 
                      birth_range, risk){
  
  #first pull appropriate columns from death_hazard and onset_hazard for given 
  # RV status and disease status
  risk.lambda  <- ifelse(RV_status == 1, risk, 1)
  onset.lambda <- onset_hazard*risk.lambda
  death.lambda <- death_hazard[, (1 + disease_status)]
  
  #Assuming that the person is not yet affected, simulate the waiting time until onset given current age
  t.onset <- ifelse(disease_status == 0, 
                    findWaitTime(u = runif(1), last_event = current_age, hazard = onset.lambda, part), 
                    NA)
  
  #simulate the waiting time until death given current age.  
  # NOTE: choosing rate = FASLE implies that we are assuming person will die 
  t.death <- findWaitTime(u = runif(1), last_event = current_age, hazard = death.lambda, part, scale = TRUE)
  
  # Want to adjust the waiting time until birth based on current age 
  # and also ensure that birth cannot occur after the maximum birth age
  nyear.birth <- rexp(1, lambda_birth)
  t.birth <- ifelse(
    # this condition handles the case when person is less than the minimum 
    # birth age when this occurs we return current age + waiting time 
    # so long as this value is less than the maximum birth age
    ((current_age < birth_range[1]) & (nyear.birth + birth_range[1] <= birth_range[2])), 
    nyear.birth + birth_range[1],
    # next we deal with people who have reached the minimum birth age, 
    ifelse(((current_age >= birth_range[1]) & (nyear.birth + current_age <= birth_range[2])),
           nyear.birth,
           #if neither condition is met we return NA
           NA))
  
  #create list of simulated times
  times <- as.list(c(Birth = t.birth,  Onset = t.onset, Death = t.death))
  
  #find the smallest waiting time
  if (sum(which.min(times)) == 0) {
    nyears = as.matrix(0, ncol = 1)
    colnames(nyears) = c("retry")
  }else if (sum(which.min(times)) != 0) {
    nyears = as.matrix(times[[which.min(times)]], ncol = 1)
    colnames(nyears) = c(paste(names(which.min(times))))
  }
  return(nyears)
}

##---------------##
##  life_step   ##
##---------------##
## create a function that determines all life events
##
## Arguments--------------------------------------------------------------------
## RV_status      - constant; 1 if individual has RV, 0 otherwise
## lambda_birth   - constant; the individuals simulated birth rate
## onset_hazard   - vector; age specific incidence rates for the disease of interest
## death_hazard   - data.frame; age specific mortality rates
##                              column 1 = unaffected mortality rates
##                              column 2 = affected mortality rates
## part           - vector; partition over which to apply onset and death rates
## birth_range    - vector; max and min birth ages
## NB_params      - vector; size and probability parameters for NB distribution
## risk           - constant: the RR for developing disease
## Function Requirements--------------------------------------------------------
## rbirth_rate
## event_step

## Package Requirements---------------------------------------------------------
## NONE
life_step = function(RV_status, onset_hazard, death_hazard, part,
                     birth_range, NB_params, risk){
  
  #initialize data frame to hold life events
  R.life  <- data.frame(Start = 0) 
  min.age <- min(part)
  max.age <- max(part)
  #initialize disease status, start age at minumum age permissable under part
  DS <- 0; t <- min.age 
  
  #generate and store the birth rate for this individual
  B.lambda <- rbirth_rate( NB_params[1], NB_params[2] , 
                           birth_range[1], birth_range[2])
  
  while(t < max.age){
    #generate next event
    my.step <- event_step(current_age = t, disease_status = DS, RV_status, 
                          lambda_birth = B.lambda, onset_hazard, death_hazard, 
                          part, birth_range, risk)
    
    #add to previous life events
    R.life <- cbind(R.life, my.step)
    
    if (colnames(my.step) == "Death") {
      #if death occurs stop simulation by setting t = 100
      t <- max.age
    } else if (colnames(my.step) == "Onset") {
      #if onset occurs change disease status and continue
      t  <- t + as.numeric(my.step[1,1])
      DS <- 1
    } else {
      #otherwise, increment counter, handles both birth event and case when no 
      # event occurs (i.e. retry)
      t <- t + as.numeric(my.step[1,1])
    }
  }
  
  return(R.life)
}
