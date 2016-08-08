##------------------##
##  create_pedFile  ##
##------------------##
## create a function returns an empty data frame with all of the
## appropriate ped.file variables
##
## Arguments____________________________________________________________________
## family_num   - constant, family identifier
## risk         - constant, relative risk of developing disease
##
## Function Requirements________________________________________________________
## NONE
##
## Package Requirements_________________________________________________________
## NONE
##

create_pedFile = function(){
  data.frame(FamID = numeric(),
             ID = numeric(),
             gender = numeric(),
             dad_id = numeric(),
             mom_id = numeric(),
             affected = numeric(),
             DA1 = numeric(),
             DA2 = numeric(),
             Gen = numeric(),
             birth_year = numeric(),
             onset_year = numeric(),
             death_year = numeric(),
             RR = numeric(),
             do_sim = numeric(),
             available = numeric(),
             stringsAsFactors=FALSE)
}

##------------##
##  add_mate  ##
##------------##
## create a function that adds all the necessary info for a mate to the pedigree
##
## Arguments____________________________________________________________________
## partner_info   - vector; info for partner
## last_id        - constant:
##
## Function Requirements________________________________________________________
## NONE
##
## Package Requirements_________________________________________________________
## NONE
##
add_mate = function(partner_info, last_id){
  new.mate.info = data.frame(FamID = partner_info$FamID,
                             ID = last_id+1,
                             gender = abs(partner_info$gender-1),
                             dad_id = NA,
                             mom_id = NA,
                             affected = 0,
                             DA1 = 0,
                             DA2 = 0,
                             Gen = partner_info$Gen,
                             birth_year = NA,
                             onset_year = NA,
                             death_year = NA,
                             RR = NA,
                             do_sim = 0,
                             available = 0,
                             stringsAsFactors=FALSE)

  mate.return = list(new.mate.info, last_id+1)
  return(mate.return)
}

##-----------------##
##  add_offspring  ##
##-----------------##
## create a function that adds all the necessary info for a child to the pedigree
##
## Arguments____________________________________________________________________
## dad_info   - vector; info for dad
## mom_info   - vector; info for mom
## last_id    - constant:
## byear      - constant; individual's birth year
##
## Function Requirements________________________________________________________
## NONE
##
## Package Requirements_________________________________________________________
## NONE
##
add_offspring = function(dad_info, mom_info, byear, last_id){
  new.child.info = data.frame(FamID = dad_info$FamID,
                              ID = last_id+1,
                              gender = round(runif(1)),
                              dad_id = dad_info$ID,
                              mom_id = mom_info$ID,
                              affected = 0,
                              DA1 = sample(x = c(dad_info$DA1, dad_info$DA2), size = 1),
                              DA2 = sample(x = c(mom_info$DA1, mom_info$DA2), size = 1),
                              Gen = dad_info$Gen + 1,
                              birth_year = byear,
                              onset_year = NA,
                              death_year = NA,
                              RR = c(dad_info$RR, mom_info$RR)[which(!is.na(c(dad_info$RR, mom_info$RR)))],
                              do_sim = 1,
                              available = 1,
                              stringsAsFactors=FALSE)
  child.return = list(new.child.info, last_id+1)
  return(child.return)
}

##-------------##
##  date_step  ##
##-------------##
## create a function to translate the waiting times between
## life events into the years at which the life events occur
##
## Arguments____________________________________________________________________
## dad_info   - vector; info for dad
## mom_info   - vector; info for mom
## last_id    - constant:
## byear      - constant; individual's birth year
##
## Function Requirements________________________________________________________
## NONE
##
## Package Requirements_________________________________________________________
## NONE
##
date_step = function(Rlife_events, byear){
  life.years = round(cumsum(as.numeric(Rlife_events[1,]))) + byear
  names(life.years) = names(Rlife_events)
  return(life.years)
}

##-------------##
##  nfam_step  ##
##-------------##
## create a function that simulates an individuals life events and then
## creates a nuclear family when appropriate
##
## Arguments____________________________________________________________________
## found_info     - vector; info for dad
## stop_year      - constant:
## last_id        - constant:
## onset_hazard   - vector; age specific incidence rates for the disease of interest
## death_hazard   - data.frame; age specific mortality rates
## part           - vector; partition over which to apply onset and death rates
## birth_range    - vector; max and min birth ages
## NB_params      - vector; size and probability parameters for NB distribution
## risk           - constant; the RR for developing disease
##
## Function Requirements________________________________________________________
## date_step
## life_step
## add_offspring
## add_mate
##
## Package Requirements_________________________________________________________
## NONE
##
nfam_step = function(found_info, stop_year, last_id,
                     onset_hazard, death_hazard, part,
                     birth_range, NB_params, risk){

  nfam_ped = found_info

  #Simulate life steps for our founder
  sim.life = life_step(RV_status = sum(c(found_info$DA1, found_info$DA2)),
                       onset_hazard, death_hazard, part,
                       birth_range, NB_params, risk)

  #convert these to years
  sim.years = date_step(Rlife_events = sim.life, byear = found_info$birth_year)

  #set the disease status based on whether or not they experience disease onset prior to stop_year
  onset.true = is.element("Onset", names(sim.years))
  if (onset.true == TRUE) {
    o.year = as.numeric(sim.years[which(names(sim.years) == "Onset")])
    if (o.year <= stop_year) {
      nfam_ped$affected = 1
      nfam_ped$onset_year = o.year
    } else {
      nfam_ped$affected = 0
      nfam_ped$onset_year = NA
    }
  } else {
    nfam_ped$affected = 0
    nfam_ped$onset_year = NA
  }

  #set the year of death if it occurs before stop_year
  d.year = as.numeric(sim.years[which(names(sim.years) == "Death")])
  if (d.year <= stop_year) {
    nfam_ped$death_year = d.year
  }

  #store the birth years of each child
  birth.events = as.numeric(sim.years[which(names(sim.years) == "Birth" &
                                              sim.years <= stop_year)])

  if (length(birth.events) > 0) {
    #add a mate
    new.mate = add_mate(partner_info = nfam_ped[1,], last_id)
    nfam_ped = rbind(nfam_ped, new.mate[[1]])
    last_id = new.mate[[2]]

    #store info for mom and dad
    dad = nfam_ped[which(nfam_ped$gender == 0), ]
    mom = nfam_ped[which(nfam_ped$gender == 1), ]

    for (k in 1:length(birth.events)) {
      #add child
      new.child = add_offspring(dad_info = dad, mom_info = mom,
                            byear = birth.events[k], last_id)
      nfam_ped = rbind(nfam_ped, new.child[[1]])
      last_id = new.child[[2]]
    }
  }

  #swich the do_sim value to 0 for the individual whose life events we have simulated
  nfam_ped$do_sim[1] = 0

  sim.fam.return = list(nfam_ped, last_id)
  return(sim.fam.return)
}


##-------------##
##  ped_step  ##
##-------------##
## create a function that simulates an entire pedigree up to the stop year
##
## Arguments____________________________________________________________________
## onset_hazard   - vector; age specific incidence rates for the disease of interest
## death_hazard   - data.frame; age specific mortality rates
## part           - vector; partition over which to apply onset and death rates
## birth_range    - vector; max and min birth ages
## NB_params      - vector; size and probability parameters for NB distribution
## risk           - constant; the RR for developing disease
## founder_byears - vector; length 2, years to be chosen from uniformly for founder birth year
## stop_year      - constant:
##
## Function Requirements________________________________________________________
## nfam_step
##
## Package Requirements_________________________________________________________
## NONE
##
ped_step = function(onset_hazard, death_hazard, part,
                    risk, founder_byears,
                    birth_range = c(18, 45),
                    NB_params = c(2, 4/7),
                    stop_year = 2015){

  #initialize a data frame to store all the necessary info for the ped file
  fam_ped = create_pedFile()

  #randomly generate birth year and gender for the founder carrying the RV,
  #and fill in all other fields with the appropriate info
  fam_ped[1,] = c(NA,               #family ID
                  1,                #RV status
                  round(runif(1)),  #gender
                  NA, NA,           #dad_id and #mom_id
                  NA,               #affected status
                  1, 0,             #alleles 1 and 2,
                  1,                #generation no
                  round(runif(1, min = founder_byears[1], max = founder_byears[2])), #birth year
                  NA, NA,           #onset and death years
                  risk,             #RR of developing disease
                  1, 1)             # do_sim and availablilty

  last_id = 1
  last.gen = 1

  #store the ID of individuals for whom we need to simulate life events
  re.sim = fam_ped$ID[which(fam_ped$do_sim == 1 & fam_ped$Gen == last.gen)]
  while (length(re.sim) > 0) {
    for (k in 1:length(re.sim)) {
      new.kin = nfam_step(found_info = fam_ped[which(fam_ped$ID == re.sim[k]),],
                          stop_year, last_id,
                          onset_hazard, death_hazard, part,
                          birth_range, NB_params, risk)

      #replace individual by their simulated self and add family members when necessary
      fam_ped = rbind(fam_ped[-which(fam_ped$ID==re.sim[k]), ],
                      new.kin[[1]])

      last_id = new.kin[[2]]
    }
    last.gen = max(fam_ped$Gen)
    re.sim = fam_ped$ID[which(fam_ped$do_sim == 1 & fam_ped$Gen == last.gen)]
  }

  return(fam_ped[, c(1:8,10:13,15,9)])
}


##--------------##
##  sim_family  ##
##--------------##
## Create a function that will:
##  a.) simulate a pedigree
##  b.) choose a proband from the available candidates
##  c.) trim the family based on recall_probs
##  d.) check to see if the familiy has the requested number of affecteds
##       if the family does not have the requested number of affecteds
##       after steps a or c, we throw away the family and simulate a new one
##
## Arguments____________________________________________________________________
## onset_hazard   - vector; age specific incidence rates for the disease of interest
## death_hazard   - data.frame; age specific mortality rates
## part           - vector; partition over which to apply onset and death rates
## risk           - constant; the RR for developing disease
## founder_byears - vector; length 2, years to be chosen from uniformly for founder birth year
## onset_yspan    - vector, length = 2, the acertainment period of the study
##                          i.e., years in which proband became affected
## num_affected   - constant; number of affected in simulated family
## family_num     - constant; a family identificaton number
## recall_probs   - vector, length n, probability proband recalls relatives of degree n
## birth_range    - vector; max and min birth ages
## NB_params      - vector; size and probability parameters for NB distribution
## stop_year      - constant:
##
## Function Requirements________________________________________________________
## ped_step
## trim_step
##
## Package Requirements_________________________________________________________
## NONE
##
sim_family = function(onset_hazard, death_hazard, part, risk,
                      founder_byears, onset_yspan,
                      num_affected, family_num,
                      recall_probs,
                      birth_range = c(18, 45),
                      NB_params = c(2, 4/7),
                      stop_year = 2015){

  #generate the family pedigree, sheck to see that the untrimmed pedigree has
  # the appropriate number of affected individuals
  D = 0
  while(D <= 0){
    d = 0
    while( d <= 0 ){
      fam_ped = ped_step(onset_hazard, death_hazard, part, risk,
                         founder_byears, birth_range, NB_params, stop_year)
      fam_ped$FamID = family_num

      if( nrow(fam_ped) == 1 | sum(fam_ped$affected) < num_affected |
          length(fam_ped$ID[which(fam_ped$onset_year %in%
                                  onset_yspan[1]:onset_yspan[2])]) < num_affected ){
        d = 0
      } else { d = 1 }
    }

    #trim the pedigree and check to see that the trimmed pedigree has
    # the appropriate number of affecteds
    if (missing(recall_probs)) {
      trim_ped = trim_step(ped_file = fam_ped, onset_yspan)
    } else {
      trim_ped = trim_step(ped_file = fam_ped, onset_yspan, recall_probs)
    }
    #determine the number of available affected individuals
    avail.affect = trim_ped$ID[which(trim_ped$available == 1 & trim_ped$affected == 1)]
    D = ifelse(length(avail.affect) < num_affected, 0, 1)
  }

  #return original and trimmed pedigrees
  my.return = list(fam_ped, trim_ped)
  names(my.return) = c("full.ped", "trim_ped")
  return(my.return)
}

##------------##
##  Sim_Peds  ##
##------------##
## Create a function that will simulate n pedigrees, gives user a progress bar,
## but no parallization here
##
## Arguments____________________________________________________________________
## npeds          - constant; number of pedigrees to simulate
## onset_hazard   - vector; age specific incidence rates for the disease of interest
## death_hazard   - data.frame; age specific mortality rates
## part           - vector; partition over which to apply onset and death rates
## risk           - constant; the RR for developing disease
## founder_byears - vector; length 2, years to be chosen from uniformly for founder birth year
## onset_yspan    - vector, length = 2, the acertainment period of the study
##                          i.e., years in which proband became affected
## num_affected   - constant; number of affected in simulated family
## recall_probs   - vector, length n, probability proband recalls relatives of degree n
## birth_range    - vector; max and min birth ages
## NB_params      - vector; size and probability parameters for NB distribution
## stop_year      - constant:
##
## Function Requirements________________________________________________________
## sim_family
##
## Package Requirements_________________________________________________________
## NONE
##
Sim_Peds = function(npeds, risk,
                    onset_hazard, death_hazard, part,
                    founder_byears, onset_yspan,
                    num_affected, recall_probs,
                    birth_range = c(18, 45), NB_params = c(2, 4/7),
                    stop_year = 2015){


  ped_files = create_pedFile()
  ped.RRs = rep(risk, npeds)
  my.count = 1
  pb <- txtProgressBar(min = 0, max = length(ped.RRs), style = 3)
  for(m in 1:length(ped.RRs)){
    loop.fam = sim_family(onset_hazard, death_hazard, part, risk = ped.RRs[m],
                          founder_byears, onset_yspan,
                          num_affected, family_num = my.count,
                          recall_probs,
                          birth_range, NB_params, stop_year)[[2]]
    ped_files = rbind(ped_files, loop.fam)
    setTxtProgressBar(pb, m)
    my.count = my.count + 1
    }
  return(ped_files)
}
