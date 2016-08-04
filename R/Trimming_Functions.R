##-------------##
##  trim_step  ##
##-------------##
## Create a function that will:
##  a.) choose a proband
##  b.) determine which of the proband's relatives will remain in pedigree
##     based on recall.prob.  Here the entries of recall.prob are the probabilies
##     of the proband recalling relatives of various degrees.
## More specifically, let k(p, r_i) be the kinship coefficient between the proband
## and relative i.  Now, let p_k denote the probability that the proband can recall a
## relative whose kinship coefficent is k(p, r_i).
##
## Arguments____________________________________________________________________
## ped_file       - data.frame; ped_file generated using ped_step
## onset_yspan    - vector, length = 2, the acertainment period of the study
##                          i.e., years in which proband became affected
## recall_probs   - vector, length n, probability proband recalls relatives of degree n
##
## Function Requirements________________________________________________________
## nfam_step
##
## Package Requirements_________________________________________________________
## kinship2
##
trim_step = function(ped_file, onset_yspan, recall_probs){
  #First we randomly choose an individual who experienced
  # onset during the time span defined by onset_yspan

  AffIDs = ped_file[which(ped_file$affected==1),which(colnames(ped_file)%in%c("onset_year", "ID"))]
  AffIDs$poss_proband = ifelse(AffIDs$onset_year%in% onset_yspan[1]:onset_yspan[2], 1, 0)


  if (sum(AffIDs$poss_proband) == 1) {
    #only 1 available proband
    probandID = AffIDs$ID[which(AffIDs$poss_proband == 1)]
  } else if ( sum(abs(AffIDs$poss_proband - 1)) > 0 ) {
    #multiple available probands but 2 affected condition already met by start of ascertainment period
    probandID = sample(size = 1, x = AffIDs$ID[which(AffIDs$poss_proband == 1)])
  } else{
    #no affecteds before ascertainment period
    AffIDs$poss_proband[which.min(AffIDs$onset_year)] = 0
    probandID =  sample(size = 1, x = AffIDs$ID[which(AffIDs$poss_proband == 1)])
  }

  #generate a vector of Unif(0,1) RVs, these will be used to
  # determine who is trimmed and who is kept
  u = runif(length(ped_file$ID))

  #calculate the kinship matrix for our pedigree
  kin.mat = kinship(ped_file,
                    id = ped_file$ID,
                    dadid = ped_file$dad_id,
                    momid = ped_file$mom_id)

  kin.proband = kin.mat[, which(ped_file$ID == probandID)]

  if (missing(recall_probs)) {
    #Now, we keep only those individuals for whom the kinship coefficent with the proband*4
    # is greater than runif(1).  This will result in keeping parents, offspring, and siblings
    # with probability 1
    ped.trim = ped_file[kin.proband/0.25 >= u, ]
  } else {
    #This next chuck (where rprobs is created) associates the recall probabilities
    # specified by the user with the appropriate individuals in the pedigree.
    rprobs = rep(NA, length(ped_file$ID))
    #set recall probability for proband to 1
    rprobs[which(kin.proband == 0.5)] = 1
    #set recall probability to non-related individuals to 0
    # NOTE: those who need to be added back in will be later
    rprobs[which(kin.proband == 0)] = 0

    #create a vector of kinship coefficient with the same length as the recall.prob
    # vector specified by the user
    kin.list = 2^{-seq(from = 2, to = (length(recall_probs)+1), by = 1)}

    for (i in 1:(length(kin.list)-1)) {
      rprobs[which(kin.proband == kin.list[i])] = recall_probs[i]
    }
    rprobs[is.na(rprobs)] = recall_probs[length(recall_probs)]

    ped.trim = ped_file[rprobs >= u, ]
  }

  d = 0
  while (d <= 0) {
    #find the dad IDs that are required but have been removed
    readd.dad = ped.trim$dad_id[!is.element(ped.trim$dad_id,
                                            ped.trim$ID[which(ped.trim$gender == 0)])]

    #find the mom IDs that are required but have been removed
    readd.mom = ped.trim$mom_id[!is.element(ped.trim$mom_id,
                                            ped.trim$ID[which(ped.trim$gender == 1)])]

    #Now pull the rows containing the required parents from the original ped_file
    readd = ped_file[which(ped_file$ID %in%
                             c(readd.dad[!is.na(readd.dad)],
                               readd.mom[!is.na(readd.mom)])), ]
    #and remove all of their simulated data, and mark unavailable
    if (nrow(readd) >= 1) {
      readd$birth_year = NA
      readd$onset_year = NA
      readd$death_year = NA
      readd$available  = 0
    }

    ped.trim = rbind(ped.trim, readd)

    #see if we are still missing any moms or dads
    readd.dad = ped.trim$dad_id[!is.element(ped.trim$dad_id,
                                            ped.trim$ID[which(ped.trim$gender == 0)])]
    readd.mom = ped.trim$mom_id[!is.element(ped.trim$mom_id,
                                            ped.trim$ID[which(ped.trim$gender == 1)])]


    if (sum(!is.na(readd.dad)) > 0 | sum(!is.na(readd.mom)) > 0) {
      d = 0
    } else {d = 1}
  }

  ped.trim$is_proband = ifelse(ped.trim$ID == probandID, 1, 0)
  return(ped.trim)
}
