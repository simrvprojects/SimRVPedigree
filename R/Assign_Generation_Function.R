##--------------------##
##  assignGeneration  ##
##--------------------##
## Create a function that will assign generation number to affecteds
## Generation 1 represents the first affected or first suspected carrier status
##
## Arguments____________________________________________________________________
## ped_file       - data.frame; ped_file generated using ped_step
## 
## Function Requirements________________________________________________________
## NONE
##
## Package Requirements_________________________________________________________
## kinship2
##
assignGeneration = function(ped_file){
  
  #pull affecteds from the pedigree
  ped_afds = ped_file[which(ped_file$affected == 1), ]
  
  if (nrow(ped_afds) == 0) {
    return(ped_afds)
  } else {
    d = 0
    while (d <= 0) {
      #find the dad IDs that are required but have been removed
      readd_dad = ped_afds$dad_id[!is.element(ped_afds$dad_id,
                                              ped_afds$ID[which(ped_afds$gender == 0)])]
      
      #find the mom IDs that are required but have been removed
      readd_mom = ped_afds$mom_id[!is.element(ped_afds$mom_id, 
                                              ped_afds$ID[which(ped_afds$gender == 1)])]
      
      #Now pull the rows containing the required parents from the original ped_file
      readd = ped_file[which(ped_file$ID %in% 
                               c(readd_dad[!is.na(readd_dad)], 
                                 readd_mom[!is.na(readd_mom)])),]
      #and remove all of their simulated data, and mark unavailable
      if (nrow(readd) >= 1) {
        readd$birth_year = NA
        readd$onset_year = NA
        readd$death_year = NA
        readd$available = 0
      }
      
      ped_afds = rbind(ped_afds, readd)
      
      #see if we are still missing any moms or dads
      readd_dad = ped_afds$dad_id[!is.element(ped_afds$dad_id,
                                              ped_afds$ID[which(ped_afds$gender == 0)])]
      readd_mom = ped_afds$mom_id[!is.element(ped_afds$mom_id, 
                                              ped_afds$ID[which(ped_afds$gender == 1)])]
      
      
      if (sum(!is.na(readd_dad)) > 0 | sum(!is.na(readd_mom)) > 0) { 
        d = 0 
      } else {d = 1}
    }
    
    
    ped_afds$Gen = ifelse((ped_afds$affected == 1 &  ped_afds$available == 1),
                          ped_afds$Gen, NA)
    
    
    Gen.tab = table(ped_afds$Gen)
    min.gen = as.numeric(names(Gen.tab[1]))
    min.gen2 = as.numeric(names(Gen.tab[2]))
    gen.diff = min.gen2-min.gen
    num.in.min.gen = as.numeric(Gen.tab[1])
    numcols.in.Gen.tab = length(names(Gen.tab))
    
    
    kin.mat = kinship(ped_afds, 
                      id = ped_afds$ID, 
                      dadid = ped_afds$dad_id, 
                      momid = ped_afds$mom_id)
    
    kin.distance = -log(kin.mat)/log(2)
    
    if (min.gen == 1 | (min.gen == 2 & num.in.min.gen >= 2)) {
      return(ped_afds)
    } else if (num.in.min.gen >= 2) {
      #find the distance between only those in the lowest generation
      kin.distance = -log(kin.mat[which(ped_afds$Gen == min.gen), which(ped_afds$Gen == min.gen)])/log(2)
      #find the gen difference
      new.gen.diff = min.gen - (max(kin.distance)/2 + 1)
      #assign new generation
      ped_afds$Gen[!is.na(ped_afds$Gen)] = ped_afds$Gen[!is.na(ped_afds$Gen)] - new.gen.diff
      return(ped_afds)
    } else if (num.in.min.gen == 1) {
      #find the distance between only those in the lowest 2 generations
      kin.distance = -log(kin.mat[which(ped_afds$Gen == min.gen), which(ped_afds$Gen %in% c(min.gen,min.gen2))])/log(2)
      #find the difference between the maximum distance in kin.distance and 
      # the value of the second smallest generation - this is the value we will use to adjust everyone
      # at or below the second smallest generation
      new.gen.diff = min.gen2 - max(kin.distance)
      #find the difference between the maximum distance in kin.distance 
      # and the difference between the smallest 2 generations 
      # This will be the new gen no for the 1 individual in the lowest generation
      new.gen.oldest = max(kin.distance) - gen.diff
      
      #assign new generation
      ped_afds$Gen[!is.na(ped_afds$Gen)] = ifelse(ped_afds$Gen[!is.na(ped_afds$Gen)] == min.gen,
                                                  new.gen.oldest,
                                                  ped_afds$Gen[!is.na(ped_afds$Gen)] - new.gen.diff)
      
      return(ped_afds)
    } 
  }
}

##--------------##
##  ReGen_Peds  ##
##--------------##
## Create a function that will assign generation number to affecteds
## for a simulated dataset
##
## Arguments____________________________________________________________________
## ped_file       - data.frame; ped_file generated using ped_step
## 
## Function Requirements________________________________________________________
## assignGeneration
##
## Package Requirements_________________________________________________________
## kinship2
##
ReGen_Peds = function(Ped_Sample) {
  fam.nums = unique(Ped_Sample$FamID)
  ReGen = assignGeneration(ped_file = Ped_Sample[which(Ped_Sample$FamID == fam.nums[1]), ])
  for(k in 2:length(fam.nums)) {
    ReGen = rbind(ReGen,
                  assignGeneration(ped_file = Ped_Sample[which(Ped_Sample$FamID==fam.nums[k]), ]))
  }
  return(ReGen)
}


##----------------------##
##  est_NumAffect_Dist  ##
##----------------------##
## Create a function that will estimate the distribution of the number of affecteds
## for a simulated dataset of families
##
## Arguments____________________________________________________________________
## Ped_Sample     - data.frame; sample dataset of ped_files generated using Sim_Peds
## 
## Function Requirements________________________________________________________
## NONE
##
## Package Requirements_________________________________________________________
## dplyr
##
library(dplyr)
est_NADist = function(Ped_Sample){
  NumAff_data = Ped_Sample %>% group_by(FamID) %>% summarize(NumAffected = sum(affected))
  
  fam.nums = unique(Ped_Sample$FamID)
  Number_Affected_tab = table(NumAff_data$NumAffected)/length(fam.nums)
  
  ecdf_NumAffected = ecdf(NumAff_data$NumAffected)
  my.return = list(NumAff_data, ecdf_NumAffected)
  return(my.return)
}


##---------------------##
##  get_Affected_info  ##
##---------------------##
## Create a function that will get info such as DOB, DX age and reassigned generation
##  from affecteds
##
## Arguments____________________________________________________________________
## Ped_Sample     - data.frame; sample dataset of ped_files generated using Sim_Peds
## 
## Function Requirements________________________________________________________
## ReGen_Peds
##
## Package Requirements_________________________________________________________
## dplyr
##
get_Affected_info = function(Ped_Sample){
  ReGen = ReGen_Peds(Ped_Sample)
  keep.members = which(ReGen$affected==1 & ReGen$available==1)
  AffectedData = data.frame(Affected_Gen = ReGen$Gen[keep.members],
                            Onset_Year = ReGen$onset_year[keep.members],
                            Birth_Year = ReGen$birth_year[keep.members],
                            Onset_Age = ReGen$onset_year[keep.members]-ReGen$birth_year[keep.members],
                            Fam_NO = ReGen$FamID[keep.members],
                            RR = ReGen$RR[keep.members])
  
  Ctest = cor.test(AffectedData$Onset_Age, AffectedData$Birth_Year, method = "kendall")
  my.return = list(AffectedData, Ctest$estimate, Ctest$conf.int, Ctest$p.value)
  return(my.return)
}

##--------------##
##  get_AveIBD  ##
##--------------##
## Create a function that will return average IBD for familiy
##
## Arguments____________________________________________________________________
## ped_file     - data.frame; single ped_file
## 
## Function Requirements________________________________________________________
## NONE
##
## Package Requirements_________________________________________________________
## kinship2
##
get_AveIBD = function(ped_file) {
  
  if (0 %in% unique(ped_file$gender)) {
    tped = with(ped_file, pedigree(id = ID, dad = dad_id, mom = mom_id, sex = gender+1))
  } else {
    tped = with(ped_file, pedigree(id = ID, dad = dad_id, mom = mom_id, sex = gender))
  }
  
  tkin = kinship(tped)[which(ped_file$affected ==1 & ped_file$available == 1),
                       which(ped_file$affected ==1 & ped_file$available == 1)]
  AveIBD = mean(2*tkin[lower.tri(tkin, diag = F)])
  #changed to 2* to get IBD
  return(AveIBD)
}

get_kinDist = function(ped_file) {
  
  if (0 %in% unique(ped_file$gender)) {
    tped = with(ped_file, pedigree(id = ID, dad = dad_id, mom = mom_id, sex = gender+1))
  } else {
    tped = with(ped_file, pedigree(id = ID, dad = dad_id, mom = mom_id, sex = gender))
  }
  
  tkin = log(kinship(tped)[which(ped_file$affected ==1 & ped_file$available == 1),
                       which(ped_file$affected ==1 & ped_file$available == 1)])/log(2)
  AveKdist = mean(tkin[lower.tri(tkin, diag = F)])
  #changed to 2* to get IBD
  return(AveKdist)
}


##----------------##
##  get_FamStats  ##
##----------------##
## Create a function that will compute family statistics for a simulated dataset of families
##
## Arguments____________________________________________________________________
## Ped_Sample     - data.frame; sample dataset of ped_files generated using Sim_Peds
## 
## Function Requirements________________________________________________________
## get_Affected_info
##
## Package Requirements_________________________________________________________
## NONE
##
get_FamStats = function(Ped_Sample){
  #reassign generation based on affected status
  ReGen_Sample = ReGen_Peds(Ped_Sample) 
  
  FamAff_data = ReGen_Sample %>% group_by(FamID) %>% summarize(NumAffected = sum(affected),
                                                               NumGens = max(Gen[which(affected == 1 & available == 1)]),
                                                               AveAoo = mean(onset_year[which(affected == 1 & available == 1)] - birth_year[which(affected == 1 & available == 1)]))
  FamAff_data$AveIBD = NA
  for (i in 1:nrow(FamAff_data)) {
    FamAff_data$AveIBD[i] = get_AveIBD(ped_file = ReGen_Sample[which(ReGen_Sample$FamID == FamAff_data$FamID[i]), ])
  }
  
  return(FamAff_data)
}
