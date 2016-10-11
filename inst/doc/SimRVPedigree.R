## ------------------------------------------------------------------------
library(SimRVPedigree)

## ------------------------------------------------------------------------
#load example hazards 
data("AgeSpecific_Hazards")
colnames(AgeSpecific_Hazards)

#the first column of the AgeSpecific_Hazards dataset provides the age-specific
#onset hazard for the population
my_onset_hazard = AgeSpecific_Hazards[,1]

#We must specify a partition of ages over which to apply the age-specifc
#hazards.  Note that a valid partition must contain 1 element more that the
#age-specific hazards. Assuming that the age-specific hazards are specified 
#in 1 year increments starting at birth we can specify the partition of ages 
#as follows. 
age_part = seq(0, length(my_onset_hazard), by = 1)  

#Let's assume that the individual has inherited a rare variant which has an
#associated relative risk of disease onset of 10, and that the indiviual is
#currently 45 years old.  To simulate the waiting time to onset we would type
#the following command:
set.seed(1)
Time_to_onset <- get_WaitTime(p = runif(1), last_event = 45, 
                 hazard = my_onset_hazard*10,
                 part = age_part)
Time_to_onset


#Note that in some instances, no waiting time to onset is simulated  
Time_to_onset <- get_WaitTime(p = 0.99, last_event = 45, 
                 hazard = my_onset_hazard*10,
                 part = age_part)
Time_to_onset

## ------------------------------------------------------------------------
#the second and third columns of the AgeSpecific_Hazards dataset, respectively,
#provide the age-specific death hazards for the unaffected and affected 
#populations.  Note: if the disease of interest is sufficiently rare the 
#population death hazard may be used as an estimate for the unaffected death 
#hazard.
my_death_hazard <- AgeSpecific_Hazards[, c(2,3)]

#Note the effect of the different death hazards: 
#First we simulate the waiting time to death for a 45 year old who has NOT 
#experienced onset.
set.seed(123)
Time_to_death_unaffected <- get_WaitTime(p = runif(1), last_event = 45, 
                                         hazard = my_death_hazard[,1],
                                         part = age_part, scale = TRUE)
Time_to_death_unaffected

#Next, we simulate the waiting time to death for a 45 year old who HAS 
#experienced disease onset.
set.seed(123)
Time_to_death_affected <- get_WaitTime(p = runif(1), last_event = 45, 
                                       hazard = my_death_hazard[,2],
                                       part = age_part, scale = TRUE)
Time_to_death_affected

## ------------------------------------------------------------------------
# First, we illustrate how all life events are simulated for an individual who 
# has inherited the rare variant 
set.seed(2)
Life_Events <- get_lifeEvents(RV_status = 1, 
                              onset_hazard = my_onset_hazard,
                              death_hazard = my_death_hazard,
                              part = age_part,
                              birth_range = c(18, 45),
                              NB_params = c(2, 4/7),
                              RR = 10, YOB = 1900)

Life_Events

# Note the effect of keeping all factors constant except for rare variant #status.  
set.seed(2)
Life_Events <- get_lifeEvents(RV_status = 0, 
                              onset_hazard = my_onset_hazard,
                              death_hazard = my_death_hazard,
                              part = age_part,
                              birth_range = c(18, 45),
                              NB_params = c(2, 4/7),
                              RR = 10, YOB = 1900)

Life_Events

## ------------------------------------------------------------------------
 #Simulate a random pedigree
 set.seed(6)
 ex_ped <- sim_ped(onset_hazard = my_onset_hazard,
                   death_hazard = my_death_hazard,
                   part = age_part, stop_year = 2015,
                   RR = 5, FamID = 1,
                   founder_byears = c(1900, 1910))

## ---- results = FALSE, tidy = TRUE, message = FALSE----------------------
# Pedigrees may be plotted usign the kinship2 package.  Assuming that 
# this package has been installed it is loaded by executing the command:
 library(kinship2)

## ---- fig.height=4, fig.width=6------------------------------------------
# Define pedigree to use kinship2's plot function
 ex_pedigree <- pedigree(id = ex_ped$ID,
                         dadid = ex_ped$dad_id,
                         momid = ex_ped$mom_id,
                         sex = (ex_ped$gender + 1),
                         affected = cbind(Affected = ex_ped$affected,
                                          RV_status = ex_ped$DA1 + 
                                                      ex_ped$DA2),
                         famid = ex_ped$FamID)['1']
 plot(ex_pedigree, cex = 0.75)
 pedigree.legend(ex_pedigree, location = "topleft",  radius = 0.25)

## ---- fig.height = 4.55, fig.width=7, echo = FALSE-----------------------
 #Read in example pedigree to trim
 data(exp_peds)

 #plot first pedigree in example_ped using kinship2
 library(kinship2)
 ex_pedigree <- pedigree(id = exp_peds$ID,
                         dadid = exp_peds$dad_id,
                         momid = exp_peds$mom_id,
                         sex = (exp_peds$gender + 1),
                         affected = cbind(Affected = exp_peds$affected,
                                          RV_status = exp_peds$DA1 +
                                                      exp_peds$DA2),
                         famid = exp_peds$FamID)['1']
 plot(ex_pedigree, cex = 0.85)
 pedigree.legend(ex_pedigree, location = "topleft",  radius = 0.25)
 mtext("Original Pedigree", side = 3, line = 2)

## ---- echo = FALSE, results = FALSE, eval = TRUE, fig.height=4.5, fig.width=7----
set.seed(2)
#trim pedigree
TrimPed <- trim_pedigree(ped_file = exp_peds[which(exp_peds$FamID == 1), ],
                             ascertain_span = c(2005, 2015),
                             num_affected = 2,
                             recall_probs = c(1))

#plot trimmed pedigree
Tped <- pedigree(id = TrimPed$ID,
                 dadid = TrimPed$dad_id,
                 momid = TrimPed$mom_id,
                 sex = (TrimPed$gender + 1),
                 affected = cbind(Affected = TrimPed$affected,
                                  Proband = TrimPed$is_proband,
                                  RV_status = TrimPed$DA1 + TrimPed$DA2),
                 famid = TrimPed$FamID)['1']

plot(Tped, cex = 0.85)
pedigree.legend(Tped, location = "topleft",  radius = 0.25)
mtext("recall_probs = c(1), \n ascertain_span = (2005, 2015), \n num_affected = 2",
      side = 3 )

## ---- echo = FALSE, results = FALSE, eval = TRUE, fig.height=4.55, fig.width=7----
set.seed(2)
#trim pedigree
TrimPed <- trim_pedigree(ped_file = exp_peds[which(exp_peds$FamID == 1), ],
                             ascertain_span = c(2005, 2015),
                             num_affected = 2,
                             recall_probs = c(1, 0.5, 0.25, 0.125))

#plot trimmed pedigree
Tped <- pedigree(id = TrimPed$ID,
                 dadid = TrimPed$dad_id,
                 momid = TrimPed$mom_id,
                 sex = (TrimPed$gender + 1),
                 affected = cbind(Affected = TrimPed$affected,
                                  Proband = TrimPed$is_proband,
                                  RV_status = TrimPed$DA1 + TrimPed$DA2),
                 famid = TrimPed$FamID)['1']

plot(Tped, cex = 0.85)
pedigree.legend(Tped, location = "topleft",  radius = 0.25)
mtext("recall_probs = c(1, 0.5, 0.25, 0.125), \n ascertain_span = (2005, 2015), \n num_affected = 2",
                 side = 3 )

## ------------------------------------------------------------------------
#Load hazard data
data(AgeSpecific_Hazards)

#specify onset hazard, death hazard, and age parition over which to apply 
#hazards.
my_onset_hazard <- AgeSpecific_Hazards[,1]
my_death_hazard <- AgeSpecific_Hazards[,c(2,3)]
age_part <- seq(0, 100, by = 1)

 #Simulate a random pedigree
 set.seed(6)
 ex_RVped <- sim_RVpedigree(onset_hazard = my_onset_hazard,
                            death_hazard = my_death_hazard,
                            part = age_part,
                            RR = 5, FamID = 1, stop_year = 2015,
                            founder_byears = c(1900, 1910),
                            ascertain_span = c(1900, 2015),
                            num_affected = 2,
                            recall_probs = c(1, 0.5, 0.25, 0.125))

## ---- results = FALSE, tidy = TRUE, echo = FALSE, message=FALSE----------
#Pedigrees may be plotted usign the kinship2 package.  Assuming that this 
#package has been installed it is loaded by executing the command:
 library(kinship2)

## ---- fig.height=5, fig.width=7, echo = FALSE----------------------------
 #Define pedigree to use kinship2's plot function
 ex_RVpedigree_FULL <- pedigree(id = ex_RVped[[1]]$ID,
                           dadid = ex_RVped[[1]]$dad_id,
                           momid = ex_RVped[[1]]$mom_id,
                           sex = (ex_RVped[[1]]$gender + 1),
                           affected = cbind(Affected = ex_RVped[[1]]$affected,
                                            RV_status = ex_RVped[[1]]$DA1 + 
                                                        ex_RVped[[1]]$DA2),
                         famid = ex_RVped[[1]]$FamID)['1']
 plot(ex_RVpedigree_FULL)
 pedigree.legend(ex_RVpedigree_FULL, location = "topleft",  radius = 0.25)
 mtext("Pedigree Prior to Trimming and Proband Selection", side = 3)

 #Define pedigree to use kinship2's plot function
 ex_RVpedigree_TRIM <- pedigree(id = ex_RVped[[2]]$ID,
                           dadid = ex_RVped[[2]]$dad_id,
                           momid = ex_RVped[[2]]$mom_id,
                           sex = (ex_RVped[[2]]$gender + 1),
                           affected = cbind(Affected = ex_RVped[[2]]$affected,
                                            Proband = ex_RVped[[2]]$is_proband,
                                            RV_status = ex_RVped[[2]]$DA1 + 
                                                        ex_RVped[[2]]$DA2),
                         famid = ex_RVped[[2]]$FamID)['1']
 plot(ex_RVpedigree_TRIM)
 pedigree.legend(ex_RVpedigree_TRIM, location = "bottomleft",  radius = 0.25)
 mtext("Ascertained Pedigree", side = 3)

## ---- eval = FALSE-------------------------------------------------------
#  #SITATION FOR DOPARALLEL VIGNETTE?
#  #Note that: while parallel processsing may acheived using only the doParallel
#  #package, to ensure that simulations are reproducible we must also incorporate
#  #the doRNG package.
#  
#  #assuming they have been installed the required packages are loaded using the
#  #commands:
#  library(doParallel)
#  library(doRNG)
#  
#  # Before we create our cluster, let's determine how many processors are
#  #currently in use using the getDoParWorkers function.  Since we have not created
#  #a cluster yet, this function should return 1.
#  getDoParWorkers()
#  
#  #The number of cores available for parallel processing will depend on the the
#  #computer  being used.  To determine how many cores are available for parallel
#  #processing on your computer, simply execute the following command:
#  detectCores()
#  
#  #To run simulations in parallel we must create a cluster and then register the
#  #cluster. The following code illustrates how to create and register a cluster
#  #that will simulate pedigrees in parallel on 2 cores.
#  cl <- makeCluster(2)       # create cluster
#  registerDoParallel(cl)     # register cluster
#  
#  #Note that getDoParWorkers() should now return 2 instead of 1
#  getDoParWorkers()
#  
#  #To avoid problems, after you are finished using the cluster, you will want to
#  #stop it. This can be acheived by executing:
#  on.exit(stopCluster(cl))
#  #which will stop the cluster when you end the R session, or by executing
#  #stopCluster(cl) after the simulation is complete.
#  
#  #to ensure reproducibility, we make use of the %dorng% operator provided by the
#  #doRNG package, in the foreach loop, by specifing a seed after .option.RNG.
#  npeds <- 8    #set the number of pedigrees to generate
#  
#  RV_peds = foreach(i = seq(npeds),
#                    .combine = rbind,
#                    .packages = c("kinship2", "SimRVPedigree"),
#                    .options.RNG = 1984
#                    ) %dorng% {
#                      sim_RVpedigree(onset_hazard = my_onset_hazard,
#                                     death_hazard = my_death_hazard,
#                                     part = age_part, RR = 1.5,
#                                     FamID = i, stop_year = 2015,
#                                     founder_byears = c(1900, 1980),
#                                     ascertain_span = c(2000, 2015),
#                                     recall_probs = c(1),
#                                     num_affected = 2)[[2]]
#                      }

## ---- echo = FALSE, fig.height=5, fig.width=7----------------------------
#Read in example pedigree to trim
data(exp_peds)

library(kinship2)
#assign to pedigree object to show before and after behavior of
#the assign_affectedGen function
ex_pedigree <- pedigree(id = exp_peds$ID,
                        dadid = exp_peds$dad_id,
                        momid = exp_peds$mom_id,
                        sex = (exp_peds$gender + 1),
                        affected = exp_peds$affected,
                        famid = exp_peds$FamID)


#create df to store peds with reassigned generation number
RAG_peds <- exp_peds[1,]
RAG_peds <- RAG_peds[-1,]

for(i in 1:4){
  RAG_peds <- rbind(RAG_peds,
                    assign_affectedGen(exp_peds[which(exp_peds$FamID == i), ]))
}

RAG_pedigrees <-  pedigree(id = RAG_peds$ID,
                           dadid = RAG_peds$dad_id,
                           momid = RAG_peds$mom_id,
                           sex = (RAG_peds$gender + 1),
                           affected = RAG_peds$affected,
                           famid = RAG_peds$FamID)

par(mfrow = c(1, 2))
k = 1
  ID1 = paste0("ID", sep = ":",
               exp_peds[which(exp_peds$FamID == k), 2],
               sep = "\n Gen:", exp_peds[which(exp_peds$FamID == k), 14])
  ID2 = paste0("ID", sep = ":",
               RAG_peds[which(RAG_peds$FamID == k), 2],
               sep = "\n Gen:", RAG_peds[which(RAG_peds$FamID == k), 14])
  plot(ex_pedigree[paste0(k)], id = ID1, cex = 0.5)
  mtext("Pedigree before generation reassignment", side = 3, line = 1)
  plot(RAG_pedigrees[paste0(k)], id = ID2, cex = 0.5)
  mtext("Pedigree after generation reassignment", side = 3, line = 1)

## ---- echo = FALSE, fig.height=5, fig.width=7----------------------------
par(mfrow = c(1, 2))
k = 3
  ID1 = paste0("ID", sep = ":",
               exp_peds[which(exp_peds$FamID == k), 2],
               sep = "\n Gen:", exp_peds[which(exp_peds$FamID == k), 14])
  ID2 = paste0("ID", sep = ":",
               RAG_peds[which(RAG_peds$FamID == k), 2],
               sep = "\n Gen:", RAG_peds[which(RAG_peds$FamID == k), 14])
  plot(ex_pedigree[paste0(k)], id = ID1, cex = 0.5)
  mtext("Pedigree before generation reassignment", side = 3, line = 1)
  plot(RAG_pedigrees[paste0(k)], id = ID2, cex = 0.5)
  mtext("Pedigree after generation reassignment", side = 3, line = 1)

