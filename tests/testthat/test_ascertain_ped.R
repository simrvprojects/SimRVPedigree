library(SimRVPedigree)
context("ascertain_ped")

data("AgeSpecific_Hazards")
data("SubtypeHazards")

EXPed <- sim_RVped(hazard_rates = hazard(AgeSpecific_Hazards),
                   GRR = 35, carrier_prob = 0.002,
                   RVfounder = TRUE,
                   FamID = 1,
                   num_affected = 2,
                   recall_probs = c(1),
                   founder_byears = c(1900, 1910),
                   ascertain_span = c(1970, 2015))


test_that("ascertain_ped always returns TRUE for sim_RVped pedigrees", {
  #one subtype
  EXPed <- sim_RVped(hazard_rates = hazard(AgeSpecific_Hazards),
                     GRR = 35, carrier_prob = 0.002,
                     RVfounder = TRUE,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1),
                     founder_byears = c(1900, 1910),
                     ascertain_span = c(1970, 2015))

  expect_true(ascertain_ped(EXPed[[1]], num_affected = 2,
                            ascertain_span = c(1970, 2015),
                            recall_probs = c(1))[[1]])

  #two subtypes
  EXPed <- sim_RVped(hazard_rates = hazard(SubtypeHazards, subtype_ID = c("HL", "NHL")),
                     GRR = c(1, 50), carrier_prob = 0.002,
                     RVfounder = TRUE,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1),
                     founder_byears = c(1900, 1910),
                     ascertain_span = c(1970, 2015))

  expect_true(ascertain_ped(EXPed[[1]], num_affected = 2,
                            ascertain_span = c(1970, 2015),
                            recall_probs = c(1))[[1]])

})

test_that("ascertain_ped always a pedigree of size less than or equal to the original pedigree", {
  EXPed <- sim_RVped(hazard_rates = hazard(AgeSpecific_Hazards),
                     GRR = 35, carrier_prob = 0.002,
                     RVfounder = TRUE,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1),
                     founder_byears = c(1900, 1910),
                     ascertain_span = c(1970, 2015))

  expect_lte(nrow(ascertain_ped(EXPed[[1]],
                                num_affected = 2,
                                ascertain_span = c(1970, 2015))[[2]]),
             nrow(EXPed[[1]]))

  #two_subtypes
  EXPed <- sim_RVped(hazard_rates = hazard(SubtypeHazards, subtype_ID = c("HL", "NHL")),
                     GRR = c(50, 50), carrier_prob = 0.002,
                     RVfounder = TRUE,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1),
                     founder_byears = c(1900, 1910),
                     ascertain_span = c(1970, 2015))

  expect_lte(nrow(ascertain_ped(EXPed[[1]],
                                num_affected = 2,
                                ascertain_span = c(1970, 2015))[[2]]),
             nrow(EXPed[[1]]))

})


test_that("ascertained pedigrees meet number affected criteria", {

  re_sim <- TRUE
  while (re_sim){
    exPed <- sim_ped(hazard_rates = hazard(AgeSpecific_Hazards),
                     GRR = 50, carrier_prob = 0.002,
                     RVfounder = TRUE, FamID = 1,
                     founder_byears = c(1900, 1910))
    my_n <- sample(x = c(2, 3), size = 1)
    as <- sort(round(runif(2, 1970, 2010)))
    as[1] <- as[1] - 5

    aped <- ascertain_ped(exPed,
                          num_affected = my_n,
                          ascertain_span = as,
                          recall_probs = c(1, 1, 0.5))


    #collect data on available affecteds
    aff_dat <- aped[[2]][aped[[2]]$affected & aped[[2]]$available, ]

    #if pedigree is not ascertained, keep trying
    re_sim <- !aped[[1]]

  }


  #pedigree contains at least n affected relatives
  expect_true(sum(aped[[2]]$affected[aped[[2]]$available == 1]) >= my_n)

  # at least one affected experienced onset during the ascertainment span
  expect_true(sum(aff_dat$onsetYr >= as[1] & aff_dat$onsetYr <= as[2]) > 0)

  # the proband experienced onset during the ascertainment span
  expect_true(aff_dat$onsetYr[aff_dat$proband == 1] >= as[1] & aff_dat$onsetYr[aff_dat$proband == 1] <= as[2])

})


test_that("pedigrees with the correct number of affected fail if onset at wrong time", {

  re_sim <- TRUE
  while (re_sim){
    exPed <- sim_ped(hazard_rates = hazard(AgeSpecific_Hazards),
                     GRR = 50, carrier_prob = 0.002,
                     RVfounder = TRUE, FamID = 1,
                     founder_byears = c(1900, 1910))

    #number of affected required for ascertainment
    my_n <- sample(x = c(2, 3), size = 1)

    #ascertainment span
    as <- sort(round(runif(2, 2005, 2015)))
    as[1] <- as[1] - 5

    #check to see if pedigree is ascertained
    aped <- ascertain_ped(exPed, num_affected = my_n, ascertain_span = as, recall_probs = c(1, 1, 0.5))

    #collect data on affected relatives
    aff_dat <- aped[[2]][aped[[2]]$affected & aped[[2]]$available, ]

    #check to see if the pedigree was disqualfied, but had the correct
    #number of affecteds
    if (nrow(aff_dat) >= my_n & !aped[[1]]) {
      re_sim <- FALSE
    }

  }

  expect_true(sum(aff_dat$onsetYr >= as[1] & aff_dat$onsetYr <= as[2]) == 0)
})




test_that("pedigrees with the incorrect number of affecteds fail", {

  re_sim <- TRUE
  while (re_sim){
    exPed <- sim_ped(hazard_rates = hazard(AgeSpecific_Hazards),
                     GRR = 50, carrier_prob = 0.002,
                     RVfounder = TRUE, FamID = 1,
                     founder_byears = c(1900, 1910))
    my_n <- sample(x = c(2, 3), size = 1)
    as <- sort(round(runif(2, 1970, 2010)))
    as[1] <- as[1] - 5


    aped <- ascertain_ped(exPed, num_affected = my_n, ascertain_span = as, recall_probs = c(1))
    my_ped <- aped[[2]][aped[[2]]$available == 1, ]

    re_sim <- aped[[1]]

  }

  expect_true(sum(my_ped$affected) < my_n)
})


test_that("pedigrees have the correct number of indivdiuals with a particular subtype", {

  re_sim <- FALSE
  while (!re_sim){
    exPed <- sim_ped(hazard_rates = hazard(SubtypeHazards, subtype_ID = c("HL", "NHL")),
                     GRR = c(100, 100), carrier_prob = 0.002,
                     RVfounder = TRUE, FamID = 1,
                     founder_byears = c(1900, 1910),
                     stop_year = 2019)

    my_n <- sample(x = c(2, 3), size = 1)
    my_n2 <- my_n - 1
    as <- c(1900, 2019)


    aped <- ascertain_ped(exPed, num_affected = my_n,
                          ascertain_span = as, recall_probs = c(1),
                          sub_criteria = list("HL", my_n2))

    my_ped <- aped[[2]][aped[[2]]$available == 1, ]

    re_sim <- aped[[1]]

  }

  expect_true(sum(my_ped$affected) >= my_n & sum(my_ped$affected[which(my_ped$subtype == "HL")]) >= my_n2)
})
