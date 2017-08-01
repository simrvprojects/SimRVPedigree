context("censor_ped")
test_that("censor_ped returns an error when no proband or censor year provided", {
  RVped <- sim_RVped(hazard_rates = new.hazard(AgeSpecific_Hazards),
                     GRR = 50, carrier_prob = 0.002,
                     RVfounder = TRUE,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1, 1, 1, 0.5, 0.25),
                     founder_byears = c(1900, 1980),
                     ascertain_span = c(1980, 2015))[[1]][, c(1, 14)]

  expect_error(censor_ped(RVped))
})

test_that("censor_ped returns a smaller or equally sized pedfile", {
  RVped <- sim_RVped(hazard_rates = new.hazard(AgeSpecific_Hazards),
                     GRR = 50, carrier_prob = 0.002,
                     RVfounder = TRUE,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1, 1, 1, 0.5, 0.25),
                     founder_byears = c(1900, 1980),
                     ascertain_span = c(1980, 2015))[[2]]

  expect_gte(nrow(RVped), nrow(censor_ped(RVped)))
})

test_that("censor_ped does not return any info after the censor year", {
  RVped <- sim_RVped(hazard_rates = new.hazard(AgeSpecific_Hazards),
                     GRR = 50, carrier_prob = 0.002,
                     RVfounder = TRUE,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1, 1, 1, 0.5, 0.25),
                     founder_byears = c(1900, 1980),
                     ascertain_span = c(1980, 2015))[[2]]

  my_Cyear <- RVped$onsetYr[which(RVped$proband == 1)]
  C_ped <- censor_ped(ped_file = RVped, censor_year = my_Cyear)
  expect_gte(my_Cyear, max(C_ped$birthYr, na.rm = T))
  expect_gte(my_Cyear, max(C_ped$onsetYr, na.rm = T))
  if(!all(is.na(C_ped$deathYr))){
    expect_gte(my_Cyear, max(C_ped$deathYr, na.rm = T))
  }
})
