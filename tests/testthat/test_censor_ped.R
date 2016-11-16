library(testthat)
context("censor_ped")
test_that("censor_ped returns an error when no proabnd or censor year provided", {
  RVped <- sim_RVped(onset_hazard = AgeSpecific_Hazards[,1],
                     death_hazard = AgeSpecific_Hazards[,c(2,3)],
                     part = seq(0, 100, by = 1),
                     RR = 50, FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1, 1, 1, 0.5, 0.25),
                     founder_byears = c(1900, 1980),
                     ascertain_span = c(1980, 2015))[[1]]

  expect_error(censor_ped(RVped))
})

test_that("censor_ped returns a smaller or equally sized pedfile", {
  RVped <- sim_RVped(onset_hazard = AgeSpecific_Hazards[,1],
                     death_hazard = AgeSpecific_Hazards[,c(2,3)],
                     part = seq(0, 100, by = 1),
                     RR = 50, FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1, 1, 1, 0.5, 0.25),
                     founder_byears = c(1900, 1980),
                     ascertain_span = c(1980, 2015))[[2]]

  expect_gte(nrow(RVped), nrow(censor_ped(RVped)))
})

test_that("censor_ped does not return any info after the censor year", {
  RVped <- sim_RVped(onset_hazard = AgeSpecific_Hazards[,1],
                     death_hazard = AgeSpecific_Hazards[,c(2,3)],
                     part = seq(0, 100, by = 1),
                     RR = 50, FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1, 1, 1, 0.5, 0.25),
                     founder_byears = c(1900, 1980),
                     ascertain_span = c(1980, 2015))[[2]]

  my_Cyear <- RVped$onset_year[which(RVped$proband == 1)]
  C_ped <- censor_ped(ped_file = RVped, censor_year = my_Cyear)
  expect_gte(my_Cyear, max(C_ped$birth_year, na.rm = T))
  expect_gte(my_Cyear, max(C_ped$onset_year, na.rm = T))
  expect_gte(my_Cyear, max(C_ped$death_year, na.rm = T))
})
