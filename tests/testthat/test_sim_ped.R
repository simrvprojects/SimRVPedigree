library(testthat)
context("sim_ped")
test_that("returns a single ped file dataframe", {
  expect_true(is.data.frame(sim_ped(onset_hazard = AgeSpecific_Hazards[,1],
                                    death_hazard = AgeSpecific_Hazards[,c(2,3)],
                                    part = seq(0, 100, by = 1),
                                    RR = 10, FamID = 1, stop_year = 2015,
                                    founder_byears = c(1900, 1980))))
})

test_that("pedigree always contains at least 1 person", {
  expect_true(nrow(sim_ped(onset_hazard = AgeSpecific_Hazards[,1],
                           death_hazard = AgeSpecific_Hazards[,c(2,3)],
                           part = seq(0, 100, by = 1),
                           RR = 10, FamID = 1, stop_year = 2015,
                           founder_byears = c(1900, 1980))) >= 1)
})

