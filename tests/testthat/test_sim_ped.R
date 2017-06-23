library(testthat)
context("sim_ped")
test_that("returns a single ped file dataframe", {
  expect_true(is.data.frame(sim_ped(hazard_rates = AgeSpecific_Hazards,
                                    part = seq(0, 100, by = 1),
                                    GRR = 10, FamID = 1, stop_year = 2015,
                                    founder_byears = c(1900, 1980))))
})

test_that("pedigree always contains at least 1 person", {
  expect_true(nrow(sim_ped(hazard_rates = AgeSpecific_Hazards,
                           part = seq(0, 100, by = 1),
                           GRR = 10, FamID = 1, stop_year = 2015,
                           founder_byears = c(1900, 1980))) >= 1)
})

