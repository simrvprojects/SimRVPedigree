library(testthat)
context("choose_proband")
test_that("returns a dataframe with the proband variable with a single proband", {
  RVped <- sim_RVped(onset_hazard = AgeSpecific_Hazards[,1],
                     death_hazard = AgeSpecific_Hazards[,c(2,3)],
                     part = seq(0, 100, by = 1),
                     RR = 50, FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1),
                     founder_byears = c(1900, 1980),
                     ascertain_span = c(1970, 2015))[[1]][, c(1:14)]

  Tped = choose_proband(ped = RVped, num_affected = 2,
                        ascertain_span = c(1970, 2015))

  expect_true("proband" %in% colnames(Tped))
  expect_equal(sum(Tped$proband), 1)
})

test_that("The proband is an affected and experienced onset during the ascertainment span", {
  RVped <- sim_RVped(onset_hazard = AgeSpecific_Hazards[,1],
                    death_hazard = AgeSpecific_Hazards[,c(2,3)],
                    part = seq(0, 100, by = 1),
                    RR = 50, FamID = 1,
                    num_affected = 2,
                    recall_probs = c(1),
                    founder_byears = c(1900, 1970),
                    ascertain_span = c(1980, 2015))[[1]][, c(1:14)]

  Tped = choose_proband(ped = RVped, num_affected = 2,
                        ascertain_span = c(1970, 2015))


  expect_true(Tped$affected[which(Tped$proband == 1)] == 1)
  expect_true(Tped$onset_year[which(Tped$proband == 1)] >= 1970 &
                Tped$onset_year[which(Tped$proband == 1)] <= 2015)

})


