library(testthat)
context("sim_ped")
test_that("returns a single ped file dataframe", {
  expect_true(is.data.frame(sim_ped(hazard_rates = new.hazard(seq(0, 100, by = 1),
                                                              AgeSpecific_Hazards),
                                    GRR = 10, prob_causalRV = 1,
                                    FamID = 1, stop_year = 2015,
                                    founder_byears = c(1900, 1905))))
})

test_that("pedigree always contains at least 1 person", {
  expect_true(nrow(sim_ped(hazard_rates = new.hazard(seq(0, 100, by = 1),
                                                     AgeSpecific_Hazards),
                           GRR = 10, prob_causalRV = 1,
                           FamID = 1, stop_year = 2015,
                           founder_byears = c(1900, 1905))) >= 1)
})

test_that("Effects of prop_causalRV = 1 and single_founderEntry = T", {
  exPed <- sim_ped(hazard_rates = new.hazard(seq(0, 100, by = 1),
                                             AgeSpecific_Hazards),
                   GRR = 10, prob_causalRV = 1,
                   FamID = 1, stop_year = 2015,
                   founder_byears = c(1900, 1905),
                   single_founderEntry = T)

  #expect that first founder introduces causal variant
  expect_true(1 %in% exPed[1, c(7,8)])
  #expect that only one founder introduces causal variant
  expect_true(sum(exPed[which(is.na(exPed$dadID)), c(7, 8)]) == 1)
  #expect that no children are heterozygous at the disease locus
  expect_true(!any(exPed$DA1 + exPed$DA2 == 2))
})

test_that("Effects of prop_causalRV != 1 and single_founderEntry = T", {
  exPed <- sim_ped(hazard_rates = new.hazard(seq(0, 100, by = 1),
                                             AgeSpecific_Hazards),
                   GRR = 10, prob_causalRV = 0.1,
                   FamID = 1, stop_year = 2015,
                   founder_byears = c(1900, 1905),
                   single_founderEntry = T)

  #expect that only one founder introduces causal variant
  expect_true(sum(exPed[which(is.na(exPed$dadID)), c(7, 8)]) <= 1)
  #expect that no children are heterozygous at the disease locus
  expect_true(!any(exPed$DA1 + exPed$DA2 == 2))
})


test_that("Effects of prop_causalRV == 1 and single_founderEntry = F", {
  exPed <- sim_ped(hazard_rates = new.hazard(seq(0, 100, by = 1),
                                             AgeSpecific_Hazards),
                   GRR = 10, prob_causalRV = 1,
                   FamID = 1, stop_year = 2015,
                   founder_byears = c(1900, 1905),
                   single_founderEntry = F)

  #expect that every founder introduces causal variant
  expect_true(sum(exPed[which(is.na(exPed$dadID)), c(7, 8)]) == nrow(exPed[which(is.na(exPed$dadID)), ]))
})

test_that("If GRR = 1 and single_founderEntry = T no one should have RV", {
  exPed <- sim_ped(hazard_rates = new.hazard(seq(0, 100, by = 1),
                                             AgeSpecific_Hazards),
                   GRR = 1, prob_causalRV = 1,
                   FamID = 1, stop_year = 2015,
                   founder_byears = c(1900, 1905),
                   single_founderEntry = T)

  #expect that every founder introduces causal variant
  expect_true(!any(exPed[, c(7,8)] == 1))
})
