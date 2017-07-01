library(testthat)
context("sim_RVped")

EXPed <- sim_RVped(hazard_rates = AgeSpecific_Hazards,
                   part = seq(0, 100, by = 1),
                   GRR = 35, prob_causalRV = 1,
                   FamID = 1,
                   num_affected = 2,
                   recall_probs = c(1),
                   founder_byears = c(1900, 1910),
                   ascertain_span = c(1970, 2015))


test_that("returns a list containing two pedfiles and two numeric constants", {
  expect_true(is.list(EXPed))
  expect_true(is.data.frame(EXPed[[1]]))
  expect_true(is.data.frame(EXPed[[2]]))
  expect_true(is.numeric(EXPed[[3]]))
  expect_true(is.numeric(EXPed[[4]]))
})

test_that("pedigree always conatains more than 1 person", {
  expect_gt(nrow(EXPed[[2]]), 1)
})


test_that("both pedigrees contains at least 2 affecteds when num_affected = 2", {
  RVped1 <- EXPed[[1]]
  RVped2 <- EXPed[[2]]
  expect_gte(sum(RVped1$affected[which(RVped1$available == 1)]), 2)
  expect_gte(sum(RVped2$affected[which(RVped2$available == 1)]), 2)
})


test_that("proband in trimmed pedigree had 1 affected relative before onset, when num_affected = 2", {
  RVped <- sim_RVped(hazard_rates = AgeSpecific_Hazards,
                     part = seq(0, 100, by = 1),
                     GRR = 35, prob_causalRV = 1,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1),
                     founder_byears = c(1900, 1910),
                     ascertain_span = c(1970, 2015))[[2]]

  Oyears <- RVped$onsetYr[which(RVped$affected == 1 &
                                     RVped$available == 1 &
                                     RVped$proband == 0)]

  expect_gte(length(Oyears[which(Oyears <= 2015 & Oyears >= 1970)]), 1)
})

test_that("issues errors when invalid partition supplied", {
  expect_error(sim_RVped(sim_RVped(hazard_rates = AgeSpecific_Hazards,
                                   part = seq(1, 100, by = 1),
                                   GRR = 35, prob_causalRV = 1, FamID = 1,
                                   num_affected = 2,
                                   recall_probs = c(1),
                                   founder_byears = c(1900, 1910),
                                   ascertain_span = c(1970, 2015))))
  })


test_that("issues error when hazard_rates contains fewer than 3 columns", {
  expect_error(sim_RVped(hazard_rates = AgeSpecific_Hazards[, c(1,2)],
                         part = seq(0, 100, by = 1),
                         GRR = 35, prob_causalRV = 1, FamID = 1,
                         num_affected = 2,
                         founder_byears = c(1900, 1980),
                         ascertain_span = c(2000, 2015)))
  })

test_that("issues error when RR < 0", {
  expect_error(sim_RVped(hazard_rates = AgeSpecific_Hazards,
                         part = seq(0, 100, by = 1),
                         GRR = -1, prob_causalRV = 1, FamID = 1,
                         num_affected = 2,
                         founder_byears = c(1900, 1980),
                         ascertain_span = c(2000, 2015)))
})

test_that("issues error when hazard contains NA values", {
  haz_rates = rbind(AgeSpecific_Hazards[c(1:99),], c(NA, 1, 1))
  expect_error(sim_RVped(hazard_rates = haz_rates,
                         part = seq(0, 100, by = 1),
                         GRR = 35, prob_causalRV = 1, FamID = 1,
                         num_affected = 2,
                         founder_byears = c(1900, 1980),
                         ascertain_span = c(2000, 2015)))
  })

test_that("issues error when part contains NA values", {
  expect_error(sim_RVped(hazard_rates = AgeSpecific_Hazards,
                         part = c(NA, seq(1, 100, by = 1)),
                         GRR = 35, prob_causalRV = 1, FamID = 1,
                         num_affected = 2,
                         founder_byears = c(1900, 1980),
                         ascertain_span = c(2000, 2015)))
  })

test_that("issues error when ascertain_span not properly specified", {
    expect_error(sim_RVped(hazard_rates = AgeSpecific_Hazards,
                           part = seq(0, 100, by = 1),
                           GRR = 35, prob_causalRV = 1, FamID = 1,
                           num_affected = 2,
                           founder_byears = c(1900, 1980),
                           ascertain_span = c(2017, 2015)))
  })

test_that("issues error when part doesn't start at zero", {
  expect_error(sim_RVped(hazard_rates = AgeSpecific_Hazards,
                         part = seq(50, 100, by = 0.5),
                         GRR = 35, prob_causalRV = 1, FamID = 1,
                         num_affected = 2,
                         founder_byears = c(1900, 1980),
                         ascertain_span = c(2000, 2015)))
  })

test_that("issues error when birth_range not properly specified", {
  expect_error(sim_RVped(hazard_rates = AgeSpecific_Hazards,
                         part = seq(0, 100, by = 1),
                         GRR = 35, prob_causalRV = 1, FamID = 1,
                         num_affected = 2,
                         birth_range = c(10, 5),
                         founder_byears = c(1900, 1980),
                         ascertain_span = c(2000, 2015)))
  })

test_that("issues error when recall_probs not properly specified", {
  expect_error(sim_RVped(hazard_rates = AgeSpecific_Hazards,
                         part = seq(0, 100, by = 1),
                         GRR = 35, prob_causalRV = 1, FamID = 1,
                         num_affected = 2,
                         recall_probs = c(10, 5),
                         founder_byears = c(1900, 1980),
                         ascertain_span = c(2000, 2015)))
  })
