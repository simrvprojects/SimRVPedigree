context("sim_RVped")

EXPed <- sim_RVped(hazard_rates = new.hazard(AgeSpecific_Hazards),
                   GRR = 35, carrier_prob = 0.002,
                   RVfounder = "first",
                   FamID = 1,
                   num_affected = 2,
                   recall_probs = c(1),
                   founder_byears = c(1900, 1910),
                   ascertain_span = c(1970, 2015))


test_that("returns a list containing two pedfiles", {
  expect_true(is.list(EXPed))
  expect_true(is.data.frame(EXPed[[1]]))
  expect_true(is.data.frame(EXPed[[2]]))
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
  RVped <- sim_RVped(hazard_rates = new.hazard(AgeSpecific_Hazards),
                     GRR = 35, carrier_prob = 0.002,
                     RVfounder = "first",
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


test_that("issues error when RR < 0", {
  expect_error(sim_RVped(hazard_rates = new.hazard(AgeSpecific_Hazards),
                         GRR = -1, carrier_prob = 0.002, FamID = 1,
                         num_affected = 2,
                         founder_byears = c(1900, 1980),
                         ascertain_span = c(2000, 2015)))
})

test_that("issues error when ascertain_span not properly specified", {
    expect_error(sim_RVped(hazard_rates = new.hazard(AgeSpecific_Hazards),
                           GRR = 35, carrier_prob = 0.002, FamID = 1,
                           num_affected = 2,
                           founder_byears = c(1900, 1980),
                           ascertain_span = c(2017, 2015)))
  })

test_that("issues error when birth_range not properly specified", {
  expect_error(sim_RVped(hazard_rates = new.hazard(AgeSpecific_Hazards),
                         GRR = 35, carrier_prob = 0.002, FamID = 1,
                         num_affected = 2,
                         birth_range = c(10, 5),
                         founder_byears = c(1900, 1980),
                         ascertain_span = c(2000, 2015)))
  })

test_that("issues error when recall_probs not properly specified", {
  expect_error(sim_RVped(hazard_rates = new.hazard(AgeSpecific_Hazards),
                         GRR = 35, carrier_prob = 0.002, FamID = 1,
                         num_affected = 2,
                         recall_probs = c(10, 5),
                         founder_byears = c(1900, 1980),
                         ascertain_span = c(2000, 2015)))
  })

test_that("issues error when hazard_rates is not a hazard object", {
  expect_error(sim_RVped(hazard_rates = list(AgeSpecific_Hazards, seq(0, 100, by = 1)),
                         GRR = 35, carrier_prob = 0.002, FamID = 1,
                         num_affected = 2,
                         founder_byears = c(1900, 1980),
                         ascertain_span = c(2000, 2015)))
})
