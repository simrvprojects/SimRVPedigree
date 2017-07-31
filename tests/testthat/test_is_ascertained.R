context("is_ascertained")

EXPed <- sim_RVped(hazard_rates = new.hazard(AgeSpecific_Hazards),
                   GRR = 35, carrier_prob = 0.002,
                   RVfounder = "first",
                   FamID = 1,
                   num_affected = 2,
                   recall_probs = c(1),
                   founder_byears = c(1900, 1910),
                   ascertain_span = c(1970, 2015))


test_that("is_ascertained always returns TRUE for sim_RVped pedigrees", {
  expect_true(is_ascertained(EXPed[[1]],
              num_affected = 2,
              ascertain_span = c(1970, 2015),
              recall_probs = c(1))[[1]])
})

test_that("is_ascertained always a pedigree of size less than or equal to the original pedigree", {
  expect_lte(nrow(is_ascertained(EXPed[[1]],
                                 num_affected = 2,
                                 ascertain_span = c(1970, 2015))[[2]]),
             nrow(EXPed[[1]]))
})
