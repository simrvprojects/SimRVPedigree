context("reassign_gen")
test_that("reassign_gen returns a smaller or equally sized pedfile", {
  RVped <- sim_RVped(hazard_rates = hazard(AgeSpecific_Hazards),
                     GRR = 50, carrier_prob = 0.002,
                     RVfounder = TRUE,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1),
                     founder_byears = c(1900, 1905),
                     ascertain_span = c(1980, 2015))[[2]]

  expect_gte(nrow(RVped), nrow(reassign_gen(RVped)))
})


test_that("maximum re-assigned gen is at most maximum gen from original pedigree ", {
  RVped <- sim_RVped(hazard_rates = hazard(AgeSpecific_Hazards),
                     GRR = 50, carrier_prob = 0.002,
                     RVfounder = TRUE,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1),
                     founder_byears = c(1900, 1905),
                     ascertain_span = c(1980, 2015))[[2]]

  expect_gte(max(RVped$Gen), max(reassign_gen(RVped)$Gen, na.rm = TRUE))
})


test_that("never have two affected individuals with reassigned gen = 1", {
  RVped <- sim_RVped(hazard_rates = hazard(AgeSpecific_Hazards),
                     GRR = 50, carrier_prob = 0.002, RVfounder = TRUE,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1),
                     founder_byears = c(1900, 1905),
                     ascertain_span = c(1980, 2015))[[2]]

  GenTab <- table(reassign_gen(RVped)$Gen)
  GenTab
  if ("1" %in% names(GenTab)) {
    expect_equal(as.numeric(GenTab[which(names(GenTab) == 1)]), 1)
  }
})