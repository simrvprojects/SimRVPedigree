library(testthat)
context("assign_affectedGen")
test_that("assign_affectedGen returns a smaller or equally sized pedfile", {
  RVped <- sim_RVped(onset_hazard = AgeSpecific_Hazards[,1],
                     death_hazard = AgeSpecific_Hazards[,c(2,3)],
                     part = seq(0, 100, by = 1),
                     RR = 50, FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1, 1, 1, 0.5, 0.25),
                     founder_byears = c(1900, 1980),
                     ascertain_span = c(1980, 2015))[[2]]

  expect_gte(nrow(RVped), nrow(assign_affectedGen(RVped)))
})


test_that("maximum re-assigned gen is at most maximum gen from original pedigree ", {
  RVped <- sim_RVped(onset_hazard = AgeSpecific_Hazards[,1],
                     death_hazard = AgeSpecific_Hazards[,c(2,3)],
                     part = seq(0, 100, by = 1),
                     RR = 50, FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1, 1, 1, 0.5, 0.25),
                     founder_byears = c(1900, 1980),
                     ascertain_span = c(1980, 2015))[[2]]

  expect_gte(max(RVped$Gen), max(assign_affectedGen(RVped)$Gen, na.rm = TRUE))
})


test_that("never have two affected individuals with reassigned gen = 1", {
  RVped <- sim_RVped(onset_hazard = AgeSpecific_Hazards[,1],
                     death_hazard = AgeSpecific_Hazards[,c(2,3)],
                     part = seq(0, 100, by = 1),
                     RR = 50, FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1, 1, 1, 0.5, 0.25),
                     founder_byears = c(1900, 1980),
                     ascertain_span = c(1980, 2015))[[2]]

  GenTab <- table(assign_affectedGen(RVped)$Gen)
  GenTab
  if ("1" %in% names(GenTab)) {
    expect_equal(as.numeric(GenTab[which(names(GenTab) == 1)]), 1)
  }
})
