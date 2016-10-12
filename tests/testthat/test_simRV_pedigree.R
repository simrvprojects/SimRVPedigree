library(testthat)
context("sim_RVpedigree")
test_that("returns a list containing two pedfiles", {
  expect_true(is.list(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                                           death_hazard = AgeSpecific_Hazards[,c(2,3)],
                                           part = seq(0, 100, by = 1),
                                           RR = 20, FamID = 1,
                                           num_affected = 2,
                                           founder_byears = c(1900, 1980),
                                           ascertain_span = c(2000, 2015))))

  expect_true(is.data.frame(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                                           death_hazard = AgeSpecific_Hazards[,c(2,3)],
                                           part = seq(0, 100, by = 1),
                                           RR = 20, FamID = 1,
                                           num_affected = 2,
                                           founder_byears = c(1900, 1980),
                                           ascertain_span = c(2000, 2015))[[1]]))

  expect_true(is.data.frame(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                                           death_hazard = AgeSpecific_Hazards[,c(2,3)],
                                           part = seq(0, 100, by = 1),
                                           RR = 20, FamID = 1,
                                           num_affected = 2,
                                           founder_byears = c(1900, 1980),
                                           ascertain_span = c(2000, 2015))[[2]]))
  })

test_that("pedigree always conatains more than 1 person", {
  expect_gt(nrow(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                                  death_hazard = AgeSpecific_Hazards[,c(2,3)],
                                  part = seq(0, 100, by = 1),
                                  RR = 20, FamID = 1,
                                  num_affected = 2,
                                  founder_byears = c(1900, 1980),
                                  ascertain_span = c(2000, 2015))[[2]]), 1)
})


test_that("pedigree always conatains more than 2 affecteds when num_affected = 2", {
  RVped <- sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                          death_hazard = AgeSpecific_Hazards[,c(2,3)],
                          part = seq(0, 100, by = 1),
                          RR = 20, FamID = 1,
                          num_affected = 2,
                          founder_byears = c(1900, 1980),
                          ascertain_span = c(2000, 2015))
  RVped1 <- RVped[[1]]
  RVped2 <- RVped[[2]]
  expect_gte(sum(RVped1$affected[which(RVped1$available == 1)]), 2)
  expect_gte(sum(RVped2$affected[which(RVped2$available == 1)]), 2)
})


test_that("proband in trimmed pedigree had 1 affected relative before onset, when num_affected = 2", {
  RVped <- sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                          death_hazard = AgeSpecific_Hazards[,c(2,3)],
                          part = seq(0, 100, by = 1),
                          RR = 20, FamID = 1,
                          num_affected = 2,
                          founder_byears = c(1900, 1980),
                          ascertain_span = c(2000, 2015))[[2]]

  Oyears <- RVped$onset_year[which(RVped$affected == 1 &
                                     RVped$available == 1 &
                                     RVped$is_proband == 0)]

  expect_gte(length(Oyears[which(Oyears <= 2015 & Oyears >= 2000)]), 1)
})

test_that("issues errors when invalid partition supplied", {
  expect_error(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                              death_hazard = AgeSpecific_Hazards[,c(2,3)],
                              part = seq(1, 100, by = 1),
                              RR = 20, FamID = 1,
                              num_affected = 2,
                              founder_byears = c(1900, 1980),
                              ascertain_span = c(2000, 2015)))
  })


test_that("issues error when death_hazard contains only 1 column", {
  expect_error(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                              death_hazard = AgeSpecific_Hazards[,c(2)],
                              part = seq(0, 100, by = 1),
                              RR = 20, FamID = 1,
                              num_affected = 2,
                              founder_byears = c(1900, 1980),
                              ascertain_span = c(2000, 2015)))
  })

test_that("issues error when hazard contains NA values", {
  expect_error(sim_RVpedigree(onset_hazard = c(AgeSpecific_Hazards[1:80,1], NA,
                                               AgeSpecific_Hazards[82:100,1]),
                              death_hazard = AgeSpecific_Hazards[,c(2,3)],
                              part = seq(0, 100, by = 1),
                              RR = 20, FamID = 1,
                              num_affected = 2,
                              founder_byears = c(1900, 1980),
                              ascertain_span = c(2000, 2015)))
  })

test_that("issues error when part contains NA values", {
  expect_error(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                              death_hazard = AgeSpecific_Hazards[,c(2,3)],
                              part = c(NA, seq(1, 100, by = 1)),
                              RR = 20, FamID = 1,
                              num_affected = 2,
                              founder_byears = c(1900, 1980),
                              ascertain_span = c(2000, 2015)))
  })

test_that("issues error when ascertain_span not properly specified", {
    expect_error(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                              death_hazard = AgeSpecific_Hazards[,c(2,3)],
                              part = seq(0, 100, by = 1),
                              RR = 20, FamID = 1,
                              num_affected = 2,
                              founder_byears = c(1900, 1980),
                              ascertain_span = c(2017, 2015)))
  })

test_that("issues error when part doesn't start at zero", {
  expect_error(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                                death_hazard = AgeSpecific_Hazards[,c(2,3)],
                                part = seq(50, 100, by = 0.5),
                                RR = 20, FamID = 1,
                                num_affected = 2,
                                founder_byears = c(1900, 1980),
                                ascertain_span = c(2000, 2015)))
  })

test_that("issues error when birth_range not properly specified", {
  expect_error(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                              death_hazard = AgeSpecific_Hazards[,c(2,3)],
                              part = seq(0, 100, by = 1),
                              RR = 20, FamID = 1,
                              num_affected = 2,
                              birth_range = c(10, 5),
                              founder_byears = c(1900, 1980),
                              ascertain_span = c(2000, 2015)))
  })

test_that("issues error when recall_probs not properly specified", {
  expect_error(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                              death_hazard = AgeSpecific_Hazards[,c(2,3)],
                              part = seq(0, 100, by = 1),
                              RR = 20, FamID = 1,
                              num_affected = 2,
                              recall_probs = c(10, 5),
                              founder_byears = c(1900, 1980),
                              ascertain_span = c(2000, 2015)))
  })

test_that("issues warnings when partition doesn't span appropriate ages", {
  expect_warning(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                                death_hazard = AgeSpecific_Hazards[,c(2,3)],
                                part = seq(0, 50, by = 0.5),
                                RR = 20, FamID = 1,
                                num_affected = 2,
                                founder_byears = c(1900, 1980),
                                ascertain_span = c(2000, 2015)))
})

test_that("issues warnings when suspected that affected death hazard listed first", {
  expect_warning(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                                death_hazard = AgeSpecific_Hazards[,c(3,2)],
                                part = seq(0, 100, by = 1),
                                RR = 20, FamID = 1,
                                num_affected = 2,
                                founder_byears = c(1900, 1980),
                                ascertain_span = c(2000, 2015)))
})
