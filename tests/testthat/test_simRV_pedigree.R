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

test_that("issues errors when check functions fail", {
  expect_error(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                              death_hazard = AgeSpecific_Hazards[,c(2,3)],
                              part = seq(1, 100, by = 1),
                              RR = 20, FamID = 1,
                              num_affected = 2,
                              founder_byears = c(1900, 1980),
                              ascertain_span = c(2000, 2015)))

  expect_error(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                              death_hazard = AgeSpecific_Hazards[,c(2,3)],
                              part = seq(0, 100, by = 1),
                              RR = 20, FamID = 1,
                              num_affected = 2,
                              founder_byears = c(1980, 1980),
                              ascertain_span = c(2000, 2015)))

  expect_error(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                              death_hazard = AgeSpecific_Hazards[,c(2,3)],
                              part = seq(0, 100, by = 1),
                              RR = 20, FamID = 1,
                              num_affected = 2,
                              founder_byears = c(1900, 1980),
                              ascertain_span = c(2017, 2015)))

  expect_error(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                              death_hazard = AgeSpecific_Hazards[,c(2,3)],
                              part = seq(0, 100, by = 1),
                              RR = 20, FamID = 1,
                              num_affected = 2,
                              founder_byears = c(1900, 1980),
                              ascertain_span = c(2017, 2015)))

  expect_error(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                                death_hazard = AgeSpecific_Hazards[,c(2,3)],
                                part = seq(50, 100, by = 0.5),
                                RR = 20, FamID = 1,
                                num_affected = 2,
                                founder_byears = c(1900, 1980),
                                ascertain_span = c(2000, 2015)))
  })

test_that("issues warning for partitions that don't span appropriate ages", {
  expect_warning(sim_RVpedigree(onset_hazard = AgeSpecific_Hazards[,1],
                                death_hazard = AgeSpecific_Hazards[,c(2,3)],
                                part = seq(0, 50, by = 0.5),
                                RR = 20, FamID = 1,
                                num_affected = 2,
                                founder_byears = c(1900, 1980),
                                ascertain_span = c(2000, 2015)))
})
