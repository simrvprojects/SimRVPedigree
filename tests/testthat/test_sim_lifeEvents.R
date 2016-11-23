context("sim_lifeEvents")
test_that("sim_lifeEvents should always begin at start and end at death", {
  Levents <- sim_lifeEvents(onset_hazard = AgeSpecific_Hazards[,1],
                            death_hazard = AgeSpecific_Hazards[,c(2:3)],
                            part = seq(0, 100, by = 1),
                            birth_range = c(18, 45), NB_params = c(2, 4/7),
                            RR = 25, YOB = 1900)

  expect_equal(names(Levents)[1], "Start")
  expect_equal(names(Levents)[length(Levents)], "Death")
})

test_that("sim_lifeEvents never returns onset more than once", {
  Levents <- sim_lifeEvents(onset_hazard = AgeSpecific_Hazards[,1],
                            death_hazard = AgeSpecific_Hazards[,c(2:3)],
                            part = seq(0, 100, by = 1),
                            birth_range = c(18, 45), NB_params = c(2, 4/7),
                            RR = 50, YOB = 1900)
  if("Onset" %in% names(table(names(Levents)))){
    expect_equal(as.numeric(table(names(Levents))[which(names(table(names(Levents))) ==
                                                          "Onset")]),
                 1)
  }
})

test_that("sim_lifeEvents always returns death event after all other events", {
  Levents <- sim_lifeEvents(onset_hazard = AgeSpecific_Hazards[,1],
                            death_hazard = AgeSpecific_Hazards[,c(2:3)],
                            part = seq(0, 100, by = 1),
                            birth_range = c(18, 45), NB_params = c(2, 4/7),
                            RR = 50, YOB = 1900)
  Levents <- as.numeric(Levents)
  num_events <- length(Levents)
  expect_equal(sum(Levents[num_events] >= Levents[-num_events])/(num_events - 1),
               1)
})
