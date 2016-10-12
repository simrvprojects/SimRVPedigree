context("get_lifeEvents")
test_that("get_lifeEvents should always begin at start and end at death", {
  Levents <- get_lifeEvents(RV_status = 0,
                            onset_hazard = AgeSpecific_Hazards[,1],
                            death_hazard = AgeSpecific_Hazards[,c(2:3)],
                            part = seq(0, 100, by = 1),
                            birth_range = c(18, 45), NB_params = c(2, 4/7),
                            RR = 1, YOB = 1900)

  expect_equal(names(Levents)[1], "Start")
  expect_equal(names(Levents)[length(Levents)], "Death")
})

test_that("get_lifeEvents never returns onset more than once", {
  Levents <- get_lifeEvents(RV_status = 1,
                            onset_hazard = AgeSpecific_Hazards[,1],
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

test_that("get_lifeEvents always returns death event after all other events", {
  Levents <- get_lifeEvents(RV_status = 0,
                            onset_hazard = AgeSpecific_Hazards[,1],
                            death_hazard = AgeSpecific_Hazards[,c(2:3)],
                            part = seq(0, 100, by = 1),
                            birth_range = c(18, 45), NB_params = c(2, 4/7),
                            RR = 50, YOB = 1900)
  Levents <- as.numeric(Levents)
  num_events <- length(Levents)
  expect_equal(sum(Levents[num_events] >= Levents[-num_events])/(num_events - 1),
               1)
})
