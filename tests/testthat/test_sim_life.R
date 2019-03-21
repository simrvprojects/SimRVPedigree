context("sim_life")
test_that("sim_life should always begin at start and end at death when stop_year is sufficiently large", {
  Levents <- sim_life(hazard_rates = hazard(AgeSpecific_Hazards),
                      GRR = 25, carrier_prob = 0.02, RV_status = T,
                      YOB = 1900, stop_year = 2001,
                      NB_params = c(2, 4/7))

  expect_equal(names(Levents[[1]])[1], "Start")
  expect_equal(names(Levents[[1]])[length(Levents[[1]])], "Death")

  expect_true(!is.na(Levents$death_event))
})

test_that("sim_life never returns more than one onset event", {
  onset_occured = FALSE

  while (onset_occured == FALSE) {
    Levents <- sim_life(hazard_rates = hazard(AgeSpecific_Hazards),
                        GRR = 500, carrier_prob = 0.02, RV_status = T,
                        YOB = 1900, stop_year = 2001,
                        NB_params = c(2, 4/7))

    if (!is.na(Levents$onset_event)) {
      onset_occured = TRUE
    }
  }

  expect_equal(as.numeric(table(names(Levents[[1]]))[names(table(names(Levents[[1]]))) == "Onset"]), 1)

})

test_that("sim_life always returns death event after all other events", {
  Levents <- sim_life(hazard_rates = hazard(AgeSpecific_Hazards),
                      GRR = 50, carrier_prob = 0.02, RV_status = T,
                      YOB = 1900, stop_year = 2001,
                      NB_params = c(2, 4/7))

  num_events <- length(Levents[[1]])
  expect_equal(sum(Levents[[1]][num_events] >= Levents[[1]][-num_events])/(num_events - 1),
               1)

})

test_that("recorded events always agree", {
  all_events_occur = FALSE

  while (all_events_occur == FALSE) {
    Levents <- sim_life(hazard_rates = hazard(AgeSpecific_Hazards),
                        GRR = 500, carrier_prob = 0.02, RV_status = T,
                        YOB = 1900, stop_year = 2001,
                        NB_params = c(2, 4/7))

    if (!is.na(Levents$onset_event) & !is.null(Levents$repro_event)) {
      all_events_occur = TRUE
    }
  }

  num_events <- length(Levents[[1]])

  #check that death event agrees
  expect_equal(as.numeric(Levents[[1]][num_events]),
               Levents$death_event)

  #check that reproduction events agree
  expect_equal(as.numeric(Levents[[1]][names(Levents[[1]]) == "Child"]),
               Levents$repro_events)

  #check that onset events agree
  expect_equal(as.numeric(Levents[[1]][names(Levents[[1]]) == "Onset"]),
               Levents$onset_event)

  #check that subtype is "no_subtypes"
  expect_equal("no_subtypes",
               Levents$subtype)

})

test_that("sim_life doesn't return any events after the stop year", {
  my_stopY <- 1900 + round(runif(1, min = 30, max = 60))
  Levents <- sim_life(hazard_rates = hazard(AgeSpecific_Hazards),
                      GRR = 50, carrier_prob = 0.02, RV_status = T,
                      YOB = 1900, stop_year = my_stopY,
                      NB_params = c(2, 4/7))

  expect_gte(my_stopY, as.numeric(Levents[[1]][length(Levents[[1]])]))
})


#---------------------------------#
# repeat for 2-subtype simulation #
#---------------------------------#
test_that("sim_life should always begin at start and end at death when stop_year is sufficiently large", {
  Levents <- sim_life(hazard_rates = hazard(SubtypeHazards,
                                            subtype_ID = c("HL", "NHL")),
                      GRR = c(25, 25), carrier_prob = 0.02, RV_status = T,
                      YOB = 1900, stop_year = 2001,
                      NB_params = c(2, 4/7))

  expect_equal(names(Levents[[1]])[1], "Start")
  expect_equal(names(Levents[[1]])[length(Levents[[1]])], "Death")

  expect_true(!is.na(Levents$death_event))
})

test_that("sim_life never returns more than one onset event", {
  onset_occured = FALSE

  while (onset_occured == FALSE) {
    Levents <- sim_life(hazard_rates = hazard(SubtypeHazards,
                                              subtype_ID = c("HL", "NHL")),
                        GRR = c(50, 50), carrier_prob = 0.02, RV_status = T,
                        YOB = 1900, stop_year = 2001,
                        NB_params = c(2, 4/7))

    if (!is.na(Levents$onset_event)) {
      onset_occured = TRUE
    }
  }

  expect_equal(as.numeric(table(names(Levents[[1]]))[names(table(names(Levents[[1]]))) == "Onset"]), 1)

})

test_that("sim_life always returns death event after all other events", {
  Levents <- sim_life(hazard_rates = hazard(SubtypeHazards,
                                            subtype_ID = c("HL", "NHL")),
                      GRR = c(50, 50), carrier_prob = 0.02, RV_status = T,
                      YOB = 1900, stop_year = 2001,
                      NB_params = c(2, 4/7))

  num_events <- length(Levents[[1]])
  expect_equal(sum(Levents[[1]][num_events] >= Levents[[1]][-num_events])/(num_events - 1),
               1)

})

test_that("recorded events always agree", {
  all_events_occur = FALSE

  while (all_events_occur == FALSE) {
    Levents <- sim_life(hazard_rates = hazard(SubtypeHazards,
                                              subtype_ID = c("HL", "NHL")),
                        GRR = c(50, 50), carrier_prob = 0.02, RV_status = T,
                        YOB = 1900, stop_year = 2001,
                        NB_params = c(2, 4/7))

    if (!is.na(Levents$onset_event) & !is.null(Levents$repro_event)) {
      all_events_occur = TRUE
    }
  }

  num_events <- length(Levents[[1]])

  #check that death event agrees
  expect_equal(as.numeric(Levents[[1]][num_events]),
               Levents$death_event)

  #check that reproduction events agree
  expect_equal(as.numeric(Levents[[1]][names(Levents[[1]]) == "Child"]),
               Levents$repro_events)

  #check that reproduction events agree
  expect_equal(as.numeric(Levents[[1]][names(Levents[[1]]) == "Onset"]),
               Levents$onset_event)

  #check that subtype is "no_subtypes"
  expect_true(Levents$subtype %in% c("HL", "NHL"))

})

test_that("sim_life doesn't return any events after the stop year", {
  my_stopY <- 1900 + round(runif(1, min = 30, max = 60))
  Levents <- sim_life(hazard_rates = hazard(SubtypeHazards,
                                            subtype_ID = c("HL", "NHL")),
                      GRR = c(50, 50), carrier_prob = 0.02, RV_status = T,
                      YOB = 1900, stop_year = my_stopY,
                      NB_params = c(2, 4/7))

  expect_gte(my_stopY, as.numeric(Levents[[1]][length(Levents[[1]])]))
})

