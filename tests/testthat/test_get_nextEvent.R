context("get_nextEvent")
test_that("If current age > max birth age and disease status = 1, next event is death", {
  expect_equal(colnames(get_nextEvent(current_age = 46,
                                      disease_status = 1,
                                      lambda_birth = 0.05,
                                      onset_hazard = AgeSpecific_Hazards[,1],
                                      death_hazard = AgeSpecific_Hazards[,c(2:3)],
                                      part = seq(0, 100, by = 1),
                                      birth_range = c(18, 45), RR = 1))
               , "Death")
})
