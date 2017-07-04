context("get_nextEvent")
test_that("If current age > max birth age and disease status = 1, next event is death", {
  expect_equal(colnames(get_nextEvent(current_age = 46,
                                      disease_status = 1,
                                      lambda_birth = 0.05,
                                      hazard_rates = new.hazard(AgeSpecific_Hazards),
                                      birth_range = c(18, 45), RR = 1))
               , "Death")
})
