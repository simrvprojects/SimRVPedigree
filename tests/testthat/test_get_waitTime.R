context("get_waitTime")
test_that("wait_time return is in [0, max(part) - last_event]", {
  expect_gte(get_WaitTime(p = 0, last_event = 0,
                          hazard = AgeSpecific_Hazards[,2],
                          part = seq(0, 100, by = 1),
                          scale = TRUE), 0)
  expect_lte(get_WaitTime(p = 0, last_event = 0,
                          hazard = AgeSpecific_Hazards[,2],
                          part = seq(0, 100, by = 1),
                          scale = TRUE),
             max(seq(0, 100, by = 1)))
  expect_gte(get_WaitTime(p = 0, last_event = 0,
                          hazard = AgeSpecific_Hazards[,2],
                          part = seq(0, 100, by = 1)), 0)
  expect_lte(get_WaitTime(p = 0, last_event = 0,
                          hazard = AgeSpecific_Hazards[,2],
                          part = seq(0, 100, by = 1)),
             max(seq(0, 100, by = 1)))
  expect_gte(get_WaitTime(p = 1, last_event = 0,
                          hazard = AgeSpecific_Hazards[,2],
                          part = seq(0, 100, by = 1),
                          scale = TRUE), 0)
  expect_lte(get_WaitTime(p = 1, last_event = 0,
                          hazard = AgeSpecific_Hazards[,2],
                          part = seq(0, 100, by = 1),
                          scale = TRUE),
             max(seq(0, 100, by = 1)))
})

test_that("when scale = TRUE, and p = 1, wait_time != NA", {
  expect_equal(get_WaitTime(p = 1, last_event = 0,
                            hazard = AgeSpecific_Hazards[,1],
                            part = seq(0, 100, by = 1),
                            scale = TRUE), 100)
  expect_equal(get_WaitTime(p = 1, last_event = 0,
                            hazard = AgeSpecific_Hazards[,2],
                            part = seq(0, 100, by = 1),
                            scale = TRUE), 100)
  expect_equal(get_WaitTime(p = 1, last_event = 0,
                            hazard = AgeSpecific_Hazards[,3],
                            part = seq(0, 100, by = 1),
                            scale = TRUE), 100)
})

test_that("when scale = FALSE, and p = 1, wait_time == NA", {
  expect_equal(get_WaitTime(p = 1, last_event = 0,
                            hazard = AgeSpecific_Hazards[,1],
                            part = seq(0, 100, by = 1)), NA)
})

