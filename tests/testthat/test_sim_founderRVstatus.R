context("sim_founderRVstatus")
test_that("If GRR = 1, always returns d = (0, 0) and RR = 1", {
  expect_equal(sim_founderRVstatus(GRR = 1, carrier_prob = 0.02,
                                   RVfounder = TRUE),
               list(c(0, 0), 1))
  expect_equal(sim_founderRVstatus(GRR = 1, carrier_prob = 0.02,
                                   RVfounder = FALSE),
               list(c(0, 0), 1))

  expect_equal(sim_founderRVstatus(GRR = c(1, 1), carrier_prob = 0.02,
                                   RVfounder = TRUE),
               list(c(0, 0), c(1, 1)))

  expect_equal(sim_founderRVstatus(GRR = c(1, 1), carrier_prob = 0.02,
                                   RVfounder = FALSE),
               list(c(0, 0), c(1, 1)))
})

test_that("RVfounder = TRUE settings", {
  first_ex <- sim_founderRVstatus(GRR = 10, carrier_prob = 0.02,
                                  RVfounder = TRUE)
  expect_true(1 %in% first_ex[[1]])
  expect_equal(first_ex[[2]], 10)

  second_ex <- sim_founderRVstatus(GRR = c(10, 1), carrier_prob = 0.02,
                                   RVfounder = TRUE)
  expect_true(1 %in% second_ex[[1]])
  expect_equal(second_ex[[2]], c(10, 1))
})


test_that("RVfounder = FALSE settings", {
  #Note that here we set the carrier prob to 1
  #so that the founder must carry a cRV even when
  #RVfounder is False

  first_ex <- sim_founderRVstatus(GRR = 10, carrier_prob = 1,
                                  RVfounder = FALSE)
  expect_true(0 %in% first_ex[[1]])
  expect_true(1 %in% first_ex[[1]])
  expect_equal(first_ex[[2]], 10)

  second_ex <- sim_founderRVstatus(GRR = c(10, 10), carrier_prob = 1,
                                   RVfounder = FALSE)
  expect_true(0 %in% second_ex[[1]])
  expect_true(1 %in% second_ex[[1]])
  expect_equal(second_ex[[2]], c(10, 10))

})
