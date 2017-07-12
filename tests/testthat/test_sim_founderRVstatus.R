context("sim_founderRVstatus")
test_that("If GRR = 1, always returns d = (0, 0) and RR = 1, and leaves intro_RV unchanged", {
  expect_equal(sim_founderRVstatus(GRR = 1, allele_freq = 0.02,
                                   RVfounder = "first", intro_RV = TRUE),
               list(c(0, 0), 1, T))
  expect_equal(sim_founderRVstatus(GRR = 1, allele_freq = 0.02,
                                   RVfounder = "first", intro_RV = FALSE),
               list(c(0, 0), 1, T))
  expect_equal(sim_founderRVstatus(GRR = 1, allele_freq = 0.02,
                                   RVfounder = "multiple", intro_RV = TRUE),
               list(c(0, 0), 1, T))
  expect_equal(sim_founderRVstatus(GRR = 1, allele_freq = 0.02,
                                   RVfounder = "multiple", intro_RV = FALSE),
               list(c(0, 0), 1, T))
})

test_that("RVfounder = first settings", {
  expect_equal(sim_founderRVstatus(GRR = 10, allele_freq = 0.02,
                                   RVfounder = "first", intro_RV = TRUE),
               list(c(0, 0), 1, T))
  first_ex <- sim_founderRVstatus(GRR = 10, allele_freq = 0.02,
                                  RVfounder = "first", intro_RV = FALSE)
  expect_true(1 %in% first_ex[[1]])
  expect_equal(first_ex[2:3],
               list(10, T))

})

test_that("RVfounder = multiple settings", {
  expect_equal(sim_founderRVstatus(GRR = 10, allele_freq = 1,
                                   RVfounder = "multiple", intro_RV = FALSE),
               list(c(1, 1), 10, F))
})
