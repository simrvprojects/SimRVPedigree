library(SimRVPedigree)
context("get_onsetHazard")
test_that("Always returns a vector of the appropriate length", {

  t_haz <- pgamma(seq(0, 99), shape = 15, scale = 4)/1000 +
    abs(rnorm(100, mean = 0, sd = 0.00005))

  expect_equal(length(get_onsetHazard(sub_hazard = t_haz,
                                      sub_GRR = 1, carrier_prob = 0.002,
                                      RV_status = T)),
               length(t_haz))

  expect_equal(length(get_onsetHazard(sub_hazard = t_haz,
                                      sub_GRR = 1, carrier_prob = 0.002,
                                      RV_status = F)),
               length(t_haz))

  expect_equal(length(get_onsetHazard(sub_hazard = t_haz,
                                      sub_GRR = 50, carrier_prob = 0.002,
                                      RV_status = T)),
               length(t_haz))


  expect_equal(length(get_onsetHazard(sub_hazard = t_haz,
                                      sub_GRR = 50, carrier_prob = 0.002,
                                      RV_status = F)),
               length(t_haz))

  expect_equal(length(get_onsetHazard(sub_hazard = t_haz,
                                      sub_GRR = 1, carrier_prob = 0,
                                      RV_status = T)),
               length(t_haz))

})


test_that("If GRR = 1 then hazard is unchanged", {

  t_haz <- pgamma(seq(0, 99), shape = 15, scale = 4)/1000 +
    abs(rnorm(100, mean = 0, sd = 0.00005))

  #for carriers
  expect_equal(get_onsetHazard(sub_hazard = t_haz,
                               sub_GRR = 1, carrier_prob = 0.002,
                               RV_status = T),
               t_haz)

  #for non-carriers
  expect_equal(get_onsetHazard(sub_hazard = t_haz,
                               sub_GRR = 1, carrier_prob = 0.002,
                               RV_status = F),
               t_haz)

})

test_that("If carrier_prob = 0 and GRR = 1 then hazard is unchanged", {

  t_haz <- pgamma(seq(0, 99), shape = 15, scale = 4)/1000 +
    abs(rnorm(100, mean = 0, sd = 0.00005))

  expect_equal(get_onsetHazard(sub_hazard = t_haz,
                               sub_GRR = 1, carrier_prob = 0,
                               RV_status = T),
               t_haz)

  expect_equal(get_onsetHazard(sub_hazard = t_haz,
                               sub_GRR = 1, carrier_prob = 0,
                               RV_status = F),
               t_haz)


})

test_that("If carrier_prob = 0, GRR != 1, and RV_status = FALSE then hazard is unchanged", {

  t_haz <- pgamma(seq(0, 99), shape = 15, scale = 4)/1000 +
    abs(rnorm(100, mean = 0, sd = 0.00005))

  expect_equal(get_onsetHazard(sub_hazard = t_haz,
                               sub_GRR = 50, carrier_prob = 0,
                               RV_status = F),
               t_haz)
})


test_that("GRR > 1 tests for carriers and non-carriers", {

  t_haz <- pgamma(seq(0, 99), shape = 15, scale = 4)/1000 +
    abs(rnorm(100, mean = 0, sd = 0.00005))

  base_rate <- t_haz/(1 + 0.002*(50 - 1))

  #non-carriers
  expect_equal(get_onsetHazard(sub_hazard = t_haz,
                               sub_GRR = 50, carrier_prob = 0.002,
                               RV_status = F),
               base_rate)

  #carriers
  expect_equal(get_onsetHazard(sub_hazard = t_haz,
                               sub_GRR = 50, carrier_prob = 0.002,
                               RV_status = T),
               base_rate*50)

})


test_that("If GRR > 1, hazard_rates increase for carriers and decrease for non-carriers", {

  t_haz <- pgamma(seq(0, 99), shape = 15, scale = 4)/1000 +
    abs(rnorm(100, mean = 0, sd = 0.00005))

  #carriers
  expect_true(all(get_onsetHazard(sub_hazard = t_haz,
                                  sub_GRR = 50, carrier_prob = 0.002,
                                  RV_status = T) > t_haz))

  #non-carriers
  expect_true(all(get_onsetHazard(sub_hazard = t_haz,
                                  sub_GRR = 50, carrier_prob = 0.002,
                                  RV_status = F) < t_haz))
})
