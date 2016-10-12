library(testthat)
data("exp_peds")
context("trim_pedigree")
test_that("returns a pedfile with the is_proband variable", {
  expect_true(is.data.frame(trim_pedigree(ped_file = exp_peds[which(exp_peds$FamID == 1), ],
                                          ascertain_span = c(2000, 2015),
                                          num_affected = 2)))

  expect_true("is_proband" %in% colnames(trim_pedigree(ped_file = exp_peds[which(exp_peds$FamID == 1), ],
                                          ascertain_span = c(2000, 2015),
                                          num_affected = 2)))
})

test_that("The proband is an affected and experienced onset during the ascertainment span", {
  Tped <- trim_pedigree(ped_file = exp_peds[which(exp_peds$FamID == 1), ],
                        ascertain_span = c(2000, 2015),
                        num_affected = 2)
  expect_true(Tped$affected[which(Tped$is_proband == 1)] == 1)
  expect_true(Tped$onset_year[which(Tped$is_proband == 1)] >= 2000 &
                Tped$onset_year[which(Tped$is_proband == 1)] <= 2015)
})

