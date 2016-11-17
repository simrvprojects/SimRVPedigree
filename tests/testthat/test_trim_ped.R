library(testthat)
data("ExamplePedigrees")
context("trim_ped")
test_that("returns an error when no proband selected", {
  expect_error(trim_ped(ped_file = ExamplePedigrees[which(ExamplePedigrees$FamID == 1), c(1,14)]))

})

test_that("The proband is an affected and experienced onset during the ascertainment span", {
  Tped <- trim_ped(ped_file = ExamplePedigrees[which(ExamplePedigrees$FamID == 1), ])
  expect_true(Tped$affected[which(Tped$proband == 1)] == 1)
  expect_true(Tped$onset_year[which(Tped$proband == 1)] >= 2000 &
                Tped$onset_year[which(Tped$proband == 1)] <= 2015)
})
