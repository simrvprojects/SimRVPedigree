data("EgPeds")
context("trim.ped")

eggped <- ped(EgPeds)
test_that("returns an error when no proband selected", {
  expect_error(trim.ped(ped_file = eggped[eggped$FamID == 1, c(1,14)]))

})

test_that("The proband is an affected and experienced onset during the ascertainment span", {
  Tped <- trim.ped(ped_file = eggped[eggped$FamID == 1, ])
  expect_true(Tped$affected[Tped$proband == 1] == 1)
  expect_true(Tped$onsetYr[Tped$proband == 1] >= 2000 &
                Tped$onsetYr[Tped$proband == 1] <= 2015)
})
