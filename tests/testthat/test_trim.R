data("EgPeds")
context("trim.ped")

eggped <- new.ped(EgPeds)
test_that("returns an error when no proband selected", {
  expect_error(trim.ped(eggped[eggped$FamID == 1, c(1,14)]))
})

test_that("returns a pedigree no larger than the pedigree that was supplied", {
  expect_lte(nrow(trim.ped(eggped[eggped$FamID == 1, ])),
             nrow(eggped[eggped$FamID == 1, ]))
})
