data("EgPeds")
context("trim.ped")

eggped <- ped(EgPeds)
test_that("returns an error when no proband selected", {
  expect_error(trim.ped(ped_file = eggped[eggped$FamID == 1, c(1,14)]))
})
