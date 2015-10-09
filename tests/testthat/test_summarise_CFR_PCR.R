context("Checking Model Summary Sanity")
preds <- CFR_PCR(summarised=TRUE)

test_that("Number of rows = conditions * items", {
  expect_equal(nrow(preds),45)

})

test_that("Accuracy range is sane", {
  r <- range(preds$acc)
  expect_true(r[1] >=0)
  expect_true(r[2] <=1)
})
