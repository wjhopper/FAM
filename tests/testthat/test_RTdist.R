library(FAM)
context("Checking Likelihood Function Sanity")

preds <- CFR_PCL(summarised=FALSE)
dist <- preds %>% group_by(class,order) %>% RTdist()

test_that("Likelihood sums to the number of groups/classes", {
  expect_equal(sum(dist$y), length(unique(dist$class)))
  }
)
